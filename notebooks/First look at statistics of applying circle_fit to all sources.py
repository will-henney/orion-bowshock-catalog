# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import sys
from pathlib import Path
from astropy.table import Table
import pandas as pd
import numpy as np
import json

ROOT_PATH = Path("..") # necessary since we are in the `notebooks/` sub-folder
# -

# Read the table of sources and convert to pandas dtaframe, using a query to exclude problematic and interproplyd sources:

qstring = "~((Group == '*P') | (Group == '*I'))"
df = Table.read(
    ROOT_PATH / "luis-programas/arcs-summary-classify.ecsv",
    format="ascii.ecsv").to_pandas().query(qstring)


df.columns

# Make a new dataframe that just keeps the planitude and alatude from the old calculations, and then adds columns for the new `circle_fit` results.

df2 = df[["Object", "Group", "D", "PA_star", "h"]].assign(
    old_Pi_in=df.Rc_in/df.R_in,
    old_Pi_out=df.Rc_out/df.R_out,
)

# Look at the mean and standard deviation by spatial group.

df2.groupby(by="Group").mean().sort_values(by="D")

df2.groupby(by="Group").std().sort_values(by="D")

# +
# Keys in JSON file for each variable
VARS = "Pi", "Lambda", "d Lambda", "R0"#, "d theta"
# And a parallel version for the table column names (no spaces!)
VVARS = [_.replace(" ", "_") for _ in VARS]


def mad(x):
    """Median absolute deviation rescaled to sigma for Gaussian"""
    return 1.4826*np.nanmedian(np.abs(x - np.nanmedian(x)))


def collate_circle_fit_one_source(source_id):
    rslt = {"Object": source_id}
    
    for arc_long, arc in [["inner", "in"], ["outer", "out"]]:
        json_paths = ROOT_PATH.glob(f"*/{source_id}-{arc_long}-*.json")

        # Initialize empty lists for each variable
        var_lists = {f"{vv}_{arc}": [] for vv in VVARS}

        # Now fill in those lists with data from JSON files
        for json_path in json_paths:
            data = json.load(json_path.open())
            if not 60 <= data["d theta"] <= 75:
                continue # skip anything outside this range of d theta
            for v, vv in zip(VARS, VVARS):
                var_lists[f"{vv}_{arc}"].append(data.get(v))

        # Finally, calculate median and dispersion
        for vv in VVARS:
            key = f"{vv}_{arc}"
            vmedian = np.nanmedian(var_lists[key])
            vsig = mad(np.array(var_lists[key]))
            rslt[key] = vmedian
            rslt[f"e_{key}"] = vsig
            rslt[f"n_{key}"] = len(var_lists[key])
            # rslt[f"list_{key}"] = var_lists[key]

    return rslt



# -

data = collate_circle_fit_one_source("LL2")
data

df3 = pd.DataFrame([collate_circle_fit_one_source(source_id) for source_id in df2["Object"]])
df3

dff = df2.merge(df3, how="outer", on="Object")
dff

columns_to_drop = [_ for _ in dff.columns if _.startswith("n_")]
dff.drop(columns_to_drop, axis=1, inplace=True)

dff.groupby(by="Group").mean().sort_values(by="D")

import seaborn as sn

dff_log = pd.concat([
    dff[["Object", "Group"]],
    dff[["old_Pi_in", "old_Pi_out", "Pi_in", "Pi_out", "Lambda_in", "Lambda_out"]].clip(lower=1.0, upper=15.0).apply(np.log10),
], axis=1)


dff_log

fig = sn.pairplot(dff_log, hue="Group", palette="magma")
for i, row in enumerate(fig.axes):
    for j, ax in enumerate(row):
        if i != j:
            ax.set(xlim=[-0.2, 1.2], ylim=[-0.2, 1.2])
            ax.plot([-0.2, 1.2], [-0.2, 1.2], "--", c="k")

ax = sn.jointplot(dff["Pi_out"], dff["Lambda_out"], kind="reg")
ax.ax_joint.set(xlim=[0.9, 15.0], ylim=[0.9, 15.0], xscale="log", yscale="log")
ax.ax_joint.set_aspect("equal")

ax = sn.jointplot(dff["Pi_in"], dff["Lambda_in"], kind="reg", color="g")
ax.ax_joint.set(xlim=[0.9, 15.0], ylim=[0.9, 15.0], xscale="log", yscale="log")
ax.ax_joint.set_aspect("equal")

from matplotlib import pyplot as plt

sn.set_context("poster")

from scipy.stats import gaussian_kde 
from matplotlib.colors import PowerNorm

# Since we are plotting $\Pi$–$\Lambda$ in logarithmic space, we want to calculate the KDE in log space also.  Otherwise, the smoothing will seem to be larger for lower values of $\Pi$, $\Lambda$. 

# +
import seaborn.distributions 

# Save original version of KDE function
smkde_original = seaborn.distributions._statsmodels_bivariate_kde
def smkde_logxy(x, y, bw, gridsize, cut, clip):
    """Calculate KDE on logarithmic grid for x and y"""
    xx, yy, z = smkde_original(np.log10(x), np.log10(y), bw, gridsize, cut, clip)
    xx = 10**xx
    yy = 10**yy
    return xx, yy, z
# Monkey patch the function that sns.kdeplot uses to calculate the KDE
seaborn.distributions._statsmodels_bivariate_kde = smkde_logxy


# +
fig, ax = plt.subplots(figsize=(12, 12))

Rc_grid = np.linspace(0.0, 15.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)


m = np.isfinite(dff["Pi_out"]) & np.isfinite(dff["Lambda_out"])
sn.kdeplot(dff[m]["Pi_out"], dff[m]["Lambda_out"], ax=ax, 
           bw=(0.07, 0.035), levels=[0.5, 1.0, 2.0, 4.0, 6.0], norm=PowerNorm(0.5),
           cmap="Blues",
          )

ax.errorbar(x=dff["Pi_out"], y=dff["Lambda_out"], xerr=dff["e_Pi_out"], yerr=dff["d_Lambda_out"], 
            fmt='none', alpha=0.3)

groups = [ "LV", "SE", "N", "NW", "W", "S", ]
colors = [ "light brown", "light orange", "cerulean", "dark pink", "purple", 
           "forest green", ]
colors = sn.xkcd_palette(colors)
groups_and_colors = list(zip(groups, colors))


for group, color in groups_and_colors:
    m = dff["Group"] == group
    ax.scatter(x=dff[m]["Pi_out"], y=dff[m]["Lambda_out"], 
               marker="o", c=[color], label=group, edgecolors="w")
ax.legend(ncol=2).set_title("Group")
ax.set(xlim=[0.35, 15.0], ylim=[0.35, 15.0], xscale="log", yscale="log")
ax.set_aspect("equal")
ax.set(xlabel=r"$\Pi_\mathrm{out}$", ylabel=r"$\Lambda_\mathrm{out}$")

fig.savefig(ROOT_PATH / "all-sources-Pi_out-Lambda_out.pdf")

# +
fig, ax = plt.subplots(figsize=(12, 12))


Rc_grid = np.linspace(0.0, 15.0, 2000)
R90_T0_grid = np.sqrt(2*Rc_grid)
R90_T1_grid = np.sqrt(2*Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0 

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color='k', alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color='k', alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c='k', lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color='k', zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color='k', zorder=-1)


m = np.isfinite(dff["Pi_in"]) & np.isfinite(dff["Lambda_in"])
sn.kdeplot(dff[m]["Pi_in"], dff[m]["Lambda_in"], ax=ax, 
           bw=(0.07, 0.035), levels=[0.5, 1.0, 2.0, 4.0, 8.0, 12.0], norm=PowerNorm(0.5),
           cmap="Purples",
          )



ax.errorbar(x=dff["Pi_in"], y=dff["Lambda_in"], xerr=dff["e_Pi_in"], yerr=dff["d_Lambda_in"], 
            fmt='none', alpha=0.3)

for group, color in groups_and_colors:
    m = dff["Group"] == group
    ax.scatter(x=dff[m]["Pi_in"], y=dff[m]["Lambda_in"], 
               marker="o", c=[color], label=group, edgecolors="w")
ax.legend(ncol=2).set_title("Group")    
ax.set(xlim=[0.35, 15.0], ylim=[0.35, 15.0], xscale="log", yscale="log")
ax.set_aspect("equal")
ax.set(xlabel=r"$\Pi_\mathrm{in}$", ylabel=r"$\Lambda_\mathrm{in}$")
fig.savefig(ROOT_PATH / "all-sources-Pi_in-Lambda_in.pdf")
# -

groups_and_colors


