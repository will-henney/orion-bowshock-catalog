# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
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
import seaborn as sn

ROOT_PATH = Path("..")  # necessary since we are in the `notebooks/` sub-folder
# -

# Read the table of sources and convert to pandas dtaframe, using a query to exclude problematic and interproplyd sources:

qstring = "~((Group == '*P') | (Group == '*I'))"
df = (
    Table.read(
        ROOT_PATH / "luis-programas/arcs-summary-classify.ecsv", format="ascii.ecsv"
    )
    .to_pandas()
    .query(qstring)
)


df.columns

# Make a new dataframe that just keeps the planitude and alatude from the old calculations, and then adds columns for the new `circle_fit` results.

df2 = df[["Object", "Group", "D", "PA_star", "h"]].assign(
    old_Pi_in=df.Rc_in / df.R_in, old_Pi_out=df.Rc_out / df.R_out,
)

# Look at the mean and standard deviation by spatial group.

df2.groupby(by="Group").mean().sort_values(by="D")

df2.groupby(by="Group").std().sort_values(by="D")

# +
# Keys in JSON file for each variable
VARS = "Pi", "Lambda", "d Lambda", "R0"  # , "d theta"
# And a parallel version for the table column names (no spaces!)
VVARS = [_.replace(" ", "_") for _ in VARS]


def mad(x):
    """Median absolute deviation rescaled to sigma for Gaussian"""
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def collate_circle_fit_one_source(source_id):
    rslt = {"Object": source_id}

    for arc_long, arc in [["inner", "in"], ["outer", "out"]]:
        json_paths = sorted(ROOT_PATH.glob(f"*/{source_id}-{arc_long}-*.json"))
        try:
            rslt["Folder"] = str(json_paths[0].parent.name)
        except IndexError:
            rslt["Folder"] = ""

        # Initialize empty lists for each variable
        var_lists = {f"{vv}_{arc}": [] for vv in VVARS}

        # Now fill in those lists with data from JSON files
        # print(json_paths)
        for json_path in json_paths:
            data = json.load(json_path.open())
            # print(data)
            if not 60 <= data["d theta"] <= 75:
                continue  # skip anything outside this range of d theta
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

data = collate_circle_fit_one_source("LL7")
data

df3 = pd.DataFrame(
    [collate_circle_fit_one_source(source_id) for source_id in df2["Object"]]
)
df3

dff = df2.merge(df3, how="outer", on="Object")
dff

columns_to_drop = [_ for _ in dff.columns if _.startswith("n_")]
dff.drop(columns_to_drop, axis=1, inplace=True)

dff.groupby(by="Group").mean().sort_values(by="D")

dff_log = pd.concat(
    [
        dff[["Object", "Group"]],
        dff[["old_Pi_in", "old_Pi_out", "Pi_in", "Pi_out", "Lambda_in", "Lambda_out"]]
        .clip(lower=1.0, upper=15.0)
        .apply(np.log10),
    ],
    axis=1,
)


dff_log

fig = sn.pairplot(dff_log.dropna(), hue="Group", palette="magma")
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
    xx = 10 ** xx
    yy = 10 ** yy
    return xx, yy, z


# Monkey patch the function that sns.kdeplot uses to calculate the KDE
seaborn.distributions._statsmodels_bivariate_kde = smkde_logxy


# +
fig, ax = plt.subplots(figsize=(12, 12))

Rc_grid = np.linspace(0.0, 15.0, 2000)
R90_T0_grid = np.sqrt(2 * Rc_grid)
R90_T1_grid = np.sqrt(2 * Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color="k", alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color="k", alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c="k", lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color="k", zorder=-1)


m = np.isfinite(dff["Pi_out"]) & np.isfinite(dff["Lambda_out"])
sn.kdeplot(
    dff[m]["Pi_out"],
    dff[m]["Lambda_out"],
    ax=ax,
    bw=(0.07, 0.035),
    levels=[0.5, 1.0, 2.0, 4.0, 6.0],
    norm=PowerNorm(0.5),
    cmap="Blues",
)

ax.errorbar(
    x=dff["Pi_out"],
    y=dff["Lambda_out"],
    xerr=dff["e_Pi_out"],
    yerr=dff["d_Lambda_out"],
    fmt="none",
    alpha=0.3,
)

groups = [
    "LV",
    "SE",
    "N",
    "NW",
    "W",
    "S",
]
colors = [
    "light brown",
    "light orange",
    "cerulean",
    "dark pink",
    "purple",
    "forest green",
]
colors = sn.xkcd_palette(colors)
groups_and_colors = list(zip(groups, colors))


for group, color in groups_and_colors:
    m = dff["Group"] == group
    ax.scatter(
        x=dff[m]["Pi_out"],
        y=dff[m]["Lambda_out"],
        marker="o",
        c=[color],
        label=group,
        edgecolors="w",
    )
ax.legend(ncol=2).set_title("Group")
ax.set(xlim=[0.35, 15.0], ylim=[0.35, 15.0], xscale="log", yscale="log")
ax.set_aspect("equal")
ax.set(xlabel=r"$\Pi_\mathrm{out}$", ylabel=r"$\Lambda_\mathrm{out}$")

fig.savefig(ROOT_PATH / "all-sources-new-Pi_out-Lambda_out.pdf")

# +
fig, ax = plt.subplots(figsize=(12, 12))


Rc_grid = np.linspace(0.0, 15.0, 2000)
R90_T0_grid = np.sqrt(2 * Rc_grid)
R90_T1_grid = np.sqrt(2 * Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color="k", alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color="k", alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c="k", lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color="k", zorder=-1)


m = np.isfinite(dff["Pi_in"]) & np.isfinite(dff["Lambda_in"])
sn.kdeplot(
    dff[m]["Pi_in"],
    dff[m]["Lambda_in"],
    ax=ax,
    bw=(0.07, 0.035),
    levels=[0.5, 1.0, 2.0, 4.0, 8.0, 12.0],
    norm=PowerNorm(0.5),
    cmap="Purples",
)


ax.errorbar(
    x=dff["Pi_in"],
    y=dff["Lambda_in"],
    xerr=dff["e_Pi_in"],
    yerr=dff["d_Lambda_in"],
    fmt="none",
    alpha=0.3,
)

for group, color in groups_and_colors:
    m = dff["Group"] == group
    ax.scatter(
        x=dff[m]["Pi_in"],
        y=dff[m]["Lambda_in"],
        marker="o",
        c=[color],
        label=group,
        edgecolors="w",
    )
ax.legend(ncol=2).set_title("Group")
ax.set(xlim=[0.35, 15.0], ylim=[0.35, 15.0], xscale="log", yscale="log")
ax.set_aspect("equal")
ax.set(xlabel=r"$\Pi_\mathrm{in}$", ylabel=r"$\Lambda_\mathrm{in}$")
fig.savefig(ROOT_PATH / "all-sources-new-Pi_in-Lambda_in.pdf")
# -

groups_and_colors

# Look at which sources to not have valid outer alatude, $\Lambda_\mathrm{out}$

" ".join(dff[~np.isfinite(dff["Lambda_out"])]["Object"].to_list())

" ".join(dff[~np.isfinite(dff["Lambda_in"])]["Object"].to_list())

" ".join(
    dff[(np.isfinite(dff["Lambda_out"])) & (~np.isfinite(dff["d_Lambda_out"]))][
        "Object"
    ].to_list()
)

" ".join(
    dff[(np.isfinite(dff["Lambda_in"])) & (~np.isfinite(dff["d_Lambda_in"]))][
        "Object"
    ].to_list()
)


def flag(x):
    return np.char.array(np.where(np.isfinite(x), "🟨", "⬜️"))


flag(dff.Lambda_out) + flag(dff.d_Lambda_out) + flag(dff.Lambda_in) + flag(
    dff.d_Lambda_in
)

dff["Lam_flags"] = (
    flag(dff.Lambda_out)
    + flag(dff.d_Lambda_out)
    + flag(dff.Lambda_in)
    + flag(dff.d_Lambda_in)
)

pd.options.display.max_rows = 999

dff[["Object", "Group", "Folder", "Lam_flags"]]

dff[["Object", "Group", "Folder", "Lam_flags"]].to_csv(
    ROOT_PATH / "missing_lambdas.csv", index=False
)

Path("xx/yy/zz.txt").parent

# ## Check for significant differences between inner/outer shapes
#
# Use K-S or Anderson-Darling

m = dff["Pi_in"] < 20.0
dff["log_Pi_in"] = np.log10(dff["Pi_in"])
dff["log_Pi_out"] = np.log10(dff["Pi_out"])
ax = sn.jointplot(dff[m]["log_Pi_out"], dff[m]["log_Pi_in"], kind="reg")
ax.ax_joint.set(
    xticks=[0, 0.5, 1.0],
    yticks=[0, 0.5, 1.0],
    #    xlim=[0.9, 15.0], ylim=[0.9, 15.0],
    #    xscale="log", yscale="log",
)
ax.ax_joint.set_aspect("equal")

ax = sn.jointplot(dff["Lambda_out"], dff["Lambda_in"], kind="reg")
# ax.ax_joint.set(xlim=[0.9, 15.0], ylim=[0.9, 15.0], xscale="log", yscale="log")
ax.ax_joint.set(
    xticks=[0, 1, 2, 3, 4, 5], yticks=[0, 1, 2, 3, 4, 5],
)
ax.ax_joint.set_aspect("equal")

# Import all the stat tests ...

from astropy.stats import kuiper_two
from scipy.stats import ks_2samp, anderson_ksamp, pearsonr, levene


def print_stats(var1, var2):
    m = np.isfinite(dff[var1]) & np.isfinite(dff[var2])
    x, y = dff[m][var1], dff[m][var2]
    _, p_kuiper = kuiper_two(x, y)
    _, p_KS = ks_2samp(x, y)
    _, _, p_AD = anderson_ksamp([x, y], midrank=False)
    _, p_BF = levene(x, y, center="median")
    r_P, p_P = pearsonr(x, y)
    print(
        f"Pearson rank correlation between {var1} and {var2}: {r_P:.5f} (p = {p_P:.5f})"
    )
    print(f"P values for tests of difference between {var1} and {var2}:")
    print("Kuiper test:", np.round(p_kuiper, 5))
    print("Kolmogorov–Smirnov test:", np.round(p_KS, 5))
    print("Anderson–Darling test:", np.round(p_AD, 5))
    print("Brown–Forsyth test (variance):", np.round(p_BF, 5))


print_stats("Pi_out", "Pi_in")

print_stats("Pi_out", "Lambda_out")

print_stats("Lambda_out", "Lambda_in")

print_stats("Pi_in", "Lambda_in")

print_stats("Pi_out", "old_Pi_out")

print_stats("Pi_in", "old_Pi_in")

print_stats("log_Pi_out", "log_Pi_in")

# Next one is a but silly.  Of course, the difference tests will be significant between one value that is logged and one that isn't ...

print_stats("log_Pi_out", "Lambda_out")
print()
print_stats("log_Pi_in", "Lambda_in")

# ### Summary of statistical tests
#
# We use two different types of tests:
#
# 1. Pearson rank correlation, which ignores the absolute values, but looks at between-source relative change in the two variables.
#
# 2. The difference tests, which *are* sensitive to the absolute values but ignore only look at the two distributions, ignoring the correlations between sources.  Kuiper, K–S, and A–D are mainly sensitive to central tendency, while B–F is sensitive to width.
#
# Prinipal findings are:
#
# * Strong correlation between `Lambda_out` and `Lambda_in`: $r = 0.61$, $p = 0.0005$.  All difference tests non-significant.
# * Weaker (but still significant) correlation between `log_Pi_out` and `log_Pi_in`: $r = 0.39$, $p = 0.002$.  Again, all difference tests non-significant.
# * `Lambda` versus `Pi` has strong correlation (both out and in).  The Brown–Forsyth test is the only one that shows a significant difference, presumably because `Lambda` does not have the fat tail towards high values that `Pi` has.
#
# Conclusions are that it is fine to merge the outer and inner arc values, since there is no significant difference between the two.

# ## Merging the inner and outer arcs
#
# We will do this by the following steps
#
# 1. Separate out the inner and outer arc columns into two separate data frames
# 2. Rename `Pi`, `Lambda` (and `d_...`) columns in both to remove `_in` and `_out` suffixes.
# 3. Add an `Arc` column to each (value `in` or `out`)
# 4. Concatenate the tables.

VVARS


# Steps 1, 2, 3, for outer arcs:

# +
def strip_suffix(s):
    return s.replace("_out", "").replace("_in", "")


COLS = (
    ["Object", "Group", "Folder"]
    + [f"{v}_out" for v in VVARS]
    + [f"e_{v}_out" for v in VVARS]
)
dff_out = dff[COLS].rename(strip_suffix, axis="columns")
dff_out["Arc"] = "out"
dff_out
# -

# Steps 1, 2, 3, for inner arcs:

COLS = (
    ["Object", "Group", "Folder"]
    + [f"{v}_in" for v in VVARS]
    + [f"e_{v}_in" for v in VVARS]
)
dff_in = dff[COLS].rename(strip_suffix, axis="columns")
dff_in["Arc"] = "in"
dff_in

# Step 4 to put them both together.

dff_merge = pd.concat([dff_out, dff_in], ignore_index=True)
dff_merge

# Now we can plot the merged dataset as before

# +
fig, ax = plt.subplots(figsize=(12, 12))


Rc_grid = np.linspace(0.0, 15.0, 2000)
R90_T0_grid = np.sqrt(2 * Rc_grid)
R90_T1_grid = np.sqrt(2 * Rc_grid - 1.0)
R90_T1_grid[~np.isfinite(R90_T1_grid)] = 0.0

ax.fill_between(Rc_grid, R90_T1_grid, R90_T0_grid, color="k", alpha=0.2)
ax.fill_between(Rc_grid, R90_T0_grid, color="k", alpha=0.1)
ax.plot(Rc_grid, R90_T0_grid, c="k", lw=0.5)
ax.axhline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.axvline(1.0, lw=0.5, alpha=0.5, color="k", zorder=-1)
ax.plot([0.0, 10.0], [0.0, 10.0], lw=0.5, alpha=0.5, color="k", zorder=-1)


m = np.isfinite(dff_merge["Pi"]) & np.isfinite(dff_merge["Lambda"])
sn.kdeplot(
    dff_merge[m]["Pi"],
    dff_merge[m]["Lambda"],
    ax=ax,
    bw=(0.05, 0.025),
    levels=[0.5, 1.0, 2.0, 4.0, 8.0, 12.0],
    norm=PowerNorm(0.5),
    cmap="Purples",
    zorder=1,
)


ax.errorbar(
    x=dff_merge["Pi"],
    y=dff_merge["Lambda"],
    xerr=dff_merge["e_Pi"],
    yerr=dff_merge["d_Lambda"],
    fmt="none",
    alpha=0.3,
    zorder=2,
)

sizes = np.empty_like(dff_merge["Pi"])
m = dff_merge["Arc"] == "out"
sizes[m] = 200.0
sizes[~m] = 100.0

for group, color in groups_and_colors:
    m = dff_merge["Group"] == group
    ax.scatter(
        x=dff_merge[m]["Pi"],
        y=dff_merge[m]["Lambda"],
        marker="o",
        c=[color],
        s=sizes[m],
        label=group,
        edgecolors="w",
        zorder=3,
    )
ax.legend(ncol=2).set_title("Group")
ax.set(xlim=[0.35, 15.0], ylim=[0.35, 15.0], xscale="log", yscale="log")
ax.set_aspect("equal")
ax.set(xlabel=r"$\Pi$", ylabel=r"$\Lambda$")
fig.savefig(ROOT_PATH / "all-sources-merge-Pi-Lambda.pdf")
# -

# This looks very promising.  We now have enough points that we can reduce the bandwidth (smoothing) of the KDE.
#
# It appears that the shape distribution is bimodal (at least).  There is a cluster at $(\Pi, \Lambda) = (2, 2)$ and another at $(\Pi, \Lambda) \approx (3, 3)$.  Then some outliers around $(1.4, 1.4)$ and other outliers with $\Pi > 5$.
#
# *We should model this with a Gaussian Mixture Model when we have more data.*  (See Jake VdP book)

# ## Fit Gaussian Cluster Model to shape data

from sklearn.mixture import GaussianMixture as GMM

# Put data in form expected by scikit-learn (and without NaNs)

m = np.isfinite(dff_merge["Pi"]) & np.isfinite(dff_merge["Lambda"])
X = np.array(list(zip(np.log10(dff_merge[m]["Pi"]), np.log10(dff_merge[m]["Lambda"]))))

gmm = GMM(n_components=4).fit(X)

labels = gmm.predict(X)

labels

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(10 ** X[:, 0], 10 ** X[:, 1], c=labels, s=40, cmap="rainbow")
ax.set(
    xscale="log",
    yscale="log",
    xlim=[0.75, 15],
    ylim=[0.75, 15],
    xlabel=r"Planitude, $\Pi$",
    ylabel=r"Alatude, $\Lambda$",
)
ax.set_aspect("equal")

# So, using 4 components, it prefers to not separate the (1.4, 1.4) clump, instead fitting one component to the spur, which is probably spurious. It also fits a component to the high-$\Pi$ tail, which is more justified.

# ### Bayesian Gaussian Mixture to automatically decide how many components

from sklearn.mixture import BayesianGaussianMixture

bgmm = BayesianGaussianMixture(n_components=4, verbose=1).fit(X)

bgmm.weights_

blabels = bgmm.predict(X)

blabels

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(10 ** X[:, 0], 10 ** X[:, 1], c=blabels, s=40, cmap="rainbow")
ax.set(
    xscale="log",
    yscale="log",
    xlim=[0.75, 15],
    ylim=[0.75, 15],
    xlabel=r"Planitude, $\Pi$",
    ylabel=r"Alatude, $\Lambda$",
)
ax.set_aspect("equal")

# So, using the default parameters means that it decides that two components are enough.  One (purple) gets the main ridge with \(\Lambda \approx \Pi\), while the other gets the outliers.

bgmm = BayesianGaussianMixture(
    n_components=30, verbose=1, weight_concentration_prior=1000.0
).fit(X)

bgmm.weights_

blabels = bgmm.predict(X)
blabels

fig, ax = plt.subplots(figsize=(8, 8))
ax.scatter(10 ** X[:, 0], 10 ** X[:, 1], c=blabels, s=40, cmap="rainbow")
ax.set(
    xscale="log",
    yscale="log",
    xlim=[0.75, 15],
    ylim=[0.75, 15],
    xlabel=r"Planitude, $\Pi$",
    ylabel=r"Alatude, $\Lambda$",
)
ax.set_aspect("equal")

# However much we fiddle with the parameters to try and favor multiple components, it only wants to fit one component to the main diagonal.  Although it can be persuaded to split the outliers up into different components.

# ### Use pomegranate library to incorporate per-sample weights
#
# With scikit.learn, we can't have different weights for each source (sample in ML parlance).  But we can do that with the pomegranate library.

from pomegranate import GeneralMixtureModel, MultivariateGaussianDistribution

rel_err = np.hypot(
    dff_merge[m]["e_Pi"] / dff_merge[m]["Pi"],
    dff_merge[m]["e_Lambda"] / dff_merge[m]["Lambda"],
)
weights = 1.0 / rel_err ** 2
weights[rel_err == 0.0] = np.nanmedian(weights)
max_weight = 4.0 * np.nanmedian(weights)
weights[weights > max_weight] = max_weight
weights[::10]

N_COMPONENTS = 4

model = GeneralMixtureModel.from_samples(
    distributions=MultivariateGaussianDistribution,
    n_components=N_COMPONENTS,
    X=X,
    weights=weights,
)

plabels = model.predict(X)
plabels

xmin, xmax = 0.75, 15
ymin, ymax = 0.75, 15
NP = 1000
xx, yy = np.meshgrid(np.linspace(xmin, xmax, NP), np.linspace(ymin, ymax, NP),)
x_ = np.array(list(zip(xx.flatten(), yy.flatten())))
p_components = model.predict_proba(np.log10(x_)).reshape((NP, NP, N_COMPONENTS))
p_tot = model.probability(np.log10(x_)).reshape((NP, NP))
p_tot /= p_tot.max()

fig, ax = plt.subplots(figsize=(8, 8))
# for k in range(1, N_COMPONENTS):
#    ax.contour(xx, yy, p_components[:, :, k], levels=[0.5], cmap="Reds_r", zorder=-2)
levels = np.array([0.001, 0.003, 0.01, 0.03, 0.1, 1.0])
GAMMA = 0.05
ax.contourf(
    xx, yy, p_tot ** GAMMA, levels=levels ** GAMMA, cmap="Greens", zorder=-1, alpha=0.5
)
scatter = ax.scatter(
    10 ** X[:, 0],
    10 ** X[:, 1],
    c=plabels,
    s=5 * np.sqrt(weights),
    edgecolors="w",
    linewidths=1,
    vmin=-1,
    vmax=4,
    cmap="seismic",
)
ax.set(
    xscale="log",
    yscale="log",
    xlim=[xmin, xmax],
    ylim=[ymin, ymax],
    xlabel=r"Planitude, $\Pi$",
    ylabel=r"Alatude, $\Lambda$",
)
ax.legend(*scatter.legend_elements(), ncol=2, title="Component")
ax.set_aspect("equal")

# So it turns out that we *can* get the components that we could intuitively see, by means of using the weights.
#
# However, this is very sensitive to initial conditions.  Re-running the fit can give wildly different decompositions.

for p in range(N_COMPONENTS):
    print("COMPONENT", p)
    print(weights[plabels == p].describe())
    print()

p_components[:, :, 0].mean()

# ### Correlation between spatial groups and GMM components
#
#

# Select just the sources that have $\Pi$ and $\Lambda$ measurements.

df = dff_merge[m].copy()

df["GMM"] = plabels
df["Weight"] = weights
df

sn.countplot(x="Group", hue="GMM", data=df, palette="seismic")
plt.gcf().set_size_inches(8, 6)
sn.despine()

sn.countplot(x="GMM", hue="Group", data=df, palette="magma")
plt.gcf().set_size_inches(8, 6)
sn.despine()

# *Note that the labels of the GMM components may have changed since I wrote the following cell*

# So there are some slight differences.
#
# The NW group has the majority of component 3, which is the high-$\Lambda$ tail. This is only 3 objects though, so may not be significant.
#
# Component 1, which is the low-$(\Pi, \Lambda)$ clump, has a slight over-representation of the S group, but this is omly 2 objects, so even less significant.
#
# Component 2, which is the compact clump at $(2, 2)$,  also has over-representation of the S group, with 6 sources, which is more significant.
#
# Conponent 0, which is the main diagonal ridge, is dominated by the W group.

# Token change to test jupytext


