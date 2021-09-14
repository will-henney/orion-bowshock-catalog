import sys
from pathlib import Path
from astropy.table import Table
import pandas as pd
import numpy as np
import json
import seaborn as sn

ROOT_PATH = Path("..")

qstring = "~((Group == '*P') | (Group == '*I'))"
df = (
    Table.read(
         "arcs-summary-classify.ecsv", format="ascii.ecsv"
    )
    .to_pandas())
    #.query(qstring)
    #)

print(df.columns)

# Make a new dataframe that just keeps the planitude and alatude from the old calculations, and then adds columns for the new circle_fit results.

df2 = df[["Object", "Group", "D", "PA_star", "h"]].assign(
    old_Pi_in=df.Rc_in / df.R_in, old_Pi_out=df.Rc_out / df.R_out,
)

# Look at the mean and standard deviation by spatial group.

print(df2.groupby(by="Group").mean().sort_values(by="D"))
print(df2.groupby(by="Group").std().sort_values(by="D"))

# Keys in JSON file for each variable
VARS = "Pi", "Lambda", "d Lambda", "R0"  # , "d theta"
# And a parallel version for the table column names (no spaces!)
VVARS = [_.replace(" ", "_") for _ in VARS]


def mad(x):
    """Median absolute deviation rescaled to sigma for Gaussian"""
    return 1.4826 * np.nanmedian(np.abs(x - np.nanmedian(x)))


def collate_circle_fit_one_source(source_id):
    rslt = {"Object": source_id}

    for arc_long, arc in [["inner", "in"], ["outer", "out"], ["ridge", "out"]]:
        json_paths = sorted(ROOT_PATH.glob(f"*/{source_id}-{arc_long}-*.json"))
        try:
            rslt["Folder"] = str(json_paths[0].parent.name)
        except IndexError:
            rslt["Folder"] = ""
        # print(arc_long)
        # print(json_paths)
        # Initialize empty lists for each variable
        var_lists = {f"{vv}_{arc}": [] for vv in VVARS}

        # Now fill in those lists with data from JSON files
        # print(json_paths)
        if not json_paths:
            # But skip the rest if there is no data for this arc type
            continue

        rslt["Outer arc type"] = arc_long
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


data = collate_circle_fit_one_source("167-317")
print(data)

df3 = pd.DataFrame(
    [collate_circle_fit_one_source(source_id) for source_id in df2["Object"]]
)

print(df3)

dff = df2.merge(df3, how="outer", on="Object")

columns_to_drop = [_ for _ in dff.columns if _.startswith("n_")]
other_columns_to_drop = ['D', 'PA_star', 'h', 'R0_in', 'e_R0_in', 'R0_out', 'e_R0_out']
total_columns_to_drop = columns_to_drop + other_columns_to_drop

dff.drop(total_columns_to_drop, axis=1, inplace=True)

print(type(dff))
print(dff.columns)

# Write table with the results
Table.from_pandas(dff).write('arcs-sources-Pi-Lambda.ecsv', format='ascii.ecsv', overwrite=True)
