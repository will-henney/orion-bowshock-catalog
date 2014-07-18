from __future__ import print_function
import glob
import json

pattern = "../../JorgeBowshocks/Jorge_prop/PC-correct/*/*-arcdata.json"

file_list = glob.glob(pattern)


def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table
        

# Initialize a list for each column in output table
col_names = ["Object", "RA", "Dec", "D", 
             "PA_star", "h", "PA_out", "PA_in", 
             "R_out", "R_in", "Rc_out", "Rc_in"]
table = {cn: [] for cn in col_names}


def PA_fmt(pa):
    """Write position angles to accuracy of 0.1 deg"""
    return "{:.1f}".format(pa)


def arcsec_fmt(r):
    """ Write distances to accuracy of 0.001 arcsec"""
    return "{:.3f}".format(r)

def bright_fmt(r):
    """ Write brightnesses to accuracy of 0.001"""
    return "{:.3f}".format(r)


for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)

    # Add this object's data to each output column
    table["Object"].append(data["star"]["id"])
    table["RA"].append(data["star"]["RA"])
    table["Dec"].append(data["star"]["Dec"])
    table["D"].append(arcsec_fmt(data["star"]["D"]))
    table["PA_star"].append(PA_fmt(data["star"]["PA"]))
    table["h"].append(arcsec_fmt(data["thickness"]["h0"]))
    
    if "outer" in data:
        table["PA_out"].append(PA_fmt(data["outer"]["PA0"]))
        table["R_out"].append(arcsec_fmt(data["outer"]["R0"]))
        table["Rc_out"].append(arcsec_fmt(data["outer"]["Rc"]))
    else:
        table["PA_out"].append("-")
        table["R_out"].append("-")
        table["Rc_out"].append("-")
        
    if "inner" in data:
        table["PA_in"].append(PA_fmt(data["inner"]["PA0"]))
        table["R_in"].append(arcsec_fmt(data["inner"]["R0"]))
        table["Rc_in"].append(arcsec_fmt(data["inner"]["Rc"]))
    else:
        table["PA_in"].append("-")
        table["R_in"].append("-")
        table["Rc_in"].append("-")
        


# Write output table to a file
with open("arcs-summary-proplyds.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))
    
