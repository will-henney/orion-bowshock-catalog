from __future__ import print_function
import glob
import json
import numpy as np
import ast

pattern = '../PC-will/*-arcdata.json'

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
col_names = ["Object", "RA", "Dec", "D", "D_binary",
             "PA_star", "h", "PA_out", "PA_in", 
             "R_out", "R_in", "Rc_out", "Rc_in",
             "Value_Bally", "Value_bg_Bally", 
             "Delta", "Dif_Bally"]
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
    table["PA_star"].append(PA_fmt(data["star"]["PA"]))
    table["h"].append(arcsec_fmt(data["thickness"]["h0"]))
    table["D"].append(arcsec_fmt(data["star"]["D"]))
    
    try:
        table["D_binary"].append(arcsec_fmt(data["star"]["D_binary"]))
    except KeyError:
        table["D_binary"].append("--")
    
    if "outer" in data:
        table["PA_out"].append(PA_fmt(data["outer"]["PA0"]))
        table["R_out"].append(arcsec_fmt(data["outer"]["R0"]))
        table["Rc_out"].append(arcsec_fmt(data["outer"]["Rc"]))
    else:
        table["PA_out"].append("--")
        table["R_out"].append("--")
        table["Rc_out"].append("--")
        
    if "inner" in data:
        table["PA_in"].append(PA_fmt(data["inner"]["PA0"]))
        table["R_in"].append(arcsec_fmt(data["inner"]["R0"]))
        table["Rc_in"].append(arcsec_fmt(data["inner"]["Rc"]))
    else:
        table["PA_in"].append("--")
        table["R_in"].append("--")
        table["Rc_in"].append("--")
        
    for k in data.keys():
        if k.startswith('Bally'):
            imagename=k
            
            try:
                shell = data[imagename]["shell"]["value"]
                bg = data[imagename]["background"]["value"]
                dbg = data[imagename]["background"]["delta"] 
            except KeyError:
                shell = 0.0
                bg = 0.0
                dbg = 0.0
            #diference between shell and backgraund for Bally 
            difference = shell - bg
            table["Value_Bally"].append(bright_fmt(shell)) 
            table["Dif_Bally"].append(bright_fmt(difference)) 
            table["Value_bg_Bally"].append(bright_fmt(bg))
            table["Delta"].append(bright_fmt(dbg))
            break               # Ensure we do not include more than one value
    else:
        # If no image was found, we still need to append to the lists
        table["Value_Bally"].append("--")
        table["Dif_Bally"].append("--") 
        table["Value_bg_Bally"].append("--")
        table["Delta"].append("--")
       
#Write output table to a file
with open("arcs-summary-PC-will.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))
       