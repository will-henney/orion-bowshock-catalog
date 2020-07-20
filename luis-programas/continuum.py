'''
Estmate the contribution of the continuum
'''

import glob
import json
import numpy as np

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

# The normalization is to the exposure times in the program GO 5085  
T_f656 = 200 #s
T_f658 = 500
T_f547 = 50
T_acs_f658 = 1
T_robb_f656 = 1

# Rectangular width
w_f656 = 28.34 #Angstrom
w_f658 = 39.23
w_f547 = 638.10
w_acs_f658 = 74.9
w_robb_f656 = 28.34

shell = {}
bg = {}
label_dict = {"acs":"Bally", "wfpc2":"Robberto_WFPC2", "wfpc2_f658":"WFC_mosaic_f658", 
"wfpc2_f656":"WFC_mosaic_f656", "wfpc2_f547":"WFC_mosaic_f547"}
label_s = []
label_ = {}

nsources = len(file_list)

# allocate empty arrays to store the mean and sigma
shell = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,)), 'wfpc2_f658': np.empty((nsources,)), 
'wfpc2_f656': np.empty((nsources,)), "wfpc2_f547": np.empty((nsources,))}
bg = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,)), 'wfpc2_f658': np.empty((nsources,)), 
'wfpc2_f656': np.empty((nsources,)), "wfpc2_f547": np.empty((nsources,))}

for isource, file_name in enumerate(file_list):
    with open(file_name) as f:
        data = json.load(f) 

    label_s.append(data["star"]["id"])
    label_['star'] = label_s

    for k in data.keys():
        for camera, label in label_dict.items():
            if k.startswith(label):
                imagename = k
                
                try:
                    shell_ = np.array(data[imagename]["shell"]["value"])
                    bg_ = np.array(data[imagename]["background"]["value"])
                except KeyError:
                    shell_ = 0.0
                    bg_ = 0.0

                shell[camera][isource]=float(shell_)
                bg[camera][isource]=float(bg_)

# Divide by the exposure times 
S_bg_f656 = (bg["wfpc2_f656"])/T_f656
S_bg_f658 = (bg["wfpc2_f658"])/T_f658
S_bg_f547 = (bg["wfpc2_f547"])/T_f547

# contribution of the continuum
Sf547_bg_f656 = (S_bg_f547*w_f656)/w_f547
Sf547_bg_f658 = (S_bg_f547*w_f658)/w_f547
Sf547_bg_f547 = (S_bg_f547*w_f547)/w_f547
Sf547_bg_acs_f658 = (S_bg_f547*w_acs_f658)/w_f547
Sf547_bg_wfpc2_f656 = (S_bg_f547*T_robb_f656)/w_f547

# Percentage
por_Sf547_bg_f656 = Sf547_bg_f656/S_bg_f656
por_Sf547_bg_f658 = Sf547_bg_f658/S_bg_f658
por_Sf547_bg_f547 = Sf547_bg_f547/S_bg_f547
por_Sf547_bg_acs_f658 = Sf547_bg_acs_f658/S_bg_f658
por_Sf547_bg_wfpc2_f656 = Sf547_bg_wfpc2_f656/S_bg_f656

# Table
def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table

# list for each column in  table
col_names = [ "Object", "S(Bg_f656)", "S(Bg_f658)", "S(Bg_f547)", "S(Bg_acs_f658)", "S(Bg_wfpc2_f656)", "S(cont_Bg_f656)", "S(cont_Bg_f658)", "S(cont_Bg_f547)", "S(cont_Bg_acs_f658)", "S(cont_Bg_wfpc2_f656)", "%(cont_Bg_f656)", "%(cont_Bg_f658)", "%(cont_Bg_f547)", "%(cont_Bg_acs_f658)", "%(cont_Bg_wfpc2_f656)"]

table = {cn: [] for cn in col_names}

def arcsec_fmt(r):
    """ Write distances to accuracy of 0.001 arcsec"""
    return "{:.3f}".format(r)

def bright_fmt(r):
    """ Write brightnesses to accuracy of 0.001"""
    return "{:.3f}".format(r)

for s, a, b, c, d, f, g, h, i, j, k, l, m, n, o, p in zip(label_s, S_bg_f656, S_bg_f658, S_bg_f547, S_bg_f658, S_bg_f656, Sf547_bg_f656, Sf547_bg_f658, Sf547_bg_f547, Sf547_bg_acs_f658, Sf547_bg_wfpc2_f656, por_Sf547_bg_f656, por_Sf547_bg_f658, por_Sf547_bg_f547, por_Sf547_bg_acs_f658, por_Sf547_bg_wfpc2_f656):
    table["Object"].append(s)
    table["S(Bg_f656)"].append(str(bright_fmt(a)))
    table["S(Bg_f658)"].append(str(bright_fmt(b)))
    table["S(Bg_f547)"].append(str(bright_fmt(c)))
    table["S(Bg_acs_f658)"].append(str(bright_fmt(d)))
    table["S(Bg_wfpc2_f656)"].append(str(bright_fmt(f)))
    table["S(cont_Bg_f656)"].append(str(bright_fmt(g)))
    table["S(cont_Bg_f658)"].append(str(bright_fmt(h)))
    table["S(cont_Bg_f547)"].append(str(bright_fmt(i)))
    table["S(cont_Bg_acs_f658)"].append(str(bright_fmt(j)))
    table["S(cont_Bg_wfpc2_f656)"].append(str(bright_fmt(k)))
    table["%(cont_Bg_f656)"].append(str(bright_fmt(l)))
    table["%(cont_Bg_f658)"].append(str(bright_fmt(m)))
    table["%(cont_Bg_f547)"].append(str(bright_fmt(n)))
    table["%(cont_Bg_acs_f658)"].append(str(bright_fmt(o)))
    table["%(cont_Bg_wfpc2_f656)"].append(str(bright_fmt(p)))

with open("arc-summary-continuum.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))
