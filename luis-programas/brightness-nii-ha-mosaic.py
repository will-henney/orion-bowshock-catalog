'''
Estimate correlation of brightness and error
'''
import glob
import json
import numpy as np
import matplotlib.pyplot as plt

# photometric keywords from the mosaic-WFPC2 header
Sfactor_f656 = 9.341771211049384e-05
Sfactor_f658 = 3.771922832375e-05

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

shell_mean = {}
sigma_mean = {}
bg_mean = {}
dbg_mean = {}
distance = {}
label_dict = {"WFC_f658":"WFC_mosaic_f658", "WFC_f656":"WFC_mosaic_f656"}
label_s = []
label_ = {}

nsources = len(file_list)

# allocate empty arrays to store the mean and sigma
shell_mean = {'WFC_f658': np.empty((nsources,)), 'WFC_f656': np.empty((nsources,))}
sigma_mean = {'WFC_f658': np.empty((nsources,)), 'WFC_f656': np.empty((nsources,))}
bg_mean = {'WFC_f658': np.empty((nsources,)), 'WFC_f656': np.empty((nsources,))}
dbg_mean = {'WFC_f658': np.empty((nsources,)), 'WFC_f656': np.empty((nsources,))}
distance = {'star': np.empty((nsources,))}

for isource, file_name in enumerate(file_list):
    with open(file_name) as f:
        data = json.load(f)
    
    #Distance to ionizing star
    distance['star'][isource] = np.array(data["star"]["D"])    
    label_s.append(data["star"]["id"])
    label_['star'] = label_s
   
    for camera, label in label_dict.items():
        imagename = label
         
        try:
            theta = np.array(data[imagename]["binned"]["theta"])  
            shell = np.array(data[imagename]["binned"]["shell"])
            bg = np.array(data[imagename]["binned"]["background"])
            dbg = np.array(data[imagename]["binned"]["background sigma"])
        except KeyError:
            theta = 0.0
            shell = 0.0  
            bg = 0.0
            dbg = 0.0
        
        #mask for filter the values and number total of points "n"
        m = (np.abs(theta) < 35.0)&(np.isfinite(shell))
        m = m & (np.isfinite(bg)&np.isfinite(dbg))
        n = m.sum()
       
        #estimate value of the shell, bg and its mean
        try:
            shell_dif_mean = np.mean(shell[m] - bg[m])
            Bg = np.mean(bg[m]) 
            dbg_dif = np.sqrt(2)*dbg[m]
            sigma_shell_mean = np.sqrt(np.sum(dbg_dif**2))/n
            sigma_bg_mean =  np.sqrt(np.sum(dbg[m]**2))/n 
        except TypeError:
            shell_dif_mean = 0.0
            Bg = 0.0
            dbg_dif = 0.0
            sigma_shell_mean = 0.0
            sigma_bg_mean =  0.0
        
        shell_mean[camera][isource] = float(shell_dif_mean)
        sigma_mean[camera][isource] = float(sigma_shell_mean)
        bg_mean[camera][isource] = float(Bg)
        dbg_mean[camera][isource] = float(sigma_bg_mean)
                       
#ratio NII and H alpha for shell and Bg
#Mask  = (shell_mean['WFC_f658']!=0.0)&(shell_mean['WFC_f656']!=0.0)
D_arcmin = distance['star']/60.0
Nii_ha = (Sfactor_f658*shell_mean['WFC_f658'])/(Sfactor_f656*shell_mean['WFC_f656'])
Bg_Nii_ha = (Sfactor_f658*bg_mean['WFC_f658'])/(Sfactor_f656*bg_mean['WFC_f656'])

#uncertainty of shell and Bg
sigma_ratio_Niiha = np.abs((Sfactor_f658*shell_mean['WFC_f658'])/(Sfactor_f656*shell_mean['WFC_f656']))*(np.sqrt((sigma_mean['WFC_f658']/shell_mean['WFC_f658'])**2 + (sigma_mean['WFC_f656']/shell_mean['WFC_f656'])**2))
Bg_sigma_ratio_Niiha = np.abs((Sfactor_f658*bg_mean['WFC_f658'])/(Sfactor_f656*bg_mean['WFC_f656']))*(
np.sqrt((dbg_mean['WFC_f658']/bg_mean['WFC_f658'])**2 + (dbg_mean['WFC_f656']/bg_mean['WFC_f656'])**2))

#Tabla
def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table

# list for each column in  table
col_names = [ "Object",  "D", "S(shell_mosaicf658)", "sigma(shell_mosaicf658)", "S(shell_mosaicf656)", 
             "sigma(shell_mosaicf656)", "S(Bg_mosaicf658)","sigma(Bg_mosaicf658)", "S(Bg_mosaicf656)",  
             "sigma(Bg_mosaicf656)", "S(shell_NII/Halpha)", "sigma(shell_NII/Halpha)", "S(Bg_NII/Halpha)", "sigma(Bg_NII/Halpha)"] 

table = {cn: [] for cn in col_names}

def arcsec_fmt(r):
    """ Write distances to accuracy of 0.001 arcsec"""
    return "{:.3f}".format(r)

def bright_fmt(r):
    """ Write brightnesses to accuracy of 0.001"""
    return "{:.5f}".format(r)

    
for s, a, b, c , d, e , f, g, h, i, j, k, l, n in zip(label_['star'], distance['star'], Sfactor_f658*shell_mean['WFC_f658'], Sfactor_f658*sigma_mean['WFC_f658'], Sfactor_f656*shell_mean['WFC_f656'], Sfactor_f656*sigma_mean['WFC_f656'], Sfactor_f658*bg_mean['WFC_f658'], Sfactor_f658*dbg_mean['WFC_f658'], Sfactor_f656*bg_mean['WFC_f656'], Sfactor_f656*dbg_mean['WFC_f656'], Nii_ha, sigma_ratio_Niiha, Bg_Nii_ha, Bg_sigma_ratio_Niiha): 
    table["Object"].append(s)
    table["D"].append(str(arcsec_fmt(a)))
    table["S(shell_mosaicf658)"].append(str(bright_fmt(b)))
    table["sigma(shell_mosaicf658)"].append(str(bright_fmt(c)))
    table["S(shell_mosaicf656)"].append(str(bright_fmt(d)))
    table["sigma(shell_mosaicf656)"].append(str(bright_fmt(e)))
    table["S(Bg_mosaicf658)"].append(str(bright_fmt(f)))  
    table["sigma(Bg_mosaicf658)"].append(str(bright_fmt(g)))
    table["S(Bg_mosaicf656)"].append(str(bright_fmt(h)))
    table["sigma(Bg_mosaicf656)"].append(str(bright_fmt(i)))  
    table["S(shell_NII/Halpha)"].append(str(bright_fmt(j)))  
    table["sigma(shell_NII/Halpha)"].append(str(bright_fmt(k)))
    table["S(Bg_NII/Halpha)"].append(str(bright_fmt(l)))
    table["sigma(Bg_NII/Halpha)"].append(str(bright_fmt(n))) 
    
with open("arc-summary-mosaic.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(xmin=0.06, xmax=10.0)
ax1.set_ylim(ymin=0.001, ymax=20.0)
ax1.set_xlabel(r'D, arcmin')
ax1.set_ylabel(r'Brightness fraction NII and Halpha, S = S(NII)/S(Halpha)')
ax1.plot(D_arcmin, Nii_ha, 'ro', label = 'Shell-brightness')
ax1.errorbar(D_arcmin, Nii_ha,  yerr = sigma_ratio_Niiha, marker='o', fmt='ro')
ax1.plot(D_arcmin, Bg_Nii_ha, 'ko', alpha=0.4, label = 'Background-brightness')
ax1.errorbar(D_arcmin, Bg_Nii_ha,  yerr = Bg_sigma_ratio_Niiha, alpha=0.4, marker='o', fmt='ko')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True)
ax1.legend()
plt.savefig('mosaic_ratio-Nii_ha-vs-D_error.pdf')
















