'''
Estimate correlation of brightness and error
'''
import glob
import json
import numpy as np
import matplotlib.pyplot as plt

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

Sfactor_ACS = 0.0025030687604156482
Sfactor_WFPC2 = 0.16895605909108011

shell_mean = {}
sigma_mean = {}
bg_mean = {}
dbg_mean = {}
distance = {}
label_dict = {"acs":"Bally", "wfpc2":"Robberto_WFPC2"}
label_s = []

nsources = len(file_list)

# allocate empty arrays to store the mean and sigma
shell_mean = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,))}
sigma_mean = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,))}
bg_mean = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,))}
dbg_mean = {'acs': np.empty((nsources,)), 'wfpc2': np.empty((nsources,))}
distance = {'star': np.empty((nsources,))}

for isource, file_name in enumerate(file_list):
    with open(file_name) as f:
        data = json.load(f)
    
    #Distance to ionizing star
    distance['star'][isource] = np.array(data["star"]["D"])    
    label_s.append(data["star"]["id"])
    #print label_s
    for k in data.keys():
        for camera, label in label_dict.items():
            if k.startswith(label):
                imagename = k
        
                #print(file_name, "Found image ", imagename
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
                m = np.abs(theta) < 45.0
                n = m.sum() 
               
                #error
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
                
                if ((np.isfinite(shell_dif_mean) and shell_dif_mean > 0.0)
                       and (np.isfinite(sigma_shell_mean) and sigma_shell_mean > 0.0)): 
                    shell_mean[camera][isource] = float(shell_dif_mean)
                    sigma_mean[camera][isource] = float(sigma_shell_mean)
                else:
                    shell_mean[camera][isource] = 0.0
                    sigma_mean[camera][isource] = 0.0
                
                
                bg_mean[camera][isource] = float(Bg)
                dbg_mean[camera][isource] = float(sigma_bg_mean)

#ratio NII and H alpha for shell
Mask  = (shell_mean['acs']!=0.0)&(shell_mean['wfpc2']!=0.0)
D_arcmin = distance['star'][Mask]/60.0
ratio_haNii =  ((Sfactor_ACS*shell_mean['acs'][Mask])/(Sfactor_WFPC2*shell_mean['wfpc2'][Mask])) - 1      

#bg
ratio_haNii_bg = ((Sfactor_ACS*bg_mean['acs'][Mask])/(Sfactor_WFPC2*bg_mean['wfpc2'][Mask])) - 1
           
#propagation the uncertainty for ratio                 
sigma_ratio_haNii = np.abs((Sfactor_ACS*shell_mean['acs'][Mask])/(Sfactor_WFPC2*shell_mean['wfpc2'][Mask]))*(
np.sqrt((sigma_mean['acs'][Mask]/shell_mean['acs'][Mask])**2 + (sigma_mean['wfpc2'][Mask]/shell_mean['wfpc2'][Mask])**2))

# sigma bg
sigma_ratio_haNiibg = np.abs((Sfactor_ACS*bg_mean['acs'][Mask])/(Sfactor_WFPC2*bg_mean['wfpc2'][Mask]))*(
np.sqrt((dbg_mean['acs'][Mask]/bg_mean['acs'][Mask])**2 + (dbg_mean['wfpc2'][Mask]/bg_mean['wfpc2'][Mask])**2))

# table
def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table

# list for each column in  table
col_names = [ "Object", "S(shell_mean_acs)", "S(shell_mean_wfpc2)", "sigma(acs)", 
             "sigma(wfpc2)"] #, "S(NII/Halpha),shell", "sigma(NII/Halpha), shell"]

table = {cn: [] for cn in col_names}

for s, a, b, c , d in zip(label_s, Sfactor_ACS*shell_mean['acs'], Sfactor_WFPC2*shell_mean['wfpc2'], Sfactor_ACS*sigma_mean['acs'], Sfactor_WFPC2*sigma_mean['wfpc2']): #, ratio_haNii, sigma_ratio_haNii)
    table["Object"].append(s)
    table["S(shell_mean_acs)"].append(str(a))
    table["S(shell_mean_wfpc2)"].append(str(b))
    table["sigma(acs)"].append(str(c))
    table["sigma(wfpc2)"].append(str(d))
    #table["S(NII/Halpha),shell"].append(str(e))
    #table["sigma(NII/Halpha), shell"].append(str(f))    

with open("arc-uncertainty.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))
                
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(xmin=0.2e-4, xmax=1.0)
ax1.set_ylim(ymin=0.1e-5, ymax=1e-1)
ax1.set_xlabel(r'H alpha and NII surface brightness, $S(\mathrm{H\alpha+NII})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$')
ax1.set_ylabel(r'H alpha surface brightness, $S(\mathrm{H\alpha})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$')
#ax1.plot(Sfactor_ACS*shell_mean['acs'], Sfactor_WFPC2*shell_mean['wfpc2'], 'r.')
ax1.errorbar(Sfactor_ACS*shell_mean['acs'], Sfactor_WFPC2*shell_mean['wfpc2'], xerr = Sfactor_ACS*sigma_mean['acs'], yerr = Sfactor_WFPC2*sigma_mean['wfpc2'], marker='o', fmt='ro')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True)
plt.savefig('acs_wfpc2-correlation-brightness_error.pdf')

fig = plt.figure()
ax2 = fig.add_subplot(1,1,1)
ax2.set_xlim(xmin=0.06, xmax=15.0)
ax2.set_ylim(ymin=0.004, ymax=250.0)
ax2.set_xlabel(r'D, arcmin')
ax2.set_ylabel(r'Brightness fraction NII and Halpha, S = S(NII)/S(Halpha)')
ax2.plot(D_arcmin, ratio_haNii, 'ro', label = 'Shell-brightness')
ax2.errorbar(D_arcmin, ratio_haNii,  yerr = sigma_ratio_haNii, marker='o', fmt='ro')
ax2.plot(D_arcmin, ratio_haNii_bg, 'bo', label = 'Background-brightness')
ax2.errorbar(D_arcmin, ratio_haNii_bg,  yerr = sigma_ratio_haNiibg, marker='o', fmt='bo')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.legend()
plt.savefig('acs_wfpc2-ratio-Nii_ha-vs-D_-mean-error.pdf')


