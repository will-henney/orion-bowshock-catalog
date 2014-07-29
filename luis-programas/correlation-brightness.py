'''
Estimate correlation of brightness and error
'''
import glob
import json
import numpy as np
import matplotlib.pyplot as plt

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

no_wfcp2 = ["117-421", "042-628", "074-421"]

Sfactor_ACS = 0.0025030687604156482
Sfactor_WFPC2 = 0.16895605909108011

shell_mean = {}
sigma_mean = {}
bg_mean = {}
dbg_mean = {}
distance = {}
label_dict = {"acs":"Bally", "wfpc2":"Robberto_WFPC2"}
label_s = []
label_ = {}
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
    label_['star'] = label_s
    
    for k in data.keys():
        for camera, label in label_dict.items():
            if k.startswith(label):
                imagename = k
                
                try:
                    theta = np.array(data[imagename]["binned"]["theta"])  
                    shell = np.array(data[imagename]["binned"]["shell"])
                    bg = np.array(data[imagename]["binned"]["background"])
                    dbg = np.array(data[imagename]["binned"]["background sigma"])
                except KeyError:
                    None
  
                #mask for filter the values and number total of points "n"
                m = (np.abs(theta) < 45.0)&(np.isfinite(shell))
                m = m & (np.isfinite(bg)&np.isfinite(dbg)) 
                n = m.sum() 
               
                #error
                try:
                    shell_dif_mean = np.mean(shell[m] - bg[m])
                    Bg = np.mean(bg[m]) 
                    dbg_dif = np.sqrt(2)*dbg[m]
                    sigma_shell_mean = np.sqrt(np.sum(dbg_dif**2))/n
                    sigma_bg_mean =  np.sqrt(np.sum(dbg[m]**2))/n 
                except TypeError:
                    None
                print Bg
                shell_mean[camera][isource] = float(shell_dif_mean)
                sigma_mean[camera][isource] = float(sigma_shell_mean)
                bg_mean[camera][isource] = float(Bg)
                dbg_mean[camera][isource] = float(sigma_bg_mean)
                
             
#ratio NII and H alpha for shell
#Mask  = (shell_mean['acs']!=0.0)&(shell_mean['wfpc2']!=0.0)
D_arcmin = distance['star']/60.0
ratio_haNii =  ((Sfactor_ACS*shell_mean['acs'])/(Sfactor_WFPC2*shell_mean['wfpc2'])) - 1      

#bg
ratio_haNii_bg = ((Sfactor_ACS*bg_mean['acs'])/(Sfactor_WFPC2*bg_mean['wfpc2'])) - 1
           
#propagation the uncertainty for ratio                 
sigma_ratio_haNii = np.abs((Sfactor_ACS*shell_mean['acs'])/(Sfactor_WFPC2*shell_mean['wfpc2']))*(
np.sqrt((sigma_mean['acs']/shell_mean['acs'])**2 + (sigma_mean['wfpc2']/shell_mean['wfpc2'])**2))

# sigma bg
sigma_ratio_haNiibg = np.abs((Sfactor_ACS*bg_mean['acs'])/(Sfactor_WFPC2*bg_mean['wfpc2']))*(
np.sqrt((dbg_mean['acs']/bg_mean['acs'])**2 + (dbg_mean['wfpc2']/bg_mean['wfpc2'])**2))

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
col_names = [ "Object", "D", "S(shell_mean_acs)", "sigma(acs)", "S(shell_mean_wfpc2)",
             "sigma(wfpc2)","S(Bg_acs)","sigma(Bg_acs)", "S(Bg_wfpc2)",  "sigma(Bg_wfpc2)", 
             "S(shell_NII/Halpha)", "sigma(shell_NII/Halpha)", "S(Bg_NII/Halpha)", "sigma(Bg_NII/Halpha)"]

table = {cn: [] for cn in col_names}

def arcsec_fmt(r):
    """ Write distances to accuracy of 0.001 arcsec"""
    return "{:.3f}".format(r)

def bright_fmt(r):
    """ Write brightnesses to accuracy of 0.001"""
    return "{:.5f}".format(r)

for s, a, b, c, d, e, f, g, h, i, j, k, l, n in zip(label_['star'], distance['star'], 
      Sfactor_ACS*shell_mean['acs'], Sfactor_ACS*sigma_mean['acs'], Sfactor_WFPC2*shell_mean['wfpc2'],
       Sfactor_WFPC2*sigma_mean['wfpc2'], Sfactor_ACS*bg_mean['acs'], Sfactor_ACS*dbg_mean['acs'], 
       Sfactor_WFPC2*bg_mean['wfpc2'], Sfactor_WFPC2*dbg_mean['wfpc2'], ratio_haNii, sigma_ratio_haNii, 
       ratio_haNii_bg, sigma_ratio_haNiibg ):
    table["Object"].append(s)
    table["D"].append(str(arcsec_fmt(a)))
    table["S(shell_mean_acs)"].append(str(bright_fmt(b)))
    table["sigma(acs)"].append(str(bright_fmt(c)))
    table["S(shell_mean_wfpc2)"].append(str(bright_fmt(d)))
    table["sigma(wfpc2)"].append(str(bright_fmt(e)))
    table["S(Bg_acs)"].append(str(bright_fmt(f)))
    table["sigma(Bg_acs)"].append(str(bright_fmt(g)))
    table["S(Bg_wfpc2)"].append(str(bright_fmt(h)))
    table["sigma(Bg_wfpc2)"].append(str(bright_fmt(i)))
    table["S(shell_NII/Halpha)"].append(str(bright_fmt(j)))
    table["sigma(shell_NII/Halpha)"].append(str(bright_fmt(k)))   
    table["S(Bg_NII/Halpha)"].append(str(bright_fmt(l)))
    table["sigma(Bg_NII/Halpha)"].append(str(bright_fmt(n))) 

with open("arc-summary-acs-wfpc.tab", "w") as f:
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


