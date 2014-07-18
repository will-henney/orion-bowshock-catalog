'''
Estimate correlation of brightness and error
'''
import glob
import json
import numpy as np
import matplotlib.pyplot as plt

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

shell_mean = {}
sigma_mean = {}
bg_mean = {}
dbg_mean = {}
distance = {}
label_dict = {"WFC_f658":"WFC_mosaic_f658", "WFC_f656":"WFC_mosaic_f656"}


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
        m = (np.abs(theta) < 35.0)
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
Mask  = (shell_mean['WFC_f658']!=0.0)&(shell_mean['WFC_f656']!=0.0)
D_arcmin = distance['star'][Mask]/60.0
Nii_ha = (shell_mean['WFC_f658'][Mask])/(shell_mean['WFC_f656'][Mask])
Bg_Nii_ha = (bg_mean['WFC_f658'][Mask])/(bg_mean['WFC_f656'][Mask])

#uncertainty of shell and Bg
sigma_ratio_Niiha = np.abs((shell_mean['WFC_f658'][Mask])/(shell_mean['WFC_f656'][Mask]))*(np.sqrt((sigma_mean['WFC_f658'][Mask]/shell_mean['WFC_f658'][Mask])**2 + (sigma_mean['WFC_f656'][Mask]/shell_mean['WFC_f656'][Mask])**2))
Bg_sigma_ratio_Niiha = np.abs((bg_mean['WFC_f658'][Mask])/(bg_mean['WFC_f656'][Mask]))*(
np.sqrt((dbg_mean['WFC_f658'][Mask]/bg_mean['WFC_f658'][Mask])**2 + (dbg_mean['WFC_f656'][Mask]/bg_mean['WFC_f656'][Mask])**2))

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(xmin=0.06, xmax=10.0)
ax1.set_ylim(ymin=0.05, ymax=20.0)
ax1.set_xlabel(r'D, arcmin')
ax1.set_ylabel(r'Vlaue Brightness fraction NII and Halpha, Value = Value(NII)/Value(Halpha)')
ax1.plot(D_arcmin, Nii_ha, 'ro', label = 'Shell-brightness')
ax1.errorbar(D_arcmin, Nii_ha,  yerr = sigma_ratio_Niiha, marker='o', fmt='ro')
ax1.plot(D_arcmin, Bg_Nii_ha, 'ko', alpha=0.4, label = 'Background-brightness')
ax1.errorbar(D_arcmin, Bg_Nii_ha,  yerr = Bg_sigma_ratio_Niiha, alpha=0.4, marker='o', fmt='ko')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True)
ax1.legend()
plt.savefig('mosaic_ratio-Nii_ha-vs-D_error.pdf')
















