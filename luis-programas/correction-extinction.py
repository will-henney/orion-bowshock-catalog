'''
Extinction correction
'''
import numpy as np
import json
from  astropy.table import Table, Column
import matplotlib.pyplot as plt


#Sfactor_WFPC2_three = 0.27849736582447426
#Sfactor_WFPC2 = 0.2734642808999356
#Sfactor_Odell = 0.1790544428870778 #multiplied by 2
Sfactor_ACS = 0.0025030687604156482

file_list = "extinction.json"

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t")
 
with open(file_list) as f:
    extinction_data = json.load(f)

tab.add_column(Column(name='Correction',
                      data=[10**(0.78*extinction_data.get(source, 0.0))
                            for source in tab['Object']]))


# Observed and corrected background brightness values
#m = (tab['D']!='-') and (tab['Value_bg_Bally']!='-')
D_arcmin = tab['D']/60.0
SHa_obs = Sfactor_ACS*tab['Value_bg_Bally']
SHa = SHa_obs*tab['Correction']
print D_arcmin
# Limits for plot
xmin, xmax = 0.06, 20.0
ymin, ymax = 5e-4, 2.0

# Comparison with O'Dell & Harris D^{-2} line
ODH_value = 2.7e-12 # Value @ 1 arcmin for S(Hb) in units of arcsec^{-2}
Ha_Hb = 2.87
radian_arcsec = np.radians(1.0/3600.0)
Dgrid = np.logspace(np.log10(xmin), np.log10(xmax))
ODH_S = ODH_value*Ha_Hb/radian_arcsec**2/Dgrid**2



#for obj, file_list in data.items():
    #for s, x, y in zip(tab['Object'], tab['D'], tab['Value_bg_Bally']):
        #if obj == s:
            #F = y*10**(0.78*file_list)  
            #d.append(x)
            #S.append(F)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r'$S(\mathrm{H\alpha})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$')
ax1.set_title(r'Surface brightness of background nebula (ACS F658N)')
ax1.plot(D_arcmin, SHa_obs, 'go', alpha=0.3, label='Observed')
ax1.plot(D_arcmin, SHa, 'bo', alpha=0.8, label='Extinction-corrected')
ax1.plot(Dgrid, ODH_S, 'k--', label=r"O'Dell & Harris: $D^{-2}$")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()
#ax1.legend(fontsize='small')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)
ax1.grid(True)

plt.savefig('S(alpha)_bg_correction_ACSluis.pdf')
plt.show()
               



               

