import numpy as np
import json
from  astropy.table import Table, Column
import matplotlib.pyplot as plt

# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0012357713647687021

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t")
 
with open("extinction.json") as f:
    extinction_data = json.load(f)
    

### Verbose version
# red_corr_ha = []
# for source in tab['Object']:
#     chb = extinction_data.get(source, 0.0)
#     red_corr_ha.append(10**(0.78*chb))
# new_column = Column(name='Correction', data=red_corr_ha)
# tab.add_column(new_column)

### Compact version
tab.add_column(Column(name='Correction',
                      data=[10**(0.78*extinction_data.get(source, 0.0))
                            for source in tab['Object']]))

# Observed and corrected background brightness values
D_arcmin = tab['D']/60.0
SHa_obs = Sfactor_ACS*tab['Value_bg_Bally']
SHa = SHa_obs*tab['Correction']

# Limits for plot
xmin, xmax = 0.06, 20.0
ymin, ymax = 5e-4, 2.0

# Comparison with O'Dell & Harris D^{-2} line
ODH_value = 2.7e-12 # Value @ 1 arcmin for S(Hb) in units of arcsec^{-2}
Ha_Hb = 2.87
radian_arcsec = np.radians(1.0/3600.0)
Dgrid = np.logspace(np.log10(xmin), np.log10(xmax))
ODH_S = ODH_value*Ha_Hb/radian_arcsec**2/Dgrid**2

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r'$S(\mathrm{H\alpha})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$')
ax1.set_title(r'Surface brightness of background nebula (ACS F658N)')
ax1.plot(D_arcmin, SHa_obs, 'ro', alpha=0.2, label='Observed')
ax1.plot(D_arcmin, SHa, 'bo', alpha=0.8, label='Extinction-corrected')
ax1.plot(Dgrid, ODH_S, 'k--', label=r"O'Dell & Harris: $D^{-2}$")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend(fontsize='small')
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(ymin, ymax)
ax1.grid(True)

plt.savefig('S(alpha)_bg_correction_extlog_will.pdf')
plt.show()
               

