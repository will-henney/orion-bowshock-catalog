'''
brightness of shell 
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table
import json

fha = 0.78

tab_acs = Table.read("arc-summary-acs-wfpc.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

mask = np.isfinite(tab_acs['S(shell_mean_acs)']) & np.isfinite(tab_acs['sigma(acs)'])

max_error = 1.16
mask = mask & (tab_acs['sigma(acs)'] <
               max_error*tab_acs['S(shell_mean_acs)'])

with open("extinction.json") as f:
    extinction_data = json.load(f)

Chb = np.array([extinction_data.get(source, 0.0) for source in tab_acs['Object'][mask]])
D_arcmin = tab_acs['D'][mask]/60.0
niiha_acs = tab_acs['S(shell_mean_acs)'][mask]
sigma_niiha_acs = tab_acs['sigma(acs)'][mask]
niiha_acs *= 10**(fha*Chb) 
sigma_niiha_acs *= 10**(fha*Chb)

# Correct for [N II] contamination of Ha filter
# This comes from the fit done in luis-programas/ratio-brightness.py
# Combined fit: Ratio = 0.28 D**0.43
Rnii_ha = 0.28*D_arcmin**0.43
niiha_acs /= 1.0 + Rnii_ha
sigma_niiha_acs /= 1.0 + Rnii_ha

#brightness WFPC2
ha_wfpc2 = tab_acs['S(shell_mean_wfpc2)']
sigma_ha_wfpc2 = tab_acs['sigma(wfpc2)']

#Ajuste
m = (niiha_acs>0.0)
Distance_log = np.log10(D_arcmin[m])
niiha_acs_log = np.log10(niiha_acs[m])
a = np.polyfit(Distance_log, niiha_acs_log, 1, w = 1./np.sqrt(sigma_niiha_acs))
p = np.poly1d(a)
x_grid_log = np.linspace(Distance_log.min(), Distance_log.max())
print('ACS brightness fit:  = {:.2f} D**{:.2f}'.format(10**a[1], 10**a[0]))

fig = plt.figure()
ax = fig.add_subplot(1,1,1, axisbg="#eeeeee")
#ax.set_ylim(ymin=-0.01, ymax=0.1)
ax.set_xlim(xmin=0.06, xmax=12.0)
#ax.set_ylim(ymin=-0.3, ymax=1.8)
ax.set_xlabel(r'Projected distance from Trapezium, D / arcmin')
ax.set_ylabel(r'H alpha surface brightness, erg/s/cm2/sr')
ax.plot(D_arcmin, niiha_acs, 'gs',  alpha = 0.4)
ax.errorbar(D_arcmin, niiha_acs,  yerr = sigma_niiha_acs, marker='s', alpha = 0.5, fmt='gs', capthick=2.0, linewidth=2.0)
ax.plot(10**(x_grid_log), 10**(p(x_grid_log)), 'k-')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.legend(fontsize='small', loc=2)
#ax.grid(True)
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, color='k', hold=None, zorder=-5.0)
plt.savefig('brightness-shell_acs_Vs-D_new.pdf')
   
