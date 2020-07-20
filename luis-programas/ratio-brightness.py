'''
ratio of NII and Halpha of the brightness of the mosaic and acs + WFPC2
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab_acs = Table.read("arc-summary-acs-wfpc.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))
tab_mos = Table.read("arc-summary-mosaic.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

# ratios brightness between acs and wfpc2
mask = np.isfinite(tab_acs['S(shell_NII/Halpha)']) & np.isfinite(tab_acs['sigma(shell_NII/Halpha)'])
mask = mask &  np.isfinite(tab_acs['S(Bg_NII/Halpha)']) &  np.isfinite(tab_acs['sigma(Bg_NII/Halpha)'])

mask_mos = np.isfinite(tab_mos['S(shell_NII/Halpha)']) & np.isfinite(tab_mos['sigma(shell_NII/Halpha)'])
mask_mos = mask_mos & np.isfinite(tab_mos['S(Bg_NII/Halpha)']) & np.isfinite(tab_mos['sigma(Bg_NII/Halpha)']) 

# Only plot points where errors are not too large
max_error = 1.16
mask = mask & (tab_acs['sigma(shell_NII/Halpha)'] <
               max_error*tab_acs['S(shell_NII/Halpha)'])
mask = mask & (tab_acs['sigma(Bg_NII/Halpha)'] <
                     max_error*tab_acs['S(Bg_NII/Halpha)'])

max_error = 1.87
mask_mos = mask_mos & (tab_mos['sigma(shell_NII/Halpha)'] <
               max_error*tab_mos['S(shell_NII/Halpha)'])
mask_mos = mask_mos & (tab_mos['sigma(Bg_NII/Halpha)'] <
                     max_error*tab_mos['S(Bg_NII/Halpha)'])


D_arcmin = tab_acs['D'][mask]/60.0
niiha_acswfpc2 = tab_acs['S(shell_NII/Halpha)'][mask]
sigma_niiha_acswfpc2 = tab_acs['sigma(shell_NII/Halpha)'][mask]
Bg_niiha_acswfpc2 = tab_acs['S(Bg_NII/Halpha)'][mask]
Bg_sigma_niiha_acswfpc2 = tab_acs['sigma(Bg_NII/Halpha)'][mask]

# ratios brightness between mosaics
D_arcmin1 = tab_mos['D'][mask_mos]/60.0
niiha_mosaic = tab_mos['S(shell_NII/Halpha)'][mask_mos]
sigma_niiha_mosaic = tab_mos['sigma(shell_NII/Halpha)'][mask_mos]    
Bg_niiha_mosaic = tab_mos['S(Bg_NII/Halpha)'][mask_mos]
Bg_sigma_niiha_mosaic = tab_mos['sigma(Bg_NII/Halpha)'][mask_mos]

#fit line
m = (niiha_acswfpc2>0.0)&(sigma_niiha_acswfpc2>0.0)
Distance_log = np.log10(D_arcmin[m])
a = np.polyfit(Distance_log, np.log10(niiha_acswfpc2[m]), 1, w=1./sigma_niiha_acswfpc2[m])
p = np.poly1d(a)
x_grid = np.linspace(0.07, 11.0, 200)
yfit = 10**p(np.log10(x_grid))
print('ACS-WFPC2 fit: Ratio = {:.2f} D**{:.2f}'.format(10**a[1], a[0]))

#fit line mosaic
m1 = (niiha_mosaic>0.0) & (niiha_mosaic < 1.0) &(sigma_niiha_mosaic>0.0)
Distance1_log = np.log10(D_arcmin1[m1])
a = np.polyfit(Distance1_log, np.log10(niiha_mosaic[m1]), 1, w=1./sigma_niiha_mosaic[m1])
p = np.poly1d(a)
yfit1 = 10**p(np.log10(x_grid))
print('Pure WFPC2 fit: Ratio = {:.2f} D**{:.2f}'.format(10**a[1], a[0]))

# fit combined dataset
x = np.concatenate((Distance_log, Distance1_log))
y = np.concatenate((niiha_acswfpc2[m], niiha_mosaic[m1]))
w = np.concatenate((3./sigma_niiha_acswfpc2[m], 1./sigma_niiha_mosaic[m1]))
a = np.polyfit(x, np.log10(y), 1)
p = np.poly1d(a)
yfit2 = 10**p(np.log10(x_grid))
print('Combined fit: Ratio = {:.2f} D**{:.2f}'.format(10**a[1], a[0]))

fig = plt.figure()
ax3 = fig.add_subplot(1,1,1)
ax3.set_xlim(xmin=0.07, xmax=11.0)
ax3.set_ylim(ymin=-0.3, ymax=1.8)
#ax3.set_xlim(xmin=0.06, xmax=11.0)
#ax3.set_ylim(ymin=0.004, ymax=250.0)
#ax3.set_ylim(ymin=0.0001, ymax=200.0)
ax3.set_xlabel(r'Projected distance from Trapezium, D / arcmin')
ax3.set_ylabel(r'Brightness ratio ACS-WFPC2, S = S(NII)/S(H alpha)')
ax3.plot(D_arcmin, niiha_acswfpc2, 'ko',  alpha = 0.3, label = 'Shell-brightness')
ax3.errorbar(D_arcmin, niiha_acswfpc2,  yerr = sigma_niiha_acswfpc2, marker='o', alpha = 0.3, fmt='ko')
ax3.plot(D_arcmin, Bg_niiha_acswfpc2, 'ro', alpha = 0.4, label = 'Background-brightness')
ax3.errorbar(D_arcmin, Bg_niiha_acswfpc2,  yerr = Bg_sigma_niiha_acswfpc2, marker='o', alpha = 0.4, fmt='ro')
ax3.plot(x_grid, yfit, 'k--')
ax3.plot(x_grid, yfit2, 'k-.')
# ax3.set_xscale('log')
#ax3.set_yscale('log')
ax3.legend(fontsize='small', loc=2)
#ax3.grid(True)
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, color='k', hold=None, zorder=-5.0)
plt.savefig('acs_wfpc2-ratio-Nii_ha-vs-D_-mean-error_new.pdf')

fig = plt.figure()
ax4 = fig.add_subplot(1,1,1)
ax4.set_xlim(xmin=0.07, xmax=11.0)
ax4.set_ylim(ymin=-0.3, ymax=1.8)
ax4.set_xlabel(r'Projected distance from Trapezium, D / arcmin')
ax4.set_ylabel(r'Brightness ratio mosaic-WFPC2 , S = S(NII)/S(H alpha)')
ax4.plot(D_arcmin1, niiha_mosaic, 'ko', alpha=0.3, label = 'Shell-brigthness')
ax4.errorbar(D_arcmin1, niiha_mosaic, yerr = sigma_niiha_mosaic, marker='o', alpha =0.3, fmt='ko')
ax4.plot(D_arcmin1, Bg_niiha_mosaic, 'ro', alpha=0.4, label = 'background-brightness')
ax4.errorbar(D_arcmin1, Bg_niiha_mosaic,  yerr = Bg_sigma_niiha_mosaic, marker='o', alpha=0.4, fmt='ro')
ax4.plot(x_grid, yfit1, 'k--')
ax4.plot(x_grid, yfit2, 'k-.')
# ax4.set_xscale('log')
#ax4.set_yscale('log')
ax4.legend(fontsize='small', loc=2)
#ax4.grid(True)
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, color='k', hold=None, zorder=-5.0)
plt.savefig('wfpc2-mosaic-ratio-Nii_ha-vs-D_-mean-error_new.pdf')
