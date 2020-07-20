'''
Compare the ratio of the brightness of the mosaic and acs + WFPC2
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab_acs = Table.read("arc-summary-acs-wfpc.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))
tab_mos = Table.read("arc-summary-mosaic.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

# ratios brightness between acs and wfpc2
mask = np.isfinite(tab_acs['S(shell_NII/Halpha)']) & np.isfinite(tab_mos['S(shell_NII/Halpha)'])
mask = mask & np.isfinite(tab_acs['sigma(shell_NII/Halpha)']) & np.isfinite(tab_mos['sigma(shell_NII/Halpha)'])
mask_Bg = np.isfinite(tab_acs['S(Bg_NII/Halpha)']) & np.isfinite(tab_mos['S(Bg_NII/Halpha)']) 
mask_Bg = mask_Bg & np.isfinite(tab_acs['sigma(Bg_NII/Halpha)']) & np.isfinite(tab_mos['sigma(Bg_NII/Halpha)'])

# Only plot points where errors are not too large
max_error = 1.3
mask = mask & (tab_acs['sigma(shell_NII/Halpha)'] <
               max_error*tab_acs['S(shell_NII/Halpha)'])
mask = mask & (tab_mos['sigma(shell_NII/Halpha)'] <
               max_error*tab_mos['S(shell_NII/Halpha)'])

max_error = 0.5
mask_Bg = mask_Bg & (tab_acs['sigma(Bg_NII/Halpha)'] <
                     max_error*tab_acs['S(Bg_NII/Halpha)'])
mask_Bg = mask_Bg & (tab_mos['sigma(Bg_NII/Halpha)'] <
                     max_error*tab_mos['S(Bg_NII/Halpha)'])

D_arcmin = tab_acs['D'][mask]/60.0
niiha_acswfpc2 = tab_acs['S(shell_NII/Halpha)'][mask]
sigma_niiha_acswfpc2 = tab_acs['sigma(shell_NII/Halpha)'][mask]
Bg_niiha_acswfpc2 = tab_acs['S(Bg_NII/Halpha)'][mask_Bg]
Bg_sigma_niiha_acswfpc2 = tab_acs['sigma(Bg_NII/Halpha)'][mask_Bg]

# ratios brightness between mosaics
D_arcmin1 = tab_mos['D'][mask]/60.0
niiha_mosaic = tab_mos['S(shell_NII/Halpha)'][mask]
sigma_niiha_mosaic = tab_mos['sigma(shell_NII/Halpha)'][mask]    
Bg_niiha_mosaic = tab_mos['S(Bg_NII/Halpha)'][mask_Bg]
Bg_sigma_niiha_mosaic = tab_mos['sigma(Bg_NII/Halpha)'][mask_Bg]

#Fitting a straight line shell
#m = (niiha_acswfpc2>0.0)&(niiha_mosaic>0.0)
#x_log = np.log10(niiha_acswfpc2[m])
#y_log = np.log10(niiha_mosaic[m])
#a = np.polyfit(x_log, x_log, 1)
#a1 = np.polyfit(niiha_acswfpc2, niiha_mosaic, 1)
#x_grid_log = np.linspace(x_log.min(), x_log.max())
a = np.polyfit(niiha_acswfpc2, niiha_acswfpc2, 1)
p = np.poly1d(a)
x_grid_log = np.linspace(niiha_acswfpc2.min(), niiha_mosaic.max())

#Fitting a straight line background
#mm = (Bg_niiha_acswfpc2>0.0)&(Bg_niiha_mosaic>0.0)
#xx_log = np.log10(Bg_niiha_acswfpc2[m1])
#yy_log = np.log10(Bg_niiha_mosaic[m1])
aa = np.polyfit(Bg_niiha_acswfpc2, Bg_niiha_acswfpc2, 1)
pp = np.poly1d(aa)
xx_grid_log = np.linspace(Bg_niiha_acswfpc2.min(), Bg_niiha_acswfpc2.max())

fig = plt.figure(figsize=(8, 6.8))
ax1 = fig.add_subplot(1,1,1)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(xmin=-0.03, xmax=1.5)
ax1.set_ylim(ymin=-0.03, ymax=1.5)
ax1.set_xlabel(r'Brightness acs and wfpc2, S = S(NII)/S(Halpha)')
ax1.set_ylabel(r'Brightness mosaic, S = S(NII)/S(Halpha)')
#ax1.plot(D_arcmin, niiha_acswfpc2, 'ro', label = 'ACS-WFPC2-Shell')
ax1.errorbar(niiha_acswfpc2, niiha_mosaic, xerr = sigma_niiha_acswfpc2, yerr = sigma_niiha_mosaic, fmt=None, zorder=-100)
scat = ax1.scatter(niiha_acswfpc2, niiha_mosaic, s=30, c=tab_mos['D'][mask])
for x, y, label in zip(niiha_acswfpc2, niiha_mosaic, tab_acs['Object'][mask]):
    ax1.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
ax1.plot(x_grid_log, p(x_grid_log), 'k-', linewidth=1.4)
cb = fig.colorbar(scat, ax=ax1)
cb.set_label('Distance from Trapezium, arcsec')
#ax1.plot(niiha_acswfpc2[m], np.polyval(a, niiha_acswfpc2[m]),'r-')
#ax1.plot(niiha_acswfpc2[m], a[0]*niiha_acswfpc2[m]+a[1],'b-')
#ax1.plot(D_arcmin1, niiha_mosaic, 'bo', label = 'Mosaic-shell')
#ax1.errorbar(D_arcmin1, niiha_mosaic,  yerr = sigma_niiha_mosaic, marker='o', fmt='bo')
#ax1.set_xscale('log')
#ax1.set_yscale('log')
#fig.tight_layout()
#plt.savefig('brightness-mosaic-vs-acswfpc2-shell.pdf')

fig = plt.figure(figsize=(8, 6.8))
ax2 = fig.add_subplot(1,1,1)
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlim(xmin=0.0, xmax=0.4)
ax2.set_ylim(ymin=0.0, ymax=0.4)
ax2.set_xlabel(r'Brightness acs and wfpc2, S = S(NII)/S(Halpha)')
ax2.set_ylabel(r'Brightness mosaic, S = S(NII)/S(Halpha)')
ax2.errorbar(Bg_niiha_acswfpc2, Bg_niiha_mosaic, xerr = Bg_sigma_niiha_acswfpc2, yerr = Bg_sigma_niiha_mosaic, fmt=None, zorder=-100 )
scatt = ax2.scatter(Bg_niiha_acswfpc2, Bg_niiha_mosaic, s=30, c=tab_mos['D'][mask_Bg])
for x, y, label in zip(Bg_niiha_acswfpc2, Bg_niiha_mosaic, tab_mos['Object'][mask_Bg]):
    ax2.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
#ax2.plot(Bg_niiha_acswfpc2, np.polyval(aa, Bg_niiha_acswfpc2),'b-')
ax2.plot(xx_grid_log, pp(xx_grid_log),'k-', linewidth=1.4)
cb = fig.colorbar(scatt, ax=ax2)
cb.set_label('Distance from Trapezium, arcsec')
#ax2.set_xscale('log')
#ax2.set_yscale('log')
#fig.tight_layout()
#plt.savefig('brightness-mosaic-vs-acswfpc2-bg.pdf')

fig = plt.figure()
ax3 = fig.add_subplot(1,1,1)
#ax3.set_xlim(xmin=0.6, xmax=8.5)
#ax3.set_ylim(ymin=-1.0, ymax=2.7)
ax3.set_xlim(xmin=0.06, xmax=15.0)
ax3.set_ylim(ymin=0.004, ymax=250.0)
#ax3.set_ylim(ymin=0.0001, ymax=200.0)
ax3.set_xlabel(r'D, arcmin')
ax3.set_ylabel(r'Brightness ratios NII and NII mosaic-WFPC2, S = S(NII)/S(Halpha)')
ax3.plot(D_arcmin, niiha_acswfpc2, 'ro', label = 'Shell-brightness')
ax3.errorbar(D_arcmin, niiha_acswfpc2,  yerr = sigma_niiha_acswfpc2, marker='o', fmt='ro')
ax3.plot(D_arcmin_Bg, Bg_niiha_acswfpc2, 'bo', label = 'Background-brightness')
ax3.errorbar(D_arcmin_Bg, Bg_niiha_acswfpc2,  yerr = Bg_sigma_niiha_acswfpc2, marker='o', fmt='bo')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.legend()
ax3.grid(True)
plt.savefig('acs_wfpc2-ratio-Nii_ha-vs-D_-mean-error_new.pdf')

fig = plt.figure()
ax4 = fig.add_subplot(1,1,1)
ax4.set_xlim(xmin=0.06, xmax=8.5)
ax4.set_ylim(ymin=0.0, ymax=0.7)
ax4.set_xlabel(r'D, arcmin')
ax4.set_ylabel(r'Brightness fraction NII and Halpha ACS-WFPC2 , S = S(NII)/S(Halpha)')
ax4.plot(D_arcmin1, niiha_mosaic, 'ro', label = 'Shell-brigthness')
ax4.errorbar(D_arcmin1, niiha_mosaic, yerr = sigma_niiha_mosaic, marker='o', fmt='ro')
ax4.plot(D_arcmin1_Bg, Bg_niiha_mosaic, 'bo', label = 'background-brightness')
ax4.errorbar(D_arcmin1_Bg, Bg_niiha_mosaic,  yerr = Bg_sigma_niiha_mosaic, marker='o', fmt='bo')
ax4.set_xscale('log')
#ax4.set_yscale('log')
ax4.legend(loc=2)
ax4.grid(True)
#plt.savefig('wfpc2-mosaic-ratio-Nii_ha-vs-D_-mean-error_new.pdf')

