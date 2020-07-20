'''
Graphic of nii vs ha of the mosaic
''' 
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab_w = Table.read("arc-summary-acs-wfpc.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))
tab_mos = Table.read("arc-summary-mosaic.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

mask = np.isfinite(tab_w['S(shell_mean_wfpc2)']) & np.isfinite(tab_mos['S(shell_mosaicf656)'])
mask = mask & (tab_w['S(shell_mean_wfpc2)']>0.0) & (tab_mos['S(shell_mosaicf656)']>0.0)
mask_Bg = np.isfinite(tab_w['S(Bg_wfpc2)']) & np.isfinite(tab_mos['S(Bg_mosaicf656)'])
mask_Bg = mask & (tab_w['S(Bg_wfpc2)']>0.0) & (tab_mos['S(Bg_mosaicf656)']>0.0)

# Error
max_error = 0.8
mask = mask & (tab_w['sigma(wfpc2)'] <
               max_error*tab_w['S(shell_mean_wfpc2)'])
mask = mask & (tab_mos['sigma(shell_mosaicf656)'] <
               max_error*tab_mos['S(shell_mosaicf656)'])

# Shell
S_ha_wfpc = tab_w['S(shell_mean_wfpc2)'][mask]
S_ha_mos = tab_mos['S(shell_mosaicf656)'][mask]
sigma_S_ha_wfpc = tab_w['sigma(wfpc2)'][mask]
sigma_S_ha_mos = tab_mos['sigma(shell_mosaicf656)'][mask]

# Background
S_ha_wfpc_bg = tab_w['S(Bg_wfpc2)'][mask_Bg]
S_ha_mos_bg = tab_mos['S(Bg_mosaicf656)'][mask_Bg]
sigma_S_ha_wfpc_bg = tab_w['sigma(Bg_wfpc2)'][mask_Bg]
sigma_S_ha_mos_bg = tab_mos['sigma(Bg_mosaicf656)'][mask_Bg]

#relation line 1:1
S_ha_wfpc_grid = np.linspace(S_ha_wfpc.min(), S_ha_wfpc.max())
S_ha_wfpc_bg_grid = np.linspace(S_ha_wfpc_bg.min(), S_ha_wfpc_bg.max()) 

#Ratios of Halpha between mosaic f656 and wfpc2 f656
F = S_ha_mos_bg/S_ha_wfpc_bg

fig = plt.figure(figsize=(8, 6.8))
ax1 = fig.add_subplot(1,1,1)
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(xmin=-0.001, xmax=0.0337)
ax1.set_ylim(ymin=-0.001, ymax=0.0337)
ax1.set_xlabel(r'Brightness H alpha WFPC2 f656n , erg/s/cm2/sr')
ax1.set_ylabel(r'Brightness H alpha mosaic f656n, erg/s/cm2/sr')
#ax1.plot(D_arcmin, niiha_acswfpc2, 'ro', label = 'ACS-WFPC2-Shell')
ax1.errorbar(S_ha_wfpc, S_ha_mos, xerr = sigma_S_ha_wfpc, yerr = sigma_S_ha_mos, fmt=None, zorder=-100)
scat = ax1.scatter(S_ha_wfpc, S_ha_mos, s=30, c=tab_mos['D'][mask])
for x, y, label in zip(S_ha_wfpc, S_ha_mos, tab_w['Object'][mask]):
    ax1.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
ax1.plot(S_ha_wfpc_grid, S_ha_wfpc_grid, 'k-')
cb = fig.colorbar(scat, ax=ax1)
cb.set_label('Distance from Trapezium, arcsec')
#ax1.set_xscale('log')
#ax1.set_yscale('log')
plt.savefig('brightness-hamosaic-vs-harobberto-shell.pdf')

fig = plt.figure(figsize=(8, 6.8))
ax2 = fig.add_subplot(1,1,1)
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlim(xmin=-0.0, xmax=0.20)
ax2.set_ylim(ymin=-0.0, ymax=0.20)
ax2.set_xlabel(r'Brightness H alpha WFPC2 f656n, erg/s/cm2/sr')
ax2.set_ylabel(r'Brightness H alpha mosaic f656n, erg/s/cm2/sr')
ax2.errorbar(S_ha_wfpc_bg, S_ha_mos_bg, xerr = sigma_S_ha_wfpc_bg, yerr = sigma_S_ha_mos_bg, fmt=None, zorder=-100)
scatt = ax2.scatter(S_ha_wfpc_bg, S_ha_mos_bg, s=30, c=tab_mos['D'][mask_Bg])
for x, y, label in zip(S_ha_wfpc_bg, S_ha_mos_bg, tab_w['Object'][mask]):
    ax2.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
ax2.plot(S_ha_wfpc_bg_grid, S_ha_wfpc_bg_grid, 'k-')
cb = fig.colorbar(scatt, ax=ax2)
cb.set_label('Distance from Trapezium, arcsec')
#ax2.set_xscale('log')
#ax2.set_yscale('log')
plt.savefig('brightness-hamosaic-vs-harobberto-bg.pdf')

fig = plt.figure()
ax3 = fig.add_subplot(1,1,1)
#ax3.set_aspect('equal', adjustable='box')
ax3.set_xlim(xmin=-0.02, xmax=0.20)
ax3.set_ylim(ymin=0.85, ymax=1.50)
ax3.set_xlabel(r'Brightness H alpha Robberto f656n, erg/s/cm2/sr')
ax3.set_ylabel(r'Ratio  H alpha  , S =H alpha(mosaic f656)/H alpha(Robberto f656)')
#ax3.errorbar(S_ha_wfpc_bg, S_ha_mos_bg, xerr = sigma_S_ha_wfpc_bg, yerr = sigma_S_ha_mos_bg, fmt=None, zorder=-100)
scattt = ax3.scatter(S_ha_wfpc_bg, F, s=30, c=tab_mos['D'][mask_Bg])
for x, y, label in zip(S_ha_wfpc_bg, F, tab_w['Object'][mask]):
    ax3.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
#ax3.plot(S_ha_wfpc_bg_grid, 'k-')
cb = fig.colorbar(scattt, ax=ax3)
cb.set_label('Distance from Trapezium, arcsec')
plt.axhline(y=1.0, xmin=0, xmax=1, hold=None, zorder=-5)
#ax2.set_xscale('log')
#ax2.set_yscale('log')
plt.savefig('ratio-hamosaic-haRobber-vs-hawfpc2-bg.pdf')


