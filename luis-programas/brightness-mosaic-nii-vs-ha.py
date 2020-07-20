'''
Graphic of nii vs ha of the mosaic
''' 
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab_mos = Table.read("arc-summary-mosaic.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

mask = np.isfinite(tab_mos['S(shell_mean_wfpc2)']) & np.isfinite(tab_mos['S(shell_mosaicf656)'])
mask_Bg = np.isfinite(tab_mos['S(Bg_wfpc2)']) & np.isfinite(tab_mos['S(Bg_mosaicf656)'])
S(Bg_wfpc2)'
#Shell
S_ha=tab_mos['S(shell_mosaicf658)'][mask]
S_ha=tab_mos['S(shell_mosaicf656)'] [mask]

# Background
S_nii_Bg=tab_mos['S(Bg_mosaicf658)'][mask_Bg]
S_ha_Bg=tab_mos['S(Bg_mosaicf656)'][mask_Bg]

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(xmin=-0.0027, xmax=0.008)
ax1.set_ylim(ymin=-0.008, ymax=0.025)
ax1.set_xlabel(r'Brightness NII f658n, erg/s/cm2/sr')
ax1.set_ylabel(r'Brightness Halpha f656n, erg/s/cm2/sr')
#ax1.plot(D_arcmin, niiha_acswfpc2, 'ro', label = 'ACS-WFPC2-Shell')
scat = ax1.scatter(S_nii, S_ha, s=30, c=tab_mos['D'][mask])
for x, y, label in zip(S_nii, S_ha, tab_mos['Object'][mask]):
    ax1.annotate(label, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
cb = fig.colorbar(scat, ax=ax1)
cb.set_label('Distance from Trapezium, arcsec')
#ax1.set_xscale('log')
#ax1.set_yscale('log')
plt.grid(True)
plt.savefig('brightness-mosaic-nii-vs-ha-shell.pdf')

fig = plt.figure()
ax2 = fig.add_subplot(1,1,1)
ax2.set_aspect('equal', adjustable='box')
ax2.set_xlim(xmin=-0.008, xmax=0.046)
ax2.set_ylim(ymin=-0.03, ymax=0.20)
ax2.set_xlabel(r'Brightness NII f658n, erg/s/cm2/sr')
ax2.set_ylabel(r'Brightness Halpha f656n, erg/s/cm2/sr')
scatt = ax2.scatter(S_nii_Bg, S_ha_Bg, s=30, c=tab_mos['D'][mask_Bg])
cb = fig.colorbar(scatt, ax=ax2)
cb.set_label('Distance from Trapezium, arcsec')
ax2.set_xscale('log')
ax2.set_yscale('log')
plt.grid(True)
plt.savefig('brightness-mosaic-nii-vs-ha-bg.pdf')


