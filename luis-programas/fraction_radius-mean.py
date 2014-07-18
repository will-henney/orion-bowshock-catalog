'''
Fraction of radius (R_out/R_in)
'''
import sys
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

tab = Table.read("arcs-summary.tab",  format="ascii.commented_header", delimiter="\t")

m = (tab['R_out']!='-')&(tab['R_in']!='-')
D_arcmin = tab['D'][m]/60.0
R_out = np.array(tab['R_out'][m], dtype = float)
R_in = np.array(tab['R_in'][m], dtype= float)

m1 = (tab['Rc_out']!='-')&(tab['Rc_in']!='-')
Rc_out = np.array(tab['Rc_out'][m1], dtype = float)
Rc_in = np.array(tab['Rc_in'][m1], dtype = float)
# fraction between the radius
f = R_out/R_in
#f_mean = np.mean(f)
bins = np.linspace(min(D_arcmin), max(D_arcmin), 6)
digitized = np.digitize(D_arcmin, bins)-1
bin_centers = bins + (bins[1]-bins[0])/2
bin_means = [np.mean(f[digitized == j]) for j in range(len(bins))]

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r'$\mathrm{R_{0}(out)}/\mathrm{ R_{0}(in)}$')
#ax1.set_title(r'Fraction of radius and mean per bins')
ax1.plot(D_arcmin, f, 'ro')
#ax1.plot(bin_centers, bin_means)
#ax1.plot([0.01, 10], [f_mean, f_mean], 'r-', linewidth=2)
#ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.grid(True)

plt.savefig("fraction-radius_mean2.pdf")

