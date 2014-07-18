'''
fraction between h and R0 versus D
'''
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t")

label = tab['Object']
m = (tab['R_out']!='-')
h = tab['h'][m]
R0 = np.array(tab['R_out'][m], dtype=float)        
D_arcmin = tab['D'][m]/60.0

#fraction betwen the width (h) and out radius (R0) and mean
f = h/R0
f_mean = np.mean(f)
f_av = np.median(f)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
for x,y,s in zip(D_arcmin, f, label):
    ax1.plot(x,y,'bo')
    ax1.annotate(s, (x, y), alpha=0.8, size=5,
                   xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r'$\mathrm{h}/\mathrm{ R_{out}}$')
ax1.set_title(r'Width (h) of the shells and its mean')
ax1.plot(D_arcmin, f, 'bo')
ax1.plot([0.01,15], [f_mean, f_mean], 'r-', linewidth=2)
ax1.plot([0.01,15], [f_av, f_av], 'g-', linewidth=2)
ax1.set_xscale('log')
#ax1.set_yscale('log')
ax1.grid(True)

plt.savefig('h(width)_shell_vs_D.pdf')
plt.show()
