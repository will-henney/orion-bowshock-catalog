import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan) )

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')
with open("interproplyd.txt") as f:
    problem_sources = f.read().split('\n')

print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

Distance = tab['D']
Dif_Bally = tab['Dif_Bally']
#Dif_Robberto = np.array(tab['Dif_Robberto_f656'][m], dtype=float)
  
#print label
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_xlim(xmin=100,xmax=600)
#ax1.set_ylim(ymin=0,ymax=10)
for x, y, s, e in zip(Distance[m], Dif_Bally[m], label[m], tab['Delta'][m]):
    if e < y:
        ax1.plot(x,y,'bo')
        ax1.errorbar(x, y, yerr=e, c='b')
        size = 3
    else:
        ax1.plot(x,y,'r.')
        size = 2
    ax1.annotate(s, (x, y), alpha=0.5, size=size,
                   xytext=(-3,3), textcoords='offset points', ha='right', va='bottom',)

ax1.plot(Distance, tab["Value_bg_Bally"], 'k.', alpha=0.4)

ax1.set_xlabel(r'$D$, arcsec')
ax1.set_ylabel(r'$\mathrm{Value(shell - background)}$, $\mathrm{electrones\ s^{-1}}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_title(r'fraction Difference shell   y   background vs D')
ax1.grid(True)

#fig.savefig("brightness_ACSshellT.pdf")
plt.show()
