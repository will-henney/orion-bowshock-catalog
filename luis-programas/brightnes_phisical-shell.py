'''
Estimate of surface brightness erg/s/cm^2/sr for shell ACS
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan) )

# photometric keywords from the ACS header
Sfactor_ACS = 0.0025030687604156482

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')
with open("interproplyd.txt") as f:
    problem_sources += f.read().split('\n')

print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

Distance = tab['D']
Dif_Bally = tab['Dif_Bally']

SBally_physical = Sfactor_ACS*Dif_Bally

#print label
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_xlim(xmin=100,xmax=600)
#ax1.set_ylim(ymin=0,ymax=10)
for x, y, s, e in zip(Distance[m], SBally_physical[m], label[m], Sfactor_ACS*tab['Delta'][m]):
    if e < y:
        ax1.plot(x,y,'bo')
        ax1.errorbar(x, y, yerr=e, c='b')
        size = 3
    else:
        ax1.plot(x,y,'r.')
        size = 2
    ax1.annotate(s, (x, y), alpha=0.5, size=size,
                   xytext=(-3,3), textcoords='offset points', ha='right', va='bottom',)

ax1.plot(Distance, tab["Value_bg_Bally"]*Sfactor_ACS, 'k.', alpha=0.4)

ax1.set_xlabel(r'$D$, arcsec')
ax1.set_ylabel(r'$S(\mathrm{H\alpha+NII})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_title(r'fraction Difference shell   y   background vs D')
ax1.grid(True)

fig.savefig("S(Ha+NII)_ACSshell.pdf")
