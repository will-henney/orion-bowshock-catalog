'''
Ratio between ACS and WFPC2 (histogram)
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

#Fsctores of calibration
Sfactor_ACS = 0.0025030687604156482
Sfactor_WFPC2 = 0.16895605909108011

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')

print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

Bally = Sfactor_ACS*tab['Value_bg_Bally']
Robberto = Sfactor_WFPC2*tab['Value_bg_Robberto_f656']

Ratio = Bally/Robberto

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$S(\mathrm{H\alpha + NII})/ S(\mathrm{H\alpha})$')
ax1.set_ylabel(r'$N$')
ax1.set_title(r'Histogram of ratio brightnes')
ax1.hist(Ratio[m], 50, color='blue')
ax1.grid(True)

plt.savefig("histogram-ratio_S_ACS-WFPC2_babgraund.pdf")

