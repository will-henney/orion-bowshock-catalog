'''
Estimate of the fration between Halpha and NII
''' 

import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0025030687604156482
Sfactor_WFPC2 = 0.16895605909108011

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')

with open("interproplyd.txt") as l:
    problem_sources += l.read().split('\n')

#print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

# Observed and corrected background brightness values
D_arcsec = tab['D']
SHaNII = Sfactor_ACS*tab['Value_bg_Bally']
SHa = Sfactor_WFPC2*tab['Value_bg_Robberto_f656']
ShaNII_shell = Sfactor_ACS*tab['Dif_Bally']
Sha_shell = Sfactor_WFPC2*tab['Dif_Robberto_f656']

# Ratio between S(Ha + NII) and S(Ha)

NIIsHa = (SHaNII/SHa) - 1
NIIsHa_shell = (ShaNII_shell/Sha_shell) - 1
#print(NIIsHa, SRobberto_value)
Bg = np.where(NIIsHa > 0.0, NIIsHa, 0.03)
Shell = np.where(NIIsHa_shell > 0.0, NIIsHa_shell, 0.03)
# for NIIsHa_a in NIIsHa[m]:
#     if NIIsHa_a > 0.0:
#         Bg_NIIsHa = NIIsHa_a
#     else:
#         Bg_NIIsHa = 0.03         
#     Bg.append(Bg_NIIsHa)  

# for NIIsHa_shell_a in NIIsHa_shell[m]:
#     if NIIsHa_shell_a > 0.0:
#         s_NIIsHa_shell = NIIsHa_shell_a
#     else:
#         s_NIIsHa_shell = 0.03         
#     Shell.append(s_NIIsHa_shell
   
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(xmin=0.06, xmax=15.0)
ax1.set_ylim(ymin=0.03, ymax=50.0)
for x, y, z, s in zip(D_arcsec[m]/60.0, Bg[m], Shell[m], label[m]):
    #ax1.plot(x_grid_log, p(x_grid_log),'r-')
    ax1.annotate(s, (x, y), alpha=0.5, size=3,
                   xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)
    ax1.annotate(s, (x, z), alpha=0.5, size=3,
                   xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)
ax1.plot(D_arcsec[m]/60.0, Bg[m],'ro', alpha=0.4, label='background-brightness')
ax1.plot(D_arcsec[m]/60.0, Shell[m] ,'ko', alpha=0.6, label='Shell-brightness')
ax1.set_xlabel(r'D, arcmin')
ax1.set_ylabel(r'Brightness fraction NII and Halpha, S = S(NII)/S(Halpha)')
ax1.set_xscale('log')
#ax1.set_yscale('symlog', linthreshy=0.1)
ax1.set_yscale('log')
ax1.legend()
ax1.grid(True)

plt.savefig('S-NII_Ha_vs_D.pdf')


