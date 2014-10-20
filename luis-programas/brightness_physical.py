'''
Estimate of surface brightness erg/s/cm^2/sr
'''
import pylab
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.text as txt
from  astropy.table import Table

# photometric keywords from the ACS header
Sfactor_ACS = 0.0025030687604156482

# ''            ''         ''     WFPC2 

Sfactor_WFPC2 = 0.16895605909108011

# O'Dell Calibration constant
#Sfactor_Odell = 0.1790544428870778


tab = Table.read("arcs-summary-merge.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')
with open("interproplyd.txt") as f:
    problem_sources += f.read().split('\n')

#m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) & (tab['R_out'] > tab['R_in'])
#m = m & np.array([not source in problem_sources for source in tab['Object']])

label = tab['Object']

#Value and Distance 
Distance = tab['D']
Bally_value = tab['Value_bg_Bally']
Robberto_value = np.array(tab['Value_bg_Robberto_f656'], dtype =float)

#multiply image pixel values by this to give surface brightness
#Distance_arcmin = (float(1)/60)*Distance 
SBally_value = Sfactor_ACS*Bally_value
SRobberto_value =  Sfactor_WFPC2*Robberto_value

#Ajuste
m1 = (SBally_value>0.0)&(SRobberto_value>0.0)
Distance_log = np.log10(Distance[m1])
Bally_value_log = np.log10(SBally_value[m1])
Robberto_value_log = np.log10(SRobberto_value[m1])
a = np.polyfit( Bally_value_log, Robberto_value_log, 1)
p = np.poly1d(a)
x_grid_log = np.linspace(Bally_value_log.min(), Bally_value_log.max())
print a

#relation line 1:1
S_ha_wfpc = np.linspace(SBally_value.min(), SBally_value.max())


fig = plt.figure(figsize=(8, 6.8))
ax1 = fig.add_subplot(1,1,1, axisbg="#f5f5dc")
ax1.set_aspect('equal', adjustable='box')
ax1.set_xlim(xmin=0.001, xmax=10**-1)
ax1.set_ylim(ymin=0.001, ymax=10**-1)
for x,y,s in zip(SBally_value, SRobberto_value, label):
    ax1.annotate(s, (x, y), fontsize='x-small', xytext=(5, 5), textcoords='offset points', bbox={"boxstyle": "round", "fc": "white", "ec": "none", "alpha": 0.5}, alpha=0.7)
ax1.set_xlabel(r'$\log$ $S(\mathrm{H\alpha + NII})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$ ACS-F658N')
ax1.set_ylabel(r'$\log$ $S(\mathrm{H\alpha})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$   WFPC2-F656N')
#ax1.plot(Distance, SBally_value, 'bo')
scat = ax1.scatter(SBally_value, SRobberto_value, s=30, c=tab['D'])
cb = fig.colorbar(scat, ax=ax1)
cb.set_label('Distance from Trapezium, arcsec')
ax1.plot(10**(x_grid_log), 10**(p(x_grid_log)),'k-')
ax1.plot(S_ha_wfpc, S_ha_wfpc, 'k--')
#ax1.set_title(r"Surface brightness background of WFPC2 (F656N) in function of ACS (F658N), with ODell's calibration constant")
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.grid(True)

plt.savefig('S(alpha)_bg_WFPC2_ACS_calibrationlog-f.pdf')

#plt.show()



