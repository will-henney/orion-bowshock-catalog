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

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t")

# photometric keywords from the ACS header
Sfactor_ACS = 0.0025030687604156482

# ''            ''         ''     WFPC2 

Sfactor_WFPC2 = 0.16895605909108011

# O'Dell Calibration constant
#Sfactor_Odell = 0.1790544428870778

label = tab['Object']

#Value and Distance 
m = (tab['Value_bg_Bally']!='-') and (tab['Value_bg_Robberto_f656']!='-')
Distance = tab['D'][m]
Bally_value = tab['Value_bg_Bally'][m]
Robberto_value = np.array(tab['Value_bg_Robberto_f656'][m], dtype =float)

#multiply image pixel values by this to give surface brightness
#Distance_arcmin = (float(1)/60)*Distance 
SBally_value = 0.8*Sfactor_ACS*Bally_value
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
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
for x,y,s in zip(SBally_value, SRobberto_value, label):
    ax1.plot(x,y,'bo')
    #ax1.plot(x_grid_log, p(x_grid_log),'r-')
    ax1.annotate(s, (x, y), alpha=0.8, size=6,
                   xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)
ax1.set_xlabel(r'$\log$ $S(\mathrm{H\alpha + NII})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$ ACS')
ax1.set_ylabel(r'$\log$ $S(\mathrm{H\alpha})$, $\mathrm{erg\ s^{-1}\ cm^{-2}\ sr^{-1}}$   WFPC2')
#ax1.plot(Distance, SBally_value, 'bo')
ax1.plot(10**(x_grid_log), 10**(p(x_grid_log)),'r-')
#ax1.set_title(r"Surface brightness background of WFPC2 (F656N) in function of ACS (F658N), with ODell's calibration constant")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True)

#plt.savefig('S(alpha)_bg_WFPC2_ACS_calibrationlogT.pdf')

plt.show()



