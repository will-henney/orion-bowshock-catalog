'''	  	
Estimation of the pressure
'''
import numpy as np
import json
from  astropy.table import Table, Column
import matplotlib.pyplot as plt

# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0025030687604156482

# Value of alpha in cm^3/s
alpha = 1.27e-13 

# Energy to 3 to 2... erg
E = 3.0267338723714944e-12 

# Constant of Bolzman in erg/K
k = 1.38e-16 

# Temperatura
T = 1e4
M_o = 3.5e-7 #Msun/yr
v = 1200     #km/s

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

with open("extinction.json") as f:
    extinction_data = json.load(f)

tab.add_column(Column(name='Correction',
                      data=[10**(0.78*extinction_data.get(source, 0.0))
                            for source in tab['Object']]))

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')

with open("interproplyd.txt") as l:
    problem_sources += l.read().split('\n')

print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

# Observed and corrected background brightness values

D_arcsec = tab['D']
SHa_obs = Sfactor_ACS*tab['Dif_Bally']
SHa = SHa_obs*tab['Correction']

# Thickness and radius of the shell for measurements of delta l
h0 = tab['h']
rc = tab['Rc_out']
deltal = 2*np.sqrt(h0*rc)

n = np.sqrt((4*np.pi*SHa)/(deltal*alpha*E*436*1.49597870691e13))

# Pressure termal shell
Pc = 2*n*k*T

# pressure external
Pex = (M_o*v*1e5*3.0e33)/(4*np.pi*3.1557e7*(436*D_arcsec*1.49597870691e13)**2)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlim(xmin=0.06, xmax=20.0)
ax1.set_ylim(ymin=2e-11, ymax=2e-7) 
ax1.plot(D_arcsec[m]/60.0, Pc[m],'ro', label='Pressure-shell')
ax1.plot(D_arcsec/60.0, Pex, 'k-', alpha=0.6, label='Pressure-wind')
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r' $\mathrm{Pressure}$, $\mathrm{g\ cm^{-1}\ s^{-2}}$')
#ax1.set_title(r"Surface brightness background of WFPC2 (F656N) in function of ACS (F658N), with ODell's calibration constant")
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend()
ax1.grid(True)

plt.savefig('pressure_correted_extintion_ACS_vs_D.pdf')

