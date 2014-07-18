'''
Estimate density
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

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('-', np.nan))

with open("extinction.json") as f:
    extinction_data = json.load(f)

tab.add_column(Column(name='Correction',
                      data=[10**(0.78*extinction_data.get(source, 0.0))
                            for source in tab['Object']]))

with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')

print(problem_sources)
label = tab['Object']
m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) 
m = m & np.array([not source in problem_sources for source in tab['Object']])

# Observed and corrected background brightness values

D_arcmin = tab['D'][m]/60.0
SHa_obs = Sfactor_ACS*tab['Dif_Bally'][m]
SHa = SHa_obs*tab['Correction'][m]

# Limits for plot
#xmin, xmax = 0.06, 20.0

# Thickness and radius of the shell for measurements of delta l
h0 = tab['h'][m]
rc = tab['Rc_out'][m]
deltal = 2*np.sqrt(h0*rc)
#for a, b in zip(tab['Object'], deltal):
    #print a, b 
n = np.sqrt((4*np.pi*SHa)/(deltal*alpha*E*436*1.49597870691e13))
for x, y in zip(tab['Object'], n):
    print x, y

# fit of line
m1 = (D_arcmin>0.0)&(np.isfinite(n))&(n>0.0)
Distance_log = np.log10(D_arcmin[m1])
n_log = np.log10(n[m1])
a = np.polyfit(Distance_log, n_log, 1)
p = np.poly1d(a)
x_grid_log = np.linspace(Distance_log.min(), Distance_log.max())


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel(r'$D$, arcmin')
ax1.set_ylabel(r'Density, $\mathrm{ cm^{-3}}$')
ax1.set_title(r'Density of the shell (ACS  f658n)')
ax1.plot(D_arcmin, n, 'bo')
ax1.plot(10**(x_grid_log), 10**(p(x_grid_log)),'r-')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_xlim(xmin, xmax)
ax1.grid(True)

#plt.savefig('density_shell_acs_corrected.pdf')

               
