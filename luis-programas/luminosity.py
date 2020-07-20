'''
estimate rate of luminosity
'''
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table
import json

# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0025030687604156482

# Value of recombination coefficient in cm^3/s
alpha_B = 2.6e-13 
alpha_Ha = 1.27e-13

# Energy to 3 to 2... erg
Eha = 3.0267338723714944e-12 

# Relative extinction at Ha from Blagrave
fha = 0.78

# Distances 
D_orion_pc = 436.0
AU = 1.49597870691e13
cm_per_arcsec = D_orion_pc*AU

# Boltzmann in erg/K
k = 1.3806503e-16

# Temperatura
T = 1e4

# Mas of proton in g
mh = 1.67e-24
m_p = 1.3*mh

#velocidad del sonido 

c = np.sqrt(k*T/mh)
print c
#number of mac
M = 2.

# cool
gamma = 3.3e-24

# factor for estimate ratio between lengh
d_cool_h = ((M**2+3)/4)**2
print d_cool_h 

famous = ['177-341', '167-317', '168-326', '161-324', '142-301']
other_interesting = ['w073-227', 'w069-601', 'w266-558', 'w000-400', 'w005-514']

def label_sources(labels, x, y, extramask=None, allmask=None):
    """Add labels at points (x, y) for selected sources"""
    mask = np.array(['LL' in source
                     or source in famous + other_interesting
                     for source in labels])
    if extramask is not None:
        mask = mask | extramask
    if allmask is not None:
        mask = mask & allmask
    radius = 5.0
    theta_max = np.pi/3.0
    for i, txt in enumerate(labels[mask]):
        theta = (2*np.random.random_sample() - 1.0)*theta_max
        dx, dy = radius*np.cos(theta), radius*np.sin(theta)
        ax.annotate(txt, (x[mask][i], y[mask][i]), (dx, dy),
                    textcoords='offset points', verticalalignment='center',
                    fontsize=5, alpha=0.6)

tab = Table.read("arcs-summary-merge.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )
with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')
with open("interproplyd.txt") as f:
    problem_sources += f.read().split('\n')

m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) & (tab['R_out'] > tab['R_in'])
m = m & np.array([not source in problem_sources for source in tab['Object']])

D60 = tab['D']/60
H = (tab['R_out'] - tab['R_in'])/tab['R_out']

with open("extinction.json") as f:
    extinction_data = json.load(f)

# Ha surface brightness, corrected for extinction
Sha = Sfactor_ACS*tab['Dif_Bally']
Chb = np.array([extinction_data.get(source, 0.0) for source in tab['Object']])
Sha *= 10**(fha*Chb)

# Combined fit: Ratio = 0.28 D**0.43
Rnii_ha = 0.28*D60**0.43
Sha /= 1.0 + Rnii_ha

# Thickness and radius of the shell for measurements of delta l
h0 = tab['h']*cm_per_arcsec
rc = tab['Rc_out']*cm_per_arcsec
deltal = 2*np.sqrt(h0*rc)


nshell = np.sqrt(4.*np.pi*Sha/(alpha_Ha*deltal*Eha))

#ratio luninosity 
ratio_L = 0.5*m_p*M*c**3/(nshell*gamma*h0)

#estimate ratio between d_cool and h

Lambda_ratio = 0.25  # This value correspond to (T_1/T_0)^-a because Lambda is proporcional to T^a, with a = 2 

d_hratio = d_cool_h*ratio_L*Lambda_ratio

pltfile = 'luminosity-ratio.pdf'
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], ratio_L[m], s=100*H[m], c='blue', cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, ratio_L, allmask=m)
#cb = plt.colorbar()
#cb.set_label('H alpha surface brightness, erg/s/cm2/sr')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Ratio luminosity. L_shock/L_shell')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlim(0.06, 20)
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, color='r', hold=None, zorder=-5.0)
plt.text(0.05, 0.9, 'Symbol size indicates shell relative thickness, H',
         transform=ax.transAxes, fontsize='small')
#ax.set_ylim(-0.03e24, 0.25e24)
fig.savefig(pltfile)

pltfile = 'thickness-ratio.pdf'
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], d_hratio[m], s=100*H[m], c='red', cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, d_hratio, allmask=m)
#cb = plt.colorbar()
#cb.set_label('H alpha surface brightness, erg/s/cm2/sr')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Ratio thickness. Lcool/h')
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlim(0.06, 20)
plt.axhline(y=0.0, xmin=0.0, xmax=1.0, color='b', hold=None, zorder=-5.0)
plt.text(0.05, 0.9, 'Symbol size indicates shell relative thickness, H',
         transform=ax.transAxes, fontsize='small')
#ax.set_ylim(-0.03e24, 0.25e24)
fig.savefig(pltfile)

