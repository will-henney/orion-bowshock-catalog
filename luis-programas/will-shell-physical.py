
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table
import json

# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0025030687604156482

# Value of recombination coefficient in cm^3/s
alpha = 2.6e-13 

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

# Stellar wind
yr = 3.15576e7
Msun = 1.989e33
km = 1.0e5
Mdot_wind = 3.5e-7*Msun/yr # g/s 
Vwind = 1200*km     # cm/s
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
def draw_inclination_arrow(x0, y0, inc=45):
    """Show effects of inclination changes on x-axis"""
    x1 = x0/np.cos(np.radians(inc))
    plt.plot([x0, x1], [y0, y0], '-k')
    plt.plot(x0, y0, 'ok', ms=3.5)
    plt.plot(x1, y0, '>k', ms=3.5)
    plt.text(x1, y0, '   i = {}'.format(inc), va='center', fontsize='x-small')
tab = Table.read("arcs-summary-merge.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )
with open("problem-sources.txt") as f:
    problem_sources = f.read().split('\n')
with open("interproplyd.txt") as f:
    problem_sources += f.read().split('\n')

m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) & (tab['R_out'] > tab['R_in'])
m = m & np.array([not source in problem_sources for source in tab['Object']])
A = tab['Rc_out']/tab['R_out']
Ain = tab['Rc_in']/tab['R_in']
H = (tab['R_out'] - tab['R_in'])/tab['R_out']
D60 = tab['D']/60
contrast = np.log10(tab['Dif_Bally']/tab['Value_bg_Bally'])
q = tab['R_out'].data/tab['D'].data
PA_star = (tab['PA_star'] - 180.0) % 360.0
dPA = ((tab ['PA_out'] - tab ['PA_star'] + 180.0) % 360.0) - 180.0
with open("../ll-data.json") as f:
    db = json.load(f)

is_proplyd = {}
for data in db.values():
    if "LL" in data["bname"]:
        key = data["bname"]
    elif data["oname"]:
        key = data["oname"].split()[-1]
    else:
        continue
    if "proplyd?" in data["Notes"]:
        is_proplyd[key] = 0
    elif "proplyd" in data["Notes"]:
        is_proplyd[key] = 1
    else:
        is_proplyd[key] = -1

with open("extinction.json") as f:
    extinction_data = json.load(f)

# Ha surface brightness, corrected for extinction
Sha = Sfactor_ACS*tab['Dif_Bally']
Chb = np.array([extinction_data.get(source, 0.0) for source in tab['Object']])
Sha *= 10**(fha*Chb)

# Thickness and radius of the shell for measurements of delta l
h0 = tab['h']*cm_per_arcsec
rc = tab['Rc_out']*cm_per_arcsec
deltal = 2*np.sqrt(h0*rc)

nshell = np.sqrt(4.*np.pi*Sha/(alpha*deltal*Eha))
pshell = 2.0*nshell*k*T
MdotV_in = pshell*4.*np.pi*(tab['R_in']*cm_per_arcsec)**2 *yr/Msun/km

D60_grid = np.linspace(D60.min(), D60.max(), 2)
Dcm_grid = 60*D60_grid*cm_per_arcsec
Pram = Mdot_wind*Vwind/(4.*np.pi*Dcm_grid**2)

proplyd_mask = np.array([is_proplyd[source] == 1 for source in tab['Object']])
not_proplyd_mask = np.array([is_proplyd[source] == -1 for source in tab['Object']])
maybe_proplyd_mask = np.array([is_proplyd[source] == 0 for source in tab['Object']])

figlist = []

pltfile = 'will-nshell-vs-D.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], nshell[m], s=10*deltal[m]/cm_per_arcsec, c=np.log10(Sha[m]), cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, nshell, (nshell > 3500.0/D60) | (nshell < 1000.0/D60), allmask=m)
cb = plt.colorbar()
draw_inclination_arrow(0.1, 100*1.3**3, 30)
draw_inclination_arrow(0.1, 100*1.3**2, 45)
draw_inclination_arrow(0.1, 100*1.3, 60)
draw_inclination_arrow(0.1, 100, 75)
plt.text(0.1, 100/1.4, 'True distance correction', fontsize='x-small')
cb.set_label('H alpha surface brightness, erg/s/cm2/sr')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Shell electron density, ne / pcc ')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.05, 20.0)
ax.set_ylim(50.0, 5.e4)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-Pshell-vs-D.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], pshell[m], s=10*deltal[m]/cm_per_arcsec, c=np.log10(Sha[m]), cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, pshell, allmask=m)
cb = plt.colorbar()
cb.set_label('H alpha surface brightness, erg/s/cm2/sr')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Shell pressure, P / dynes/cm^2 ')
plt.plot(D60_grid, Pram, '-k')
plt.plot(D60_grid*np.cos(np.radians(30)), Pram, '--k')
plt.plot(D60_grid*np.cos(np.radians(60)), Pram, ':k')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.05, 20.0)
ax.set_ylim(2e-10, 1.5e-7)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-MdotV-vs-D.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
mm = m & proplyd_mask
plt.scatter(D60[mm], MdotV_in[mm], s=10*deltal[mm]/cm_per_arcsec, c='red', alpha=0.6)
mm = m & maybe_proplyd_mask
plt.scatter(D60[mm], MdotV_in[mm], s=10*deltal[mm]/cm_per_arcsec, c='orange', alpha=0.6)
mm = m & not_proplyd_mask
plt.scatter(D60[mm], MdotV_in[mm], s=10*deltal[mm]/cm_per_arcsec, c='black', alpha=0.6)
label_sources(tab['Object'], D60, MdotV_in, not_proplyd_mask | (MdotV_in > 1.e-6/D60), allmask=m)
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Inner flow Mdot V, Msun/yr km/s')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.05, 20.0)
ax.set_ylim(6e-9, 3e-5)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))
