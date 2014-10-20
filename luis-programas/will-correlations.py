
import numpy as np
import matplotlib.pyplot as plt
from  astropy.table import Table

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
A = tab['Rc_out']/tab['R_out']
Ain = tab['Rc_in']/tab['R_in']
H = (tab['R_out'] - tab['R_in'])/tab['R_out']
D60 = tab['D']/60
contrast = np.log10(tab['Dif_Bally']/tab['Value_bg_Bally'])
q = tab['R_out'].data/tab['D'].data
PA_star = (tab['PA_star'] - 180.0) % 360.0
dPA = ((tab ['PA_out'] - tab ['PA_star'] + 180.0) % 360.0) - 180.0
 
figlist = []

pltfile = 'will-A-vs-q.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(q[m], A[m], s=20*tab['R_out'][m], c=D60[m], vmin=0.0, cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], q, A, allmask=m)
cb = plt.colorbar()
cb.set_label('Projected distance from Trapezium, D / arcmin')
plt.xlabel('Bowshock fractional size, q = r0/D')
plt.ylabel('Bowshock bluntness, A = Rc/r0')
ax.set_xscale('log')
ax.set_xlim(0.001, 1.0)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-H-vs-A.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(A[m], H[m], s=20*tab['R_out'][m], c=D60[m], vmin=0.0, cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], A, H, (A > 5.0) | (A < 1.0) | (H <= 0.2) | (H >= 0.6), allmask=m)
cb = plt.colorbar()
cb.set_label('Projected distance from Trapezium, D / arcmin')
plt.xlabel('Bowshock bluntness, A = Rc/r0')
plt.ylabel('Shell relative thickness, H = h/r0')
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-H-vs-q.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(q[m], H[m], s=20*tab['R_out'][m], c=np.log10(D60[m]), cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], q, H, allmask=m)
cb = plt.colorbar()
cb.set_label('Projected distance from Trapezium, D / arcmin')
plt.xlabel('Bowshock fractional size, q = r0/D')
plt.ylabel('Shell relative thickness, H = h/r0')
ax.set_xscale('log')
ax.set_xlim(0.001, 1.0)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-q-vs-D.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], q[m], s=20*tab['R_out'][m], c=contrast[m], cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, q, allmask=m)
cb = plt.colorbar()
cb.set_label('Shell/background brightness contrast')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Bowshock fractional size, q = r0/D')
ax.set_xlim(0.05, 20.0)
ax.set_ylim(0.001, 1.0)
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-r0-vs-D.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.scatter(D60[m], tab['R_out'][m], s=100*H[m], c=contrast[m], cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], D60, tab['R_out'],
              (tab['R_out'] >= 5.0) | (tab['R_out'] <= 0.8), allmask=m)
Darray = np.linspace(D60.min(), D60.max())
r0norm = 0.8
plt.plot(Darray, r0norm*Darray/Darray.min(), 'k--', zorder=-100)
plt.plot(Darray, r0norm*(Darray/Darray.min ())**0.25, 'k-', zorder=-100)
cb = plt.colorbar()
cb.set_label('Brightness contrast, log10(Shell / BG)')
plt.xlabel('Projected distance from Trapezium, D / arcmin')
plt.ylabel('Bowshock radius, r0 / arcsec')
plt.text(0.05, 0.05, 'Symbol size indicates shell relative thickness, H',
         transform=ax.transAxes, fontsize='x-small')
ax.set_xlim(0.05, 20.0)
ax.set_ylim(0.3, 11.0)
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))

pltfile = 'will-PA-vs-PA.pdf'
fig = plt.figure(figsize=(7,6))
ax = fig.add_subplot(111, axisbg="#eeeeee")
plt.fill_betweenx([-90.0, 90.0], [0.0, 0.0], [90.0, 90.0], zorder=-10, alpha=0.05)
plt.fill_betweenx([-90.0, 90.0], [180.0, 180.0], [270.0, 270.0], zorder=-10, alpha=0.05)
plt.fill_betweenx([-90.0, 90.0], [360.0, 360.0], [450.0, 450.0], zorder=-10, alpha=0.05)
plt.text(45.0, -80.0, 'NE\nquadrant',  ha='center', fontsize='x-small')
plt.text(135.0, -80.0, 'SE\nquadrant', ha='center', fontsize='x-small')
plt.text(225.0, -80.0, 'SW\nquadrant', ha='center', fontsize='x-small')
plt.text(315.0, -80.0, 'NW\nquadrant', ha='center', fontsize='x-small')
plt.axhline(zorder=-5)
plt.scatter(PA_star[m], dPA[m], s=20*tab['R_out'][m], c=D60[m], cmap=plt.cm.hot, alpha=0.6)
label_sources(tab['Object'], PA_star, dPA, np.abs(dPA) > 45.0, allmask=m)
cb = plt.colorbar()
cb.set_label('Projected distance from Trapezium, D / arcmin')
plt.xlabel('PA of source from Trapezium, deg')
plt.ylabel('Angle between bowshock axis and radial direction, deg')
ax.set_xlim(-30.0, 375.0)
ax.set_ylim(-90.0, 90.0)
fig.savefig(pltfile)
figlist.append('[[file:luis-programas/{0}][{0}]]'.format(pltfile))
