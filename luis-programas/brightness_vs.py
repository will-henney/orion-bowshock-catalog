'''
Graphic of brightness of Robberto in function of brightness of Bally
'''
import pylab
import math
import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.text as txt
from  astropy.table import Table

tab = Table.read("arcs-summary.tab", format="ascii.commented_header", delimiter="\t")

label = tab['Object']

#MASK
m = (tab['Dif_Bally']!='-') and (tab['Dif_Robberto_f656']!='-')
x = np.array(tab['Dif_Bally'][m], dtype = float)
y = np.array(tab['Dif_Robberto_f656'][m], dtype = float)

#Fraction of brightness
#fract = y/x
#fract_av = np.median(fract) #mediana
print x, y

#Fitting a straight line
m1 = (x>0.0)&(y>0.0)
x_log = np.log10(x[m1])
y_log = np.log10(y[m1])
a = np.polyfit(x_log, y_log, 1)
p = np.poly1d(a)
x_grid_log = np.linspace(x_log.min(), x_log.max())

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_xlim(xmin=100,xmax=600)
#ax1.set_ylim(ymin=0,ymax=10)
#for x_t,y_t,s in zip(x, y, label):
    #ax1.plot(x_t,y_t,'bo')
    #ax1.annotate(s, (x_t, y_t), alpha=0.8, size=8,
                   #xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)
ax1.set_xlabel(r'$\mathrm{Value(shell - background)}$, $\mathrm{electrons}\ s^{-1}$, (ACS F658N)')
ax1.set_ylabel(r'$\mathrm{Value(shell - background)}$, $\mathrm{electrons}\ s^{-1}$, (WFPC2 F656N)')
ax1.plot(x, y,'bo')
ax1.plot(10**(x_grid_log), 10**(p(x_grid_log)),'r-')
#ax1.plot(x,fract, 'bo')
#ax1.plot([0.01,40], [fract_av, fract_av], 'r-', linewidth=2)
#ax1.set_title(r'Mediana. Brightness of source of Bally and Robberto')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True)

plt.savefig('dif_brightness_vsT.pdf')

plt.show()

#fig.savefig("fit_dif_brightness_vs.pdf")
