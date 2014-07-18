from __future__ import print_function
import glob
import json
import pylab
import random
import numpy
import matplotlib.pyplot as plt
import matplotlib.text as txt

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)
D = []
Bg = []
f_diffe = []
label = []  
  
        label.append(data["star"]["id"])
        D.append(data["star"]["D"])
        print(D)

        b = data.keys()
        for k in b:
        if k.startswith('Bally'):
            imagename=k
            try:
                shell = data[imagename]["shell"]["value"]
                bg = data[imagename]["background"]["value"]
                dbg = data[imagename]["background"]["delta"]
            except KeyError:
                shell = 0.0
                bg = 0.0
                dbg = 0.0
            # difference between sheel and background 
            try:
                difference = float(shell-bg)/bg
            except ZeroDivisionError:
                difference = 0.0            
                f_diffe.append(difference) 
                Bg.append(bg) 
#a = numpy.polyfit(D,f_diffe,1)
#print('Ajuste lineal:', a)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_xlim(xmin=100,xmax=600)
#ax1.set_ylim(ymin=0,ymax=10)
for x,y,s in zip(D,Bg,label):
    ax1.plot(x,y,'bo')
    ax1.annotate(s, (x, y), alpha=0.8, size=5,
                   xytext=(-2,2), textcoords='offset points', ha='right', va='bottom',)

ax1.set_xlabel(r'$D$, arcsec')
ax1.set_ylabel(r'Value, $background$')

#ax1.set_title(r'Difference shell and  background vs D')
ax1.grid(True)

fig.savefig("brightness-bgT.pdf")
