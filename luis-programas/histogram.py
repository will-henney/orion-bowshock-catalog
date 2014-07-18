import numpy as np
import matplotlib.pyplot as plt
import math
import random

infile = 'arcs-summary.tab'

with open(infile,'r') as inputFile:
    for i,line in enumerate(inputFile):
        pass
nLines = i

r_out = []
rc_out = []
F = []
Rc =['0.5','1','1.5','2', '2.5', '3', '3.5', '4', '4.5','5', '5.5', '6', '6.5', '7', '7.5','8', '8.5', '9']
Number = []

with open(infile,'r') as inputFile:
    for j in range(int(nLines)):
        line = inputFile.readline()
        linelist = line.split()
        R_out = linelist[7].strip()
        Rc_out = linelist[9].strip()
        
        if (R_out != '-' and Rc_out !='-' ):
            r_out.append(R_out)
            rc_out.append(Rc_out)

r_out_p = r_out[1:]
rc_out_p = rc_out[1:]
for r_out, rc_out in zip(r_out_p, rc_out_p):
    r_out=float(r_out)
    rc_out=float(rc_out)

#fraction of the  radius
    fracc=float(rc_out)/r_out
    F.append(fracc)

cnat = len(F)
a = 0
b = 0
c = 0
d = 0
e = 0
f = 0
g = 0
h = 0
k = 0
l = 0
m = 0
n = 0
o = 0
p = 0
q = 0
r = 0
s = 0
t = 0
for i in F:
    if (0<=i<=0.5):
        a=a+1
    if (0.5<=i<=1):
        b=b+1
    if (1<=i<=1.5):
        c=c+1
    if (1.5<=i<=2):
        d=d+1
    if (2<=i<=2.5):
        e=e+1
    if (2.5<=i<=3):
        f=f+1
    if (3<=i<=3.5):
        g=g+1
    if (3.5<=i<=4):
        h=h+1
    if (4<=i<=4.5):
        k=k+1
    if (4.5<=i<=5):
        l=l+1
    if (5<=i<=5.5):
        m=m+1
    if (5.5<=i<=6):
        n=n+1
    if (6<=i<=6.5):
        o=o+1
    if (6.5<=i<=7):
        p=p+1
    if (7<=i<=7.5):
        q=q+1
    if (7.5<=i<=8):
        r=r+1
    if (8<=i<=8.5):
        s=s+1
    if (8.5<=i<=9):
        t=t+1
    if (9<=i<=9.5):
        u=u+1
Number.append(a)
Number.append(b)
Number.append(c)
Number.append(d)
Number.append(e)
Number.append(f)
Number.append(g)
Number.append(h)
Number.append(k)
Number.append(l)
Number.append(m)
Number.append(n)
Number.append(o)
Number.append(p)
Number.append(q)
Number.append(r)
Number.append(s)
Number.append(t)

print cnat ,Number

pos = np.arange(len(Rc))
width = 1.0     

ax = plt.axes()
ax.set_xticks(pos + (width / 2))
ax.set_xticklabels(Rc)
ax.set_xlabel(r'$Rc_{out}/R_{out}$')
ax.set_ylabel(r'$N$')
ax.set_title(r'$Histograma$')
plt.bar(pos, Number, width, color='b')
plt.show()
