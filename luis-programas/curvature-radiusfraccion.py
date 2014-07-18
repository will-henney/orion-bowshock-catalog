import math
import matplotlib.pyplot as plt
import matplotlib.text as txt

infile = 'arcs-summary.tab' 

with open(infile,'r') as inputFile:
    for i,line in enumerate(inputFile):
        pass
nLines = i+1
#print nLines

d = []
r_out = []
rc_out = []
F = []
with open(infile,'r') as inputFile:
    for j in range(int(nLines)):
        line = inputFile.readline()
        linelist = line.split()
        D = linelist[3].strip()
        R_out = linelist[8].strip()
        Rc_out = linelist[10].strip()
        
        if (R_out != '-' and Rc_out !='-' ):
            d.append(D)
            r_out.append(R_out)
            rc_out.append(Rc_out)
dp=d[1:]
r_out_p = r_out[1:]
rc_out_p = rc_out[1:]
for d, r_out, rc_out in zip(dp, r_out_p, rc_out_p):
    d=float(d)
    r_out=float(r_out)
    rc_out=float(rc_out)

#cocientes de los radios
    fracc=float(rc_out)/r_out
    F.append(fracc)

print( dp, F)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
#ax1.set_xlim(xmin=100,xmax=600)
#ax1.set_ylim(ymin=0,ymax=10)
ax1.set_xlabel(r'$D$')
ax1.set_ylabel(r'$Rc_{in}/R_{in}$')
ax1.plot(dp,F,'bo')
ax1.set_title(r'$Rc_{in}/R_{in} vs D$')
ax1.grid(True)

plt.savefig('Rc-fraction.pdf')
