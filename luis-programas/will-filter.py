from __future__ import print_function
import pysynphot
from matplotlib import pyplot as plt

AIR_REFRACTIVE_INDEX = 1.000277 # @STP according to Wikipedia
VEL0 = +25.0
LIGHTSPEED = 2.99792458e5

doppler = (1.0 + VEL0/LIGHTSPEED)

for fn in 'wfpc2,F656N', 'wfpc2,F658N', 'wfpc2,F547M', 'acs,hrc,f658n':
    bandpass = pysynphot.ObsBandpass(fn)
    print(fn, 'Rectangular width =', bandpass.rectwidth())
    plt.plot(bandpass.wave/AIR_REFRACTIVE_INDEX, bandpass.throughput, label=fn)

plt.xlim(6450, 6700)
for wav0, height in (6548.05, 0.33), (6562.79, 2), (6583.45, 1):
    plt.axvline(wav0, 0.0, 0.8*height, color='r', lw=2)
plt.legend(fontsize='small')
plt.xlabel("Air wavelength, Angstrom")
plt.ylabel("Filter throughput")

plt.savefig('will-filter.pdf')
