from pathlib import Path
import sys
from astropy.io import fits
import numpy as np

field = sys.argv[1]
filter = sys.argv[2]

files = Path(".").glob(f"{field}-*-{filter}.fits")

imlist = []
for file in files:
    if "-best-" in str(file): continue
    if "-Bally-" in str(file): continue
    hdu = fits.open(file)[0]
    imlist.append(hdu.data)
imstack = np.stack(imlist)
imbest = np.nanmedian(imstack, axis=0)
fits.PrimaryHDU(
    header=hdu.header,
    data=imbest,
).writeto(f"{field}-best-{filter}.fits", overwrite=True)
