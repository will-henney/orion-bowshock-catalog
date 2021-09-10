from pathlib import Path
import sys
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import regions


field = sys.argv[1]
filter = sys.argv[2]

files = Path(".").glob(f"{field}-*-{filter}.fits")

imlist = []
for file in files:
    if "-best-" in str(file): continue
    hdu = fits.open(file)[0]
    regfile = str(file).replace(".fits", "-exclude.reg")
    try:
        regs = regions.read_ds9(regfile)
        print("Masking:", regs)
        w = WCS(hdu.header)
        shape = hdu.data.shape
        for reg in regs:
            m = reg.to_pixel(w).to_mask().to_image(shape).astype(bool)
            hdu.data[m] = np.nan
    except:
        pass
    imlist.append(hdu.data)
imstack = np.stack(imlist)
imbest = np.nanmin(imstack, axis=0)
fits.PrimaryHDU(
    header=hdu.header,
    data=imbest,
).writeto(f"{field}-best-{filter}.fits", overwrite=True)
