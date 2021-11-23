from pathlib import Path
import sys
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import regions


field = sys.argv[1]
filter1 = sys.argv[2]
filter2 = sys.argv[3]

files1 = list(Path(".").glob(f"{field}-*-{filter1}.fits"))
files2 = list(Path(".").glob(f"{field}-*-{filter2}.fits"))

imlist = []
for file in files1 + files2:
    if "-best-" in str(file):
        continue
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
    if "f435w" in str(file):
        # Don't let f435w be the minimum unless we absolutely have to
        hdu.data *= 2
    imlist.append(hdu.data)
imstack = np.stack(imlist)
imbest = np.nanmin(imstack, axis=0)
fits.PrimaryHDU(
    header=hdu.header,
    data=imbest,
).writeto(f"{field}-best-{filter1}-and-{filter2}.fits", overwrite=True)
