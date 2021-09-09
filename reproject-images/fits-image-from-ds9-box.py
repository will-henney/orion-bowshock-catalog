import sys
import numpy as np
from astropy.wcs import WCS
import astropy.units as u
from astropy.io import fits
import regions


def get_fits_image(box, pix_arcsec=0.05):
    w = WCS(naxis=2)
    w.wcs.crval = box.center.ra.deg, box.center.dec.deg
    w.wcs.ctype = "RA---TAN", "DEC--TAN"
    w.wcs.cdelt = -pix_arcsec / 3600, pix_arcsec / 3600
    bb = box.to_pixel(w).bounding_box
    w.wcs.crpix = -bb.ixmin, -bb.iymin
    bb = box.to_pixel(w).bounding_box
    nx, ny = bb.ixmax + 1, bb.iymax + 1
    return fits.PrimaryHDU(header=w.to_header(), data=np.ones((ny, nx)))


if __name__ == "__main__":
    regfile = sys.argv[1]
    regionlist = regions.read_ds9(regfile)
    for region in regionlist:
        if type(region) == regions.shapes.rectangle.RectangleSkyRegion:
            hdu = get_fits_image(region)
            fname = "-".join(["blank"] + region.meta["label"].split()) + ".fits"
            hdu.writeto(fname, overwrite=True)
