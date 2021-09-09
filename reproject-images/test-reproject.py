INFILES=[["Robb-13", "5r"], ["Robb-13", "0l"], ["Bally-06", 1], ["Bally-06", 9]]
NEWFIELD="NW-Field"
from pathlib import Path
import numpy as np
from astropy.io import fits
from reproject import reproject_interp

base_path = Path("/Users/will/Work/OrionTreasury/")

def get_filename(dataset, field, filter):
    if dataset == "Robb-13":
        return base_path / "acs" / f"hlsp_orion_hst_acs_strip{field}_{filter}_v1_drz.fits"
    elif dataset == "Bally-06":
        if filter == "f658n":
            return base_path / "Bally-ACS" / f"j8oc{field:02d}010_wcs.fits"
        else:
            return None
    else:
        return None

hdu0 = fits.open(f"blank-{NEWFIELD}.fits")[0]
FILTERS = "f658n", "f555w", "f435w"
for dataset, field in INFILES:
    for filter in FILTERS:
        fn = get_filename(dataset, field, filter)
        if fn is not None:
            hdulist = fits.open(fn)
            hdu = hdulist["SCI"]
            # Remove any unwanted SIP keywords
            del hdu.header["A_*"]
            del hdu.header["B_*"]
            # switch zeros for NaN
            hdu.data = np.where(
                hdu.data > 0.0,
                hdu.data,
                np.nan,
            )
            newdata, _ = reproject_interp(hdu, hdu0.header, order="nearest-neighbor")
            newhdu = fits.PrimaryHDU(header=hdu0.header, data=newdata)
            try:
                field = f"{field:02d}"
            except:
                pass
            newfn = f"{NEWFIELD}-{dataset}-{field}-{filter}.fits"
            newhdu.writeto(newfn, overwrite=True)
        else:
            pass
