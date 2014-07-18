"""
Find the extinction at the position of each source

Read data from C(Hbeta) map 
"""

from astropy.io import fits
import json
import numpy as np
import pyregion

hdu = fits.open("chbeta-radec.fits")[0]

regions = pyregion.open("ll-boxes-new.reg")
results = {}
for region in regions:
    # Extract name of source
    source = region.comment.split('{')[-1].split('}')[0]
    
    # Find C(Hb) in source box

    # Make a ShapeList that contains only this source
    shape_list =  pyregion.ShapeList([region])
    mask = shape_list.get_mask(hdu)
    chb = hdu.data[mask].mean()
  
    # Check that the C(Hb) value is sensible (finite and non-negative)
   
    if np.isfinite(chb) and chb > 0.0:
        results[source] = float(chb)
    else:
        results[source] = 0.0

# Save results in JSON file
jsonfile = "extinction.json"
with open(jsonfile, "w") as f:
    json.dump(results, f, indent=4)


