"""
Find the WFPC2 field that corresponds to each source
"""
from __future__ import print_function
import math
import numpy as np
import json

file_name = "all-images.json"

with open(file_name) as f:
    data = json.load(f)

#ximage="$LARGE_FITS_DIR/wfpc2/hlsp_orion_hst_wfpc2_17_f656n_v1_sci.fits"

wfpc2_pattern = "wfpc2/hlsp_orion_hst_wfpc2"
filter_pattern = "f656n"

#b= data.keys()
source ="SRC=/fs/posgrado01/other0/angel/Dropbox/LuisBowshocks/programas"
print(source)
for obj, file_list in data.items():
    nfound = 0
    for filename in file_list:
        if wfpc2_pattern in filename and filter_pattern in filename:
            print("python", "$SRC/extract-image.py", obj, "--fitsfile", filename.split("/")[-1])
            nfound += 1
    
    if nfound == 0:
        print("#"+obj, "NO IMAGES FOUND")


      

