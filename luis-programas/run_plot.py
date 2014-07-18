"""
Find the argument for each source to run "plot-image.py"
"""
from __future__ import print_function
import glob
import json

pattern = "j8oc??010_wcs/*-arcdata.json"

file_list = glob.glob(pattern)

source ="SRC=/fs/posgrado01/other0/angel/Dropbox/LuisBowshocks/programas"
print(source)
for file_name in file_list:
    with open(file_name) as f:
        data = json.load(f)

    for k in data.keys():
        if k.startswith('Robberto_WFPC2'):
            imagename=k
            print("python", "$SRC/plot-image.py", data["star"]["id"], "--image", imagename)    
           
     
