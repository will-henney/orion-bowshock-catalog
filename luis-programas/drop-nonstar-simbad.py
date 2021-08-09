import numpy as np
from astropy.table import Table, vstack, Column, MaskedColumn
import astropy.coordinates as coord

#Using ll-stars-arcs-simbad-047.ecsv on which was used a searched radi of 0.47 arcsec
tabsim = Table.read("ll-stars-arcs-simbad-053.ecsv", format="ascii")

#Dropping sources classified as ISM, unknown and HH on SIMBAD
mask = ((tabsim["OTYPE"] != "ISM") &
                (tabsim["OTYPE"] != "Unknown") &
                             (tabsim["OTYPE"] != "HH")) 

tabsim_result = tabsim[mask]
print("Number of sources:", len(tabsim_result))

# Read the table with the arcs
tab = Table.read('../luis-programas/arcs-summary-unique.ecsv', 
                      format='ascii.ecsv')
tab['coord'] = coord.SkyCoord(ra=tab['RA'], dec=tab['Dec'],
                                   unit=('hourangle', 'deg'))

cols = ['Object', 'RA', 'Dec', 'MAIN_ID', 'RA_star', 'DEC_star', 'SP_TYPE', 'OTYPE', 'Sep'] 

tabsim_result['coord'] = coord.SkyCoord(ra=tabsim_result['RA'], dec=tabsim_result['DEC'],
                                       unit=('hourangle', 'deg'))

id_, ra, dec, sp_type, otype, seps = [], [], [], [], [], []
for arcs in tab:
    sep = arcs['coord'].separation(tabsim_result['coord']).arcsec
    sepmin_i = sep.argmin()
    id_.append(tabsim_result['MAIN_ID'][sepmin_i])
    ra.append(tabsim_result['RA'][sepmin_i])
    dec.append(tabsim_result['DEC'][sepmin_i])
    sp_type.append(tabsim_result['SP_TYPE'][sepmin_i])
    otype.append(tabsim_result['OTYPE'][sepmin_i])
    seps.append(sep[sepmin_i])

m = np.array(seps) > 1.0
tab.add_column(
          MaskedColumn(name='MAIN_ID', data=id_, mask=m))
tab.add_column(
          MaskedColumn(name='RA_star', data=ra, mask=m))
tab.add_column(
          MaskedColumn(name='DEC_star', data=dec, mask=m))
tab.add_column(
          MaskedColumn(name='SP_TYPE', data=sp_type, mask=m))
tab.add_column(
          MaskedColumn(name='OTYPE', data=otype, mask=m))
tab.add_column(
          MaskedColumn(name='Sep', data=seps, mask=m, format='{:.3f}'))   
# Write the final table
tab[cols].write("ll-stars-arcs-simbad-clean.ecsv", format='ascii.ecsv', overwrite=True)
