'''
Estimated physical quantities
'''
import json
import numpy as np
from  astropy.table import Table


# Conversion counts/pixel -> erg/s/cm2/sr
Sfactor_ACS = 0.0025030687604156482

# Value of recombination coefficient in cm^3/s
alpha_B = 2.6e-13 
alpha_Ha = 1.27e-13

# Energy to 3 to 2... erg
Eha = 3.0267338723714944e-12 

# Relative extinction at Ha from Blagrave
fha = 0.78

# Distances 
D_orion_pc = 436.0
AU = 1.49597870691e13
cm_per_arcsec = D_orion_pc*AU

# Boltzmann in erg/K
k = 1.3806503e-16

# Temperatura
T = 1e4

# Stellar wind
yr = 3.15576e7
Msun = 1.989e33
km = 1.0e5

tab = Table.read("arcs-summary-merge.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )
tab.sort('RA')
#m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) & (tab['R_out'] > tab['R_in'])

D60 = tab['D']/60

with open("extinction.json") as f:
    extinction_data = json.load(f)

# Ha surface brightness, corrected for extinction
Sha = Sfactor_ACS*tab['Dif_Bally']
Chb = np.array([extinction_data.get(source, 0.0) for source in tab['Object']])
Sha *= 10**(fha*Chb)
# Correct for [N II] contamination of Ha filter
# This comes from the fit done in luis-programas/ratio-brightness.py
# Combined fit: Ratio = 0.28 D**0.43
Rnii_ha = 0.28*D60**0.43
Sha /= 1.0 + Rnii_ha

# Thickness and radius of the shell for measurements of delta l
h0 = tab['h']*cm_per_arcsec
rc = tab['Rc_out']*cm_per_arcsec
deltal = 2*np.sqrt(2*h0*rc)

nshell = np.sqrt(4.*np.pi*Sha/(alpha_Ha*deltal*Eha))
pshell = 2.0*nshell*k*T
MdotV_in = pshell*4.*np.pi*(tab['R_in']*cm_per_arcsec)**2*yr/Msun/km
MdotV_out = pshell*4.*np.pi*(60*D60*cm_per_arcsec)**2*yr/Msun/km

def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table
        
# Initialize a list for each column in output table
col_names = ["Object", "S(Halpha)", "Density", "Presure Shell", "Flow inner", "Flow outer"]
table = {cn: [] for cn in col_names}

for label, SHa, Nshell, Pshell, MV_in, MV_outer in zip(tab['Object'], Sha, nshell, pshell, MdotV_in, MdotV_out):
    table["Object"].append(label)
    if SHa!='--' and Nshell!='--' and Pshell!='--' and MV_in!='--' and MV_outer!='--':
        table["S(Halpha)"].append(str(SHa))
        table["Density"].append(str(Nshell))
        table["Presure Shell"].append(str(Pshell))
        table["Flow inner"].append(str(MV_in))
        table["Flow outer"].append(str(MV_outer))
    else:
        table["S(Halpha)"].append(str(SHa))
        table["Density"].append(str(Nshell))
        table["Presure Shell"].append(str(Pshell))
        table["Flow inner"].append(str(MV_in))
        table["Flow outer"].append(str(MV_outer))
#print table["S(Halpha)"]

with open("summary-physics.tab", "w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))

