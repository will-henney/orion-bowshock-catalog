'''
Estimated physical quantities
'''
import json
import numpy as np
from  astropy.table import Table

def write_table(columns, col_names):
    """
    Write an ascii table of columns (sequence of sequences), using col_names as the header
    """
    table = "# " + "\t".join(col_names) + "\n"
    for row in zip(*columns):
        table += "\t".join(row) + "\n"
    return table
        
# Initialize a list for each column in output table
col_names = ["Object", "RA", "Dec", "S(Halpha)", "Density", "Pressure_Shell", "Flow_inner", "Flow_outer"]
table = {cn: [] for cn in col_names}

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

sources = ['066-652', '160-350', '162-456', '168-326N', '173-342', '204-330'] 
tab_merge = Table.read("arcs-summary-merge.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )
tab_pc = Table.read("arcs-summary-PC-will.tab", format="ascii.commented_header", delimiter="\t",
                 fill_values=('--', np.nan) )

mask = np.array([not a in sources for a in tab_merge['Object']])
m =  np.array([a in sources for a in tab_pc['Object']])
#m = np.isfinite(tab['R_out']) & np.isfinite(tab['R_in']) & (tab['R_out'] > tab['R_in'])

D60 = tab_merge['D'][mask]/60
D60_binary = tab_pc['D_binary'][m]/60

with open("extinction.json") as f:
    extinction_data = json.load(f)

# Ha surface brightness, corrected for extinction
Sha = Sfactor_ACS*tab_merge['Dif_Bally'][mask]
Chb = np.array([extinction_data.get(source, 0.0) for source in tab_merge['Object'][mask]])
Sha *= 10**(fha*Chb)
# Correct for [N II] contamination of Ha filter
# This comes from the fit done in luis-programas/ratio-brightness.py
# Combined fit: Ratio = 0.28 D**0.43
Rnii_ha = 0.28*D60**0.43
Sha /= 1.0 + Rnii_ha

# Thickness and radius of the shell for measurements of delta l
h0 = tab_merge['h'][mask]*cm_per_arcsec
rc = tab_merge['Rc_out'][mask]*cm_per_arcsec
deltal = 2*np.sqrt(2*h0*rc)

nshell = np.sqrt(4.*np.pi*Sha/(alpha_Ha*deltal*Eha))
pshell = 2.0*nshell*k*T
MdotV_in = pshell*4.*np.pi*(tab_merge['R_in'][mask]*cm_per_arcsec)**2*yr/Msun/km
MdotV_out = pshell*4.*np.pi*(60*D60*cm_per_arcsec)**2*yr/Msun/km

# from interproplyd

# Ha surface brightness, corrected for extinction
Sha_pc = Sfactor_ACS*np.array(tab_pc['Dif_Bally'][m], dtype = float)
Chb_pc = np.array([extinction_data.get(source, 0.0) for source in tab_pc['Object'][m]])
Sha_pc *= 10**(fha*Chb_pc)
# Correct for [N II] contamination of Ha filter
# This comes from the fit done in luis-programas/ratio-brightness.py
# Combined fit: Ratio = 0.28 D**0.43
Rnii_ha_pc = 0.28*D60_binary**0.43
Sha_pc /= 1.0 + Rnii_ha_pc

# Thickness and radius of the shell for measurements of delta l
h0_pc = tab_pc['h'][m]*cm_per_arcsec
rc_pc = tab_pc['Rc_out'][m]*cm_per_arcsec
deltal_pc = 2*np.sqrt(2*h0_pc*rc_pc)

nshell_pc = np.sqrt(4.*np.pi*Sha_pc/(alpha_Ha*deltal_pc*Eha))
pshell_pc = 2.0*nshell_pc*k*T
MdotV_in_pc = pshell_pc*4.*np.pi*(tab_pc['R_in'][m]*cm_per_arcsec)**2*yr/Msun/km
MdotV_out_pc = pshell_pc*4.*np.pi*(60*D60_binary*cm_per_arcsec)**2*yr/Msun/km

#for a, b in zip(tab_pc['Object'][m], nshell_pc):
   # print a, b

for label, ra, dec, SHa, Nshell, Pshell, MV_in, MV_outer in zip(tab_merge['Object'][mask], tab_merge['RA'][mask], tab_merge['Dec'][mask], Sha, nshell, pshell, MdotV_in, MdotV_out):
    table["Object"].append(label)
    table["RA"].append(ra)
    table["Dec"].append(dec)
    if SHa!='--' and Nshell!='--' and Pshell!='--' and MV_in!='--' and MV_outer!='--':
        table["S(Halpha)"].append("{:.5f}".format(SHa))
        table["Density"].append("{:.3f}".format(Nshell))  
        table["Pressure_Shell"].append(str(Pshell))
        table["Flow_inner"].append(str(Pshell))
        table["Flow_outer"].append(str(Pshell))
    else:
        table["S(Halpha)"].append('--')
        table["Density"].append('--')  
        table["Pressure_Shell"].append('--')
        table["Flow_inner"].append('--')
        table["Flow_outer"].append('--') 

for label_pc, ra_pc, dec_pc, SHa_pc, Nshell_pc, Pshell_pc, MV_in_pc, MV_outer_pc in zip(tab_pc['Object'][m], tab_pc['RA'][m], tab_pc['Dec'][m], Sha_pc, nshell_pc, pshell_pc, MdotV_in_pc, MdotV_out_pc):
    table["Object"].append(label_pc)
    table["RA"].append(ra_pc)
    table["Dec"].append(dec_pc)
    if SHa_pc!='--' and Nshell_pc!='-' and Pshell_pc!='-' and MV_in_pc!='-' and MV_outer_pc!='-':
        table["S(Halpha)"].append("{:.5f}".format(SHa_pc))
        table["Density"].append("{:.3f}".format(Nshell_pc))  
        table["Pressure_Shell"].append(str(Pshell_pc))
        table["Flow_inner"].append(str(Pshell_pc))
        table["Flow_outer"].append(str(Pshell_pc))
    else:
        table["S(Halpha)"].append('--')
        table["Density"].append('--')  
        table["Pressure_Shell"].append('--')
        table["Flow_inner"].append('--')
        table["Flow_outer"].append('--') 

with open("summary-physical.tab","w") as f:
    f.write(write_table([table[cn] for cn in col_names], col_names))
