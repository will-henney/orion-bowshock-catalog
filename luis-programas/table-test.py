'''
Automatically creates the table for latex
'''
import numpy as np
from  astropy.table import Table

arctab = Table.read("arcs-summary-merge.tab", format = "ascii")
Col = ['Object', 'RA', 'Dec', 'D','h','R_out','R_in','Rc_out','Rc_in']
arctab[Col].write('table-test-luis.tex', format = "ascii.latex") 
