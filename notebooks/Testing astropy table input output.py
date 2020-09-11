# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

from pathlib import Path
from astropy.table import Table
import numpy as np

datadir = Path("../luis-programas")

tab = Table.read(datadir / "arcs-summary-merge.ecsv")

tab.show_in_notebook()

t2 = Table.read(datadir / "arcs-summary-merge.ecsv", format="ascii.ecsv")

t2.show_in_notebook()

t2 == tab

t2[0] == tab[0]

t2[0]

tab[0]

# So it looks like `ascii.ecsv` is the best format for saving and loading tabular data.
#
# * The other ascii formats are too fragile, especially if there are missing values
# * FITS Table would work, but has the disadvantage that files are not human-readable, and also the string columns end up being `bytes`, not real strings
# * Also, it is explicitly called out in the astropy documentation that ECSV is preferred

tab = Table.read(datadir / "arcs-summary-classify.ecsv", format="ascii.ecsv")

tab.show_in_notebook()


