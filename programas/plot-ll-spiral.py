from __future__ import print_function
import json
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
import aplpy
from misc_utils import expand_fits_path
import seaborn as sns
sns.set_style('ticks')

groups = [ "LV", "SE", "N", "NW", "SW", "W", "S", ]
colors = [ "light brown", "light orange", "pale yellow",
           "dark pink", "purple", "cerulean", "apple green", ]
colors = sns.xkcd_palette(colors)
groups_and_colors = zip(groups, colors)

def find(name, path):
    """
    Original from http://stackoverflow.com/questions/1724693/find-a-file-in-python
    """
    for root, dirs, files in os.walk(path):
        for realname in name, "w" + name:
            if realname in files and not "_drz" in root:
                return os.path.join(root, realname)
    return None
            

def plot_map(limits, figname, canvas_size,
             fitsfile='$LARGE_FITS_DIR/wfi/Orion_H_A_shallow.fits',
             north=False,
             vmin=0.0, vmax=None, stretch='linear',
             innerbox=None, arrowscale=1.0):
    # Use an image as a backdrop
    fig = aplpy.FITSFigure(expand_fits_path(fitsfile),
                           figsize=canvas_size, north=north)
    # Set the viewport
    xc, yc = (limits[0] + limits[1])/2, (limits[2] + limits[3])/2
    w, h = limits[1] - limits[0], limits[3] - limits[2]
    fig.recenter(c0.ra.deg - xc/3600, c0.dec.deg + yc/3600,
                 width=w/3600, height=h/3600)
    fig.show_grayscale(vmin=vmin, vmax=vmax, invert=True,
                       stretch=stretch, interpolation='none')
    ax = fig._ax1


    c = coord.SkyCoord(RAs, Decs, unit=(u.hourangle, u.degree))
    # Cartesian pixel coordinates of each source
    x, y = fig.world2pixel(c.ra.deg, c.dec.deg)

    # Pixel size in degrees
    pix_scale = aplpy.wcs_util.celestial_pixel_scale(fig._wcs)
    # Convert to arcsec
    pix_scale *= 3600

    for group, color in groups_and_colors:
        mm = groups == group
        ax.plot(x[mm], y[mm], "o", mfc=color, mec='black', mew=1.0, alpha=0.8)


    c = coord.SkyCoord(pRAs, pDecs, unit=(u.hourangle, u.degree))
    x, y = fig.world2pixel(c.ra.deg, c.dec.deg)
    ax.scatter(x, y, c=pColors, s=pSizes, edgecolors='none', alpha=0.3, zorder=100)
    
    # Now plot the spiral shape
    # Radius in arcmin
    Rspiral = np.linspace(0.1, 20.0, 200)
    # Angle in degrees
    PAspiral = -180.0*np.log10(Rspiral)

    def polar2pixel(r, pa):
        ra = c0.ra.deg + r*np.sin(np.radians(pa))
        dec = c0.dec.deg + r*np.cos(np.radians(pa))
        return fig.world2pixel(ra, dec)
    
    x, y = polar2pixel(Rspiral/60.0, PAspiral)
    x1, y1 = polar2pixel(Rspiral/60.0, PAspiral-90.0)
    x2, y2 = polar2pixel(Rspiral/60.0, PAspiral+90.0)
    ax.plot(x, y, '-')
    ax.plot(x1, y1, '-')
    ax.plot(x2, y2, '-')
    
    fig.save(figname)


if __name__ == "__main__":

    #
    # Set up arc data
    #
    table = Table.read("luis-programas/arcs-summary-classify.tab", 
                     format="ascii.commented_header", delimiter="\t",
                     fill_values=('--', np.nan) ).filled(np.nan)
    names = table["Object"].data
    with open("luis-programas/problem-sources.txt") as f:
        problem_sources = f.read().split('\n')
    with open("luis-programas/interproplyd.txt") as f:
        interprop_sources = f.read().split('\n')
    problem_mask = np.array([name in problem_sources for name in table['Object']])
    interprop_mask = np.array([name in interprop_sources for name in table['Object']])
    m = (~problem_mask) & (~interprop_mask)
    RAs = table["RA"][m].data
    Decs = table["Dec"][m].data
    groups = table["Group"][m].data

    #
    # Set up proplyd data 
    #
    pColor_from_Type = {
        "i": "r", "d": "k", "rn": "c", "j": "g",
    }
    pSize_from_Type = {
        "i": 1.0, "d": 2.0, "rn": 0.7, "j": 0.7,
    }

    PROPLYD_TABLE_FILE = "ricci-data.json"
    ptable = json.load(open(PROPLYD_TABLE_FILE))

    pRAs = [v["RA"] for v in ptable.values()]
    pDecs = [v["Dec"] for v in ptable.values()]
    pColors = [pColor_from_Type[v["Type"]] for v in ptable.values()]
    pSizes = [pSize_from_Type[v["Type"]] for v in ptable.values()]
    pSizes = np.array(pSizes)*8.0


        
    c0 = coord.SkyCoord("05:35:16.463", "-05:23:23.18",
                        unit=(u.hourangle, u.degree))

    fullbox = [-425, 875, -825, 375]
    plot_map(fullbox, "ll-pos-spiral.pdf", (10, 10),
             vmin=5.0, vmax=2000.0, stretch='sqrt')










