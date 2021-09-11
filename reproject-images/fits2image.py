from matplotlib import pyplot as plt
from astropy.io import fits
import click

@click.command()
@click.argument("fitsname")
@click.option("--vmin", default=0.5)
@click.option("--vmax", default=1.5)
@click.option("--format", default="png")
@click.option("--verbose/--no-verbose", default=False)
def fits2image(fitsname, vmin, vmax, format, verbose):
    """
    Convert FITSNAME.fits to an image
    """
    hdu = fits.open(fitsname)[0]
    imname = fitsname.replace(".fits", f".{format}")
    plt.imsave(imname, hdu.data, vmin=vmin, vmax=vmax, cmap="gray", origin="lower")
    if verbose:
        print("Saved", imname)

if __name__ == "__main__":
    fits2image()
