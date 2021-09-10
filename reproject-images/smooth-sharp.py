import click
import numpy as np
from astropy.io import fits
from astropy.convolution import convolve_fft, Gaussian2DKernel

@click.command()
@click.argument("infile")
@click.option("--width", default=20, show_default=True,
              help="Width of Gaussian smoothing kernel in pixels")
@click.option("--suffix", default="",
              help="Extra id string to add to output filename")
@click.option("--twopass/--no-twopass", default=False, show_default=True,
              help="Whether to do a second pass for eliminating negative halos")
@click.option("--threshold", default=1.5, show_default=True,
              help="Mask out relative brightnesses higher than this in second pass")
def main(infile, width, suffix, twopass, threshold):
    """
    Do high-pass filtering of Fits file INFILE
    """
    hdu = fits.open(infile)[0]

    im = hdu.data

    gauss_kernel = Gaussian2DKernel(width)
    smoothim = convolve_fft(im, gauss_kernel, allow_huge=True)
    sharpim = im/smoothim

    if twopass:
        mask = (sharpim > threshold) | (im == 0.0)
        im[mask] = np.nan
        print('Starting second pass: N(masked) =', mask.sum())
        smoothim = convolve_fft(im, gauss_kernel, normalize_kernel=True, allow_huge=True)
        sharpim = im/smoothim

    outhdu = fits.PrimaryHDU(data=smoothim, header=hdu.header)

    outfile = infile.replace(".fits", f"_smooth_{width}{suffix}.fits")
    outhdu.writeto(outfile, overwrite=True)

    outhdu.data = sharpim
    outfile = outfile.replace("_smooth_", "_sharp_")
    outhdu.writeto(outfile, overwrite=True)

if __name__ == "__main__":
    main()
