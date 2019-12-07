from astropy.io import fits
import numpy as np
import argparse
import ntpath

parser = argparse.ArgumentParser(description='WCS for error images')
parser.add_argument("input_file", help="input fits file")
parser.add_argument("err_file", help="error image")
args = parser.parse_args()

filename = args.err_file
fnm = ntpath.basename(filename).split('.fits')[0]
img_frame = fits.open(args.input_file)
hdr_inp = img_frame[1].header
err_hdul = fits.open(filename)
hdr = err_hdul[1].header

hdr['CD1_1'] = hdr_inp['CD1_1']
hdr['CD1_2'] = hdr_inp['CD1_2']
hdr['CD2_1'] = hdr_inp['CD2_1']
hdr['CD2_2'] = hdr_inp['CD2_2']
hdr['CRVAL1'] = hdr_inp['CRVAL1']
hdr['CRVAL2'] = hdr_inp['CRVAL2']

err_hdul[0].verify('fix')
err_hdul[1].verify('fix')

err_hdul.writeto(fnm+'_corrected.fits', overwrite=True)
