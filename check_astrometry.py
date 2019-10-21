import matplotlib
matplotlib.use('Agg')
import pyvo as vo
from astropy.io import fits
from astropy import wcs
import os
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
import ntpath
import argparse
from numpy.linalg import lstsq
from lmfit import Minimizer, Parameters, report_fit

parser = argparse.ArgumentParser(description='Photometric zeropoint calculation')
parser.add_argument("input_file", help="input fits file")
args = parser.parse_args()

filter_2mass_converter = {'K': 'k_m', 'J': 'j_m'}

filename = args.input_file #'/data1/Data/MMIRS/ComaA_K/Coma1_519_K_mean.fits'
fnm = ntpath.basename(filename).split('.fits')[0]
service = vo.dal.TAPService("http://vao.stsci.edu/PS1DR2/tapservice.aspx")
hdul = fits.open(filename)
frame = hdul[1]
w = wcs.WCS(frame.header)
center_coords = w.wcs_pix2world(frame.data.shape[0]/2., frame.data.shape[1]/2., 1)
query_text = "SELECT RAMean, DecMean, yKronMag  FROM dbo.StackObjectView WHERE CONTAINS(POINT('ICRS',RAMean, DecMean),CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND nDetections > 5 AND yKronMag < 24 AND yKronMag > -999" % (center_coords[0], center_coords[1], 0.15)
print query_text
#ir_cat = service.search(query_text)
#ir_cat = ir_cat.to_table()
#ir_cat.write('pannstarrs_cat.fits', format='fits', overwrite=True)
sex_command = "sex %s -CATALOG_NAME %s" % (filename, fnm+'_astrometry_scatalog.fits')
os.system(sex_command)
stilts_command = "java -jar ~/stilts.jar tskymatch2 ifmt1=fits ifmt2=fits find=best join=1and2 in1=%s#2 in2=pannstarrs_cat.fits ra1=ALPHA_J2000 dec1=DELTA_J2000 ra2=RAMean dec2=DecMean out=%s ofmt=fits error=4.5" % (fnm+'_astrometry_scatalog.fits', fnm+'_astrometry_matched.fits')
os.system(stilts_command)
final_tbl = Table.read(fnm+'_astrometry_matched.fits')
final_tbl = final_tbl[final_tbl['FLUX_RADIUS'] < 5.0]
final_tbl_ = final_tbl.copy()

params = Parameters()
params.add('CD1_1', value=1e-6)
params.add('CD1_2', value=1e-6)
params.add('CD2_1', value=1e-6)
params.add('CD2_2', value=-1e-6)
params.add('CRVAL1', value=195)
params.add('CRVAL2', value=27)


#xy_matrix = np.array([(final_tbl_['X_IMAGE'] - hdul[1].header['CRPIX1']).T, (final_tbl_['Y_IMAGE'] - hdul[1].header['CRPIX2']).T, np.ones(final_tbl_['DecMean'].T.shape)])
#rd_matrix = np.array([final_tbl_['RAMean'].T, final_tbl_['DecMean'].T])
#print xy_matrix.shape, rd_matrix.shape
#res = lstsq(xy_matrix.T, rd_matrix.T)
#print res
for i in range(5):
	xy_matrix = np.array([(final_tbl_['X_IMAGE'] - hdul[1].header['CRPIX2']).T, (final_tbl_['Y_IMAGE'] - hdul[1].header['CRPIX1']).T, np.ones(final_tbl_['DecMean'].T.shape)])
	rd_matrix = np.array([final_tbl_['RAMean'].T, final_tbl_['DecMean'].T])
	#print xy_matrix.shape, rd_matrix.shape
	res = lstsq(xy_matrix.T, rd_matrix.T)
	wcs_matrix = res[0]
	residuals = rd_matrix.T - np.dot(xy_matrix.T, wcs_matrix)
	print residuals.shape
	stds = np.std(residuals, axis=0)
	final_tbl_ = final_tbl_[((np.abs(residuals.T[0])<2*stds[0]) & (np.abs(residuals.T[1])<2*stds[1]))]
	print wcs_matrix
	ra_diff = final_tbl_['RAMean'] - final_tbl_['ALPHA_J2000']
	dec_diff = final_tbl_['DecMean'] - final_tbl_['DELTA_J2000']
	mean_ra_shift = np.mean(ra_diff)
	mean_ra_shift_std = np.std(ra_diff)
	mean_dec_shift = np.mean(dec_diff)
	mean_dec_shift_std = np.std(dec_diff)
	NUM = len(final_tbl_)
	print "Iter %i: d_RA=%.5f SIG_d_RA=%.5f d_DEC=%.5f SIG_d_DEC=%.5f %i" % (i, mean_ra_shift, mean_ra_shift_std, mean_dec_shift, mean_dec_shift_std, NUM)
	#final_tbl_ = final_tbl_[((np.abs(ra_diff - mean_ra_shift)<2*mean_ra_shift_std) & (np.abs(dec_diff - mean_dec_shift)<2*mean_dec_shift_std))]

print mean_ra_shift, mean_ra_shift_std, mean_dec_shift, mean_dec_shift_std
plt.subplot(121)
n, bins, patches = plt.hist(dec_diff, 10, density=True, facecolor='g', alpha=0.75)
plt.subplot(132)
n, bins, patches = plt.hist(ra_diff, 10, density=True, facecolor='g', alpha=0.75)
plt.savefig(fnm + '_shifts_hist.png')
plt.clf()
plt.quiver(final_tbl_['RAMean'], final_tbl_['DecMean'], ra_diff - mean_ra_shift, dec_diff - mean_dec_shift)
plt.savefig(fnm + '_shifts_vecfield.png')
hdr = hdul[1].header
hdr['CD1_1'] = wcs_matrix[0][0]
hdr['CD1_2'] = wcs_matrix[1][0]
hdr['CD2_1'] = wcs_matrix[0][1]
hdr['CD2_2'] = wcs_matrix[1][1]
hdr['CRVAL1'] = wcs_matrix[2][0]
hdr['CRVAL2'] = wcs_matrix[2][1]
print hdr
hdul[0].verify('fix')
hdul[1].verify('fix')
hdul.writeto(fnm + '_corrected.fits', overwrite=True)

