import matplotlib
matplotlib.use('Agg')
import pyvo as vo
from astropy.io import fits
from astropy import wcs
import os
import matplotlib.pyplot as plt
from astropy.table import Table, hstack
import numpy as np
import ntpath
import argparse
from numpy.linalg import lstsq
from lmfit import Minimizer, Parameters, report_fit
from astropy.coordinates import SkyCoord
from astropy import units as u

parser = argparse.ArgumentParser(description='WCS correction for images obtained with MMIRS. Kirill Grishin, email: kirillg6@gmail.com')
parser.add_argument("input_file", help="input fits file")
parser.add_argument("-d_ra", help="initial guess for shift in ra (deg)", default=0, type=float)
parser.add_argument("-d_dec", help="initial guess for shift in dec (deg)", default=0, type=float)
parser.add_argument("-bgsub", help="substract background (keyword)",nargs='?', const=True, default=True)
parser.add_argument("-mkplt", help="make plot",nargs='?', const=True, default=False)
parser.add_argument("-stilts_path", help="path to stilts .jar file", default='~/stilts.jar', type=str)
parser.add_argument("-sex_path", help="path to sextractor", default='sex', type=str)
parser.add_argument("-ext", help="extension", default=1, type=int)
parser.add_argument("-sep", help="minimun separation", default=12.5, type=float)
cd_pref = 'CD'
args = parser.parse_args()

min_sep = args.sep
filter_2mass_converter = {'K': 'k_m', 'J': 'j_m'}

filename = args.input_file #'/data1/Data/MMIRS/ComaA_K/Coma1_519_K_mean.fits'
fnm = ntpath.basename(filename).split('.fits')[0]
bg_framename = 'check.fits'
#service = vo.dal.TAPService("http://vao.stsci.edu/PS1DR2/tapservice.aspx")
service = vo.dal.TAPService("http://gea.esac.esa.int/tap-server/tap")
hdul = fits.open(filename)
frame = hdul[args.ext]
frame.header['CRVAl1'] = frame.header['CRVAl1'] + args.d_ra
frame.header['CRVAl2'] = frame.header['CRVAl2'] + args.d_dec
hdul[0].verify('fix')
hdul[args.ext].verify('fix')
tmp_filename = fnm + '_tmp.fits'
hdul.writeto(tmp_filename, overwrite=True)
w = wcs.WCS(frame.header)
center_coords = w.wcs_pix2world(frame.data.shape[0]/2., frame.data.shape[1]/2., 1)
query_text = "SELECT RAMean, DecMean, yKronMag  FROM dbo.StackObjectView WHERE CONTAINS(POINT('ICRS',RAMean, DecMean),CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND nDetections > 5 AND yKronMag < 24 AND yKronMag > -999" % (center_coords[0], center_coords[1], 0.15)
query_text = "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS',%.5f,%.5f,%.5f))=1" % (center_coords[0], center_coords[1], 0.15)
ir_cat = service.search(query_text)
ir_cat = ir_cat.to_table()
rf = open("test.reg", "a+")
rf.write("global color=green\nfk5\n")
for row in ir_cat:
	rf.write("circle(%.5f,%.5f,1.0\")\n" % (row['ra'], row['dec']))
rf.close()
#print ir_cat
#ir_cat.write('pannstarrs_cat.fits', format='fits', overwrite=True)
sex_command = "%s %s -CATALOG_NAME %s" % (args.sex_path, tmp_filename, fnm+'_astrometry_scatalog.fits')
os.system(sex_command)
stilts_command = "java -jar %s tskymatch2 ifmt1=fits ifmt2=fits find=best join=1and2 in1=%s#2 in2=pannstarrs_cat.fits ra1=ALPHA_J2000 dec1=DELTA_J2000 ra2=ra dec2=dec out=%s ofmt=fits error=14.5" % (args.stilts_path, fnm+'_astrometry_scatalog.fits', fnm+'_astrometry_matched.fits')
#os.system(stilts_command)
sx_cat = Table(fits.open(fnm+'_astrometry_scatalog.fits')[2].data)
sf = open("test_sex.reg", "a+")
sf.write("global color=red\n fk5\n")
xy_c = []
for row in sx_cat:
        sf.write("circle(%.5f,%.5f,1.0\")\n" % (row['ALPHA_J2000'], row['DELTA_J2000']))
	xy_c.append([row['X_IMAGE'], row['Y_IMAGE']])
sf.close()
sx_c = w.all_pix2world(xy_c,1)
c = SkyCoord(sx_c, unit=u.degree)
#c = SkyCoord(ra=sx_c[:][0]*u.degree,
#                dec=sx_c[:][1]*u.degree)
twomass_cat = SkyCoord(ra=ir_cat['ra'], 
                        dec=ir_cat['dec'])
idx, d2d, d3d = c.match_to_catalog_sky(twomass_cat)
inp_ind = np.arange(idx.shape[0])
idx_sld = idx[d2d.value<min_sep/3600.0]
inp_sld = inp_ind[d2d.value<min_sep/3600.0]
final_tbl = hstack([sx_cat[inp_sld], ir_cat[idx_sld]], join_type='inner')
print w.all_pix2world([[500,500]],1)
#final_tbl = Table.read(fnm+'_astrometry_matched.fits')
final_tbl = final_tbl[final_tbl['FLUX_RADIUS'] < 5.0]
final_tbl_ = final_tbl.copy()

ff = open("final_cat.reg", "a+")
ff.write("global color=blue\n fk5\n")
for row in final_tbl_:
        ff.write("circle(%.5f,%.5f,1.0\")\n" % (row['ra'], row['dec']))
ff.close()

params = Parameters()
params.add('CD1_1', value=w.pixel_scale_matrix[0,0])
params.add('CD1_2', value=w.pixel_scale_matrix[0,1])
params.add('CRVAL1', value=frame.header['CRVAL1'])
params.add('CRVAL2', value=frame.header['CRVAL2'])


#xy_matrix = np.array([(final_tbl_['X_IMAGE'] - hdul[1].header['CRPIX1']).T, (final_tbl_['Y_IMAGE'] - hdul[1].header['CRPIX2']).T, np.ones(final_tbl_['DecMean'].T.shape)])
#rd_matrix = np.array([final_tbl_['RAMean'].T, final_tbl_['DecMean'].T])
#print xy_matrix.shape, rd_matrix.shape
#res = lstsq(xy_matrix.T, rd_matrix.T)
#print res
hdr_c = frame.header.copy()

def fcn2min(params, res=False, diff=False):
	hdr_c[cd_pref+'1_1'] = params['CD1_1'].value
	hdr_c[cd_pref+'1_2'] = params['CD1_2'].value
	hdr_c[cd_pref+'2_1'] = params['CD1_2'].value
	hdr_c[cd_pref+'2_2'] = -params['CD1_1'].value
	hdr_c['CRVAL1'] = params['CRVAL1'].value
        hdr_c['CRVAL2'] = params['CRVAL2'].value
	w = wcs.WCS(hdr_c) 
	w2p = np.array([final_tbl_['ra'],final_tbl_['dec']])
	print w2p
	pix_c = w.all_world2pix(w2p.T,1)
	pix_fr = np.array([final_tbl_['X_IMAGE'],final_tbl_['Y_IMAGE']]).T
	resid = pix_c - pix_fr
	chi = np.concatenate((resid[:,0], resid[:,1]))
	if res:
		if diff:
			return resid
		else:
			return resid[:,0]**2/np.std(resid[:,0]) + resid[:,1]**2/np.std(resid[:,1])
	return chi

for i in range(5):
	minner = Minimizer(fcn2min, params)
	result = minner.minimize()
	chi2 = fcn2min(result.params, res=True)
	if i < 4:
		final_tbl_ = final_tbl_[chi2 < 6.0]
	print "Iter %i: CHI2=%.3f N_obj=%i" % (i, np.sum(chi2)/len(chi2), len(chi2))

#for i in range(5):
if args.mkplt:
	diffs = fcn2min(result.params, res=True, diff=True)
	plt.subplot(121)
	n, bins, patches = plt.hist(diffs[:,0], 10, density=True, facecolor='g', alpha=0.75)
	plt.xlabel('x offset, pix')
	plt.subplot(122)
	n, bins, patches = plt.hist(diffs[:,1], 10, density=True, facecolor='g', alpha=0.75)
	plt.xlabel('y offset, pix')
	plt.savefig(fnm + '_shifts_hist.png')
	plt.clf()
	q = plt.quiver(final_tbl_['X_IMAGE'], final_tbl_['Y_IMAGE'], diffs[:,0], diffs[:,1])
	plt.quiverkey(q, X=0.3, Y=-0.1, U=1, label='1 pix', labelpos='E')
	plt.xlim(0, frame.data.shape[1])
	plt.ylim(0, frame.data.shape[0])
	plt.gca().set_aspect('equal', adjustable='box')
	plt.title(fnm)
	plt.savefig(fnm + '_shifts_vecfield.png')
hdr = hdul[args.ext].header
hdr[cd_pref+'1_1'] = result.params['CD1_1'].value
hdr[cd_pref+'1_2'] = result.params['CD1_2'].value
hdr[cd_pref+'2_1'] = result.params['CD1_2'].value
hdr[cd_pref+'2_2'] = -result.params['CD1_1'].value
hdr['CRVAL1'] = result.params['CRVAL1'].value
hdr['CRVAL2'] = result.params['CRVAL2'].value
hdul[0].verify('fix')
hdul[args.ext].verify('fix')
hdul.writeto(fnm + '_corrected.fits', overwrite=True)
if args.bgsub:
	bg_img = fits.open(bg_framename)[0].data
	bg_level = np.median(bg_img)
	hdul[args.ext].data = hdul[args.ext].data - bg_level
hdul.writeto(fnm + '_corr_bgsubstr.fits', overwrite=True)



