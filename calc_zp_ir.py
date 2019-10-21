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
from astropy.coordinates import SkyCoord
from astropy import units as u

parser = argparse.ArgumentParser(description='Photometric zeropoint calculation for images obtained with MMIRS. Kirill Grishin, email: kirillg6@gmail.com')
parser.add_argument("input_file", help="input fits file")
parser.add_argument("-sc", help="source catalog (default 2MASS), options(keys): 2MASS(2M), UKIDSS(UK) ", default='2M', choices=['2M', 'UK'])
parser.add_argument("-sep", help="minimum separation (in arcsec) default 0.8 arcsec", default=0.8, type=float)
parser.add_argument("-log", help="put zerpoint into log file (keyword)",nargs='?', const=True, default=False)
parser.add_argument("-logfile", help="name of log-file", default='LOG_ZP.csv', type=str)
parser.add_argument("-savecat", help="save input catalog",nargs='?', const=True, default=False)
parser.add_argument("-sex_path", help="path to sextractor", default='sex', type=str)
parser.add_argument("-savesc", help="save Sextractor catalog (keyword)",nargs='?', const=True, default=False)
parser.add_argument("-plot", help="make plots (keyword)",nargs='?', const=True, default=False)
args = parser.parse_args()

data_source = args.sc
min_sep = args.sep
filter_2mass_converter = {'K': {'2M':'Kmag', 'UK':'Kmag'},
			'J': {'2M':'Jmag', 'UK':'Jmag1'}}
s_design = {'2M':'2MASS',
	'UK':'UKIDSS DR9'}
filename = args.input_file #'/data1/Data/MMIRS/ComaA_K/Coma1_519_K_mean.fits'
fnm = ntpath.basename(filename).split('.fits')[0]
service = vo.dal.TAPService("http://tapvizier.u-strasbg.fr/TAPVizieR/tap")
hdul = fits.open(filename)
frame = hdul[1]
w = wcs.WCS(frame.header)
center_coords = w.wcs_pix2world(frame.data.shape[0]/2., frame.data.shape[1]/2., 1)
query_2mass = "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT('ICRS',RAJ2000,DEJ2000), CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND Kmag>0 AND Jmag>0" % (center_coords[0], center_coords[1], 0.3)
query_ukidss = "SELECT * FROM \"II/319/las9\" WHERE CONTAINS(POINT('ICRS',RAJ2000,DEJ2000), CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND Kmag>0 AND Jmag1>0" % (center_coords[0], center_coords[1], 0.3)
query_text = query_2mass if data_source == '2M' else query_ukidss
ir_cat = service.search(query_text)
ir_cat = ir_cat.to_table()
for c_name in ir_cat.colnames:
	if ir_cat[c_name].dtype == 'O':
		ir_cat.remove_column(c_name)
if args.savecat:
	catnm = fnm + '_stellar_cat.fits'
	ir_cat.write(catnm, format='fits', overwrite=True)
sex_command = "%s %s -CATALOG_NAME %s -CHECKIMAGE_TYPE NONE" % (args.sex_path, filename, fnm+'_zpcalc_scatalog.fits')
os.system(sex_command)
sx_cat = Table(fits.open(fnm+'_zpcalc_scatalog.fits')[2].data)
c = SkyCoord(ra=sx_cat['ALPHA_J2000']*u.degree,
		dec=sx_cat['DELTA_J2000']*u.degree)
twomass_cat = SkyCoord(ra=ir_cat['RAJ2000'], 
			dec=ir_cat['DEJ2000'])
idx, d2d, d3d = c.match_to_catalog_sky(twomass_cat)
inp_ind = np.arange(idx.shape[0])
idx_sld = idx[d2d.value<min_sep/3600.0]
inp_sld = inp_ind[d2d.value<min_sep/3600.0]
final_tbl = hstack([sx_cat[inp_sld], ir_cat[idx_sld]], join_type='inner')
final_tbl = final_tbl[final_tbl['FLUX_RADIUS'] < 5.0]
y_plt = -2.5*np.log10(final_tbl['FLUX_AUTO'])
print filter_2mass_converter[frame.header['FILTID1']][data_source]
x_plt = final_tbl[filter_2mass_converter[frame.header['FILTID1']][data_source]]

zp = x_plt - y_plt
zp_c = np.array(zp.copy())
for i in range(10):
	zp_mean = np.nanmean(zp_c)
	zp_std = np.nanstd(zp_c)
	NUM = len(zp_c)
	zp_c = zp_c[np.abs(zp_c - zp_mean)<2*zp_std]
	print "Iter %i: ZP=%.3f SIG=%.3f NUM=%i" % (i, zp_mean, zp_std, NUM)
if not args.savesc:
	os.system("rm %s" % (fnm+'_zpcalc_scatalog.fits'))
if args.plot:
	plt.plot(x_plt, y_plt, 'k+')
	line_x = np.array([min(x_plt), max(x_plt)])
	plt.plot(line_x, line_x - zp_mean, 'b')
	plt.xlabel(s_design[data_source]+ ' '+ filter_2mass_converter[frame.header['FILTID1']][data_source])
	plt.ylabel('-2.5lg(F)')
	plt.title(fnm)
	plt.savefig(fnm+'_zpcalc_plot.png')
	plt.clf()
	y_plt_c = -x_plt + y_plt
	x_plt_c = final_tbl['Kmag'] - final_tbl[filter_2mass_converter['J'][data_source]]
	plt.plot(x_plt_c, y_plt_c, 'b+')
	plt.savefig(fnm+'_zpcalc_plot_color.png')

if args.log:
	f = open(args.logfile,"a+")
	f.write("%s, %s, %s, %.3f, %.3f, %i\r\n" % (fnm+'.fits',frame.header['FILTID1'],s_design[data_source],zp_mean, zp_std, NUM))
	f.close()


