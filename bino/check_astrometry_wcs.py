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
from .utils import *


def corr_wcs(input_file, d_ra=0, d_dec=0, bgsub=True,
             mkplt=True, stilts_path='~/stilts.jar', sex_path='sex',
             ext=1, sep=12.5, SE_method='PU', cd_pref = 'CD', 
             sc='GDR2', savecat=False):
    """
    WCS correction for images obtained with MMIRS. Kirill Grishin, email: grishin@voxastro.org
    input_file -- input fits file with images
    d_ra -- initial guess in ra shift (deg)
    d_dec -- initial guess in dec shift (deg)
    bgsub -- substract background
    mkplt -- make plots
    stilts_path -- path to stilts .jar file
    sex_path -- path to sextractor
    ext -- extension
    sep -- minimun separation
    SE_method - method for source extraction. Options: SE - SExtractor, PU (default) - PhotUtils
    """
    min_sep = sep
    filter_2mass_converter = {'K': 'k_m', 'J': 'j_m'}

    filename = input_file
    fnm = ntpath.basename(filename).split('.fits')[0]
    bg_framename = 'check.fits'
    hdul = fits.open(filename)
    frame = hdul[ext]
    frame.header['CRVAl1'] = frame.header['CRVAl1'] + d_ra
    frame.header['CRVAl2'] = frame.header['CRVAl2'] + d_dec
    hdul[0].verify('fix')
    hdul[ext].verify('fix')
    tmp_filename = fnm + '_tmp.fits'
    hdul.writeto(tmp_filename, overwrite=True)
    w = wcs.WCS(frame.header)
    center_coords = w.wcs_pix2world(frame.data.shape[0]/2., frame.data.shape[1]/2., 1)
    query_text = make_query(center_coords[0], center_coords[1], source=sc, r=0.2)
    ir_cat = get_table(query_text, source=sc, savecat=savecat, catnm=fnm + '_stellar_cat.fits')
    rf = open("test.reg", "a+")
    rf.write("global color=green\nfk5\n")
    for row in ir_cat:
            rf.write("circle(%.5f,%.5f,1.0\")\n" % (row['RAJ2000'], row['DEJ2000']))
    rf.close()
    source_extraction(tmp_filename, fnm+'_astrometry_scatalog.fits', method=SE_method, sex_path=sex_path)
    sctbl_ext = 2 if SE_method=='SE' else 1
    sx_cat = Table(fits.open(fnm+'_astrometry_scatalog.fits')[sctbl_ext].data)
    sf = open("test_sex.reg", "a+")
    sf.write("global color=red\n fk5\n")
    xy_c = []
    for row in sx_cat:
            sf.write("circle(%.5f,%.5f,1.0\")\n" % (row['ALPHA_J2000'], row['DELTA_J2000']))
            xy_c.append([row['X_IMAGE'], row['Y_IMAGE']])
    sf.close()
    sx_c = w.all_pix2world(xy_c,1)
    final_tbl = cross_match(sx_cat, ir_cat, b_ra_name='RAJ2000', b_dec_name='DEJ2000')
    if SE_method == 'SE':
        final_tbl = final_tbl[final_tbl['FLUX_RADIUS'] < 5.0]
    final_tbl_ = final_tbl.copy()

    ff = open("final_cat.reg", "a+")
    ff.write("global color=blue\n fk5\n")
    for row in final_tbl_:
            ff.write("circle(%.5f,%.5f,1.0\")\n" % (row['RAJ2000'], row['DEJ2000']))
    ff.close()

    params = Parameters()
    params.add('CD1_1', value=w.pixel_scale_matrix[0,0])
    params.add('CD1_2', value=w.pixel_scale_matrix[0,1])
    params.add('CRVAL1', value=frame.header['CRVAL1'])
    params.add('CRVAL2', value=frame.header['CRVAL2'])

    hdr_c = frame.header.copy()

    def fcn2min(params, res=False, diff=False):
            hdr_c[cd_pref+'1_1'] = params['CD1_1'].value
            hdr_c[cd_pref+'1_2'] = params['CD1_2'].value
            hdr_c[cd_pref+'2_1'] = params['CD1_2'].value
            hdr_c[cd_pref+'2_2'] = -params['CD1_1'].value
            hdr_c['CRVAL1'] = params['CRVAL1'].value
            hdr_c['CRVAL2'] = params['CRVAL2'].value
            w = wcs.WCS(hdr_c) 
            w2p = np.array([final_tbl_['RAJ2000'],final_tbl_['DEJ2000']])
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
            print("Iter %i: CHI2=%.3f N_obj=%i" % (i, np.sum(chi2)/len(chi2), len(chi2)))

    if mkplt:
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
    hdr = hdul[ext].header
    hdr[cd_pref+'1_1'] = result.params['CD1_1'].value
    hdr[cd_pref+'1_2'] = result.params['CD1_2'].value
    hdr[cd_pref+'2_1'] = result.params['CD1_2'].value
    hdr[cd_pref+'2_2'] = -result.params['CD1_1'].value
    hdr['CRVAL1'] = result.params['CRVAL1'].value
    hdr['CRVAL2'] = result.params['CRVAL2'].value
    hdul[0].verify('fix')
    hdul[ext].verify('fix')
    hdul.writeto(fnm + '_corrected.fits', overwrite=True)
    if bgsub:
            bg_img = fits.open(bg_framename)[0].data
            bg_level = np.median(bg_img)
            hdul[ext].data = hdul[ext].data - bg_level
    hdul.writeto(fnm + '_corr_bgsubstr.fits', overwrite=True)