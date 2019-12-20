import matplotlib
matplotlib.use('Agg')
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
from .utils import *


def calc_zp_ir(input_file, sc='2M', sep=0.8, log=False,
               logfile='LOG_ZP.csv', savecat=False, sex_path='sex',
               savesc=True, plot=True, SE_method='PU'):
    """
    Zeropoint correction for MMIRS data. Kirill Grishin, email: grishin@voxastro.org
    Function for zeropoint calculation
    input_file - input fits file
    sc - source catalog, available options are 2M - 2MASS and UK - UKIDSS, default - 2M
    logfile - full path to log file (if file doesn't exists, it will be created), default - LOG_ZP.csv
    sex_path - path to sexstractor (or full command to execute it), default - sex
    savecat - save input catalog< optained with VO services, default - False
    savesc - save SExtractor output, default - yes
    plot - make a photometric plot
    """
    data_source = sc
    min_sep = sep
    if sc not in ['2M', 'UK']:
        raise Exception('Not allowed source catalog used. Please use 2MASS (2M keyword) or UKIDSS (UK keyword)')
    filter_2mass_converter = {'K': {'2M':'Kmag', 'UK':'Kmag'},
                                'J': {'2M':'Jmag', 'UK':'Jmag1'}}
    s_design = {'2M':'2MASS',
            'UK':'UKIDSS DR9'}
    filename = input_file
    fnm = ntpath.basename(filename).split('.fits')[0]
    hdul = fits.open(filename)
    frame = hdul[1]
    w = wcs.WCS(frame.header)
    center_coords = w.wcs_pix2world(frame.data.shape[0]/2., frame.data.shape[1]/2., 1)
    query_text = make_query(center_coords[0], center_coords[1], source=sc, r=0.2)
    ir_cat = get_table(query_text, source=sc, savecat=savecat, catnm=fnm + '_stellar_cat.fits')
    source_extraction(filename, fnm+'_zpcalc_scatalog.fits', method=SE_method, sex_path=sex_path)
    sctbl_ext = 2 if SE_method=='SE' else 1
    sx_cat = Table(fits.open(fnm+'_zpcalc_scatalog.fits')[sctbl_ext].data)
    final_tbl = cross_match(sx_cat, ir_cat)
    if SE_method == 'SE':
        final_tbl = final_tbl[final_tbl['FLUX_RADIUS'] < 5.0]
    y_plt = -2.5*np.log10(final_tbl['FLUX_AUTO'])
    x_plt = final_tbl[filter_2mass_converter[frame.header['FILTID1']][data_source]]
    zp = x_plt - y_plt
    zp_c = np.array(zp.copy())
    for i in range(10):
            zp_mean = np.nanmean(zp_c)
            zp_std = np.nanstd(zp_c)
            NUM = len(zp_c)
            zp_c = zp_c[np.abs(zp_c - zp_mean)<2*zp_std]
            print("Iter %i: ZP=%.3f SIG=%.3f NUM=%i" % (i, zp_mean, zp_std, NUM))
    if not savesc:
            os.system("rm %s" % (fnm+'_zpcalc_scatalog.fits'))
    if plot:
            plt.clf()
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
            plt.clf()
    if log:
            f = open(logfile,"a+")
            f.write("%s, %s, %s, %.3f, %.3f, %i\r\n" % (fnm+'.fits',frame.header['FILTID1'],s_design[data_source],zp_mean, zp_std, NUM))
            f.close()