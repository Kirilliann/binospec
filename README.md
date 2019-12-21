# binospec
This is package for processing imaging data obtained with MMIRS and Binospec. It can be also used for other kind of images if they  have WCS somolar to MMIRS.

Installation process:

```python3 setup.py install```

There are two general functions in this package - astrometric correction (check_astrometry_wcs) and photometric calibration (calc_zp_ir). For both these functions script uses VO services to get reference catalogs, so Internet connection is required.

After installation process we recommend you to perform a suite of tests:

```python3 setup.py test```

To make sure that all function work well, please run them on test data:

in ipython3:

```
import bino

bino.corr_wcs('Coma1_519_K_mean.fits', sex_path='sextractor', SE_method='PU', sep=12.5, sc='GDR2')
bino.calc_zp_ir('Coma1_519_K_mean_corr_bgsubstr.fits', sc='UK', SE_method='PU')

```
