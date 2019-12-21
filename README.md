# binospec
This is package for processing imaging data obtained with MMIRS and Binospec. It can be also used for other kind of images if they  have WCS somolar to MMIRS.

Installation process:

```python3 setup.py install```

There are two general functions in this package - astrometric correction (check_astrometry_wcs) and photometric calibration (calc_zp_ir). For both these functions script uses VO services to get reference catalogs, so Internet connection is required.