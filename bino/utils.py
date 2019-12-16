from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np

def cross_match(a_cat, b_cat, a_ra_name='ALPHA_J2000',
                b_ra_name='RAJ2000', a_dec_name='DELTA_J2000',
                b_dec_name='DEJ2000', sep=5.0):
    c = SkyCoord(ra=a_cat[a_ra_name]*u.degree,
                    dec=a_cat[a_dec_name]*u.degree)
    twomass_cat = SkyCoord(ra=b_cat[b_ra_name], 
                            dec=b_cat[b_dec_name])
    idx, d2d, d3d = c.match_to_catalog_sky(twomass_cat)
    inp_ind = np.arange(idx.shape[0])
    idx_sld = idx[d2d.value<sep/3600.0]
    inp_sld = inp_ind[d2d.value<sep/3600.0]
    final_tbl = hstack([a_cat[inp_sld], b_cat[idx_sld]], join_type='inner')
    return final_tbl

def make_query(ra, dec, source='2M', r=0.3):
    if source == '2M':
        s = "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT('ICRS',RAJ2000,DEJ2000), CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND Kmag>0 AND Jmag>0" % (ra, dec, r)
    if source == 'UK':
        s = "SELECT * FROM \"II/319/las9\" WHERE CONTAINS(POINT('ICRS',RAJ2000,DEJ2000), CIRCLE('ICRS',%.5f,%.5f,%.5f))=1 AND Kmag>0 AND Jmag1>0" % (ra, dec, r)
    return s