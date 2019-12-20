from __future__ import unicode_literals

import doctest
import unittest

import bino


class voQueryTestCase(unittest.TestCase):
    def test_query_2M(self):
        self.assertEqual(bino.make_query(0,0,source='2M'), "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT(\'ICRS\',RAJ2000,DEJ2000), CIRCLE(\'ICRS\',0.00000,0.00000,0.30000))=1 AND Kmag>0 AND Jmag>0", msg='2MASS query text test')
    def test_query_UK(self):
        self.assertEqual(bino.make_query(0,0,source='2M'), "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT(\'ICRS\',RAJ2000,DEJ2000), CIRCLE(\'ICRS\',0.00000,0.00000,0.30000))=1 AND Kmag>0 AND Jmag>0", msg='UKIDSS query text test')
    def test_query_GDR2(self):
        self.assertEqual(bino.make_query(0,0,source='GDR2'), "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS',0.00000,0.00000,0.30000))=1", msg='GAIA query text test')
    def test_query_PSDR2(self):
        self.assertEqual(bino.make_query(0,0,source='PSDR2'), "SELECT RAMean, DecMean, yKronMag  FROM dbo.StackObjectView WHERE CONTAINS(POINT('ICRS',RAMean, DecMean),CIRCLE('ICRS',0.00000,0.00000,0.30000))=1 AND nDetections > 5 AND yKronMag < 24 AND yKronMag > -999", msg='UKIDSS query test')
    def test_VO_2M(self):
        t = bino.get_table(bino.make_query(0,0,source='2M'), source='2M')
        self.assertEqual(len(t), 462)
    def test_VO_UK(self):
        t = bino.get_table(bino.make_query(0,0,source='UK'), source='UK')
        self.assertEqual(len(t), 3417)
    def test_VO_GDR2(self):
        t = bino.get_table(bino.make_query(0,0,source='GDR2'), source='GDR2')
        self.assertEqual(len(t), 858)
    def test_VO_PSDR2(self):
        t = bino.get_table(bino.make_query(0,0,source='PSDR2',r=0.10), source='PSDR2')
        self.assertEqual(len(t), 318)
        
def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests

if __name__ == '__main__':
    unittest.main()