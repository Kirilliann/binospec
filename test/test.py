from __future__ import unicode_literals

import doctest
import unittest

import bino


class IsSnakeTestCase(unittest.TestCase):
    def test_quey_2M(self):
        self.assertEqual(bino.make_query(0,0,source='2M'), "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT(\'ICRS\',RAJ2000,DEJ2000), CIRCLE(\'ICRS\',0.00000,0.00000,0.30000))=1 AND Kmag>0 AND Jmag>0", msg='2MASS query test')
    def test_quey_UK(self):
        self.assertEqual(bino.make_query(0,0,source='2M'), "SELECT * FROM \"II/246/out\" WHERE CONTAINS(POINT(\'ICRS\',RAJ2000,DEJ2000), CIRCLE(\'ICRS\',0.00000,0.00000,0.30000))=1 AND Kmag>0 AND Jmag>0", msg='UKIDSS query test')
    def test_quey_GDR2(self):
        self.assertEqual(bino.make_query(0,0,source='GDR2'), "SELECT * FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS',0.00000,0.00000,0.30000))=1", msg='UKIDSS query test')
        
def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests

if __name__ == '__main__':
    unittest.main()