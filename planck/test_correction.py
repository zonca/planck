from __future__ import division
import cPickle
import numpy as np
import healpy
import logging as l
import unittest
import matplotlib.pyplot as plt

from pointing import *
from correction import *
import private
from LFI import LFI

class TestCorrection(unittest.TestCase):

    def setUp(self):

        l.basicConfig(level=l.DEBUG,
                       format='%(asctime)s %(levelname)s %(message)s')

    def test_deaberration(self):
        """Test with Michele's corrections"""
        ch = LFI(instrument_db='/u/zonca/planck/data/mission/SIAM/LFI_instrumentDB_9.3.fits')['LFI28M']
        obt = np.array([1628860882.826]) # OD92
        pnt = Pointing(obt, coord='E', deaberration=False, wobble=False)
        vec = pnt.get(ch) # theta = 166 deg
        vecc = vec + deaberration(vec, obt, coord='E')
        qarray.norm_inplace(vecc)
        np.testing.assert_array_almost_equal(np.degrees(np.arccos(qarray.arraylist_dot(vec, vecc)))*60**2, np.array([[ 19.8740605]]))
        obt += 10. #30 sec after 14.6 deg
        pnt = Pointing(obt, coord='E', deaberration=False, wobble=False)
        vec = pnt.get(ch) # theta = 106 deg
        vecc = vec + deaberration(vec, obt, coord='E')
        qarray.norm_inplace(vecc)
        np.testing.assert_array_almost_equal(np.degrees(np.arccos(qarray.arraylist_dot(vec, vecc)))*60**2, np.array([[ 6.31009871]]))

    def test_get_wobble_psi2(self):
        self.assertAlmostEqual(get_wobble_psi2(1.067546522419189e+14/2**16), np.radians(-28.236426/60))

    def test_null_correction(self):
        r = wobble([0,0], wobble_psi2_model=lambda x:np.array([private.WOBBLE['psi2_ref']]*len(x))) 
        np.testing.assert_array_almost_equal(r, np.array([[-0., -0., -0., -1.], [-0., -0., -0., -1.]])) 

    def test_wobble_correction(self):
        r = wobble(np.array([1635615346.1697083])
        np.testing.assert_array_almost_equal(r, np.array([  2.08679712e-08,  -1.56235725e-05,  -0.00000000e+00, 1.00000000e+00]))


if __name__ == '__main__':
    unittest.main()
