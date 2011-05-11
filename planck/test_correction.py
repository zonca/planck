from __future__ import division
import cPickle
import numpy as np
import healpy
import logging as l
import unittest
import matplotlib.pyplot as plt

from pointing import *
from correction import *
from LFI import LFI

class TestCorrection(unittest.TestCase):

    def setUp(self):

        l.basicConfig(level=l.DEBUG,
                       format='%(asctime)s %(levelname)s %(message)s')

    def test_deaberration(self):
        """Test with Michele's corrections"""
        obt = np.array([1628860882.826]) # OD92
        pnt = Pointing(obt, coord='E', deaberration=False, wobble=False)
        vec = pnt.get('LFI28M') # theta = 166 deg
        vecc = vec + deaberration(vec, obt, coord='E')
        qarray.norm_inplace(vecc)
        np.testing.assert_array_almost_equal(np.degrees(np.arccos(qarray.arraylist_dot(vec, vecc)))*60**2, np.array([[ 19.8740605]]))
        obt += 10. #30 sec after 14.6 deg
        pnt = Pointing(obt, coord='E', deaberration=False, wobble=False)
        vec = pnt.get('LFI28M') # theta = 106 deg
        vecc = vec + deaberration(vec, obt, coord='E')
        qarray.norm_inplace(vecc)
        np.testing.assert_array_almost_equal(np.degrees(np.arccos(qarray.arraylist_dot(vec, vecc)))*60**2, np.array([[ 6.31009871]]))

if __name__ == '__main__':
    unittest.main()
