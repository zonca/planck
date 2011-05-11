from __future__ import division
import cPickle
import numpy as np
import healpy
import logging as l
import unittest
import matplotlib.pyplot as plt

from pointing import *
from LFI import LFI
from testenv.remix import read_exchange

class TestPointing(unittest.TestCase):

    def setUp(self):

        l.basicConfig(level=l.DEBUG,
                       format='%(asctime)s %(levelname)s %(message)s')
        self.chlfi = LFI()['LFI28M']
        self.TOLERANCE = 1e-7
        #from healpix
        self.ecl2gal_healpix = np.matrix([
                                        [ -5.48824860e-02,  -9.93821033e-01,  -9.64762490e-02],
                                        [  4.94116468e-01,  -1.10993846e-01,   8.62281440e-01],
                                        [ -8.67661702e-01,  -3.46354000e-04,   4.97154957e-01]
                                        ])

    def test_prepare_for_dipole(self):
        '''Saves pointing array for testing dipole generation'''
        read_exchange([self.chlfi], ods = [100], discard_flag = False,type='R')
        pnt = Pointing(self.chlfi.f.obtx,coord='G')
        vec = pnt.get(self.chlfi)
        np.save('vec_LFI28M_OD100_G',vec)

    def test_prepare_for_M3(self):
        '''Saves pointing array for comparing with M3'''
        read_exchange([self.chlfi], ods = [100], discard_flag = False,type='R')
        span = 56300
        obt = self.chlfi.f.obtx[:span]
        pnt = Pointing(obt, coord='E')
        vec = pnt.get(self.chlfi)
        np.save('vec_LFI28M_OD100_E',vec)
        np.save('obt_LFI28M_OD100_E',obt)

    def test_pointing_vs_dpc(self):
        '''Check pointing against LFI DPC'''
        # pointing extracted from DPC trieste
        dpc=cPickle.load(open('/u/zonca/p/issues/pointing/pointing100_DPC.pkl'))
        i_dpc = dpc['sampleOBT'].searchsorted(106793429004442.0)
        thetadpc=dpc['theta'][i_dpc]
        phidpc=dpc['phi'][i_dpc]
        psidpc=dpc['psi'][i_dpc]
        dpc['psi'] = dpc['psi'][i_dpc:]
        dpc['theta'] = dpc['theta'][i_dpc:]
        dpc['phi'] = dpc['phi'][i_dpc:]
        dpc['sampleOBT'] = dpc['sampleOBT'][i_dpc:]/2**16
        dpc['name'] = 'dpc'

        read_exchange([self.chlfi], ods = [100], discard_flag = False,type='R')
        i_te = 697
        te = {'sampleOBT':self.chlfi.f.obtx[i_te:],'name':'te'}
        #pnt = Pointing(self.chlfi.f.obtx[697:698],coord='E')
        pnt = Pointing(te['sampleOBT'],coord='E')
        vec = pnt.get(self.chlfi)
        thetav, phiv, psiv = pnt.get_3ang(self.chlfi)
        te['theta'] = thetav
        te['phi'] = phiv
        theta, phi, psi = thetav[0], phiv[0], psiv[0]
        print('Theta DPC %f testenv %f' % (thetadpc, theta))
        print('Phi DPC %f testenv %f' % (phidpc, phi))
        span = 32.5 * 60 * 10
        for angle in ['theta','phi','diff_theta','diff_phi']:
            plt.figure()
            if angle.startswith('diff'):
                a = angle.split('_')[1]
                plt.plot(d['sampleOBT'][:span], dpc[a][:span]-te[a][:span],label=a.capitalize() + ' difference')
            else:
                for d in [dpc,te]:
                    plt.plot(d['sampleOBT'][:span], d[angle][:span], label=d['name'])
            plt.legend();plt.grid()
            plt.title(angle.capitalize().replace('_',' '))
            plt.xlabel('OBT')
            plt.ylabel('%s [rad]' % angle)
            plt.savefig('%s_dpc_te.png' % angle)
        assert abs(thetadpc - theta) < self.TOLERANCE
        assert abs(phidpc - phi) < self.TOLERANCE
        assert abs(psidpc - psi) < self.TOLERANCE

    def test_quaternion_ecl2gal(self):
        q = np.array([0, 0, 0, 1])
        vecl = np.array([ 0.29192658,  0.45464871,  0.84147098])
        qgal = quaternion_ecl2gal(q)
        vgal_ecl2gal = qarray.rotate(qgal,vecl).flatten()
        vgal_matrix = np.asarray(self.ecl2gal_healpix * vecl[:,np.newaxis]).flatten()
        np.testing.assert_array_almost_equal(vgal_ecl2gal , vgal_matrix)

    def test_vector_ecl2gal(self):
        vecl = np.array([ 0.29192658,  0.45464871,  0.84147098])
        vgal_ecl2gal = vector_ecl2gal(vecl).flatten()
        vgal_matrix = np.asarray(self.ecl2gal_healpix * vecl[:,np.newaxis]).flatten()
        np.testing.assert_array_almost_equal(vgal_ecl2gal , vgal_matrix)

    def test_iqu_vs_m3(self):
        obt = np.array([1629538385.0881653, 1629538385.3650208])
        pnt = Pointing(obt)
        pix, qw, uw = pnt.get_pix_iqu(self.chlfi)
        np.testing.assert_array_almost_equal(pix, [6100717, 6102734])
        np.testing.assert_array_almost_equal(qw, [-8.252502679824829e-01, -8.330312371253967e-01])
        np.testing.assert_array_almost_equal(uw, [5.647672414779663e-01, 5.532259941101074e-01])

if __name__ == '__main__':
    unittest.main()
