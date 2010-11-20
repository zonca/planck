import cPickle
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


    def test_prepare_for_dipole(self):
        '''Saves pointing array for testing dipole generation'''
        read_exchange([self.chlfi], ods = [100], discard_flag = False,type='R')
        pnt = Pointing(self.chlfi.f.obtx,coord='G')
        vec = pnt.get(self.chlfi)
        np.save('vec_LFI28M_OD100_G',vec)
        assert True

    def test_prepare_for_M3(self):
        '''Saves pointing array for comparing with M3'''
        read_exchange([self.chlfi], ods = [100], discard_flag = False,type='R')
        span = 56300
        obt = self.chlfi.f.obtx[:span]
        pnt = Pointing(obt, coord='E')
        vec = pnt.get(self.chlfi)
        np.save('vec_LFI28M_OD100_E',vec)
        np.save('obt_LFI28M_OD100_E',obt)
        assert True

    def test_pointing_vs_dpc(self):
        '''Check pointing against LFI DPC'''
        # pointing extracted from DPC trieste
        dpc=cPickle.load(open('/u/zonca/p/issues/pointing/pointing100_DPC.pkl'))
        i_dpc = dpc['sampleOBT'].searchsorted(106793429004442.0)
        thetadpc=dpc['theta'][i_dpc]
        phidpc=dpc['phi'][i_dpc]
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
        thetav, phiv = healpy.vec2ang(vec)
        te['theta'] = thetav
        te['phi'] = phiv
        theta, phi = thetav[0], phiv[0]
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

    def test_quaternion_ecl2gal(self):
        #from Quaternion module .transform
        ecl2gal_matrix = np.matrix([[ -5.48755398e-02,  -9.93821384e-01,  -9.64765918e-02],
                                    [  4.94109453e-01,  -1.10990693e-01,   8.62285866e-01],
                                    [ -8.67666136e-01,  -3.51593739e-04,   4.97147215e-01]])
        q = np.array([0, 0, 0, 1])
        vecl = np.array([ 0.29192658,  0.45464871,  0.84147098])
        qgal = quaternion_ecl2gal(q)
        vgal_ecl2gal = qarray.rotate(qgal,vecl)
        vgal_matrix = np.asarray(ecl2gal_matrix * vecl[:,np.newaxis]).flatten()
        print(vgal_ecl2gal)
        print(vgal_matrix)
        assert (vgal_ecl2gal - vgal_matrix).std() < 1e-8

    def test_vector_ecl2gal(self):
        #from Quaternion module .transform
        ecl2gal_matrix = np.matrix([[ -5.48755398e-02,  -9.93821384e-01,  -9.64765918e-02],
                                    [  4.94109453e-01,  -1.10990693e-01,   8.62285866e-01],
                                    [ -8.67666136e-01,  -3.51593739e-04,   4.97147215e-01]])
        vecl = np.array([ 0.29192658,  0.45464871,  0.84147098])
        vgal_ecl2gal = vector_ecl2gal(vecl)
        vgal_matrix = np.asarray(ecl2gal_matrix * vecl[:,np.newaxis]).flatten()
        print(vgal_ecl2gal)
        print(vgal_matrix)
        assert (vgal_ecl2gal - vgal_matrix).std() < 1e-8
