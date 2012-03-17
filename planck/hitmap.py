from __future__ import division
import healpy
import cPickle
import logging as l
import numpy as np
from LFI import LFI
from pointing import Pointing, DiskPointing
from testenv.remix import read_exchange
import glob

def testBit(int_type, offset):
    """returns a nonzero result, 2**offset, if the bit at 'offset' is one."""
    mask = 1 << offset
    return(int_type & mask)

def concat_hitmaps(folder = 'pkl/'):
    files = glob.iglob(folder + '*pkl')
    for f in files:
        print(f)
        odhitmap = cPickle.load(open(f,'rb'))
        try:
            hitmap += odhitmap
        except:
            hitmap = odhitmap
    return hitmap 
    

class HitMap(object):

    def __init__(self, freq, od, nside=1024, use_flag=True,folder=None):
        LOG_FILENAME = '/project/projectdirs/planck/user/zonca/issues/hitmap/full.log'
        l.basicConfig(filename=LOG_FILENAME,level=l.DEBUG)
        self.freq = freq
        self.od = od
        self.use_flag = use_flag
        self.lfi = LFI()
        self.f = self.lfi.f[self.freq]
        self.nside = nside
        #self.ecl2gal = healpy.Rotator(coord=['E','G'])
        self.folder = folder
        l.info('%s ready' % self)

    def __repr__(self):
        return 'HitMap %d GHz, od %d' % (self.freq, self.od)

    def run(self):
        print('Reading data')
        read_exchange(self.f.ch, ods = [self.od], discard_flag = False,type='C')
        obt_good = testBit(self.f.commonflag,0)==0
        print('Preparing pointing')
        self.pnt = DiskPointing(self.od, self.freq, folder=self.folder)
        self.hitmap = np.zeros(healpy.nside2npix(self.nside))
        for ch in self.f.ch: 
            print('Processing rad %s' % ch)
            #theta, phi = self.ecl2gal(*self.pnt.get_ang(ch))
            theta, phi = self.pnt.get_ang(ch)
            ch_good = obt_good & (ch.flagx == 0) & (ch.pair.flagx == 0)
            print('Flagged %.2f perc' % (100*(len(ch_good)-ch_good.sum())/len(ch_good)))
            if len(theta[ch_good]) > 0:
                ids = np.bincount(healpy.ang2pix(self.nside, theta[ch_good], phi[ch_good], nest=True))
                self.hitmap[:len(ids)] += ids
            else:
                print('Skip channel')
        print('Writing to file')
        cPickle.dump(self.hitmap, open('/project/projectdirs/planck/user/zonca/issues/hitmap/pkl/%d_%d.pkl' % (self.freq,self.od),'wb'),protocol=-1)

def pix2map(pix, nside, tod=None):
    """Pixel array to hitmap, if TOD with same lenght of PIX is provided, 
    it is binned to a map"""
    #TODO test case
    pix = pix.astype(np.int)
    ids = np.bincount(pix, weights=None)
    hitmap = np.ones(healpy.nside2npix(nside)) * healpy.UNSEEN
    hitmap[:len(ids)] = ids
    hitmap = healpy.ma(hitmap)
    if tod is None:
        return hitmap
    else:
        ids_binned = np.bincount(pix, weights=tod)
        binned = np.ones(healpy.nside2npix(nside)) * healpy.UNSEEN
        binned[:len(ids_binned)] = ids_binned
        binned = healpy.ma(binned)/hitmap
        return hitmap, binned

    
