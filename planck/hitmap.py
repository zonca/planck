import matplotlib
matplotlib.use('Agg')
import healpy
import cPickle
import logging as l
import numpy as np
from LFI import LFI
from pointing import Pointing
from testenv.remix import read_exchange
import glob

def testBit(int_type, offset):
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

    def __init__(self, freq, od, nside=1024, use_flag=True):
        LOG_FILENAME = '/u/zonca/p/issues/hitmap/full.log'
        l.basicConfig(filename=LOG_FILENAME,level=l.DEBUG)
        self.freq = freq
        self.od = od
        self.use_flag = use_flag
        self.lfi = LFI()
        self.f = self.lfi.f[self.freq]
        self.nside = nside
        l.info('%s ready' % self)

    def __repr__(self):
        return 'HitMap %d GHz, od %d' % (self.freq, self.od)

    def run(self):
        l.debug('Reading data')
        read_exchange_obt_flag(self.freq, [self.f.r[0]], ods = [self.od], discard_flag = False,type='R')
        obt = np.ma.masked_where(testBit(self.f.commonflag,0)!=0,self.f.obtx)
        l.debug('Preparing pointing')
        self.pnt = Pointing(obt.compressed(),coord='G')
        self.hitmap = np.zeros(healpy.nside2npix(self.nside))
        for rad in self.f.r: 
            l.debug('Processing rad %s' % rad)
            vec = self.pnt.get(rad)
            ids = np.bincount(healpy.vec2pix(self.nside, vec[:,0], vec[:,1], vec[:,2]))
            self.hitmap[:len(ids)] += ids
        l.debug('Writing to file')
        cPickle.dump(self.hitmap, open('/u/zonca/p/issues/hitmap/pkl/%d_%d.pkl' % (self.freq,self.od),'wb'),protocol=-1)

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

    
