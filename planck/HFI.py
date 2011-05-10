#!/usr/bin/env python
#
# Generic python class for dealing with Planck HFI
# by zonca@deepspace.ucsb.edu

import numpy as np
import Planck
import private

class HFIChannel(Planck.Channel):

    MS = { 0 : 'a', 1 : 'b' }
    fromMS = { 'a' : 0, 'b' : 1}

    @property
    def centralfreq(self):
        return self.f.freq

    @property
    def calibdiff(self):
        return self.diff / self.inst.cal['cal'][self.inst.cal['ch']==self.tag]

    @property
    def horn(self):
        return int(self.tag.replace('a','').replace('b','').replace('-',''))

    def Planck_to_RJ(self, data):
        return data / private.mKRJ_2_mKcmb[self.f.freq]

class HFI(Planck.Instrument):
    
    uncal = 'R'

    Channel = HFIChannel

    def load_cal(self):
        self.cal = np.loadtxt(private.HFI_calibfile,dtype=[('ch','S8'),('cal',np.double)])
    
    def __init__(self, name = 'HFI', rimo =private.HFI_rimo):
        super(HFI, self).__init__(name,rimo)
        self.load_cal()

    @staticmethod
    def freq_from_tag(tag):
        return int(tag[:3])
