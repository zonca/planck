#!/usr/bin/env python
#
# Generic python class for dealing with Planck LFI
# by zonca@deepspace.ucsb.edu

import numpy as np
import Planck
import private

def flatten_d(chlist):
    return [d for ch in chlist for d in ch.d]

class LFIChannel(Planck.Channel):

    MS = { 0 : 'M', 1 : 'S' }
    fromMS = { 'M' : 0, 'S' : 1}

    def __init__(self, data, inst=None):
        super(LFIChannel, self).__init__(data, inst)
        self.d = [Detector(self, 0), Detector(self, 1)]

    @property
    def RCA(self):
        return LFI.RCA_from_tag(self.tag)

    @property
    def arm(self):
        return self.tag[-1]

    @property
    def n(self):
        return self.fromMS[self.arm]

    @property
    def centralfreq(self):
        return self.get_instrument_db_field('nu_cen')

    @property
    def wn(self): 
        return self.get_instrument_db_field('net_KCMB')

    def get_instrument_db_field(self, field): 
        return self.inst.instrument_db(self)[field][0]

    def __getitem__(self, n):
        return self.d[n]

    def Planck_to_RJ(self, data):
        import dipole
        return dipole.Planck_to_RJ(data, self.centralfreq)

    @property
    def pair(self):
        return self.inst[self.tag.replace(self.MS[self.n], self.MS[not self.n])]

class LFIFrequencySet(Planck.FrequencySet):
    
    @property
    def d(self):
        return flatten_d(self.ch)
class LFI(Planck.Instrument):

    uncal = 'C'

    Channel = LFIChannel
    FrequencySet = LFIFrequencySet

    def __init__(self, name = 'LFI', rimo = private.LFI_rimo):
        super(LFI, self).__init__(name,rimo)

    @classmethod
    def freq_from_tag(cls, tag):
        RCA = cls.RCA_from_tag(tag)
        if RCA <= 23:
            return 70
        elif RCA <= 26:
            return 44
        elif RCA <= 28:
            return 30
        else:
            return None
        
    @staticmethod
    def RCA_from_tag(tag):
        return int(tag[3:5])
        
    def instrument_db(self,ch):
        if not hasattr(self,'_instrument_db') or self._instrument_db is None:
            import pyfits
            self._instrument_db = np.array(pyfits.open(private.instrument_db,ignore_missing_end=True)[1].data)
        det_index, = np.where([rad.strip().endswith(ch.tag) for rad in self._instrument_db['Radiometer']])
        return self._instrument_db[det_index]

    @property
    def d(self):
        return flatten_d(self.ch)
class Detector(Planck.ChannelBase):
    def __init__(self, ch, n):
        self.n = n
        self.ch = ch
        self.tag = '%s-%s%s' % (ch.tag, self.ch.n, self.n)

    def savfilename(self, od):
        return private.savfilename % (od, self.ch.RCA, self.ch.n, self.n)

    @property
    def cdstag(self):
        return 'RCA%s%s%s' % (self.ch.RCA,self.ch.n, self.n)
