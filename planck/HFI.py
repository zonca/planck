#!/usr/bin/env python
#
# Generic python class for dealing with Planck HFI
# by zonca@deepspace.ucsb.edu

import Planck
import private

class HFIChannel(Planck.Channel):

    @property
    def centralfreq(self):
        return self.f.freq

    def Planck_to_RJ(self, data):
        return data / private.mKRJ_2_mKcmb[self.f.freq]

class HFI(Planck.Instrument):
    

    uncal = 'R'

    Channel = HFIChannel
    
    def __init__(self, name = 'HFI', rimo =private.HFI_rimo):
        super(HFI, self).__init__(name,rimo)

    @staticmethod
    def freq_from_tag(tag):
        return int(tag[:3])
