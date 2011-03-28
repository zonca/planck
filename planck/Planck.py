import pyfits
import numpy as np
from exceptions import KeyError
import itertools

class ChannelBase(object):
    '''Base for Channel, frequencyset and detector'''

    def __repr__(self):
        return self.tag

class Channel(ChannelBase):
    '''Abstract channel class for LFI and HFI channels'''

    def __init__(self, data, inst=None):
        self.tag = data[0].strip()
        self.rimo = np.array(data)
        self.inst = inst

    @property
    def sampling_freq(self):
        return self.rimo['F_SAMP']

class FrequencySet(ChannelBase):
    def __init__(self, freq, ch, inst=None):
        self.freq = freq
        self.ch = ch
        self.inst = inst
        for ch in self.ch:
            ch.f = self
        self.tag = '%d' % self.freq 

    def __repr__(self):
        return '%d GHz' % self.freq

    @property
    def sampling_freq(self):
        return self.ch[0].sampling_freq

    @property
    def wn(self):
        return np.mean([ch.wn for ch in self.ch])

class Instrument(object):
    '''Common base class for LFI and HFI'''

    Channel = Channel
    FrequencySet = FrequencySet

    def __init__(self, name, rimo):
        '''Rimo is full path to Reduced Instrument Model FITS file'''
        self.name = name
        self.rimo = rimo
        rimo_file = np.array(pyfits.open(rimo)[1].data)
        rimo_file.sort()
        self.rimo_fields = rimo_file.dtype.names
        self.ch = map(self.Channel, rimo_file, [self]*len(rimo_file))
        self.chdict = dict( (ch.tag, ch) for ch in self.ch)
        self.f = self.create_frequency_sets()

    def create_frequency_sets(self):
        freqs = [self.freq_from_tag(ch.tag) for ch in self.ch]
        f = {}
        for freq in set(freqs):
            chlist = [self.ch[i] for i,chfreq in enumerate(freqs) if chfreq == freq]
            f[freq] = self.FrequencySet(freq, chlist, self)
        return f
            
    def __getitem__(self, key):
        return self.chdict[key]
        
    
import LFI
import HFI

class Planck(object):
    '''Planck class, gives an iterator .ch for all LFI and HFI channels'''
    
    def __init__(self):
        self.inst = {'LFI':LFI.LFI(), 'HFI':HFI.HFI()}
        self.ch = [ch for inst in self.inst.values() for ch in inst.ch]
        self.f = dict((freq,f) for inst in self.inst.values() for freq,f in inst.f.iteritems())

    def __getitem__(self, key):
        try:
            return self.inst['LFI'][key]
        except KeyError:
            return self.inst['HFI'][key]
