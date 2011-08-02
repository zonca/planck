import pyfits
import numpy as np
from exceptions import KeyError
import itertools
import operator

def group_by_horn(chlist):
    return itertools.groupby(chlist, operator.attrgetter('horn'))
    

def parse_channels(chfreq):
    if chfreq is None:
        return None
    if isinstance(chfreq, list) and isinstance(chfreq[0], Channel):
        return chfreq
    else:
        pl = Planck()
        if isinstance(chfreq, int):
            return pl.f[chfreq].ch
        elif isinstance(chfreq, str):
            return [pl[chfreq]]
        elif isinstance(chfreq, list) and isinstance(chfreq[0], str):
            return [pl[tag] for tag in chfreq]

def parse_channel(chfreq):
    if isinstance(chfreq, Channel):
        return chfreq
    else:
        pl = Planck()
        if isinstance(chfreq, str):
            return pl[chfreq]

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
    def arm(self):
        return self.tag[-1]

    @property
    def sampling_freq(self):
        return float(self.rimo['F_SAMP'])

    @property
    def white_noise(self): 
        return float(self.rimo['NET'])

    @property
    def pair(self):
        return self.inst[self.tag.replace(self.MS[self.n], self.MS[not self.n])]

    @property
    def n(self):
        return self.fromMS[self.arm]

    def get_beam_real(self, m_b, component='main'):
        """Return real"""
        if m_b == -1:
            return 2 * -1 * private.BEAM[component][self.tag][-m_b]
        else:
            return 2 * private.BEAM[component][self.tag][m_b]

    def get_beam_imag(self, m_b, component='main'):
        """Return imaginary"""
        if m_b == -1:
            return private.BEAM[component][self.tag][2]
        elif m_b == 0:
            return 0
        elif m_b == 1:
            return private.BEAM[component][self.tag][2]

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

def freq2inst(freq):
    return ['LFI','HFI'][freq>=100]
EXCLUDED_CH = ['143-8', '545-3', '857-4']
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
        self.ch = [ch for ch in self.ch if ch.tag not in EXCLUDED_CH]
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
import private


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

