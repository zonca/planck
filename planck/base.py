from astropy.io import fits as pyfits
import logging as l
import numpy as np
try:
    from exceptions import KeyError, ValueError
except:
    pass
import itertools
import operator
import collections
from . import private

def group_by_horn(chlist):
    return itertools.groupby(chlist, operator.attrgetter('horn'))
    
class ChannelBase(object):
    '''Base for Channel, frequencyset and detector'''

    def __repr__(self):
        return self.tag

    def get_gaussian_beam(self, lmax=1024, pol=False, beam_eff=False):
        """Equivalent gaussian beam from RIMO FWHM

        Returns the transfer function of a gaussian beam until lmax,
        either polarized or not"""
        import healpy as hp
        beam = hp.gauss_beam( self.fwhm, lmax, pol) 

        if beam_eff:
            beam *= self.beam_efficiency
        return beam

class Channel(ChannelBase):
    '''Abstract channel class for LFI and HFI channels'''

    def __init__(self, data, inst=None):
        try:
            self.tag = data[0].strip().decode('unicode_escape')
        except:
            self.tag = data[0].strip()
        self.rimo = np.array(data)
        self.inst = inst

    @property
    def num(self):
        return self.f.ch.index(self)

    @property
    def arm(self):
        return self.tag[-1]

    @property
    def sampling_freq(self):
        return float(self.rimo['F_SAMP'])

    @property
    def white_noise(self): 
        return float(self.rimo[self.inst.white_noise_field])

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

    @property
    def fwhm(self):
        fwhm = self.get_instrument_db_field("FWHM")
        l.info("Channel %s: FWHM %.2f arcmin" % (self.tag, fwhm))
        return np.radians(fwhm/60.),

    @property
    def beam_efficiency(self):
        import private
        try:
            return private.beam_efficiency[self.tag] / 100.
        except KeyError:
            l.warning("Missing beam efficiency for channel %s" % self.tag)
            return 1.

    def get_instrument_db_field(self, field): 
        try:
            return self.inst.instrument_db(self)[field][0]
        except ValueError: 
            return self.inst.instrument_db(self)[field.upper()][0]
        
class FrequencySet(ChannelBase):
    def __init__(self, freq, ch, inst=None):
        self.freq = freq
        self.ch = ch
        self.inst = inst
        for ch in self.ch:
            ch.f = self
        self.tag = '%d' % self.freq 
        self.f = self

    @property
    def horns(self):
        return group_by_horn(self.ch)

    def __repr__(self):
        return '%d GHz' % self.freq

    @property
    def sampling_freq(self):
        return self.ch[0].sampling_freq

    def get_aggregated_property(self, name):
        name_attrgetter = operator.attrgetter(name)
        return np.mean([name_attrgetter(ch) for ch in self.ch])

    def __getattr__(self, name):
        return self.get_aggregated_property(name)

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
        f = collections.OrderedDict() 
        for freq in sorted(set(freqs)):
            chlist = [self.ch[i] for i,chfreq in enumerate(freqs) if chfreq == freq]
            f[freq] = self.FrequencySet(freq, chlist, self)
        return f
        
    def instrument_db(self,ch):
        if not hasattr(self,'_instrument_db') or self._instrument_db is None:
            if isinstance(self.instrument_db_file, list):
                self.instrument_db_file = self.instrument_db_file[0]
            if isinstance(self.instrument_db_file, dict):
                self.instrument_db_file = self.instrument_db_file[self.name]
            self._instrument_db = np.array(pyfits.open(self.instrument_db_file,ignore_missing_end=True)[1].data)
            l.warning('Loading instrumentdb %s' % self.instrument_db_file)
        try:
            det_index, = np.where([rad.strip().endswith(ch.tag.encode()) for rad in self._instrument_db['Radiometer']])
        except ValueError:
            det_index, = np.where([rad.strip().endswith(ch.tag.encode()) for rad in self._instrument_db['DETECTOR']])
        return self._instrument_db[det_index]
            
    def __getitem__(self, key):
        return self.chdict[key]
