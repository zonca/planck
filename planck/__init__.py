from .lfi import LFI
from .hfi import HFI
from .base import Channel

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


class Planck(object):
    '''Planck class, gives an iterator .ch for all LFI and HFI channels'''
    
    def __init__(self):
        self.inst = {'LFI':LFI(), 'HFI':HFI()}
        self.ch = [ch for inst in self.inst.values() for ch in inst.ch]
        self.f = dict((freq,f) for inst in self.inst.values() for freq,f in inst.f.items())

    def __getitem__(self, key):
        try:
            return self.inst['LFI'][key]
        except KeyError:
            return self.inst['HFI'][key]
