import glob
import os
import sqlite3
import exceptions
import logging as l
from collections import namedtuple

from Planck import parse_channels, freq2inst

try:
    import private
except exceptions.ImportError:
    print('private.py is needed to use the planck module, this is available only to members of the Planck collaboration')
    raise

Period = namedtuple('Period', ['number','start','stop'])
Observation = namedtuple('Observation', ['od','start','stop','PP','EFF'])

def obt2od(obt, freq=30):
    """Get precise OD from obt stamp2"""
    conn = sqlite3.connect(private.database)
    c = conn.cursor()
    query = c.execute('select startOBT, endOBT, od from efdd_od_ranges where startOBT<%d and endOBT>%d' % (obt,obt))
    q = query.fetchone()
    od = int(q[2]) + float((obt-q[0]))/(q[1]-q[0])
    c.close()
    return od

class DataSelector(object):
    """Planck data selector
    channels can be integer frequency, list of channel names (same frequency) or a single channel name string
    efftype is R for reduced, C for converted [uncal with dip]
    ods are the AHF ODs
    eff_ods are the EFF ODs"""

    def __init__(self, channels=None, efftype='R'):
        self.channels = parse_channels(channels)
        self.f = self.channels[0].f
        self.config = {}
        self.config['database'] = private.database
        self.config['exchangefolder'] = private.exchangefolder[self.f.inst.name]
        self.config['ahf_folder'] = private.AHF
        self.efftype = efftype

    @property
    def ods(self):
        try:
            return self._ods
        except exceptions.AttributeError:
            self._ods = set()
            for obt_range in self.obt_ranges:
                self._ods.update(self.get_ODs(obt_range))
            self._ods = sorted(list(self._ods))
            return self.ods

    @property
    def obt_ranges(self):
        try:
            return self._obt_ranges
        except exceptions.AttributeError:
            self._obt_ranges = map(get_obt_range_from_od, self.ods)
            return self._obt_ranges

    def by_od_range(self, od_range):
        """ods is a list of start and stop OD (INCLUDED)"""
        self.od_range = od_range
        self._ods = range(od_range[0], od_range[1]+1)
    

    def by_obt(self, obt_ranges):
        """obt_ranges is a list of 2 element lists (or tuple) with start-stop obt"""
        raise exceptions.NotImplementedError()
        self._obt_ranges = obt_ranges

    def by_rings(self, rings):
        """rings is a list of ring numbers"""
        #FIXME LFI or HFI?
        raise exceptions.NotImplementedError()

    def get_AHF_ods(self, obt_range):
        conn = sqlite3.connect(self.config['database'])
        c = conn.cursor()
        values = (obt_range[0]*2**16, obt_range[-1]*2**16)
        query = c.execute('select od from ahf_files where endOBT>=? and startOBT<=?', values)
        ods = [int(q[0]) for q in query]
        c.close()
        return ods

    def get_one_AHF(self, obt_range):
        ods = self.get_AHF_ods(obt_range)
        files = [glob.glob(
            os.path.join(self.config['ahf_folder'], '%04d' % od, 'vel*')
            )[0] for od in ods]
        return files

    def get_AHF(self):
        return [self.get_one_AHF(obt_range) for obt_range in self.obt_ranges]

    def get_EFF(self):
        return self.latest_exchange(self.eff_ods)

    def get_OBS(self):
        OBS = []
        for od, obt_range in zip(self.ods, self.obt_ranges):
            OBS.append(Observation(od=od, start=obt_range[0], stop=obt_range[1], PP=self.get_PP(od), EFF=self.latest_exchange(eff_ods_from_obt_range(self.f.freq, obt_range))))
        return OBS

    def get_PP(self, od):
        conn = sqlite3.connect(self.config['database'])
        c = conn.cursor()
        query = c.execute('select pointID_unique, start_time, end_time from list_ahf_infos where od==? order by pointID_unique ASC', (str(od),))
        PP = [Period(int(q[0]),q[1]/2.**16,q[2]/2.**16) for q in query]
        c.close()
        return PP

    #def get_obt_range_from_ods(self, freq=None, ods=None):
    #    if ods is None:
    #        ods = self.ods
    #    if freq is None:
    #        freq = self.f.freq
    #    conn = sqlite3.connect(self.config['database'])
    #    c = conn.cursor()
    #    query = c.execute('select startOBT, endOBT from efdd_od_ranges where freq=? and od in ? order by od ASC', (freq, ods))
    #    c.close()
    #    obt_ranges = [[q[0]/2.**16, q[1]/2.**16] for q in query]
    #    return obt_ranges

    #def get_obt_range_from_od_range(self, freq=None, od_range=None):
    #    if od_range is None:
    #        od_range = self.od_range
    #    if freq is None:
    #        freq = self.f.freq
    #    conn = sqlite3.connect(self.config['database'])
    #    c = conn.cursor()
    #    obt_range = []
    #    #start of first pp
    #    query = c.execute('select start_time from list_ahf_infos where od==? order by pointID_unique ASC limit 1', od_range[0])
    #    obt_range.append(query.fetchone()[0]/2.**16)
    #    #end of last pp
    #    query = c.execute('select end_time from list_ahf_infos where od==? order by pointID_unique DESC limit 1', od_range[1])
    #    obt_range.append(query.fetchone()[0]/2.**16)
    #    c.close()
    #    return obt_range


    @property
    def eff_ods(self):
        """List of ODs within the obt range provided"""
        return eff_ods_from_obt_range(self.f.freq, [self.obt_ranges[0][0], self.obt_ranges[-1][-1]])

    def latest_exchange(self, od):
        return latest_exchange(self.f.freq, od, self.config['exchangefolder'], self.efftype)

def eff_ods_from_obt_range(freq, obt_range, database=None):
    conn = sqlite3.connect(database or private.database)
    c = conn.cursor()
    query = c.execute('select od from efdd_od_ranges where freq=? and endOBT>=? and startOBT<=?',(freq, obt_range[0]*2**16, obt_range[1]*2**16) )
    eff_ods = [q[0] for q in query]
    c.close()
    return eff_ods

def get_obt_range_from_od(od, database=None):
    conn = sqlite3.connect(database or private.database)
    c = conn.cursor()
    query = c.execute('select start_time,end_time from list_ahf_infos where od==? order by pointID_unique ASC', (str(od),))
    all = query.fetchall()
    obt_range = (all[0][0]/2.**16, all[-1][-1]/2.**16)
    c.close()
    return obt_range

def latest_exchange(freq, ods, exchangefolder = None, type = 'R'):
    """Returns the latest version of an exchange format file"""
    single = False
    if exchangefolder is None:
        exchangefolder = private.exchangefolder[freq2inst(freq)]

    if isinstance(ods, int):
        ods = [ods]
        single = True
    if glob.glob(exchangefolder + '/*.fits'):
        ods = [0]
    if type == 'K':
        type = ''
    EFF = []
    for od in ods:
        pattern = '/*%03d-%04d-%s*.fits' % (freq, od, type)
        if od:
            pattern = ('/%04d' % od) + pattern
        l.debug('Exchange format: %s' % (exchangefolder + pattern))
        allversions = glob.glob(exchangefolder + pattern)
        if not allversions:
            error_message = 'Cannot find file from pattern: %s' % (exchangefolder +         pattern)
            l.error(error_message)
            raise exceptions.IOError(error_message)
        EFF.append(max(allversions, key = lambda x: int(x[-13:-5])))
    if single:
        EFF = EFF[0]
    return EFF

if __name__ == '__main__':
    ds = DataSelector(channels=100)
    ds.config['exchangefolder'] = 'toast_planck_data/hfi_m3_v41'
    ds.by_od_range([97, 103])
    print(ds.get_EFF())
    print(ds.get_AHF())
    print(ds.get_OBS())
    obs=ds.get_OBS()
