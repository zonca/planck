import exceptions
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

Period = namedtuple('Period', ['number','start','stop','splitnumber'])
Observation = namedtuple('Observation', ['od','tag','start','stop','PP','EFF', 'break_startrow', 'break_stoprow'])

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
    eff_ods are the EFF ODs
    
    some configuration options can be modified after creating the object by accessing the .config dictionary"""

    def __init__(self, channels=None, efftype='R', include_preFLS=False):
        self.include_preFLS = include_preFLS
        self.channels = parse_channels(channels)

        if self.channels is None:
            self.f = None
        else:
            self.f = self.channels[0].f
        self.config = {}
        self.config['database'] = private.database
        self.config['breaks'] = private.TOAST['breaks']
        if self.channels is None:
            self.config['exchangefolder'] = None
        else:
            self.config['exchangefolder'] = private.exchangefolder[self.f.inst.name]

        self.config['ahf_folder'] = private.AHF
        self.config['exclude_454_455'] = True
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
        if self.config['exclude_454_455'] and self.f.inst.name == 'LFI':
            for OD in [454, 455]:
                try:
                    self._ods.remove(OD)
                except:
                    pass
    

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
        values = ((obt_range[0]-1)*2**16, (obt_range[-1]+1)*2**16)
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

    def get_OD_OBS(self):
        """Gets one observation for each Operational Day"""
        OBS = []
        for od, obt_range in zip(self.ods, self.obt_ranges):
            eff_ods = eff_ods_from_obt_range(self.f.freq, obt_range)
            if self.config['exclude_454_455'] and self.f.inst.name == 'LFI':
                for OD in [454, 455]:
                    try:
                        eff_ods.remove(OD)
                    except:
                        pass
            OBS.append(Observation(od=od, tag='', start=obt_range[0], stop=obt_range[1], PP=self.get_PP(od), EFF=self.latest_exchange(eff_ods), break_startrow=None, break_stoprow=None))
        return OBS

    def get_OBS(self):
        """Checks the breaks table and splits the observations accordingly"""
        OBS = self.get_OD_OBS()

        conn = sqlite3.connect(self.config['breaks'])
        c = conn.cursor()
        #obt are in clocks in the database
        query = c.execute('select od, startobt, stopobt, startrow, stoprow from eff_breaks where freq=? and startobt > ? and stopobt < ?', (self.f.freq,OBS[0].start*2**16, OBS[-1].stop*2**16))
        for od, startobt, stopobt, startrow, stoprow in query:
            #convert in seconds
            startobt /= 2.**16
            stopobt /= 2.**16
            l.warning('Break found in OD %d' % od)
            try:
                OB = [o for o in OBS if o.start < stopobt and o.stop > startobt][0]
            except exceptions.IndexError:
                l.error('Cannot identify the observation related to the break in OD %d' % od)
                sys.exit(1)
            i = OBS.index(OB)
            splitted_OBS = split_observation(OB, startobt, stopobt, startrow, stoprow)
            OBS[i] = splitted_OBS[1]
            OBS.insert(i, splitted_OBS[0])
        c.close()
        return OBS

    def get_PP(self, od):
        """Gets all the pointing periods in one Operational Day, returns a list of Period named tuples"""
        conn = sqlite3.connect(self.config['database'])
        c = conn.cursor()
        if self.include_preFLS:
            query = c.execute('select pointID_unique, start_time, end_time from list_ahf_infos where od==? AND start_time < end_time order by start_time ASC', (str(od),))
        else:
            query = c.execute('select pointID_unique, start_time, end_time from list_ahf_infos where od==? and start_time > 106743579730069 AND  start_time < end_time order by start_time ASC', (str(od),))

        PP = []
        for q in query:
            pid_numbers= map(int, q[0].split('-'))
            if len(pid_numbers) == 1:
                pid_numbers.append(0)
            PP.append( Period(pid_numbers[0],q[1]/2.**16,q[2]/2.**16, pid_numbers[1]) )
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
        eff_ods = eff_ods_from_obt_range(self.f.freq, [self.obt_ranges[0][0], self.obt_ranges[-1][-1]])
        if self.config['exclude_454_455'] and self.f.inst.name == 'LFI':
            for OD in [454, 455]:
                try:
                    eff_ods.remove(OD)
                except:
                    pass
        return eff_ods

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
    query = c.execute('select start_time,end_time from list_ahf_infos where od==? order by start_time ASC', (str(od),))
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
        pattern = '*%03d-%04d-%s*.fits' % (freq, od, type)
        if od:
            pattern = os.path.join('%04d' % od, pattern)
        l.debug('Exchange format: %s' % (os.path.join(exchangefolder, pattern)))
        allversions = glob.glob(os.path.join(exchangefolder, pattern))
        if not allversions:
            error_message = 'Cannot find file from pattern: %s' % os.path.join(exchangefolder, pattern)
            l.error(error_message)
            raise exceptions.IOError(error_message)
        EFF.append(max(allversions, key = lambda x: int(x[-13:-5])))
    if single:
        EFF = EFF[0]
    return EFF

def split_observation(OB, startobt, stopobt, startrow, stoprow):
    """Splits one observation in 2 observations"""
    PP1 = [p for p in OB.PP if p.start < startobt]
    PP1[-1] = Period(PP1[-1].number, PP1[-1].start, startobt, PP1[-1].splitnumber)
    OB1 = Observation(od=OB.od, tag=OB.tag + 'a', start=OB.start, stop=startobt, PP=PP1, EFF=OB.EFF, break_startrow=startrow, break_stoprow=None)

    PP2 = [p for p in OB.PP if p.stop > stopobt]
    PP2[0] = Period(PP2[0].number, stopobt, PP2[0].stop, PP2[0].splitnumber)
    OB2 = Observation(od=OB.od, tag=OB.tag + 'b', start=stopobt, stop=OB.stop, PP=PP2, EFF=OB.EFF, break_startrow=None, break_stoprow=stoprow)

    return OB1, OB2


if __name__ == '__main__':
    ds = DataSelector(channels=30)
    ds.config['exchangefolder'] = '/global/scratch/sd/planck/user/zonca/data/LFI_DX7S_conv/'
    ds.by_od_range([452, 457])
    print(ds.get_EFF())
    print(ds.get_AHF())
    print(ds.get_OBS())
    obs=ds.get_OBS()
