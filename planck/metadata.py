import exceptions
from Planck import parse_channels

try:
    import private
except exceptions.ImportError:
    print('private.py is needed to use the planck module, this is available only to members of the Planck collaboration')
    raise

class DataSelector(object):
    """Planck data selector
    channels can be integer frequency, list of channel names (same frequency) or a single channel name string"""

    def __init__(self, channels=None):
        self.config = {}
        self.config['database'] = private.database
        self.config['exchangefolder'] = private.exchangefolder
        self.channels = parse_channels(channels)
        self.f = self.channels[0].f

    def by_od(self, ods):
        """ods is a list of operational days"""
        pass

    def by_obt(self, obt_ranges):
        """obt_ranges is a list of 2 element lists (or tuple) with start-stop obt"""
        self.obt_ranges = obt_ranges

    def by_rings(self, rings):
        """rings is a list of ring numbers"""
        #FIXME LFI or HFI?
        pass

    def get_one_AHF(self, obt_range):
        conn = sqlite3.connect(self.config['database'])
        c = conn.cursor()
        values = (obt_range[0]*2**16, obt_range[-1]*2**16)
        query = c.execute('select file_path from ahf_files where endOBT>=? and startOBT<=?', values)
        files = [q[0] for q in query]
        c.close()
        return files

    def get_AHF(self):
        return [self.get_one_AHF(obt_range) for obt_range in self.obt_ranges]
        

def ods_btw_obt(freq, obtrange):
    '''List of operational days in the obtrange supplied'''
    assert obtrange[1]>=obtrange[0]
    conn = sqlite3.connect(default_config()['database'])
    c = conn.cursor()
    query = c.execute('select od from efdd_od_ranges where freq=? and endOBT>=? and startOBT<=?',(freq, obtrange[0]*2**16, obtrange[1]*2**16) )
    ods = [q[0] for q in query]
    c.close()
    return ods

def latest_exchange(freq, od, exchangefolder = None, type = 'R'):
    '''Returns the latest version of an exchange format file'''
    if exchangefolder is None:
        exchangefolder = default_config()['exchangefolder'][freq2inst(freq)]
    if glob.glob(exchangefolder + '/*.fits'):
        od = 0
    if type == 'K':
        type = ''
    pattern = '/*%03d-%04d-%s*.fits' % (freq, od, type)
    if od:
        pattern = ('/%04d' % od) + pattern
    l.debug('Exchange format: %s' % (exchangefolder + pattern))
    allversions = glob.glob(exchangefolder + pattern)
    if not allversions:
        error_message = 'Cannot find file from pattern: %s' % (exchangefolder +         pattern)
        l.error(error_message)
        raise exceptions.IOError(error_message)

    return max(allversions, key = lambda x: int(x[-13:-5]))
