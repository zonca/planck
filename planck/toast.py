import sys

import private
from Planck import parse_channels
from metadata import DataSelector

import pytoast 

def strconv(f):
    return "%.16f" % f

class ToastConfig(object):
    """Toast configuration class"""

    def __init__(self, odrange, channels, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder=None, ahf_folder=None, fpdb=None, output_xml='toastrun.xml'):
        """odrange: list of start and end OD
           channels: one of integer frequency, channel string, list of channel strings"""
        self.odrange = odrange
        self.nside = nside
        self.coord = coord
        self.ordering = ordering
        self.outmap = outmap
        self.channels = parse_channels(channels)
        self.f = self.channels[0].f
        self.output_xml = output_xml
        self.fpdb = fpdb or private.rimo[self.f.inst.name]

        self.data_selector = DataSelector(channels=self.channels)
        self.data_selector.config['exchangefolder'] = 'toast_planck_data/hfi_m3_v41'
        self.data_selector.by_od_range(self.odrange)

    def run(self):
        """Call the python-toast bindings to create the xml configuration file"""
        # create run
        self.conf = pytoast.Run()
          
        # Add sky
        sky = self.conf.sky_add ( "sky", "native", pytoast.ParMap() )

        # Add healpix IQU mapset to sky
        params = pytoast.ParMap()
          
        params[ "path" ] = self.outmap
        params[ "fullsky" ] = "TRUE"
        params[ "stokes" ] = "IQU"
        params[ "order" ] = self.ordering
        params[ "coord" ] = self.coord
        params[ "nside" ] = str(self.nside)
        params[ "units" ] = "micro-K"
          
        mapset = sky.mapset_add ( "healpix", "healpix", params )

        # Add telescope
        tele = self.conf.telescope_add ( "planck", "native", pytoast.ParMap() )

        # Add HFI focalplane
        params = pytoast.ParMap()
        params[ "path" ] = self.fpdb
        fp = tele.focalplane_add ( "FP_%s" % self.f.inst.name, "planck_rimo", params )

        self.add_pointing(tele)
        self.add_observations(tele)
        self.add_channels(tele)
        # write out XML
        self.conf.write ( self.output_xml )

    def add_pointing(self, telescope):
        # Add pointing files
        params = pytoast.ParMap()
        for i,ahf in enumerate(self.data_selector.get_AHF()):  
          params[ "path" ] = str(ahf[0]).replace('att','vel_att')
          telescope.pointing_add ( "%04d" % i, "planck_ahf", params )

    def add_observations(self, telescope):
        """Each observation is a OD as specified in the AHF files"""
        # Add streamset
        params = pytoast.ParMap()
        strset = telescope.streamset_add ( self.f.inst.name, "native", params )
          
        # Add observations

        # For Planck, each observation must cache the timestamps for all days included
        # in the observation, which is very expensive.  For this reason, the number of
        # ODs in a single observation should be the minimum necessary for the size of
        # the interval being considered.  For applications which do not care about noise
        # stationarity or which are only dealing with ring-by-ring quantities, one OD
        # per observation is fine.  If you wish to consider intervals longer than a day
        # as a single stationary interval, then each observation will need to include
        # multiple ODs.

        # For "planck_exchange" format, if no "start" or "stop" times are specified,
        # then the observation will span the time range of the EFF data specified in
        # the "times1", "times2", "times3", etc parameters.
          
        params = pytoast.ParMap()

        eff_files = self.data_selector.get_EFF()

        for observation in self.data_selector.get_OBS():  
            for i, eff in enumerate(observation.EFF):
                params[ "times%d" % (i+1) ] = eff
            params[ "start" ] = strconv(observation.start)
            params[ "stop" ] = strconv(observation.stop)
            obs = strset.observation_add ( "%04d" % observation.od , "planck_exchange", params )

            for pp in observation.PP:
                params = pytoast.ParMap()
                params[ "start" ] = strconv(pp.start)
                params[ "stop" ] = strconv(pp.stop)
                obs.interval_add( "%05d" % pp.number, "native", params )
          
        # Add streams for real data

        for ch in self.channels:
          params = pytoast.ParMap()
          rawname = "raw_" + ch.tag
          strm = strset.stream_add ( rawname, "native", params )
          
          # Add TODs for this stream
          params = pytoast.ParMap()
          params[ "flagmask" ] = "1"
          for od, eff in zip(self.data_selector.ods, eff_files):
            params[ "hdu" ] = ch.eff_tag
            params[ "path" ] = eff
            strm.tod_add ( "%s_%d" % (ch.tag, od), "planck_exchange", params )
            
    def add_channels(self, telescope):
        params = pytoast.ParMap()
        params[ "focalplane" ] = self.conf.telescopes()[0].focalplanes()[0].name()
        for ch in self.channels:
          params[ "detector" ] = ch.tag
          params[ "stream" ] = "%s/raw_%s" % (self.f.inst.name, ch.tag)
          telescope.channel_add ( ch.tag, "native", params )
          
if __name__ == '__main__':

    toast_config = ToastConfig([97, 103], 100, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder='toast_planck_data/hfi_m3_v41', ahf_folder='toast_planck_data/AHF_v1', output_xml='toastrun.xml')
    toast_config.run()
