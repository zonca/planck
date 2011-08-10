import sys

import private
from Planck import parse_channels
from metadata import DataSelector

import logging as l
import pytoast 

l.basicConfig(level=l.INFO)

def strconv(f):
    """Formatting for the xml"""
    return "%.16g" % f

def Params(dic=None):
    """Creates a Toast ParMap from a python dictionary"""
    params = pytoast.ParMap()
    if not dic is None:
        for k,v in dic.iteritems():
            if isinstance(v, float):
                v = strconv(v)
            else:
                v = str(v)
            params[k] = v
    return params

DEFAULT_OBTMASK = {'LFI':1, 'HFI':1}
DEFAULT_FLAGMASK = {'LFI':255, 'HFI':1}

class ToastConfig(object):
    """Toast configuration class"""

    def __init__(self, odrange, channels, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder=None, fpdb=None, output_xml='toastrun.xml', ahf_folder=None, components='IQU', obtmask=None, flagmask=None, log_level=l.INFO):
        """odrange: list of start and end OD, AHF ODS, i.e. with whole pointing periods as the DPC is using
           channels: one of integer frequency, channel string, list of channel strings
           obtmask and flagmask: default LFI 1,255 HFI 1,1
           """
        l.root.level = log_level
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
        if not exchange_folder is None:
            self.data_selector.config['exchangefolder'] = exchange_folder
        if not ahf_folder is None:
            self.data_selector.config['ahf_folder'] = ahf_folder
        self.data_selector.by_od_range(self.odrange)

        self.wobble = private.WOBBLE
        self.components = components

        self.obtmask = obtmask or DEFAULT_OBTMASK[self.f.inst.name]
        self.flagmask = flagmask or DEFAULT_FLAGMASK[self.f.inst.name]

    def run(self):
        """Call the python-toast bindings to create the xml configuration file"""
        self.conf = pytoast.Run()
          
        sky = self.conf.sky_add ( "sky", "native", pytoast.ParMap() )

        mapset = sky.mapset_add ( '_'.join(['healpix',self.components, self.ordering]), "healpix", 
            Params({
                "path"  : self.outmap,
                "fullsky"  : "TRUE",
                "stokes"  : self.components,
                "order"  : self.ordering,
                "coord"  : self.coord,
                "nside"  : str(self.nside),
                "units"  : "micro-K"
            }))

        tele = self.conf.telescope_add ( "planck", "planck", 
            Params({  
                "wobblepsi2dir":self.wobble["psi2_dir"],
                "wobblepsi2_ref":self.wobble["psi2_ref"],
                "wobblepsi1_ref":self.wobble["psi1_ref"],
            }))

        fp = tele.focalplane_add ( "FP_%s" % self.f.inst.name, "planck_rimo", Params({"path":self.fpdb}) )

        self.add_pointing(tele)
        self.add_observations(tele)
        self.add_channels(tele)
        # write out XML
        self.conf.write ( self.output_xml )

    def add_pointing(self, telescope):
        # Add pointing files
        for i,ahf in enumerate(self.data_selector.get_AHF()):  
          telescope.pointing_add ( "%04d" % i, "planck_ahf", Params({"path": str(ahf[0])}))

    def add_observations(self, telescope):
        """Each observation is a OD as specified in the AHF files"""
        # Add streamset
        strset = telescope.streamset_add ( self.f.inst.name, "native", Params() )
          
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
          
        eff_files = self.data_selector.get_EFF()

        for observation in self.data_selector.get_OBS():  
            params = {"start":observation.start, "stop":observation.stop}
            for i, eff in enumerate(observation.EFF):
                params[ "times%d" % (i+1) ] = eff
            obs = strset.observation_add ( "%04d%s" % (observation.od, observation.tag) , "planck_exchange", Params(params) )

            for pp in observation.PP:
                obs.interval_add( "%05d" % pp.number, "native", Params({"start":pp.start, "stop":pp.stop}) )
          
        # Add streams for real data

        for ch in self.channels:
          rawname = "raw_" + ch.tag
          strm = strset.stream_add ( rawname, "native", Params() )
          
          # Add TODs for this stream
          params = {}
          params[ "flagmask" ] = self.flagmask
          params[ "obtmask" ] = self.obtmask
          for od, eff in zip(self.data_selector.ods, eff_files):
            params[ "hdu" ] = ch.eff_tag
            params[ "path" ] = eff
            strm.tod_add ( "%s_%d" % (ch.tag, od), "planck_exchange", Params(params) )

        # Add white-noise 1 PSD per mission
        for ch in self.channels:
            noise_name = "noise_" + ch.tag
            noise_stream = strset.noise_add ( noise_name, "native", Params() )
            noise_stream.psd_add ( "white", "native", Params({
                "start" : strset.observations()[0].start(),
                "stop" : strset.observations()[-1].stop(),
                "rate": ch.sampling_freq,
                "rms": ch.white_noise 
            }))
            
    def add_channels(self, telescope):
        params = pytoast.ParMap()
        params[ "focalplane" ] = self.conf.telescopes()[0].focalplanes()[0].name()
        for ch in self.channels:
            params[ "detector" ] = ch.tag
            params[ "stream" ] = "%s/raw_%s" % (self.f.inst.name, ch.tag)
            params[ "noise" ] = "%s/noise_%s" % (self.f.inst.name, ch.tag)
            telescope.channel_add ( ch.tag, "native", params )
          
if __name__ == '__main__':

    toast_config = ToastConfig([451, 458], 30, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder='/global/scratch/sd/planck/user/zonca/data/LFI_DX7S_conv/', output_xml='30_break.xml')
    toast_config.run()
