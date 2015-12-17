import sys
import numpy as np
import copy
import os
import exceptions
import glob
import sqlite3

import private
from Planck import parse_channels, EXCLUDED_CH
from metadata import DataSelector
from collections import defaultdict

import logging as l
from pytoast.core import Run, ParMap

from astropy.io import fits as pyfits

l.basicConfig(level=l.INFO)


class PPBoundaries:
    def __init__(self, freq, dbfile):
        """Load the start and stop OBT timestamps extracted from exchange format files"""
        self.dbfile = dbfile

        conn = sqlite3.connect( dbfile )
        c = conn.cursor()

        self.ppf = []

        tabname = 'ring_times_hfi'

        if ( freq == 30 ):
            tabname = 'ring_times_lfi30'
        if ( freq == 44 ):
            tabname = 'ring_times_lfi44'
        if ( freq == 70 ):
            tabname = 'ring_times_lfi70'

        execstr = "select id, pointID_unique, start, stop from %s" % ( tabname )
        query = c.execute( execstr )
        for id, pointID_unique, start, stop in query:
            self.ppf += [ ( start, stop ) ]

        c.close()

    def get(self, PID):
        # LFI PID 3 is row 1
        filerow = PID - 2
        return self.ppf[filerow]

def get_eff_od(file_path):
    return int(os.path.basename(file_path).split('-')[1])

def strconv(f):
    """Formatting for the xml"""
    return "%.16g" % f

def Params(dic=None):
    """Creates a Toast ParMap from a python dictionary"""
    params = ParMap()
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

    def __init__(self, odrange=None, channels=None, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder=None, fpdb=None, output_xml='toastrun.xml', ahf_folder=None, components='IQU', obtmask=None, flagmask=None, log_level=l.INFO, remote_exchange_folder=None, remote_ahf_folder=None, calibration_file=None, dipole_removal=False, noise_tod=False, noise_tod_weight=None, efftype=None, flag_HFI_bad_rings=None, include_preFLS=False, ptcorfile=None, include_repointings=False, psd=None, deaberrate=True, extend_857=False, no_wobble=False, eff_is_for_flags=False, exchange_weights=None, beamsky=None, beamsky_weight=None, interp_order=5, horn_noise_tod=None, horn_noise_weight=None, horn_noise_psd=None, observation_is_interval=False, lfi_ring_range=None, hfi_ring_range=None, wobble_high=False):
        """TOAST configuration:

            odrange: list of start and end OD, AHF ODS, i.e. with whole pointing periods as the DPC is using
            lfi_ring_range: first and last LFI pointing ID to include
            hfi_ring_range: first and last HFI ring number to include
            channels: one of integer frequency, channel string, list of channel strings
            obtmask and flagmask: default LFI 1,255 HFI 1,1
            exchange_folder: either a string or a container or strings listing locations for Exchange Format data
            exchange_weights: list of scaling factors to apply to Exchange data prior to coadding
            eff_is_for_flags : only use Exchange data for flags
            remote_exchange_folder, remote_ahf_folder: they allow to run toast.py in one environment using ahf_folder and exchange_folder and then replace the path with the remote folders
            calibration_file: path to a fits calibration file, with first extension OBT, then one extension per channel with the calibration factors
            dipole_removal: dipole removal is performed ONLY if calibration is specified
            noise_tod: Add simulated noise TODs
            noise_tod_weight: scaling factor to apply to noise_tod
            flag_HFI_bad_rings: if None, flagged just for HFI. If a valid file, use that as input.
            include_preFLS : if None, True for LFI
            include_repointings : Construct intervals with the 4 minute repointing maneuvers: 10% dataset increase, continuous time stamps
            psd : templated name of an ASCII PSD files. Tag CHANNEL will be replaced with the appropriate channel identifier.
            horn_noise_tod : False
            horn_noise_weight : None
            horn_noise_psd : templated name of an ASCII PSD files. Tag HORN will be replaced with the appropriate identifier.
            deaberrate : Correct pointing for aberration
            extend_857 : Whether or not to include the RTS bolometer, 857-4 in processing
            no_wobble : Disable all flavors of wobble angle correction
            wobble_high : Use 8Hz wobble correction rather than per ring
            beamsky : templated name of the beamsky files for OTF sky convolution. Tag CHANNEL will be replaced with the appropriate channel identifier.
            beamsky_weight : scaling factor to apply to the beamsky
            interp_order : beamsky interpolation order defines number of cartesian pixels to interpolate over
            observation_is_interval : If True, do not split operational days into pointing period intervals

            additional configuration options are available modifying:
            .config
            and Data Selector configuration:
            .data_selector.config
            dictionaries before running .run()
        """
        l.root.level = log_level
        self.extend_857 = extend_857
        if self.extend_857:
            if '857-4' in EXCLUDED_CH: EXCLUDED_CH.remove('857-4')
        if sum([odrange!=None, lfi_ring_range!=None, hfi_ring_range!=None]) != 1:
            raise Exception('Must specify exactly one type of data span: OD, LFI PID or HFI ring')
        self.odrange = odrange
        self.lfi_ring_range = lfi_ring_range
        self.hfi_ring_range = hfi_ring_range
        self.nside = nside
        self.coord = coord
        self.ordering = ordering
        self.outmap = outmap
        if channels == None:
            raise Exception('Must define which channels to include')
        self.channels = parse_channels(channels)
        self.f = self.channels[0].f
        self.output_xml = output_xml
        self.ptcorfile = ptcorfile
        if no_wobble and wobble_high:
            raise Exception('no_wobble and wobble_high are mutually exclusive')
        self.no_wobble = no_wobble
        self.wobble_high = wobble_high
        self.include_repointings = include_repointings
        self.deaberrate = deaberrate
        self.noise_tod = noise_tod
        self.noise_tod_weight = noise_tod_weight
        self.psd = psd
        self.horn_noise_tod = horn_noise_tod
        self.horn_noise_weight = horn_noise_weight
        self.horn_noise_psd = horn_noise_psd
        if self.horn_noise_tod and not self.horn_noise_psd:
            raise Exception('Must specify horn_noise_psd template name when enabling horn_noise')
        self.fpdb = fpdb or private.rimo[self.f.inst.name]
        self.beamsky = beamsky
        self.beamsky_weight = beamsky_weight
        self.interp_order = interp_order
        self.observation_is_interval = observation_is_interval
        self.rngorder = {
            'LFI18M' : 0,
            'LFI18S' : 1,
            'LFI19M' : 2,
            'LFI19S' : 3,
            'LFI20M' : 4,
            'LFI20S' : 5,
            'LFI21M' : 6,
            'LFI21S' : 7,
            'LFI22M' : 8,
            'LFI22S' : 9,
            'LFI23M' : 10,
            'LFI23S' : 11,
            'LFI24M' : 12,
            'LFI24S' : 13,
            'LFI25M' : 14,
            'LFI25S' : 15,
            'LFI26M' : 16,
            'LFI26S' : 17,
            'LFI27M' : 18,
            'LFI27S' : 19,
            'LFI28M' : 20,
            'LFI28S' : 21,
            '100-1a' : 22,
            '100-1b' : 23,
            '100-2a' : 24,
            '100-2b' : 25,
            '100-3a' : 26,
            '100-3b' : 27,
            '100-4a' : 28,
            '100-4b' : 29,
            '143-1a' : 30,
            '143-1b' : 31,
            '143-2a' : 32,
            '143-2b' : 33,
            '143-3a' : 34,
            '143-3b' : 35,
            '143-4a' : 36,
            '143-4b' : 37,
            '143-5'  : 38,
            '143-6'  : 39,
            '143-7'  : 40,
            '143-8'  : 41,
            '217-5a' : 42,
            '217-5b' : 43,
            '217-6a' : 44,
            '217-6b' : 45,
            '217-7a' : 46,
            '217-7b' : 47,
            '217-8a' : 48,
            '217-8b' : 49,
            '217-1'  : 50,
            '217-2'  : 51,
            '217-3'  : 52,
            '217-4'  : 53,
            '353-3a' : 54,
            '353-3b' : 55,
            '353-4a' : 56,
            '353-4b' : 57,
            '353-5a' : 58,
            '353-5b' : 59,
            '353-6a' : 60,
            '353-6b' : 61,
            '353-1'  : 62,
            '353-2'  : 63,
            '353-7'  : 64,
            '353-8'  : 65,
            '545-1'  : 66,
            '545-2'  : 67,
            '545-3'  : 68,
            '545-4'  : 69,
            '857-1'  : 70,
            '857-2'  : 71,
            '857-3'  : 72,
            '857-4'  : 73,
            'LFI18' : 74,
            'LFI19' : 75,
            'LFI20' : 76,
            'LFI21' : 77,
            'LFI22' : 78,
            'LFI23' : 79,
            'LFI24' : 80,
            'LFI25' : 81,
            'LFI26' : 82,
            'LFI27' : 83,
            'LFI28' : 84,
            '100-1' : 85,
            '100-2' : 86,
            '100-3' : 87,
            '100-4' : 88,
            '143-1' : 89,
            '143-2' : 90,
            '143-3' : 91,
            '143-4' : 92,
            '217-5' : 93,
            '217-6' : 94,
            '217-7' : 95,
            '217-8' : 96,
            '353-3' : 97,
            '353-4' : 98,
            '353-5' : 99,
            '353-6' : 100,
            }

        self.config = {}
        if self.f.inst.name == 'LFI':
            self.config['pairflags'] = True
        else:
            self.config['pairflags'] = False

        if efftype is None:
            efftype ='R'
            if self.f.inst.name == 'LFI' and (not calibration_file is None):
                efftype ='C'
        self.data_selector = DataSelector(channels=self.channels, efftype=efftype, include_preFLS=include_preFLS)

        self.exchange_weights = exchange_weights

        if remote_exchange_folder:
            if isinstance(remote_exchange_folder, str):
                remote_exchange_folder = [remote_exchange_folder]
            if remote_exchange_folder[-1] != '/':
                remote_exchange_folder += '/'

        self.remote_exchange_folder = remote_exchange_folder
        if remote_ahf_folder:
            if remote_ahf_folder[-1] != '/':
                remote_ahf_folder += '/'
        self.remote_ahf_folder = remote_ahf_folder
        if exchange_folder:
            if isinstance(exchange_folder, str):
                exchange_folder = [exchange_folder]
            self.exchange_folder = exchange_folder

            self.data_selector.config[ 'exchangefolder' ] = exchange_folder[0]
            
        if ahf_folder:
            self.data_selector.config['ahf_folder'] = ahf_folder

        if self.odrange:
            self.data_selector.by_od_range(self.odrange)
        elif self.lfi_ring_range:
            self.data_selector.by_lfi_rings(self.lfi_ring_range)
        elif self.hfi_ring_range:
            self.data_selector.by_hfi_rings(self.hfi_ring_range)
        else:
            raise Exception('Must specify one type of data span')
        

        self.wobble = private.WOBBLE
        self.components = components

        self.obtmask = obtmask or DEFAULT_OBTMASK[self.f.inst.name]
        self.flagmask = flagmask or DEFAULT_FLAGMASK[self.f.inst.name]

        self.calibration_file = calibration_file
        self.dipole_removal = dipole_removal

        if flag_HFI_bad_rings is None:
            if self.f.inst.name == 'HFI':
                flag_HFI_bad_rings = True
            else:
                flag_HFI_bad_rings = False
        if flag_HFI_bad_rings:
            if ( os.path.isfile(str(flag_HFI_bad_rings)) ):
                self.bad_rings = flag_HFI_bad_rings
            else:
                self.bad_rings = private.HFI_badrings
        else:
            self.bad_rings = None

        self.eff_is_for_flags = eff_is_for_flags

        self.fptab = None

    def parse_fpdb(self, channel, key):
        if self.fptab == None:
            fptab = pyfits.getdata(self.fpdb, 1)
        ind = (fptab.field('DETECTOR') == channel).ravel()
        if np.sum(ind) != 1:
            raise Exception('Error matching {} in {}: {} matches.'.format(channel, self.fpdb, np.sum(ind)))
        return fptab.field(key)[ind][0]

    def run(self, write=True):
        """Call the python-toast bindings to create the xml configuration file"""
        self.conf = Run()

        if self.noise_tod or self.horn_noise_tod:
            self.conf.variable_add ( "rngbase", "native", Params({"default":"0"}) )
          
        sky = self.conf.sky_add ( "sky", "native", ParMap() )

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

        #if self.f.inst.name == 'LFI':
        #    wobble_offset = 0;
        #else:
        #    wobble_offset = self.wobble["psi2_offset"]

        if self.no_wobble:
            teleparams = {
                "wobble_ahf_high":"FALSE",
                "wobble_ahf_obs":"FALSE",
                "wobblepsi2dir":""
                }
        else:
            if self.wobble_high:
                wobble_obs = 'FALSE'
                wobble_high = 'TRUE'
            else:
                wobble_obs = 'TRUE'
                wobble_high = 'FALSE'
            
            teleparams = {  
                #"wobblepsi2dir":self.wobble["psi2_dir"],
                "wobblepsi2_ref":self.wobble["psi2_ref"],
                "wobblepsi1_ref":self.wobble["psi1_ref"],
                "wobble_ahf_obs":wobble_obs,
                "wobble_ahf_high":wobble_high
                #"wobblepsi2_offset":wobble_offset
                }

        if self.ptcorfile:
            teleparams['ptcorfile'] = self.ptcorfile

        if self.deaberrate != None:
            if self.deaberrate:
                teleparams['deaberrate'] = 'TRUE'
            else:
                teleparams['deaberrate'] = 'FALSE'

        tele = self.conf.telescope_add ( "planck", "planck", 
            Params(teleparams))

        fp = tele.focalplane_add ( "FP_%s" % self.f.inst.name, "planck_rimo", Params({"path":self.fpdb}) )

        self.add_pointing(tele)
        self.add_observations(tele)
        self.add_streams()
        self.add_eff_tods()
        self.add_noise()
        self.add_channels(tele)
        if write:
            self.write()


    def write(self):
        # write out XML
        self.conf.write ( self.output_xml )


    def add_pointing(self, telescope):
        # Add pointing files
        for i,ahf in enumerate(self.data_selector.get_AHF()):  
            path = str(ahf[0])
            if self.remote_ahf_folder:
                path = path.replace(self.data_selector.config['ahf_folder'], self.remote_ahf_folder)
            telescope.pointing_add ( "%04d" % i, "planck_ahf", Params({"path": path}))


    @property
    def observations(self):
        try:
            return self._observations
        except:
            self._observations = self.data_selector.get_OBS()
            return self.observations


    def add_observations(self, telescope):
        """Each observation is an OD as specified in the AHF files"""
        # Add streamset
        self.strset = telescope.streamset_add ( self.f.inst.name, "native", Params() )
          
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

        broken_od = defaultdict(None)

        # Observations are same for all datasets but now there may be
        # several exchange folders and raw streams to add. Therefore
        # tod_name_list and tod_par_list are lists of dictionaries.
        
        self.tod_name_list = []
        self.tod_par_list = []
        for i in range(len(self.exchange_folder)):
            self.tod_name_list.append( defaultdict(list) )
            self.tod_par_list.append( defaultdict(list) )
            
        # Add observations

        for iobs, observation in enumerate(self.observations):
            params = {"start":observation.start, "stop":observation.stop}
            for i, eff in enumerate(observation.EFF):
                if self.remote_exchange_folder:
                    params[ "times%d" % (i+1) ] = eff.replace(self.exchange_folder[0], self.remote_exchange_folder[0])
                else:
                    params[ "times%d" % (i+1) ] = eff
            obs = self.strset.observation_add ( "%04d%s" % (observation.od, observation.tag) , "planck_exchange", Params(params) )

            if self.observation_is_interval:
                # intervals are observations
                obs.interval_add( '{:04}'.format(iobs), "native", Params({"start":observation.start, "stop":observation.stop}) )
            else:
                pointing_periods = observation.PP
                if self.include_repointings:
                    # First interval begins with the operational day, last one ends with it. New interval begins when the old stable pointing ends.
                    # First and last included samples remain the same.
                    for ipp, pp in enumerate(pointing_periods):
                        if ipp == 0:
                            # Remove the comments in this section to include period between OD and ring start
                            #if iobs == 0:
                            obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":pp.start, "stop":pointing_periods[ipp+1].start}) )
                            #else:
                            #    obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":observation.start, "stop":pointing_periods[ipp+1].start}) )
                        elif ipp == len(pointing_periods) - 1:
                            if iobs == len(self.observations) - 1:
                                obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":pp.start, "stop":pp.stop}) )
                            else:
                                obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":pp.start, "stop":observation.stop}) )
                        else:
                            obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":pp.start, "stop":pointing_periods[ipp+1].start}) )
                else:
                    # Intervals are the stable pointing periods
                    for pp in pointing_periods:
                        obs.interval_add( "%05d-%d" % (pp.number, pp.splitnumber), "native", Params({"start":pp.start, "stop":pp.stop}) )

            for ch in self.channels:
                print("Observation %d%s, EFF ODs:%s" % (observation.od, observation.tag, str(map(get_eff_od, observation.EFF))))
                for ix, fx in enumerate(self.exchange_folder):
                    for i, file_path in enumerate(observation.EFF):
                        eff_od = get_eff_od(file_path)
                        # Add TODs for this stream
                        params = {}
                        params[ "flagmask" ] = self.flagmask
                        params[ "obtmask" ] = self.obtmask
                        params[ "hdu" ] = ch.eff_tag
                        if self.config['pairflags']:
                            params[ "pairflags" ] = 'TRUE'
                        if self.remote_exchange_folder:
                            params[ "path" ] = file_path.replace(self.exchange_folder[0], self.remote_exchange_folder[ix])
                        else:
                            params[ "path" ] = file_path.replace(self.exchange_folder[0], self.exchange_folder[ix])
                        tag = ''
                        if i==(len(observation.EFF)-1) and not observation.break_startrow is None:
                            params['rows'] = observation.break_startrow + 1
                            tag = 'a'
                            broken_od[ch.tag] = eff_od
                        if not observation.break_stoprow is None and broken_od[ch.tag]==eff_od:
                            params['startrow'] = observation.break_stoprow
                            tag = 'b'
                            broken_od[ch.tag] = None
                        name = "%s_%d%s" % (ch.tag, eff_od, tag)
                        if name not in self.tod_name_list[ix][ch.tag]:
                            print('add ' + name)
                            self.tod_name_list[ix][ch.tag].append(name)
                            self.tod_par_list[ix][ch.tag].append(params)
                        else:
                            print("skip " + name)


    def add_streams(self):
        # Add streams for data components
        basename = "@rngbase@"

        self.strm = {}

        for ch in self.channels:
            stack_elements = []

            # add simulated noise stream
            if self.noise_tod:
                rngstream = self.rngorder[ ch.tag ] * 100000

                noisename = "/planck/" + self.f.inst.name + "/noise_" + ch.tag
                self.strm["simnoise_" + ch.tag] = self.strset.stream_add( "simnoise_" + ch.tag, "native", Params( ) )
                suffix = ''
                if self.noise_tod_weight:
                    suffix += ',PUSH:$' + strconv(self.noise_tod_weight) + ',MUL'
                stack_elements.append( "PUSH:simnoise_" + ch.tag + suffix)
                if self.observation_is_interval:
                    # one noise tod per observation.
                    #
                    # FIXME: this code will not work.  the start and stop times
                    #    must be the *detector* start and stop times inside
                    #    the observation boundaries (which are set by the AHF).
                    #    so you must iterate over self.pp_boundaries and find
                    #    the first and last valid detector times inside this observation.
                    #
                    for iobs, observation in enumerate(self.observations):
                        self.strm["simnoise_" + ch.tag].tod_add ( "nse_%s_%05d" % (ch.tag, iobs), "sim_noise", Params({
                               "noise" : noisename,
                               "base" : basename,
                               "start" : observation.start,
                               "stop" : observation.stop,
                               "offset" : rngstream + iobs
                        }))                        
                else:
                    # one noise tod per pointing period
                    self.pp_boundaries = PPBoundaries(self.f.freq, self.data_selector.config['database'])
                    for row, pp_boundaries in enumerate(self.pp_boundaries.ppf):
                        if pp_boundaries[1] < self.observations[0].start or pp_boundaries[0] > self.observations[-1].stop:
                            continue
                        self.strm["simnoise_" + ch.tag].tod_add ( "nse_%s_%05d" % (ch.tag, row), "sim_noise", Params({
                               "noise" : noisename,
                               "base" : basename,
                               "start" : pp_boundaries[0],
                               "stop" : pp_boundaries[1],
                               "offset" : rngstream + row
                        }))

            # add simulated noise stream common to each horn
            if self.horn_noise_tod and ch.tag[-1] in 'MSab':
                horn = ch.tag[:-1]
                rngstream = self.rngorder[ horn ] * 100000

                noisename = "/planck/" + self.f.inst.name + "/noise_" + horn
                self.strm["simnoise_" + horn] = self.strset.stream_add( "simnoise_" + horn, "native", Params( ) )
                suffix = ''
                if self.horn_noise_weight:
                    suffix += ',PUSH:$' + strconv(self.horn_noise_weight) + ',MUL'
                if len(stack_elements) != 0:
                    suffix += ',ADD'
                stack_elements.append( "PUSH:simnoise_" + horn + suffix)
                if self.observation_is_interval:
                    # one noise tod per observation
                    for iobs, observation in enumerate(self.observations):
                        self.strm["simnoise_" + horn].tod_add ( "nse_%s_%05d" % (horn, iobs), "sim_noise", Params({
                               "noise" : noisename,
                               "base" : basename,
                               "start" : observation.start,
                               "stop" : observation.stop,
                               "offset" : rngstream + iobs
                        }))                        
                else:
                    # one noise tod per pointing period
                    self.pp_boundaries = PPBoundaries(self.f.freq)
                    for row, pp_boundaries in enumerate(self.pp_boundaries.ppf):
                        if pp_boundaries[1] < self.observations[0].start or pp_boundaries[0] > self.observations[-1].stop:
                            continue
                        self.strm["simnoise_" + horn].tod_add ( "nse_%s_%05d" % (horn, row), "sim_noise", Params({
                               "noise" : noisename,
                               "base" : basename,
                               "start" : pp_boundaries[0],
                               "stop" : pp_boundaries[1],
                               "offset" : rngstream + row
                        }))

            # add the beam sky
            if self.beamsky:
                epsilon = self.parse_fpdb(ch.tag, 'epsilon')
                #beamskyname = '/planck/' + self.f.inst.name + '/' + ch.tag
                beamskyname = '/planck/' + ch.tag
                self.strm["beamsky_" + ch.tag] = self.strset.stream_add( "beamsky_" + ch.tag, "planck_beamsky", Params({
                    'path' : self.beamsky.replace('CHANNEL', ch.tag),
                    'epsilon' : str(epsilon),
                    'interp_order' : str(self.interp_order),
                    'coord' : self.coord,
                    'channel' : beamskyname
                    }) )
                suffix = ''
                if self.beamsky_weight:
                    suffix += ',PUSH:$' + strconv(self.beamsky_weight) + ',MUL'
                if len(stack_elements) != 0:
                    suffix += ',ADD'
                stack_elements.append( "PUSH:beamsky_" + ch.tag + suffix)

            # add real data streams, either for flags or data plus flags
            for i, x in enumerate(self.exchange_folder):
                self.strm[ 'raw{}_{}'.format(i, ch.tag) ] = self.strset.stream_add( "raw{}_{}".format(i, ch.tag), "native", Params() )
                if self.eff_is_for_flags:
                    suffix = ',FLG'
                else:
                    suffix = ''
                    if self.exchange_weights:
                        suffix = ',PUSH:$' + strconv(self.exchange_weights[i]) + ',MUL'
                    if len(stack_elements) != 0:
                        suffix += ',ADD'
                stack_elements.append( "PUSH:raw{}_{}{}".format(i, ch.tag, suffix) )

            # add calibration stream
            if (not self.calibration_file is None):
                self.strm["cal_" + ch.tag] = self.strset.stream_add( "cal_" + ch.tag, "planck_cal", Params( {"hdu":ch.tag, "path":self.calibration_file } ) )
                stack_elements.append("PUSH:cal_" + ch.tag + ",MUL")

            # dipole subtract
            if self.dipole_removal:
                self.strm["dipole_" + ch.tag] = self.strset.stream_add( "dipole_" + ch.tag, "dipole", Params( {"channel":ch.tag, "coord":"E"} ) )
                stack_elements.append("PUSH:dipole_" + ch.tag + ",SUB")

            # bad rings
            if not self.bad_rings is None:
                self.strm["bad_" + ch.tag] = self.strset.stream_add ( "bad_" + ch.tag, "planck_bad", Params({'detector':ch.tag, 'path':self.bad_rings}) )
                stack_elements.append("PUSH:bad_" + ch.tag + ",FLG")

            # stack
            expr = ','.join([el for el in stack_elements])
            self.strm["stack_" + ch.tag] = self.strset.stream_add ( "stack_" + ch.tag, "stack", Params( {"expr":expr} ) )


    def add_eff_tods(self):
        """Add TOD files already included in tod_name_list and tod_par_list to the streamset"""
        for ix in range(len(self.exchange_folder)):
            # remove duplicate files on breaks
            for ch in self.channels:
                for name in [n for n in self.tod_name_list[ix][ch.tag] if n.endswith('a')]:
                    try:
                        delindex = self.tod_name_list[ix][ch.tag].index(name.rstrip('a'))
                        print('Removing %s because of breaks' % self.tod_name_list[ix][ch.tag][delindex])
                        del self.tod_name_list[ix][ch.tag][delindex]
                        del self.tod_par_list[ix][ch.tag][delindex]
                    except exceptions.ValueError:
                        pass

            # add EFF to stream
            for ch in self.channels:
                for name, par in zip(self.tod_name_list[ix][ch.tag], self.tod_par_list[ix][ch.tag]):
                    self.strm[ "raw{}_{}".format(ix, ch.tag) ].tod_add ( name, "planck_exchange", Params(par) )


    def add_noise(self):
        # Add RIMO noise model (LFI) or mission average PSD (HFI)

        for ch in self.channels:
            noise = self.strset.noise_add ( "noise_" + ch.tag, "native", Params() )

            # add PSD
            if self.psd:                
                psdname = self.psd.replace('CHANNEL', ch.tag)
                noise.psd_add ( "psd", "ascii", Params({
                    "start" : self.strset.observations()[0].start(),
                    "stop" : self.strset.observations()[-1].stop(),
                    "path": psdname
                    }))
            else:
                noise.psd_add ( "psd", "planck_rimo", Params({
                    "start" : self.strset.observations()[0].start(),
                    "stop" : self.strset.observations()[-1].stop(),
                    "path": self.fpdb,
                    "detector": ch.tag
                    }))
            # add horn noise psd
            if self.horn_noise_tod and ch.tag[-1] in 'aM':
                horn = ch.tag[:-1]
                noise = self.strset.noise_add ( "noise_" + horn, "native", Params() )
                psdname = self.horn_noise_psd.replace('HORN', horn)
                noise.psd_add ( "psd", "ascii", Params({
                    "start" : self.strset.observations()[0].start(),
                    "stop" : self.strset.observations()[-1].stop(),
                    "path": psdname
                    }))

            
    def add_channels(self, telescope):
        params = ParMap()
        params[ "focalplane" ] = self.conf.telescopes()[0].focalplanes()[0].name()
        for ch in self.channels:
            params[ "detector" ] = ch.tag
            params[ "stream" ] = "/planck/%s/stack_%s" % (self.f.inst.name, ch.tag)
            params[ "noise" ] = "/planck/%s/noise_%s" % (self.f.inst.name, ch.tag)
            telescope.channel_add ( ch.tag, "native", params )


if __name__ == '__main__':

    # Default LFI run

    conf = ToastConfig([740, 760], 30, nside=1024, ordering='RING', coord='E', efftype='R', output_xml='test_30_default.xml')
    conf.run()

#    # LFI run with noise simulation and real data flags
#    
    #conf = ToastConfig([450, 460], 30, nside=1024, ordering='RING', coord='E', efftype='C', output_xml='test_30_simnoise.xml', noise_tod=True)
    #conf.run()
#
#    # LFI real data run with calibration and dipole subtraction
#
#    conf = ToastConfig([97, 101], 30, nside=1024, ordering='RING', coord='E', efftype='C', output_xml='test_30_dical.xml', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/lfi/369S/C030-0000-369S-20110713.fits', dipole_removal=True)
#    conf.run()
#
#    # LFI noise simulation with real flags, calibration and dipole subtraction
#
#    conf = ToastConfig([97, 101], 30, nside=1024, ordering='RING', coord='E', efftype='C', output_xml='test_30_simdical.xml', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/lfi/369S/C030-0000-369S-20110713.fits', dipole_removal=True, noise_tod=True)
#    conf.run()
#
#    # HFI default run
#
#    conf = ToastConfig([97, 101], 100, nside=2048, ordering='NEST', coord='G', output_xml='test_100_default.xml')
#    conf.run()
#
#    # HFI noise simulation
#
#    conf = ToastConfig([97, 101], 100, nside=2048, ordering='NEST', coord='G', output_xml='test_100_simnoise.xml', noise_tod=True)
#    conf.run()
#
#    # HFI real data run with calibration and dipole subtraction
#
#    conf = ToastConfig([97, 101], 100, nside=2048, ordering='NEST', coord='G', output_xml='test_100_dical.xml', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/hfi/variable_gains_100GHz.fits', dipole_removal=True)
#    conf.run()
#
#    # HFI noise simulation with real flags, calibration and dipole subtraction
#
#    conf = ToastConfig([97, 101], 100, nside=2048, ordering='NEST', coord='G', output_xml='test_100_simdical.xml', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/hfi/variable_gains_100GHz.fits', dipole_removal=True, noise_tod=True)
#    conf.run()
#
