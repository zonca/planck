import sys
import numpy as np
import copy
import os
import exceptions
import glob

import private
from Planck import parse_channels
from metadata import DataSelector
from collections import defaultdict

import logging as l
import pytoast 

l.basicConfig(level=l.INFO)

class PPBoundaries:
    def __init__(self, freq):
        """Load the start and stop OBT timestamps extracted from exchange format files"""
        ppfile = '/global/scratch/sd/colliera/LS/%d/ts.txt' % freq
        print('Loading ' + ppfile)
        self.ppf = np.loadtxt(ppfile)

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

    def __init__(self, odrange, channels, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder=None, fpdb=None, output_xml='toastrun.xml', ahf_folder=None, components='IQU', obtmask=None, flagmask=None, log_level=l.INFO, remote_exchange_folder=None, remote_ahf_folder=None, calibration_file=None, dipole_removal=True, efftype=None, flag_HFI_bad_rings=None, include_preFLS=None):
        """TOAST configuration:

            odrange: list of start and end OD, AHF ODS, i.e. with whole pointing periods as the DPC is using
            channels: one of integer frequency, channel string, list of channel strings
            obtmask and flagmask: default LFI 1,255 HFI 1,1
            remote_exchange_folder, remote_ahf_folder: they allow to run toast.py in one environment using ahf_folder and exchange_folder and then replace the path with the remote folders
            calibration_file: path to a fits calibration file, with first extension OBT, then one extension per channel with the calibration factors
            dipole_removal: dipole removal is performed ONLY if calibration is specified
            flag_HFI_bad_rings: if None, flagged just for HFI
            include_preFLS : if None, True for LFI

            additional configuration options are available modifying:
            .config
            and Data Selector configuration:
            .ds.config
            dictionaries before running .run()
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

        self.config = {}
        if self.f.inst.name == 'LFI':
            self.config['pairflags'] = True
            if include_preFLS is None:
                include_preFLS = True
        else:
            self.config['pairflags'] = False
            if include_preFLS is None:
                include_preFLS = False

        if efftype is None:
            efftype ='R'
            if self.f.inst.name == 'LFI' and (not calibration_file is None):
                efftype ='C'
        self.data_selector = DataSelector(channels=self.channels, efftype=efftype, include_preFLS=include_preFLS)
        if remote_exchange_folder:
            if remote_exchange_folder[-1] != '/':
                remote_exchange_folder += '/'
        self.remote_exchange_folder = remote_exchange_folder
        if remote_ahf_folder:
            if remote_ahf_folder[-1] != '/':
                remote_ahf_folder += '/'
        self.remote_ahf_folder = remote_ahf_folder
        if not exchange_folder is None:
            self.data_selector.config['exchangefolder'] = exchange_folder
        if not ahf_folder is None:
            self.data_selector.config['ahf_folder'] = ahf_folder
        self.data_selector.by_od_range(self.odrange)

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
            self.bad_rings = private.HFI_badrings
        else:
            self.bad_rings = None

    def run(self, write=True):
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

        if self.f.inst.name == 'LFI':
            wobble_offset = 0;
        else:
            wobble_offset = self.wobble["psi2_offset"]

        tele = self.conf.telescope_add ( "planck", "planck", 
            Params({  
                "wobblepsi2dir":self.wobble["psi2_dir"],
                "wobblepsi2_ref":self.wobble["psi2_ref"],
                "wobblepsi1_ref":self.wobble["psi1_ref"],
                "wobblepsi2_offset":wobble_offset
            }))

        fp = tele.focalplane_add ( "FP_%s" % self.f.inst.name, "planck_rimo", Params({"path":self.fpdb}) )

        self.add_pointing(tele)
        self.add_observations(tele)
        self.add_tods()
        self.add_noise()
        if not self.bad_rings is None:
            self.add_bad_rings()
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
        """Each observation is a OD as specified in the AHF files"""
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

        # Add streams for real data

        strm = {}

        for ch in self.channels:
          rawname = "raw_" + ch.tag
          strm[ch.tag] = self.strset.stream_add ( rawname, "native", Params() )
        if (not self.calibration_file is None) and self.dipole_removal:
            for ch in self.channels:
                strm["dipole_" + ch.tag] = self.strset.stream_add( "dipole_" + ch.tag, "dipole", Params( {"channel":ch.tag, "coord":"E"} ) )

        if not self.calibration_file is None:
            for ch in self.channels:
                strm["cal_" + ch.tag] = self.strset.stream_add( "cal_" + ch.tag, "planck_cal", Params( {"hdu":ch.tag, "path":self.calibration_file } ) )

            #stack
            for ch in self.channels:
                stack_elements = ["raw_" + ch.tag]
                if (not self.calibration_file is None):
                    stack_elements.append("cal_" + ch.tag + ",MUL")
                    if self.dipole_removal:
                        stack_elements.append("dipole_" + ch.tag + ",SUB")
                expr = ','.join(['PUSH:' + el for el in stack_elements])
                if not self.bad_rings is None:
                    expr += ',PUSH:bad_%s,ADD' % ch.tag
                calname = "cal_" + ch.tag
                strm["stack_" + ch.tag] = self.strset.stream_add ( "stack_" + ch.tag, "stack", Params( {"expr":expr} ) )

        broken_od = defaultdict(None)
        # Add observations
        self.tod_name_list = defaultdict(list) 
        self.tod_par_list = defaultdict(list) 
        for observation in self.observations:  
            params = {"start":observation.start, "stop":observation.stop}
            for i, eff in enumerate(observation.EFF):
                if self.remote_exchange_folder:
                    params[ "times%d" % (i+1) ] = eff.replace(self.data_selector.config['exchangefolder'], self.remote_exchange_folder)
                else:
                    params[ "times%d" % (i+1) ] = eff
            obs = self.strset.observation_add ( "%04d%s" % (observation.od, observation.tag) , "planck_exchange", Params(params) )

            pointing_periods = observation.PP
            for pp in pointing_periods:
                obs.interval_add( "%05d" % pp.number, "native", Params({"start":pp.start, "stop":pp.stop}) )

            for ch in self.channels:
              print("Observation %d%s, EFF ODs:%s" % (observation.od, observation.tag, str(map(get_eff_od, observation.EFF))))
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
                      params[ "path" ] = file_path.replace(self.data_selector.config['exchangefolder'], self.remote_exchange_folder)
                  else:
                      params[ "path" ] = file_path
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
                  if name not in self.tod_name_list[ch.tag]:
                      print('add ' + name)
                      self.tod_name_list[ch.tag].append(name)
                      self.tod_par_list[ch.tag].append(params)
                  else:
                      print("skip " + name)

    def add_tods(self):
        """Add TOD files already included in tod_name_list and tod_par_list to the streamset"""
        # remove duplicate files on breaks
        for ch in self.channels:
            for name in [n for n in self.tod_name_list[ch.tag] if n.endswith('a')]:
                try:
                    delindex = self.tod_name_list[ch.tag].index(name.rstrip('a'))
                    print('Removing %s because of breaks' % self.tod_name_list[ch.tag][delindex])
                    del self.tod_name_list[ch.tag][delindex]
                    del self.tod_par_list[ch.tag][delindex]
                except exceptions.ValueError:
                    pass
    
        # add EFF to stream
        for ch, stream in zip(self.channels, self.strset.streams()):
             for name, par in zip(self.tod_name_list[ch.tag], self.tod_par_list[ch.tag]):
                  stream.tod_add ( name, "planck_exchange", Params(par) ) 

    def add_noise(self):
        # Add white-noise 1 PSD per mission
        for ch in self.channels:
            noise_name = "noise_" + ch.tag
            noise_stream = self.strset.noise_add ( noise_name, "native", Params() )
            noise_stream.psd_add ( "white", "native", Params({
                "start" : self.strset.observations()[0].start(),
                "stop" : self.strset.observations()[-1].stop(),
                "rate": ch.sampling_freq,
                "rms": ch.white_noise 
            }))

    def add_bad_rings(self):
        for ch in self.channels:
           self.strset.stream_add ( '_'.join(['bad', ch.tag]), "planck_bad", Params({'detector':ch.tag, 'path':self.bad_rings}) )
            
    def add_channels(self, telescope):
        params = pytoast.ParMap()
        params[ "focalplane" ] = self.conf.telescopes()[0].focalplanes()[0].name()
        for ch in self.channels:
            params[ "detector" ] = ch.tag
            if (self.calibration_file is None):
                params[ "stream" ] = "/planck/%s/raw_%s" % (self.f.inst.name, ch.tag)
            else:
                params[ "stream" ] = "/planck/%s/stack_%s" % (self.f.inst.name, ch.tag)
            params[ "noise" ] = "/planck/%s/noise_%s" % (self.f.inst.name, ch.tag)
            telescope.channel_add ( ch.tag, "native", params )

class ToastNoiseMC(ToastConfig):

    def run(self, write=True):
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

        # current RNG stream offset for the entire realization.  This is set to zero initially, and then can be incremented by Madam
        # NOT IMPLEMENTED IN TOAST
        self.conf.variable_add ( "rngbase", "native", Params({"default":"0"}) )

        if self.f.inst.name == 'LFI':
            wobble_offset = 0;
        else:
            wobble_offset = self.wobble["psi2_offset"]

        tele = self.conf.telescope_add ( "planck", "planck", 
            Params({  
                "wobblepsi2dir":self.wobble["psi2_dir"],
                "wobblepsi2_ref":self.wobble["psi2_ref"],
                "wobblepsi1_ref":self.wobble["psi1_ref"],
                "wobblepsi2_offset":wobble_offset
            }))

        fp = tele.focalplane_add ( "FP_%s" % self.f.inst.name, "planck_rimo", Params({"path":self.fpdb}) )

        self.add_pointing(tele)
        self.add_observations(tele)
        self.add_tods()
        self.add_sim_noise()
        if not self.bad_rings is None:
            self.add_bad_rings()

        params = pytoast.ParMap()
        params[ "focalplane" ] = self.conf.telescopes()[0].focalplanes()[0].name()

        for ch in self.channels:
# create stack
            stack_elements = ["nse_" + ch.tag, "raw_" + ch.tag]
            expr = ','.join(['PUSH:' + el for el in stack_elements])
            expr += ',FLG'

            if not self.bad_rings is None:
                expr += ',PUSH:bad_%s,ADD' % ch.tag

            self.strset.stream_add ( "stack_" + ch.tag, "stack", Params( {"expr":expr} ) )

            params[ "detector" ] = ch.tag
            params[ "stream" ] = "/planck/%s/stack_%s" % (self.f.inst.name, ch.tag)
            params[ "noise" ] = "/planck/%s/noise_%s" % (self.f.inst.name, ch.tag)
            tele.channel_add ( ch.tag, "native", params )

        if write:
            self.write()

    def add_sim_noise(self):
        """Add noise streams"""
# This code simply assigns the RNG stream to each TOD incrementing by one each time.
# If we require reproducibility regardless of which channel and ring combinations
# are selected, then we need to agree upon the per-channel maximum number of rings
# and on the absolute order of each channel within the list of streams used for a
# single realization.

        self.pp_boundaries = PPBoundaries(self.f.freq)

        rngstream = 0

# this is the TOAST variable holding the per-realization offset
        basename = "@rngbase@"
        noisestrm = {}

        for ch in self.channels:
               noisestrm[ch.tag] = self.strset.stream_add ( "nse_" + ch.tag, "native", Params() )
               noise = self.strset.noise_add ( "noise_" + ch.tag, "native", Params() )
               noisename = "/planck/" + self.f.inst.name + "/noise_" + ch.tag

               # add PSD
               if self.f.inst.name == "LFI":
                   noise.psd_add ( "psd", "planck_rimo", Params({
                               "start" : self.strset.observations()[0].start(),
                               "stop" : self.strset.observations()[-1].stop(),
                               "path": self.fpdb,
                               "detector": ch.tag
                   }))
               elif self.f.inst.name == "HFI":
                   noise.psd_add ( "psd", "ascii", Params({
                   "start" : self.strset.observations()[0].start(),
                   "stop" : self.strset.observations()[-1].stop(),
                   "path": private.hfi_psd
               }))

               #for pp in observation.hfiPP:
               for observation in self.observations:
                   for pp in observation.PP:
                       pp_boundaries = self.pp_boundaries.get(pp.number)
                       noisestrm[ch.tag].tod_add ( "%05d" % pp.number, "sim_noise", Params({
                           "noise" : noisename,
                           "base" : basename,
                           "start" : pp_boundaries[0],
                           "stop" : pp_boundaries[1],
                           "offset" : rngstream
                   }))
                       rngstream += 1

          
if __name__ == '__main__':
    self = ToastNoiseMC([96, 103], 30, nside=1024, ordering='RING', coord='E', efftype='C', output_xml='30_noise.xml')
    self.run()

    #toast_config = ToastConfig([95, 102], 30, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder='/project/projectdirs/planck/data/mission/lfi_ops_dx7', output_xml='30_break.xml', remote_exchange_folder='/scratch/scratchdirs/planck/data/mission/lfi_dx7s_conv/', remote_ahf_folder='/scratch/scratchdirs/planck/data/mission/AHF_v2/', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/lfi/369S/C030-0000-369S-20110713.fits')
    #toast_config = ToastConfig([91, 200], 30, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder='/project/projectdirs/planck/data/mission/lfi_ops_dx7', output_xml='30_break.xml', calibration_file='/project/projectdirs/planck/data/mission/calibration/dx7/lfi/369S/C030-0000-369S-20110713.fits', flag_HFI_bad_rings=False)
    #toast_config.run()
    #toast_config = ToastConfigCal([91, 95], 30, nside=1024, ordering='RING', coord='E', outmap='outmap.fits', exchange_folder='/project/projectdirs/planck/data/mission/lfi_ops_dx7', output_xml='30_cal_91-95.xml', efftype='C')
    #toast_config.run(False)
    #toast_config.add_exchange_format('solar_system_dipole', '/global/scratch/sd/planck/user/zonca/data/LFI_dipole_solar_system/')
    #toast_config.add_exchange_format('orbital_dipole', '/global/scratch/sd/planck/user/zonca/data/LFI_dipole_orbital/')
    #toast_config.write()
