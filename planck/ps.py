import subprocess
import pkgutil
import numpy as np
import os
import pyfits
from configobj import ConfigObj
import math

import healpy

def apply_mask(mask, m):
    m[mask == 0] = healpy.UNSEEN
    return m

def sum_diff_maps(a, b):
    '''Returns half-sum and half-difference of input maps on common pixels
    '''

    valid_pixels = np.logical_and(a != healpy.UNSEEN, b != healpy.UNSEEN)
    jackmaps = np.zeros_like([a,b])
    jackmaps[:] = healpy.UNSEEN
    jackmaps[0][valid_pixels] = ( a[valid_pixels] + b[valid_pixels]) / 2
    jackmaps[1][valid_pixels] = ( a[valid_pixels] - b[valid_pixels]) / 2

    return jackmaps

def smooth(m, arcmin, lmax = None):
    '''Utility to smooth a map with smooting by Healpix
    '''
    unseen = m == healpy.UNSEEN
    healpy.write_map('tempmap.fits', m, nest = False)
    config_filename = 'config_smooth.txt'
    config = ConfigObj()
    config.filename = config_filename 
    config['simul_type'] = 1
    if lmax:
        config['nlmax'] =  lmax
    config['infile'] = 'tempmap.fits' 
    config['outfile'] = 'tempmap_smoothed.fits'
    config['fwhm_arcmin'] = arcmin
    config.write()
    if os.path.exists('tempmap_smoothed.fits'):
        os.remove('tempmap_smoothed.fits')
    callstring = 'smoothing --double %s' % config_filename
    subprocess.call(callstring, shell=True)
    smoothed_m = healpy.read_map('tempmap_smoothed.fits')
    smoothed_m[unseen] = healpy.UNSEEN
    return smoothed_m

def remove_dipole(m, gal_cut = 30):
    #module abs path
    abspath = os.path.dirname(__file__)
    if os.path.exists('tempmap.fits'):
        os.remove('tempmap.fits')
    healpy.write_map('tempmap.fits',m,nest = False)
    callstring = 'idl %s/fixmap.pro  -IDL_QUIET 1 -quiet  -args tempmap.fits %d' % (abspath,gal_cut)
    subprocess.call(callstring, shell=True)
    out = healpy.read_map('no_dipole_tempmap.fits', nest=True)
    os.remove('no_dipole_tempmap.fits')
    return out

def anafast(m, m2=None, gal_cut = 30, lmax = None):
    '''Utility to run anafast by Healpix'''
    healpy.write_map('tempmap.fits', m, nest = False)
    config_filename = 'anafastconfig.txt'
    config = ConfigObj()
    config.filename = config_filename 
    config['simul_type'] = 1
    if gal_cut:
        config['theta_cut_deg'] = gal_cut
    if lmax:
        config['nlmax'] =  lmax
    config['infile'] = 'tempmap.fits' 
    if not m2 is None:
        healpy.write_map('tempmap2.fits', m2, nest = False)
        config['infile2'] = 'tempmap2.fits' 
    config['outfile'] = 'tempcl.fits'
    config['won'] = 0
    config.write()
    if os.path.exists('tempcl.fits'):
        os.remove('tempcl.fits')
    callstring = 'anafast --double %s' % config_filename
    subprocess.call(callstring, shell=True)
    cl = pyfits.open('tempcl.fits')[1].data.field('TEMPERATURE')
    os.remove('tempcl.fits')
    return cl

def cut_planet(lat, lon, radius = 1.5, nside = 512):
    '''Cut planet

    Example:
    Jupiter [lat -40.33,lon 33.75] radius = 1.5 deg'''

    callstring = 'idl cut_planet.pro -IDL_QUIET 1 -quiet -args %f %f %f %f' % (lat, lon, radius, nside)
    popen = subprocess.Popen(callstring, shell=True, stdout=subprocess.PIPE, cwd = os.path.dirname(__file__))
    output = popen.communicate()[0]
    return map(float, output.strip().split())
