from __future__ import division

#TODO remove useless imports
import math
import pyfits
import logging as l 
import csv
import glob
from itertools import *
import numpy as np
#from IPython.Debugger import Tracer; debug_here = Tracer()
import re
import quaternionarray as qarray
from utils import grouper
import private
from pointingtools import *
import exceptions

SPIN2BORESIGHT = np.radians(85.0)

def angles2siam(theta, phi, psi):
    mat_spin2boresight=qarray.rotation([0,1,0], np.pi/2-SPIN2BORESIGHT)

    mat_theta_phi = qarray.rotation([-math.sin(phi),math.cos(phi),0], theta)
    mat_psi = qarray.rotation([0,0,1], psi)
    # detector points to X axis
    total = qarray.mult(mat_spin2boresight, qarray.mult(mat_theta_phi, mat_psi))
    # siam is defined as pointing to Z axis
    return np.dot(qarray.to_rotmat(total[0]), np.array([[0,0,1],[0,1,0],[1,0,0]]))

try:
    import pysqlite2.dbapi2 as sqlite3
except:
    import sqlite3

QECL2GAL_HYDRA = np.array((-0.37382079227204573, 0.33419217216073838, 0.64478939348298625, 0.57690575088960561))
QECL2GAL_HEALPIX = np.array([-0.37381693504678937, 0.33419069514234978, 0.64479285220138716, 0.57690524015582401])
QECL2GAL = QECL2GAL_HEALPIX

class Siam(object):

    def __init__(self, horn_pointing=False, siamfile=None):
        self.horn_pointing = horn_pointing
        if siamfile is None:
            siamfile = private.siam
        f = open(siamfile)
        lines = f.readlines()
        self.siam = {}
        for line in grouper(4,lines[1:]):
            chtag = line[0].split()[0]
            m = np.array(np.matrix(';'.join(line[1:])))
            self.siam[chtag] = m
        self.siamfile = siamfile
    def get(self, ch):
        if ch.inst.name == 'HFI':
            l.debug('using SIAM %s' % self.siamfile)
            return self.siam[ch.tag].T
        else:
            l.warning('For LFI using instrument DB angles')
            return SiamAngles(self.horn_pointing).get(ch)

class SiamAngles(object):

    
    def __init__(self, horn_pointing):
        self.horn_pointing = horn_pointing

    def get_angles(self, ch):
        if self.horn_pointing and ch.arm == 'M':
            l.warning('USING HORN POINTING')
            S_ch = ch.inst[ch.tag.replace('M','S')]
            theta = np.radians(S_ch.get_instrument_db_field('theta_uv'))
            phi = np.radians(S_ch.get_instrument_db_field('phi_uv'))
        else:
            theta = np.radians(ch.get_instrument_db_field('theta_uv'))
            phi = np.radians(ch.get_instrument_db_field('phi_uv'))
        psi = np.radians(ch.get_instrument_db_field('psi_uv')+ch.get_instrument_db_field('psi_pol'))
        return theta, phi, psi

    def get(self, ch):
        return angles2siam(*self.get_angles(ch))

class SiamForcedAngles(SiamAngles):

    def __init__(self, theta, phi, psi):
        self.theta = theta
        self.phi = phi
        self.psi = psi

    def get_angles(self, ch):
        return self.theta, self.phi, self.psi

def AHF_btw_OBT(obt):

    conn = sqlite3.connect(private.database)
    c = conn.cursor()
    values = ((obt[0]-60*5)*2**16, (obt[-1]+60*5)*2**16)
    query = c.execute('select file_path from ahf_files where endOBT>=? and startOBT<=?', values)
    files = [q[0] for q in query]
    c.close()
    return files

def generate_repointing_flag(obt):
    flag = np.zeros_like(obt)
    files = [pyfits.open(file)[1].data for file in AHF_btw_OBT(obt)]

    files[-1] = files[-1][:(files[-1].field('OBT_SPL')/2.**16).searchsorted(obt[-1])+1]
    files[0] = files[0][(files[0].field('OBT_SPL')/2.**16).searchsorted(obt[0])-1:]
    AHF = np.concatenate(files)

    i_start_repointing, = np.nonzero(np.diff(AHF['OBT_BEG']))
    start_repointing = AHF['OBT_SPL'][i_start_repointing+1]/2.**16
    end_repointing = AHF['OBT_BEG'][i_start_repointing+1]/2.**16
    for start, end in zip(start_repointing,end_repointing):
        flag[obt.searchsorted(start):obt.searchsorted(end)] = 1
    return flag

def quaternion_ecl2gal(qsat):
    '''Convert array of quaternions from Ecliptic to Galactic'''
    l.info('Rotating to Galactic frame')
    qsatgal = qarray.mult(QECL2GAL ,qsat)
    # renormalizing to unity
    qarray.norm_inplace(qsatgal)
    return qsatgal

def vector_ecl2gal(vecl):
    '''Convert arrays from Ecliptic to Galactic'''
    l.info('Rotating to Galactic frame')
    return qarray.rotate(QECL2GAL ,vecl)

def vector_gal2ecl(vecl):
    '''Convert arrays from Ecliptic to Galactic'''
    l.info('Rotating to Galactic frame')
    return qarray.rotate(qarray.inv(QECL2GAL) ,vecl)
