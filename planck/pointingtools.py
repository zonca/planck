from __future__ import division

#TODO remove useless imports
import math
import pyfits
import logging as l 
import csv
import glob
from itertools import *
from Quaternion import Quat as quat
import numpy as np
from IPython.Debugger import Tracer; debug_here = Tracer()
import re
import quaternionarray as qarray
from utils import grouper
import private
from cgkit.cgtypes import *
from pointingtools import *
import exceptions

try:
    import pysqlite2.dbapi2 as sqlite3
except:
    import sqlite3

QECL2GAL_HYDRA = np.array((-0.37382079227204573, 0.33419217216073838, 0.64478939348298625, 0.57690575088960561))
QECL2GAL_HEALPIX = np.array([-0.37381693504678937, 0.33419069514234978, 0.64479285220138716, 0.57690524015582401])
QECL2GAL = QECL2GAL_HEALPIX

class Siam(object):

    def __init__(self, horn_pointing=False):
        self.horn_pointing = horn_pointing
        siamfile = private.siam
        l.debug('using SIAM %s' % siamfile)
        f = open(siamfile)
        lines = f.readlines()
        self.siam = {}
        for line in grouper(4,lines[1:]):
            chtag = line[0].split()[0]
            m = np.array(np.matrix(';'.join(line[1:])))
            self.siam[chtag] = m
    def get(self, ch):
        if ch.inst.name == 'HFI':
            return self.siam[ch.tag].T
        else:
            l.warning('For LFI using instrument DB angles')
            return SiamAngles(self.horn_pointing).get(ch)

class SiamAngles(object):

    SPIN2BORESIGHT = 85.0
    
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
        mat_spin2boresight=mat3.rotation(np.pi/2-self.SPIN2BORESIGHT/180.*np.pi,vec3(0,1,0))
        theta, phi, psi = self.get_angles(ch)
        #debug_here()

        mat_theta_phi = mat3.rotation(theta,vec3(-math.sin(phi),math.cos(phi),0))
        mat_psi = mat3.rotation(psi,vec3(0,0,1))
        # detector points to X axis
        total = mat_spin2boresight * (mat_theta_phi * mat_psi)
        total_mat = np.matrix(np.split(np.array(total.toList(rowmajor=True)),3))
        # siam is defined as pointing to Z axis
        return np.array(total_mat * np.matrix([[0,0,1],[0,1,0],[1,0,0]]))

class SiamForcedAngles(SiamAngles):

    def __init__(self, theta, phi, psi):
        self.theta = theta
        self.phi = phi
        self.psi = psi

    def get_angles(self, ch):
        return self.theta, self.phi, self.psi

def AHF_btw_OBT(obt):

    conn = sqlite3.connect('/project/projectdirs/planck/user/zonca/remix/database.db')
    c = conn.cursor()
    values = (obt[0]*2**16, obt[-1]*2**16)
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
    qsatgal = qarray.norm(qsatgal)
    return qsatgal

def vector_ecl2gal(vecl):
    '''Convert arrays from Ecliptic to Galactic'''
    l.info('Rotating to Galactic frame')
    return qarray.rotate(QECL2GAL ,vecl)

def vector_gal2ecl(vecl):
    '''Convert arrays from Ecliptic to Galactic'''
    l.info('Rotating to Galactic frame')
    return qarray.rotate(qarray.inv(QECL2GAL) ,vecl)
