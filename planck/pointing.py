from __future__ import division

import logging as l 
import numpy as np
from LFI import LFI
from IPython.Debugger import Tracer; debug_here = Tracer()
import re
import quaternionarray as qarray
from utils import grouper
import Planck
import private
from pointingtools import *
import pycfitsio
import correction
import glob

class Pointing(object):
    '''Pointing interpolation and rotation class
    
    usage:
    >>> ch= Planck()['100-1a']
    >>> pnt = Pointing(obt, coord='G') #interpolates AHF to obt
    >>> vec = pnt.get(ch) #rotates to detector frame and gives x,y,z vector
    >>> pix = pnt.get_pix(ch, 2048, nest=True) #healpix pixel number nside 2048
    '''
    comp =  ['X','Y','Z','S']

    def __init__(self,obt,coord='G', horn_pointing=False, deaberration=True, wobble=True, interp='slerp', siamfile=None, wobble_offset=0):
        '''
        nointerp to use the AHF OBT stamps'''
        l.warning('Pointing setup, coord:%s, deab:%s, wobble:%s' % (coord, deaberration, wobble))
        #get ahf limits
        self.deaberration = deaberration
        self.wobble = wobble

        filenames = AHF_btw_OBT(obt)
        files = [pycfitsio.open(f, False) for f in filenames]
        l.debug('reading files %s' % str(files))
        AHF_data_iter = [f[0] for f in files]

        l.debug('reading files')

        ahf_obt = np.concatenate([h.read_column('OBT_SPL') for h in AHF_data_iter])
        ahf_obt /= 2.**16
        i_start = max(ahf_obt.searchsorted(obt[0])-1,0)
        i_end = min(ahf_obt.searchsorted(obt[-1])+1,len(ahf_obt)-1)
        ahf_obt = ahf_obt[i_start:i_end]

        ahf_quat = np.empty((len(ahf_obt),4))
        for i,c in enumerate(self.comp):
            ahf_quat[:,i] = np.concatenate([h.read_column('QUATERNION_'+c) for h in AHF_data_iter])[i_start:i_end]

        if coord == 'G':
            ahf_quat = quaternion_ecl2gal(ahf_quat)

        #debug_here()
        if self.wobble:
           ahf_quat = qarray.mult(ahf_quat, correction.wobble(ahf_obt,offset=wobble_offset))
           qarray.norm_inplace(ahf_quat)

        if interp is None:
            self.qsatgal_interp = ahf_quat 
            # save AHF obt for later interpolation
            self.ahf_obt = ahf_obt
        else:
            l.info('Interpolating quaternions with %s' % interp)
            interpfunc = getattr(qarray, interp)
            self.qsatgal_interp = interpfunc(obt, ahf_obt, ahf_quat)

        #if self.wobble:
        #    self.qsatgal_interp = qarray.mult(self.qsatgal_interp, correction.wobble(obt))
        #    qarray.norm_inplace(self.qsatgal_interp)

        l.info('Quaternions interpolated')
        self.siam = Siam(horn_pointing, siamfile=siamfile)

        self.obt = obt
        self.coord = coord

        l.debug('Closing AHF files')
        for f in files:
            f.close()

    def interp_get(self, rad):
        '''Interpolation after rotation to gal frame'''
        from Quaternion import Quat
        l.info('Rotating to detector %s' % rad)
        siam_quat = Quat(self.siam.get(rad)).q
        totquat = qarray.mult(self.qsatgal_interp, siam_quat)
        totquat_interp = qarray.nlerp(self.obt, self.ahfobt, totquat)
        x = np.array([1, 0, 0])
        vec = qarray.rotate(totquat_interp, x)
        l.info('Rotated to detector %s' % rad)
        return vec

    def get(self, rad):
        rad = Planck.parse_channel(rad)
        l.info('Rotating to detector %s' % rad)
        x = np.dot(self.siam.get(rad),[1, 0, 0])
        vec = qarray.rotate(self.qsatgal_interp, x)
        qarray.norm_inplace(vec)
        if self.deaberration:
            l.warning('Applying deaberration correction')
            vec += correction.simple_deaberration(vec, self.obt, self.coord)
            qarray.norm_inplace(vec)
        l.info('Rotated to detector %s' % rad)
        return vec

    def inv(self, rad, vec):
        rad = Planck.parse_channel(rad)
        l.info('Rotating to detector %s' % rad)
        if self.deaberration:
            l.warning('Applying deaberration correction')
            vec -= correction.simple_deaberration(vec, self.obt, self.coord)
            qarray.norm_inplace(vec)
        vec_rad = qarray.rotate(qarray.inv(self.qsatgal_interp), vec)
        invsiam = np.linalg.inv(self.siam.get(rad))
        #invsiamquat = qarray.inv(qarray.norm(qarray.from_rotmat(self.siam.get(rad))))
        #qarray.rotate(invsiamquat, vec_rad)
        return np.array([np.dot(invsiam , row) for row in vec_rad])

    def compute_psi(self, theta, phi, rad):
        z = np.dot(self.siam.get(rad),[0, 0, 1])
        vecz = qarray.norm(qarray.rotate(self.qsatgal_interp, z))
        e_phi = np.hstack([-np.sin(phi)[:,np.newaxis], np.cos(phi)[:,np.newaxis], np.zeros([len(phi),1])])
        e_theta = np.hstack([(np.cos(theta)*np.cos(phi))[:,np.newaxis], (np.cos(theta)*np.sin(phi))[:,np.newaxis], -np.sin(theta)[:,np.newaxis]])
        psi = np.arctan2(-qarray.arraylist_dot(vecz, e_phi), -qarray.arraylist_dot(vecz, e_theta))
        return psi.flatten()

    def get_vecpsi(self, rad):
        from healpy import vec2ang
        vec = self.get(rad)
        theta, phi = vec2ang(vec)
        psi = self.compute_psi(theta, phi, rad)
        return vec, psi

    def get_3ang(self, rad):
        l.info('Rotating to detector %s' % rad)
        theta, phi = self.get_ang(rad)
        psi = self.compute_psi(theta, phi, rad)
        #psi = np.arcsin(-qarray.arraylist_dot(vecz, e_phi))
        l.info('Rotated to detector %s' % rad)
        return theta, phi, psi

    def get_pix_iqu(self, rad, nside=1024, nest=True):
        from healpy import vec2pix, vec2ang
        vec = self.get(rad)
        theta, phi = vec2ang(vec)
        psi = self.compute_psi(theta, phi, rad)
        #return vec2pix(nside, vec[:,0], vec[:,1], vec[:,2], nest), np.cos(2*psi), np.sin(2*psi)
        spsi = np.sin(psi)
        cpsi = np.cos(psi)
        cf = 1./(cpsi**2 + spsi**2)
        return vec2pix(nside, vec[:,0], vec[:,1], vec[:,2], nest), (cpsi**2 - spsi**2)*cf, 2*cpsi*spsi*cf

    def get_pix(self, rad, nside=1024, nest=True):
        from healpy import vec2pix
        vec = self.get(rad)
        return vec2pix(nside, vec[:,0], vec[:,1], vec[:,2], nest)

    def get_ang(self, rad, degrees=False):
        from healpy import vec2ang
        vec = self.get(rad)
        ang = vec2ang(vec)
        if degrees:
            return map(np.rad2deg, ang)
        else:
            return ang

class DiskPointing(Pointing):
    '''Read pointing from disk'''
    def __init__(self, od, freq, folder=None):
        self.folder = folder or private.pointingfolder
        self.filename = glob.glob(self.folder + '/%04d/?%03d-*.fits' % (od,freq))[0]

    def get_3ang(self, ch):
        l.debug('Reading %s' % self.filename)
        with pycfitsio.open(self.filename) as f:
            h = f[ch.tag]
            return h.read_column('THETA'), h.read_column('PHI'), h.read_column('PSI')

    def get(self, ch):
        import healpy
        theta, phi, psi = self.get_3ang(ch)
        return healpy.ang2vec(theta, phi)
