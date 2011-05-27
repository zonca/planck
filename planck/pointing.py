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
from cgkit.cgtypes import *
import physcon
import pfits
import correction

class Pointing(object):
    '''Pointing interpolation and rotation class
    
    usage:
    >>> ch= Planck()['100-1a']
    >>> pnt = Pointing(obt, coord='G') #interpolates AHF to obt
    >>> vec = pnt.get(ch) #rotates to detector frame and gives x,y,z vector
    >>> pix = pnt.get_pix(ch, 2048, nest=True) #healpix pixel number nside 2048
    '''

    fields= ['OBT_SPL','QUATERNION_X','QUATERNION_Y','QUATERNION_Z','QUATERNION_S']

    def __init__(self,obt,coord='G', AHF_d=None, nointerp=False, horn_pointing=False, deaberration=True, wobble=True):
        '''AHF_d is the pfits AHF data if already loaded in the main file
        nointerp to use the AHF OBT stamps'''
        l.warning('Pointing setup, coord:%s, deab:%s, wobble:%s' % (coord, deaberration, wobble))
        #get ahf limits
        self.deaberration = deaberration
        self.wobble = wobble

        if  AHF_d is None:
            files = AHF_btw_OBT(obt)
            l.debug('reading files %s' % str(files))
            AHF_data_iter = [pfits.FITS(f+'[ATT-HIST-HGH][col %s]' % (';'.join(self.fields))).get_hdus()[1].get_data() for f in files]
        else:
            AHF_data_iter = AHF_d

        ahfobt = np.array([])
        qsat = None
        l.debug('concatenating quaternions')
        for AHF_data in AHF_data_iter:
            AHF_data['OBT_SPL'] /= 2.**16

        i_start = max(AHF_data_iter[0]['OBT_SPL'].searchsorted(obt[0])-1,0)
        i_end = min(AHF_data_iter[-1]['OBT_SPL'].searchsorted(obt[-1])+1,len(AHF_data_iter[-1]['OBT_SPL'])-1)
        for field in self.fields:
            AHF_data_iter[0][field]=AHF_data_iter[0][field][i_start:]
            AHF_data_iter[-1][field]=AHF_data_iter[-1][field][:i_end]
        
        qsat = np.hstack([np.concatenate([AHF_data['QUATERNION_%s' % comp] for AHF_data in AHF_data_iter]).reshape((-1,1)) for comp in ['X','Y','Z','S']])
        self.ahfobt = np.concatenate([AHF_data['OBT_SPL'] for AHF_data in AHF_data_iter])

        if coord == 'E':
            qsatgal = qsat
        elif coord == 'G':
            qsatgal = quaternion_ecl2gal(qsat)

        #debug_here()
        if self.wobble and len(self.ahfobt)<len(obt):
            qsatgal = qarray.mult(qsatgal, correction.wobble(self.ahfobt))

        if nointerp:
            self.qsatgal_interp = qsatgal 
        else:
            l.info('Interpolating quaternions')
            #nlerp
            self.qsatgal_interp = qarray.nlerp(obt, self.ahfobt, qsatgal)

        if self.wobble and len(self.ahfobt)>=len(obt):
            self.qsatgal_interp = qarray.mult(self.qsatgal_interp, correction.wobble(obt))

        l.info('Quaternions interpolated')
        self.siam = Siam(horn_pointing)

        self.obt = obt
        self.coord = coord

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
        #if self.wobble:
        #    x = qarray.rotate(correction.wobble(self.obt), x)
        #    qarray.norm_inplace(x)
        vec = qarray.rotate(self.qsatgal_interp, x)
        qarray.norm_inplace(vec)
        if self.deaberration:
            l.warning('Applying deaberration correction')
            vec += correction.simple_deaberration(vec, self.obt, self.coord)
            qarray.norm_inplace(vec)
        l.info('Rotated to detector %s' % rad)
        return vec

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
