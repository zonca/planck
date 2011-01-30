from __future__ import division

import pyfits
import logging as l 
import numpy as np
from LFI import LFI
#from IPython.Debugger import Tracer; debug_here = Tracer()
import re
import quaternionarray as qarray
from utils import grouper
import Planck
import private
from pointingtools import *

class Pointing(object):
    '''Pointing interpolation and rotation class
    
    usage:
    >>> ch= Planck()['100-1a']
    >>> pnt = Pointing(obt, coord='G') #interpolates AHF to obt
    >>> vec = pnt.get(ch) #rotates to detector frame and gives x,y,z vector
    >>> pix = pnt.get_pix(ch, 2048, nest=True) #healpix pixel number nside 2048
    '''

    def __init__(self,obt,coord='G', AHF_d=None, nointerp=False):
        '''AHF_d is the pyfits AHF data if already loaded in the main file
        nointerp to use the AHF OBT stamps'''
        l.warning('Pointing setup, coord:%s' % coord)
        #get ahf limits

        if  AHF_d is None:
            AHF_data_iter = (pyfits.open(file)[1].data for file in AHF_btw_OBT(obt))
        else:
            AHF_data_iter = [AHF_d]

        ahfobt = np.array([])
        qsat = None
        for AHF_data in AHF_data_iter:

            obt_spl = AHF_data.field('OBT_SPL')/2.**16
            i_start = max(obt_spl.searchsorted(obt[0])-1,0)
            i_end = min(obt_spl.searchsorted(obt[-1])+1,len(obt_spl)-1)
            AHF = AHF_data[i_start:i_end]

            allquat = np.hstack([AHF.field('QUATERNION_X')[:,np.newaxis], AHF.field('QUATERNION_Y')[:,np.newaxis], AHF.field('QUATERNION_Z')[:,np.newaxis], AHF.field('QUATERNION_S')[:,np.newaxis]])

            if qsat is None:
                qsat = allquat
            else:
                qsat = np.vstack([qsat,allquat])
            ahfobt = np.concatenate([ahfobt, AHF.field('OBT_SPL')/2.**16])

        if coord == 'E':
            qsatgal = qsat
        elif coord == 'G':
            hfobt = np.array([])
            qsatgal = quaternion_ecl2gal(qsat)

        if nointerp:
            self.qsatgal_interp = qsatgal 
        else:
            l.info('Interpolating quaternions')
            #nlerp
            self.qsatgal_interp = qarray.nlerp(obt, ahfobt, qsatgal)

        l.info('Quaternions interpolated')
        self.siam = Siam()

        self.ahfobt = ahfobt
        self.obt = obt

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
        l.info('Rotating to detector %s' % rad)
        x = np.dot(self.siam.get(rad),[1, 0, 0])
        vec = qarray.norm(qarray.rotate(self.qsatgal_interp, x))
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
        return vec2pix(nside, vec[:,0], vec[:,1], vec[:,2], nest), np.ones(len(psi)), np.cos(2*psi), np.sin(2*psi)

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
