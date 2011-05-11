import numpy as np

from cgkit.cgtypes import *
import physcon

from dipole import SatelliteVelocity

def deaberration(vec, obt, coord):
    satvel = SatelliteVelocity(coord).orbital_v(obt)
    return np.cross(vec, np.cross(vec, satvel/physcon.c))

def wobble(od):

    R_psi1 = mat3.rotation(private.WOBBLE['psi1_ref'], vec3(0,0,1)).transpose()
    R_psi2 = mat3.rotation(private.WOBBLE['psi2_ref'], vec3(0,1,0)).transpose()

    psi2 = get_wobble_psi2(od)
    R_psi2T = mat3.rotation(psi2, vec3(0,1,0))

    wobble_rotation = R_psi1.transpose() * (R_psi2T * (R_psi2 * R_psi1))

    return wobble_rotation
