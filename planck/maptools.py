import numpy as np
import healpy as hp

def degrade_mask(m, nside):
    return np.ceil(hp.ud_grade(m.astype(np.float), nside)).astype(np.bool)
