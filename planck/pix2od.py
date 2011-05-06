#!/usr/bin/env python
import datetime
import sys
import numpy as np
from optparse import OptionParser
import sys
import os.path
import exceptions

try:
    from bitstring import ConstBitArray
except exceptions.ImportError:
    print('pix2od requires bitstring: easy_install bitstring')
    sys.exit(1)

PIX2ODPATH = 'data'
OBTSTARTDATE = datetime.datetime(1958,1,1,0,0,0)
LAUNCH = datetime.datetime(2009, 5, 13, 13, 11, 57, 565826)
SECONDSPERDAY = 3600 * 24

def pix2od(ch, pixels):
    """ nside 512 NEST pixel number to OD [91-563 excluding 454-455] for Planck channels
    it returns a list of sets.
    Each set contains the ODs hit by each pixel
    """
    filename = os.path.join(PIX2ODPATH, 'od_by_pixel_%s.bin' % ch.replace('M','S'))
    pixels_by_od = ConstBitArray(filename = filename)
    odrange = [91, 563+1]
    num_ods = odrange[1] - odrange[0]
    NSIDE = 512
    NPIX = 12 * NSIDE**2
    tot_ods = list()
    if np.any(np.array(pixels) >= NPIX):
        raise exceptions.ValueError('ERROR: input pixels must be NSIDE 512, RTFM!')
    for pix in pixels:
        ods = set(np.array(list(pixels_by_od[num_ods*pix:num_ods*pix+num_ods].findall([True]))) + odrange[0])
        tot_ods.append(ods)
    return tot_ods


if __name__ == '__main__':
    usage = '''nside 512 NEST pixel number to OD [91-563 excluding 454-455] for Planck channels
               pix2od -c LFI27M 1000 [1001 1002...]
               returns list of operational days
               pointing generated with TestEnv using HORN pointing (M has same pointing as S)
               Using Galactic coordinates:
               pix2od -c LFI28S -a -- Lat[+-90 deg] Long[+-180 deg]
                 '''
    parser = OptionParser (usage = usage)
    parser.add_option('-c','--channel', dest = 'ch', 
                    help = 'channel tag (def %default)', action = 'store', type = 'string',
                    default = 'LFI28M')
    parser.add_option('-a','--angles', dest = 'angles', 
                    help = 'Use galactic latitude [+-90] and longitude [+-180] in degrees instead of pixel numbers (def %default)', action = 'store_true', 
                    default = False)

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)
    NSIDE=512
    if options.angles:
        tot_ods = []
        from healpy import ang2pix
        args = map(float, args)
        pixels = []
        for latitude, longitude in zip(args[::2], args[1::2]):
            pixels.append(ang2pix(NSIDE, np.radians(90-latitude), np.radians(longitude),nest=True))
        tot_ods = set.intersection(*utils.pix2od(options.ch, pixels))
    else:
        pixels = map(int, args)
        tot_ods = set.union(*utils.pix2od(options.ch, pixels))
        
    all_ods = np.unique(tot_ods)
    print(' '.join(map(str,list(all_ods))))
    sys.exit(0)
