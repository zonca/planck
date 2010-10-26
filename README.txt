Python package for dealing with Planck satellite data
part of the US Planck Test Environment

license: GPL v3
author: Andrea Zonca 
website: http://andreazonca.com

Access to Planck Data is needed to make use of this package, this
software does not include any performance number or any other 
information covered by the Planck Data Agreement.

If you are in the Planck collaboration and interested in using 
and contributing to the software please contact me.

Includes:

* planck, LFI, HFI: metadata classes for LFI and HFI channels
created dynamically from the Reduced Instrument Model (RIMO),
not publicly available, not even channel names are available in
this package.

* pointing: pointing library which builds detector pointing from
satellite quaternions, it is based on quaternionarray
[http://github.com/zonca/quaternionarray]

* utils, ps, hitmap: utilities for date conversion, angular power spectra
and hitmaps
