
from __future__ import absolute_import, division, print_function

import numpy
import cosmolopy.perturbation as cp
import cosmolopy.parameters as cparam
import pylab

# Get WMAP5 cosmological parameters.
cosmo = cparam.WMAP5_BAO_SN_mean() 

# Set up an array of z values.
z = numpy.arange(6.0, 20.0, 0.1) 

# Calculate collapse fraction.
f_col = cp.collapse_fraction(*cp.sig_del(1e4, z, **cosmo))

pylab.figure(figsize=(8,6))
pylab.plot(z, f_col)
pylab.xlabel('z')
pylab.ylabel(r'$f_\mathrm{col}$')
pylab.show()
