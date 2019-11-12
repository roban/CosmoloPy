"""CosmoloPy is a package of cosmology routines built on NumPy/SciPy.

Capabilities include
--------------------

`cosmolopy.density`
  Various cosmological densities.

`cosmolopy.distance`
  Various cosmological distance measures.

`cosmolopy.luminosityfunction`
  Galaxy luminosity functions (Schecter functions).

`cosmolopy.magnitudes`
  Conversion in and out of the AB magnitude system.

`cosmolopy.parameters`
  Pre-defined sets of cosmological parameters (e.g. from WMAP).

`cosmolopy.perturbation`
  Perturbation theory and the power spectrum.

`cosmolopy.reionization`
  The reionization of the IGM.

Functions take cosmological parameters (which can be numpy arrays)
as keywords, and ignore any extra keywords. This means you can make a
dictionary of all of your cosmological parameters and pass it to any
function.  

The `parameters` module supplies some convenient pre-defined parameter sets.

Usage
-----

The easiest way to use CosmoloPy is to create a dictionary of the
cosmology parameters and pass it to each function using the ** syntax.

>>> import cosmolopy.distance as cd
>>> cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
>>> d_co = cd.comoving_distance(6., **cosmo)
>>> print "Comoving distance to z=6 is %.1f Mpc" % (d_co)
Comoving distance to z=6 is 8017.8 Mpc

The cosmolopy package also defines some convenient shortcuts,
including a fiducial cosmology (currently the WMAP7+BAO+H0 mean), so
you can just do this:

>>> from cosmolopy import *
>>> d_a = cd.angular_diameter_distance(6, **fidcosmo)
>>> print "Angluar-diameter distance to z=6 is %.1f Mpc" % (d_a)
Angluar-diameter distance to z=6 is 1209.9 Mpc
>>> d_light = cd.light_travel_distance(6, **fidcosmo)
>>> print "Light-travel distance to z=6 is %.1f Mpc" % (d_light)
Light-travel distance to z=6 is 3922.9 Mpc

Calculate the mass of a halo with Virial temperature of 10^4 kelvin,
then verify the Virial temperature for a halo of that mass:

>>> import cosmolopy.perturbation as cp
>>> cosmo = {'omega_M_0' : 0.27, 
... 'omega_lambda_0' : 1-0.27, 
... 'omega_b_0' : 0.045, 
... 'omega_n_0' : 0.0,
... 'N_nu' : 0,
... 'h' : 0.72,
... 'n' : 1.0,
... 'sigma_8' : 0.9
... } 
>>> mass = cp.virial_mass(1e4, 6.0, **cosmo)
>>> temp = cp.virial_temp(mass, 6.0, **cosmo)
>>> print "Mass = %.3g M_sun" % mass
Mass = 1.68e+08 M_sun
>>> print round(temp, 4)
10000.0

Calculate the critical and matter densities:

>>> from cosmolopy import *
>>> 'rho_crit=%.3g Msun/Mpc^3, rho_0=%.3g Msun/Mpc^3' % cden.cosmo_densities(**fidcosmo)
'rho_crit=1.38e+11 Msun/Mpc^3, rho_0=3.75e+10 Msun/Mpc^3'

Look in the tests/ and examples/ directories for more examples. 

"""

from __future__ import absolute_import, division, print_function

from .__version__ import __version__
from . import constants as cc
from . import density as cden
from . import distance as cd
from . import perturbation as cp
from . import reionization as cr
from . import parameters as cparam
from . import magnitudes as cmag
from . import luminosityfunction as cl

fidcosmo = cparam.WMAP7_BAO_H0_mean(flat=True, extras=True)
