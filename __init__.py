"""CosmoloPy: A package of cosmology routines built on NumPy/SciPy.

Most functions in the modules of this package are designed to take
cosmological parameters as keywords and will ignore any extraneous
keywords. This means you can make a dictionary of all of your cosmological
parameters and pass it to any function. 

The parameters module supplies some convenient pre-defined parameter sets.

Examples:
--------

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

"""
