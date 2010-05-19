"""Calculate various cosmological distance measures. 

Mostly follows David Hogg's pedagogical paper arXiv:astro-ph/9905116v4 .

"""

import math

import numpy
import scipy
import scipy.integrate as si
import scipy.interpolate
import scipy.optimize

import constants as cc

def get_omega_k_0(**cosmo):
    """Calculate omega_k_0 for a cosmology, if needed.

    If omega_k_0 is specified, return it. Otherwise return 1.0 -
    omega_M_0 - omega_lambda_0

    """
 
    if 'omega_k_0' in cosmo:
        omega_k_0 = cosmo['omega_k_0']
    else:
        omega_k_0 = 1. - cosmo['omega_M_0'] - cosmo['omega_lambda_0']
    return omega_k_0

def set_omega_k_0(cosmo):
    """Returns the cosmo dictionary with omega_k_0 set.
    See get_omega_k_0.
    
    Note that cosmo is not passed as **cosmo for once. This function
    modifies the dictionary in place and returns the result.

    """
    if 'omega_k_0' in cosmo:
        return cosmo
    else:
        cosmo['omega_k_0'] = get_omega_k_0(**cosmo)
        return cosmo


### Distance measures ###

def e_z(z, **cosmo):
    """Calculate the unitless Hubble expansion rate at redshift z.

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to E(z), defined in his eq. 14.

    """

    return (cosmo['omega_M_0'] * (1+z)**3. + 
            cosmo['omega_k_0'] * (1+z)**2. + 
            cosmo['omega_lambda_0'])**0.5

def hubble_z(z, **cosmo):
    """Calculate the value of the Hubble constant at redshift z.

    Units are s^-1

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to H_0 * E(z) (see his eq. 14).

    """
    H_0 = cosmo['h'] * cc.H100_s

    return H_0 * e_z(z, **cosmo)

def hubble_distance_z(z, **cosmo):
    """Calculate the value of the Hubble distance at redshift z.

    Units are Mpc.

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to D_H / E(z) = c / (H_0 E(z)) [see his eq. 14], which
    appears in the definitions of many other distance measures.

    """
    H_0 = cosmo['h'] * cc.H100_s
    
    return cc.c_light_Mpc_s / (H_0 * e_z(z, **cosmo))

def _comoving_integral(z, omega_M_0, omega_lambda_0, omega_k_0, h):

    e_z = (omega_M_0 * (1+z)**3. + 
           omega_k_0 * (1+z)**2. + 
           omega_lambda_0)**0.5
    
    H_0 = h * cc.H100_s
    
    H_z =  H_0 * e_z

    return cc.c_light_Mpc_s / (H_z)

def comoving_distance(z, z0 = 0, **cosmo):
    """Calculate the line-of-sight comoving distance (in Mpc) to redshift z.

    See equation 15 of David Hogg's arXiv:astro-ph/9905116v4

    Optionally calculate the integral from z0 to z.

    Returns:
    -------
    
    d_co: ndarray
       Comoving distance in Mpc.

    err: ndarray
       Extimate of numerical integration error from scipy.integrate.quad.

    Examples:
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> d_co, err = cd.comoving_distance(6., **cosmo)
    >>> print "Comoving distance = %.1f Mpc" % (d_co)
    Comoving distance = 8017.8 Mpc

    """

    #cosmo = set_omega_k_0(cosmo)

    dc_func = \
        numpy.vectorize(lambda z, z0, omega_M_0, omega_lambda_0, omega_k_0, h: 
                        si.quad(_comoving_integral, z0, z, limit=1000,
                                args=(omega_M_0, omega_lambda_0, omega_k_0, h)))
    d_co, err = dc_func(z, z0, 
                        cosmo['omega_M_0'],
                        cosmo['omega_lambda_0'],
                        cosmo['omega_k_0'],
                        cosmo['h']
                        )
    return d_co, err

def propper_motion_distance(**args):
    """Returns comoving_distance_transverse."""
    return comoving_distance_transverse(z, **cosmo)

def comoving_distance_transverse(z, **cosmo):
    """Calculate the transverse comoving distance (in Mpc) to redshift z.

    This is also called the proper motion distance, D_M.

    See equation 16 of David Hogg's arXiv:astro-ph/9905116v4

    This is the distance d_m, such that the comoving distance between
    two events at the same redshift, but separated on the sky by some
    angle delta_theta is d_m * delta_theta.

    Warning: currently returns the error on the line-of-sight comoving
    distance D_C, not the error on the transverse comoving distance
    D_M.

    """

    d_c, d_c_err = comoving_distance(z, 0.0, **cosmo)

    # We use atleast_1d to allow us to write code for arrays of omega
    # values.
    #omega_k_0 = numpy.atleast_1d(get_omega_k_0(**cosmo))
    omega_k_0 = get_omega_k_0(**cosmo)
    

    if numpy.all(omega_k_0 == 0.0):
        return d_c, d_c_err
    
    d_h_0 = hubble_distance_z(0.0, **cosmo)
    sqrt_ok0 = numpy.sqrt(numpy.abs(omega_k_0))
    sqrt_ok0[omega_k_0 == 0.0] = 1.0
    argument = sqrt_ok0 * d_c / d_h_0
    factor = d_h_0 * (1./sqrt_ok0)
    
    d_m = ((omega_k_0 > 0.0) * (factor * numpy.sinh(argument)) +
           (omega_k_0 == 0.0) * d_c +
           (omega_k_0 < 0.0) * (factor * numpy.sin(argument)))

    return d_m, d_c_err

def angular_diameter_distance(z, z0 = 0, **cosmo):
    """Calculate the angular-diameter distance (Mpc) to redshift z.
    
    Optionally find the angular diameter distante between objects at
    z0 and z.

    See equations 18-19 of David Hogg's arXiv:astro-ph/9905116v4

    Warning: returns two error estimates, one from each invocation of
    comoving_distance_transverse (see the docstring for that function
    for explanation of the error estimate).

    Examples:
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> d_a, err1, err2 = cd.angular_diameter_distance(6., **cosmo)
    >>> print "Angular diameter distance = %.1f Mpc" % (d_a)
    Angular diameter distance = 1145.4 Mpc

    """

    omega_k = numpy.atleast_1d(get_omega_k_0(**cosmo))
    if (numpy.any(omega_k < 0)):
        raise ValueError("Not implemented for Omega_k < 0")

    dm2, err2 = comoving_distance_transverse(z, **cosmo)
    if z0 is 0:
        return dm2 / (1. + z), 0, err2

    dm1, err1 = comoving_distance_transverse(z0, **cosmo)

    d_h_0 = hubble_distance_z(0.0, **cosmo)

    term1 = dm1 * numpy.sqrt(1. + omega_k * (dm2/d_h_0)**2.)
    term2 = dm2 * numpy.sqrt(1. + omega_k * (dm1/d_h_0)**2.)

    da12 = (term2 - term1)/(1+z) # only for Omega_k > 0

    return da12, err1, err2

def luminosity_distance(z, **cosmo):
    """Calculate the luminosity distance to redshift z.
    
    Optionally calculate the integral from z0 to z.

    See, for example, David Hogg's arXiv:astro-ph/9905116v4

    """
    da, err1, err2 = angular_diameter_distance(z, **cosmo)
    return da * (1+z)**2., err1, err2


def _lookback_integral(z, omega_M_0, omega_lambda_0, omega_k_0, h):

    e_z = (omega_M_0 * (1+z)**3. + 
           omega_k_0 * (1+z)**2. + 
           omega_lambda_0)**0.5

    H_0 = h * cc.H100_s

    H_z =  H_0 * e_z
    
    return 1./((1. + z) * H_z)

def lookback_time(z, z0 = 0.0, **cosmo):
    """Calculate the lookback time (in s) to redshift z.

    See equation 30 of David Hogg's arXiv:astro-ph/9905116v4

    Optionally calculate the integral from z0 to z.

    Returns:
    -------

    t_look: ndarray
       Lookback time in seconds.

    err: ndarray 
       Estimate of numerical integration error in lookback time (in s)
       from scipy.integrate.quad.

    """

    #cosmo = set_omega_k_0(cosmo)

    lt_func = \
        numpy.vectorize(lambda z, z0, omega_M_0, omega_lambda_0, omega_k_0, h: 
                        si.quad(_lookback_integral, z0, z, limit=1000,
                                args=(omega_M_0, omega_lambda_0, omega_k_0, h)))
    t_look, err = lt_func(z, z0, 
                          cosmo['omega_M_0'],
                          cosmo['omega_lambda_0'],
                          cosmo['omega_k_0'],
                          cosmo['h']
                          )
    return t_look, err

def age(z, use_flat=True, **cosmo):
    """Calculate the age of the universe as seen at redshift z.

    Age at z is lookback time at z'->Infinity minus lookback time at z.
    
    See also: lookback_time.

    Returns two error estimates: one for the full age of the universe,
    the other for the lookback time.
    """
    if use_flat and get_omega_k_0(**cosmo) == 0:
        return age_flat(z, **cosmo), float('nan'), float('nan')
    fullage, err_f = lookback_time(numpy.Inf, **cosmo)
    tl, err_t = lookback_time(z, **cosmo)
    age = fullage - tl
    return age, err_f, err_t

def age_flat(z, **cosmo):
    """Calculate the age of the universe assuming a flat cosmology.

    Analytical formula from Peebles, p. 317, eq. 13.2.
    """

    omega_k = get_omega_k_0(**cosmo)
    if (numpy.any(omega_k != 0)):
        #raise ValueError("Not implemented for Omega_k != 0")
        print "Warning: using lambda = 1 - omega_M for non-flat cosmology!"

    om = cosmo['omega_M_0']
    lam = 1. - cosmo['omega_M_0']
    t_z = (2. * 
           numpy.arcsinh(numpy.sqrt(lam/om) * (1. + z)**(-3./2.)) /
           (cc.H100_s * cosmo['h'] * 3. * numpy.sqrt(lam)))

    return t_z

def quick_age_function(zmax = 20., zmin = 0., zstep = 0.1,
                       return_inverse=False,
                       **cosmo):
    """Return an interpolation function that will give age as a funtion of z

    If return_inverse is True, will also return a function giving z as
    a function of age.

    Returns
    -------

    agefunc, err_f, err_t

    or
    
    agefunc, redfunc, err_f, err_t
    
    """
    z = numpy.arange(zmin, zmax, zstep)
    ages, err_f, err_t = age(z, **cosmo)
    agefunc = scipy.interpolate.interp1d(z, ages)
    if return_inverse:
        redfunc = scipy.interpolate.interp1d(ages[::-1], z[::-1])
        return agefunc, redfunc, err_f, err_t
    else:
        return agefunc, err_f, err_t

def quick_redshift_age_function(zmax = 20., zmin = 0., zstep = 0.1, **cosmo):
    """Return an interpolation function that will give z as a funtion of
       the age of the universe.
    """
    z = numpy.arange(zmin, zmax, zstep)
    z = z[::-1]
    ages, err_f, err_t = age(z, **cosmo)
    return scipy.interpolate.interp1d(ages, z), err_f, err_t
    
def light_travel_distance(z, z0 = 0, **cosmo):
    """Calculate the light travel distance to redshift z.

    Units are Mpc.
    
    """
    t_look, err = lookback_time(z, z0, **cosmo)
    return cc.c_light_Mpc_s * t_look, cc.c_light_Mpc_s * err

def redshift_d_light(dl, z_guess = 6.0, **cosmo):
    """Calculate the redshift corresponding to a given light travel
    distance.

    Units are the same as light_travel_distance (Mpc).

    """
    
    dl_diff = lambda z: abs(dl - light_travel_distance(z, **cosmo)[0])
    z = scipy.optimize.fmin(dl_diff, z_guess)
    return z

if __name__ == "__main__":
    import doctest
    doctest.testmod()
