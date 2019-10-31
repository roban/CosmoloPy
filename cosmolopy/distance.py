"""Cosmological distance measures. 

Mostly follows David Hogg's pedagogical paper arXiv:astro-ph/9905116v4 .

Distance units are Mpc, time units are seconds.

"""

from __future__ import absolute_import, division, print_function

import math

import numpy
import scipy
import scipy.integrate as si
import scipy.interpolate
import scipy.optimize

from . import constants as cc

def get_omega_k_0(**cosmo):
    """'Spatial curvature density' omega_k_0 for a cosmology (if needed).

    If omega_k_0 is specified, return it. Otherwise return:

      1.0 - omega_M_0 - omega_lambda_0

    """
 
    if 'omega_k_0' in cosmo:
        omega_k_0 = cosmo['omega_k_0']
    else:
        omega_k_0 = 1. - cosmo['omega_M_0'] - cosmo['omega_lambda_0']
    return omega_k_0

def set_omega_k_0(cosmo):
    """Returns the cosmo dictionary with omega_k_0 set.
    See get_omega_k_0.
    
    Note that cosmo is not passed as \*\*cosmo for once. This function
    modifies the dictionary in place and returns the result.

    """
    if 'omega_k_0' in cosmo:
        return cosmo
    else:
        cosmo['omega_k_0'] = get_omega_k_0(**cosmo)
        return cosmo


### Distance measures ###

def e_z(z, **cosmo):
    """The unitless Hubble expansion rate at redshift z.

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to E(z), defined in his eq. 14.

    Modified (JBJ, 29-Feb-2012) to include scalar w parameter

    """

    if 'w' in cosmo:
        #ow = exp( si.quad( lambda zp: (1.+cosmo['w']) / (1.+zp), 0., z, limit=1000 ) )
        return (cosmo['omega_M_0'] * (1+z)**3. + 
                cosmo['omega_k_0'] * (1+z)**2. + 
                cosmo['omega_lambda_0'] * (1+z)**(1+cosmo['w']) )**0.5
    else:
        return (cosmo['omega_M_0'] * (1+z)**3. + 
                cosmo['omega_k_0'] * (1+z)**2. + 
                cosmo['omega_lambda_0'])**0.5

def hubble_z(z, **cosmo):
    """The value of the Hubble constant at redshift z.

    Units are s^-1

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to H_0 * E(z) (see his eq. 14).

    """
    H_0 = cosmo['h'] * cc.H100_s

    return H_0 * e_z(z, **cosmo)

def hubble_distance_z(z, **cosmo):
    """The value of the Hubble distance at redshift z.

    Units are Mpc.

    In David Hogg's (arXiv:astro-ph/9905116v4) formalism, this is
    equivalent to D_H / E(z) = c / (H_0 E(z)) [see his eq. 14], which
    appears in the definitions of many other distance measures.

    """
    H_0 = cosmo['h'] * cc.H100_s
    
    return cc.c_light_Mpc_s / (H_0 * e_z(z, **cosmo))

def _comoving_integrand(z, omega_M_0, omega_lambda_0, omega_k_0, h, w=-1.):

    e_z = (omega_M_0 * (1+z)**3. + 
           omega_k_0 * (1+z)**2. + 
           omega_lambda_0 * (1+z)**(1.+w))**0.5
    
    H_0 = h * cc.H100_s
    
    H_z =  H_0 * e_z

    return cc.c_light_Mpc_s / (H_z)

def comoving_integrand(z, **cosmo):
    """The derivative of the comoving distance with redshift: dd_c/dz.

    See equation 15 of David Hogg's arXiv:astro-ph/9905116v4

    Units are Mpc.
    
    """
    if 'w' in cosmo:
        w = cosmo['w']
    else:
        w = -1.

    return _comoving_integrand(z,
                               cosmo['omega_M_0'],
                               cosmo['omega_lambda_0'],
                               cosmo['omega_k_0'],
                               cosmo['h'], w)

def comoving_distance(z, z0 = 0, **cosmo):
    """The line-of-sight comoving distance (in Mpc) to redshift z.

    See equation 15 of David Hogg's arXiv:astro-ph/9905116v4

    Units are Mpc.

    Optionally calculate the integral from z0 to z.

    Returns
    -------
    
    d_co: ndarray
       Comoving distance in Mpc.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> d_co = cd.comoving_distance(6., **cosmo)
    >>> print "Comoving distance to z=6 is %.1f Mpc" % (d_co)
    Comoving distance to z=6 is 8017.8 Mpc

    """

    #cosmo = set_omega_k_0(cosmo)

    if 'w' in cosmo:
        w = cosmo['w']
    else:
        w = -1.

    dc_func = \
        numpy.vectorize(lambda z, z0, omega_M_0, omega_lambda_0, omega_k_0, h, w: 
                        si.quad(_comoving_integrand, z0, z, limit=1000,
                                args=(omega_M_0, omega_lambda_0, omega_k_0, h, w)))
    d_co, err = dc_func(z, z0, 
                        cosmo['omega_M_0'],
                        cosmo['omega_lambda_0'],
                        cosmo['omega_k_0'],
                        cosmo['h'],
                        w
                        )
    return d_co

def proper_motion_distance(z, **cosmo):
    """Returns comoving_distance_transverse.

    Units are Mpc.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> d_M = cd.proper_motion_distance(6., **cosmo)
    >>> print "Transverse comoving distance to z=6 is %.1f Mpc" % (d_M)
    Transverse comoving distance to z=6 is 8017.8 Mpc

    """
    return comoving_distance_transverse(z, **cosmo)

def comoving_distance_transverse(z, **cosmo):
    """The transverse comoving distance (in Mpc) to redshift z.

    This is also called the proper motion distance, D_M.

    See equation 16 of David Hogg's arXiv:astro-ph/9905116v4

    Units are Mpc.

    This is the distance d_m, such that the comoving distance between
    two events at the same redshift, but separated on the sky by some
    angle delta_theta is d_m * delta_theta.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> d_M = cd.comoving_distance_transverse(6., **cosmo)
    >>> print "Transverse comoving distance to z=6 is %.1f Mpc" % (d_M)
    Transverse comoving distance to z=6 is 8017.8 Mpc

    """

    d_c = comoving_distance(z, 0.0, **cosmo)

    omega_k_0 = get_omega_k_0(**cosmo)

    if numpy.all(omega_k_0 == 0.0):
        return d_c
    
    d_h_0 = hubble_distance_z(0.0, **cosmo)
    sqrt_ok0 = numpy.sqrt(numpy.abs(omega_k_0))
    if not numpy.isscalar(omega_k_0):
        sqrt_ok0[omega_k_0 == 0.0] = 1.0
    argument = sqrt_ok0 * d_c / d_h_0
    factor = d_h_0 * (1./sqrt_ok0)
    d_m = ((omega_k_0 > 0.0) * (factor * numpy.sinh(argument)) +
           (omega_k_0 == 0.0) * d_c +
           (omega_k_0 < 0.0) * (factor * numpy.sin(argument)))

    return d_m

def angular_diameter_distance(z, z0 = 0, **cosmo):
    """The angular-diameter distance (Mpc) to redshift z.
    
    Optionally find the angular diameter distance between objects at
    z0 and z (only implemented for omega_k_0 >= 0).

    See equations 18-19 of David Hogg's arXiv:astro-ph/9905116v4

    Units are Mpc.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> d_a = cd.angular_diameter_distance(6., **cosmo)
    >>> print "Angular diameter distance = %.1f Mpc" % (d_a)
    Angular diameter distance = 1145.4 Mpc

    """

    omega_k = numpy.atleast_1d(get_omega_k_0(**cosmo))
    if (numpy.any(omega_k < 0) and not(z0 == 0)):
        raise ValueError("Not implemented for Omega_k < 0 and z0 > 0")

    dm2  = comoving_distance_transverse(z, **cosmo)
    if z0 == 0:
        return dm2 / (1. + z)

    dm1 = comoving_distance_transverse(z0, **cosmo)

    d_h_0 = hubble_distance_z(0.0, **cosmo)

    term1 = dm1 * numpy.sqrt(1. + omega_k * (dm2/d_h_0)**2.)
    term2 = dm2 * numpy.sqrt(1. + omega_k * (dm1/d_h_0)**2.)

    da12 = (term2 - term1)/(1+z) # only for Omega_k > 0

    return da12

def luminosity_distance(z, **cosmo):
    """The luminosity distance to redshift z.
    
    Units are Mpc.

    See, for example, David Hogg's arXiv:astro-ph/9905116v4

    """
    da = angular_diameter_distance(z, **cosmo)
    return da * (1+z)**2.

def diff_comoving_volume(z, **cosmo):
    """The differential comoving volume element dV_c/dz/dSolidAngle.

    Dimensions are volume per unit redshift per unit solid angle.

    Units are Mpc**3 Steradians^-1.

    See David Hogg's arXiv:astro-ph/9905116v4, equation 28.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> dVc = cd.diff_comoving_volume(6.0, **cosmo)
    >>> print "dV/dz/dSolidAngle at z=6 is %.3g Mpc**3" % (dVc)
    dV/dz/dSolidAngle at z=6 is 2.63e+10 Mpc**3
    """
    
    d_h_0 = hubble_distance_z(0.0, **cosmo)
    d_m = comoving_distance_transverse(z, **cosmo)
    ez = e_z(z, **cosmo)
    return d_h_0 * d_m**2. / ez

def comoving_volume(z, **cosmo):
    """The comoving volume out to redshift z.

    See David Hogg's arXiv:astro-ph/9905116v4, equation 29.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> Vc = cd.comoving_volume(6.0, **cosmo)
    >>> print "Vc = %.3g Mpc**3" % (Vc)
    Vc = 2.16e+12 Mpc**3


    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.0, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> Vc = cd.comoving_volume(6.0, **cosmo)
    >>> print "Vc = %.3g Mpc**3" % (Vc)
    Vc = 1.68e+12 Mpc**3

    """

    dm = comoving_distance_transverse(z, **cosmo)

    omega_k_0 = get_omega_k_0(**cosmo)

    flat_volume = 4. * numpy.pi * dm**3. / 3.

    if numpy.all(omega_k_0 == 0.0):
        return flat_volume
    
    d_h_0 = hubble_distance_z(0.0, **cosmo)

    sqrt_ok0 = numpy.sqrt(numpy.abs(omega_k_0))
    dmdh = dm/d_h_0
    argument = sqrt_ok0 * dmdh
    f1 = 4. * numpy.pi * d_h_0**3. / (2. * omega_k_0)
    f2 = dmdh * numpy.sqrt(1. + omega_k_0 * (dmdh)**2.)
    f3 = 1./sqrt_ok0

    if numpy.isscalar(omega_k_0):
        if omega_k_0 > 0.0:
            return f1 * (f2 - f3 * numpy.arcsinh(argument))
        elif omega_k_0 == 0.0:
            return flat_volume
        elif omega_k_0 < 0.0:
            return f1 * (f2 - f3 * numpy.arcsin(argument))
    else:
        b = numpy.broadcast(omega_k_0,z,dm)
        Vc = numpy.zeros(b.shape)
        m1 = numpy.ones(b.shape, dtype=bool) * (omega_k_0 > 0.0)
        Vc[m1] = (f1 * (f2 - f3 * numpy.arcsinh(argument)))[m1]

        m1 = numpy.ones(b.shape, dtype=bool) * (omega_k_0 == 0.0)
        Vc[m1] = flat_volume[m1]

        m1 = numpy.ones(b.shape, dtype=bool) * (omega_k_0 < 0.0)
        Vc[m1] = (f1 * (f2 - f3 * numpy.arcsin(argument)))[m1]
        return Vc

def _lookback_integrand(z, omega_M_0, omega_lambda_0, omega_k_0, h):

    e_z = (omega_M_0 * (1+z)**3. + 
           omega_k_0 * (1+z)**2. + 
           omega_lambda_0)**0.5

    H_0 = h * cc.H100_s

    H_z =  H_0 * e_z
    
    return 1./((1. + z) * H_z)

def lookback_integrand(z, **cosmo):
    """The derivative of the lookback time with redshift: dt_L/dz.

    See equation 30 of David Hogg's arXiv:astro-ph/9905116v4

    Units are seconds.

    """
    return _lookback_integrand(z,
                               cosmo['omega_M_0'],
                               cosmo['omega_lambda_0'],
                               cosmo['omega_k_0'],
                               cosmo['h'])

def lookback_time(z, z0 = 0.0, **cosmo):
    """The lookback time (in s) to redshift z.

    See equation 30 of David Hogg's arXiv:astro-ph/9905116v4

    Units are s.

    Optionally calculate the integral from z0 to z.

    Returns
    -------

    t_look: ndarray
       Lookback time in seconds.

    """

    #cosmo = set_omega_k_0(cosmo)

    lt_func = \
        numpy.vectorize(lambda z, z0, omega_M_0, omega_lambda_0, omega_k_0, h: 
                        si.quad(_lookback_integrand, z0, z, limit=1000,
                                args=(omega_M_0, omega_lambda_0, omega_k_0, h)))
    t_look, err = lt_func(z, z0, 
                          cosmo['omega_M_0'],
                          cosmo['omega_lambda_0'],
                          cosmo['omega_k_0'],
                          cosmo['h']
                          )
    return t_look

def age(z, use_flat=True, **cosmo):
    """The age of the universe as seen at redshift z.

    Age at z is lookback time at z'->Infinity minus lookback time at z.
    
    See also: lookback_time.

    Units are s.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> t = cd.age(6.0, **cosmo)
    >>> print "age at z=6.0 = %.3g Gyr" % (t/cc.Gyr_s)
    age at z=6.0 = 0.892 Gyr

    """
    if use_flat and numpy.all(get_omega_k_0(**cosmo) == 0):
        return age_flat(z, **cosmo)
    fullage = lookback_time(numpy.Inf, **cosmo)
    tl = lookback_time(z, **cosmo)
    age = fullage - tl
    return age

def age_flat(z, **cosmo):
    """The age of the universe assuming a flat cosmology.
    
    Units are s.

    Analytical formula from Peebles, p. 317, eq. 13.2.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> t = cd.age_flat(6.0, **cosmo)
    >>> print "age at z=6.0 is %.3g Gyr" % (t/cc.Gyr_s)
    age at z=6.0 is 0.892 Gyr

    """

    omega_k = get_omega_k_0(**cosmo)
    if (numpy.any(omega_k != 0)):
        #raise ValueError("Not implemented for Omega_k != 0")
        print("Warning: using lambda = 1 - omega_M for non-flat cosmology!")

    om = cosmo['omega_M_0']
    lam = 1. - cosmo['omega_M_0']
    t_z = (2. * 
           numpy.arcsinh(numpy.sqrt(lam/om) * (1. + z)**(-3./2.)) /
           (cc.H100_s * cosmo['h'] * 3. * numpy.sqrt(lam)))

    return t_z

def quick_distance_function(function, zmax = 20., zmin = 0., zstep = 0.001,
                            return_inverse=False, k=3,
                            **cosmo):
    """Return an interpolation function that will give distance as a
    funtion of z

    If return_inverse is True, will also return a function giving z as
    a function of distance.

    Inputs
    ------

    function -- the distance function to interpolate (can be any
    callable that takes a redshift argument plus cosmology keywords).

    k -- spline order (`scipy.interpolate.InterpolatedUnivariateSpline`)

    Returns
    -------

    distfunc

    or
    
    distfunc, zfunc

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> distfunc, redfunc = cd.quick_distance_function(cd.luminosity_distance, return_inverse=True, **cosmo)
    >>> d = distfunc(6.3333)
    >>> z = redfunc(d)
    >>> "%.1g" % (distfunc(6.3333)/cd.luminosity_distance(6.3333, **cosmo) - 1.0)
    '-2e-16'
    >>> "%.1g" % (z/6.3333 - 1.0)
    '0'

    """
    z = numpy.linspace(zmin, zmax, math.ceil((zmax-zmin)/zstep))
    dists = function(z, **cosmo)
    distfunc = scipy.interpolate.InterpolatedUnivariateSpline(z, dists, k=k)
    if return_inverse:
        redfunc = scipy.interpolate.InterpolatedUnivariateSpline(dists, z, k=k)
        return distfunc, redfunc
    else:
        return distfunc

def quick_age_function(zmax = 20., zmin = 0., zstep = 0.001,
                       return_inverse=False,
                       **cosmo):
    """Return an interpolation function that will give age as a funtion of z

    Units are s.

    If return_inverse is True, will also return a function giving z as
    a function of age.

    Returns
    -------

    agefunc

    or
    
    agefunc, redfunc

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> agefunc = cd.quick_age_function(**cosmo)
    >>> t = agefunc(6.0)
    >>> print "age at z=6.0 is %.3g Gyr" % (t/cc.Gyr_s)
    age at z=6.0 is 0.892 Gyr

    
    """
    z = numpy.arange(zmin, zmax, zstep)
    ages = age(z, **cosmo)
    agefunc = scipy.interpolate.interp1d(z, ages)
    if return_inverse:
        redfunc = scipy.interpolate.interp1d(ages[::-1], z[::-1])
        return agefunc, redfunc
    else:
        return agefunc

def quick_redshift_age_function(zmax = 20., zmin = 0., zstep = 0.001, **cosmo):
    """Return an interpolation function giving z as a funtion of age
    of the universe.

    Units of time are s.

    Returns
    -------

    redfunc

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> redfunc = cd.quick_redshift_age_function(**cosmo)
    >>> z = redfunc(1.0 * cc.Gyr_s)
    >>> print "When age=1.0Gyr z=%.2f" % (z)
    When age=1.0Gyr z=5.49

    """
    z = numpy.arange(zmin, zmax, zstep)
    z = z[::-1]
    ages = age(z, **cosmo)
    return scipy.interpolate.interp1d(ages, z)
    
def light_travel_distance(z, z0 = 0, **cosmo):
    """The light travel distance to redshift z.

    Units are Mpc.

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> dlookback = cd.light_travel_distance(3.0, 2.0, **cosmo)
    >>> print "Lookback distance from z=2 to 3 is %.2g Mpc" % (dlookback)
    Lookback distance from z=2 to 3 is 3.3e+02 Mpc
    
    """
    t_look = lookback_time(z, z0, **cosmo)
    return cc.c_light_Mpc_s * t_look

def redshift_d_light(dl, z_guess = 6.0, fmin_args={}, **cosmo):
    """The redshift corresponding to a given light travel distance.

    Units are the same as light_travel_distance (Mpc).

    Examples
    --------

    >>> import cosmolopy.distance as cd
    >>> import cosmolopy.constants as cc
    >>> cosmo = {'omega_M_0' : 0.3, 'omega_lambda_0' : 0.7, 'h' : 0.72}
    >>> cosmo = cd.set_omega_k_0(cosmo)
    >>> z = cd.redshift_d_light(10. * cc.c_light_Mpc_Gyr, **cosmo)
    Optimization terminated successfully.
             Current function value: 0.000112
             Iterations: 26
             Function evaluations: 52
    >>> print "Redshift for a lookback time of 10Gyr is z=%.3f" % (z)
    Redshift for a lookback time of 10Gyr is z=2.025

    """
    
    dl_diff = lambda z: abs(dl - light_travel_distance(z, **cosmo)[0])
    z = scipy.optimize.fmin(dl_diff, z_guess, **fmin_args)
    return z

if __name__ == "__main__":
    import doctest
    doctest.testmod()
