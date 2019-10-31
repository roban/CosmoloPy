"""Cosmological densities like matter density, baryon density, etc.
"""

from __future__ import absolute_import, division, print_function

import math

import numpy
import scipy
import scipy.special
import scipy.integrate as si

from . import constants as cc
from . import distance as cd
from .distance import get_omega_k_0, set_omega_k_0

def omega_M_z(z, **cosmo):
    """Matter density omega_M as a function of redshift z.

    Notes
    -----

    From Lahav et al. (1991, MNRAS 251, 128) equations 11b-c. This is
    equivalent to equation 10 of Eisenstein & Hu (1999 ApJ 511 5).

    """
    if get_omega_k_0(**cosmo) == 0:
        return 1.0 / (1. + (1. - cosmo['omega_M_0'])/
                      (cosmo['omega_M_0'] * (1. + z)**3.))
    else:
        return (cosmo['omega_M_0'] * (1. + z)**3. / 
                cd.e_z(z, **cosmo)**2.)

def cosmo_densities(**cosmo):
    """The critical and mean densities of the universe.

    Returns
    -------
    rho_crit and rho_0 in solar masses per cubic Megaparsec.

    """

    omega_M_0 = cosmo['omega_M_0']
    h = cosmo['h']

    rho_crit = 3. * (h * cc.H100_s)**2. / (8. * math.pi * cc.G_const_Mpc_Msun_s)
    rho_0 = omega_M_0 * rho_crit
    
    #print " Critical density rho_crit = %.3g Msun/Mpc^3" % rho_crit
    #print " Matter density      rho_0 = %.3g Msun/Mpc^3" % rho_0

    return rho_crit, rho_0

def get_X_Y(**cosmo):
    """The fraction of baryonic mass in hydrogen and helium.

    Assumes X_H + Y_He = 1.

    You must specify either 'X_H', or 'Y_He', or both.
    """
    if 'X_H' in cosmo and 'Y_He' not in cosmo:
        X_H = cosmo['X_H']
        Y_He = 1. - X_H
    elif 'Y_He' in cosmo and 'X_H' not in cosmo:
        Y_He = cosmo['Y_He']
        X_H = 1. - Y_He
    else:
        X_H = cosmo['X_H']
        Y_He = cosmo['Y_He']
    return X_H, Y_He
    
def baryon_densities(**cosmo):
    """Hydrogen number density at z=0.

    Parameters
    ----------

    cosmo: cosmological parameters

    parameters used: 'omega_b_0', 'X_H' and/or 'Y_He', plus those
    needed by cosmo_densities.
       

    Returns
    -------

    rho_crit, rho_0, n_He_0, n_H_0

    The first two are in units of solar masses per cubic
    Megaparsec. The later two are in number per cubic Megaparsec.
    
    """
    
    X_H, Y_He = get_X_Y(**cosmo)

    rho_crit, rho_0 = cosmo_densities(**cosmo)

    n_H_0 = (rho_crit * cosmo['omega_b_0'] * X_H * cc.M_sun_g / 
             cc.m_H_g)
    n_He_0 = (rho_crit * cosmo['omega_b_0'] * Y_He * cc.M_sun_g / 
              cc.m_He_g)
    #    print " Hydrogen number density n_H_0 = %.4g (Mpc^-3)" % n_H_0
    return rho_crit, rho_0, n_He_0, n_H_0
