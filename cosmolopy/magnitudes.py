"""Conversions between fluxes, luminosities and AB magnitudes.
"""
import math

import numpy
import cosmolopy.distance as cd
import cosmolopy.constants as cc

def f_nu_from_magAB(magAB):
    """Convert apparent magnitude into flux (erg s^-1 cm^-2 Hz^-1)."""
    f_nu = 10.**((magAB + 48.6)/(-2.5))
    return f_nu

def L_nu_from_magAB(magAB):
    """Convert absolute magnitude into luminosity (erg s^-1 Hz^-1)."""
    const = 4. * math.pi * (10. * cc.pc_cm)**2.
    L_nu =  const * 10.**((magAB + 48.6)/(-2.5))
    return L_nu

def magnitude_AB_from_L_nu(luminosity_nu):
    """Convert luminosity (erg s^-1 Hz^-1) into absolute magnitude."""

    const = 4. * math.pi * (10. * cc.pc_cm)**2.
    magAB = -2.5 * numpy.log10(luminosity_nu/const) - 48.6
    return magAB

def magnitude_AB(z, f_lambda, wavelength, **cosmo):
    """The apparent and absolute AB magnitude given a flux.

    Inputs
    ------

    z: array or scalar
        the redshift of the source. Set to None to get absolute
        magnitude from a luminosity.

    f_lambda: array or scalar
        observed flux from the source in units of erg s^-1 cm^-2 Ang^-1

    wavelength: array or scalar
        the observed wavelength of the flux measurement(s) in angstroms

    Returns
    -------

    Returns ab (apparent), and AB (absolute) magnitudes.

    """
    
    # Correction to the flux due to redshifted differential wavelength.
    f_rest = f_lambda * (1+z)
    
    # Rest wavelength and frequency.
    lambda_rest = wavelength/(1+z)
    nu_rest = cc.c_light_cm_s / (lambda_rest * cc.angstrom_cm)

    # Observed frequency.
    nu_0 = cc.c_light_cm_s / (wavelength * cc.angstrom_cm)

    ab_app_nok = -2.5 * numpy.log10(f_lambda * (wavelength / nu_0)) - 48.60
    #print "apparent AB mag without K correction = %.4g" % ab_app_nok

    # Apparent AB magnitude.
    ab_app = -2.5 * numpy.log10(f_rest * (lambda_rest / nu_rest)) - 48.60

    # Distance modulus mu = m-M
    if z is None:
        mu = -5.0
    else:
        dl = cd.luminosity_distance(z, **cosmo)[0]
        mu = 5 * numpy.log10(dl/(10*cc.pc_cm))

    # Absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app, ab_abs

def magnitude_AB1450(z, f_lambda, wavelength, nu_power=-0.5, **cosmo):
    """Extrapolate to the AB magnitude at 1450 cc.angstrom_cms.

    Inputs
    ------

    z: array or scalar
        the redshift of the source

    f_lambda: array or scalar
        observed flux from the source in units of erg s^-1 cm^-2 Ang^-1

    wavelength: array or scalar
        the observed wavelength of the flux measurement(s) in cc.angstrom_cms

    nu_power:
        the powerlaw index (f_nu ~ nu^nu_power) used to extrapolate
        the flux to 1450 anstroms.

    Returns
    -------

    Apparent and absolute magnitudes extrapolated to 1450 angstroms.


    Notes
    -----
    
    Follows Fan et al. 2003:

        We extrapolate the continuum to rest-frame 1450A, assuming a
        continuum shape f_nu ~ nu^-0.5 to calculate AB_1450.

    """
    
    # correction to the flux due to redshifted differential wavelength
    f_rest = f_lambda * (1+z)

    # rest wavelength and frequency
    lambda_rest = wavelength/(1+z)
    nu_rest = cc.c_light_cm_s / (lambda_rest * cc.angstrom_cm)

    # rest flux per unit freq.
    f_nu_rest = f_rest * (lambda_rest / nu_rest)

    nu_1450 = cc.c_light_cm_s / (1450 * cc.angstrom_cm)
    f_nu_1450 = f_nu_rest * (nu_1450/nu_rest)**nu_power
    
    # apparent AB magnitude
    ab_app = -2.5 * numpy.log10(f_nu_1450) - 48.59

    # distance modulus mu = m-M
    dl = cd.luminosity_distance(z, **cosmo)[0]
    mu = 5 * numpy.log10(dl/(10*cc.pc_cm))

    # absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app, ab_abs
