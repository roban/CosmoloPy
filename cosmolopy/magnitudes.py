"""Simple routines for conversion in and out of the AB magnitude system.

Also contains some functions for luminosity and flux conversions.
"""
import math

import numpy
import cosmolopy.distance as cd
import cosmolopy.constants as cc

def distance_modulus(z, **cosmo):
    """The distance modulus mu = m-M to redshift z.
    
    """
    dl = cd.luminosity_distance(z, **cosmo)[0]
    mu = 5. * numpy.log10((dl * 1e6)/10)
    return mu

def distance_from_modulus(mu, **cosmo):
    """Calculate the distance (in Mpc) given the distance modulus mu = m-M.
    
    """
    return (10**(0.2 * mu + 1))/1e6

def Kfactor(z, spec_unit):
    """Factor to correct for redshift distorting the spectral dimension.

    This is the equivalent of the 'K-correction' for a dispersed
    spectrum that is always observed at the same source-frame frequency
    or wavelength.

    If frequency units are used, the observed spectrum is squashed, so
    the source-frame spectrum would be lower by a factor of
    1/(1+z). If wavelength units are used, the observed spectrum is
    strechted, so the source-frame spectrum would be higher by a
    factor of 1+z.

    If the width of the integrated region is adjusted to correct for
    the redshift distortion, then no correction is needed. This may be
    the case for a measurement of the total flux in an emission line,
    for instance (though care must be taken with the definition of the
    total flux and the integration range).

    Parameters
    ----------
    
    spec_unit: string 'freq', 'wave', or None to specify flux per unit
        frequency, per unit wavelength, or total.
    
    """
    Kfactor = None
    if spec_unit == 'freq':
        Kfactor = 1./(1. + z)
    elif spec_unit == 'wave':
        Kfactor = 1. + z
    elif spec_unit is None:
        Kfactor = 1.
    else:
        raise ValueError, "Kfactor must be \'freq\', \'wave\' or None."
    return Kfactor

def luminosity_from_flux(z, flux=1.0, spec_unit='wave', **cosmo):
    """Convert a flux at redshift z to a luminosity.

    Parameters
    ----------
    
    spec_unit: see `Kfactor` function.
    
    """
    dl = cd.luminosity_distance(z, **cosmo)[0]
    luminosity = flux * Kfactor(z, spec_unit) * 4. * math.pi * (dl * cc.Mpc_cm)**2.
    return luminosity

def flux_from_luminosity(z, luminosity=1.0, spec_unit='wave', **cosmo):
    """Convert a luminosity at redshift z to a flux.

    Parameters
    ----------
    
    spec_unit: see `Kfactor` function.
    """
    dl = cd.luminosity_distance(z, **cosmo)[0]
    flux = (luminosity /
            (Kfactor(z, spec_unit) * 4. * math.pi * (dl * cc.Mpc_cm)**2.))
    return flux

def f_lambda_from_f_nu(flux_nu, frequency):
    """Convert from flux (or luminosity) in Hz^-1 to angstroms^-1.
    """
    return flux_nu * frequency**2. / (cc.c_light_cm_s / cc.angstrom_cm)

def f_nu_from_f_lambda(flux_lambda, wavelength):
    """Convert from flux (or luminosity) in angstroms^-1 to Hz^-1.
    """
    return flux_lambda * wavelength**2. / (cc.c_light_cm_s / cc.angstrom_cm)

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
    """Calculate the apparent and absolute AB magnitude given a flux.

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

    # Apparent AB magnitude.
    ab_app = -2.5 * numpy.log10(f_rest * (lambda_rest / nu_rest)) - 48.60

    # Distance modulus mu = m-M
    if z is None:
        mu = -5.0
    else:
        mu = distance_modulus(z, **cosmo)

    # Absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app_nok, ab_abs

def magnitude_AB1450(z, f_lambda, wavelength, nu_power=-0.5, **cosmo):
    """Extrapolate to the AB magnitude at 1450 angstroms.

    Assumes a powerlaw spectrum with index nu_power.

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
    f_rest = f_lambda * (1.+z)

    # rest wavelength and frequency
    lambda_rest = wavelength/(1+z)
    nu_rest = cc.c_light_cm_s / (lambda_rest * cc.angstrom_cm)

    # rest flux per unit freq.
    f_nu_rest = f_rest * (lambda_rest / nu_rest)

    nu_1450 = cc.c_light_cm_s / (1450 * cc.angstrom_cm)
    f_nu_1450 = f_nu_rest * (nu_1450/nu_rest)**nu_power

    # apparent AB magnitude
    ab_app = -2.5 * numpy.log10(f_nu_1450) - 48.6
    ab_app_nok = -2.5 * numpy.log10(f_nu_1450 * (1. + z)) - 48.6

    # distance modulus mu = m-M
    mu = distance_modulus(z, **cosmo)

    # absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app_nok, ab_abs
