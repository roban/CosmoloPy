"""Conversions between fluxes, luminosities and AB magnitudes.
"""
from __future__ import absolute_import, division, print_function

import math

import numpy
import cosmolopy.distance as cd
import cosmolopy.constants as cc

"""AB Magnitude zero point."""
MAB0 = -2.5 * numpy.log10(3631.e-23)

def nu_lambda(coordinate):
    """Convert between frequency and wavelength, nu to lambda or
    lambda to nu.

    Either:
     given `lambda` returns 'nu' or
     given `nu` returns `lambda`.

    Units are:
     `Hz` for nu and `Ang` for `lambda`.

    Works because `nu = c/lambda` and `lambda = c/nu`, and I use `c`
    in units of `Angs/s`.

    Usage
    -----

     >>> from cosmolopy import magnitudes
     >>> nu = magnitudes.nu_lambda(1216.)
     >>> lam = magnitudes.nu_lambda(nu)
     >>> lam
     1216.0
     """
    return (cc.c_light_cm_s / cc.angstrom_cm) / coordinate

def f_nu_lambda(flux, coordinate):
    """Convert f_nu to f_lambda or f_lambda to f_nu.

    Either:
     given `f_lambda` and `lambda` returns `f_nu` and 'nu' or
     given `f_nu` and `nu` returns `f_lambda` and `lambda`.

    Units are:
     `erg s^-1 cm^-2 Hz^-1` for f_nu and 
     `erg s^-1 cm^-2 Ang^-1` for `f_lambda`.

    Works because `f_nu = f_lambda * lambda**2/c` and `f_lambda = f_nu
    * nu**2/c`, and I use `c` in units of `Angs/s`.

    Usage
    -----

     >>> from cosmolopy import magnitudes
     >>> fnu, nu = magnitudes.f_nu_lambda(2.0, 1216.)
     >>> flam, lam = magnitudes.f_nu_lambda(fnu, nu)
     >>> flam, lam
     (2.0, 1216.0)
    """
    
    return (flux * coordinate**2. / (cc.c_light_cm_s / cc.angstrom_cm),
            (cc.c_light_cm_s / cc.angstrom_cm) / coordinate)


def f_nu_from_magAB(magAB):
    """Convert apparent magnitude into flux (erg s^-1 cm^-2 Hz^-1).

    Usage
    -----

    Check that the AB magnitude zero point is 3631 Jy:

     >>> from cosmolopy import magnitudes
     >>> "%.4g" % (magnitudes.f_nu_from_magAB(0.0)/1e-23)
     '3631'

    """
    f_nu = 10.**((magAB + MAB0)/(-2.5))
    return f_nu

def L_nu_from_magAB(magAB):
    """Convert absolute magnitude into luminosity (erg s^-1 Hz^-1).

    Usage
    -----

    Check that the AB magnitude zero point is 3631 Jy:

     >>> from cosmolopy import magnitudes
     >>> import math
     >>> L_nu = magnitudes.L_nu_from_magAB(0.0)
     >>> "%.4g" % (L_nu/(1e-23 * 4. * math.pi * (10*cc.pc_cm)**2))
     '3631'

    """
    const = 4. * math.pi * (10. * cc.pc_cm)**2.
    L_nu =  const * 10.**((magAB + MAB0)/(-2.5))
    return L_nu

def magnitude_AB_from_L_nu(luminosity_nu):
    """Convert luminosity (erg s^-1 Hz^-1) into absolute magnitude.

    Usage
    -----

    Check that the AB magnitude zero point is 3631 Jy:

     >>> import numpy, math
     >>> from cosmolopy import magnitudes, cc
     >>> L_nu = 3631e-23 * (4. * math.pi * (10*cc.pc_cm)**2)
     >>> "%.3f" % numpy.abs(magnitudes.magnitude_AB_from_L_nu(L_nu))
     '0.000'

    """

    const = 4. * math.pi * (10. * cc.pc_cm)**2.
    magAB = -2.5 * numpy.log10(luminosity_nu/const) - MAB0
    return magAB

def distance_modulus(z, **cosmo):
    """Distance modulus mu = m-M.

    The distance modulus is the difference between the apparent and
    absolute magnitudes,

      mu = 5 log(d/10 pc)

    Usage
    -----

    >>> from cosmolopy import fidcosmo, magnitudes
    >>> "mu(z=6) = %.4g" % magnitudes.distance_modulus(6.0, **fidcosmo)
    'mu(z=6) = 48.86'
    
    """
    dl = cd.luminosity_distance(z, **cosmo)
    mu = 5 * numpy.log10(dl/(10e-6))
    return mu

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
        the observed wavelength of the flux measurement(s) in Angstroms

    Returns
    -------

    Returns ab (apparent), and AB (absolute) magnitudes.

    Notes
    -----

    Note that here you pass fluxes that are per unit wavelength, not
    per unit frequency. To get the absolute magnitude for a
    *luminosity* specified in units of erg s^-1 Ang^-1, set z=None.

    Usage
    -----

    Check that the AB magnitude zero point is 3631 Jy:

     >>> from cosmolopy import fidcosmo, magnitudes, cc, cd
     >>> import numpy, math
     >>> L_nu = 3631e-23 * (4. * math.pi * (10*cc.pc_cm)**2)
     >>> nu = magnitudes.nu_lambda(1216.)
     >>> L_lambda, lamb = magnitudes.f_nu_lambda(L_nu, nu)
     >>> mAB, MAB = magnitudes.magnitude_AB(None, L_lambda, 1216., **fidcosmo)
     >>> "%.3f" % numpy.abs(MAB)
     '0.000'

    Find the apparent (and absolute, which should be zero) magnitudes
    of a 3631 Jy source at z=6.0:

     >>> from cosmolopy import fidcosmo, magnitudes, cc, cd
     >>> import numpy, math
     >>> L_nu = 3631e-23 * (4. * math.pi * (10*cc.pc_cm)**2)
     >>> nu = magnitudes.nu_lambda(1216.)
     >>> L_lambda, lamb = magnitudes.f_nu_lambda(L_nu, nu)
     >>> dl = cd.luminosity_distance(6.0, **fidcosmo)
     >>> f_lambda = L_lambda/(4. * math.pi * (dl*cc.Mpc_cm)**2 * (1. + 6.0))
     >>> mAB, MAB = magnitudes.magnitude_AB(6.0, f_lambda, 7.*1216., **fidcosmo)
     >>> "%.3f, %.3f" % (mAB, MAB)
     '48.865, 0.000'

    """

    # Distance modulus mu = m-M
    if z is None:
        mu = 0.0
        z = 0
        f_lambda = f_lambda / (4. * math.pi * (10. * cc.pc_cm)**2.)
    else:
        mu = distance_modulus(z, **cosmo)

    # Correction to the flux due to redshifted differential wavelength.
    f_rest = f_lambda * (1+z)
    
    # Rest wavelength and frequency.
    lambda_rest = wavelength/(1+z)
    nu_rest = cc.c_light_cm_s / (lambda_rest * cc.angstrom_cm)

    # Observed frequency.
    nu_0 = cc.c_light_cm_s / (wavelength * cc.angstrom_cm)

    # Apparent AB magnitude.
    ab_app = -2.5 * numpy.log10(f_rest * (lambda_rest / nu_rest)) - MAB0

    # Absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app, ab_abs

def magnitude_AB1450(z, f_lambda, wavelength, nu_power=-0.5, **cosmo):
    """Extrapolate to the AB magnitude at 1450 Angstroms.

    Inputs
    ------

    z: array or scalar
        the redshift of the source

    f_lambda: array or scalar
        observed flux from the source in units of erg s^-1 cm^-2 Ang^-1

    wavelength: array or scalar
        the observed wavelength of the flux measurement(s) in Angstroms.

    nu_power:
        the powerlaw index (f_nu ~ nu^nu_power) used to extrapolate
        the flux to 1450 Angstroms.

    Returns
    -------

    Apparent and absolute magnitudes extrapolated to 1450 Angstroms.


    Notes
    -----
    
    Follows Fan et al. 2003:

        We extrapolate the continuum to rest-frame 1450A, assuming a
        continuum shape f_nu ~ nu^-0.5 to calculate AB_1450.

    Usage
    -----

    Find the apparent and absolute rest-frame 1450 Angstrom magnitudes
    of source with a flux of 3631 Jy at rest-frame 1216 Angstroms at
    z=6.0:


     >>> from cosmolopy import fidcosmo, magnitudes, cc, cd
     >>> import numpy, math
     >>> L_nu = 3631e-23 * (4. * math.pi * (10*cc.pc_cm)**2)
     >>> nu = magnitudes.nu_lambda(1216.)
     >>> L_lambda, lamb = magnitudes.f_nu_lambda(L_nu, nu)
     >>> dl = cd.luminosity_distance(6.0, **fidcosmo)
     >>> f_lambda = L_lambda/(4. * math.pi * (dl*cc.Mpc_cm)**2 * (1. + 6.0))
     >>> mAB, MAB = magnitudes.magnitude_AB1450(6.0, f_lambda, 7.*1216., 
     ...                                        **fidcosmo)
     >>> "%.3f, %.3f" % (mAB, MAB)
     '48.769, -0.096'

    And is that offset from an absolute magnitude of zero consisten
    with our assumed powerlaw index?  

     >>> "%.3f" %(-2.5 * numpy.log10((1216./1450)**0.5))
     '0.096'

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
    ab_app = -2.5 * numpy.log10(f_nu_1450) - MAB0

    # distance modulus mu = m-M
    mu = distance_modulus(z, **cosmo)

    # absolute magnitude
    ab_abs = ab_app - mu 

    return ab_app, ab_abs

if __name__ == "__main__":
    import doctest
    doctest.testmod()
