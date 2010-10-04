"""Test magnitude functions.
"""

import numpy

import cosmolopy.parameters as cparam
import cosmolopy.magnitudes as magnitudes
from cosmolopy.magnitudes import *
import cosmolopy.constants as cc

def test_magnitudes():
    cosmo = cparam.WMAP7_BAO_H0_mean(flat=True)

    z0 = 6.0
    mag0 = -18
    lambda0 = 1450.
    lambda_obs = lambda0 * (1.+z0)
    nu0 = cc.c_light_cm_s / (lambda0 * cc.angstrom_cm)
    nu_obs = nu0 / (1.+z0)

    ###
    
    lum_nu = L_nu_from_magAB(mag0)
    mag_abs3 = magnitude_AB_from_L_nu(lum_nu)
    assert abs(mag_abs3 - mag0) < abs(1e-10 * mag0)
    
    flux_nu = flux_from_luminosity(z0, luminosity=lum_nu, spec_unit='freq',
                                   **cosmo)
    flux_lambda = f_lambda_from_f_nu(flux_nu, nu_obs)
    flux_nu2 = f_nu_from_f_lambda(flux_lambda, lambda_obs)
    assert abs(flux_nu2 - flux_nu) < 1e-10 * flux_nu

    lum_nu2 = luminosity_from_flux(z0, flux=flux_nu2,
                                   spec_unit='freq', **cosmo)
    assert abs(lum_nu2 - lum_nu) < 1e-10 * lum_nu

    mag_obs, mag_abs = magnitude_AB(z0, flux_lambda,
                                    wavelength=lambda_obs, **cosmo)
    assert abs(mag_abs - mag0) < abs(1e-10 * mag0)

    mag_obs2, mag_abs2 = magnitude_AB1450(z0, flux_lambda,
                                          wavelength=lambda_obs,
                                          nu_power=-0.5, **cosmo)
    assert abs(mag_abs2 - mag0) < abs(1e-10 * mag0)
    
    flux_nu3 = f_nu_from_magAB(mag_obs2)
    assert abs(flux_nu3 - flux_nu) < 1e-10 * flux_nu

    flux_lambda3 = f_lambda_from_f_nu(flux_nu3, nu_obs)
    lum_lambda3 = luminosity_from_flux(z0, flux=flux_lambda3,
                                       spec_unit='wave', **cosmo)
    lum_nu3 = f_nu_from_f_lambda(lum_lambda3, lambda0)
    assert abs(lum_nu3 - lum_nu) < 1e-10 * lum_nu

    mu = distance_modulus(z0, **cosmo)
    lum_nu4 = L_nu_from_magAB(mag_obs - mu) / (1. + z0)
    assert abs(lum_nu4 - lum_nu) < 1e-10 * lum_nu

if __name__ == '__main__':
    test_magnitudes()

    
