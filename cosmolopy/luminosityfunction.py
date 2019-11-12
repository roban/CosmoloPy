"""Galaxy luminosity functions (Schechter functions).

The `LFHistory` class implements a luminosity function history,
encapsulating the changes in the galaxy luminosity distribution
parameters as a function of redshift.

"""

from __future__ import absolute_import, division, print_function

import os
import optparse

import numpy
import scipy.special

import cosmolopy.reionization as cr
import cosmolopy.distance as cd
import cosmolopy.parameters as cp
import cosmolopy.constants as cc
import cosmolopy.utils as utils
from cosmolopy.saveable import Saveable
import cosmolopy.magnitudes as magnitudes

def mass_from_sfr(sfr):
    """Use Labbe et al. (2009) relation between stellar mass and SFR.
    
    See arXiv:0911.1356v4

    Parameters
    ----------

    sfr:
        the star formation rate in M_Sun / year.

    Returns
    -------

    mass:
        stellar mass in units of M_Sun
    """
    return 10**(8.7 + 1.06 * numpy.log10(sfr))

def sfr_from_mass(mass):
    """Use Labbe et al. (2009) relation between stellar mass and SFR.
    
    See arXiv:0911.1356v4

    Parameters
    ----------

    mass:
        stellar mass in units of M_Sun

    Returns
    -------

    sfr:
        the star formation rate in M_Sun / year.

    """
    return 10**((numpy.log10(mass) - 8.7)/1.06)

def sfr_from_L_nu(luminosity):
    """Use Kennicutt (1998) conversion from UV luminosity to star formation rate.

    Parameters
    ----------

    luminosity:
        the luminosity in units of ergs s^-1 Hz^-1 anywhere between
        1500-2800 Angstroms.

    Returns
    -------

    The SFR in Msun/year.

    Notes
    -----

    Kennicutt (1998ARA&A..36..189K) says:

       SFR/(MSun/year) = 1.4 * 10^-28 (L_nu/ergs s^-1 Hz^-1)

    where L_nu is the UV luminosity anywhere between 1500-2800 Angstroms.
    """
    return 1.4e-28 * luminosity

def L_nu_from_sfr(sfr):
    """Use Kennicutt (1998) conversion from UV luminosity to star formation rate.

    Parameters
    ----------

    sfr:
        The SFR in Msun/year.

    Returns
    -------

    luminosity:
        the luminosity in units of ergs s^-1 Hz^-1 anywhere between
        1500-2800 Angstroms.

    Notes
    -----

    Kennicutt (1998ARA&A..36..189K) says:

       SFR/(MSun/year) = 1.4 * 10^-28 (L_nu/ergs s^-1 Hz^-1)

    where L_nu is the UV luminosity anywhere between 1500-2800 Angstroms.
    """
    return sfr/(1.4e-28)

def magnitudeAB_from_sfr(sfr):
    """Use Kennicutt (1998) conversion from UV luminosity to AB magnitude.

    Convenience function: uses L_nu_from_sfr and
    magnitudes.magnitude_AB_from_L_nu.
    """
    lnu = L_nu_from_sfr(sfr)
    return magnitudes.magnitude_AB_from_L_nu(lnu)

def schechterL(luminosity, phiStar, alpha, LStar):
    """Schechter luminosity function."""
    LOverLStar = (luminosity/LStar)
    return (phiStar/LStar) * LOverLStar**alpha * numpy.exp(- LOverLStar)

def schechterM(magnitude, phiStar, alpha, MStar):
    """Schechter luminosity function by magnitudes."""
    MStarMinM = 0.4 * (MStar - magnitude)
    return (0.4 * numpy.log(10) * phiStar *
            10.0**(MStarMinM * (alpha + 1.)) * numpy.exp(-10.**MStarMinM))

def schechterCumuLL(luminosity, phiStar, alpha, LStar):
    """Integrate luminosity in galaxies above luminosity=L.

    Uses an analytical formula.
    """
    # Note that the scipy.special definition of incomplete gamma is
    # normalized to one and is the lower incomplete gamma function, so
    # we have to use the complement and multiply by the unnormalized
    # integral (which is computed in schechterTotLL).
    return (schechterTotLL(phiStar, alpha, LStar) *
            (1. - scipy.special.gammainc(alpha+2., luminosity/LStar)))

def schechterCumuLM(magnitudeAB, phiStar, alpha, MStar):
    """Integrate luminosity in galaxies brighter than magnitudeAB.

    Uses an analytical formula.
    """
    LStar = magnitudes.L_nu_from_magAB(MStar)
    lum = magnitudes.L_nu_from_magAB(magnitudeAB)
    return schechterCumuLL(lum, phiStar, alpha, LStar)

def schechterTotLL(phiStar, alpha, LStar):
    """Integrate total luminosity in galaxies.

    Uses an analytical formula.
    """
    return phiStar * LStar * scipy.special.gamma(alpha + 2.)

def schechterTotLM(phiStar, alpha, MStar):
    """Integrate total luminosity in galaxies.

    Uses an analytical formula.
    """
    LStar = magnitudes.L_nu_from_magAB(MStar)
    return schechterTotLL(phiStar, alpha, LStar)

def iPhotonRateDensity(schechterParams,
                       maglim=None,
                       sedParams={},
                       wavelength=1500.):
    """Ionizing photon rate density from a luminosity function.

    in units of photons s^-1.

    Given schecterParams, the parameters of a Schechter luminosity
    function (in terms of AB Magnitudes), sedParams, the parameters of
    the galactic Spectral Energy Distribution, and the wavelength of
    the AB Magnitudes, calculate the emission rate density of ionizing
    photons.

    See Also
    --------

    BrokenPowerlawSED

    schechterTotLM

    """
    if maglim is None:
        lum = schechterTotLM(**schechterParams)
    else:
        lum = schechterCumuLM(maglim, **schechterParams)
    rQL = BrokenPowerlawSED(**sedParams).iPhotonRateRatio(wavelength)
    return lum * rQL

# From Bouwens et al. (2007)
B2007_z4 = {'MStar': -20.98,
            'phiStar': 1.3e-3,
            'alpha': -1.73,
            'z':3.8}

B2007_z5 = {'MStar': -20.64,
            'phiStar': 1.0e-3,
            'alpha': -1.66,
            'z':5.0}

B2007_z6 = {'MStar': -20.24,
            'phiStar': 1.4e-3,
            'alpha': -1.74,
            'z':5.9}

B2007_z7 = {'MStar': -19.3,
            'phiStar': 1.4e-3,
            'alpha': -1.74,
            'z':7.4
            }

# From Oesch et al. (2010)
O2010_z7 = {'MStar': -19.9,
            'phiStar': 1.4e-3,
            'alpha': -1.77,
            'z':6.8 # Not sure about this
            }

# From Bouwens 2008ApJ...686..230B
B2008 = {'z' : [3.8, 5.0, 5.9, 7.3, 9.0],
         'MStar': [-20.98, -20.64, -20.24, -19.8, -19.6],
         'phiStar': [1.3e-3, 1.0e-3, 1.4e-3, 1.1e-3, 1.1e-3],
         'alpha': [-1.73, -1.66, -1.74, -1.74, -1.74]}

# Like B2008, but with fixed phiStar, alpha
B2008_fixed = {'z' : [3.8, 5.0, 5.9, 7.3, 9.0],
               'MStar': [-20.98, -20.64, -20.24, -19.8, -19.6],
               'phiStar': [1.1e-3, 1.1e-3, 1.1e-3, 1.1e-3, 1.1e-3],
               'alpha': [-1.74, -1.74, -1.74, -1.74, -1.74]}

# Like B2008, but with fixed phiStar, alpha
zlinear = numpy.array([3.8, 5.0, 5.9, 7.3, 9.0])
B2008_linear_z5 = {'z' : zlinear,
                   'MStar': 0.36 * (zlinear - 5.0) - 20.64,
                   'phiStar': [1.1e-3, 1.1e-3, 1.1e-3, 1.1e-3, 1.1e-3],
                   'alpha': [-1.74, -1.74, -1.74, -1.74, -1.74]}

# From Trenti et al. 2010, EXCEPT the z=9 alpha value, which I lowered.
T2010ICLF = {'z' : [4., 5., 7., 8., 9.],
             'MStar' : [-20.90, -20.55, -20.00, -19.70, -19.55],
             'phiStar' : [1.3e-3, 1.5e-3, 1.0e-3, 0.60e-3, 0.22e-3],
             'alpha' : [-1.57, -1.63, -1.84, -1.90, -1.999]}

class LFHistory(Saveable):
    """Interpolate / extrapolate the Schechter parameters.
    
    By default, above the observed redshift range:
    
    MStar is linearly extrapolated as a function of time (not z) to high z.

    phiStar and alpha are constant at high z.
    
    """

    def __init__(self, params=B2008,
                 MStar_bounds = ['extrapolate', float('NaN')],
                 phiStar_bounds = ['constant', float('NaN')],
                 alpha_bounds = ['constant', float('NaN')],
                 extrap_args = {},
                 extrap_var = 'z',
                 sedParams = {},
                 wavelength = 1500.,
                 **cosmo):
        
        for (k, v) in params.items():
            params[k] = numpy.asarray(v)
            
        self.params = params
        self.MStar_bounds = MStar_bounds
        self.phiStar_bounds = phiStar_bounds
        self.alpha_bounds = alpha_bounds
        self.extrap_args = extrap_args
        self.extrap_var = extrap_var
        self.sedParams = sedParams
        self.wavelength = wavelength
        self.cosmo = cosmo

        self.zobs = params['z']
        self.tobs = cd.age(self.zobs, **cosmo)[0]
        self.MStar = params['MStar']
        self.phiStar = params['phiStar']
        self.alpha = params['alpha']

        if extrap_var == 't':
            self.xobs = self.tobs
            self._iPhotFunc = \
                       lambda t1, mag: (self.iPhotonRateDensity_t(t1,
                                                                  maglim=mag))
            
        elif extrap_var == 'z':
            self.xobs = self.zobs
            MStar_bounds = MStar_bounds[::-1]
            phiStar_bounds = phiStar_bounds[::-1]
            alpha_bounds = alpha_bounds[::-1]
            self._iPhotFunc = \
                        lambda z1, mag: (self.iPhotonRateDensity_z(z1,
                                                                   maglim=mag))
            
        self._MStarfunc = utils.Extrapolate1d(self.xobs, self.MStar,
                                              bounds_behavior=MStar_bounds,
                                              **extrap_args
                                              )
        print("M*:", end=' ')
        print(self._MStarfunc.extrap_string())

        self._phiStarfunc = utils.Extrapolate1d(self.xobs, self.phiStar,
                                                bounds_behavior=phiStar_bounds,
                                                **extrap_args
                                                )
        print("phi*:", end=' ')
        print(self._phiStarfunc.extrap_string())
        self._alphafunc = utils.Extrapolate1d(self.xobs, self.alpha,
                                              bounds_behavior=alpha_bounds,
                                              **extrap_args
                                              )
        print("alpha:", end=' ')
        print(self._alphafunc.extrap_string())
        
        self._SED = BrokenPowerlawSED(**sedParams)
        self._rQL = self._SED.iPhotonRateRatio(wavelength)

        for name, func in list(globals().items()):
            if not name.startswith('schechter'):
                continue
            def newfunc(z, _name=name, _func=func, **kwargs):
                params = self.params_z(z)
                M = params['MStar']
                phi = params['phiStar']
                alpha = params['alpha']
                return _func(MStar=M, phiStar=phi, alpha=alpha, **kwargs)
            newfunc.__name__ = func.__name__
            newfunc.__doc__ = func.__doc__
            self.__dict__[name] = newfunc

    def iPhotonRateDensity_z(self, z, maglim=None):
        """Ionizing photon rate density from a luminosity function.

        See the iPhotonRateRatio function.
        """
        params = self.params_z(z)
        if maglim is None:
            lum = schechterTotLM(**params)
        else:
            lum = schechterCumuLM(maglim, **params)
        return lum * self._rQL

    def iPhotonRateDensity_t(self, t, maglim=None):
        """Ionizing photon rate density from a luminosity function.
        
        See the iPhotonRateRatio function.
        """
        params = self.params_t(t)
        if maglim is None:
            lum = schechterTotLM(**params)
        else:
            lum = schechterCumuLM(maglim, **params)
        return lum * self._rQL

    def ionization(self, z, maglim=None):
        xH = cr.ionization_from_luminosity(z,
                                           self._iPhotFunc,
                                           rate_is_tfunc =
                                           self.extrap_var == 't',
                                           ratedensityfunc_args=[maglim],
                                           **self.cosmo)
        return xH

    def params_t(self, t):
        """Return interp/extrapolated Schechter function parameters."""
        if self.extrap_var == 't':
            return {'MStar':self._MStarfunc(t),
                    'phiStar':self._phiStarfunc(t),
                    'alpha':self._alphafunc(t)}
        elif self.extrap_var == 'z':
            ### FIX THIS ###
            raise NotImplementedError("params_t not implemented for z interps!")

    def params_z(self, z):
        """Return interp/extrapolated Schechter function parameters."""
        if self.extrap_var == 'z':
            return {'MStar':self._MStarfunc(z),
                    'phiStar':self._phiStarfunc(z),
                    'alpha':self._alphafunc(z)}
        elif self.extrap_var == 't':
            z = numpy.atleast_1d(z)
            t = cd.age(z, **self.cosmo)[0]
            return self.params_t(t)

def plotLFevo(hist=None,
              params=B2008,
              extrap_var='t',
              #maglim = -21.07 - 2.5 * numpy.log10(0.2)):
              maglim = -21. - 2.5 * numpy.log10(0.2),
              z_max=20.0,
              skipIon=True
              ):

    """Plot evolution of luminosity function params and total luminsity.

    Schechter function at each redshift is integrated up to maglim to
    find total luminsity.
    """

    for (k, v) in params.items():
        params[k] = numpy.asarray(v)

    if hist is None:
        cosmo = cp.WMAP7_BAO_H0_mean(flat=True)
        hist = LFHistory(params, extrap_var=extrap_var, **cosmo)
    else:
        params = hist.params
        extrap_var = hist.extrap_var
        cosmo = hist.cosmo

    z = params['z']
    t = cd.age(z, **cosmo)[0] / cc.yr_s 
    MStar = params['MStar']
    phiStar = params['phiStar']
    alpha = params['alpha']

    if maglim is None:
        ltot = schechterTotLM(MStar=MStar, phiStar=phiStar, alpha=alpha)
    else:
        ltot = schechterCumuLM(magnitudeAB=maglim,
                               MStar=MStar, phiStar=phiStar, alpha=alpha)

    print(hist._MStarfunc.extrap_string())

    zPlot = numpy.arange(z.min()-0.1, z_max, 0.1)
    tPlot = cd.age(zPlot, **cosmo)[0] / cc.yr_s
    newparams = hist.params_z(zPlot)
    MPlot = newparams['MStar']
    phiPlot = newparams['phiStar']
    alphaPlot = newparams['alpha']

    if maglim is None:
        ltotPlot = hist.schechterTotLM(zPlot)
    else:
        ltotPlot = hist.schechterCumuLM(zPlot, magnitudeAB=maglim)
    # From Table 6 of 2008ApJ...686..230B
    lB2008 =10.** numpy.array([26.18, 25.85, 25.72, 25.32, 25.14])

    iPhot = hist.iPhotonRateDensity_z(zPlot, maglim=maglim)
#    iPhotFunc = lambda t1: cc.yr_s * hist.iPhotonRateDensity_t(t1,
#                                                               maglim=maglim)

    if not skipIon:
        xH = hist.ionization(zPlot, maglim)
    import pylab
    pylab.figure(1)
    pylab.gcf().set_label('LFion_vs_z')
    if skipIon:
        pylab.subplot(111)
    else:
        pylab.subplot(211)
    pylab.plot(zPlot, iPhot)
    if not skipIon:
        pylab.subplot(212)
        pylab.plot(zPlot, xH)
        pylab.ylim(0.0, 1.5)
        
    pylab.figure(2)
    pylab.gcf().set_label('LFparams_vs_z')
    pylab.subplot(311)
    pylab.plot(z, MStar, '-')
    pylab.plot(z, MStar, 'o')
    pylab.plot(zPlot, MPlot, ':')

    pylab.subplot(312)
    pylab.plot(z, phiStar, '-')
    pylab.plot(z, phiStar, 'o')
    pylab.plot(zPlot, phiPlot, ':')

    pylab.subplot(313)
    pylab.plot(z, alpha, '-')
    pylab.plot(z, alpha, 'o')
    pylab.plot(zPlot, alphaPlot, ':')

    pylab.figure(3)
    pylab.gcf().set_label('LFparams_vs_t')
    pylab.subplot(311)
    pylab.plot(t, MStar, '-')
    pylab.plot(t, MStar, 'o')
    pylab.plot(tPlot, MPlot, ':')

    pylab.subplot(312)
    pylab.plot(t, phiStar, '-')
    pylab.plot(t, phiStar, '.')
    pylab.plot(tPlot, phiPlot, ':')

    pylab.subplot(313)
    pylab.plot(t, alpha, '-')
    pylab.plot(t, alpha, 'o')
    pylab.plot(tPlot, alphaPlot, ':')

    pylab.figure(4)
    pylab.gcf().set_label('LFlum_vs_z')
    pylab.subplot(121)
    pylab.plot(z, ltot, 'o')
    pylab.plot(z, lB2008, 'x')
    pylab.plot(zPlot, ltotPlot)

    pylab.subplot(122)
    pylab.plot(t, ltot, 'o')
    pylab.plot(t, lB2008, 'x')
    pylab.plot(tPlot, ltotPlot)

def test_plot_schechter():

    phiStar = 1.8e-3
    MStar = -20.04
    alpha = -1.71
    
    LStar = magnitudes.L_nu_from_magAB(MStar)

    mags = numpy.arange(-22, -11.0, 0.5)
    lums = magnitudes.L_nu_from_magAB(mags)

    phi_L = schechterL(lums, phiStar, alpha, LStar)
    phi_M = schechterM(mags, phiStar, alpha, MStar)

    L_L = schechterCumuLL(lums, phiStar, alpha, LStar)
    L_M = schechterCumuLM(mags, phiStar, alpha, MStar)

    phi_L_func = lambda l: l * schechterL(l, phiStar, alpha, LStar)
    L_L_num = utils.logquad(phi_L_func, lums, 1e35)[0]
    L_L_num2 = utils.vecquad(phi_L_func, lums, 1e29)[0]
    
    phi_M_func = lambda m: (magnitudes.L_nu_from_magAB(m) *
                            schechterM(m, phiStar, alpha, MStar))
    L_M_num2 = utils.vecquad(phi_M_func, -25, mags)[0]

    Ltot_L = schechterTotLL(phiStar, alpha, LStar)
    Ltot_M = schechterTotLM(phiStar, alpha, MStar)

    pylab.figure()
    pylab.subplot(221)
    pylab.plot(lums, lums * lums * phi_L)
    pylab.xscale('log')
    pylab.yscale('log')
    pylab.ylabel(r'$ L^2 \Phi_L$')

    pylab.subplot(222)
    pylab.plot(mags, -mags * lums * phi_M)
    pylab.yscale('log')
    pylab.ylabel(r'$ -M L \Phi_M$')

    pylab.subplot(223)
    pylab.plot(lums, Ltot_L - L_L)
    pylab.plot(lums, L_M)
    pylab.plot(lums, L_L_num, '--')
    pylab.plot(lums, L_L_num2, ':')
    pylab.plot(lums, L_M_num2, 'x')
    pylab.axhline(y=Ltot_L)
    pylab.axhline(y=Ltot_M)
    pylab.xscale('log')
    pylab.yscale('log')

    pylab.subplot(224)
    pylab.plot(mags, Ltot_M - L_M)
    pylab.plot(mags, L_L)
    pylab.plot(mags, L_L_num, '--')
    pylab.plot(mags, L_L_num2, ':')
    pylab.plot(mags, L_M_num2, 'x')
    pylab.axhline(y=Ltot_L)
    pylab.axhline(y=Ltot_M)
    pylab.yscale('log')

    #pylab.show()

class BrokenPowerlawSED(object):
    """Define an SED with a break at 912 Angstroms and different
    slopes above and below the break.
    """
    
    # From http://www.astro.wisc.edu/~dolan/constants.html
    planck_erg_s = 6.6260755e-27 #erg s

    def lambdanu(self, wavelength):
        """Convert between wavelength (Ang) and frequency (Hz) or vice-versa.
        """
        if numpy.all(wavelength == 0):
            wavelength = numpy.array(wavelength) 
        return cc.c_light_cm_s / (wavelength * cc.angstrom_cm)


    def luminosity_wavelength(wavelength0, wavelength1):
        """Normalized luminosity in the defined band.

        wavelength0 and wavelength1 in Angstroms.

        This is the fraction of total luminosity in the band.
        Multiply the result by the total luminosity (energy per
        second) to get physical units.
        """
        return self.sed.integrate(self.lambdanu(wavelength0),
                                  self.lambdanu(wavelength1))

    def photonRate_wavelength(self, wavelength0, wavelength1):
        """Normalized photon emission rate in the band between
        wavelength0 and wavelength1.

        Units are erg^-1 (which could also be expressed as s^-1 per
        (erg/s)).

        Multiply the result by the total luminosity (in ergs per unit
        time), to get the actual photon emission rate.

        Example
        -------

        To get the ionizing photon emission rate (per unit luminosity):
        >>> BrokenPowerlawSED().photonRate_wavelength(0., 912.)
        3272819078.0292048
        """
        return ((1./self.planck_erg_s) *
                self.sed.integrate(self.lambdanu(wavelength0),
                                   self.lambdanu(wavelength1),
                                   weight_power=-1.))
    
    def __init__(self, s_ion=-3.0, s_red=0.0, break_factor=6.0):
        """Return a model SED for a galaxy.

        Parameters
        ----------

        s_ion:
            spectral index (f_nu ~ nu^s_ion) at lambda < 912 Ang

        s_red:
            spectral index (f_nu ~ nu^s_red) at lambda > 912 Ang

        Notes
        -----

        Bolton and Haehnelt (2007MNRAS.382..325B) use an SED with

            eps_nu ~ v^0  for 912 < lambda < 3000 Ang.
                   ~ v^-3 for       labmda < 912 Ang.
    
        'with an additional break in the spectrum at the Lyman limit'
        
            eps_L = eps(1500)/6
        
        """

        limits_Ang=numpy.array((3000., 912., 0.))
        limits_nu = self.lambdanu(limits_Ang)

        nu_912 = self.lambdanu(912.)

        # Ratio of flux at 1500 Ang to flux at 912 Ang (redward of the break).
        ratio_1500_912 = (912./1500.)**s_red
        
        # Ratio of the ionizing to red coefficients.
        coeff_ratio = ratio_1500_912 * nu_912**(s_red - s_ion) / break_factor

        coefficients = numpy.array((1., coeff_ratio))

        # Define an SED:
        powers = numpy.array((s_red, s_ion))
        sed = utils.PiecewisePowerlaw(limits_nu, powers, coefficients)
        self.sed = sed
        self.break_factor = break_factor
        self.s_ion = s_ion
        self.s_red = s_red

    def __call__(self, nu):
        """The spectrum at the given frequency/frequencies.

        Multiply by the total luminosity to get the luminosity per
        unit frequency.

        Units are (fraction of total luminosity) per Hz.
        """
        return self.sed(nu)

    def iPhotonRateRatio(self, wavelength=1500.):
        """The ratio of ionizing photon emission rate to luminosity.

        Paramters
        ---------
        
        lambda:
            wavelength in Angstroms

        Returns
        -------

        The ratio of ionizing photon emission rate to luminosity at
        the given wavelength Q/L_nu(lambda) in units of photons s^-1
        (erg s^-1 Hz^-1)^-1.

        Notes
        -----

        While this function takes an argument in angstroms, the ratio
        is calculated using the luminosity per unit frequence, so as
        to be easily comensurate with luminosities inferred from AB
        magnitudes.
        
        """

        return (self.photonRate_wavelength(0.,912.) /
                self(self.lambdanu(wavelength)))

if __name__ == '__main__':

    usage = """Run with a filename argument to produce image files, e.g.:
    python luminosityfunction.py lum.png
    """
    
    ### Argument parsing. ###
    parser = optparse.OptionParser(usage)
    parser.add_option("-f", "--file", action='store_true', dest="filename",
                      default=None)
    (options, args) = parser.parse_args()
    if (len(args) > 0) and (options.filename is None):
        options.filename = args[0]

    if options.filename is None:
        print("No filename given.")
        print(usage)
    else:
        prefix, extension = os.path.splitext(options.filename)
    
    import pylab

    plotLFevo() # Canonical observed LF
    plotLFevo(extrap_var='z')
    plotLFevo(params=B2008_fixed)
    plotLFevo(params=B2008_fixed, extrap_var='z')

    test_plot_schechter()
    
#    plotLFevo(maglim=None) # Integrate LF all the way to zero luminosity
#    plotLFevo(B2008_fixed) # Use fixed alpha, phiStar
#    plotLFevo(B2008_fixed, maglim=None) # ... integrated to L=0

    # ... with alpha raised to -1.4 (from -1.74)
#    newparams = B2008_fixed
#    newparams['alpha'] = numpy.array(newparams['alpha'])
#    newparams['alpha'][:] = -1.4

#    plotLFevo(newparams, maglim=None)

#    test_plot_schechter()

    ### Plot output code. ###
    if options.filename is None:
        pylab.show()
    else:
        from matplotlib import _pylab_helpers
        for manager in _pylab_helpers.Gcf.get_all_fig_managers():
            fig = manager.canvas.figure
            if len(fig.get_label()) > 0:
                label = fig.get_label()
            else:
                label = '_Fig' + str(manager.num)
            newfilename = prefix + '_' + label + extension
            fig.savefig(newfilename, dpi=75)#, bbox_inches="tight")
	
