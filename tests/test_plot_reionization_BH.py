"""Reproduce some results from Bolton and Haehnelt (2007MNRAS.382..325B).
"""

import pylab

import numpy
import scipy.interpolate

import cosmolopy.distance as cd
import cosmolopy.reionization as cr
import cosmolopy.perturbation as cp

def nIonDotLowz(z):
    """Ionizations per second per Mpc^3 from BH2007 for z < 6
    """
    return 10.**(50.5 - 0.06 * (z - 6.0))

def nIonDotBHmodel2(z):
    """Ionization model 2 from BH2007: constant above z=6.
    """
    return ((z < 6) * nIonDotLowz(z) +
            (z >= 6) * nIonDotLowz(6))

def test_plot_nIonDotBH():
    """Compare with BH2007 figure 7."""
    z = numpy.linspace(1.75, 15., 100)
    n2 = nIonDotBHmodel2(z)

    pylab.figure()
    pylab.plot(z, numpy.log10(n2), label='Model 2')
    pylab.xlim(1.75, 15.)
    pylab.ylim(48.9, 51.8)
    pylab.legend(loc='best')

def test_plot_integrate_ionization_recomb_BH(xHe=0.0):

    cosmo = {'omega_b_0' : 0.0463,
             'omega_M_0' : 0.26,
             'omega_lambda_0' : 1.-0.26,
             'omega_k_0' : 0.0,
             'h' : 0.72,
             'n' : 0.95,
             'sigma_8' : 0.85,
             #'tau' : 0.09,
             #'tau_err' : 0.03,
             'omega_n_0' : 0.0,
             'N_nu': 0,
             'Y_He': 0.24,
             'baryonic_effects' : False,
             }

    nIonDot2 = lambda z1: cr.ionization_from_luminosity(z1,
                                                        nIonDotBHmodel2,
                                                        xHe=xHe,
                                                        **cosmo)
    
    zi = numpy.linspace(3.,15.,500.)
    w2i = nIonDot2(zi)
    ifunc2 = scipy.interpolate.interp1d(zi, w2i)
    ifunc2b = scipy.interpolate.interp1d(zi, w2i * 2.0)

    clump_fact_func = lambda z1: 2.0

    z = numpy.linspace(15,4,100.)
    fIon2, w2, t2 = cr.integrate_ion_recomb(z,
                                            ifunc2,
                                            clump_fact_func,
                                            xHe=xHe,
                                            **cosmo)
    fIon2b, w2b, t2b = cr.integrate_ion_recomb(z,
                                               ifunc2b,
                                               clump_fact_func,
                                               xHe=xHe,
                                               **cosmo)

    pylab.plot(z, fIon2, "--", label='Model 2')
    pylab.plot(z, fIon2b, "-", label='Model 2b')
    pylab.plot(z, w2, "0.75", label='M2 no rec.')

    pylab.xlim(4,15)
    pylab.ylim(0,1.0)
    pylab.legend(loc='best')

if __name__ == '__main__':
    test_plot_nIonDotBH()
    pylab.figure()
    test_plot_integrate_ionization_recomb_BH()
    pylab.show()
    
