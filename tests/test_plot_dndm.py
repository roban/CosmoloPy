"""Reproduce halo mass function plot from Reed et al. (2007 MNRAS 374 2) 
"""
import inspect

import numpy
import matplotlib.pyplot as pylab

import cosmolopy.distance as cd
import cosmolopy.constants as cc
import cosmolopy.perturbation as cp
import cosmolopy.reionization as cr
import cosmolopy.parameters as cparams

def reed_cosmo(flat=True, extras=True):
    """Cosmological parameters used in Reed et al. (2007 MNRAS 374 2)."""
    cosmo = {'omega_b_0' : 0.045,
             'omega_M_0' : 0.25,
             'omega_lambda_0' : 0.75,
             'h' : 0.73,
             'n' : 1.0,
             'sigma_8' : 0.9
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        cparams.add_extras(cosmo)
    return cosmo

def test_figure1():
    """Plot Reed et al. fig. 1: halo mass functions (no quant. tests).
    """

    cosmo = reed_cosmo()

    m = numpy.logspace(5, 12, 200) / cosmo['h']
        
    z = numpy.array([[30., 20., 10.]]).transpose()
             
    # P-S halo mass function:
    dndmPS = cp.dndm_PS(m, z, **cosmo)
    dndlogmPS_plot = dndmPS * m / cosmo['h']**3. / numpy.log(10.)

    # S-T halo mass function:
    dndmST = cp.dndm_ST(m, z, **cosmo)
    dndlogmST_plot = dndmST * m / cosmo['h']**3. / numpy.log(10.) 

    pylab.figure(figsize=(6,6))

    for iz in range(len(z)):
        print m.shape
        print dndlogmPS_plot.shape
        pylab.plot(m * cosmo['h'], dndlogmST_plot[iz], 'k--',
                   label='Sheth-Tormen')
        pylab.plot(m * cosmo['h'], dndlogmPS_plot[iz], 'r:',
                   label='Press-Schechter')

    pylab.xlabel(r"$m$ ($h^{-1} M_\mathrm{sun}$)")
    pylab.ylabel(r"$dn/d\log_{10}(m)$ ($h^{3} Mpc^{-3}$)")

    pylab.yscale('log')
    pylab.xscale('log')
    pylab.ylim(1e-7, 3e4)
    pylab.xlim(1e5, 1e12)

    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " Reed et al. (2007 MNRAS 374 2).")

if __name__ == '__main__':
    test_figure1()
    pylab.show()


