"""Reproduce collapse fraction plot from Haiman and Holder (2003)

Haiman and Holder 2003ApJ...595....1H.

"""
import inspect

import numpy
import matplotlib.pyplot as pylab

import cosmolopy.distance as cd
import cosmolopy.constants as cc
import cosmolopy.perturbation as cp
import cosmolopy.reionization as cr
import cosmolopy.parameters as cparams

def hh_cosmo(flat=True, extras=True):
    """Cosmological parameters used in HH2003."""
    cosmo = {'omega_b_0' : 0.047,
             'omega_M_0' : 0.29,
             'omega_lambda_0' : 0.71,
             'h' : 0.72,
             'n' : 0.99,
             'sigma_8' : 0.9
             }
    if flat:
        cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
        cosmo['omega_k_0'] = 0.0
    if extras:
        cparams.add_extras(cosmo)
    return cosmo

def test_figure1():
    """Plot HH fig. 1: Evolution of the collapsed fraction (no quant. tests).

    "Evolution of the collapsed fraction of all baryons in halos in
    three different virial temperature ranges."

    Haiman and Holder 2003ApJ...595....1H.

    """

    dz = 1.
    z = numpy.arange(45., 0.0, -1. * dz)

    cosmo = hh_cosmo()
    
    T_min = numpy.array([[1.0e2],
                         [1.0e4],
                         [2.0e5]])
             
    linestyle = ['-', ':', '--']

    # Collapse fraction with the various minimum Virial temps.
    fc = cp.collapse_fraction(*cp.sig_del(T_min, z, **cosmo))

    # Calculate fraction of mass in halos in the ranges of min T.
    fc_diff = numpy.empty((T_min.shape[0], z.shape[0]))
    fc_diff[0:-1] = fc[0:-1] - fc[1:]
    fc_diff[-1] = fc[-1]

    pylab.figure(figsize=(6,6))    
    for i in range(len(linestyle)):
        pylab.plot(z, fc_diff[i], ls=linestyle[i])

    pylab.xlabel(r"$\mathrm{z}$")
    pylab.ylabel(r"$\mathrm{F_{coll}}$")

    pylab.yscale('log')
    pylab.ylim(1e-5, 1.0)
    pylab.xlim(45.,0.)

    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " H&H (2003ApJ...595....1H).")

if __name__ == '__main__':
    test_figure1()
    pylab.show()


