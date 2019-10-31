"""Reionization consistency check with BK&P 2009MNRAS.397..971B"""

from __future__ import absolute_import, division, print_function

import numpy
import numpy.testing.utils as ntest
import matplotlib.pyplot as pylab

import cosmolopy.perturbation as cp
import cosmolopy.distance as cd
import cosmolopy.density as cden
import cosmolopy.reionization as cr
import cosmolopy.constants as cc
import cosmolopy.parameters as cparam

def test_tau_BKP():
    """Reionization consistency check with BK&P 2009MNRAS.397..971B

    Bagla J.~S., Kulkarni G., Padmanabhan T., 2009 (2009MNRAS.397..971B).
    
    Take their cannonical set of cosmological parameters and ionization
    coefficients and make sure I get the same optical depth value.

    The first figure demonstrates the calculation of the ionization
    fraction and the optical depth (not shown in the paper).

    The second figure shows dependence of the f_* f_esc,gamma values on
    tau (actually calculated the other way around, f_* f_esc,gamma ->
    tau). The plotted point (x) marks the value from the paper. The
    agreement between the calculated optical depth and the WMAP tau value
    is printed in the text output.
    
    """
    N_gamma = 6804.

    # f_* f_esc N_gamma
    f_ion = numpy.transpose(numpy.atleast_2d(numpy.arange(10., 200., 10.)))
    
    alpha_B = 1e-13 # cm^3 s^-1
    m_min = 1e8 # M_sun

    x_ionHe = 2.0

    dz = 0.1
    z = numpy.arange(20., 6. - 1.5 * dz, -1. * dz)
    #z = numpy.arange(6., 20. + 1.5 * dz, dz)

    cosmos = [cparam.WMAP5_ML(flat=True),
              cparam.WMAP5_mean(flat=True)]    

    pylab.figure(figsize=(8,8))
    colors=['r', 'k']
    names = ["WMAP5 ML", "WMAP5 mean"]
    i = -1
    for cosmo in cosmos:
        i += 1
        print
        print names[i]
        # calculate ionized fraction, including recombinations
        x_rec, w_rec, t = cr.integrate_ion_recomb_collapse(z, f_ion,
                                                           m_min,
                                                           passed_min_mass=True,
                                                           alpha_B=alpha_B,
                                                           **cosmo)
        # calculate the optical depth in this scenario
        tau_z0 = cr.optical_depth_instant(z[-1],
                                          x_ionH=1.0,
                                          x_ionHe=2.0,
                                          **cosmo)

        tau_later = cr.integrate_optical_depth(x_rec, x_ionHe * x_rec, 
                                               z, **cosmo)
        tau_0 = tau_later[:, -1] + tau_z0

        tau = tau_later

        print "tau(WMAP) = %.3f" % cosmo['tau']
        for j in range(len(f_ion.flat)):
            if round(f_ion[j],1) != 50.0:
                continue
            print "with f_* f_esc_gamma N_gamma = %.1f:" % f_ion[j]
            pylab.plot(z, x_rec[j], ls='-', color=colors[i])
            pylab.plot(z, w_rec[j], ls=':', color=colors[i])
            #pylab.plot(z, 10. * tau[j], ls='--', color=colors[i])
            pylab.plot(z, 10. * (tau_0[j] - tau[j]), ls='--', color=colors[i])
            pylab.axhline(y=10. * tau_0[j], ls='--', color=colors[i])
            print "tau(z=0)  = %.4f" % tau_0[j]
            print " fractional diff. = %.3g" % ((tau_0[j] - cosmo['tau']) / 
                                                cosmo['tau'])
            if i==1:
                # Make sure we recover the WMAP value.
                ntest.assert_approx_equal(tau_0[j], cosmo['tau'], 2)

    pylab.ylim(0,1.01)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"ionized fraction or optical depth $\tau \times 10$")

    pylab.figure(figsize=(8,8))
    pylab.plot(tau_0, f_ion / N_gamma, '-', color='k')
    pylab.plot([cosmo['tau']], [50.0 / N_gamma], 'x', color='b')
    pylab.xlim(0.06, 0.12)
    pylab.ylim(0., 0.022)
    pylab.xlabel(r'$\tau$')
    pylab.ylabel(r'$f_* f_{esc,\gamma}$')
        
if __name__ == "__main__":
    test_tau_BKP()
    pylab.show()

