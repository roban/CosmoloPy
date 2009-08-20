"""Test that analytical and numerical optical depth calculations agree.

Also reproduces a figure from Griffiths et al. astro-ph/9812125v3
(don't use the published version, which contains errors).
"""

import math

import numpy
import numpy.testing.utils as ntest
import matplotlib.pyplot as pylab

import cosmolopy.distance as cd
import cosmolopy.constants as cc
import cosmolopy.density as cden
import cosmolopy.reionization as cr


def test_GBL_tau_star():
    """Test tau_* against GB&L astro-ph/9812125v3.

    tau_* is a quantity used in optical_depth_instant.
    """
    z = 1.0

    # Fully ionized H and He
    x_ionH = 1.0
    x_ionHe = 2.0

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[0.3],[0.6],[1.0]])
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.65
    cosmo['omega_b_0'] = 0.02 / cosmo['h']**2.
    cosmo['Y_He'] = 0.24

    tau_inst, tau_star = cr.optical_depth_instant(z, 
                                                  x_ionH=x_ionH, 
                                                  x_ionHe=x_ionHe, 
                                                  return_tau_star=True,
                                                  **cosmo)
    print "tau_star = %.7f" % tau_star
    print ("tau_star/(h Omega_b) = %.7f =? 0.061" % 
           (tau_star / (cosmo['h'] * cosmo['omega_b_0'])))

    ntest.assert_approx_equal(tau_star / (cosmo['h'] * cosmo['omega_b_0']),
                              0.061,
                              2)

    print "(1 - Y_He/2) = %.3f =? 0.88" % (1. - (cosmo['Y_He']/2.))
    ntest.assert_approx_equal((1. - (cosmo['Y_He']/2.)),
                              0.88,
                              7)

    H_0 = cc.H100_s * cosmo['h']

    # s^-1 * Mpc s^-1 * Mpc^2 / Mpc^3 msun^-1 s^-2 / Msun -> 
    tau_star_explicit =  ((1. - (cosmo['Y_He']/2.)) * 
                          ((3. * H_0 * cosmo['omega_b_0'] * cc.c_light_Mpc_s *
                            cc.sigma_T_Mpc) / 
                           (8. * math.pi * cc.G_const_Mpc_Msun_s * 
                            (cc.m_p_g/cc.M_sun_g))))

    print "tau_star_explicit = %.7f =? tau_star" % tau_star_explicit
    ntest.assert_approx_equal(tau_star, tau_star_explicit, 3)

def test_GBL_tau_inst():
    """Test match between analytical and numerical tau with instant
    reionization.

    Also makes a plot reproducing figure 1 of arXiv:astro-ph/9812125v3.
    """
    dz = 0.05
    z = numpy.arange(0., 80. + 1.5*dz, dz)

    # Fully ionized H and He
    x_ionH = 1.0
    x_ionHe = 2.0

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[0.3],[0.6],[1.0]])
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.65
    cosmo['omega_b_0'] = 0.02 / cosmo['h']**2.
    cosmo['Y_He'] = 0.24

    tau_inst = cr.optical_depth_instant(z, x_ionH=x_ionH, x_ionHe=x_ionHe, 
                                        **cosmo)
    tau_int = cr.integrate_optical_depth(x_ionH, x_ionHe, z, **cosmo)

    linestyle = ['-', ':', '--']
    
    pylab.figure()
    pylab.subplot(2,1,1)
    pylab.title("Compare to GB&L fig. 1 (astro-ph/9812125v3.)")
    for i in range(len(linestyle)):
        pylab.plot(z, tau_inst[i], ls=linestyle[i], color='b')
        pylab.plot(z[1:], tau_int[i], ls=linestyle[i], color='r')

    pylab.xlim(0,80)
    pylab.ylim(0,1)
    pylab.xlabel(r"$\mathrm{z_{ion}}$")
    pylab.ylabel(r"$\tau$")
    
    pylab.subplot(2,1,2)
    for i in range(len(linestyle)):
        pylab.plot(z[1:], 
                   1.e4 * (tau_int[i,:] - tau_inst[i,1:])/tau_inst[i,1:], 
                   ls=linestyle[i], color='k')
        diff = (tau_int[i,:] - tau_inst[i,1:]) / tau_inst[i,1:]

        print ("max fractional error in num. int. = %.3g" % 
               numpy.max(numpy.abs(diff))
               )
        ntest.assert_array_less(numpy.abs(diff), 2.e-4)

    pylab.xlim(0,40)
    pylab.xlabel(r"$\mathrm{z_{ion}}$")
    pylab.ylabel(r"$\mathrm{10^4 \times (num.\tau - ana.\tau)/ana.\tau}$")

if __name__ == "__main__":
    test_GBL_tau_star()
    test_GBL_tau_inst()
    pylab.show()

