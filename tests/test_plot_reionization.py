"""Plot IGM ionized fraction.

Note that there are no actual quantitative tests here. Just checks
that it runs without errors.

"""

import numpy
import numpy.testing.utils as ntest
import matplotlib.pyplot as pylab

import cosmolopy.perturbation as cp
import cosmolopy.distance as cd
import cosmolopy.density as cden
import cosmolopy.reionization as cr
import cosmolopy.constants as cc
import cosmolopy.parameters as cparam

import testingutils as tu

def test_plot_FZH():
    """Plot figure 3 from FZH (2004ApJ...613....1F) (no quantitative tests).
    """

    cosmo = {}
    cosmo['omega_M_0'] = 0.3
    cosmo['omega_k_0'] = 0.0
    cosmo['omega_lambda_0'] = 0.7
    cosmo['h'] = 0.7
    cosmo['n'] = 1.0
    cosmo['sigma_8'] = 0.9

    cosmo['omega_b_0'] = 0.0462 # not specified in paper
    cosmo['omega_n_0'] = 0.0  # not specified in paper
    cosmo['N_nu'] = 0  # not specified in paper
    cosmo['Y_He'] = 0.24  # not specified in paper

    dz = 0.1
    z = numpy.arange(25., 8. - 1.5 * dz, -1. * dz)
    
    T_min = 1e4 #K
    c_ion = numpy.array([[500.], [40.], [12.]])

    # calculate ionized fraction from collapse fraction
    x_fcol = cr.ionization_from_collapse(z, 
                                         c_ion, 
                                         T_min,
                                         **cosmo)
    #linestyle = ['-', ':', '--']
    color = ['r', 'g', 'b']
    pylab.figure()
    pylab.subplot(2,1,1)
    for i in range(len(color)):
        pylab.plot(z, x_fcol[i], ls='--', color=color[i])
    pylab.axhline(y=0.75)
    pylab.yscale('log')
    pylab.xlim(8,25)
    pylab.ylim(1e-4, 1)
    pylab.title("Compare to figure 3 from FZH (2004ApJ...613....1F)")

    pylab.subplot(2,1,2)
    for i in range(len(color)):
        pylab.plot(z, x_fcol[i], ls='--', color=color[i])
    pylab.axhline(y=0.75)
    pylab.xlim(8,25)
    pylab.ylim(0, 1)

def test_plot_integrate_ion_recomb_collapse():
    """Plot results of integrate_ion_recomb_collapse (no quantitative tests).
    """

    cosmo = cparam.WMAP5_mean(flat=True)

    dz = 0.5
    z = numpy.arange(25., 8. - 1.5 * dz, -1. * dz)
    
    T_min = 1e4 #K
    c_ion = numpy.array([[500.], [40.], [12.]])

    # calculate ionized fraction from collapse fraction
    x_fcol = cr.ionization_from_collapse(z, 
                                         c_ion, 
                                         T_min,
                                         **cosmo)

    x_rec = numpy.empty(x_fcol.shape)
    w_rec = numpy.empty(x_fcol.shape)
    for i in range(x_fcol.shape[0]):
        # calculate ionized fraction, including recombinations
        x_rec[i], w_rec[i], t = cr.integrate_ion_recomb_collapse(z, c_ion[i,0],
                                                               temp_min = T_min,
                                                               **cosmo)
    
    #linestyle = ['-', ':', '--']
    color = ['r', 'g', 'b']
    pylab.figure()
    pylab.subplot(2,1,1)
    for i in range(len(color)):
        pylab.plot(z, x_fcol[i], ls='--', color=color[i])
        pylab.plot(z, x_rec[i], ls='-', color=color[i])
        pylab.plot(z, w_rec[i], ls=':', color=color[i])
    pylab.axhline(y=0.75)
    pylab.yscale('log')
    pylab.xlim(8,25)
    pylab.ylim(1e-4, 1)

    pylab.subplot(2,1,2)
    for i in range(len(color)):
        pylab.plot(z, x_fcol[i], ls='--', color=color[i])
        pylab.plot(z, x_rec[i], ls='-', color=color[i])
        pylab.plot(z, w_rec[i], ls=':', color=color[i])
    pylab.axhline(y=0.75)
    pylab.xlim(8,25)
    pylab.ylim(0, 1)

if __name__ == "__main__":
    test_plot_FZH()
    test_plot_integrate_ion_recomb_collapse()
    pylab.show()

