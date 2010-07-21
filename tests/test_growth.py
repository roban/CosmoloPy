import os
import numpy
import numpy.testing.utils as ntest
import cosmolopy.perturbation as cp
import cosmolopy.constants as cc

import matplotlib.pyplot as pylab

def cosmo_wmap_5():
    """Cosmology from Komatsu et al. 2008 (arXiv:0803.0547v1)"""
    cosmo = {}
    omega_c = 0.233 #+-0.013, 
    cosmo['omega_b_0'] = 0.0462 #+-0.0015, 
    cosmo['omega_M_0'] = omega_c + cosmo['omega_b_0'] # 0.2792
    cosmo['omega_k_0'] = 0.0
    cosmo['omega_lambda_0'] = (1. - cosmo['omega_M_0'] - cosmo['omega_k_0'])
    #cosmo['omega_lambda_0'] = 0.721 #+-0.015
    cosmo['h'] = 0.701 #+-0.0013 
    cosmo['n'] = 0.960 #+0.014-0.013
    cosmo['sigma_8'] = 0.817 #+-0.026
    y_p = 0.240 # +- 0.006: the observed helium abundance
    cosmo['X_H'] = 1. / (1. + 1. / (4. * ((1. / y_p) - 1.)))

    cosmo['omega_n_0'] = 0.0
    cosmo['N_nu'] = 0

    return cosmo

def test_growth(cosmo=None):
    if cosmo is None:
        cosmo = cosmo_wmap_5()
    print "Comparing growth factor with calculations from http://icosmo.org/"

    # load external distance calculations
    # z D
    growth_file = os.path.dirname(os.path.abspath(__file__))
    growth_file = os.path.join(growth_file, 'icosmo_testdata', 'growth.txt')
    ic_growth = numpy.loadtxt(growth_file)

    z = ic_growth[:,0]

    cp_growth = cp.fgrowth(z, cosmo['omega_M_0'])

    label = r"$D(z)$"
    pylab.figure()
    pylab.plot(z, ic_growth[:,1], 
               label=label + ' IC', ls=':')
    pylab.plot(z, cp_growth, label=label + ' perturbation.py', ls='-')
    pylab.legend(loc='best')

    pylab.figure()
    diff = (ic_growth[:,1] - cp_growth) / ic_growth[:,1]

    maxdiff = numpy.max(numpy.abs(diff))
    print "Maximum fraction difference in %s is %e." % (label,
                                                        maxdiff)
    pylab.plot(z,
               diff, 
               label=label, ls='-')
    pylab.legend(loc='best')

    ntest.assert_array_less(numpy.abs(diff), 5e-3)

if __name__ == "__main__":
    test_growth()
    pylab.show()



