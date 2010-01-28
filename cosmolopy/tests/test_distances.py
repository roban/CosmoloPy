import os
import numpy
import numpy.testing.utils as ntest
import cosmolopy.distance as cd
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

def test_distances(threshold = 1e-3):
    """Compare distance measures with calculations from http://icosmo.org/"""
    cosmo = cosmo_wmap_5()

    print "Comparing distances with calculations from http://icosmo.org/"

    # load external distance calculations
    # z  DA(z)   DT(z)   DL(z) 
    distance_file = os.path.dirname(os.path.abspath(__file__))
    distance_file = os.path.join(distance_file, 'icosmo_testdata', 
                                 'distances.txt')
    ic_dists = numpy.loadtxt(distance_file)

    z = ic_dists[:,0]

    cd_da, err1, err2 = cd.angular_diameter_distance(z, **cosmo)
    cd_dt, err3 = cd.light_travel_distance(z, **cosmo)
    cd_dl, err4, err5 = cd.luminosity_distance(z, **cosmo)
    cd_dm, err6 = cd.comoving_distance_transverse(z, **cosmo)
    
    cd_dists = numpy.vstack((cd_dm, cd_da, cd_dt, cd_dl))

    labels = [ r'$D_M$', r'$D_A$', r'$D_T$', r'$D_L$']
    threshold = [1e-7,    1e-7,     1e-3,    1e-7    ]  

    # print ic_dists[-1].transpose()
    # print cd_dists[:,-1]
    pylab.figure()
    for i in range(len(labels)):
        pylab.plot(ic_dists[:,0], ic_dists[:,i+1], 
                   label=labels[i] +' IC', ls=':')
        pylab.plot(z, cd_dists[i], label=labels[i] +' distance.py', ls='-')
        pylab.legend(loc='best')

    pylab.figure()

    for i in range(len(labels)):
        diff = (cd_dists[i] - ic_dists[:,i+1]) / ic_dists[:,i+1]
        maxdiff = numpy.max(numpy.abs(diff[ic_dists[:,i+1] > 0]))
        print "Maximum fraction difference in %s is %e." % (labels[i],
                                                            maxdiff)
        assert(numpy.all(maxdiff < threshold[i]))
        pylab.plot(ic_dists[:,0], 
                   diff, 
                   label=labels[i], ls='-')

    pylab.plot(z, err2, label=labels[1] + ' err.', ls=':')
    pylab.plot(z, err3, label=labels[2] + ' err.', ls=':')
    pylab.plot(z, err5, label=labels[3] + ' err.', ls=':')
    pylab.plot(z, err6, label=labels[0] + ' err.', ls=':')

    pylab.legend(loc='best')

def test_hubble(threshold = 1e-7):
    cosmo = cosmo_wmap_5()
    print "Comparing hubble constant with calculations from http://icosmo.org/"

    # load external distance calculations
    # z H(z)
    hz_file = os.path.dirname(os.path.abspath(__file__))
    hz_file = os.path.join(hz_file, 'icosmo_testdata', 
                           'h_z.txt')
    ic_hz = numpy.loadtxt(hz_file)

    z = ic_hz[:,0]

    cd_hz = cd.hubble_z(z, **cosmo) * cc.Mpc_km

    label = r"$H(z)$"
    pylab.figure()
    pylab.plot(z, ic_hz[:,1], 
               label=label + ' IC', ls=':')
    pylab.plot(z, cd_hz, label=label + ' distance.py', ls='-')
    pylab.legend(loc='best')

    pylab.figure()
    diff = (ic_hz[:,1] - cd_hz) / ic_hz[:,1]
    
    maxdiff = numpy.max(numpy.abs(diff))
    print "Maximum fraction difference in %s is %e." % (label,
                                                        maxdiff)
    if maxdiff > threshold:
        print "Warning: difference exceeds threshold %e !!!" % threshold
    assert(maxdiff < threshold)

    pylab.plot(z,
               diff, 
               label=label, ls='-')
    pylab.legend(loc='best')


if __name__ == "__main__":
    test_distances()
    test_hubble()
    pylab.show()



