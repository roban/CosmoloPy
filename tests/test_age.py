"""Test integrated age against analytical age."""

import numpy
import numpy.testing.utils as ntest
import matplotlib.pyplot as pylab

import cosmolopy.distance as cd
import cosmolopy.constants as cc


def test_age():
    """Test integrated age against analytical age."""
    z = numpy.arange(0, 10.0, 0.05)

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[0.99],[0.01],[0.3]])
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.7
    cd.set_omega_k_0(cosmo)

    linestyle = ['-', ':', '--']

    gyr = 1e9 * cc.yr_s

    tl = cd.lookback_time(z, **cosmo)
    age = cd.age(z, **cosmo)
    age_ana = cd.age_flat(z, **cosmo)

    pylab.figure(figsize=(6,6))
    for i in range(len(linestyle)):
        pylab.plot(z, (tl/gyr)[i], ls=linestyle[i], color='0.5')
        pylab.plot(z, (age/gyr)[i], ls=linestyle[i], color='r')
        pylab.plot(z, (age_ana/gyr)[i], ls=linestyle[i], color='k')
    pylab.xlabel("redshift z")
    pylab.ylabel(r"age $t_L/$Gyr")

    pylab.figure(figsize=(6,6))
    for i in range(len(linestyle)):
        pylab.plot(z, ((age - age_ana)/age_ana)[i], ls=linestyle[i], 
                   color='k')
        # Make sure errors are small:
        ntest.assert_array_less((numpy.abs((age - age_ana)/age_ana)[i]),
                                3e-13)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"age: (integral - analytical)/analytical")

if __name__ == "__main__":
    test_age()
    pylab.show()
