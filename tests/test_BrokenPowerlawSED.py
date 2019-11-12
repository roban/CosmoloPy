
from __future__ import absolute_import, division, print_function

import numpy

import cosmolopy.utils as utils
from cosmolopy.luminosityfunction import BrokenPowerlawSED

def test_BrokenPowerlawSED(plot=False):
    """Test that the BrokenPowerlawSED has the right break and normalization.
    """
    if plot: import pylab
    runtest_BrokenPowerlawSED(BrokenPowerlawSED(), plot=plot)
    runtest_BrokenPowerlawSED(BrokenPowerlawSED(s_ion=-2., s_red=-0.5), plot=plot)
    if plot: pylab.show()

def runtest_BrokenPowerlawSED(sed, plot=False):
    
    nu_912 = sed.lambdanu(912.)
    
    ratio_1500_L_num = (sed(sed.lambdanu(1500.))/
                        sed(sed.lambdanu(912. - 1e-10)))[0]

    # Test the magnitude of the break.
    assert numpy.abs(ratio_1500_L_num - sed.break_factor) < 1e-9

    # Test the normalization of the integral.
    numintegral = utils.logquad(sed, sed.lambdanu(3000.), sed.lambdanu(1e-6))[0]
    assert numpy.abs(numintegral - 1.) < 1e-6

    if plot:
        import pylab
        wav = numpy.arange(100.,2000.,1.0)
        nu = sed.lambdanu(wav)
        pylab.subplot(121)
        #norm = sed(nu_912 + 1e10)
        norm = 1.0
        pylab.plot(wav, sed(nu)/norm)
        #pylab.axhline(y=6, ls=':')
        #pylab.axvline(x=1500, ls=':')
        pylab.subplot(122)
        pylab.plot(nu, sed(nu)/norm)

if __name__ == '__main__':
    test_BrokenPowerlawSED(True)
