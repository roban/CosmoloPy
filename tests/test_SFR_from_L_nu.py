"""Test SFR - UV luminosity conversion.
"""

import numpy

import cosmolopy.parameters as cparam

import cosmolopy.magnitudes as magnitudes
import cosmolopy.luminosityfunction as luminosityfunction


def test_sfr():
    """Check that M = -18 -> SFR ~ 1 Msun/yr.

    From Oesch et al (2009 ApJ 690 1350).
    """

    mag = -18
    lum = magnitudes.L_nu_from_magAB(mag)
    sfr = luminosityfunction.sfr_from_L_nu(lum)

    print "M = %.3g -> L_nu = %.3g erg/s/Hz -> SFR = %.3g Msun/yr" % (mag,
                                                                      lum,
                                                                      sfr)
    assert numpy.round(sfr) == 1

def test_sfrd():
    """ Check conversion of Lmin -> SFRD.

    At z=7:
    Lmin = 0.2 L*z=6 -> Mmin=-18.5 -> log(rho/erg/s/Hz/Mpc^3)=25.5 ->
    SFRD = 10^-2.32 Msun/yr^-1/Mpc^-3.

    From Oesch et al (2009 ApJ 690 1350).
    """

    cosmo = cparam.WMAP7_BAO_H0_mean(flat=True)


#     # From Bouwens 2008ApJ...686..230B
#     lfhist = luminosityfunction.LFHistory(params=luminosityfunction.B2008,
#                                           **cosmo)
    
#     mStarz6 = lfhist.params_z(6.0)['MStar']

    alpha = -1.74
    mStarz6 = -20.24
    mStarz7 = -19.7
    phiStar = 1.4e-3
    
    lStarz6 = magnitudes.L_nu_from_magAB(mStarz6)
    lmin = 0.2 * lStarz6
    magmin = magnitudes.magnitude_AB_from_L_nu(lmin)
    ltot = luminosityfunction.schechterCumuLM(magnitudeAB=magmin,
                                              MStar=mStarz7,
                                              phiStar=phiStar,
                                              alpha=alpha)
    sfrd = luminosityfunction.sfr_from_L_nu(ltot)

    print """Lmin = 0.2 L*z=6 -> Lmin/erg/s/Hz = %.3g -> Mmin = %.3g ->
    log(rho/(erg/s/Hz/Mpc^3)) = %.3g -> log(SFRD/(MSun/yr)) = %.3g"""\
    % (lmin, magmin, numpy.log10(ltot), numpy.log10(sfrd))

    ltotz6 = luminosityfunction.schechterCumuLM(magnitudeAB=magmin,
                                                MStar=mStarz6,
                                                phiStar=phiStar,
                                                alpha=alpha)

    print "luminosity density increase from z=7 to z=6 is %.2g percent." \
          % (100 * (1. - ltot/ltotz6))


    assert numpy.abs(numpy.round(1. - ltot/ltotz6, 1) - 0.5) < 0.1
    assert numpy.round(magmin, 1) == -18.5
    assert numpy.abs(numpy.round(numpy.log10(ltot), 1) - 25.5) < 0.2
    assert numpy.abs(numpy.round(numpy.log10(sfrd), 1) - -2.32) < 0.05

if __name__ == '__main__':
    test_sfr()
    test_sfrd()

    
