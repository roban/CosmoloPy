import numpy

from cosmolopy.utils import *
from cosmolopy.saveable import Saveable


######### Testing ############

def test_integrate_piecewise(pieces=2, method='quad'):
    # Test modified from scipy test_quadrature.py.
    n = 2
    z = 1.8
    def func(x):       # Bessel function integrand
        return numpy.cos(n*x-z*numpy.sin(x))/numpy.pi
    x = numpy.linspace(0, numpy.pi, pieces)
    val = integrate_piecewise(func,x, method=method)
    table_val = 0.30614353532540296487
    diff = val[-1] - table_val
    print "Error with %i %s pieces = %.3g" % (pieces, method, diff)
    numpy.testing.assert_almost_equal(val[-1], table_val, decimal=7)

def test_PiecewisePowerlaw(n=4, plot=False):

    powers = 1.5 * (0.5 - numpy.random.rand(n))
    limits = 10. * numpy.cumsum(numpy.random.rand(n+1))

    pfunc = PiecewisePowerlaw(limits, powers, externalval=0)
    x = numpy.linspace(limits[0] - 0.1, limits[-1] + 0.1, 20.)
    #x = numpy.linspace(limits[0], limits[-1], 20.)
    y = pfunc(x)

    integral = pfunc.integrate(0, x)
    numintegral = integrate_piecewise(pfunc, x, method='quad')

    integral2 = pfunc.integrate(x, x[-1])
    numintegral2 = numintegral[-1] - numintegral

    # Weighted integral
    integral3 = pfunc.integrate(0, x, weight_power=1.5)
    weightedfunc = lambda x: (x**1.5) * pfunc(x)
    numintegral3 = integrate_piecewise(weightedfunc, x, method='quad')

    if plot:
        import pylab

        pylab.subplot(221)
        pylab.title('x vs. y')
        pylab.plot(x,y)
        pylab.xlim(min(x), max(x))
        for xlim in pfunc._limits.flat:
            pylab.axvline(x=xlim)
        
        pylab.subplot(223)
        pylab.title('x vs. forward/reverse integral')
        pylab.plot(x, integral)
        pylab.plot(x, numintegral)
        pylab.plot(x, integral2)
        pylab.plot(x, numintegral2)
        pylab.plot(pfunc._limits.flat[1:], numpy.cumsum(pfunc._integrals), '.')
        pylab.xlim(min(x), max(x))

        pylab.subplot(222)
        pylab.title('frac. diffs (num - ana)/ana')
        pylab.plot(x, (integral - numintegral)/integral)
        pylab.plot(x, (integral2 - numintegral2)/integral2)
        pylab.plot(x, (integral3 - numintegral3)/integral3)

        pylab.subplot(224)
        pylab.title('weighted integrals')
        pylab.plot(x, integral3)
        pylab.plot(x, numintegral3)
        
        pylab.show()
    assert numpy.all(numpy.abs(integral - numintegral) < 1e-4)
    assert numpy.all(numpy.abs(integral2 - numintegral2) < 1e-4)
    assert numpy.all(numpy.abs(integral3[integral3>0] -
                               numintegral3[integral3>0]) /
                     integral3[integral3>0] < 1e-3)
    assert numpy.abs(integral[-1] - 1.) < 1e-4
    assert numpy.abs(integral2[0] - 1.) < 1e-4


def test_Extrapolate1d():

    slope = 2. * (0.5 - numpy.random.rand(1))
    intercept = 20. * (0.5 - numpy.random.rand(1))
    x = numpy.array([3., 4.,
                     5., 6., 7.,
                     8., 9.])
    y = slope * x + intercept

    # Random deviations in the midde.
    y[2:5] = y[2:5] + (5. * (0.5 - numpy.random.rand(3)))
    x1 = numpy.linspace(0., 15., 100.)

    extrap = Extrapolate1d(x,y)
    print extrap.extrap_string()
    y1 = extrap(x1)
    ytrue = slope * x1 + intercept

    # Test extrapolation with a fixed slope.
    newslopes = [3.0, 2.0]
    extrap2 = Extrapolate1d(x, y, slopes=newslopes)
    print extrap2.extrap_string()
    y2 = extrap2(x1)

    mask = numpy.logical_or(x1 >= x[5],
                            x1 <= x[1])

    assert numpy.all(numpy.abs((y1 - ytrue)[mask]) < 1e-10)

    import pylab
    pylab.plot(x, y, 'o')
    pylab.plot(x1, y1, '-')
    pylab.plot(x1, ytrue, ':')
    pylab.plot(x1, y2, '--')

if __name__ == '__main__':

    import pylab
    pylab.figure()
    for i in range(4):
        test_Extrapolate1d()
    pylab.show()


    test_PiecewisePowerlaw(plot=True)
    
    import numpy.testing
    for method in ['quad', 'romberg']:
        for pieces in xrange(2,1000,300):
            test_integrate_piecewise(pieces, method)
