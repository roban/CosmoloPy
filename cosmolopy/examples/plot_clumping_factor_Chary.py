import numpy

import cosmolopy.reionization as cr

def plot_clumping_fractor_Chary():

    z = numpy.linspace(5,15,500)
    cf = cr.clumping_factor_Chary(z)

    import pylab
    pylab.plot(z,cf)
    pylab.yscale('log')
    pylab.xlim(5,15)
    pylab.ylim(1,100)
    pylab.show()
