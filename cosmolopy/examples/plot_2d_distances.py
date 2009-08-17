"""Plot some distance measures versus redshift and omega_M.

"""
import sys
import getopt

import numpy
import matplotlib.pyplot as pylab
import matplotlib.cm as cm

import cosmolopy.distance as cd
import cosmolopy.constants as cc

def plot_DM(filename):
    """The dimensionless proper motion distance DM/DH. 
    """

    dz = 0.1
    z = numpy.arange(0., 5. + 1.1 * dz, dz)

    cosmo = {}
    dom = 0.1
    om = numpy.atleast_2d(numpy.arange(0.1,
                                       1. + 1.1 * dom,
                                       dom)).transpose()
    cosmo['omega_M_0'] = om
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.701
    
    dh = cd.hubble_distance_z(0, **cosmo)
    dm, dm_err = cd.comoving_distance_transverse(z, **cosmo)

    plot_dist(z, dz, om, dom, dm, dh, 'proper motion distance', r'D_M', 
              filename)

def plot_DA(filename):
    """The dimensionless angular diameter distance DA/DH. 
    """

    dz = 0.1
    z = numpy.arange(0., 5. + dz, dz)

    cosmo = {}
    dom = 0.05
    om = numpy.atleast_2d(numpy.arange(0.1,
                                       1. + dom,
                                       dom)).transpose()
    cosmo['omega_M_0'] = om
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.701
    
    dh = cd.hubble_distance_z(0, **cosmo)
    da, da_err1, da_err2 = cd.angular_diameter_distance(z, **cosmo)

    plot_dist(z, dz, om, dom, da, dh, 'angular diameter distance', r'D_A',
              filename)

def make_and_plot_DA():
    """The dimensionless angular diameter distance DA/DH. 
    """

    dz = 0.1
    z = numpy.arange(0., 5. + dz, dz)

    cosmo = {}
    dom = 0.05
    om = numpy.atleast_2d(numpy.arange(0.1,
                                       1. + dom,
                                       dom)).transpose()
    cosmo['omega_M_0'] = om
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.701
    
    dh = cd.hubble_distance_z(0, **cosmo)
    da, da_err1, da_err2 = cd.angular_diameter_distance(z, **cosmo)

    plot_dist(z, dz, om, dom, da, dh, 'angular diameter distance', r'D_A')

def plot_dist(z, dz, om, dom, dist, dh, name, mathname, filename=None):
    x, y = numpy.meshgrid(z, om)
    pylab.figure(figsize=(5.5,4.5))    
    pylab.imshow(dist/dh, 
                 extent=(z.min() - dz/2., 
                         z.max() + dz/2.,
                         om.max() + dom/2.,
                         om.min() - dom/2., 
                         ),
                 interpolation='nearest',
                 aspect = z.max()/om.max(),
                 cmap = cm.Spectral,
                 )
    cb = pylab.colorbar()
    cb.ax.set_ylabel(r'$' + mathname + '/D_H$')

    pylab.contour(x, y, dist/dh, 10, colors='k')
    pylab.xlim(z.min(), z.max())
    pylab.ylim(om.min(), om.max()) 
    pylab.xlabel("redshift z")
    pylab.ylabel(r"$\Omega_M = 1 - \Omega_\lambda$")
    pylab.title(name)
    if filename is not None:
        prefix, extension = filename.split('.')
        pylab.savefig(prefix + '_' + mathname + '.' + extension,
                      bbox_inches="tight")

if __name__ == "__main__":
    if len(sys.argv)==1:
        print "Run with a filename argument to produce image files, e.g.:"
        print " python plot_2d_distances.py dist2d.png"
        print " python plot_2d_distances.py dist2d.eps"
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = None
        
    plot_DM(filename)
    plot_DA(filename)

    if filename is None:
        pylab.show()

