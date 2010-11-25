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

    # Set up an array of redshift values.
    dz = 0.1
    z = numpy.arange(0., 10. + 1.1 * dz, dz)

    # Set up a cosmology dictionary, with an array of matter density values.
    cosmo = {}
    dom = 0.01
    om = numpy.atleast_2d(numpy.linspace(0.1, 1.0, (1.-0.1)/dom)).transpose()
    cosmo['omega_M_0'] = om
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.701
    cosmo['omega_k_0'] = 0.0

    # Calculate the hubble distance.
    dh = cd.hubble_distance_z(0, **cosmo)
    # Calculate the comoving distance.
    dm = cd.comoving_distance_transverse(z, **cosmo)

    # Make plots.
    plot_dist(z, dz, om, dom, dm, dh, 'proper motion distance', r'D_M', 
              filename)
    plot_dist_ony(z, dz, om, dom, dm, dh, 'proper motion distance', r'D_M', 
              filename)

def plot_DA(filename):
    """The dimensionless angular diameter distance DA/DH. 
    """

    # Set up an array of redshift values.
    dz = 0.1
    z = numpy.arange(0., 10. + dz, dz)

    # Set up a cosmology dictionary, with an array of matter density values.
    cosmo = {}
    dom = 0.01
    om = numpy.atleast_2d(numpy.linspace(0.1, 1.0, (1.-0.1)/dom)).transpose()
    cosmo['omega_M_0'] = om
    cosmo['omega_lambda_0'] = 1. - cosmo['omega_M_0']
    cosmo['h'] = 0.701
    cosmo['omega_k_0'] = 0.0

    # Calculate the hubble distance.
    dh = cd.hubble_distance_z(0, **cosmo)
    # Calculate the angular diameter distance.
    da = cd.angular_diameter_distance(z, **cosmo)

    # Make plots.
    plot_dist(z, dz, om, dom, da, dh, 'angular diameter distance', r'D_A',
              filename)
    plot_dist_ony(z, dz, om, dom, da, dh, 'angular diameter distance', r'D_A',
                  filename)

def plot_dist(z, dz, om, dom, dist, dh, name, mathname, filename=None):
    """Make a 2-D plot of a distance versus redshift (x) and matter density (y).
    """
    # Grid of redshift and matter density values.
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


def plot_dist_ony(z, dz, om, dom, dist, dh, name, mathname, filename=None):
    """Make a 2-D plot of matter density versus redshift (x) and distance (y)
    """


    dist = dist/dh
    z = z * numpy.ones(dist.shape)
    om = om * numpy.ones(dist.shape)

    pylab.figure(figsize=(5.5,4.5))    


    pylab.contour(z, dist, om, 50)
    cb = pylab.colorbar()
    cb.ax.set_ylabel(r'$\Omega_M = 1 - \Omega_\lambda$')
    
    pylab.xlim(z.min(), z.max())
    pylab.ylim(dist.min(), dist.max()) 
    pylab.xlabel("redshift z")
    pylab.ylabel(name + r': $'+mathname+'/D_H$')
    pylab.title(name)
    if filename is not None:
        prefix, extension = filename.split('.')
        pylab.savefig(prefix + '_' + mathname + '_ony.' + extension,
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

