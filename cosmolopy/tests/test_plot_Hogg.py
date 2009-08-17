"""Reproduce some plots from David Hogg's arXiv:astro-ph/9905116v4

"""
import inspect

import numpy
import matplotlib.pyplot as pylab

import cosmolopy.distance as cd
import cosmolopy.constants as cc

def test_figure1():
    """Plot Hogg fig. 1: The dimensionless proper motion distance DM/DH. 

    The three curves are for the three world models, Einstein-de
    Sitter (omega_M, omega_lambda) = (1, 0), solid; low-density,
    (0.05, 0), dotted; and high lambda, (0.2, 0.8), dashed.

    Hubble distance DH = c / H0

    z from 0--5
    DM / DH from 0--3

    """

    z = numpy.arange(0, 5.05, 0.05)

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[1.0],[0.05],[0.2]])
    cosmo['omega_lambda_0'] = numpy.array([[0.0],[0.0],[0.8]])
    cosmo['h'] = 0.5
    
    linestyle = ['-', ':', '--']

    dh = cd.hubble_distance_z(0, **cosmo)
    dm, dm_err = cd.comoving_distance_transverse(z, **cosmo)

    pylab.figure(figsize=(6,6))    
    for i in range(len(linestyle)):
        pylab.plot(z, (dm/dh)[i], ls=linestyle[i])
        #pylab.plot(z, (dm_err/dh)[i], ls=linestyle[i])
    pylab.xlim(0,5)
    pylab.ylim(0,3)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"proper motion distance $D_M/D_H$")
    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " (astro-ph/9905116v4)")

def test_figure2():
    """Plot Hogg fig. 2: The dimensionless angular diameter distance DA/DH.

    The three curves are for the three world models, 

    - Einstein-de Sitter (omega_M, omega_lambda) = (1, 0) [solid]
    
    : Low-density (0.05, 0) [dotted]

    -- High lambda, (0.2, 0.8) [dashed]

    Hubble distance DH = c / H0

    z from 0--5
    DA / DH from 0--0.5

    """

    z = numpy.arange(0, 5.05, 0.05)

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[1.0],[0.05],[0.2]])
    cosmo['omega_lambda_0'] = numpy.array([[0.0],[0.0],[0.8]])
    cosmo['h'] = 0.5
    
    linestyle = ['-', ':', '--']

    dh = cd.hubble_distance_z(0, **cosmo)
    da, err1, err2 = cd.angular_diameter_distance(z, **cosmo)

    # Also test the pathway with non-zero z0
    da2, err12, err22 = cd.angular_diameter_distance(z, z0=1e-8, **cosmo)

    pylab.figure(figsize=(6,6))
    for i in range(len(linestyle)):
        pylab.plot(z, (da/dh)[i], ls=linestyle[i])
        pylab.plot(z, (da2/dh)[i], ls=linestyle[i])
    pylab.xlim(0,5)
    pylab.ylim(0,0.5)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"angular diameter distance $D_A/D_H$")
    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " (astro-ph/9905116v4)")

def test_figure3():
    """Plot Hogg fig. 3: The dimensionless luminosity distance DL/DH

    The three curves are for the three world models, 

    - Einstein-de Sitter (omega_M, omega_lambda) = (1, 0) [solid]
    
    : Low-density (0.05, 0) [dotted]

    -- High lambda, (0.2, 0.8) [dashed]

    Hubble distance DH = c / H0

    z from 0--5
    DL / DH from 0--16

    """

    z = numpy.arange(0, 5.05, 0.05)

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[1.0],[0.05],[0.2]])
    cosmo['omega_lambda_0'] = numpy.array([[0.0],[0.0],[0.8]])
    cosmo['h'] = 0.5
    
    linestyle = ['-', ':', '--']

    dh = cd.hubble_distance_z(0, **cosmo)
    dl, err1, err2 = cd.luminosity_distance(z, **cosmo)

    pylab.figure(figsize=(6,6))
    for i in range(len(linestyle)):
        pylab.plot(z, (dl/dh)[i], ls=linestyle[i])
    pylab.xlim(0,5)
    pylab.ylim(0,16)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"luminosity distance $D_L/D_H$")
    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " (astro-ph/9905116v4)")

def test_figure6():
    """Plot Hogg fig. 6: The dimensionless lookback time t_L/t_H and age t/t_H.

    The three curves are for the three world models, 

    - Einstein-de Sitter (omega_M, omega_lambda) = (1, 0) [solid]
    
    : Low-density (0.05, 0) [dotted]

    -- High lambda, (0.2, 0.8) [dashed]

    Hubble distance DH = c / H0

    z from 0--5
    t/th from 0--1.2

    """

    z = numpy.arange(0, 5.05, 0.05)

    cosmo = {}
    cosmo['omega_M_0'] = numpy.array([[1.0],[0.05],[0.2]])
    cosmo['omega_lambda_0'] = numpy.array([[0.0],[0.0],[0.8]])
    cosmo['h'] = 0.5
    
    linestyle = ['-', ':', '--']

    th = 1/ cd.hubble_z(0, **cosmo)

    tl, err = cd.lookback_time(z, **cosmo)
    age, err_f, err_t = cd.age(z, **cosmo)

    pylab.figure(figsize=(6,6))
    for i in range(len(linestyle)):
        pylab.plot(z, (tl/th)[i], ls=linestyle[i])
        pylab.plot(z, (age/th)[i], ls=linestyle[i])
    pylab.xlim(0,5)
    pylab.ylim(0,1.2)
    pylab.xlabel("redshift z")
    pylab.ylabel(r"lookback timne $t_L/t_H$")
    pylab.title("compare to " + inspect.stack()[0][3].replace('test_', '') + 
                " (astro-ph/9905116v4)")

if __name__ == "__main__":
    test_figure1()
    test_figure2()
    test_figure3()
    test_figure6()
    pylab.show()

