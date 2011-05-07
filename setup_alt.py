#!/usr/bin/env python

######
# This setup.py file does not compile or install the EH perturbation
# theory modole.
######

from setuptools import setup, find_packages
import os
#import nose

packages = find_packages()
setup(
    name = "CosmoloPy",
    version = "0.1.103",
    packages = packages,
    install_requires = ['numpy', 'scipy',],
    tests_require = ['nose', 'matplotlib'],
    test_suite = 'nose.collector',

    # metadata for upload to PyPI
    author = "Roban Hultman Kramer",
    author_email = "robanhk@gmail.com",
    description = "a cosmology package for Python.",
    url = "http://roban.github.com/CosmoloPy/",   # project home page
    keywords = ("astronomy cosmology cosmological distance density galaxy" +
                "luminosity magnitude reionization Press-Schechter Schecter"),
    license = "MIT",
    long_description = \
"""CosmoloPy is a package of cosmology routines built on NumPy/SciPy.

Capabilities include
--------------------

`cosmolopy.density`
  Various cosmological densities.

`cosmolopy.distance`
  Various cosmological distance measures.

`cosmolopy.luminosityfunction`
  Galaxy luminosity functions (Schecter functions).

`cosmolopy.magnitudes`
  Conversion in and out of the AB magnitude system.

`cosmolopy.parameters`
  Pre-defined sets of cosmological parameters (e.g. from WMAP).

`cosmolopy.perturbation`
  Perturbation theory and the power spectrum.

`cosmolopy.reionization`
  The reionization of the IGM.
  
""",
    classifiers = ['License :: OSI Approved :: MIT License',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.6',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Operating System :: OS Independent'
                   ]
)
