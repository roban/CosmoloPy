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
    version = "0.1.101",
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
  Calculate various cosmological densities.

`cosmolopy.distance`
  Calculate various cosmological distance measures.

`cosmolopy.luminosityfunction`
  Routines related to galaxy luminosity functions (Schecter functions).

`cosmolopy.magnitudes`
  Simple routines for conversion in and out of the AB magnitude system.

`cosmolopy.parameters`
  Provides some pre-defined sets of cosmological parameters (e.g. from WMAP).

`cosmolopy.perturbation`
  Routines related to perturbation theory and the power spectrum.

`cosmolopy.reionization`
  Routines related to the reionization of the IGM.
  
""",
    classifiers = ['License :: OSI Approved :: MIT License',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 2.6',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Operating System :: OS Independent'
                   ]
)
