#!/usr/bin/env python
from setuptools import setup, find_packages, Extension
import os
import nose

eh_dir = os.path.join('.','cosmolopy','EH')

# Stuff used to build the cosmolopy.EH._power module:
power_module = Extension('cosmolopy.EH._power',
                         sources=[os.path.join(eh_dir, 'power.i'),
                                  os.path.join(eh_dir, 'power.c')]
                         )

setup(
    name = "CosmoloPy",
    version = "0.1.100",
    packages = find_packages(),
    package_data = {
        # If any package contains *.so files, include them:
        '': ['*.so'],
        },
    install_requires = ['numpy', 'scipy',],

    ext_modules = [power_module],

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
