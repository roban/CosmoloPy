#!/usr/bin/env python
# Future imports

from __future__ import absolute_import, division, print_function

import re
import subprocess
from setuptools import setup, find_packages, Extension
from os import path
import sys

dirpath = path.abspath(path.dirname(__file__))

eh_dir = path.join(dirpath,'cosmolopy','EH')

ext_names = ['power', 'tf_fit']

def generate_swig():
    print("Swigging sources")
    for d in ext_names:
        filename = path.join(eh_dir, d+'.i')        
        p = subprocess.call(['swig', '-python', filename],
                            cwd=dirpath)
        if p != 0:
            raise RuntimeError("Running swig failed!")

# Generate swig if build is requested
for command in ['develop', 'build', 'bdist_wheel', 'build_ext', 'build_src']:
    if command in sys.argv[1:]:
        generate_swig()
        break

# Generate dict of all data files to include
EH_files = []
if 'sdist' in sys.argv[1:]:
    for name in ext_names:
        EH_files.extend(['%s.c' % (name), '%s.i' % (name)])
package_data = {
    '': ['*.txt', '*.rst', 'LISCENSE'],
    'cosmolopy.EH': EH_files}

# Create extension modules
ext_mods = []
for name in ext_names:
    mod = Extension('cosmolopy.EH._%s' % (name),
                    sources=[path.join(eh_dir, '%s_wrap.c' % (name)),
                             path.join(eh_dir, '%s.c' % (name))])
    ext_mods.append(mod)

# Get the requirements list
with open(path.join(dirpath, 'requirements.txt'), 'r') as f:
    requirements = f.read().splitlines()

# Read the __version__.py file
with open(path.join(dirpath, 'cosmolopy/__version__.py'), 'r') as f:
    vf = f.read()

# Obtain version from read-in __version__.py file
version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", vf, re.M).group(1)

setup(
    name = "cosmolopy",
    version = version,
    packages = find_packages(),
    install_requires = requirements,

    ext_modules = ext_mods,

    tests_require = ['nose', 'matplotlib'],
    test_suite = 'nose.collector',
    platforms=["Windows", "Linux", "Unix"],

    # metadata for upload to PyPI
    author = "Roban Hultman Kramer",
    author_email = "robanhk@gmail.com",
    description = "a cosmology package for Python.",
    url = "http://roban.github.com/CosmoloPy/",   # project home page
    keywords = ("astronomy cosmology cosmological distance density galaxy" +
                "luminosity magnitude reionization Press-Schechter Schecter"),
    license = "MIT",
    python_requires = ">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4",
    package_data=package_data,
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
                   'Programming Language :: Python :: 2.7',
                   'Programming Language :: Python :: 3.5',
                   'Programming Language :: Python :: 3.6',
                   'Programming Language :: Python :: 3.7',                   
                   'Programming Language :: Python :: 2',                   
                   'Programming Language :: Python :: 3',                   
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Operating System :: OS Independent'
                   ]
)
