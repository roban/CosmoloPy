#!/usr/bin/env python
from setuptools import setup, find_packages, Extension
import os
import nose

eh_dir = os.path.join('.','cosmolopy','EH')

power_module = Extension('cosmolopy.EH._power',
                         sources=[os.path.join(eh_dir, 'power.i'),
                                  os.path.join(eh_dir, 'power.c')]
                         )

setup(
    name = "CosmoloPy",
    version = "0.0.alpha",
    packages = find_packages(),
    package_data = {
        # If any package contains *.txt or *.rst files, include them:
        '': ['*.so'],
        },
    install_requires = ['numpy', 'scipy',],

    ext_modules = [power_module],

    tests_require = ['nose',],
    test_suite = 'nose.collector',

    # metadata for upload to PyPI
    author = "Roban Hultman Kramer",
    author_email = "robanhk@gmail.com",
    description = "a basic NumPy/SciPy-based cosmology package for Python.",
    license = "MIT",
    url = "http://roban.github.com/CosmoloPy/",   # project home page
)
