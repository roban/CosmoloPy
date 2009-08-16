#!/usr/bin/env python
from setuptools import setup, find_packages
setup(
    name = "CosmoloPy",
    version = "0.0.alpha",
    packages = find_packages(),
    install_requires = ['numpy', 'scipy', 'nose'],

    # metadata for upload to PyPI
    author = "Roban Hultman Kramer",
    author_email = "robanhk@gmail.com",
    description = "a basic NumPy/SciPy-based cosmology package for Python.",
    license = "MIT",
    url = "http://roban.github.com/CosmoloPy/",   # project home page
)
