=========
CosmoloPy 
=========

A cosmology package for Python.

For documentation and installation instructions, see
http://roban.github.com/CosmoloPy/

CosmoloPy is released under the MIT software liscense (see LISCENSE).

Example
-------

Calculate the comoving distance to redshift 6.

 >>> import cosmolopy.distance as cd
 >>> cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.72}
 >>> d_co = cd.comoving_distance(6., **cosmo)
 >>> print "Comoving distance to z=6 is %.1f Mpc" % (d_co)
 Comoving distance to z=6 is 8017.8 Mpc


Prerequisites
=============

  Python
  NumPy
  SciPy

For tests (optional):
  nose
  matplotlib

For power spectrum calculation (needed for most of perturbation module):
  python-dev 

Installation from PyPI
======================

If you use Python 3.5–3.7 on Linux, you can easily install the package directly from the Python Package
Index with pip.

Run with:

    > pip install cosmolopy

Installation from Source
========================

If you use a different operating system or a more recent version of Python,
you first need to install `SWIG <https://github.com/swig/swig/>`_ v3 or later.
Then, download the CosmoloPy source and
install it by running (in CosmoloPy folder)

    > pip install . 

Testing
=======

Note that currently the _udot integration appears to be failing in nosetests.
The prefered way to run all tests is:

    > python setup.py nosetests --with-doctest

Or if you don't want to run the doctests:

    > python setup.py nosetests

If you don't have nose:

    > python setup.py test
    > python -m doctest cosmolopy/*.py

Contributors
============

- Python 3 implementation by @JohannesBuchner and @1313e 
- Automated PyPI deployment on Travis by @1313e
