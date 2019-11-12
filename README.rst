=========
CosmoloPy 
=========

a cosmology package for Python

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

You can easily install the pacakge directly from the Python Package
Index with pip.

Run with:

    > pip install cosmolopy

Installation from Source
========================

If you've downloaded the source, first install Swig-3 or later and then 
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
