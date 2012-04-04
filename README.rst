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

Installation of Prerequisites under Ubuntu
------------------------------------------

    > sudo apt-get install python-numpy python-scipy python-matplotlib 
    > sudo apt-get install python-dev

Installation from PYPI
======================

You can easily install the pacakge directly from the Python Package
Index with easy_install or pip.

Run either:

    > sudo easy_install cosmolopy

Or:

    > sudo pip install cosmolopy

Installation from Source
========================

If you've downloaded the source, install it by running

    > sudo python setup.py install

If you have trouble compiling from source, it's probably the EH.power
module (everything else is pure python). You can install without it by
running:

    > sudo python setup_alt.py install

Testing
=======

The prefered way to run all tests is:

    > python setup.py nosetests --with-doctest

Or if you don't want to run the doctests:

    > python setup.py nosetests

If you don't have nose:

    > python setup.py test
    > python -m doctest cosmolopy/*.py
