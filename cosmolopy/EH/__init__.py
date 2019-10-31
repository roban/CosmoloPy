"""See `cosmolopy.perturbation.transfer_function_EH`

Eisenstein & Hu power spectrum transfer function used in
`cosmolopy.perturbation.transfer_function_EH`.

This package contains two modules:

i. `power` (or `cosmolopy.EH.power`), which is a python wrapper of power.c 
    from Eisenstein & Hu (1999 ApJ 511 5)

ii. `tf_fit` (or `cosmolopy.EH.tf_fit`), which is a python wrapper of
     tf_fit.c from Eisenstein & Hu (1998 ApJ 496 605)

The first module contains no baryonic features in the transfer function, but 
applies to wider range of CDM variants (e.g. including neutrinos). The second
module applies more narrowly but possesses effects including the baryonic
acoustic osciallations.

See also: http://background.uchicago.edu/~whu/transfer/transferpage.html

power.c is redistributed under the MIT License with the permission of
Wayne Hu, and is described in Eisenstein, D. J., & Hu, W. "Power
Spectra for Cold Dark Matter and Its Variants," 1999, ApJ, 511, 5
[astro-ph/9710252]. Please cite this paper to acknowledge use of power
spectrum calculations.

See LISCENSE for additional information.

"""

from __future__ import absolute_import, division, print_function
