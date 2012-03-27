"""Check that the tf_fit module is working."""

import numpy.testing.utils as ntest
def test_tf_fit():
    """Check that the tf_fit module is working."""
    import cosmolopy.EH.tf_fit as tf_fit
    tf_fit.TFset_parameters(0.136, 0.2, 2.728)
    ntest.assert_approx_equal( tf_fit.TFfit_onek(10.0),
                              3.69328954548e-05)
