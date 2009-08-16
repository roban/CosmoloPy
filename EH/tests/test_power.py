"""Check that the power module is working."""

import numpy.testing.utils as ntest
def test_power():
    """Check that the power module is working."""
    import cosmolopy.EH.power as power
    power.TFmdm_set_cosm(0.3, 0.04, 0.0, 0, 0.7, 0.72, 0.0)
    ntest.assert_approx_equal(power.TFmdm_onek_hmpc(100.0),
                              1.48236347286e-06)
