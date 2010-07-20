import numpy

def fractional_diff_string(actual, desired, places=3):
    formats = "% ." + ("%i" % places) + "g"
    template = "frac. diff. of "+formats+" and "+formats+" = "+formats
    return template % (actual, desired, (desired - actual)/desired)


