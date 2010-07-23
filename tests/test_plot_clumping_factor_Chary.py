import sys

import numpy
import pylab

import cosmolopy.reionization as cr

def plot_clumping_factor_Chary():
    """Plot clumping factor from Chary paper.
    """
    z = numpy.linspace(5,15,500)
    cf = cr.clumping_factor_Chary(z)

    pylab.plot(z,cf)
    pylab.yscale('log')
    pylab.xlim(5,15)
    pylab.ylim(1,100)

if __name__ == "__main__":
    if len(sys.argv)==1:
        print "Run with a filename argument to produce image files, e.g.:"
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = None
        
    plot_clumping_factor_Chary()

    ### Plot output code. ###
    if filename is None:
        pylab.show()
    else:
        from matplotlib import _pylab_helpers
        for manager in _pylab_helpers.Gcf.get_all_fig_managers():
            fig = manager.canvas.figure
            if len(fig.get_label()) > 0:
                label = fig.get_label()
            else:
                label = 'clumpFig' + str(manager.num)
            newfilename = prefix + '_' + label + extension
            fig.savefig(newfilename, bbox_inches="tight")
