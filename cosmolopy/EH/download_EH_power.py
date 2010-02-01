"""Download the power.c power spectrum code from Wayne Hu's website.

For more information see:

 http://background.uchicago.edu/~whu/transfer/transferpage.html

 Eisenstein & Hu (1999ApJ...511....5E, or astro-ph/9710252)

The url for the code is:

 http://background.uchicago.edu/~whu/transfer/power.c

"""


import urllib
import os

def download_power():
    # Find the directory this script resides in.
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Name for saved file.
    target = os.path.join(current_dir, 'power.c')

    url = r'http://background.uchicago.edu/~whu/transfer/power.c'
    
    print __doc__

    urllib.urlretrieve(url, filename=target)

    print "Downloaded to %s" % target

if __name__ == "__main__":
    download_power()
