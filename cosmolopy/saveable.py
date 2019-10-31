"""A Saveable class with methods to save and restore.

Saveable is designed to be subclassed to create new types of objects
that can easily be pickled and reloaded.
"""

from __future__ import absolute_import, division, print_function

import pickle

class NullWriter:
    """Dummy file-like object that does nothing.

    From
    http://stackoverflow.com/questions/1809958/hide-stderr-output-in-unit-tests

    Used in Saveable.
    """
    def write(self, s):
        pass


def loadSaveable(filename):
    """Return an instance of an object unpickled from a file.
    """
    picfile = open(filename)
    loaded = pickle.load(picfile)
    picfile.close()
    return loaded

class Saveable(object):
    """An object with methods to save and restore.

    Unpickleable attributes will simply be deleted from the object.
    """

    def __init__(self, filename=None):
        if filename is not None:
            self.load(filename)
                  
    def save(self, filename):
        """Save object to a file."""
        picfile = open(filename, 'w')
        pickle.dump(self, picfile)
        picfile.close()

    def load(self, filename):
        """Return an instance of an object unpickled from a file.
        """
        picfile = open(filename)
        loaded = pickle.load(picfile)
        picfile.close()
        self.__dict__.update(loaded.__dict__)
        return self

    def dumb_pickle_filter(self):
        """Filter out attributes that can't be pickled.

        Returns a copy of dict with 
        """
        picfile = NullWriter()
        sdict = self.__dict__.copy()
        for k, v in list(sdict.items()):
            # Avoid self references (and thereby infinite recurion).
            if v is self:
                del sdict[k]
                continue
            # Remove any attributes that can't be pickled.
            try:
                pickle.dump(v, picfile)
            except (TypeError, pickle.PicklingError) as err:
                if hasattr(self, 'verbose') and self.verbose:
                    print("Won't pickle", k, type(v), ": ")
                    print("'", err, "'")
                del sdict[k]
        return sdict
    def __getstate__(self):
        """Prepare a state of pickling."""
        return self.dumb_pickle_filter()

    def __setstate__(self, dict):
        """Unpickle."""
        self.__dict__.update(dict)
