# -*- coding: utf-8 -*-

"""Low level IO routines.

The `io` sub-package provides functions to read/write various formats
(e.g., SFIT4 I/O ascii files, GEOMS hdf4) into/from very basic,
common data structures (i.e., Python dictionaries and Numpy arrays).

The returned dictionaries have the following, xray-friendly structure,
which follow the Common Data Model (CDM)::

    d = {
        'variables': {
            'var1': (('x', 'y'), data, {'units': 'kg',}),
        },
        'coords': {
            'x': (xdata, {'units': 'cm'}),
            'y': (ydata, {'units': 'cm'}),
        },
        'attrs': {}
    }

"""

from .sfit4ascii import *
