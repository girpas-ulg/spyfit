# -*- coding: utf-8 -*-

"""
The `io` sub-package provides functions to read (write) SFIT4
input and output files into (from) very basic, common data structures
(i.e., Numpy arrays and Python dictionaries).

TODO: all functions should use a consistent data model here.
Maybe follow a 'minimal' netCDF data model (coordinates, variables and
attributes). See how `xray` datasets dictionary arguments must be arranged.
proposition: return dict structured as
d = {
    'coordinates': [
        {'name': '', ''points': array, 'bounds': array,
         'attributes': {'name'...}},
    ],
    'variables': [
        {'name': '', 'data': array, 'dim_coords': ((0, 'name'), (1, 'name')...),
         'attributes' : {...}},
    ],
    'global_attributes': {}
}

"""

from .inputs import *
from .outputs import *
