# -*- coding: utf-8 -*-

"""
The `io` sub-package provides functions to read (write) SFIT4
input and output files into (from) very basic, common data structures
(i.e., Numpy arrays and Python dictionaries).

TODO: all functions should use a consistent data model here.
Maybe follow a 'minimal' netCDF data model (coordinates, variables and
attributes). See how `xray` datasets dictionary arguments must be arranged.

"""

from .inputs import *
from .outputs import *
