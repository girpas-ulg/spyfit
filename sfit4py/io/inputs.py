# -*- coding: utf-8 -*-

"""
Functions for reading/writing the SFIT4 input text files from/into
Python dictionaries and Numpy arrays.

"""

import numpy as np


__all__ = ['read_layers',]


def read_layers(filename):
    """
    Read profile layers.

    Use this function to load 'in.stalayers'.

    """
    with open(filename, 'r') as f:
        outputd = dict()
        outputd['header'] = f.readline().strip()
        n_layers = int(f.readline())
        dummy = f.readline()

        data = np.loadtxt(f)
        assert data.shape[0] == n_layers
        index, lbound, thick, growth, points = data.transpose()

        outputd['index'] = index[:-1].astype('i')
        outputd['altitude'] = points[:-1]
        outputd['altitude_lbound'] = lbound

    return outputd
