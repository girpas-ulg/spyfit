# -*- coding: utf-8 -*-

"""
Functions for reading/writing the SFIT4 input text files from/into
Python dictionaries and Numpy arrays.

"""

import re

import numpy as np


__all__ = ['read_layers', 'read_reference_profiles']


REF_GAZ_PATTERN = r"\s*(?P<index>\d+)\s+(?P<name>\w+)\s+(?P<description>.+)"


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
        # TODO: follow this way of storing coordinates in all read functions
        # see in `xray` how dictionnary arguments should be arranged for
        # coordinates and variables
        outputd['altitude'] = {
            'points': points[:-1],
            'lower_bounds': lbound
        }

    return outputd


def read_reference_profiles(filename):
    """
    Read coordinates, atmosphere and gas reference profiles.

    Use this function to load 'in.refprofile'.

    """
    with open(filename, 'r') as f:
        outputd = dict()
        sorting, n_levels, n_gases = map(int, f.readline().split())
        outputd['sort_descending'] = bool(sorting)

        # altitude, pressure and temperature profiles
        for k in ('altitude', 'pressure', 'temperature'):
            description = f.readline().strip()
            data = np.fromfile(f, count=n_levels, sep=" ")
            outputd[k] = {'points': data}
            outputd[k]['points'] = data
            outputd[k]['attributes'] = {'description': description}

        # gas profiles
        gases = []
        for i in range(n_gases):
            header = re.match(REF_GAZ_PATTERN, f.readline()).groupdict()
            data = np.fromfile(f, count=n_levels, sep=" ")
            if header['name'] == 'OTHER':
                data = np.array([])
            gas = {
                'name': '{}_reference_profile'.format(header['name']),
                'data': data,
                'attributes': {
                    'gas': header['name'],
                    'gas_index': int(header['index']),
                    'description': header['description'].strip(),
                }
            }
            gases.append(gas)

        outputd['gases'] = gases

    return outputd
