# -*- coding: utf-8 -*-

"""
Functions for reading/writing the SFIT4 input text files from/into
Python dictionaries and Numpy arrays.

"""

import re
import datetime

import numpy as np


__all__ = ['read_layers', 'read_reference_profiles', 'read_spectrum']


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
        order_desc, n_levels, n_gases = map(int, f.readline().split())
        if order_desc:
            order = 'descending'
        else:
            order = 'ascending'

        # altitude, pressure and temperature profiles
        for profile in ('altitude', 'pressure', 'temperature'):
            description = f.readline().strip()
            data = np.fromfile(f, count=n_levels, sep=" ")
            outputd[profile] = {
                'points': data,
                'attributes': {
                    'description': description,
                    'order': order
                }
            }

        # gas profiles
        gases = []
        for iprofile in range(n_gases):
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
                    'order': order
                }
            }
            gases.append(gas)

        outputd['gases'] = gases

    return outputd


def read_spectrum(filename):
    """
    Read a spectrum (bands and scans data).

    Use this function to load 'in.spectrum'.

    """
    with open(filename, 'r') as f:
        spectra = []

        while True:
            line = f.readline()
            if not line:
                break
            sza, earth_radius, lat, lon, snr = map(float, line.split())

            dt_items = list(map(int, re.split('[\s\.]+', f.readline())[:-1]))
            dt_items[-1] *= 10     # convert decimal seconds to milliseconds
            dt = datetime.datetime(*dt_items)

            title = f.readline().strip()

            wn_min, wn_max, wn_step, size = map(eval, f.readline().split())

            data = np.fromfile(f, count=size, sep=" ")
            wavenumber = np.arange(wn_min, wn_max + wn_step / 2., wn_step)

            spectra.append({
                'data': data,
                'wavenumber': wavenumber,  # TODO: treat it as a coordinate
                'attributes': {
                    'title': title,
                    'solar_zenith_angle': sza,
                    'earth_radius': earth_radius,
                    'latitude': lat,
                    'longitude': lon,
                    'snr': snr,
                    'datetime': dt
                }
            })

    return {'spectra': spectra}
