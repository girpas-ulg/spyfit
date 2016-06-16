# -*- coding: utf-8 -*-

"""
Read (write) SFIT4 input ascii files into/from xarray-compliant
dictionnaries.
"""

import datetime
import re
import os
import math
from itertools import groupby

import numpy as np


__all__ = ['read_layers', 'read_ref_profiles', 'read_spectrum']


REF_GAZ_PATTERN = r"\s*(?P<index>\d+)\s+(?P<name>\w+)\s+(?P<description>.+)"


def read_layers(filename, ldim='level'):
    """
    Read profile layers.

    Use this function to load 'in.stalayers'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    ldim : str
        Name of the dimension of the profiles (default: 'level').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.


    """
    global_attrs = {'source': os.path.abspath(filename)}

    with open(filename, 'r') as f:
        attrs = {'header': f.readline().strip()}
        n_layers = int(f.readline())
        dummy = f.readline()

        data = np.loadtxt(f)
        assert data.shape[0] == n_layers
        index, lbound, thick, growth, points = data.transpose()

        coords = {
            ldim: (ldim, index[:-1].astype('i')),
        }
        variables = {
            'station_altitude': (ldim, points[:-1], attrs),
            'station_altitude_lbound': (ldim + '_lbound', lbound)
        }

    dataset = {'data_vars': variables, 'coords': coords, 'attrs': global_attrs}
    return dataset


def read_ref_profiles(filename, rdim='rlevel'):
    """
    Read coordinates, atmosphere and gas reference profiles.

    Use this function to load 'in.refprofile'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    rdim : str
        Name of the dimension of the input reference profiles
        (default: 'rlevel').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    global_attrs = {'source': os.path.abspath(filename)}

    with open(filename, 'r') as f:
        order_desc, n_levels, n_gases = map(int, f.readline().split())
        if order_desc:
            order = 'descending'
        else:
            order = 'ascending'

        data_vars = {}

        # altitude, pressure and temperature profiles
        for profile in ('altitude', 'pressure', 'temperature'):
            description = f.readline().strip()
            data = np.fromfile(f, count=n_levels, sep=" ")
            attrs = {'description': description,
                     'order': order}
            data_vars['reference__' + profile] = (rdim, data, attrs)

        # gas profiles
        for iprofile in range(n_gases):
            header = re.match(REF_GAZ_PATTERN, f.readline()).groupdict()
            data = np.fromfile(f, count=n_levels, sep=" ")
            if header['name'] == 'OTHER':
                continue     # ignore unused gas indexes
            attrs = {
                'gas_index': int(header['index']),
                'description': header['description'].strip(),
                'order': order
            }
            vname = 'reference__{}'.format(header['name'])
            data_vars[vname] = (rdim, data, attrs)

    dataset = {'data_vars': data_vars, 'attrs': global_attrs}
    return dataset


def read_spectrum(filename, spdim='spectrum', bdim='band', sdim='scan',
                  wcoord='spec_wn', scoord='spec_scan', bcoord='spec_band'):
    """
    Read a spectrum (bands and scans data).

    Use this function to load 'in.spectrum'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    spdim : str
        Name of the dimension of the spectral data (default: 'spectra').
        spectral data for all micro-windows (i.e., bands and scans)
        will be flattened as a 1-d array.
    wcoord : str
        Name of the wavenumber coordinate for spectral data
        (default: 'spec_wn').
    scoord : str
        Name of the coordinate for spectral scans (default: 'spec_scan').
    bcoord : str
        Name of the coordinate for spectral bands (default: 'spec_band').
    bdim : str
        Name of the dimension of the band (default: 'band').
    sdim : str
        Name of the dimension of the scan (default: 'scan').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    global_attrs = {'source': os.path.abspath(filename)}

    with open(filename, 'r') as f:
        mwindows = []

        while True:
            line = f.readline()
            if not line:
                break
            sza, earth_radius, lat, lon, snr = map(float, line.split())

            dt_items = list(map(int, re.split('[\s\.]+', f.readline())[:-1]))
            dt_items[-1] *= 10     # convert decimal seconds to milliseconds
            dt = datetime.datetime(*dt_items)
            # TODO: represent datetime as string for serialization

            title = f.readline().strip()

            wn_min, wn_max, wn_step, wn_size = map(eval, f.readline().split())

            data = np.fromfile(f, count=wn_size, sep=" ")
            wn = np.arange(wn_min, wn_max + wn_step / 2., wn_step)

            mwindows.append({
                'data': data,
                'wavenumber': wn,
                'initial_snr': snr,  # SNR different for each band / scan
                'title': title,
                'solar_zenith_angle': sza,
                'earth_radius': earth_radius,
                'latitude': lat,
                'longitude': lon,
                'datetime': dt
            })

    # no band or scan id information given here.
    # assume micro windows order in file = for s in scans for b in bands
    # and assume bands and scans are sorted by id.
    # flatten scans and bands as 1-d spectral data.
    n_scan, n_band = 1, 1
    i_scan, i_band = 1, 1
    prev_datetime = None
    spec_data, wavenumber, scan, band, snr_ind_vals = [], [], [], [], []

    for mw in mwindows:
        if (i_band > 1 or i_scan > 1) and mw['datetime'] != prev_datetime:
            # new scan
            i_scan += 1
            n_band = max([n_band, i_band - 1])
            i_band = 1
        snr_ind_vals.append(((i_band, i_scan), mw['initial_snr']))
        spec_data.append(mw['data'])
        wavenumber.append(mw['wavenumber'])
        band.append(np.repeat(i_band, mw['data'].size))
        scan.append(np.repeat(i_scan, mw['data'].size))
        i_band += 1
        prev_datetime = mw['datetime']
    n_scan = i_scan
    n_band = max([n_band, i_band - 1])

    initial_snr = np.full((n_band, n_scan), np.nan)
    for (i, j), v in snr_ind_vals:
        initial_snr[i - 1][j - 1] = v

    coords = {
        wcoord: (spdim, np.concatenate(wavenumber)),
        scoord: (spdim, np.concatenate(scan)),
        bcoord: (spdim, np.concatenate(band)),
        bdim: np.arange(1, n_band + 1),
        sdim: np.arange(1, n_scan + 1)
    }
    data_vars = {
        'spec_observed': (spdim, np.concatenate(spec_data)),
        'initial_snr': ((bdim, sdim), initial_snr)
    }

    dataset = {
        'coords': coords,
        'data_vars': data_vars,
        'attrs': global_attrs
    }

    return dataset
