# -*- coding: utf-8 -*-

"""
Read SFIT4 input ascii files into (from) xarray-compliant dictionnaries.
"""

import datetime
import re
import os
from ast import literal_eval
from collections import OrderedDict

import numpy as np


REF_GAZ_PATTERN = r"\s*(?P<index>\d+)\s+(?P<name>\w+)\s+(?P<description>.*)"


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
        _ = f.readline()

        data = np.loadtxt(f)
        assert data.shape[0] == n_layers
        _, lbound, _, _, points = data.transpose()

        coords = {
            # start level at 1
            #ldim: (ldim, index[:-1].astype('i')),
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

    def get_data(f, count):
        # separators can be space or comma
        current_position = f.tell()
        data = np.fromfile(f, count=count, sep=" ")
        if data.size != count:
            f.seek(current_position)
            data = np.fromfile(f, count=count, sep=",")
            current_position = f.tell()
            if f.read(1) != '\n':
                # no trailing comma
                f.seek(current_position)
        return data

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
            data = get_data(f, n_levels)
            attrs = {'description': description,
                     'ordering': order}
            data_vars['reference__' + profile] = (rdim, data, attrs)

        # gas profiles
        for i in range(n_gases):
            header = re.match(REF_GAZ_PATTERN, f.readline()).groupdict()
            data = get_data(f, n_levels)
            if header['name'] == 'OTHER':
                continue     # ignore unused gas indexes
            attrs = {
                'gas_index': int(header['index']),
                'description': header['description'].strip(),
                'ordering': order
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

            line = f.readline().strip()
            dt_items = list(map(int, re.split(r'[\s\.]+', line[:-1])))
            dt_items[-1] *= 10     # convert decimal seconds to milliseconds
            dt = datetime.datetime(*dt_items)
            # TODO: represent datetime as string for serialization?
            # TODO: address potential time zone issues

            title = f.readline().strip()
            wn_min, wn_max, wn_step, wn_size = np.fromfile(f, count=4, sep=" ")
            data = np.fromfile(f, count=int(wn_size), sep=" ")
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
    # and assume bands and scans are sorted by id in ascending order.
    # scan and band ids begin at 1.
    # flatten scans and bands as 1-d spectral data.
    n_scan, n_band = 1, 1
    i_scan, i_band = 0, 0
    prev_datetime = None
    spec_data, wavenumber, scan, band, snr_ind_vals = [], [], [], [], []
    spec_sza, spec_radius, spec_lat, spec_lon, spec_dt = [], [], [], [], []
    spec_attrs = {}

    for mw in mwindows:
        if mw['datetime'] != prev_datetime:
            # new scan
            i_scan += 1
            n_band = max([n_band, i_band - 1])
            i_band = 1
            spec_attrs['spec_header__scan{}'.format(i_scan)] = mw['title']
            spec_sza.append(mw['solar_zenith_angle'])
            spec_radius.append(mw['earth_radius'])
            spec_lat.append(mw['latitude'])
            spec_lon.append(mw['longitude'])
            spec_dt.append(mw['datetime'])
        snr_ind_vals.append(((i_band - 1, i_scan - 1), mw['initial_snr']))
        spec_data.append(mw['data'])
        wavenumber.append(mw['wavenumber'])
        band.append(np.repeat(i_band, mw['data'].size))
        scan.append(np.repeat(i_scan, mw['data'].size))
        i_band += 1
        prev_datetime = mw['datetime']
    n_band = max([n_band, i_band - 1])
    n_scan = i_scan

    initial_snr = np.full((n_band, n_scan), np.nan)
    for (i, j), v in snr_ind_vals:
        initial_snr[i][j] = v

    coords = {
        wcoord: (spdim, np.concatenate(wavenumber)),
        scoord: (spdim, np.concatenate(scan)),
        bcoord: (spdim, np.concatenate(band)),
        bdim: np.arange(1, n_band + 1),
        sdim: np.arange(1, n_scan + 1)
    }
    data_vars = {
        'spec_observed': (spdim, np.concatenate(spec_data), spec_attrs),
        'spec_snr_initial': ((bdim, sdim), initial_snr),
        'spec_datetime': (sdim, np.asarray(spec_dt)),
        'spec_sza': (sdim, np.asarray(spec_sza)),
        'spec_earth_radius': (sdim, np.asarray(spec_radius)),
        'spec_latitude': (sdim, np.asarray(spec_lat)),
        'spec_longitude': (sdim, np.asarray(spec_lon))
    }

    dataset = {
        'coords': coords,
        'data_vars': data_vars,
        'attrs': global_attrs
    }

    return dataset


def read_ctl(filename, ldim='level'):
    """
    Read control file.

    Use this function to load 'sfit4.ctl'.

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
    variables = OrderedDict()
    attrs = OrderedDict()

    def sanatize_input_name(name):
        return name.replace('.', '__').strip()

    def eval_value(string):
        s = string.strip()
        try:
            return literal_eval(s)
        except:
            if s == 'T':
                return True
            elif s == 'F':
                return False
            re_double_fortran = re.compile(r'(\d*\.\d+)[dD]([-+]?\d+)')
            if re_double_fortran.match(s):
                return float(re_double_fortran.sub(r'\1E\2', s))
            return s

    with open(filename, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            splitted_line = line.split('=')
            if len(splitted_line) < 2:
                splitted_line.append('')
            name, val = splitted_line
            name = sanatize_input_name(name)
            if 'profile' in name and 'sigma' in name:
                # sigma profile
                # assume that 'gas.layers' is already parsed
                if not val.strip():
                    val = np.fromfile(f, count=attrs['gas__layers'],
                                      sep=" ")
                else:
                    val = np.fromstring(val, count=attrs['gas__layers'],
                                        sep=" ")
                variables[name] = (ldim, val)
            else:
                attrs[name] = eval_value(val)

    variables['sfit4_ctl'] = ([], np.nan, attrs)

    dataset = {
        'data_vars': variables,
        'attrs': global_attrs
    }

    return dataset
