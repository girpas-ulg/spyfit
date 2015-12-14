# -*- coding: utf-8 -*-

"""Read (write) SFIT4 input and output ascii files."""

import datetime
import re
import os
import math

import numpy as np


__all__ = ['read_matrix', 'read_table', 'read_profiles',
           'read_state_vector', 'read_param_iterations',
           'read_spectra', 'read_single_spectrum', 'read_solar_spectrum',
           'read_summary', 'read_layers', 'read_reference_profiles',
           'read_spectrum']


# ---------------------------------------
# output files
# ---------------------------------------

HEADER_PATTERN = (r"\s*SFIT4:V(?P<sfit4_version>[0-9.]+)"
                  r"[\w\s:-]*RUNTIME:(?P<runtime>[0-9:\-]+)"
                  r"\s*(?P<description>.+)")

# TODO: check resolution ?
MW_HEADER_PATTERN = (r"\s*(?P<date>[0-9/]+),"
                     r"\s*(?P<time>[0-9:]+),"
                     r"\s*SCAN TIME:(?P<scan_time>[0-9.]*)\s*SEC\."
                     r"\s*RES=\s*(?P<resolution>[0-9.]*).+")


def parse_header(line):
    """
    Parse the header line of an output file.

    """
    m = re.match(HEADER_PATTERN, line)
    header = m.groupdict()

    header['runtime'] = datetime.datetime.strptime(
        header['runtime'], "%Y%m%d-%H:%M:%S"
    )
    header['description'] = header['description'].strip().lower()

    return header


def parse_mw_header(line):
    """
    Parse the header line of a micro-window.

    """
    m = re.match(MW_HEADER_PATTERN, line)
    if m is None:
        return {'failed to parse header': None}
    header = m.groupdict()

    header['datetime'] = datetime.datetime.strptime(
        header.pop('date') + header.pop('time'), "%d/%m/%Y%H:%M:%S"
    )
    for k in ('scan_time', 'resolution'):
        header[k] = float(header[k])

    return header


def read_matrix(filename):
    """
    Read single matrix (or single vector) output files.

    Use this function to load 'out.ak_matrix' and 'out.seinv_vector'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        data_shape = tuple([int(d) for d in f.readline().split() if int(d) > 1])
        data = np.loadtxt(f)
        assert data.shape == data_shape
        outputd['data'] = data

    return outputd


def read_table(filename):
    """
    Read (labeled) tabular output data.

    Use this function to load 'out.k_matrix', 'out.g_matrix', 'out.kb_matrix',
    'out.sa_matrix', 'out.sainv_matrix' and 'out.shat_matrix'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        nrows, ncols = [int(n) for n in f.readline().split()[:2]]
        outputd['column_names'] = [c.strip() for c in f.readline().split()]

        data = np.loadtxt(f)
        assert data.shape == (nrows, ncols)

        if len(outputd['column_names']) != ncols:
            assert len(outputd['column_names']) == nrows
            data = data.transpose()

        outputd['data'] = data

    return outputd


def read_profiles(filename):
    """
    Read a-priori or retrieved profiles.

    Use this function to load 'out.retprofiles' and 'out.aprprofiles'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        meta = f.readline().split()
        nrows = int(meta[1])
        outputd['retrieved_gas'] = [g.strip() for g in meta[3:]]
        outputd['gas_index'] = list(map(int, f.readline().split()))
        outputd['column_names'] = [c.strip() for c in f.readline().split()]
        # TODO: redefine (and translate) the first 5 column names (coordinates)

        data = np.loadtxt(f)
        assert data.shape == (nrows, len(outputd['column_names']))
        outputd['data'] = data

    return outputd


def read_state_vector(filename):
    """
    Read the state vector (a-priori and retrieved profiles
    - with calculated total columns - and extra parameters).

    Use this function to load 'out.statevec'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        meta = f.readline().split()
        nlevels, niter, nitermax = list(map(int, meta[:3]))
        is_temp, has_converged, has_divwarn = list(map(
            lambda s: True if s == "T" else False, meta[3:]
        ))
        outputd['n_iteration'] = niter
        outputd['n_iteration_max'] = nitermax
        outputd['is_temp'] = is_temp      # TODO: refactor ! what is 'istemp'?
        outputd['has_converged'] = has_converged
        outputd['has_division_warnings'] = has_divwarn

        # coords
        for c in ('altitude', 'pressure', 'temperature'):
            dummy = f.readline()
            outputd[c] = np.fromfile(f, count=nlevels, sep=" ")

        # gas profiles / columns
        for k in ('apriori_profiles', 'apriori_total_columns',
                  'retrieved_profiles', 'retrieved_total_columns'):
            outputd[k] = dict()

        for i in range(int(f.readline())):
            dummy = f.readline()
            gas = f.readline().strip()
            outputd['apriori_total_columns'][gas] = float(f.readline())
            outputd['apriori_profiles'][gas] = np.fromfile(f, count=nlevels,
                                                           sep=" ")
            dummy = f.readline()
            dummy = f.readline()
            outputd['retrieved_total_columns'][gas] = float(f.readline())
            outputd['retrieved_profiles'][gas] = np.fromfile(f, count=nlevels,
                                                             sep=" ")

        # other parameters
        nparams = int(f.readline())

        #    np.fromfile raises a SystemError for parameter names
        pnames = []
        while 1:
            if len(pnames) >= nparams:
                break
            pnames += [n.strip() for n in f.readline().split()]
        outputd['param_names'] = pnames

        outputd['apriori_params'] = np.fromfile(f, count=nparams, sep=" ")
        outputd['retrieved_params'] = np.fromfile(f, count=nparams, sep=" ")

    return outputd


def read_param_iterations(filename):
    """
    Read the state vector factors (parameters) for each iteration.

    Use this function to load 'out.parm_vectors'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        n_params = int(f.readline())
        outputd['param_index'] = list(map(int, f.readline().split()))
        outputd['param_names'] = [n.strip() for n in f.readline().split()]
        assert len(outputd['param_index']) == n_params
        assert len(outputd['param_names']) == n_params

        raw_data = np.loadtxt(f)
        outputd['iterations'] = raw_data[:, 0].astype('i')
        outputd['data'] = raw_data[:, 1:]

    return outputd


def read_spectra(filename):
    """
    Read observed, fitted and difference spectra.

    Use this function to load 'out.pbpfile'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        # n_fits = nb. of micro-windows * nb. of spectra
        # n_mw = nb. of micro windows
        n_fits, n_mw = list(map(int, f.readline().split()))
        micro_windows = []

        # loop over each individual micro-windows (n_fits)
        for i in range(n_fits):
            mw = dict()

            # micro-window header line
            line = f.readline()
            mw['header'] = line.strip()
            mw.update(parse_mw_header(line))

            # micro-window metadata line
            line = f.readline()
            metadata = list(map(eval, line.split()))
            spec_code, wn_step, size, wn_min, wn_max = metadata[:5]
            # TODO: refactor! what is u value???
            u, band_id, scan_id, n_ret_gas = metadata[5:]
            mw['spectrum_code'] = spec_code
            mw['n_retrieval_gas'] = n_ret_gas
            mw['band_id'], mw['scan_id'] = band_id, scan_id

            # micro-window data: 3-line blocks of 12 values for each observed,
            # fitted and difference spectra.
            # 1st value of 1st line is the wavenumber (ignored, re-calculated)
            n_vals_line = 12
            labels = ('observed', 'fitted', 'difference')
            slices = [slice(1, None), slice(None), slice(None)]
            mw_data = [list(), list(), list()]

            for block in range(int(math.ceil(1. * size / n_vals_line))):
                for s, data in zip(slices, mw_data):
                    data += list(map(float, f.readline().split()))[s]

            for lbl, data in zip(labels, mw_data):
                mw[lbl] = np.array(data)

            mw['wavenumber'] = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
            assert mw['wavenumber'].size == size

            micro_windows.append(mw)

        outputd['micro_windows'] = micro_windows

    return outputd


def read_single_spectrum(filename):
    """
    Read a single spectrum (i.e., for a given gaz, band, scan and
    at a given iteration).

    Use this function to load 'out.gas_spectra' files.

    """
    with open(filename, 'r') as f:
        outputd = dict()
        outputd['filename'] = os.path.abspath(filename)

        header = f.readline().split()
        outputd['gas'] = header[1].strip()
        outputd['band_id'], outputd['scan_id'] = int(header[3]), int(header[5])
        outputd['iteration'] = int(header[-1])

        wn_min, wn_max, wn_step, size = map(eval, f.readline().split())
        outputd['wavenumber'] = np.arange(wn_min, wn_max + wn_step / 2.,
                                          wn_step)
        outputd['data'] = np.loadtxt(f).flatten()
        assert len(outputd['wavenumber']) == size
        assert len(outputd['data']) == size

    return outputd


def read_solar_spectrum(filename):
    """
    Read a calculated solar spectrum.

    Use this function to load 'out.solarspectrum'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        size, wn_min, wn_step = map(eval, f.readline().split())

        data = np.loadtxt(f)
        outputd['wavenumber'] = data[:, 0]
        outputd['data'] = data[:, 1]
        assert outputd['data'].size == size

    return outputd


def read_summary(filename):
    """
    Read retrieval summary.

    Use this function to read 'out.summary'.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(f.readline())
        outputd['filename'] = os.path.abspath(filename)

        # micro-window headers
        dummy = f.readline()
        mw_headers = []
        for i in range(int(f.readline())):
            mw_headers.append(f.readline().strip())
        outputd['micro_window_headers'] = mw_headers

        # retrieved gases
        dummy = f.readline()
        s2bool_func = lambda s: True if s == 'T' else False
        cfuncs = [int, lambda s: s, s2bool_func, float, float]
        ckeys = ['index', 'name', 'has_retrieved_profile',
                 'apriori_total_column', 'retrieved_total_column']
        ret_gases = []
        n_gases = int(f.readline())
        dummy = f.readline()
        for i in range(n_gases):
            cvals = [s.strip() for s in f.readline().split()]
            ret_gases.append({k: f(v) for k, f, v in zip(ckeys, cfuncs, cvals)})
        outputd['retrieved_gases'] = ret_gases

        # bands
        dummy = f.readline()
        icfuncs = [int, float, float, float, int, float, float, float, int]
        jcfuncs = [int, float, float]
        ickeys = ['index', 'wavenumber_start', 'wavenumber_end',
                  'wavenumber_step', 'n_points', 'pmax', 'fovdia']
        jckeys = ['index', 'initial_snr', 'calculated_snr']
        bands = []
        n_bands = int(f.readline())
        dummy = f.readline()
        for i in range(n_bands):
            icvals = [s.strip() for s in f.readline().split()]
            d = {k: f(v) for k, f, v in zip(ickeys, icfuncs, icvals[:-1])}
            scans = []
            for j in range(int(icvals[-1])):
                jcvals = [s.strip() for s in f.readline().split()]
                scans.append(
                    {k: f(v) for k, f, v in zip(jckeys, jcfuncs, jcvals)}
                )
            d['scans'] = scans
            bands.append(d)
        outputd['bands'] = bands

        # other returned values
        dummy = f.readline()
        cfuncs = [lambda s: float(s) / 100, float, float, float, float,
                  int, int, s2bool_func, s2bool_func]
        ckeys = ['fit_rms', 'chi_square_obs', 'dofs_total',
                 'dofs_trg', 'dofs_tpr', 'n_iterations', 'n_iterations_max',
                 'has_converged', 'has_division_warnings']
        dummy = f.readline()
        cvals = [s.strip() for s in f.readline().split()]
        outputd.update({k: f(v) for k, f, v in zip(ckeys, cfuncs, cvals)})

    return outputd


# ---------------------------------------
# input files
# ---------------------------------------

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
                'attrs': {
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
                'attrs': {
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
                'attrs': {
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
