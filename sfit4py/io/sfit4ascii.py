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


def parse_header(line):
    """Parse the header line of an output file."""
    m = re.match(HEADER_PATTERN, line)
    header = m.groupdict()

    runtime = datetime.datetime.strptime(header.pop('runtime'),
                                         "%Y%m%d-%H:%M:%S")
    header['sfit4_runtime'] = str(runtime)
    header['description'] = header['description'].strip().lower()

    return header


def read_matrix(filename, var_name='', dims=()):
    """Read a single matrix (or a single vector) in SFIT4 output ascii files.

    Use this function to load 'out.ak_matrix' and 'out.seinv_vector'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the matrix/vector variable. If empty (default),
        the name is set from `filename`.
    dims : tuple
        Name of the dimension(s) of the matrix/vector. If empty (default),
        `('x',)` or `('x', 'y')` is set.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    if not var_name:
        var_name = os.path.basename(filename)

    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        attrs = {'description': header.pop('description'),
                 'source': os.path.abspath(filename)}
        global_attrs = header

        data_shape = tuple([int(d) for d in f.readline().split() if int(d) > 1])
        data = np.loadtxt(f)
        assert data.shape == data_shape

        if not len(dims):
            dims = ('x', 'y')[:data.ndim]

        dataset = {
            'variables': {
                var_name: (dims, data, attrs),
            },
            'attrs': global_attrs
        }

    return dataset


def read_table(filename, var_name='', dims=()):
    """
    Read (labeled) tabular data in SFIT4 output ascii files.

    Use this function to load 'out.k_matrix', 'out.g_matrix', 'out.kb_matrix',
    'out.sa_matrix', 'out.sainv_matrix' and 'out.shat_matrix'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the tabular variable. If empty (default),
        the name is set from `filename`.
    dims : tuple
        Name of the dimension(s) of the table. If empty (default),
        `('rows', 'cols')` is set.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    if not var_name:
        var_name = os.path.basename(filename)

    if not len(dims):
            dims = ('rows', 'cols')

    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        attrs = {'description': header.pop('description'),
                 'source': os.path.abspath(filename)}
        global_attrs = header

        nrows, ncols = [int(n) for n in f.readline().split()[:2]]
        col_names = [c.strip() for c in f.readline().split()]
        coords = {dims[1]: col_names}

        data = np.loadtxt(f)
        assert data.shape == (nrows, ncols)

        if len(col_names) != ncols:
            assert len(col_names) == nrows
            data = data.transpose()

        dataset = {
            'variables': {
                var_name: (dims, data, attrs),
            },
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_profiles(filename, var_name_prefix='', ldim='levels'):
    """
    Read a-priori or retrieved profiles in SFIT4 output ascii files.

    Use this function to load 'out.retprofiles' and 'out.aprprofiles'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name_prefix : str
        Prefix to prepend to each of the profile names (e.g., 'apriori_' or
        'retrieved') (default: no prefix).
    ldim : str
        Name of the dimension of the profiles (default: 'levels').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = {'source': os.path.abspath(filename)}
        global_attrs.update(header)

        meta = f.readline().split()
        nrows = int(meta[1])
        retrieved_gases = [g.strip() for g in meta[3:]]
        gas_index = list(map(int, f.readline().split()))
        col_names = [c.strip() for c in f.readline().split()]
        # TODO: redefine (and translate) the first 5 column names (coordinates)

        data = np.loadtxt(f)
        assert data.shape == (nrows, len(col_names))

        variables = {}
        for cname, gindex, prof in zip(col_names, gas_index, data.transpose()):
            if cname == 'OTHER':
                continue
            if gindex:
                is_retrieved_gas = cname in retrieved_gases
                attrs = {'gas_index': gindex,
                         'is_retrieved_gas': is_retrieved_gas}
            else:
                attrs = {}
            variables[var_name_prefix + cname] = ((ldim,), prof, attrs)

        dataset = {
            'variables': variables,
            'attrs': global_attrs
        }

    return dataset


def read_state_vector(filename, ldim='levels', pdim='param'):
    """
    Read the state vector in SFIT4 output ascii files.

    The state vector includes a-priori and retrieved profiles
    - with calculated total columns - and extra parameters.

    Use this function to load 'out.statevec'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    ldim : str
        Name of the dimension of the profiles (default: 'levels').
    pdim : str
        Name of the dimension of the parameters (default: 'param').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = os.path.abspath(filename)

        meta = f.readline().split()
        nlevels, niter, nitermax = list(map(int, meta[:3]))
        is_temp, has_converged, has_divwarn = list(map(
            lambda s: True if s == "T" else False, meta[3:]
        ))
        global_attrs['n_iteration'] = niter
        global_attrs['n_iteration_max'] = nitermax
        global_attrs['is_temp'] = is_temp   # TODO: refactor ! what is 'istemp'?
        global_attrs['has_converged'] = has_converged
        global_attrs['has_division_warnings'] = has_divwarn

        # altitude is a coordinate
        coords = {}
        dummy = f.readline()
        coords['altitude'] = ((ldim,), np.fromfile(f, count=nlevels, sep=" "))

        # apriori profiles of pressure and temperature
        variables = {}
        for p in ('apriori_pressure', 'apriori_temperature'):
            dummy = f.readline()
            variables[p] = ((ldim,), np.fromfile(f, count=nlevels, sep=" "))

        # apriori/retrieved gas profiles (and columns)
        for i in range(int(f.readline())):
            dummy = f.readline()
            gas = f.readline().strip()
            attrs = {'total_column': float(f.readline())}
            variables['apriori_'+gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" "), attrs)

            dummy = f.readline()
            dummy = f.readline()
            attrs = {'total_column': float(f.readline())}
            variables['retrieved_'+gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" "), attrs
            )

        # parameters
        nparams = int(f.readline())

        #    np.fromfile raises a SystemError for parameter names
        pnames = []
        while 1:
            if len(pnames) >= nparams:
                break
            pnames += [n.strip() for n in f.readline().split()]
        coords[pdim] = pnames

        variables['apriori_parameters'] = (
            (pdim,), np.fromfile(f, count=nparams, sep=" ")
        )
        variables['retrieved_parameters'] = (
            (pdim,), np.fromfile(f, count=nparams, sep=" ")
        )

    dataset = {
        'variables': variables,
        'coords': coords,
        'attrs': global_attrs
    }

    return dataset


def read_param_iterations(filename, sdim='statevector', idim='iterations'):
    """
    Read the state vector for each iteration in SFIT4 ascii output files.

    Use this function to load 'out.parm_vectors'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    sdim : str
        Name of the dimension of the statevector (default: 'statevector').
    idim : str
        Name of the dimension of the iterations (default: 'iterations').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        attrs = {'description': header.pop('description')}
        global_attrs = header
        global_attrs['source'] = os.path.abspath(filename)

        n_params = int(f.readline())
        param_index = list(map(int, f.readline().split()))
        param_names = [n.strip() for n in f.readline().split()]
        assert len(param_index) == n_params
        assert len(param_names) == n_params

        raw_data = np.loadtxt(f)
        iterations = raw_data[:, 0].astype('i')
        data = raw_data[:, 1:]

        coords = {
            sdim: param_names,
            'statevector_index': (sdim, param_index),
            idim: iterations
        }
        variables = {
            'statevector_iterations': ((idim, sdim), data, attrs)
        }

        dataset = {
            'variables': variables,
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_spectra(filename, wdim='wavenumber', parse_mw_header=None):
    """
    Read observed and fitted spectra in SFIT4 ascii files.

    Use this function to load 'out.pbpfile'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    wdim : str
        Name of the dimension of the wavenumber (default: 'wavenumber'). This
        is actually the prefix to which the name of the micro-window will be
        added.
    parse_mw_header : callable or None
        A callable wich must accept the header line of a micro-window as input
        and must return a dictionary of extracted metadata that will be added in
        the attributes of the observed micro-window entry.
        If None (default), only the plain header line (string) will be added
        in the attributes.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    Notes
    -----
    To avoid duplicates, micro-window metadata will be added only in the
    observed micro-window entries and not in the fitted entries.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = os.path.abspath(filename)

        # n_fits = nb. of micro-windows * nb. of spectra
        # n_mw = nb. of micro windows
        n_fits, n_mw = list(map(int, f.readline().split()))

        # loop over each individual micro-windows (n_fits)
        coords = {}
        variables = {}
        for i in range(n_fits):

            # micro-window header line
            line = f.readline()
            attrs = {'header_line': line.strip()}
            if parse_mw_header is not None:
                attrs.update(parse_mw_header(line))

            # micro-window metadata line
            line = f.readline()
            metadata = list(map(eval, line.split()))
            spec_code, wn_step, size, wn_min, wn_max = metadata[:5]
            # TODO: refactor! what is u value???
            u, band_id, scan_id, n_ret_gas = metadata[5:]
            attrs['spectrum_code'] = spec_code
            attrs['n_retrieval_gas'] = n_ret_gas
            attrs['band_id'], attrs['scan_id'] = band_id, scan_id

            w_suffix = "_window_s{}b{}".format(scan_id, band_id)
            wdim_suffix = wdim + w_suffix

            # micro-window data: 3-line blocks of 12 values for each observed,
            # fitted and difference spectra.
            # 1st value of 1st line is the wavenumber (ignored, re-calculated)
            # difference spectra can be easily calculated, it is not returned.
            n_vals_line = 12
            labels = ('observed', 'fitted', 'difference')
            slices = [slice(1, None), slice(None), slice(None)]
            w_data = [list(), list(), list()]

            for block in range(int(math.ceil(1. * size / n_vals_line))):
                for s, data in zip(slices, w_data):
                    data += list(map(float, f.readline().split()))[s]

            for lbl, data in zip(labels, w_data):
                if lbl == 'difference':
                    continue
                a = attrs
                if lbl == 'fitted':
                    a = {}
                variables[lbl + w_suffix] = (wdim_suffix, np.array(data), a)

            wavenumber = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
            assert wavenumber.size == size
            coords[wdim_suffix] = wavenumber

        dataset = {
            'variables': variables,
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_single_spectrum(filename, var_name='', wdim='wavenumber'):
    """
    Read a single spectrum in SFIT4 output ascii files.

    a single spectrum stored in one file, e.g., for a given gaz, band, scan
    and at a given iteration.

    Use this function to load 'out.gas_spectra' files.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the spectrum variable. If empty (default),
        the name is set from `filename`.
    wdim : str
        Name of the dimension of the wavenumber (default: 'wavenumber').
        This is actually the prefix to which the name of the micro-window
        will be added.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        global_attrs = {'filename': os.path.abspath(filename)}

        if not var_name:
            # TODO: provide a function to sanatize var_name
            # e.g., replace '.' by '_' to work with xray and autocompletion
            var_name = os.path.basename(filename)

        header = f.readline().split()
        gas = header[1].strip()
        band_id, scan_id = int(header[3]), int(header[5])
        iteration = int(header[-1])

        attrs = {'gaz': gas, 'band_id': band_id, 'scan_id': scan_id,
                 'iteration': iteration}

        w_suffix = "_window_s{}b{}".format(scan_id, band_id)
        wdim += w_suffix

        wn_min, wn_max, wn_step, size = map(eval, f.readline().split())
        wavenumber = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
        data = np.loadtxt(f).flatten()
        assert len(wavenumber) == size
        assert len(data) == size

        dataset = {
            'variables': {var_name: (wdim, data, attrs)},
            'coords': {wdim: wavenumber},
            'attrs': global_attrs
        }

    return dataset


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
