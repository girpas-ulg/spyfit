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


def sanatize_name(var_name):
    """Sanatize `var_name` so that it works with autocompletion."""
    return var_name.lower().replace('.', '__')


def read_matrix(filename, var_name='', dims=''):
    """Read a single matrix (or a single vector) in SFIT4 output ascii files.

    Use this function to load 'out.ak_matrix' and 'out.seinv_vector'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the matrix/vector variable. If empty (default),
        the name is set from `filename`.
    dims : tuple or str
        Name of the dimension(s) of the matrix/vector. If empty (default),
        `'x'` or `('x', 'y')` is set.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    if not var_name:
        var_name = sanatize_name(os.path.basename(filename))

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
            'data_vars': {
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
        var_name = sanatize_name(os.path.basename(filename))

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
            'data_vars': {
                var_name: (dims, data, attrs),
            },
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_profiles(filename, var_name_prefix='', ldim='level', ret_gases=False):
    """
    Read a-priori or retrieved profiles in SFIT4 output ascii files.

    Use this function to load 'out.retprofiles' and 'out.aprprofiles'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name_prefix : str
        Prefix to prepend to each of the profile names (e.g., 'apriori' or
        'retrieved') (default: no prefix).
    ldim : str
        Name of the dimension of the profiles (default: 'level').
    ret_gases : bool
        If True, returns the profiles of only the retrieved gases
        (default: False).

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
                if ret_gases and not is_retrieved_gas:
                    continue
                attrs = {'gas_index': gindex,
                         'is_retrieved_gas': is_retrieved_gas}
            else:
                cname = cname.lower()   # lower-case name for non-gases profiles
                attrs = {}
            variables[var_name_prefix + '__' + cname] = ((ldim,), prof, attrs)

        dataset = {
            'data_vars': variables,
            'attrs': global_attrs
        }

    return dataset


def read_state_vector(filename, ldim='level', pdim='param'):
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
        Name of the dimension of the profiles (default: 'level').
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
            variables['apriori__'+gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" "), attrs)

            dummy = f.readline()
            dummy = f.readline()
            #attrs = {'total_column': float(f.readline())}
            variables['retrieved_total_column__' + gas] = (
                (), float(f.readline())
            )
            variables['retrieved__' + gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" ") #, attrs
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
        'data_vars': variables,
        'coords': coords,
        'attrs': global_attrs
    }

    return dataset


def read_param_iterations(filename, vdim='statevector', idim='iteration'):
    """
    Read the state vector for each iteration in SFIT4 ascii output files.

    Use this function to load 'out.parm_vectors'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    vdim : str
        Name of the dimension of the statevector (default: 'statevector').
    idim : str
        Name of the dimension of the iterations (default: 'iteration').

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
            vdim: param_names,
            'statevector_index': (vdim, param_index),
            idim: iterations
        }
        variables = {
            'statevector_iterations': ((idim, vdim), data, attrs)
        }

        dataset = {
            'data_vars': variables,
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_spectra(filename, wdim='wn', bdim='band', sdim='scan',
                 parse_sp_header=None):
    """
    Read observed and fitted spectra in SFIT4 ascii files.

    Use this function to load 'out.pbpfile'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    wdim : str
        Name of the dimension of the wavenumber (default: 'wn'). This
        is actually a prefix to which the band number will be added.
    bdim : str
        Name of the dimension of the bands (default: 'band').
    sdim : str
        Name of the dimension of the scans (default: 'scan').
    parse_sp_header : callable or None
        A callable wich must accept the header line of a spectrum as input
        and must return a dictionary of extracted metadata that will be added in
        the attributes of the observed spectrum entry.
        If None (default), only the plain header line (string) will be added
        in the attributes.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    Notes
    -----
    To avoid duplicates, micro-window metadata will be added only in the
    observed spectrum entries and not in the fitted entries.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = os.path.abspath(filename)

        # n_fits = nb. of bands * nb. of spectra used in each band
        # n_bands = nb. of bands
        n_fits, n_bands = list(map(int, f.readline().split()))

        bands, scans = [], []
        spectrum_header = {}
        wavenumber = {}
        sdim_vars = {
            'spectrum_sza_code': {},
        }
        bdim_vars = {
            'n_retrieved_gas': {}
        }
        spectrum_data = {'observed': {}, 'fitted': {}}

        # loop over each individual fitted spectra (n_fits)
        for i in range(n_fits):

            # fit header line
            line = f.readline()
            iheader = {'spectrum_header': line.strip()}
            if parse_sp_header is not None:
                iheader.update(parse_sp_header(line))

            # fit metadata line
            line = f.readline()
            metadata = list(map(eval, line.split()))
            spec_code, wn_step, size, wn_min, wn_max = metadata[:5]
            u, band_id, scan_id, n_ret_gas = metadata[5:]

            # TODO: refactor! what is u value???
            # TODO: make sure wn_min, wn_max, wn_step do not vary in a band!!
            #       whether the band has one or multiple scans
            # TODO: make sure that one spec_code correspond to one scan id
            # TODO: make sure n_ret_gas is the same for all scans in a band
            spectrum_header[scan_id] = iheader
            sdim_vars['spectrum_sza_code'][scan_id] = spec_code
            bdim_vars['n_retrieved_gas'][band_id] = n_ret_gas
            bands.append(band_id)
            scans.append(scan_id)

            # fit data: 3-line blocks of 12 values for each observed,
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
                dkey = '{}_{}'.format(band_id, scan_id)
                spectrum_data[lbl][dkey] = np.array(data)

            if band_id not in wavenumber.keys():
                wn = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
                assert wn.size == size
                wavenumber[band_id] = wn

        # group data and metadata by band and scan
        unique_bands = sorted(set(bands))
        unique_scans = sorted(set(scans))
        coords = {bdim: unique_bands, sdim: unique_scans}
        for b in unique_bands:
            wdim_band = '{}__band{}'.format(wdim, b)
            coords[wdim_band] = wavenumber[b]

        variables = {}
        for k, v in bdim_vars.items():
            variables[k] = (bdim, [v[b] for b in unique_bands])
        for k, v in sdim_vars.items():
            variables[k] = (sdim, [v[s] for s in unique_scans])

        header_variables = {}
        for s in unique_scans:
            for k, v in spectrum_header[s].items():
                if k in header_variables.keys():
                    header_variables[k].append(v)
                else:
                    header_variables[k] = [v]
        for k, v in header_variables.items():
            variables[k] = (sdim, v)

        for sptype in ('observed', 'fitted'):
            for b in unique_bands:
                dkey = '{}_{}'.format(b, s)
                vname = '{}_spectrum__band{}'.format(sptype, b)
                wdim_band = '{}__band{}'.format(wdim, b)
                data = [spectrum_data[sptype][dkey]
                        if dkey in spectrum_data[sptype].keys()
                        else np.ma.array(np.empty_like(wavenumber[b]),
                                         mask=True)
                        for s in unique_scans]
                variables[vname] = ((sdim, wdim_band), np.ma.row_stack(data))

        dataset = {
            'data_vars': variables,
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def read_single_spectrum(filename, var_name=None, wdim='wn'):
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
        Name chosen for the spectrum variable. If empty, the name is set
        from `filename`. If None (default), the name is defined from the
        spectrum metadata (gas, band, scan and iteration).
    wdim : str
        Name of the dimension of the wavenumber (default: 'wn').
        This is actually the prefix to which the band id will be added.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        global_attrs = {'filename': os.path.abspath(filename)}

        header = f.readline().split()
        gas = header[1].strip()
        band_id, scan_id = int(header[3]), int(header[5])
        iteration = int(header[-1])

        attrs = {'gaz': gas, 'band_id': band_id, 'scan_id': scan_id,
                 'iteration': iteration}

        if var_name is None:
            var_name = "fitted_spectrum__{}__band{}__s{}i{}".format(
                gas, band_id, scan_id, iteration
            )
        if not var_name:
            var_name = sanatize_name(os.path.basename(filename))

        wn_min, wn_max, wn_step, size = map(eval, f.readline().split())
        wavenumber = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
        wdim_band = '{}__band{}'.format(wdim, band_id)
        data = np.loadtxt(f).flatten()
        assert len(wavenumber) == size
        assert len(data) == size

        dataset = {
            'data_vars': {var_name: (wdim_band, data, attrs)},
            'coords': {wdim_band: wavenumber},
            'attrs': global_attrs
        }

    return dataset


def read_solar_spectrum(filename, var_name='', wdim='wn'):
    """
    Read a calculated solar spectrum in SFIT4 output ascii files.

    Use this function to load 'out.solarspectrum'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the spectrum variable. If empty (default),
        the name is set from `filename`.
    wdim : str
        Name of the dimension of the wavenumber (default: 'wn').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = os.path.abspath(filename)

        if not var_name:
            var_name = sanatize_name(os.path.basename(filename))

        size, wn_min, wn_step = map(eval, f.readline().split())

        data = np.loadtxt(f)
        wavenumber = data[:, 0]
        spectrum = data[:, 1]
        assert spectrum.size == size

        dataset = {
            'data_vars': {var_name: ((wdim,), spectrum)},
            'coords': {wdim: wavenumber},
            'attrs': global_attrs
        }

    return dataset


def read_summary(filename):
    """
    Read retrieval summary in SFIT4 output ascii files.

    Use this function to load 'out.summary'.

    Parameters
    ----------
    filename : str
        Name or path to the file.

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

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
