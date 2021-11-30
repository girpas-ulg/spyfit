# -*- coding: utf-8 -*-

"""Read SFIT4 output ascii files into xarray-compliant dictionnaries."""

import datetime
import re
import os
import math
import glob
import warnings

import numpy as np

from .utils import expand_path, sanatize_var_name, eval_value


HEADER_PATTERN = (r"\s*SFIT4:V(?P<sfit4_version>[0-9.]+)"
                  r"[\w\s:-]*RUNTIME:(?P<runtime>[0-9:\-]+)"
                  r"\s*(?P<description>.+)")


def parse_header(line,ignore_case=False):
    """Parse the header line of an output file."""
    if ignore_case:
        m = re.match(HEADER_PATTERN, line,re.IGNORECASE)
    else:
        m = re.match(HEADER_PATTERN, line)
    
    header = m.groupdict()

    runtime = datetime.datetime.strptime(header.pop('runtime'),
                                         "%Y%m%d-%H:%M:%S")
    header['sfit4_runtime'] = str(runtime)
    header['description'] = header['description'].strip().lower()

    return header


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
        var_name = sanatize_var_name(os.path.basename(filename))

    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        attrs = {'description': header.pop('description'),
                 'source': expand_path(filename)}
        global_attrs = header

        data_shape = tuple([int(d) for d in f.readline().split()
                            if int(d) > 1])
        data = np.fromfile(f, sep=" ").reshape(data_shape)

        if not len(dims):
            dims = ('x', 'y')[:data.ndim]

        dataset = {
            'data_vars': {
                var_name: (dims, data, attrs),
            },
            'attrs': global_attrs
        }

    return dataset


def read_ak_matrix(filename, ldim='level', kdim='kernel'):
    """Read averaging kernel matrix.

    Use this function to load 'out.ak_matrix'.

    See Also
    --------
    :func:`read_matrix`

    """
    return read_matrix(filename, var_name='averaging_kernels',
                       dims=(kdim, ldim))


def read_seinv_vector(filename, ddim='diag'):
    """Read seinv.

    Use this function to load 'out.seinv_vector'.

    See Also
    --------
    :func:`read_matrix`

    """
    return read_matrix(filename, var_name='seinv', dims=ddim)


def read_table(filename, var_name='', dims=(), index_cols=True):
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
        var_name = sanatize_var_name(os.path.basename(filename))

    if not len(dims):
        dims = ('rows', 'cols')

    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        attrs = {'description': header.pop('description'),
                 'source': expand_path(filename)}
        global_attrs = header

        nrows, ncols = [int(n) for n in f.readline().split()[:2]]
        index = [c.strip() for c in f.readline().split()]

        data = np.loadtxt(f)
        assert data.shape == (nrows, ncols)

        if index_cols:
            assert len(index) == ncols
            coord_name = dims[1]
        else:
            assert len(index) == nrows
            coord_name = dims[0]
        coords = {coord_name: index}

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
        global_attrs = {'source': expand_path(filename)}
        global_attrs.update(header)

        meta = f.readline().split()
        nrows = int(meta[1])
        retrieved_gases = [g.strip() for g in meta[3:]]
        gas_index = [int(i) for i in f.readline().split()]
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
                cname = cname.lower()  # lower-case name for non-gases profiles
                attrs = {}
            variables[var_name_prefix + '__' + cname] = ((ldim,), prof, attrs)

        dataset = {
            'data_vars': variables,
            'attrs': global_attrs
        }

    return dataset


def read_aprfs(filename, ldim='level'):
    """
    Read a-priori profiles in SFIT4 output ascii files.

    Use this function to load 'out.aprprofiles'.

    See Also
    --------
    :func:`read_profiles`

    """
    return read_profiles(filename, ldim=ldim, var_name_prefix='apriori_',
                         ret_gases=False)


def read_rprfs(filename, ldim='level'):
    """
    Read retrieved profiles in SFIT4 output ascii files.

    Use this function to load 'out.aprprofiles'.

    See Also
    --------
    :func:`read_profiles`

    """
    return read_profiles(filename, ldim=ldim, var_name_prefix='retrieved_',
                         ret_gases=True)


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
    variables = {}
    coords = {}

    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = expand_path(filename)

        meta = f.readline().split()
        nlevels, niter, nitermax = list(map(int, meta[:3]))
        is_temp, has_converged, has_divwarn = list(map(
            lambda s: True if s == "T" else False, meta[3:]
        ))
        attrs = {'n_iteration': niter,
                 'n_iteration_max': nitermax,
                 'is_temp': is_temp,   # TODO: refactor ! what is 'istemp'?
                 'has_converged': has_converged,
                 'has_division_warnings': has_divwarn}
        variables['sfit4_summary'] = ((), np.nan, attrs)

        # altitude is a coordinate
        _ = f.readline()
        variables['station_altitude'] = (
            (ldim,), np.fromfile(f, count=nlevels, sep=" ")
        )

        # apriori profiles of pressure and temperature
        for p in ('apriori_pressure', 'apriori_temperature'):
            _ = f.readline()
            variables[p] = ((ldim,), np.fromfile(f, count=nlevels, sep=" "))

        # apriori/retrieved gas profiles (and columns)
        for i in range(int(f.readline())):
            _ = f.readline()
            gas = f.readline().strip()
            variables['apriori_total_column__' + gas] = (
                (), float(f.readline())
            )
            variables['apriori__' + gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" ")
            )

            _ = f.readline()
            _ = f.readline()
            variables['retrieved_total_column__' + gas] = (
                (), float(f.readline())
            )
            variables['retrieved__' + gas] = (
                (ldim,), np.fromfile(f, count=nlevels, sep=" ")
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
        global_attrs['source'] = expand_path(filename)

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


def read_spectra(filename, spdim='spectrum', idim='iteration',
                 wcoord='spec_wn', scoord='spec_scan', bcoord='spec_band',
                 parse_sp_header=None):
    """
    Read observed and fitted spectra in SFIT4 ascii files.

    Use this function to load 'out.pbpfile'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    spdim : str
        Name of the dimension of the spectral data (default: 'spectra').
        spectral data for all fitted micro-windows (i.e., bands and scans)
        will be flattened as a 1-d array.
    idim : str
        Name of the iteration dimension (default: 'iteration'). This is not
        very useful here as 'out.pbpfile' stores spectral data for only the
        last iteration but this is to be consistent with data returned by
        the `read_single_spectrum` and `read_single_spectra` functions.
    wcoord : str
        Name of the wavenumber coordinate for spectral data
        (default: 'spec_wn').
    scoord : str
        Name of the coordinate for spectral scans (default: 'spec_scan').
    bcoord : str
        Name of the coordinate for spectral bands (default: 'spec_band').
    parse_sp_header : callable or None
        A callable wich must accept the header line of a spectrum as input
        and must return a dictionary of extracted metadata that will be added
        in the attributes of the observed spectrum entry.
        If None (default), only the plain header line (string) will be added
        in the attributes.

    Returns
    -------
    datase t : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = expand_path(filename)

        # n_fits = nb. of bands * nb. of spectra/scans used in each band
        n_fits, _ = list(map(int, f.readline().split()))

        spec_data = {'observed': [], 'fitted': []}
        coords_data = {wcoord: [], bcoord: [], scoord: []}
        spec_attrs = {}
        data_size = 0

        # loop over each individual fitted spectra (n_fits)
        for i in range(n_fits):

            # fit header line
            scan_hdr = f.readline().strip()
            # TODO: add parsed header elements as separate attrs
            # iheader = {'spectrum_header': line.strip()}
            # if parse_sp_header is not None:
            #    iheader.update(parse_sp_header(line))

            # fit metadata line
            line = f.readline()
            metadata = list(map(eval, line.split()))
            _, wn_step, size, wn_min, wn_max = metadata[:5]
            _, band_id, scan_id, n_ret_gas = metadata[5:]

            # TODO: refactor! what is u value???
            # TODO: make sure wn_min, wn_max, wn_step do not vary in a band!!
            #       whether the band has one or multiple scans
            # TODO: make sure that one spec_code correspond to one scan id
            # TODO: make sure n_ret_gas is the same for all scans in a band
            data_size += size
            spec_attrs['spec_header__scan{}'.format(scan_id)] = scan_hdr
            spec_attrs['n_retrieved_gas__band{}'.format(band_id)] = n_ret_gas

            # fit data: 3-line blocks of 12 values for each observed,
            # fitted and difference spectra.
            # 1st value of 1st line is the wavenumber (ignored, re-calculated)
            # skip difference spectrum as can be easily re-calculated.
            n_vals_line = 12
            labels = ('observed', 'fitted', 'difference')
            slices = [slice(1, None), slice(None), slice(None)]
            w_data = [list(), list(), list()]

            for block in range(int(math.ceil(1. * size / n_vals_line))):
                for s, data in zip(slices, w_data):
                    data += [float(i) for i in f.readline().split()][s]

            for lbl, data in zip(labels, w_data):
                if lbl == 'difference':
                    continue
                spec_data[lbl].append(np.array(data))

            wn = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
            assert wn.size == size
            coords_data[wcoord].append(wn)
            coords_data[bcoord].append(np.repeat(band_id, size))
            coords_data[scoord].append(np.repeat(scan_id, size))

        coords = {k: (spdim, np.concatenate(v))
                  for k, v in coords_data.items()}
        coords[idim] = (idim, np.array([-1]))

        variables = {
            'spec_observed': (spdim, np.concatenate(spec_data['observed']),
                              spec_attrs),
            'spec_fitted__ALL': (
                (idim, spdim),
                np.concatenate(spec_data['fitted'])[np.newaxis, :]
            )
        }

        dataset = {
            'data_vars': variables,
            'coords': coords,
            'attrs': global_attrs
        }

    return dataset


def _read_single_spec(filename):
    with open(filename, 'r') as f:
        header = f.readline().split()
        gas = header[1].strip()
        band_id, scan_id = int(header[3]), int(header[5])
        iteration = int(header[-1])

        wn_min, wn_max, wn_step, size = map(eval, f.readline().split())
        wavenumber = np.arange(wn_min, wn_max + wn_step / 2., wn_step)
        data = np.loadtxt(f).flatten()
        assert len(wavenumber) == size
        assert len(data) == size

    return gas, band_id, scan_id, iteration, wavenumber, data


def read_single_spectrum(filename, var_name=None, wdim='spec_wn'):
    """
    Read a single spectrum in a SFIT4 output ascii files.

    Read one spectrum that is stored in a single file, i.e.,
    for a given gas, band, scan and at a given iteration.

    Use this function to load 'out.gas_spectra' files.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str or None
        Name chosen for the spectrum variable. If empty, the name is set
        from `filename`. If None (default), the name is defined from the
        spectrum metadata (gas, band, scan and iteration).
    wdim : str
        Name of the dimension of the wavenumber (default: 'spec_wn').
        This is actually the prefix to which the band id will be added.

    Returns
    -------
    dataarray : dict
        A CDM-structured dictionary.

    """
    raw_data = _read_single_spec(filename)
    gas, band_id, scan_id, iteration, wavenumber, data = raw_data
    attrs = {'source': expand_path(filename), 'gas': gas,
             'band_id': band_id, 'scan_id': scan_id,
             'iteration': iteration}

    if var_name is None:
        var_name = "spec_fitted__{}__band{}__scan{}__iter{}".format(
            gas, band_id, scan_id, iteration
        )
    if not var_name:
        var_name = sanatize_var_name(os.path.basename(filename))

    wdim_band = '{}__band{}'.format(wdim, band_id)

    dataarray = {
        'data': data,
        'coords': {wdim_band: wavenumber},
        'dims': (wdim_band,),
        'name': var_name,
        'attrs': attrs
    }

    return dataarray


def read_single_spectra(filename, spdim='spectrum', idim='iteration',
                        wcoord='spec_wn', scoord='spec_scan',
                        bcoord='spec_band'):
    """
    Read single spectra in a SFIT4 output ascii files.

    Read one or more spectra that are each stored in a single file, i.e.,
    for a given gas, band, scan and at a given iteration.

    Use this function to load all 'out.gas_spectra' files at once.

    Parameters
    ----------
    filename : str
        Name or path to the file(s). Support UNIX-like file specification
        (e.g., using '*' to specify all files that match a given pattern).
    spdim : str
        Name of the dimension of the spectral data (default: 'spectrum').
        spectral data or coordinates for all fitted micro-windows
        (i.e., bands and scans) will be flattened as a 1-d array.
    idim : str
        Name of the dimension of the iterations (default: 'iteration').
    wcoord : str
        Name of the wavenumber coordinate for spectral data
        (default: 'spec_wn').
    scoord : str
        Name of the coordinate for spectral scans (default: 'spec_scan').
    bcoord : str
        Name of the coordinate for spectral bands (default: 'spec_band').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    global_attrs = {'source': expand_path(filename)}

    file_data = [_read_single_spec(f) for f in glob.glob(filename)]
    if not file_data:
        raise ValueError('no single spectrum file found in {}'
                         .format(expand_path(filename)))

    keys = ('gas', 'band', 'scan', 'iteration', 'wn', 'spec')
    data_vals = {k: np.asarray(v) for k, v in zip(keys, zip(*file_data))}
    unique_vals = {k: np.unique(v) for k, v in data_vals.items()
                   if k not in ('wn', 'spec')}

    def _get_sb_loc(s, b):
        indexer = ((data_vals['scan'] == s) & (data_vals['band'] == b))
        loc = indexer.nonzero()[0]
        if loc.size > 0:
            return loc[0]
        else:
            return loc

    def _get_loc(g, i, s, b):
        indexer = ((data_vals['gas'] == g) & (data_vals['iteration'] == i) &
                   (data_vals['scan'] == s) & (data_vals['band'] == b))
        loc = indexer.nonzero()[0]
        if not loc.size:
            # check if files exist for that band and scan and for
            # other gases and/or iterations
            if np.any((data_vals['scan'] == s) & (data_vals['band'] == b)):
                warnings.warn("expected output spectrum for gas {}, "
                              "scan {}, band {} & iteration {} but "
                              "no file found".format(g, s, b, i))
                return None
        if loc.size == 1:
            return loc[0]
        elif loc.size > 1:
            raise ValueError('duplicate file found for spectrum '
                             '{} (iteration: {}, scan: {}, band: {})'
                             .format(g, i, s, b))
        return loc

    band_size = {b: data_vals['spec'][data_vals['band'] == b][0].size
                 for b in unique_vals['band']}
    sb_size = [(s, b, data_vals['spec'][_get_sb_loc(s, b)].size)
               for s in unique_vals['scan'] for b in unique_vals['band']]

    def _get_flat_sorted(data, g, i, fallback_wn=False):
        sel_data = []
        for s in unique_vals['scan']:
            for b in unique_vals['band']:
                loc = _get_loc(g, i, s, b)
                if loc is None:
                    if fallback_wn:
                        # get wavenumber values from another gas / iteration
                        sel_data.append(data[_get_sb_loc(s, b)])
                    else:
                        sel_data.append(np.full(band_size[b], np.nan))
                else:
                    sel_data.append(data[loc])
        return np.concatenate(sel_data)

    spec_data = {}
    for g in unique_vals['gas']:
        spec_data[g] = np.vstack([_get_flat_sorted(data_vals['spec'], g, i)]
                                 for i in unique_vals['iteration'])

    wavenumber = _get_flat_sorted(data_vals['wn'], unique_vals['gas'][0],
                                  unique_vals['iteration'][0],
                                  fallback_wn=True)
    scan = np.concatenate([np.repeat(s, size) for s, b, size in sb_size])
    band = np.concatenate([np.repeat(b, size) for s, b, size in sb_size])

    coords = {wcoord: (spdim, wavenumber),
              scoord: (spdim, scan),
              bcoord: (spdim, band),
              idim: unique_vals['iteration']}
    variables = {'spec_fitted__{}'.format(k): ((idim, spdim), v)
                 for k, v in spec_data.items()}

    dataset = {
        'data_vars': variables,
        'coords': coords,
        'attrs': global_attrs
    }

    return dataset


def read_solar_spectrum(filename, var_name='', wdim='spec_wn_solar'):
    """
    Read a calculated sol ar spectrum in SFIT4 output ascii files.

    Use this function to load 'out.solarspectrum'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    var_name : str
        Name chosen for the spectrum variable. If empty (default),
        the name is set from `filename`.
    wdim : str
        Name of the dimension of the wavenumber (default: 'spec_wn_solar').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline())
        global_attrs = header
        global_attrs['source'] = expand_path(filename)

        if not var_name:
            var_name = sanatize_var_name(os.path.basename(filename))

        size, _, wn_step = map(eval, f.readline().split())

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


def read_summary(filename, spdim='spectrum', wcoord='spec_wn',
                 scoord='spec_scan', bcoord='spec_band',
                 bdim='band', sdim='scan',ignore_case_header=False):
    """
    Read retrieval summary in SFIT4 output ascii files.

    Use this function to load 'out.summary'.

    Parameters
    ----------
    filename : str
        Name or path to the file.
    spdim : str
        Name of the dimension of the spectral data (default: 'spectrum').
        spectral data or coordinates for all fitted micro-windows
        (i.e., bands and scans) will be flattened as a 1-d array.
    wcoord : str
        Name of the wavenumber coordinate for spectral data
        (default: 'spec_wn').
    scoord : str
        Name of the coordinate for spectral scans (default: 'spec_scan').
    bcoord : str
        Name of the coordinate for spectral bands (default: 'spec_band').
    bdim : str
        Name of the dimension of the bands (default: 'band').
    sdim : str
        Name of the dimension of the scans (default: 'scan').

    Returns
    -------
    dataset : dict
        A CDM-structured dictionary.

    """
    with open(filename, 'r') as f:
        header = parse_header(f.readline(),ignore_case=ignore_case_header)
        global_attrs = header
        global_attrs['source'] = expand_path(filename)

        variables = {}
        coords = {}

        # spectrum headers (assume same header for every band)
        _ = f.readline()
        spec_headers = []
        for i in range(int(f.readline())):
            header = f.readline().strip()
            if header not in spec_headers:
                spec_headers.append(header)

        # retrieved gases
        _ = f.readline()
        n_gases = int(f.readline())
        _ = f.readline()
        for i in range(n_gases):
            cvals = [eval_value(s) for s in f.readline().split()]
            index, name, ret_prof, atcol, rtcol = cvals
            variables['apriori_total_column__' + name] = ((), atcol)
            attrs = {'has_retrieved_profile': ret_prof, 'index': index}
            variables['retrieved_total_column__' + name] = ((), rtcol, attrs)

        # bands and scans
        _ = f.readline()
        ickeys = ['index', 'wn_start', 'wn_end', 'wn_step',
                  'n_points', 'spec_opd_max', 'spec_fov']
        jckeys = ['index', 'spec_snr_initial', 'spec_snr_calculated']
        bands = []
        n_bands = int(f.readline())
        _ = f.readline()
        for i in range(n_bands):
            icvals = [s.strip() for s in f.readline().split()]
            d = {k: eval_value(v) for k, v in zip(ickeys, icvals[:-1])}
            scans = []
            for j in range(int(icvals[-1])):
                jcvals = [s.strip() for s in f.readline().split()]
                scans.append(
                    {k: eval_value(v) for k, v in zip(jckeys, jcvals)}
                )
            d['scans'] = scans
            bands.append(d)

        wn = [None] * (len(bands) + 1)
        for b in bands:
            wn[b['index']] = np.arange(b['wn_start'],
                                       b['wn_end'] + b['wn_step'] / 2.,
                                       b['wn_step'])

        coords[sdim] = np.unique([s['index']
                                  for b in bands for s in b['scans']])
        coords[bdim] = np.unique([b['index'] for b in bands])

        wavenumber = []
        scan = []
        band = []
        for s in coords[sdim]:
            for b in coords[bdim]:
                band_scans = bands[b - 1]['scans']
                if s not in [v['index'] for v in band_scans]:
                    continue
                wavenumber.append(wn[b])
                scan.append(np.repeat(s, wn[b].size))
                band.append(np.repeat(b, wn[b].size))

        coords[wcoord] = (spdim, np.concatenate(wavenumber))
        coords[scoord] = (spdim, np.concatenate(scan))
        coords[bcoord] = (spdim, np.concatenate(band))

        for vname in ('spec_fov', 'spec_opd_max'):
            variables[vname] = (bdim, np.array([b[vname] for b in bands]))

        for vname in ('spec_snr_initial', 'spec_snr_calculated'):
            data = np.full((coords[bdim].size, coords[sdim].size), np.nan)
            for b in bands:
                for s in b['scans']:
                    data[b['index'] - 1, s['index'] - 1] = s[vname]
            variables[vname] = ((bdim, sdim), data)

        # other returned values
        _ = f.readline()
        # TODO: check names
        ckeys = ['fit_rms', 'chi_square_obs', 'dofs_total',
                 'dofs_trg', 'dofs_tpr', 'n_iteration', 'n_iteration_max',
                 'has_converged', 'has_division_warnings']
        _ = f.readline()
        cvals = [s.strip() for s in f.readline().split()]
        attrs = {k: eval_value(v) for k, v in zip(ckeys, cvals)}
        attrs['fit_rms'] /= 100   # convert value
        variables['sfit4_summary'] = ((), np.nan, attrs)

    dataset = {
        'data_vars': variables,
        'coords': coords,
        'attrs': global_attrs
    }

    return dataset
