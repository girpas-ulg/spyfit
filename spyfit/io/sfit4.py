# -*- coding: utf-8 -*-

"""Read and write SFIT4 input and output ascii files."""

from __future__ import absolute_import

import os

import xarray as xr

from ._sfit4out import *
from ._sfit4in import *


_read_g_out = lambda f: read_table(f, var_name="gain_matrix",
                                   dims=('statevector', 'diag'),
                                   index_cols=False)
_read_k_out = lambda f: read_table(f, var_name="k_matrix",
                                   dims=('diag', 'statevector'))
_read_kb_out = lambda f: read_table(f, var_name='kb_vectors',
                                    dims=('diag', 'mparam'))
_read_shat_out = lambda f: read_table(f, var_name="cov_matrix_fitted",
                                      dims=('statevector', 'statevector'))
_read_sa_out = lambda f: read_table(f, var_name="cov_matrix_initial",
                                    dims=('statevector', 'statevector'))
_read_sainv_out = lambda f: read_table(f, var_name="cov_matrix_initial_inv",
                                       dims=('statevector', 'statevector'))


def _read_summary_clean(f):
    d = read_summary(f)
    # total column values may be slightly different than in
    # `file.out.retprofiles` (not the same precision): skip it here.
    for k in list(d['data_vars'].keys()):
        if 'total_column' in k:
            del d['data_vars'][k]
    # wavenumber start/stop values are the same than in `sfit4.ctl`
    # but not than in `file.out.pbpfile`: skip it here.
    for k in ('spec_wn', 'spec_band', 'spec_scan'):
        del d['coords'][k]
    return d


def _read_spectrum_raw(f):
    # spectral, band and scan data in `file.in.spectrum` don't
    # seem to be always equal to observed spectrum in `file.out.pbpfile`.
    # use different dimension and variable names here.
    d = read_spectrum(f, spdim='spectrum_in', bdim='band_in',
                      sdim='scan_in', wcoord='spec_in_wn',
                      scoord='spec_in_scan', bcoord='spec_in_band')
    d['data_vars'] = {k.replace('spec_', 'spec_in_'): v
                      for k, v in d['data_vars'].items()}
    return d


_map_file_read_func = {
    # sfit4.ctl field name: (default sfit4 filename, read function)
    'file__in__stalayers': ('station.layers', read_layers),
    'file__in__refprofile': ('reference.prf', read_ref_profiles),
    'file__in__spectrum': ('t15asc.4', _read_spectrum_raw),
    'file__out__ak_matrix': ('ak.out', read_ak_matrix),
    'file__out__seinv_vector': ('seinv.out', read_seinv_vector),
    'file__out__g_matrix': ('g.out', _read_g_out),
    'file__out__k_matrix': ('k.out', _read_k_out),
    'file__out__kb_matrix': ('kb.out', _read_kb_out),
    'file__out__shat_matrix': ('shat.complete', _read_shat_out),
    'file__out__sa_matrix': ('sa.complete', _read_sa_out),
    'file__out__sainv_matrix': ('sainv.complete', _read_sainv_out),
    'file__out__aprprofiles': ('aprfs.table', read_aprfs),
    'file__out__retprofiles': ('rprfs.table', read_rprfs),
    'file__out__pbpfile': ('pbpfile', read_spectra),
    'file__out__statevec': ('statevec', read_state_vector),
    'file__out__summary': ('summary', _read_summary_clean),
    'file__out__parm_vectors': ('parm.vectors', read_param_iterations),
    'file__out__solarspectrum': ('solspec.dat', read_solar_spectrum),
    # TODO: raytrace.mix, chnspec1-2
}


def _read_and_merge(ds, read_func, filename, **kwargs):
    temp_ds = xr.Dataset(**read_func(filename))
    ds.merge(temp_ds, inplace=True, compat='equals', **kwargs)

    # xarray doesn't handle attribute merging
    # update attrs dicts without checking
    for k, v in temp_ds.variables.items():
        ds[k].attrs.update(v.attrs)
    ds.attrs.update(temp_ds.attrs)


def load_sfit4_rundir(dirname):
    """Load data from a single SFIT4 run directory.

    It will load data from all supported SFIT4 input and output files
    found in the directory (non recursive).

    Parameters
    ----------
    dirname : str
        name or path to the directory.

    Returns
    -------
    `xarray.Dataset`

    """
    ds = xr.Dataset(**read_ctl(os.path.join(dirname, "sfit4.ctl")))
    sfit4_inputs = ds.sfit4_ctl.attrs

    for input_name, v in _map_file_read_func.items():
        default_filename, read_func = v
        filename = sfit4_inputs.get(input_name, default_filename)
        if os.path.exists(filename):
            _read_and_merge(ds, read_func, filename)

    if sfit4_inputs.get('out__gas_spectra', False):
        # overwrite spec_fitted__ALL (not the same precision
        # than in pbpfile)
        overwrite_spec = 'spec_fitted__ALL'
        _read_and_merge(ds, read_single_spectra,
                        os.path.join(dirname, 'spc.*'),
                        overwrite_vars=overwrite_spec)

    # TODO: convert spectrum coord to multi-index
    # TODO: check 'diag' dimension, same scan/band ordering than 'spectrum' dim?

    # update/overwrite global attributes
    ds.attrs['source'] = os.path.abspath(dirname)
    ds.attrs['description'] = 'data from a single sfit4 run'

    return ds
