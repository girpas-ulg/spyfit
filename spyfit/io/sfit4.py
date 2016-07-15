# -*- coding: utf-8 -*-

"""Read and write SFIT4 input and output ascii files."""

from __future__ import absolute_import

import os

import xarray as xr

from ._sfit4_read_out import (read_ak_matrix, read_aprfs, read_matrix,
                              read_param_iterations, read_profiles, read_rprfs,
                              read_seinv_vector, read_single_spectra,
                              read_single_spectrum, read_solar_spectrum,
                              read_spectra, read_state_vector, read_summary,
                              read_table)
from ._sfit4_read_in import (read_ctl, read_layers, read_ref_profiles,
                             read_spectrum)
from ._sfit4_write_in import write_reference


# public API
read_single_spectrum = read_single_spectrum
read_profiles = read_profiles
read_matrix = read_matrix
write_reference = write_reference


def _read_g_out(fname):
    return read_table(fname, var_name="gain_matrix",
                      dims=('statevector', 'diag'), index_cols=False)


def _read_k_out(fname):
    return read_table(fname, var_name="k_matrix", dims=('diag', 'statevector'))


def _read_kb_out(fname):
    return read_table(fname, var_name='kb_vectors', dims=('diag', 'mparam'))


def _read_shat_out(fname):
    return read_table(fname, var_name="cov_matrix_fitted",
                      dims=('statevector', 'statevector'))


def _read_sa_out(fname):
    return read_table(fname, var_name="cov_matrix_initial",
                      dims=('statevector', 'statevector'))


def _read_sainv_out(fname):
    return read_table(fname, var_name="cov_matrix_initial_inv",
                      dims=('statevector', 'statevector'))


def _read_summary_clean(fname):
    d = read_summary(fname)
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


def _read_spectrum_raw(fname):
    # spectral, band and scan data in `file.in.spectrum` don't
    # seem to be always equal to observed spectrum in `file.out.pbpfile`.
    # use different dimension and variable names here.
    d = read_spectrum(fname, spdim='spectrum_in', bdim='band_in',
                      sdim='scan_in', wcoord='spec_in_wn',
                      scoord='spec_in_scan', bcoord='spec_in_band')
    d['data_vars'] = {k.replace('spec_', 'spec_in_'): v
                      for k, v in d['data_vars'].items()}
    return d


_MAP_CTL_READ_FUNC = {
    # sfit4.ctl entry: (default sfit4 filename, read function)
    'in__stalayers': ('station.layers', read_layers),
    'in__refprofile': ('reference.prf', read_ref_profiles),
    'in__spectrum': ('t15asc.4', _read_spectrum_raw),
    'out__ak_matrix': ('ak.out', read_ak_matrix),
    'out__seinv_vector': ('seinv.out', read_seinv_vector),
    'out__g_matrix': ('g.out', _read_g_out),
    'out__k_matrix': ('k.out', _read_k_out),
    'out__kb_matrix': ('kb.out', _read_kb_out),
    'out__shat_matrix': ('shat.complete', _read_shat_out),
    'out__sa_matrix': ('sa.complete', _read_sa_out),
    'out__sainv_matrix': ('sainv.complete', _read_sainv_out),
    'out__aprprofiles': ('aprfs.table', read_aprfs),
    'out__retprofiles': ('rprfs.table', read_rprfs),
    'out__pbpfile': ('pbpfile', read_spectra),
    'out__statevec': ('statevec', read_state_vector),
    'out__summary': ('summary', _read_summary_clean),
    'out__parm_vectors': ('parm.vectors', read_param_iterations),
    'out__solarspectrum': ('solspec.dat', read_solar_spectrum),
    # TODO: raytrace.mix, chnspec1-2
}


def _read_and_merge(dataset, read_func, filename, **kwargs):
    temp_ds = xr.Dataset(**read_func(filename))
    dataset.merge(temp_ds, inplace=True, compat='equals', **kwargs)

    # xarray doesn't handle attribute merging
    # update attrs dicts without checking (overwrite by default)
    for k, v in temp_ds.variables.items():
        dataset[k].attrs.update(v.attrs)
    dataset.attrs.update(temp_ds.attrs)


def load_sfit4(ctl_filename):
    """Load data from a SFIT4 run.

    It loads data from SFIT4 input and output files
    enabled / specified in the given .ctl file.

    Parameters
    ----------
    ctl_filename : str
        name or path to the sfit4.ctl file.

    Returns
    -------
    `xarray.Dataset`

    """
    ctl_path = os.path.expanduser(os.path.realpath(ctl_filename))
    ctl_dir = os.path.dirname(ctl_path)
    dataset = xr.Dataset(**read_ctl(ctl_filename))
    sfit4_inputs = dataset.sfit4_ctl.attrs

    for input_name, v in _MAP_CTL_READ_FUNC.items():
        default_filename, read_func = v
        if input_name.startswith('out'):
            if not sfit4_inputs.get(input_name, False):
                continue
        filename = sfit4_inputs.get('file__{}'.format(input_name),
                                    default_filename)
        if not os.path.isabs(filename):
            filename = os.path.join(ctl_dir, filename)
        if os.path.exists(filename):
            _read_and_merge(dataset, read_func, filename)

    if sfit4_inputs.get('out__gas_spectra', False):
        # overwrite spec_fitted__ALL (not the same precision
        # than in pbpfile)
        overwrite_spec = 'spec_fitted__ALL'
        _read_and_merge(dataset, read_single_spectra,
                        os.path.join(ctl_dir, 'spc.*'),
                        overwrite_vars=overwrite_spec)

    # TODO: convert spectrum coord to multi-index
    # TODO: check 'diag' dim, same scan/band ordering than 'spectrum' dim?

    # update/overwrite global attributes
    dataset.attrs['source'] = ctl_path
    dataset.attrs['description'] = 'data from a single sfit4 run'

    return dataset
