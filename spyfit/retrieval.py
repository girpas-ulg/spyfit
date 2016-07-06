# -*- coding: utf-8 -*-

"""Load/Save retrieval data from/into `xarray` datasets.

It allows to easily handle, concatenate and plot retrieval data from
different sources (e.g., SFIT4 IO ascii files, netCDF or HDF4 GEOMS files)
using the data model and high-level API functions implemented in `xarray`
(https://github.com/pydata/xarray).

"""

import xarray as xr

from .io.sfit4 import load_sfit4


def load_dataset(filename_or_path, fmt='netcdf'):
    """
    Load a retrieval dataset stored in a given format.

    Parameters
    ----------
    filename_or_path : str or sequence
        name or path to the file or directory. Multiple file name(s) or
        directorie(s) can be given (not yet implemented).
    fmt : {'netcdf', 'sfit4', 'geoms'}
        supported formats are netcdf (default) and sfit4 input/output
        ascii files. GEOMS format is not yet implemented.

    Returns
    -------
    `xarray.Dataset`
        if multiple filenames or paths are given, it will try to
        merge the data into a single `xarray.Dataset`.

    """
    if fmt == 'netcdf':
        return xr.open_dataset(filename_or_path)
    elif fmt == 'sfit4':
        return load_sfit4(filename_or_path)
    else:
        raise ValueError("Unrecognized format '{}'".format(fmt))
