# -*- coding: utf-8 -*-

"""
Functions for loading/saving/concatenating retrievals using the
`xray` package (https://github.com/xray/xray).

"""

import xray


def open_retrieval(path):
    """
    Load retrieval data as a :class:`xray.Dataset`.

    Parameters
    ----------
    path : str
        Either path to a netCDF file of previously exported
        (or concatenated) retrieval data or path to a directory
        containing all the input/output files of a single SFIT4 run.

    Returns
    -------
    dataset : :class:`xray.Dataset` object
        The loaded dataset.

    """
    pass   # Not yet implemented


def open_mfretrieval(paths):
    """
    Load and concatenate retrieval data from multiple SFIT4 runs.

    The data of each run must have been already exported to netCDF files.
    The data is concatenated over a new 'time' dimension, so each single run
    must be compatible (i.e., using the same bands, retrieved gases and
    parameters, etc) and can't overlap in time..

    Parameters
    ----------
    paths : str or sequence
        Either a string glob in the form “path/to/my/files/*.nc” or an
        explicit list of netCDF files to open.

    Returns
    -------
    dataset : :class:`xray.Dataset` object
        The loaded and concatenated dataset.

    """
    pass   # Not yet implemented
