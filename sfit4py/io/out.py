# -*- coding: utf-8 -*-

"""
Functions for reading the SFIT4 output text files into
Python dictionaries and Numpy arrays.

"""

import datetime
import re
import os

import numpy as np


HEADER_PATTERN = (r"\s*SFIT4:V(?P<sfit4_version>[0-9.]+)"
                  r"[\w\s:-]*RUNTIME:(?P<runtime>[0-9:\-]+)"
                  r"\s*(?P<description>.+)")


def parse_header(pattern, line):
    """
    Parse the header line of an output file
    given the regular expression `pattern`.

    """
    m = re.match(pattern, line)
    metadata = m.groupdict()

    if 'runtime' in metadata.keys():
        metadata['runtime'] = datetime.datetime.strptime(
            metadata['runtime'], "%Y%m%d-%H:%M:%S"
        )

    if 'description' in metadata.keys():
        metadata['description'] = metadata['description'].strip().lower()

    return metadata


def read_matrix(filename):
    """
    Read single matrix (or single vector) output files.

    Use this function to load the 'ak.out' and 'seinv.out'
    TODO: list here all files.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(HEADER_PATTERN, f.readline())
        outputd['filename'] = os.path.abspath(filename)

        data_shape = tuple([int(d) for d in f.readline().split() if int(d) > 1])
        data = np.loadtxt(f)
        assert data.shape == data_shape
        outputd['data'] = data

    return outputd


def read_table(filename):
    """
    Read tabular output data.

    Use this function to load the 'k.out', 'kb.out', 'sa.out', 'shat.out'...
    TODO: list here all files.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(HEADER_PATTERN, f.readline())
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

    Use this function to load the 'aprfs.out' and 'rprfs.out' files.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(HEADER_PATTERN, f.readline())
        outputd['filename'] = os.path.abspath(filename)

        meta = f.readline().split()
        nrows = int(meta[1])
        outputd['retrieved_gaz'] = [g.strip() for g in meta[3:]]
        outputd['gaz_index'] = list(map(int, f.readline().split()))
        outputd['column_names'] = [c.strip() for c in f.readline().split()]

        data = np.loadtxt(f)
        assert data.shape == (nrows, len(outputd['column_names']))
        outputd['data'] = data

    return outputd


def read_statevec(filename):
    """
    Read the state vector (a-priori and retrieved profiles / columns / params).

    Use this function to load the 'statevec' file.

    """
    with open(filename, 'r') as f:
        outputd = parse_header(HEADER_PATTERN, f.readline())
        outputd['filename'] = os.path.abspath(filename)

        meta = f.readline().split()
        nlevels, niter, nitermax = list(map(int, meta[:3]))
        is_temp, has_converged, has_divwarn = list(map(
            lambda s: True if s == "T" else False, meta[3:]
        ))
        outputd['n_iteration'] = niter
        outputd['n_iteration_max'] = nitermax
        outputd['is_temp'] = is_temp      # refactor ! what is 'istemp'?
        outputd['has_converged'] = has_converged
        outputd['has_division_warnings'] = has_divwarn

        # coords
        for c in ('altitude', 'pressure', 'temperature'):
            dummy = f.readline()
            outputd[c] = np.fromfile(f, count=nlevels, sep=" ")

        # gaz profiles / columns
        for k in ('apriori_profiles', 'apriori_total_columns',
                  'retrieved_profiles', 'retrieved_total_columns'):
            outputd[k] = dict()

        for i in range(int(f.readline())):
            dummy = f.readline()
            gaz = f.readline().strip()
            outputd['apriori_total_columns'][gaz] = float(f.readline())
            outputd['apriori_profiles'][gaz] = np.fromfile(f, count=nlevels,
                                                           sep=" ")
            dummy = f.readline()
            dummy = f.readline()
            outputd['retrieved_total_columns'][gaz] = float(f.readline())
            outputd['retrieved_profiles'][gaz] = np.fromfile(f, count=nlevels,
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
