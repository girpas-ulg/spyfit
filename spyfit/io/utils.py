# -*- coding: utf-8 -*-

"""
Various utility functions for spyfit io internals.
"""

import os
import re
from ast import literal_eval

import numpy as np


def expand_path(filename):
    """Get full path."""
    return os.path.expanduser(os.path.realpath(filename))


def get_next_line(ofile, comment='#'):
    "Read next line, skip commented or empty lines."
    while True:
        line = ofile.readline()
        if not line:
            return None
        line = line.strip()
        if line and not line.startswith('#'):
            return line


def get_1d_array(ofile, count):
    """Get 1-d array where separators can be space or comma or both."""
    current_position = ofile.tell()
    array = np.fromfile(ofile, count=count, sep=" ")
    if array.size != count:
        ofile.seek(current_position)
        array = np.fromfile(ofile, count=count, sep=",")
        current_position = ofile.tell()
        if ofile.read(1) != '\n':
            # no trailing comma
            ofile.seek(current_position)
    return array


def sanatize_var_name(var_name):
    """Sanatize `var_name` so that it works with autocompletion."""
    return var_name.lower().strip().replace('.', '__')


def eval_value(string):
    """Eval value (number, bool 'T'/'F', Fortran double)."""
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


def format_line(values, fmt):
    """fmt is either a string or a list of string of `values` length."""
    if isinstance(fmt, str):
        fmt = [fmt] * len(values)
    return ''.join([f.format(v) for f, v in zip(fmt, values)]) + '\n'


def format_array(arr, fmt='.18e', delimiter=' ', pad=0, left_pad=0,
                 newline='\n', max_line_items=5):
    """Convert array to a string."""
    sep = delimiter + " " * pad
    n_splits = np.ceil(arr.size / max_line_items)
    lines = []
    for sub_arr in np.array_split(arr, n_splits):
        str_items = sep.join(['{:{}}'.format(i, fmt) for i in sub_arr])
        lines.append(" " * left_pad + str_items)
    return newline.join(lines) + '\n'
