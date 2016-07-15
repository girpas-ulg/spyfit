# -*- coding: utf-8 -*-

"""
Write SFIT4 input ascii files.
"""

from collections import OrderedDict

import numpy as np


def _format_line(values, fmt):
    "fmt is either a string or a list of string of `values` length"
    if isinstance(fmt, str):
        fmt = [fmt] * len(values)
    return ''.join([f.format(v) for f, v in zip(fmt, values)]) + '\n'


def _format_array(arr, fmt='.18e', delimiter=' ', pad=0, left_pad=0,
                  newline='\n', max_line_items=5):
    sep = delimiter + " " * pad
    n_splits = np.ceil(arr.size / max_line_items)
    lines = []
    for sub_arr in np.array_split(arr, n_splits):
        str_items = sep.join(['{:{}}'.format(i, fmt) for i in sub_arr])
        lines.append(" " * left_pad + str_items)
    return newline.join(lines) + '\n'


def write_reference(dataset, filename, map_gas_id_name,
                    var_names='reference__{}', ordering='ascending'):
    """
    Export reference ZPT and gas profiles to SFIT4 ascii format.

    Use this function to create 'reference.prf' file.

    Parameters
    ----------
    dataset : dict or xr.Dataset
        Either a xarray dataset or a mapping that contains profile data
        for ZPT and gas profiles.
        Values of the mapping may be array-like objects or dictionaries
        that must at least contain a 'data' item.
    filename : str
        Name or path to the file to write.
    map_gas_id_name : dict
        A mapping of SFIT4 gas id (keys) and gas name or 'OTHER' (values).
    var_names : str
        Template for variable names (i.e., keys of `dataset`). Must at least
        contain one parameter '{}' that will be replaced by 'altitude',
        'pressure', 'temperature' or by the gas name.
    ordering: {'ascending', 'descending'}
        Ordering of the atmospheric profile values ('ascending': from bottom
        to top, 'descending': from top to bottom).

    Notes
    -----
    This function also looks for the 'description' and 'ordering' attributes
    or items in each variable/item of `dataset` to write profile headers
    and to set the correct data ordering (ignored if these attributes/items
    don't exist).

    """
    sfit4_ordering = {'ascending': 0, 'descending': 1}

    map_id_name = OrderedDict(
        ((-3, 'altitude'), (-2, 'pressure'), (-1, 'temperature'))
    )
    map_id_name.update(sorted(map_gas_id_name.items()))

    profiles = []
    for id, name in map_id_name.items():
        if name == 'OTHER':
            profiles.append([id, name, '', None])
        else:
            variable = dataset[var_names.format(name)]
            if hasattr(variable, 'attrs'):
                # case of xarray dataset
                description = variable.attrs.get('description', '')
                data = variable.values
                order = variable.attrs.get('ordering')
            else:
                description = variable.get('description', '')
                data = variable['data']
                order = variable.attrs.get('ordering')
            if data.ndim != 1:
                raise ValueError("profile '{}' must be 1-dimensional, found "
                                 "{} dimensions".format(name, data.ndim))
            if order is not None:
                if order not in ('ascending', 'descending', 0, 1):
                    raise ValueError("unrecognized ordering value: {} for "
                                     "profile {}".format(order, name))
                if order != ordering and order != sfit4_ordering[ordering]:
                    data = data[-1::-1]
            profiles.append([id, name, description, data])

    valid_profiles = [data for _, name, _, data in profiles
                      if name != 'OTHER']
    sizes = set([p.size for p in valid_profiles])
    if len(sizes) > 1:
        raise ValueError("all profiles must have the same size, found {}"
                         .format(sizes))
    n_levels = sizes.pop()
    for p in profiles:
        if p[1] == 'OTHER':
            p[3] = np.zeros((n_levels), dtype='f')

    with open(filename, 'w') as ref_file:
        ref_file.write(_format_line(
            (sfit4_ordering[ordering], n_levels, len(map_gas_id_name)),
            '{:>12}'
        ))
        for id, name, description, data in profiles:
            if id < 0:
                ref_file.write(' {}\n'.format(description))
            else:
                ref_file.write(_format_line((id, name, description),
                                            ['{:>5}', '{:>8}', '{:>75}']))
            ref_file.write(_format_array(data, fmt='.4E', pad=2, left_pad=2))
