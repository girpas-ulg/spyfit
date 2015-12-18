sfit4py
=======

.. image:: https://img.shields.io/travis/benbovy/sfit4py.svg
        :target: https://travis-ci.org/benbovy/sfit4py
.. image:: https://img.shields.io/pypi/v/sfit4py.svg
        :target: https://pypi.python.org/pypi/sfit4py

**This package is currently under heavy development!
It hasn't been released yet.**

**sfit4py** is a user-friendly interface to the SFIT4_ retrieval software.

It provides command-line utilities, a Python API and a collection of
helpful functions, which all aim to allow easier I/O handling, full automation
and management of SFIT4 runs.

.. _SFIT4: https://wiki.ucar.edu/display/sfit4/Infrared+Working+Group+Retrieval+Code,+SFIT

Why sfit4py?
------------

There are already Python codes related to SFIT4, e.g., those
written by Eric Nussbaumer (the SFIT4 Layer0 and Layer1), Mathias Palm
and Bavo Langerock.
**sfit4py** takes ideas from those codes and has features in common.
Additionaly, it tries to follow as much as possible the good practices
in the Python scientific ecosystem (e.g., readability, reproducibility,
using standard formats, PEP8 style guide / numpy docstrings, packaging...).
It also put emphasis on its integration with the libraries of the Python
scientific ecosystem (including emerging packages like xray_) and using
standard data models such as netCDF_.

Main Features
-------------

- Use the `Common Data Model`_ to load, export and handle retrieval data.
  The xray_ package - a required dependency - implements this data model and
  provides a powerful framework for easy inspection, merging and processing of
  retrieval data.
- Support for the SFIT4 ascii formats (reading output files, reading/writing
  input files), the GEOMS_ compliant HDF4 format (reading/writing), the
  netCDF_ format (reading/writing).
- Easy data extraction and export to various formats directly or indirectly
  supported by xray_ (e.g., hdf5, csv, excel, SQL...).
- More readable, "pythonic" names for parameters and variables (also follow
  the `CF`_ standard names when possible). Possibility to revert to the names
  defined in the SFIT4 core code or in the `sfit4.ctl` file.
- Calculation of retrieval error budgets.
- Full-automation of SFIT4 runs. Archiving and management using a
  sqlite database.
- Highly configurable and modular.

.. _Common Data Model: http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _xray: https://github.com/xray/xray
.. _CF: http://cfconventions.org/
.. _GEOMS: http://avdc.gsfc.nasa.gov/index.php?site=1178067684

Documentation
-------------

Not yet available.

.. The official documentation is hosted on ReadTheDocs: https://sfit4py.readthedocs.org.

Report Issues
-------------

Use the Github issue tracker: https://github.com/girpas-ulg/sfit4py/issues

License
-------

Copyright (C) Benoit Bovy, GIRPAS (Ulg) 2015.

Licensed under the GNU General Public License (GPLv3_). See LICENSE.

Some portions of this code have been inspired and/or modified from code
written by Bavo Langerock.

.. _GPLv3: http://www.gnu.org/licenses/gpl-3.0.fr.html
