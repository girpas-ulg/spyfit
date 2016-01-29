spyfit
=======

.. image:: https://img.shields.io/travis/benbovy/spyfit.svg
        :target: https://travis-ci.org/benbovy/spyfit
.. image:: https://img.shields.io/pypi/v/spyfit.svg
        :target: https://pypi.python.org/pypi/spyfit

**This package is currently under heavy development!
It hasn't been released yet.**

**spyfit** is a user-friendly interface to the SFIT4_ retrieval software.

It provides command-line utilities, a Python API and a collection of
helpful functions, which all aim to allow easier I/O handling, full automation
and management of SFIT4 runs.

Why spyfit?
------------

There are already Python codes related to SFIT4, e.g., those
written by Eric Nussbaumer (the SFIT4 Layer0 and Layer1), Mathias Palm
and Bavo Langerock.
**spyfit** takes ideas from those codes and has features in common.
Additionaly, it tries to follow as much as possible the good practices
in the Python scientific ecosystem (e.g., readability, reproducibility,
using standard formats, PEP8 style guide / numpy docstrings, packaging...).
It also put emphasis on its integration with the libraries of the Python
scientific ecosystem (including emerging packages like xray_) and using
standard data models such as netCDF_.

Main Features
-------------

- Use the `Common Data Model`_ to store and handle retrieval data.
  The xray_ package - a required dependency - implements this data model and
  provides a powerful framework for easy inspection, merging, processing and
  plotting of retrieval data.
- Support various formats including:
    - netCDF_ (read/write)
    - SFIT4_ ascii output files (read-only) and input files (read/write)
    - GEOMS_ compliant HDF4 format (read/write)
    - Easy export to various formats supported by xray_ and pandas_
      (e.g., hdf5, csv, excel, SQL-databases...)
- More readable, "pythonic" names for parameters and variables (also follow
  the `CF`_ standard names when possible). Possibility to revert to the names
  defined in the SFIT4 core code or in the `sfit4.ctl` file.
- Calculation of retrieval error budgets.
- Full automation of SFIT4 runs. Possibility to extract line list data directly
  from the HITRANonline_ database and to extract pressure and temperature
  profiles directly from online climate data via OpenDAP access (see, e.g.,
  the datasets available at PSD_).
- Highly configurable.

.. _SFIT4: https://wiki.ucar.edu/display/sfit4/Infrared+Working+Group+Retrieval+Code,+SFIT
.. _Common Data Model: http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf
.. _xray: https://github.com/xray/xray
.. _pandas: http://pandas.pydata.org/
.. _CF: http://cfconventions.org/
.. _GEOMS: http://avdc.gsfc.nasa.gov/index.php?site=1178067684
.. _HITRANonline: http://hitran.org/
.. _PSD: http://www.esrl.noaa.gov/psd/data/gridded/

Documentation
-------------

Not yet available.

.. The official documentation is hosted on ReadTheDocs: https://spyfit.readthedocs.org.

Report Issues
-------------

Use the Github issue tracker: https://github.com/girpas-ulg/spyfit/issues

License
-------

Copyright (C) Benoit Bovy, GIRPAS (Ulg) 2015.

Licensed under the GNU General Public License (GPLv3_). See LICENSE.

Some portions of this code have been inspired and/or modified from code
written by Bavo Langerock.

.. _GPLv3: http://www.gnu.org/licenses/gpl-3.0.fr.html
