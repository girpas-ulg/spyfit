.. _installing:

Installation
============

Required dependencies
---------------------

- Python 3.4 or 3.5.
- `numpy <http://www.numpy.org/>`__ (1.7 or later)
- `pandas <http://pandas.pydata.org/>`__ (0.15.0 or later)
- `xarray <http://dask.pydata.org>`__ (0.7.0 or later)

Optional dependencies
---------------------

For netCDF and IO
~~~~~~~~~~~~~~~~~

- `netCDF4 <https://github.com/Unidata/netcdf4-python>`__: used by xarray
  to read/write netCDF4 files
- `scipy <http://scipy.org/>`__: used by xarray for reading/writing netCDF3
- `h5netcdf <https://github.com/shoyer/h5netcdf>`__: an alternative library for
  reading and writing netCDF4 files that does not use the netCDF-C libraries
- `HDF4 <https://www.hdfgroup.org/products/hdf4/>`__ and
  `pyhdf <http://pysclint.sourceforge.net/pyhdf/>`__ or
  `python-hdf4 <http://fhs.github.io/python-hdf4/>`__: if you want to read/write
  HDF4 GEOMS files

For plotting
~~~~~~~~~~~~

- `matplotlib <http://matplotlib.org/>`__
- `cartopy <http://scitools.org.uk/cartopy/>`__

Instructions
------------

There is currently no spyfit release available at PyPi_ or Anaconda_.

Be sure you have the required dependencies (numpy, pandas and xarray)
installed first. You might consider using conda_ to install them::

    $ conda install xarray netcdf4 pip

Then you can clone the spyfit git repository and install it using pip::

    $ git clone https://github.com/girpas-ulg/spyfit
    $ cd spyfit
    $ pip install .

For development purpose, use the following command::

    $ pip install -e .

.. _PyPi: https://pypi.python.org/pypi
.. _Anaconda: https://docs.continuum.io/anaconda/index
.. _conda: http://conda.io/
