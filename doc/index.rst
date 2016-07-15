Spyfit
======

**spyfit** is a Python package that provides a set of tools for
easy handling of FTIR retrieval data and for flexible setup and
execution of retrieval processing pipelines.

Our goal with spyfit is to provide deep integration with the libraries
of the Python scientific stack and to promote the use of standard data models.
Through the use of the `xarray`_ Python package, we adopt the
`Common Data Model`_ for self-describing retrieval data (it can be viewed as
an in-memory representation of a `netCDF`_ file).

**NOTE:** This package is currently under heavy development! API is not stable.

.. _xarray: http://xarray.pydata.org
.. _Common Data Model: http://www.unidata.ucar.edu/software/thredds/current/netcdf-java/CDM
.. _netCDF: http://www.unidata.ucar.edu/software/netcdf

Documentation
-------------

.. toctree::
   :maxdepth: 1

   installing
   api
   contributing

Get in touch
------------

- Report bugs, suggest feature ideas or view the source code `on GitHub`_.

.. _on GitHub: https://github.com/girpas-ulg/spyfit

License
-------

spyfit is available under the open-source GNU General Public License (`GPLv3`_).

.. _GPLv3: http://www.gnu.org/licenses/gpl-3.0.fr.html

About
-----

spyfit is part of various Python tools developped for FTIR spectroscopy and
atmospheric chemistry modelling at the Infrared Group of Atmospheric and
Solar Physics - `GIRPAS`_ laboratory (University of Liege, Belgium).

.. _GIRPAS: http://labos.ulg.ac.be/girpas/en/
