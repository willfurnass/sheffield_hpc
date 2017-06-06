.. _gdal_sharc:

Geospatial Data Abstraction Library (GDAL)
==========================================

.. sidebar:: GDAL

   :Latest version: 2.2.0
   :URL: http://www.gdal.org/

GDAL is a library used by many Geographic Information Systems (GIS) packages for converting between many different raster and vector GIS data formats.  It also includes command-line utilities for data translation and processing.  It is released under an an X/MIT style Open Source license by the [Open Source Geospatial Foundation](http://www.osgeo.org/).

Usage
-----

By running ::

    $ module load libs/gdal/2.2.0/gcc-6.2

you

* add several GDAL programs to your ``PATH`` environment variable
* allow other programs to make use of (dynamically link against) the GDAL library
* activate version 4.9.4 of the GCC compiler suite (as its C++ standard library is required by GDAL)

You can run ``gdal-config --version`` to test that you are running the required version ::

    $ gdal-config --version
    2.2.0

The command-line programs that GDAL provides are ::
    $ gdaladdo
    $ gdalbuildvrt
    $ gdal-config
    $ gdal_contour
    $ gdaldem
    $ gdalenhance
    $ gdal_grid
    $ gdalinfo
    $ gdallocationinfo
    $ gdalmanage
    $ gdal_rasterize
    $ gdalserver
    $ gdalsrsinfo
    $ gdaltindex
    $ gdaltransform
    $ gdal_translate
    $ gdalwarp
    $ nearblack
    $ ogr2ogr
    $ ogrinfo
    $ ogrlineref
    $ ogrtindex
    $ testepsg

The ``gdal...`` utilities are mostly for processing raster GIS data formats, whilst the ``ogr...`` utilities are for processing vector GIS data formats.

Documentation
-------------
Standard ``man`` pages are available for the provided commands/functions.

These can be viewed using e.g. ::

    $ man gdal_contour

Much more information is available on the `project site <http://www.gdal.org/>`_.

Supported file formats
----------------------

GDAL has been compiled on this system with support for only a limited set of GIS data formats.  See **Installation notes** below for a list of those provided by each available version of GDAL.

Installation notes
------------------
This section is primarily for administrators of the system.

**Version 2.2.0**

GDAL 2.2.0 was compiled with v4.9.4 of the GCC compiler suite.

#. Download, configure, build and install by switching to a scratch directory and running :download:`this script </sharc/software/install_scripts/libs/gdal/2.2.0/gcc-4.9.4/install.sh>`, ensuring that all stderr and stdout is redirected to :download:`a log file </sharc/software/install_scripts/libs/gdal/2.2.0/gcc-4.9.4/install.log>`. 
#  Test by 
#. Install :download:`this modulefile </sharc/software/modulefiles/libs/gdal/2.2.0/gcc-4.9.4>` as ``/usr/local/modulefiles/libs/gdal/2.2.0/gcc-4.9.4``

The **file formats** supported by this build are listed in the compilation log file.
