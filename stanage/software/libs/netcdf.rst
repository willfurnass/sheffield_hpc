.. _netcdf_stanage:

netCDF
======

.. sidebar:: netCDF

   :URL: https://www.unidata.ucar.edu/software/netcdf/

"NetCDF (Network Common Data Form) is a set of interfaces for array-oriented data access and a freely distributed collection of data access libraries for C, Fortran, C++, Java, and other languages. The netCDF libraries support a machine-independent format for representing scientific data. Together, the interfaces, libraries, and format support the creation, access, and sharing of scientific data."

Usage
-----

All of the module files for netCDF will also load:

* the :ref:`gompi or iimpi toolchain<stanage_eb_toolchains>`
* an :ref:`HDF5 library <hdf5_stanage>`
* (and the zlib and Szip libraries)

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

To load **netCDF** run *one* of the following:

.. include:: /referenceinfo/imports/stanage/packages/netcdf-ml-el7-icelake-znver-stanage.rst

To load the **Fortran bindings** for netCDF run *one* of the following:

.. include:: /referenceinfo/imports/stanage/packages/netcdf-fortran-ml-el7-icelake-znver-stanage.rst

To load the **C++ bindings** for netCDF run *one* of the following:

.. include:: /referenceinfo/imports/stanage/packages/netcdf-c++4-ml-el7-icelake-znver-stanage.rst

Testing
-------

To verify the installation of the NetCDF components, you can use the following commands to check the general information and confirm the presence of specific bindings.
For more details, you can refer to the `NetCDF installation guide <https://jules.jchmr.org/check-netcdf>`_.

.. code-block::

   # General NetCDF installation
   nc-config --all # General information on the NetCDF installation

   # NetCDF C++ components
   nc-config --all # General information on the NetCDF installation
   ncxx4-config --all # Check whether NetCDF C++ components have been installed

   # NetCDF Fortran components
   nc-config --all # General information on the NetCDF installation
   nf-config --all # Check whether NetCDF Fortran components have been installed
