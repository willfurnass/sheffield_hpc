.. _hdf5_stanage:

HDF5
====

.. sidebar:: HDF5

   :Latest Version: 1.14.0
   :URL: https://www.hdfgroup.org/solutions/hdf5/

"HDF5 is a data model, library, and file format for storing and managing data. It supports an unlimited variety of datatypes, and is designed for flexible and efficient I/O and for high volume and complex data. HDF5 is portable and is extensible, allowing applications to evolve in their use of HDF5. The HDF5 Technology suite includes tools and applications for managing, manipulating, viewing, and analyzing data in the HDF5 format."

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

To load this library plus
the :ref:`gompi or iimpi toolchain<stanage_eb_toolchains>`
(and the zlib and Szip libraries)
run *one* of the following:
   
.. include:: /referenceinfo/imports/stanage/packages/hdf5-ml-el7-icelake-znver-stanage.rst

Installation notes
------------------

This section is primarily for administrators of the system. HDF5 has been installed using the default Easybuild config files.

Build logs and test reports can be found in ``$EBROOTHDF5/easybuild`` with a given module loaded.

