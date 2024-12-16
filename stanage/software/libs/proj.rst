.. _proj_stanage:

PROJ
====

.. sidebar:: PROJ

   :Latest Version: 9.1.1
   :URL: https://github.com/OSGeo/proj

PROJ consists of programs and a library for managing cartographic projections.

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

Running one of the following:
        
.. include:: /referenceinfo/imports/stanage/packages/proj-ml-el7-icelake-znver-stanage.rst

will

* add several PROJ programs to your ``PATH`` environment variable
* allow other programs to make use of (dynamically link against) the PROJ library
* activate a version of the GCC compiler (as its C++ standard library is required when using PROJ.X)

You can run ``proj`` to test that you are running the required version ::

    $ proj 
    Rel. 9.1.1, December 1st, 2022
    usage: proj [-bdeEfiIlmorsStTvVwW [args]] [+opt[=arg] ...] [file ...]

Documentation
-------------
Standard ``man`` pages are available for the provided commands/functions ``cs2cs``, ``geod``, ``proj``, ``geodesic`` and ``pj_init``.

These can be viewed using e.g. ::

    $ man proj

Installation notes
------------------

This section is primarily for administrators of the system. PROJ has been installed using the default Easybuild config files.

Build logs and test reports can be found in ``$EBROOTPROJ/easybuild`` with a given module loaded.
