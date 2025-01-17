.. _udunits_stanage:

udunits
=======

.. sidebar:: udunits

    :Versions: 2.2.28, 2.2.26
    :URL: https://www.unidata.ucar.edu/software/udunits


The UDUNITS package supports units of physical quantities. 
Its C library provides for arithmetic manipulation of units and for conversion 
of numeric values between compatible units. The package contains an extensive unit database, 
which is in XML format and user-extendable. The package also contains a command-line utility 
for investigating units and converting values.

.. caution::

        UDUNITS is typically loaded as an external dependency for R. Please ensure you select the matching 
        GCC compiler versions of your version of R and the UDUNITS libraries.

--------

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

To make this library available, run one of the following: 

.. include:: /referenceinfo/imports/stanage/packages/udunits-ml-el7-icelake-znver-stanage.rst

Installation notes
------------------
This section is primarily for administrators of the system. 

udunits was installed using Easybuild 4.7.0, build details can be found in ``$EBROOTGMP/easybuild`` with the module loaded.

------------------

Testing
-------

1. Load module.
2. Run “udunits2“.
3. For this test we convert 5km into miles, which produces the following results:

.. code-block::

    You have: 5km
    You want: miles
    5 km = 3.10686 miles
    x/miles = 0.621371*(x/km)

.. include:: /referenceinfo/imports/stanage/packages/udunits-dpnd-el7-icelake-znver-stanage.rst
