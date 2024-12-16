.. _gmp_stanage:

GnuMP
=====

.. sidebar:: GMP

    :Latest Version: 6.2.1
    :Dependencies: GCC compiler
    :Documentation: https://gmplib.org/manual/ 


GnuMP is a portable library written in C for arbitrary precision arithmetic on integers, rational numbers, and floating-point numbers. It aims to provide the fastest possible arithmetic for all applications that need higher precision than is directly supported by the basic C types. 


Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

GnuMP can be activated as follows:

.. include:: /referenceinfo/imports/stanage/packages/gmp-ml-el7-icelake-znver-stanage.rst


Installation notes
------------------

GnuMP was installed using Easybuild 4.7.0, build details can be found in ``$EBROOTGMP/easybuild`` with the module loaded.

