.. _openblas_stanage:

OpenBLAS
========

.. sidebar:: OpenBLAS
   
   :URL: https://www.openblas.net/
   :Documentation: https://github.com/xianyi/OpenBLAS/wiki/User-Manual

OpenBLAS is one of the :ref:`BLAS <blas_stanage>` implementations installed on Stanage.
It also provides some optimised LAPACK routines.

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

OpenBLAS can be activated using one of:

.. include:: /referenceinfo/imports/stanage/packages/openblas-ml-el7-icelake-znver-stanage.rst
   
which also loads a version of the :ref:`GCC <gcc_stanage>` compiler.
