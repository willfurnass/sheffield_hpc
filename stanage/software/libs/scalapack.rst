.. _scalapack_stanage:

ScaLAPACK
=========

.. sidebar:: ScaLAPACK
   
   :URL: http://www.netlib.org/scalapack/
   :Documentation: http://www.netlib.org/scalapack/#_documentation

ScaLAPACK is a library of high-performance linear algebra routines
for parallel distributed memory machines.
ScaLAPACK solves:

* Dense and banded linear systems
* Least squares problems
* Eigenvalue problems
* Singular value problems

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

ScaLAPACK can be activated in several ways.

To load ScaLAPACK plus

* a version of :ref:`OpenBLAS <openblas_stanage>`,
* and the :ref:`gompi or gompic toolchain <stanage_eb_toolchains>`

ScaLAPACK can be activated using one of:

.. include:: /referenceinfo/imports/stanage/packages/scalapack-ml-el7-icelake-znver-stanage.rst
   
Note that all load OpenBLAS, despite the change in the module naming convention for more recent toolchains. :ref:`See matching foss toolchains<foss-toolchain-table>`

