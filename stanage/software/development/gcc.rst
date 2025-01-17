.. _gcc_stanage:

GNU Compiler Collection (gcc)
=============================

The GNU Compiler Collection (gcc) is a widely used, free collection of compilers
for C (gcc), C++ (g++) and Fortran (gfortran).

It is possible to switch versions of the gcc compiler suite using modules.

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

then choose the version of the compiler you wish to use
by running *one* of the following lines:

.. include:: /referenceinfo/imports/stanage/packages/gcc-ml-el7-icelake-znver-stanage.rst

Confirm that you've loaded the version of gcc you wanted using ``gcc -v``. :ref:`See matching foss toolchains<foss-toolchain-table>`

Language support
----------------

* Which version(s) of GCC support which features of the `C99 <https://gcc.gnu.org/c99status.html>`__ and `C11 <https://gcc.gnu.org/wiki/C11Status>`__ standards?
* `Which version(s) of GCC support which features of the C++98, C++11, C++14 and C++17 standards? <https://gcc.gnu.org/projects/cxx-status.html>`__

Documentation
-------------

man pages are available on the system.
Once you have loaded the required version of ``gcc``, type ::

    man gcc

* `What's new in the gcc version 12 series? <https://gcc.gnu.org/gcc-12/changes.html>`__

.. include:: /referenceinfo/imports/stanage/packages/gcc-dpnd-el7-icelake-znver-stanage.rst
