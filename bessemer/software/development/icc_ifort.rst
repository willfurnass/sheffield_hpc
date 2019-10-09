.. _icc_ifort_bessemer:

Intel Compilers
===============

Intel has compilers for C (icc), C++ (icpc) and Fortran (ifort).

To activate icc/icpc, 
first connect to Bessemer, 
start an interactive session: ::

   srun --pty bash -i

then run: :: 

   module load icc/2019.1.144-GCC-8.2.0-2.31.1

which implicitly loads :ref:`GCC <gcc_bessemer>` 8.2.0, which is a dependency.

Confirm that you've loaded the version of icc/icpc you require: ::

   icc -v

To instead activate ifort you need the following ``module load`` command: ::

   module load ifort/2019.1.144-GCC-8.2.0-2.31.1

which also implicitly loads :ref:`GCC <gcc_bessemer>` 8.2.0.

To load *both* icc/icpc *and* ifort: ::

   module load iccifort/2019.1.144-GCC-8.2.0-2.31.1

Documentation
-------------

man pages are available on the system.
Once you have loaded the required version of icc/icpc/ifort, run: ::

   man icc
   man ifort
