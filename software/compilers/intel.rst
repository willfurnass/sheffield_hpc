Intel Compilers
===============
Intel Compilers help create C, C++ and Fortran applications that can take full advantage of the advanced hardware capabilities available in Intel processors and co-processors. They also simplify that development by providing high level parallel models and built-in features like explicit vectorization and optimization reports.

Making the Intel Compilers available
------------------------------------

After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` or :code:`qrsh` command. To make one of the versions of the Intel Compiler available, run one of the following module commands ::

        module load compilers/intel/14.0
        module load compilers/intel/12.1.15
        module load compilers/intel/11.0              

Compilation examples
--------------------
To compile the C++ hello world example into an executable called ``hello`` using the Intel C++ compiler ::

      icc hello.cpp -o hello

Detailed Documentation
----------------------
The following links are to Intel's website

* `User and Reference Guide for the Intel® C++ Compiler 15.0 <https://software.intel.com/en-us/compiler_15.0_ug_c>`_
* `User and Reference Guide for the Intel® Fortran Compiler 15.0 <https://software.intel.com/en-us/compiler_15.0_ug_f>`_

Installation Notes
------------------
The following notes are primarily for system sysadmins.

version 14 and below
~~~~~~~~~~~~~~~~~~~~
Installation notes are not available for these older versions of the Intel Compiler. 