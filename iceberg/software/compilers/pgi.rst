.. _`PGI Compilers`:

PGI Compilers
=============
The PGI Compiler suite offers C,C++ and Fortran Compilers. For full details of the features of this compiler suite, see `PGI's website <http://www.pgroup.com/products/pgiworkstation.htm>`_.

Making the PGI Compilers available
----------------------------------

After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qsh` or :code:`qrsh` command. To make one of the versions of the PGI Compiler Suite available, run one of the following module commands ::

    module load compilers/pgi/16.10
    module load compilers/pgi/15.10
    module load compilers/pgi/15.7
    module load compilers/pgi/14.4
    module load compilers/pgi/13.1

Once you've loaded the module, you can check the version with ::

    pgcc --version

Compilation examples
--------------------
**C**

To compile a C hello world example into an executable called ``hello`` using the PGI C compiler ::

    pgcc hello.c -o hello

**C++**

To compile a C++ hello world example into an executable called ``hello`` using the PGI C++ compiler ::

      pgc++ hello.cpp -o hello

**Fortran**

To compile a Fortran hello world example into an executable called ``hello`` using the PGI Fortran compiler ::

      pgf90 hello.f90 -o hello

Compiling on the GPU using the PGI Compiler
-------------------------------------------

**Additional permissions are needed to use GPUs on Iceberg/ShARC. See** :ref:`GPUComputing_iceberg` **for more information.**

Start an interctive GPU session (:ref:`GPUInteractive_iceberg`) and run one of the following module commands, depending on which version of the compilers you wish to use ::

  module load compilers/pgi/16.10
  module load compilers/pgi/15.10
  module load compilers/pgi/15.7
  module load compilers/pgi/14.4
  module load compilers/pgi/13.1


The PGI compilers have several features that make them interesting to users of GPU hardware:-

OpenACC Directives
^^^^^^^^^^^^^^^^^^

OpenACC is a relatively new way of programming GPUs that can be significantly simpler to use than low-level language extensions such as CUDA or OpenCL. From the `OpenACC website <http://www.openacc-standard.org/About_OpenACC>`_ :

    The OpenACC Application Program Interface describes a collection of compiler directives to specify loops and regions of code in standard C, C++ and Fortran to be offloaded from a host CPU to an attached accelerator. OpenACC is designed for portability across operating systems, host CPUs, and a wide range of accelerators, including APUs, GPUs, and many-core coprocessors.

    The directives and programming model defined in the OpenACC API document allow programmers to create high-level host+accelerator programs without the need to explicitly initialize the accelerator, manage data or program transfers between the host and accelerator, or initiate accelerator startup and shutdown.

For more details concerning OpenACC using the PGI compilers, see `The PGI OpenACC website <http://www.pgroup.com/resources/accel.htm>`_.

CUDA Fortran
^^^^^^^^^^^^

In mid 2009, PGI and NVIDIA cooperated to develop CUDA Fortran. CUDA Fortran includes a Fortran 2003 compiler and tool chain for programming NVIDIA GPUs using Fortran.

* `CUDA Fortran Programming Guide <http://www.pgroup.com/lit/whitepapers/pgicudaforug.pdf>`_.

CUDA-x86
^^^^^^^^

NVIDIA CUDA was developed to enable offloading computationally intensive kernels to massively parallel GPUs. Through API function calls and language extensions, CUDA gives developers explicit control over the mapping of general-purpose computational kernels to GPUs, as well as the placement and movement of data between an x86 processor and the GPU.

The PGI CUDA C/C++ compiler for x86 platforms allows developers using CUDA to compile and optimize their CUDA applications to run on x86-based workstations, servers and clusters with or without an NVIDIA GPU accelerator. When run on x86-based systems without a GPU, PGI CUDA C applications use multiple cores and the streaming SIMD (Single Instruction Multiple Data) capabilities of Intel and AMD CPUs for parallel execution.

* `PGI CUDA-x86 guide <http://www.pgroup.com/resources/cuda-x86.htm>`_.

Installation Notes
------------------
*Version 16.10*

The interactive installer was slightly different to that of 15.7 (below) but the questions and answers were essentially the same.

*Version 15.10*

The interactive installer was slightly different to that of 15.7 (below) but the questions and answers were essentially the same.

*Version 15.7*

The installer is interactive. Most of the questions are obvious.
Here is how I answered the rest

Installation type ::

  A network installation will save disk space by having only one copy of the
  compilers and most of the libraries for all compilers on the network, and
  the main installation needs to be done once for all systems on the network.

  1  Single system install
  2  Network install

  Please choose install option: 1

Path ::

  Please specify the directory path under which the software will be installed.
  The default directory is /opt/pgi, but you may install anywhere you wish,
  assuming you have permission to do so.

  Installation directory? [/opt/pgi] /usr/local/packages6/compilers/pgi

CUDA and AMD components ::

  Install CUDA Toolkit Components? (y/n) y
  Install AMD software components? (y/n) y

AMCL version ::

  This PGI version links with ACML 5.3.0 by default.  Also available:
    (1) ACML 5.3.0
    (2) ACML 5.3.0 using FMA4
  Enter another value to override the default (1)
  1

Other questions ::

  Install JAVA JRE [yes] yes
  Install OpenACC Unified Memory Evaluation package? (y/n) n
  Do you wish to update/create links in the 2015 directory? (y/n) y
  Do you wish to install MPICH? (y/n) y
  Do you wish to generate license keys? (y/n) n
  Do you want the files in the install directory to be read-only? (y/n) n

The license file is on the system at ``/usr/local/packages6/compilers/pgi/license.dat`` and is a 5 seat network license. Licenses are only used at compile time.

Extra install steps
-------------------
Unlike gcc, the PGI Compilers do not recognise the environment variable LIBRARY_PATH which is used by a lot of installers to specify the locations of libraries at compile time. This is fixed by creating a ``siterc`` file at ``/usr/local/packages6/compilers/pgi/linux86-64/VER/bin/siterc`` with the following contents ::

  # get the value of the environment variable LIBRARY_PATH
  variable LIBRARY_PATH is environment(LD_LIBRARY_PATH);
  variable inc_path is environment(CPATH);

  # split this value at colons, separate by -L, prepend 1st one by -L
  variable library_path is
  default($if($LIBRARY_PATH,-L$replace($LIBRARY_PATH,":", -L)));

  # add the -L arguments to the link line
  append LDLIBARGS=$library_path;
  append SITEINC=$inc_path;

Where VER is the version number in question: 15.7, 15.10 etc

At the time of writing (August 2015), this is `documented on PGI's website <https://www.pgroup.com/support/link.htm#lib_path_ldflags>`_.

Modulefile
----------
**Version 15.10**
The PGI compiler installer creates a suitable modulefile that's configured to our system. It puts it at ``/usr/local/packages6/compilers/pgi/modulefiles/pgi64/15.10`` so all that is required is to copy this to where we keep modules at ``/usr/local/modulefiles/compilers/pgi/15.10``

**Version 15.7**

The PGI compiler installer creates a suitable modulefile that's configured to our system. It puts it at ``/usr/local/packages6/compilers/pgi/modulefiles/pgi64/15.7`` so all that is required is to copy this to where we keep modules at ``/usr/local/modulefiles/compilers/pgi/15.7``
