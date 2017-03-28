.. _`PGI Compilers_sharc`:

PGI Compilers
=============
The PGI Compiler suite offers C,C++ and Fortran Compilers. For full details of the features of this compiler suite, see `PGI's website <http://www.pgroup.com/products/pgiworkstation.htm>`_.

Making the PGI Compilers available
----------------------------------

After connecting to the ShARC cluster (see :ref:`ssh`),  start an interactive session with the :ref:`qrshx` or :ref:`qrsh` command then activate a specific version of the compiler suite using one of: ::

    module load dev/PGI-compilers/16.10

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

Installation Notes
------------------
*Version 16.10*

The installer is interactive. Here is a log of the questions and answers. ::

  A network installation will save disk space by having only one copy of the
  compilers and most of the libraries for all compilers on the network, and
  the main installation needs to be done once for all systems on the network.

  1  Single system install
  2  Network install

  Please choose install option: 1

  Please specify the directory path under which the software will be installed.
  The default directory is /opt/pgi, but you may install anywhere you wish,
  assuming you have permission to do so.

  Installation directory? [/opt/pgi] /usr/local/packages/dev/pgi

  If you use the 2016 directory in your path, you may choose to
  update the links in that directory to point to the 16.10 directory.

  Do you wish to update/create links in the 2016 directory? (y/n) y
  Making symbolic links in /usr/local/packages/dev/pgi/linux86-64/2016

  Installing PGI JAVA components into /usr/local/packages/dev/pgi
  Installing PGI CUDA components into /usr/local/packages/dev/pgi
  Installing AMD GPU components into /usr/local/packages/dev/pgi
  Installing PGI OpenACC Unified Memory components into /usr/local/packages/dev/pgi ...

  ************************************************************************
  MPI
  ************************************************************************
  This release contains version 1.10.2 of the Open MPI library.

  Press enter to continue...

  Do you want to install Open MPI onto your system? (y/n) y
  Do you want to enable NVIDIA GPU support in Open MPI? (y/n) y

  Do you wish to generate license keys or configure license service? (y/n) n
  The PGI license management script is available at:
  /usr/local/packages/dev/pgi/linux86-64/16.10/bin/pgi_license_tool

  Do you want the files in the install directory to be read-only? (y/n) n

Modulefile
----------
**Version 16.10**
The PGI compiler installer creates a suitable modulefile that's configured to our system. It puts it at
``/usr/local/packages/dev/pgi/modulefiles/pgi64/16.10`` so all that is required is to copy this to where we keep modules at ``/usr/local/modulefiles/dev/PGI-compilers/16.10``
