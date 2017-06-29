.. _gsl_iceberg:

GNU Scientific Libary (GSL)
===========================

.. sidebar:: GSL

   :Latest version: 2.3
   :URL: https://www.gnu.org/software/gsl

The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License.

The library provides a wide range of mathematical routines such as random number generators, special functions and least-squares fitting. There are over 1000 functions in total with an extensive test suite.

Usage
-----

You can make the GSL available by running: ::

        $ module load libs/gcc/4.9.2/gsl/2.3

This makes the GSL's library files, header files, utility programs and documentation all accessible.

YOu can check that you are running the requested version: ::

        $ gsl-config --version
        2.3

When building software that uses the GSL you will need to *link to it*.  `From the GSL documentation <https://www.gnu.org/software/gsl/manual/html_node/Linking-programs-with-the-library.html#Linking-programs-with-the-library>`_:

    To link against the library you need to specify both the main library and a supporting CBLAS library, which provides standard basic linear algebra subroutines. A suitable CBLAS implementation is provided in the library ``libgslcblas.a`` if your system does not provide one. The following example shows how to link an application with the library: ::

        $ $CC example.o -lgsl -lgslcblas -lm

Another CBLAS implementation that you may want to use with the GSL is the :ref:`Intel Math Kernel Library (MKL) <iceberg_intel_mkl>`.

Documentation
-------------

For documentation either run: ::

        $ module load libs/gcc/4.9.2/gsl/2.3
        $ info gsl

or visit the project website.

Installation notes
------------------
This section is primarily for administrators of the system.

The GSL was installed as a user wanted to build `EIGENSOFT <https://github.com/DReichLab/EIG>`_.

Version 2.3
^^^^^^^^^^^

This was compiled using GCC on an Intel x5650 node.

#. ``cd`` to a scratch directory.
#. Run :download:`this install script (install.sh) </iceberg/software/install_scripts/libs/gcc/4.9.2/gsl/2.3/install.sh>`
#. Install :download:`this module file </iceberg/software/modulefiles/libs/gcc/4.9.2/gsl/2.3>` as ``/usr/local/modulefiles/libs/gcc/4.9.2/gsl/2.3``
