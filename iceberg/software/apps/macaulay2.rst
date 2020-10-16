.. include:: warning.rst 

.. _macaulay2:

Macaulay2
=========

.. sidebar:: Macaulay2

   :Latest version: 1.12.0.1
   :Dependancies: compilers/gcc/5.3 libs/gcc/5.2/boost/1.59 libs/gcc/lapack/3.3.0 apps/gcc/4.4.7/xzutils/5.2.4
   :URL: http://www.math.uiuc.edu/Macaulay2/

*Macaulay2* is a software system devoted to supporting research in
`algebraic geometry <http://en.wikipedia.org/wiki/Algebraic_geometry>`__
and `commutative algebra <http://en.wikipedia.org/wiki/Commutative_algebra>`__.

*Macaulay2* includes core algorithms for computing
`Gröbner bases <http://www.math.uiuc.edu/Macaulay2/Documentation/___Gröbner_spbases.html>`__
and graded or multi-graded free `resolutions <http://www.math.uiuc.edu/Macaulay2/Documentation/_resolution_lp__Module_rp.html>`__
of `modules <http://www.math.uiuc.edu/Macaulay2/Documentation/_modules.html>`__
over `quotient rings <http://www.math.uiuc.edu/Macaulay2/Documentation/>`__ of
`graded or multi-graded <http://www.math.uiuc.edu/Macaulay2/Documentation/_graded_spand_spmultigraded_sppolynomial_springs.html>`__
`polynomial rings <http://www.math.uiuc.edu/Macaulay2/Documentation/_polynomial_springs.html>`__
with a `monomial ordering <http://www.math.uiuc.edu/Macaulay2/Documentation/_monomial_sporderings.html>`__.
The core algorithms are accessible through a versatile high level interpreted user
`language <http://www.math.uiuc.edu/Macaulay2/Documentation/___The_sp__Macaulay2_splanguage.html>`__
with a powerful `debugger <http://www.math.uiuc.edu/Macaulay2/Documentation/_the_spdebugger.html>`__
supporting the `creation <http://www.math.uiuc.edu/Macaulay2/Documentation/_making_spnew_spclasses.html>`__
of new `classes <http://www.math.uiuc.edu/Macaulay2/Documentation/_what_spa_spclass_spis.html>`__
of mathematical objects and the
`installation of methods <http://www.math.uiuc.edu/Macaulay2/Documentation/_installing_spmethods.html>`__
for computing specifically with them.

*Macaulay2* can compute `Betti numbers <http://www.math.uiuc.edu/Macaulay2/Documentation/_betti_lp__Graded__Module_rp.html>`__,
`Ext <http://www.math.uiuc.edu/Macaulay2/Documentation/___Ext.html>`__,
`cohomology of coherent sheaves <http://www.math.uiuc.edu/Macaulay2/Documentation/___H__H%5E__Z__Z_sp__Coherent__Sheaf.html>`__ on projective varieties,
`primary decomposition of ideals <http://www.math.uiuc.edu/Macaulay2/Documentation/_primary_spdecomposition.html>`__,
`integral closure <http://www.math.uiuc.edu/Macaulay2/doc/Macaulay2-1.9.2/share/doc/Macaulay2/IntegralClosure/html/>`__
of `rings <http://www.math.uiuc.edu/Macaulay2/Documentation/_rings.html>`__,
and `more <http://www.math.uiuc.edu/Macaulay2/Documentation/>`__.


Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive session with the :code:`qsh` or :code:`qrsh` command.
To make Macaulay2 available in this session, run one of the following command:

.. code-block:: none

      module load apps/gcc/5.3/macaulay2/1.9.2
      module load apps/gcc/5.3/macaulay2/1.10
      module load apps/gcc/5.3/macaulay2/1.12.0.1

You should then be able to start an interactive Macaulay2 session using

.. code-block:: none

      M2

Testing
-------

Macaulay2 can perform an extensive series of self-tests at compile-time (by running :code:`make check`).  These tests were performed and passed when building Macaulay2 on Iceberg.

Installation notes and modulefile
---------------------------------

Version 1.12.0.1
################

* Built using :download:`install.sh </iceberg/software/install_scripts/apps/gcc/5.3/macaulay2/1.12.0.1/install.sh>`
* Made available to users via :download:`this modulefile </iceberg/software/modulefiles/apps/gcc/5.3/macaulay2/1.12.0.1>` located on the system at ``/usr/local/modulefiles/apps/gcc/5.3/macaulay2/1.12.0.1``

* The install script and module file enable several other module files to provide:
    * GCC's compiler and Fortran library (:code:`compilers/gcc/5.3`)
    * the Boost C++ libraries (:code:`libs/gcc/5.2/boost/1.59`)
    * LAPACK linear algebra routines (:code:`libs/gcc/lapack/3.3.0`)
    * LZMA (xz) (de)compression (:code:`apps/gcc/4.4.7/xzutils/5.2.4`)

Version 1.9.2
#############

* Built using :download:`install.sh </iceberg/software/install_scripts/apps/gcc/5.3/macaulay2/1.9.2/install.sh>`
* Built using a specific git commit (aedebfc1e6326416cb01598e09e8b4dbdc76c178) rather than the tagged 1.9.2 release as Macaulay2's self-tests failed for the latter.  There are few differences between this commit and the tagged 1.9.2 release.
* Made available to users via :download:`this modulefile </iceberg/software/modulefiles/apps/gcc/5.3/macaulay2/1.9.2>` located on the system at ``/usr/local/modulefiles/apps/gcc/5.3/macaulay2/1.9.2``

* The install script and module file enable several other module files to provide:
    * OpenMPI (:code:`mpi/intel/openmpi/1.10.0`)
    * GCC's compiler and Fortran library (:code:`compilers/gcc/5.3`)
    * the Boost C++ libraries (:code:`libs/gcc/5.2/boost/1.59`)
    * LAPACK linear algebra routines (:code:`libs/gcc/lapack/3.3.0`)
    * LZMA (xz) (de)compression (:code:`apps/gcc/4.4.7/xzutils/5.2.2`)
