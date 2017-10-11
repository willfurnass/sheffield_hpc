.. _lapack_sharc:

LAPACK
======
LAPACK is written in Fortran 90 and provides routines for solving systems of simultaneous linear equations, least-squares solutions of linear systems of equations, eigenvalue problems, and singular value problems. The associated matrix factorizations (LU, Cholesky, QR, SVD, Schur, generalized Schur) are also provided, as are related computations such as reordering of the Schur factorizations and estimating condition numbers. Dense and banded matrices are handled, but not general sparse matrices. In all areas, similar functionality is provided for real and complex matrices, in both single and double precision.

Usage
-----
ShARC runs on the CentOS 7.x flavour of Linux which provides a standard build of LAPACK via RPM files.  We have made this version of LAPACK available via the following module command ::

    $ module load libs/lapack/3.4.2-5/binary

Installation notes
------------------
This section is primarily for administrators of the system.

**Version 3.4.2-5 binary**

LAPACK was extracted from the RPMs for **lapack** and **lapack-devel**

# :download:`Install script </sharc/software/install_scripts/libs/lapack/3.4.2-5/binary/install.sh>`
# :download:`Module File </sharc/software/modulefiles/libs/lapack/3.4.2-5/binary>` as ``/usr/local/modulefiles/modulefiles/libs/lapack/3.4.2-5/binary``
