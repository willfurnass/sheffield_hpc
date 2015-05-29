R
=

.. sidebar:: R
   
   :Maintainer: Stuart Mumford
   :Email: stuart.mumford@sheffield.ac.uk
   :Version: 3.1.2
   :Support Level: bronze
   :Dependancies: BLAS
   :URL: http://www.r-project.org/ 
   :Documentation: http://www.r-project.org/  

R is a statistical computing language.

Usage
-----

R can be loaded with ::

        module load apps/R

R can then be run with ::

        $ R

Installing
----------

R was compiled from source using the following commands::

        $ module load libs/gcc/lapack
        $ module load libs/gcc/blas
        $ ./configure --prefix /usr/local/packages6/R/3.1.2 --use-blas --use-lapack
        $ make -j 6
        $ make install

It seems like it did not pick up the lapack library, but it did pick up the BLAS lib ok.
