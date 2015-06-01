R
=

.. sidebar:: R
   
   :Support Level: bronze
   :Dependancies: BLAS
   :URL: http://www.r-project.org/ 
   :Documentation: http://www.r-project.org/  

R is a statistical computing language.

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command.

The lastest version of R can be loaded with ::

        module load apps/R

Alternatively, you can load a specific version of R using one of the following ::

        module load apps/R/3.1.2
        module load apps/R/3.2.0

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
