MrBayes
=======

.. sidebar:: MrBayes

   :Version:  3.2.6
   :URL: https://sourceforge.net/projects/mrbayes/

MrBayes is a program for the Bayesian estimation of phylogeny.

Interactive Usage
-----------------

After connecting to Iceberg (see :ref:`ssh`), :ref:`start an interactive session <sched_interactive>` then

load a specific version of MrBayes with: ::

   module load apps/gcc/4.4.7/mrbayes/3.2.6

This command makes the MrBayes ``mb`` executable available to your session.

Installation notes
------------------

MrBayes was installed with gcc 4.4.7: ::

    tar -xvzf ./mrbayes-3.2.6.tar.gz
    ls
    cd mrbayes-3.2.6
    autoconf
    ./configure --with-beagle=/usr/local/packages6/libs/gcc/4.4.7/beagle/2.1.2 --prefix=/usr/local/packages6/apps/gcc/4.4.7/mrbayes/3.2.6/
    make
    make install

:download:`This modulefile </iceberg/software/modulefiles/apps/gcc/4.4.7/mrbayes/3.2.6>`
was installed as ``/usr/local/modulefiles/apps/gcc/4.4.7/mrbayes/3.2.6``
