.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _scotch:

SCOTCH
======

.. sidebar:: SCOTCH

   :Latest version: 6.0.4
   :URL: http://www.labri.fr/perso/pelegrin/scotch/
   :Location: /usr/local/packages6/libs/gcc/5.2/scotch/6.0.4

Software package and libraries for sequential and parallel graph partitioning,
static mapping and clustering, sequential mesh and hypergraph partitioning, and
sequential and parallel sparse matrix block ordering.

Usage
-----
To make this library available, run the following module command ::

        module load libs/gcc/5.2/scotch/6.0.4

Installation notes
------------------
This section is primarily for administrators of the system. SCOTCH 6.0.4 was compiled with gcc 5.2 using the bash file :download:`install_scotch.sh </iceberg/software/install_scripts/libs/gcc/5.2/install_scotch.sh>`
