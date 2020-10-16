.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 


SU2
====

.. sidebar:: SU2
   
   :Version: 4.1.0
   :Dependancies: gcc/4.4.7 mpi/gcc/openmpi/1.10.1
   :URL: http://su2.stanford.edu/
   :Documentation: https://github.com/su2code/SU2/wiki


The SU2 suite is an open-source collection of C++ based software tools for performing Partial Differential Equation (PDE) analysis and solving PDE constrained optimization problems.

Usage
-----

SU2 can be activated using the module file::

    module load apps/gcc/4.4.7/su2


Installation notes
------------------

Su2 was compiled using the :download:`install_su2.sh </iceberg/software/install_scripts/apps/gcc/4.4.7/su2/install_su2.sh>` script. 

SU2 also has it's own version of CGNS which was compiled with the script :download:`install_cgns.sh </iceberg/software/install_scripts/apps/gcc/4.4.7/su2/install_cgns.sh>`.

The module file is :download:`4.1.0 </iceberg/software/modulefiles/apps/gcc/4.4.7/su2/4.1.0>`.
