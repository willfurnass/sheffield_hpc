.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 


Phyluce
=======

.. sidebar:: Phyluce
   
   :Version: 1.5.0
   :Dependancies: apps/python/conda
   :URL: https://github.com/faircloth-lab/phyluce 
   :Documentation: http://phyluce.readthedocs.io/

Phyluce (phy-loo-chee) is a software package that was initially developed for analyzing data collected from ultraconserved elements in organismal genomes.

Usage
-----
Phyluce can be activated using the module file::

    module load apps/binapps/phyluce/1.5.0

Phyluce makes use of the ``apps/python/conda`` module, therefore this module will be loaded by loading Phyluce.
As Phyluce is a Python package your default Python interpreter will be changed by loading Phyluce.

Installation notes
------------------

As root: ::

      $ module load apps/python/conda
      $ conda create -p /usr/local/packages6/apps/binapps/conda/phyluce python=2
      $ source activate /usr/local/packages6/apps/binapps/conda/phyluce
      $ conda install -c https://conda.binstar.org/faircloth-lab phyluce

This installs Phyluce as a conda environment in the ``/usr/local/packages6/apps/binapps/conda/phyluce`` folder, 
which is then loaded by the module file :download:`phyluce/1.5.0 </iceberg/software/modulefiles/apps/binapps/phyluce/1.5.0>`, 
which is a modification of the anaconda module files.
