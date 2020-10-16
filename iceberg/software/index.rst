.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _iceberg-software:

Software on iceberg
===================

These pages list the software available on iceberg. If you notice an error or
an omission, or wish to request new software please email
the research computing team at `research-it@sheffield.ac.uk <research-it@sheffield.ac.uk>`_.

There are many versions of applications, libraries and compilers installed on the iceberg cluster. In order to avoid conflict between these software items we employ a system called modules. Having logged onto a worker nodes, users are advised to select the software (and the version of the software) they intend to use by using the ``module`` command.

.. toctree::
    :maxdepth: 2

    ../../hpc/modules
    apps/index
    libs/index
    compilers/index
    mpi/index
