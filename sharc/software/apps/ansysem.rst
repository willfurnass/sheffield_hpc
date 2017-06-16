Ansys EM
========

.. sidebar:: Ansys EM
   
   :Versions: 16.1, 17.2
   :Dependencies: For integration with Ansys Workbench requires Ansys 16.1, 17.2
   :URL: http://www.ansys.com 
   :Local URL: http://www.shef.ac.uk/cics/research/software/fluent


ANSYS Electromagnetics Suite for Linux/Unix.


Usage
-----

Ansys EM can be activated using the module files::

    module load apps/ansysem/16.1
    module load apps/ansysem/17.2

The Ansys EM exectuable is ``ansysedt``. The following is an example batch submission script which is submitted to the queue by typing ``qsub my_job.sh``::

    #!/bin/bash
    #$ -cwd
    #$ -l h_rt=00:30:00
    #$ -l rmem=1G
    #$ -pe mpi 8

    module load apps/ansysem/16.1

    ansysedt -ng -BatchSolve -Distributed -machinelist num=8 -batchoptions 'HFSS/HPCLicenseType'='pool' -batchoptions 'HFSS-IE/HPCLicenseType'='pool' Tee.aedt
	
The script requests 8 cores using the MPI parallel environment with a runtime of 30 mins and 1G of real memory per core. The Ansys EM input file is ``Tee.aedt`` and batch options ``'HFSS/HPCLicenseType'='pool'`` and ``'HFSS-IE/HPCLicenseType'='pool'`` change the HPC licencing from "pack" (the default) to "pool".	
	
Installation notes
------------------

Ansys EM 16.1 was installed using the
:download:`install_openfoam.sh </sharc/software/install_scripts/apps/ansysem/16.1/install_ansysem.sh>` script, the module
file is
:download:`16.1 </sharc/software/modulefiles/apps/ansysem/16.1>`.
Ansys EM 17.2 was installed using the
:download:`install_openfoam.sh </sharc/software/install_scripts/apps/ansysem/17.2/install_ansysem.sh>` script, the module
file is
:download:`17.2 </sharc/software/modulefiles/apps/ansysem/17.2>`. The binary installations were tested using the above batch submission script.
