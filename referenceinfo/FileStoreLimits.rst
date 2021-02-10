=========================================================
Filestore Limits / file store performance characteristics
=========================================================

Every HPC user is allocated a file-storage area of their own. Please read the section on  `filestore allocation <filestore>`_ on the clusters for further information. Any attempt to exceed this allocation during the execution of a job can have disastrous consequences. This is because any program or package writing into files will produce a fatal error and stop if the filestore limit happens to be exceeded during that operation. 

Filestore limits are not associated with jobs and can not be specified while submitting a job. Users must make sure that there is sufficient spare space in their filestore areas before submitting any job that is going to produce large amounts of output.
quota command can be used to check your current filestore allocation and usage.

Each filestore has the relevant detail and performance characteristics listed within the section on `filestore allocation <filestore>`_, this will indicate where your program is best suited to run from.

