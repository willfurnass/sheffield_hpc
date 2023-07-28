.. group-tab:: Stanage

        For a quick summary, the :ref:`seff` command can be used as follows with the job's ID to give summary of important job info including the memory usage / efficiency:

        .. code-block:: console
            :emphasize-lines: 11,12

            $ seff 64626
            Job ID: 64626
            Cluster: stanage.alces.network
            User/Group: a_user/clusterusers
            State: COMPLETED (exit code 0)
            Nodes: 2
            Cores per node: 1
            CPU Utilized: 00:02:37
            CPU Efficiency: 35.68% of 00:07:20 core-walltime
            Job Wall-clock time: 00:03:40
            Memory Utilized: 137.64 MB (estimated maximum)
            Memory Efficiency: 1.71% of 7.84 GB (3.92 GB/core)

        For more specific info, you can use the **sacct** / **sstat** commands:

        While a job is still running find out its job id by: ::

            sacct

        And check its current usage of memory by: ::

            sstat job_id --format='JobID,MaxVMSize,MaxRSS'

        If your job has already finished you can list the memory usage with sacct: ::

            sacct --format='JobID,Elapsed,MaxVMSize,MaxRSS'

        It is the **MaxVMSize** / **MaxRSS** figures that you will need to use to determine the ``--mem=`` parameter for your next job.

.. group-tab:: Bessemer
        
        For a quick summary, the :ref:`seff` command can be used as follows with the job's ID to give summary of important job info including the memory usage / efficiency:

        .. code-block:: console
            :emphasize-lines: 11,12

            $ seff 64626
            Job ID: 64626
            Cluster: bessemer
            User/Group: a_user/a_group
            State: COMPLETED (exit code 0)
            Nodes: 1
            Cores per node: 2
            CPU Utilized: 00:02:37
            CPU Efficiency: 35.68% of 00:07:20 core-walltime
            Job Wall-clock time: 00:03:40
            Memory Utilized: 137.64 MB (estimated maximum)
            Memory Efficiency: 1.71% of 7.84 GB (3.92 GB/core)

        For more specific info, you can use the **sacct** / **sstat** commands:

        While a job is still running find out its job id by: ::

            sacct

        And check its current usage of memory by: ::

            sstat job_id --format='JobID,MaxVMSize,MaxRSS'

        If your job has already finished you can list the memory usage with sacct: ::

            sacct --format='JobID,Elapsed,MaxVMSize,MaxRSS'

        It is the **MaxVMSize** / **MaxRSS** figures that you will need to use to determine the ``--mem=`` parameter for your next job.

.. group-tab:: ShARC

        You can detect the memory used by your job while it is running by using the **qstat** command for SGE as follows:
        
        While a job is still running find out its job id by: ::

            qstat

        And check its current usage of memory by: ::

            qstat -F -j job_id | grep mem

        If your job has already finished you can list the memory usage with **qacct**: ::

            qacct -j job_id | grep vmem

        The reported figures will indicate:

        * the currently used memory ( vmem ).
        * maximum memory needed since startup ( maxvmem ).

        It is the **maxvmem** figure that you will need to use to determine the ``-l rmem=`` parameter for your next job.

        For example:

        .. code-block:: console
            :emphasize-lines: 2

            $ qacct -j 628 | grep vmem
            maxvmem      4.807G
            category     -u username -l h_vmem=8G -pe smp 1 -P SHEFFIELD

        .. code-block:: console
            :emphasize-lines: 2

            $ qstat -F -j 628 | grep vmem
            usage         1:            cpu=77:57:45, mem=61962.33471 GB s, io=38.68363 GB, vmem=4.334G, maxvmem=4.807G