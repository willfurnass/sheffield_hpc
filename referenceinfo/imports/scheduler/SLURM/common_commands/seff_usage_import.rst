The ``seff`` script can be used as follows with the job's ID to give summary of important job info : 

.. code-block:: console

    $ seff job-id

For example, on the Stanage cluster:

.. code-block:: console

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
