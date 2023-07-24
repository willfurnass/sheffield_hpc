The ``seff`` script can be used as follows with the job's ID to give summary of important job info : 

.. code-block:: console

    $ seff job-id

For example:

.. code-block:: console

    $ seff 63888
      Job ID: 63888
      Cluster: stanage.alces.network
      User/Group: a_person/a_group
      State: COMPLETED (exit code 0)
      Cores: 1
      CPU Utilized: 00:00:00
      CPU Efficiency: 0.00% of 00:00:26 core-walltime
      Job Wall-clock time: 00:00:26
      Memory Utilized: 0.00 MB (estimated maximum)
      Memory Efficiency: 0.00% of 512.00 MB (512.00 MB/core)
