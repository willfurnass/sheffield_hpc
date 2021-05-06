.. toctree::
    :hidden:
    :maxdepth: 1
    :glob:

==============================
CPU Allocation Limits
==============================


.. list-table:: CPU Allocation Limits Table
   :widths: 23 23 23 30
   :header-rows: 1

   * - Scheduler Type
     - No. CPU Cores Available |br| |br| Interactive Job |br| (Default/ Min / Max )
     - No. CPU Cores Available |br| |br| Batch Job |br| (Default/ Min / Max )
     - Submission Argument

   * - SGE (ShARC)
     - 1 / 1 / 16
     - 1 / 1 / ~1536
     - ``-pe <env> <nn>``

   * - SLURM (Bessemer)
     - 1 / 1 / 40
     - 1 / 1 / 40
     - ``-c <nn>``
