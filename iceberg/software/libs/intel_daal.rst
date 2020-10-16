.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _iceberg_intel_daal:

Intel Data Analytics Acceleration Library 
=========================================

Intel's Data Analytics Acceleration Library (DAAL) provides functions for data
analysis (characterization, summarization, and transformation) and machine
learning (regression, classification, and more).

Parallel Studio Composer Edition version
----------------------------------------

DAAL can be used with and without :ref:`other Parallel Studio packages
<iceberg_intel_parallel_studio>`.
To access it: ::

        module load libs/binlibs/intel-daal/2017.0

Licensing and availability
--------------------------

See the information on :ref:`Parallel Studio licensing
<iceberg_intel_parallel_studio>`.

Installation Notes
------------------

The following notes are primarily for system administrators.

**Intel DAAL 2017.0**

Installed as part of :ref:`Parallel Studio Composer Edition 2017
<iceberg_intel_parallel_studio>`.

:download:`This modulefile
</iceberg/software/modulefiles/libs/binlibs/intel-daal/2017.0>` was installed as
``/usr/local/modulefiles/libs/binlibs/intel-daal/2017.0``.
