.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

.. _iceberg_intel_ipp:

Intel Integrated Performance Primitives
=======================================

Integrated Performance Primitives (IPP) are "high-quality, production-ready,
low-level building blocks for image processing, signal processing, and data
processing (data compression/decompression and cryptography) applications."

Parallel Studio Composer Edition version
----------------------------------------

IPP can be used with and without :ref:`other Parallel Studio packages
<iceberg_intel_parallel_studio>`.
To access it: ::

        module load libs/binlibs/intel-ipp/2017.0

Licensing and availability
--------------------------

See the information on :ref:`Parallel Studio licensing
<iceberg_intel_parallel_studio>`.

Installation Notes
------------------

The following notes are primarily for system administrators.

**Intel IPP 2017.0**

Installed as part of :ref:`Parallel Studio Composer Edition 2017
<iceberg_intel_parallel_studio>`.

:download:`This modulefile
</iceberg/software/modulefiles/libs/binlibs/intel-ipp/2017.0>` was installed as
``/usr/local/modulefiles/libs/binlibs/intel-ipp/2017.0/binary``.
