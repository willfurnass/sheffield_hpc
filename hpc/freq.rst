

Frequently Asked SGE Questions
==============================

**How many jobs can I submit at any one time**

You can submit up to 2000 jobs to the cluster, and the scheduler will allow up to 200 of your jobs to run simultaneously (we occasionally alter this value depending on the load on the cluster).

**How do I specify the processor type on Iceberg?**

Add the following line to your submission script ::

    #$ -l arch=intel-e5-2650v2

This specifies nodes that have the Ivybridge `E5-2650 CPU <http://ark.intel.com/products/75269/Intel-Xeon-Processor-E5-2650-v2-20M-Cache-2_60-GHz>`_.
All such nodes on Iceberg have 16 cores.

To only target the older, 12 core nodes that contain `X5650 CPUs <http://ark.intel.com/products/47922/Intel-Xeon-Processor-X5650-12M-Cache-2_66-GHz-6_40-GTs-Intel-QPI>`_ add the following line to your submission script ::

    #$ -l arch=intel-x5650


**How do I specify multiple email addresses for job notifications?**

Specify each additional email with it's own `-M` option ::

  #$ -M foo@example.com
  #$ -M bar@example.com

**How do you ensure that a job starts after a specified time?**

Add the following line to your submission script ::

    #$ -a time

but replace ``time`` with a time in the format MMDDhhmm

For example, for 22nd July at 14:10, you’d do ::

    #$ -a 07221410

This won’t guarantee that it will run precisely at this time since that depends on available resources. It will, however, ensure that the job runs *after* this time. If your resource requirements aren’t too heavy, it will be pretty soon after. When I tried it, it started about 10 seconds afterwards but this will vary.
