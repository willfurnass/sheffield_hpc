.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 

MOMFBD
======

.. sidebar:: MOMFBD

   :Version: 2016-04-14
   :URL: http://dubshen.astro.su.se/wiki/index.php?title=MOMFBD


MOMFBD or Multi-Object Multi-Frame Blind Deconvolution is a image processing
application for removing the effects of atmospheric seeing from solar image
data.

   
Usage
-----

The MOMFBD binaries can be added to your path with: ::

        module load apps/gcc/5.2/MOMFBD

Installation notes
------------------

MOMFBD was installed using gcc 5.2. 
Using :download:`this script </iceberg/software/install_scripts/apps/gcc/5.2/MOMFBD/install_momfbd.sh>` 
and is loaded by :download:`this modulefle </iceberg/software/modulefiles/apps/gcc/5.2/momfbd/2016.04.14>`.
