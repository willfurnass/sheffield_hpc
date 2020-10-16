.. Warning:: 
    Iceberg reaches end-of-life on **30th November 2020.**
    If you are running jobs on Iceberg then you need to take urgent action to ensure that your jobs/scripts will run on ShARC or Bessemer. 
 
    If you have never used ShARC or Bessemer then now is the time to test your scripts.
    Not all software on Iceberg is available on ShARC/Bessemer. 


IDL
===

.. sidebar:: IDL
   
   :Version: 8.5
   :Support Level: extras
   :Dependancies: java
   :URL: http://www.exelisvis.co.uk/ProductsServices/IDL.aspx 
   :Documentation: http://www.exelisvis.com/docs/using_idl_home.html

IDL is a data analysis language that first appeared in 1977.

Usage
-----
If you wish to use the IDLDE then you may need to request more memory for the interactive 
session using something like ``qsh -l rmem=8G``.

IDL versions can be activated using specific module files::

	module load apps/idl/8.5
	module load apps/idl/8.5ssw
	module load apps/idl/8.4

then run using ``idl`` or ``idlde`` for the interactive development environment. Note apps/idl/8.5ssw is the SolarSoft IDL environment.

The IDL licence is restricted - please contact us via helpdesk if you need to use this software.

Installation notes
------------------

Extract the supplied linux x86-64 tar file, run the install shell script and point it at the install directory.
