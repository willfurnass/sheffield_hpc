IDL
===

.. sidebar:: IDL

   :Version: 8.5
   :Dependencies: Java for GUI. Module loaded for Java 1.8.0_102.
   :URL: http://www.exelisvis.co.uk/ProductsServices/IDL.aspx
   :Documentation: http://www.exelisvis.com/docs/using_idl_home.html

IDL is a data analysis language that first appeared in 1977.

Usage
-----
If you wish to use the IDLDE then you may need to request more memory for the interactive session using something like ``qsh -l mem=8G``.

IDL can be activated using the module file::

    module load apps/idl/8.5/binary

Then run using ``idl`` or ``idlde`` for the interactive development environment.

Installation notes
------------------

Installation of IDL 8.5 on Sharc was a binary installation.
IDL 8.5 was installed using the
:download:`install_idl.sh </sharc/software/install_scripts/apps/idl/8.5/binary/install_idl.sh>` script; the module
file is
:download:`/usr/local/modulefiles/apps/idl/8.5/binary </sharc/software/modulefiles/apps/idl/8.5/binary>`.
