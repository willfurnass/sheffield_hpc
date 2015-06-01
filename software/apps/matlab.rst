MATLAB
======

.. sidebar:: MATLAB 
   
   :Versions:  2013a , 2013b , 2014a, 2015a
   :Support Level: FULL 
   :Dependancies: None
   :URL: http://uk.mathworks.com/products/matlab 
   :Local URL:  http://www.shef.ac.uk/wrgrid/software/matlab
   :Documentation: http://uk.mathworks.com/help/matlab

Scientific Computing and Visualisation 

Usage
-----
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` command.

The lastest version of MATLAB (currently 2015a) is made available with the command

.. code-block:: none

        module load apps/matlab

Alternatively, you can load a specific version with one of of the following commands

.. code-block:: none

       module apps/matlab/2013a
       module apps/matlab/2013b
       module apps/matlab/2014a

You can then run MATLAB by entering :code:`matlab` 

Installation notes
------------------

Requires the floating license server licserv4.shef.ac.uk to serve the licenses 
for the version of MATLAB to be installed ( or higher versions ) .
An install script and associated files are downloadable from Mathworks site along with all the required toolbox specific installation files. 


