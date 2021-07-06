Dakota
========

.. sidebar:: Dakota

   :Version: 6.14.0
   :Dependencies: https://dakota.sandia.gov/content/install-required-tools-and-libraries
   :URL: https://dakota.sandia.gov/
   :Documentation: https://dakota.sandia.gov/quickstart.html

Dakota is a general-purpose software toolkit for 
performing systems analysis and design on high performance computers. 
Dakota provides algorithms for design optimization, uncertainty quantification, 
parameter estimation, design of experiments, and sensitivity analysis, as well 
as a range of parallel computing and simulation interfacing services.

Usage
-----

Dakota 6.14.0 can be activated using the module file

.. code-block:: bash

    module load apps/dakota/6.14.0/gcc-8.2-cmake-3.17.1


.. caution::

    * Use of the ``#$ -V`` SGE option will instruct SGE to import your current terminal environment variables to be imported - **CAUTION** - this may not be desirable as it could cause conflicts between modules.


.. hint::

    * It is recommended that users submit each Dakota subtask: ``simulator_script.sh`` should submit numerous SGE jobs.
    * `Bash Heredoc scripting <https://linuxize.com/post/bash-heredoc/>`_ may be useful or required as part of submitting Dakota subtasks to the scheduler engine with ``simulator_script.sh`` to adjust resource requests or subscript variables.

Interactive Usage
-----------------


Batch Usage
------------



Installation notes
------------------

Dakota 6.14.0 was compiled from source using the
:download:`install_dakota.sh </sharc/software/install_scripts/apps/dakota/6.14.0/install_dakota.sh>` script.
