Bedtools
========

.. sidebar:: Bedtools

   :Versions:  2.31.0
   :Dependencies: GCC/12.3.0
   :URL: https://bedtools.readthedocs.org/en/latest/

Collectively, the bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks.

Interactive usage
-----------------
.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

The latest version of bedtools (currently version 2.31.0) is made available with the command:

.. code-block:: none

        module load BEDTools/2.31.0-GCC-12.3.0


After this any of the bedtools commands can be run from the prompt.



Installation notes
------------------

Bed tools was installed using EasyBuild/4.9.4, build details can be found in ``$EBROOTBEDTOOLS/easybuild`` with the module loaded.
