.. _cmake_stanage:

CMake
=====

CMake is a build tool commonly used when compiling other libraries.

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

CMake can be loaded with one of:

.. include:: /referenceinfo/imports/stanage/packages/cmake-ml-el7-icelake-znver-stanage.rst

Usage of CMake often involves: 

1. Creating and ``cd``-ing into a dedicated build directory within a source tree then
2. Running something like ``cmake -DSOME_OPTION -DANOTHER_OPTION ..``

.. include:: /referenceinfo/imports/stanage/packages/cmake-dpnd-el7-icelake-znver-stanage.rst
