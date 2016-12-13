.. _cmake_sharc:

CMake
=====

CMake is a build tool commonly used when compiling other libraries.

Usage
-----

CMake can be loaded with: ::

    module load dev/cmake/3.7.1/gcc-4.9.4

NB CMake has a run-time dependency on `libstdc++` so the above also needs to
(and does) load the :ref:`GCC compiler <gcc_sharc>` version 4.9.4.

Usage often involves: 

1. Creating and ``cd``-ing into a dedicated build directory within a source tree then
2. Running something like ``cmake -DSOME_OPTION -DANOTHER_OPTION ..``

Installation
------------

Version 3.7.1
^^^^^^^^^^^^^

1. Install using :download:`this script </sharc/software/install_scripts/dev/cmake/3.7.1/gcc-4.9.4/install.sh>`
2. Install :download:`this modulefile </sharc/software/modulefiles/dev/cmake/3.7.1/gcc-4.9.4>` as ``/usr/local/modulefiles/dev/cmake/3.7.1/gcc-4.9.4``
