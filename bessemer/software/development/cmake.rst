.. _cmake_bessemer:

CMake
=====

CMake is a build tool commonly used when compiling other libraries.

Usage
-----

CMake can be loaded with: ::

    module load CMake/3.13.3-GCCcore-8.2.0


CMake has a run-time dependency on `libstdc++` so the above also needs to
(and does) load the :ref:`GCC compiler <gcc_bessemer>` version 8.2.0,
plus versions of the ncurses, zlib, bzip2 and cURL libraries.

Usage of CMake often involves: 

1. Creating and ``cd``-ing into a dedicated build directory within a source tree then
2. Running something like ``cmake -DSOME_OPTION -DANOTHER_OPTION ..``
