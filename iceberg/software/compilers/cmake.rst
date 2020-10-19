.. include:: ../../../iceberg-eol.rst 

.. _cmake_iceberg:

CMake
=====

CMake is a build tool commonly used when compiling other libraries.

Usage
-----

CMake can be loaded with: ::

    module load compilers/cmake/3.3.0

Usage often involves: 

1. Creating and ``cd``-ing into a dedicated build directory within a source tree then
2. Running something like ``cmake -DSOME_OPTION -DANOTHER_OPTION ..``

Installation
------------

Run the following commands::

    module load apps/python/2.7

    ./bootstrap --prefix=/usr/local/packages6/cmake/3.3.0/
    --mandir=/usr/local/packages6/cmake/3.3.0/man --sphinx-man

    gmake -j8

    gmake install

Module file: ::

        module-whatis    loads the necessary `cmake-3.3.0' library paths 
        prepend-path     PATH /usr/local/packages6/cmake/3.3.0/bin 
        prepend-path     MANPATH /usr/local/packages6/cmake/3.3.0/man 
