Boost C++ Library v1.58
=======================

.. sidebar:: Boost version 1.58
   
   :Version: 1.58
   :Support Level: Bronze
   :Dependancies: libs/gcc/4.8.2/libunistring/0.9.5, libs/gcc/4.8.2/icu/55, compilers/gcc/4.8.2
   :URL: www.boost.org
   :Documentation: http://www.boost.org/doc/libs/1_58_0/
   :Location: /usr/local/packages6/libs/gcc/4.8.2/boost/1.58.0/ 

Boost provides free peer-reviewed portable C++ source libraries.

Usage
-----
This build of the Boost library requires gcc version 4.8.2. To make the compiler and library available, run the following module commands

`module load compilers/gcc/4.8.2`
`module load libs/gcc/4.8.2/boost/1.58`

Installing
----------
.. code-block:: none

    module load compilers/gcc/4.8.2
    module load libs/gcc/4.8.2/libunistring/0.9.5
    module load libs/gcc/4.8.2/icu/55
    tar -xvzf ./boost_1_58_0.tar.gz
    cd boost_1_58_0
    ./bootstrap.sh --prefix=/usr/local/packages6/libs/gcc/4.8.2/boost/1.58.0/

It complained that it could not find the icu library but when I ran

.. code-block:: none
   
    ./b2 install --prefix=/usr/local/packages6/libs/gcc/4.8.2/boost/1.58.0
    
It said that it had detected the icu library and was compiling it in

Testing
-------
TODO

Module File
-----------
Module File Location: `/usr/local/modulefiles/libs/gcc/4.8.2/boost/1.58`

.. code-block:: none

        #%Module1.0#####################################################################
        ##
        ## libunistring 0.9.5 module file
        ##

        ## Module file logging
        source /usr/local/etc/module_logging.tcl
        ##

        module load libs/gcc/4.8.2/libunistring/0.9.5
        module load libs/gcc/4.8.2/icu/55

        proc ModulesHelp { } {
                puts stderr "Makes the Boost 1.58 library available"
        }

        set BOOST_DIR /usr/local/packages6/libs/gcc/4.8.2/boost/1.58.0

        module-whatis   "Makes the Boost 1.58 library available"

        prepend-path LD_LIBRARY_PATH $BOOST_DIR/lib
        prepend-path CPLUS_INCLUDE_PATH $BOOST_DIR/include
        prepend-path LIBRARY_PATH $BOOST_DIR/lib


