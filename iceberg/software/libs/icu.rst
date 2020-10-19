.. include:: ../../../iceberg-eol.rst 

icu - International Components for Unicode
==========================================

.. sidebar:: icu

   :Latest Version: 58.1
   :URL: http://site.icu-project.org/

ICU is a mature, widely used set of C/C++ and Java libraries providing Unicode and Globalization support for software applications. ICU is widely portable and gives applications the same results on all platforms and between C/C++ and Java software.

ICU is released under a nonrestrictive open source license that is suitable for use with both commercial software and with other open source or free software.

Usage
-----
Version 58.1 of the icu library for C requires gcc version 4.9.2 (for the C++ standard library); To make the library **and this compiler** available, run the following: ::

        module load libs/gcc/4.9.2/icu/58.1

Version 55 of the icu library for C requires gcc version 4.8.2. To make the compiler and library available, run the following module commands: ::

        module load compilers/gcc/4.8.2
        module load libs/gcc/4.8.2/icu/55

Version 42 of the icu library for C uses gcc version 4.4.7 which is the default on Iceberg. To make the library available, run the following command: ::

        module load libs/gcc/4.4.7/icu/42

Installation Notes
------------------
This section is primarily for administrators of the system.

Version 58.1
^^^^^^^^^^^^

This Icu 58.1 build links against the GCC 4.9.2 C++ standard library and was installed as a dependency of `boost_iceberg` (build using the same C++ standard library); Boost in turn was installed as a dependency of Gromacs 2016.1. ::

        module load compilers/gcc/4.9.2
        tar -xvzf icu4c-58_1-src.tgz
        cd icu/source
        ./runConfigureICU Linux/gcc --prefix=/usr/local/packages6/libs/gcc/4.9.2/icu/58.1/
        make
        make check
        make install

Version 55
^^^^^^^^^^

Icu 55 is a pre-requisite for the version of `boost_iceberg` required for an experimental R module used by one of our users. It was built using gcc 4.8.2. ::

        module load compilers/gcc/4.8.2
        tar -xvzf icu4c-55_1-src.tgz
        cd icu/source
        ./runConfigureICU Linux/gcc --prefix=/usr/local/packages6/libs/gcc/4.8.2/icu/55/
        make
        make check
        make install

Version 42
^^^^^^^^^^

Icu version 42 was originally installed as a system RPM. This install moved icu to a module-based install. ::

        tar -xvzf icu4c-4_2_1-src.tgz
        cd icu/source/
        ./runConfigureICU Linux/gcc --prefix=/usr/local/packages6/libs/gcc/4.4.7/icu/42
        make
        make check
        make install

Testing
-------

Version 58.1
^^^^^^^^^^^^

Last few lines of output from `make check` were: ::

        ALL TESTS SUMMARY:
        All tests OK:  testdata intltest iotest cintltst
        make[1]: Leaving directory `/home/te1st/icu_gcc_4_9_2/icu/source/test'
        make[1]: Entering directory `/home/te1st/icu_gcc_4_9_2/icu/source'
        verifying that icu-config --selfcheck can operate
        verifying that make -f Makefile.inc selfcheck can operate
        PASS: config selfcheck OK

Version 55
^^^^^^^^^^

Last few lines of output from `make check` were: ::

        [All tests passed successfully...]
        Elapsed Time: 00:00:00.086
        make[2]: Leaving directory `/home/te1st/icu/icu/source/test/letest'
        ---------------
        ALL TESTS SUMMARY:
        All tests OK:  testdata intltest iotest cintltst letest
        make[1]: Leaving directory `/home/te1st/icu/icu/source/test'
        make[1]: Entering directory `/home/te1st/icu/icu/source'
        verifying that icu-config --selfcheck can operate
        verifying that make -f Makefile.inc selfcheck can operate
        PASS: config selfcheck OK
        rm -rf test-local.xml

Version 42
^^^^^^^^^^

Last few lines of output from `make check` were: ::

        All tests passed successfully...]
        Elapsed Time: 00:00:12.000
        make[2]: Leaving directory `/home/te1st/icu/source/test/cintltst'
        ---------------
        ALL TESTS SUMMARY:
        ok:  testdata iotest cintltst
        ===== ERRS:  intltest
        make[1]: *** [check-recursive] Error 1
        make[1]: Leaving directory `/home/te1st/icu/source/test'
        make: *** [check-recursive] Error 2

The error can be ignored since it is a `bug in the test itself <http://sourceforge.net/p/icu/mailman/message/32443311/>`__.

Module Files
------------

Version 58.1
^^^^^^^^^^^^

Module File Location: ``/usr/local/modulefiles/libs/gcc/4.9.2/icu/58.1`` ::

        #%Module1.0#####################################################################
        ##
        ## icu 58.1 module file
        ##

        ## Module file logging
        source /usr/local/etc/module_logging.tcl
        ##

        set vers 58.1
        set gccvers 4.9.2

        proc ModulesHelp { } {
            global vers
            global gccvers
            puts stderr "Makes icu library $vers (and GCC $gccvers) available"
        }
        module-whatis "Makes icu library $vers (and GCC $gccvers) available"

        # Run-time dependency on C++ std lib
        module load compilers/gcc/$gccvers

        set ICU_DIR /usr/local/packages6/libs/gcc/$gccvers/icu/$vers

        prepend-path LD_LIBRARY_PATH $ICU_DIR/lib
        prepend-path LIBRARY_PATH $ICU_DIR/lib
        prepend-path CPATH $ICU_DIR/include


Version 55
^^^^^^^^^^

Module File Location: ``/usr/local/modulefiles/libs/gcc/4.8.2/icu/55`` ::

        #%Module1.0#####################################################################
        ##
        ## icu 55 module file
        ##

        ## Module file logging
        source /usr/local/etc/module_logging.tcl
        ##

        proc ModulesHelp { } {
                puts stderr "Makes the icu library available"
        }

        set ICU_DIR /usr/local/packages6/libs/gcc/4.8.2/icu/55

        module-whatis   "Makes the icu library available"

        prepend-path LD_LIBRARY_PATH $ICU_DIR/lib
        prepend-path LIBRARY_PATH $ICU_DIR/lib
        prepend-path CPLUS_INCLUDE_PATH $ICU_DIR/include

Version 42
^^^^^^^^^^

Module File Location: ``/usr/local/modulefiles/libs/gcc/4.4.7/icu/42`` ::

        #%Module1.0#####################################################################
        ##
        ## icu 42 module file
        ##

        ## Module file logging
        source /usr/local/etc/module_logging.tcl
        ##

        proc ModulesHelp { } {
                puts stderr "Makes the icu library available"
        }

        set ICU_DIR /usr/local/packages6/libs/gcc/4.4.7/icu/42

        module-whatis   "Makes the icu library available"

        prepend-path LD_LIBRARY_PATH $ICU_DIR/lib
        prepend-path LIBRARY_PATH $ICU_DIR/lib
        prepend-path CPLUS_INCLUDE_PATH $ICU_DIR/include
