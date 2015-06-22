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
This build of the Boost library requires gcc version 4.8.2. To make the compiler and library available, run the following module commands ::

    module load compilers/gcc/4.8.2
    module load libs/gcc/4.8.2/boost/1.58

Build a simple program using Boost
----------------------------------

Many boost libraries are header-only which makes them particularly simple to compile. The following program reads a sequence of integers from standard input, uses Boost.Lambda to multiply each number by three, and writes them to standard output (taken from http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html):

.. code-block:: c++

        #include <boost/lambda/lambda.hpp>
        #include <iostream>
        #include <iterator>
        #include <algorithm>

        int main()
        {
            using namespace boost::lambda;
            typedef std::istream_iterator<int> in;

            std::for_each(
                in(std::cin), in(), std::cout << (_1 * 3) << " " );
        }

Copy this into a file called example1.cpp and compile with

:code:`g++ example1.cpp -o example`

Provided you loaded the modules given above, and you are using gcc version 4.8.2, the program should compile without error.

Linking to a Boost library
--------------------------

The following program is taken from the official Boost documentation http://www.boost.org/doc/libs/1_58_0/more/getting_started/unix-variants.html

.. code-block:: c++

        #include <boost/regex.hpp>
        #include <iostream>
        #include <string>

        int main()
        {
            std::string line;
            boost::regex pat( "^Subject: (Re: |Aw: )*(.*)" );

            while (std::cin)
            {
                std::getline(std::cin, line);
                boost::smatch matches;
                if (boost::regex_match(line, matches, pat))
                    std::cout << matches[2] << std::endl;
            }
        }

This program makes use of the Boost.Regex library, which has a separately-compiled binary component we need to link to.
Assuming that the above program is called example2.cpp, compile with the following command

:code:`g++ example2.cpp -o example2 -lboost_regex`

If you get an error message that looks like this:

:code:`example2.cpp:1:27: error: boost/regex.hpp: No such file or directory`

the most likely cause is that you forgot to load the module as detailed above.

Installation Notes
------------------
This section is primarily for administrators of the system.

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
Compiled and ran the two example files given above.

Module File
-----------
Module File Location: :code:`/usr/local/modulefiles/libs/gcc/4.8.2/boost/1.58`

.. code-block:: none

        #%Module1.0#####################################################################
        ##
        ## boost 1.58 module file
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


