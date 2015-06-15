Boost C++ Library v1.41
=======================

.. sidebar:: Boost version 1.41
   
   :Version: 1.41
   :Support Level: Bronze
   :Dependancies:  libs/gcc/4.4.7/icu/42
   :URL: www.boost.org
   :Documentation: http://www.boost.org/doc/libs/1_41_0/
   :Location: /usr/local/packages6/libs/gcc/4.4.7/boost/1.41

Boost provides free peer-reviewed portable C++ source libraries.

Usage
-----
This build of the Boost library was built with gcc 4.4.7 and so should only be used with that version of gcc. To make the library available, run the following module command.

:code:`module load libs/gcc/4.4.7/boost/1.41`

Build a simple program using Boost
----------------------------------

Many boost libraries are header-only which makes them particularly simple to compile. The following program reads a sequence of integers from standard input, uses Boost.Lambda to multiply each number by three, and writes them to standard output (taken from http://www.boost.org/doc/libs/1_41_0/more/getting_started/unix-variants.html):

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

Provided you loaded the module given above, and you are using gcc version 4.4.7, the program should compile without error.

Linking to a Boost library
--------------------------
The following program is taken from the official Boost dcoumentation http://www.boost.org/doc/libs/1_41_0/more/getting_started/unix-variants.html

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
This build of boost was built with gcc 4.4.7 and ICU version 42.

.. code-block:: none
        
        module load libs/gcc/4.4.7/icu/42
        tar -xvzf ./boost_1_41_0.tar.gz 
        cd boost_1_41_0
        ./bootstrap.sh --prefix=/usr/local/packages6/libs/gcc/4.4.7/boost/1.41
        ./bjam -sICU_PATH=/usr/local/packages6/libs/gcc/4.4.7/icu/42 install
    

Testing
-------
The two examples above were compiled and ran.

Module File
-----------
Module File Location: :code:`/usr/local/modulefiles/libs/gcc/4.4.7/boost/1.41`

.. code-block:: none

        #%Module1.0#####################################################################
        ##
        ## Boost 1.41 module file
        ##

        ## Module file logging
        source /usr/local/etc/module_logging.tcl
        ##

        module load libs/gcc/4.4.7/icu/42

        proc ModulesHelp { } {
                puts stderr "Makes the Boost 1.41 library available"
        }

        set BOOST_DIR /usr/local/packages6/libs/gcc/4.4.7/boost/1.41

        module-whatis   "Makes the Boost 1.41 library available"

        prepend-path LD_LIBRARY_PATH $BOOST_DIR/lib
        prepend-path CPLUS_INCLUDE_PATH $BOOST_DIR/include
        prepend-path LIBRARY_PATH $BOOST_DIR/lib

