.. _boost_stanage:

Boost C++ Library
==================

.. sidebar:: Boost C++ Library

   :Latest version: 1.81.0
   :URL: www.boost.org

Boost provides free, peer-reviewed and portable C++ source libraries.

Usage
-----

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

To activate the library run one of the following:
        
.. include:: /referenceinfo/imports/stanage/packages/boost-ml-el7-icelake-znver-stanage.rst

Boost has been built *without* Python support: use :ref:`conda <python_stanage>` and install the ``boost`` conda package if you want to use Boost with Python.

Build a simple program using Boost
----------------------------------

Many Boost libraries are header-only which makes them particularly simple to compile. The following program reads a sequence of integers from standard input, uses ``Boost.Lambda`` to multiply each number by three, and writes them to standard output (taken from http://www.boost.org/doc/libs/1_64_0/more/getting_started/unix-variants.html):

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

Copy this into a file called ``example1.cpp`` and compile with: ::

        g++ example1.cpp -o example

Provided you loaded the correct modules given above, the program should compile without error.

Linking to a Boost library
--------------------------
The following program is taken from the official Boost documentation http://www.boost.org/doc/libs/1_64_0/more/getting_started/unix-variants.html

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


This program makes use of the ``Boost.Regex`` library, which has a separately-compiled binary component we need to link to.
Assuming that the above program is called ``example2.cpp``, compile with the following command: ::

        g++ example2.cpp -o example2 -lboost_regex

Installation Notes
------------------

This section is primarily for administrators of the system. Boost has been installed using the default Easybuild config files.

Build logs and test reports can be found in ``$EBROOTBOOST/easybuild`` with a given module loaded.

Tested by compiling and running the two programs shown above.

.. include:: /referenceinfo/imports/stanage/packages/boost-dpnd-el7-icelake-znver-stanage.rst
