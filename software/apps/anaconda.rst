Anaconda Python
===============

.. sidebar:: Anaconda Python

   :Support Level: gold
   :Dependancies: None
   :URL: https://store.continuum.io/cshop/anaconda/
   :Version: 2.3 (Providing Python 2.7.10)

Anaconda Python is a Python distribution that contains Python itself and a large number of popular modules.

Interactive Usage
-----------------
After connecting to iceberg (see :ref:`ssh`),  start an interactive sesssion with the :code:`qsh` or `qrsh` command.

The latest version of Anaconda can be loaded with ::

        module load apps/binapps/anacondapython

Alternatively, you can load a specific version of Anaconda Python using ::

        module load apps/binapps/anacondapython/2.3

Python can then be started with the ``python`` command::

        $ python

        Python 2.7.10 |Anaconda 2.3.0 (64-bit)| (default, May 28 2015, 17:02:03)
        [GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux2
        Type "help", "copyright", "credits" or "license" for more information.
        Anaconda is brought to you by Continuum Analytics.
        Please check out: http://continuum.io/thanks and https://binstar.org
        >>>

Available Python Modules
------------------------
Anaconda Python provides a large number of Python modules including numpy, scipy, matplotlib, scikit-learn, scikit-image, ipython and many more. In addition to these standard modules, we have installed others following user requests. To see which modules are available, run the following command in an interactive Python session ::

     help('modules')

Installation Notes
------------------
These are primarily for administrators of the system. Anaconda Python is a binary install of Python. ::

  mkdir -p /usr/local/packages6/apps/binapps/anacondapython/2.3
  chmod +x ./Anaconda-2.3.0-Linux-x86_64.sh

  #Run the installer
  ./Anaconda-2.3.0-Linux-x86_64.sh

When the installer asks for location, use ::

  /usr/local/packages6/apps/binapps/anacondapython/2.3

Some modules were requested by users that were not included in the standard Anaconda installer. I installed these by doing ::

  module load apps/binapps/anacondapython/2.3
  conda install wxPython

Module file
-----------
location ``/usr/local/modulefiles/apps/binapps/anacondapython/2.3`` ::

  #%Module10.2####################################################################
  #

  ## Module file logging
  source /usr/local/etc/module_logging.tcl
  ##

  proc ModulesHelp { } {
          global ver

          puts stderr "   Adds Anaconda $ver to your environment variables."
  }

  # Anaconda version (not in the user's environment)
  set     ver     2.3

  module-whatis   "sets the necessary Anaconda $ver paths"

  set ANACONDADIR /usr/local/packages6/apps/binapps/anacondapython/
  set ANACONDAHOME $ANACONDADIR/$ver

  setenv ANACONDAHOME $ANACONDAHOME
  prepend-path PATH $ANACONDAHOME/bin
