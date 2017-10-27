.. _jh_conda: 

Jupyter on SHARC: languages, packages and environments
======================================================

.. toc:

Background
----------

Once you have a Jupyter Notebook server running 
(e.g. on a cluster worker node)
you typically want to consider your **execution environment**,
which is your choice of 

* **kernel** (code cell language)
* **software packages/libraries** 

Kernels
^^^^^^^

Jupyter Notebooks were originally known as IPython Notebooks
but now Jupyter supports Notebook code cells written in `many different languages <https://github.com/jupyter/jupyter/wiki/Jupyter-kernels>`__.
This is achieved by the (per-user) Jupyter Notebook server 
sending the contents of code cells to a **kernel** for evaluation.
The most widely used kernels are:

* the `IPython kernel`_, 
  the default kernel, 
  which can run code cells containing Python >= 3.3 and Python 2.7;
* IRKernel_, a kernel for R_ >= 3.2.

Packages
^^^^^^^^

Most notebooks make use of external software packages for e.g. fitting statistical models to data.
There are several different ways you might enable/load packages on ShARC 
(including :ref:`module files <env_modules>` and :ref:`Singularity containers <singularity_sharc>`) 
but if using Jupyter it is recommended that you install and activate software using 
the conda_ package manager if possible.  
This can be done using Jupyter's graphical interface (in most cases)
or from the command-line.

Environments
^^^^^^^^^^^^

conda_ packages are installed into **environments**.  
An environment is an isolated collection of packages.  
You might create:

* One environment for one project containing Python 3, the IPython kernel and the ``pandas`` and ``matplotlib`` Python packages (plus dependencies) for data analysis.
* Another environment for a second project containing R, the IRKernel and the ``dplyr`` and ``gplot2`` R packages.
* A third environment containing Python 2.7 plus ``pandas`` but no kernel for work that doesn't involve Jupyter Notebooks.

conda_ allows users to 

* Install and manage packages without a system adminstrator needing to be involved;
* Isolate and audit the set of packages used per project (good for reproducible research)
* Share environment definitions with others and with automated test/build systems (i.e. `continuous integration`_)

Using conda on ShARC via Jupyter
--------------------------------

From the browser tab containing the Jupyter's :ref:`file browser <jh_file_browse>` 
activate the **Conda** tab within Jupyter's interface.

You should then see something like:

.. image:: /images/jupyterhub/jupyter-conda-view.png

Here we have (in anticlockwise order) lists of:

#. All **conda environments** that have been found and can be *activated* from within Jupyter;
#. **Packages** (inc. versions) that **can be installed** into a selected environment; 
   by default only the latest version from the ``default`` conda *channel* (package repository) can be installed this way.
#. **Packages** that **are installed** in the environment selected in the above pane 
   (plus the package version and the build version (e.g. version of Python it was built for)).


These three views may take a few seconds to populate after clicking **Conda**.


CLICK THE CONDA TAB




Before we can run a Notebook we typically want to 
create a new conda_ **environment** containing:
 

**TODO**

NOTES BELOW
-----------

SOME NOTES
----------

* LIST SHOWS KERNEL NAME, TYPE, ENV NAME
* DISCOURAGED FROM USING OTHER ENVS
* The **anaconda** environments have been installed by a system administrator and cannot be changed.
* They give you access to a large number of Python packages commonly used by data scientists.
* The **jupyterhub** environments are for administrative purposes and should never be used.
* CENTRALLY MANAGED OR NOT
* NOT A CONDA GUIDE
* WHERE TO GO FOR CONDA HELP
* ANACONDA - HOW TO INSTALL


SWITCHING KERNEL

.. image:: /images/jupyterhub/sharc-jh-change-kernel-in-nb.png

Managing (Conda) environments
-----------------------------

**TODO**
WHAT CAN YOU DO
WHY WANT TO DO IT
LIMITATIONS - CAN'T MUTATE ANYTHING!

.. image:: /images/jupyterhub/sharc-jh-pkg-mgr.png


Choosing your programming language and packages
-----------------------------------------------

Before we can run a Notebook we typically want to 
create a new `Conda <https://conda.io/docs/using/envs.html>`_ **environment** 
containing:

* a set of **packages** (e.g. particular packages for creating plots);
* a Jupyter **kernel** that runs the code snippets in your Notebook.  
  This will typically be the **Python** or **R** kernel.  

The **Conda** tab *should* allow you to create a new environment containing specific packages
but this is not working at present; instead we can create a new environment by: 

#. From the **Files** tab click **New** then **Terminal** to bring up a command-line prompt in a new browser tab;
#. Create a new environment containing the Python kernel (``ipykernel``), a specific version of Python and the latest compatible versions of the ``numpy`` and ``matplotlib`` packages: ::

        conda create -n my-new-env python=3.6 numpy matplotlib ipykernel

#. Press ``y`` when asked if you want to download and install those packages and their dependencies;
#. A new environment called ``my-new-env`` will then be created
   (all files relating to it are in the ``.conda/envs/my-new-env`` directory within your home directory);
#. Save the state of your environment (i.e. all versions of installed packages) to a file (for posterity) using: ::

        conda env export -n my-new-env > /home/YOURUSERNAME/my-new-env.yml

#. Type ``exit`` then press Enter in your terminal browser tab then close the tab.

.. _conda: https://conda.io/docs/using/envs.html
.. _IPython kernel: https://github.com/ipython/ipykernel
.. _IRKernel: https://irkernel.github.io/
.. _R: https://www.r-project.org/
.. _continuous integration: https://en.wikipedia.org/wiki/Continuous_integration
