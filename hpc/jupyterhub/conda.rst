.. _jh_conda: 

Jupyter on SHARC: languages, packages and environments
======================================================

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
#. Save the state of your environment (i.e. all versions of installed packages) to a file (for posterity) using:

        conda env export -n my-new-env > /home/YOURUSERNAME/my-new-env.yml

#. Type ``exit`` then press Enter in your terminal browser tab then close the tab.
