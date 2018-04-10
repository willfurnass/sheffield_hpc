.. _jh_nb_usage:

Creating, editing and running Jupyter Notebooks
===============================================

Creating Notebooks
------------------

After :ref:`creating or deciding on a conda environment <jh_conda>` 
containing the Jupyter Kernel(s) you want to execute Notebook(s) with plus the packages want to import within Notebook(s)
you can now create a Notebook or open an existing one.

To create a Notebook:

#. Return to the Jupyter *Home* browser tab; 
#. Click the *Files* Jupyter tab;
#. Browse to the directory where you want to create your new Notebook;
#. Click **New** then (beneath **Notebooks**) the name of the Kernel/environment you wish to use 
   (e.g. ``rdkit-sharc``) - see `Selecting a Jupyter Kernel`_ for more information on selecting kernels.

   .. image:: /images/jupyterhub/sharc-jh-new-nb-w-kernel.png

#. A blank Notebook should appear in a new browser tab.

Your Notebook will have access to the packages installed in the selected environment.

Opening existing Notebooks
--------------------------

Alternatively you can click on an existing Notebook (``.ipynb``) file in Jupyter's file browser to open it.

.. image:: /images/jupyterhub/sharc-jh-example-nb.png

Selecting a Jupyter Kernel
--------------------------

After opening a Notebook, you can **change the Kernel used for executing code cells** by 
clicking *Kernel* -> *Select Kernel* from the menu bar to 
bring up a list of availble Kernels.

Some of the Kernels in this list correspond to conda environments created by the system administrator; 
others were automatically found by a Jupyter plug-in that 
searches for valid Jupyter Kernels in all conda environments visible to you.

It is **recommended that you create your own environments** (typically one per project/workflow).

**Do not use** the ``jupyterhub`` or ``jupyterhub-dev`` environments.
You are advised not to use the ``anaconda`` Kernels/environments either as these are read-only to most users
and users have little control over if/when they are updated and what packages they contain.  

Using Jupyter Notebooks
-----------------------

The basics of using Jupyter Notebooks to create self-describing, runable workflow documents 
are explained in the `Jupyter Notebook official documentation`_.

.. _Jupyter Notebook official documentation: http://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/Notebook%20Basics.html
