.. _jh_nb_usage:

Creating, editing and running Jupyter Notebooks
===============================================

Creating Notebooks
------------------

After :ref:`creating or deciding on a conda environment <jh_conda>` that you want to run a Notebook in 
you can now create a Notebook or open an existing one.

To create a Notebook:

#. Return to the Jupyter *Home* browser tab; 
#. Click the *Files* Jupyter tab;
#. Browse to the directory where you want to create your new Notebook;
#. Click **New** then (beneath **Notebooks**) the name of the conda environment you wish to use 
   (e.g. ``rdkit-sharc``).  Do **not use** the ``jupyterhub`` or ``jupyterhub-dev`` environments.
   You are advised not to use the ``anaconda`` environments either as these are read-only to most users
   and users have little control over if/when they are updated and what packages they contain.  
   It is recommended that you create your own environments.

   .. image:: /images/jupyterhub/sharc-jh-new-nb-w-kernel.png

#. A blank Notebook should appear in a new browser tab.

Your Notebook will have access to the packages installed in the selected environment.

Opening existing Notebooks
--------------------------

Alternatively you can click on an existing Notebook (``.ipynb``) file in Jupyter's file browser to open it.

.. image:: /images/jupyterhub/sharc-jh-example-nb.png

Using Jupyter Notebooks
-----------------------

The basics of using Jupyter Notebooks to create self-describing, runable workflow documents 
are explained in the `Jupyter Notebook official documentation`_.

.. _Jupyter Notebook official documentation: http://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/Notebook%20Basics.html
