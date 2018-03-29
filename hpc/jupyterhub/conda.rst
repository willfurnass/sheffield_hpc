.. _jh_conda: 

Jupyter on SHARC: preparing your environment
============================================

.. toc:

Background
----------

Once you have a Jupyter Notebook server running 
(e.g. on a cluster worker node)
you typically want to consider your **execution environment**,
which is your choice of 

* **Kernel** (code cell language)
* **software packages/libraries** 

Kernels
^^^^^^^

Jupyter Notebooks were originally known as IPython Notebooks
but now Jupyter supports Notebook code cells written in `many different languages <https://github.com/jupyter/jupyter/wiki/Jupyter-kernels>`__.
This is achieved by the (per-user) Jupyter Notebook server 
sending the contents of code cells to a **Kernel** for evaluation.
The most widely used Kernels are:

* the `IPython Kernel`_, 
  the default Kernel, 
  which can run code cells containing Python >= 3.3 and Python 2.7;
* IRKernel_, a Kernel for R_ >= 3.2.

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

* One environment for one project containing Python 3, the IPython Kernel and the ``pandas`` and ``matplotlib`` Python packages (plus dependencies) for data analysis.
* Another environment for a second project containing R, the IRKernel and the ``dplyr`` and ``gplot2`` R packages.
* A third environment containing Python 2.7 plus ``pandas`` but no Kernel for work that doesn't involve Jupyter Notebooks.

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

.. warning::

   This interface currently does **not reliably allow conda environments to be created or modified**.
   Until this issue is resolved you can use this interface to inspect your environments 
   but you should create and modify conda environments from a :ref:`terminal <jh_terminal>`.

.. comment: 
   This Jupyter extension currently yields an unhelpful 'Error: not found' error in the UI
   

Creating a new conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before we run a Notebook we typically want to 
create a new conda_ **environment** containing
the packages we are interested in 
plus a Jupyter Kernel.

See the general documentation for :ref:`using conda on ShARC <sharc-python-conda>` for 
generic instructions on how to create conda environments from the command-line.  
Note that if you are using a :ref:`Jupyter Terminal <jh_terminal>` 
then you do not need to load conda using ``module load ...``.
Make sure you **install a package containing a Jupyter Kernel** (e.g. ``ipykernel`` for Python work) into your environment.

When following that documentation you might want to use the following as starting points for creating Jupyter(Hub)-compatible environments:

Python 3:

.. code-block:: sh

   conda create -n example-python-env python=3.6 ipykernel jupyter_client

R: 

.. code-block:: sh

   conda create -n example-r-env python=3.6 r-irkernel jupyter_client libiconf

Python from the `Intel Python Distribution <https://software.intel.com/en-us/distribution-for-python>`:

.. code-block:: sh

   conda create -n example-intel-python-env -c intel intelpython3_core ipykernel jupyter_client

.. comment:
   Omit the following until 
   1. finalised a conda env location policy that works on both clusters
   2. updated batch job submission script so CONDARC set 
      (as don't think .bashrc is read when start Jupyter session via batch job)
   3. Figured out how/why nb_conda does not allow new envs to be created.
   ---

   Before that, we need to configure conda so that it stores any environments we create in ``/data/username/.conda-sharc`` as
   
   * the set of conda packages used in our environments can grow quite large; ``/home`` can fill up quickly if environments are created there;
   * ``/home`` and ``/data`` are shared between ShARC and Iceberg but 
     sharing environments between the clusters causes problems; 
     we therefore need the environments to be stored in a per-cluster directory.
   
   First, tell conda where to look for a config file by starting a :ref:`Jupyter Terminal <jh_terminal>` then running ::
   
       echo "export CONDARC=$HOME/.condarc-${SGE_CLUSTER_NAME}.yml" >> ~/.bashrc
   
   Secondly, create a config file that tells conda where to look for / create conda environments.  Again, from a Jupyter Terminal: ::
   
       nano ~/.condarc-${SGE_CLUSTER_NAME}.yml
   
   to start editing the config file in text editor, then type: ::
   
       envs_dirs:
         - /data/te1st/.conda-sharc
   
   ...making sure to replace ``te1st`` with your username.
   
   To then create an environment using the *Conda* tab in the Jupyter user interface:
   
   #. Click the **+** above the upper most pane.
   #. Enter a name for your environment.
   #. Select a language/Kernel (currently only Python 2, Python 3 and R are supported)

Capturing the state of an environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is important to track the versions of packages you used to generate your research outputs, 
primarily to allow you and others to easily repeat your workflows.  
Ideally you should manage a file detailing the packages in your environment, 
plus your Notebook and other project files, 
in a version control system such as `git <https://en.wikipedia.org/wiki/Git>`__.

After creating/modifying an environment:

* Click on the export icon (left-most icon beneath **Action**) for a given environment
  to download a conda environment definition file.
* Alternatively you can generate a definition file 
  from a :ref:`Jupyter Terminal <jh_terminal>`: ::

    source activate my-env-name
    cd /data/$USER/research-project-7
    conda env export > environment.yml

.. comment out the following until modifying envs via nb_conda has been tested
   Modifying an environment
   ^^^^^^^^^^^^^^^^^^^^^^^^
   **TODO**
   * NB some are read-only

.. commented out the following notes until nb_conda has had more testing
   Managing (Conda) environments
   -----------------------------
   WHAT CAN YOU DO
   WHY WANT TO DO IT
   * DISCOURAGED FROM USING ANACONDA ENVS
   * SHOULD NOT USE JUPYTER ENVS 
   * The **anaconda** environments have been installed by a system administrator and cannot be changed.
   * They give you access to a large number of Python packages commonly used by data scientists.
   * The **jupyterhub** environments are for administrative purposes and should never be used.
   .. image:: /images/jupyterhub/sharc-jh-pkg-mgr.png

Next
----

After you have assessed what environments you have available, 
you can start :ref:`creating, editing and running Jupyter Notebooks <jh_nb_usage>`.

.. _conda: https://conda.io/docs/using/envs.html
.. _IPython Kernel: https://github.com/ipython/ipykernel
.. _IRKernel: https://irkernel.github.io/
.. _R: https://www.r-project.org/
.. _continuous integration: https://en.wikipedia.org/wiki/Continuous_integration
