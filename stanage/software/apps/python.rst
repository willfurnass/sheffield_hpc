.. _python_stanage:
.. _anaconda3_stanage:

Python
======

.. sidebar:: Python

    :Versions: 2019.7,2020.11,2021.11,2022.05,2022.10
    :Documentation: https://www.python.org/doc/
    :URL: https://python.org


This page documents the "anaconda" installation on Stanage. This is the
recommended way of using Python, and the best way to be able to configure custom
sets of packages for your use.

"conda" a Python package manager, allows you to create "environments" which are
sets of packages that you can modify. It does this by installing them in your
home area. This page will guide you through loading conda and then creating and
modifying environments so you can install and use whatever Python packages you
need.

Using Conda Python
------------------

.. attention::
         
        **We recommend that you use the following 2022 (sub)version of Anaconda3:** ``Anaconda3/2022.05``
        
        *The latest module,* ``Anaconda3/2022.10`` *, is under investigation as this has demonstrated odd 
        behaviour on conda environment exit (for some users). We will investigate this, and advise in due course.*

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

Anaconda Python can be loaded with one of the following:

.. include:: /referenceinfo/imports/stanage/packages/anaconda3-ml-el7-icelake-znver-stanage.rst

The ``root`` conda environment (the default) provides Python 3 and no extra
modules, it is automatically updated, and not recommended for general use, just
as a base for your own environments.

.. warning::

    Due to Anaconda being installed in a module you must use the ``source`` command instead of ``conda`` 
    when activating or deactivating environments!


Creating a Conda Environment
----------------------------

Every user can create their own environments, and packages shared with the
system-wide environments will not be reinstalled or copied to your file store,
they will be *symlinked*, this reduces the space you need in your ``/users/$USER``
directory to install many different Python environments.

To create a clean environment with just Python 3.8 and numpy you can run::

    conda create -n mynumpy python=3.8 numpy

This will download the latest release of Python 3.8 and numpy, and create an
environment named ``mynumpy``.

Any version of Python or list of packages can be installed::

    conda create -n myscience python=3.5 numpy=1.15.2 scipy

If you wish to modify an existing environment, such as one of the anaconda
installations, you can ``clone`` that environment::

    conda create --clone myscience -n myexperiment

This will create an environment called ``myexperiment`` which has all the
same conda packages as the ``myscience`` environment.

How to avoid large conda environments filling up your home directory
--------------------------------------------------------------------

.. include:: /referenceinfo/imports/software/python/conda_in_fastdata.rst

Installing Packages Inside a Conda Environment
----------------------------------------------

Once you have created your own environment you can install additional packages
or different versions of packages into it. There are two methods for doing
this, ``conda`` and ``pip``, if a package is available through conda it is
strongly recommended that you use conda to install packages. You can search for
packages using conda::

    conda search pandas

then install the package using::

    conda install pandas

if you are not in your environment you will get a permission denied error
when trying to install packages, if this happens, create or activate an
environment you own.

If a package is not available through conda you can search for and install it
using pip, *i.e.*::

    pip search colormath
    pip install colormath


Using conda Environments
------------------------

Once the conda module is loaded you have to load or create the desired
conda environments. For the documentation on conda environments see
`the conda documentation <http://conda.pydata.org/docs/using/envs.html>`_.

You can load a conda environment with::

    source activate myexperiment

where ``myexperiment`` is the name of the environment, and unload one with::

    source deactivate

which will return you to the ``root`` environment.

It is possible to list all the available environments with::

    conda env list

Provided system-wide are a set of anaconda environments, these will be
installed with the anaconda version number in the environment name, and never
modified. They will therefore provide a static base for derivative environments
or for using directly.

Using Conda and Python in a batch job
-------------------------------------

Create a batch job submission script called ``myscript.slurm`` that is similar to the following:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --ntasks=1
   #SBATCH --time=10:00
   #SBATCH --mem-per-cpu=100

   export SLURM_EXPORT_ENV=ALL
   module load Anaconda3/2022.10

   # We assume that the conda environment 'myexperiment' has already been created
   source activate myexperiment
   srun python mywork.py

Then submit this to Slurm by running:

.. code-block:: bash

   sbatch myscript.slurm

Further Conda Python Learning Resources
---------------------------------------

.. include:: /referenceinfo/imports/software/python/python_learning_resources_import.rst

------------------

.. _venv_stanage: 

Python Virtual Environments
---------------------------

Please note :ref:`conda <anaconda3_stanage>` is the recommended way of using Python, and the best way to be able to configure custom
sets of packages for your use, allowing greater flexibility and ease when handling complex dependencies.
However, for simpler **Python-only** projects, venv can be a lightweight and effective alternative. 

Why Use venv with Python?
-------------------------

1. **Optimised Performance**
   - Python modules available via the module system are compiled specifically for Stanage's hardware,
   which means they are optimised for the CPU microarchitecture, leading to improved performance for compute-intensive tasks.

2. **Seamless Integration with Other module system Software**
   - By using ``--system-site-packages``, your virtual environment can inherit system-level packages like ``wheel``.   

3. **Lightweight and Fast**:
   - ``venv`` environments are lightweight compared to Conda environments, as they do not need to install a separate Python distribution.
   This means that the virtual environment only needs to store the user-specific packages and symlinks to the system-level Python libraries, reducing the amount of disk space required.

When Should You Use venv Instead of Conda?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- **Performance-critical applications**: When you want to take advantage of the CPU-specific optimisations.
- **Integration with other module system software**: Particularly when working with tightly coupled HPC tools like MPI.
- **Minimalist environments**: When you prefer lightweight management and the overhead of managing a separate Conda distribution is unnecessary.

Conda offers a wider selection of precompiled packages, but for environments tailored to the cluster hardware,
``venv`` is lighter and avoids some of the conflicts that can arise with Conda's own Python distribution.

This guide explains how to set up Python environments using ``venv`` as an alternative to Conda,
specifically for users who want to leverage the Python interpreters optimised for the HPC environment.
By following this guide, you will learn how to create, activate, and manage a virtual environment that integrates seamlessly with other module system built software.

.. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

Python can be loaded with one of the following:

.. include:: /referenceinfo/imports/stanage/packages/python-ml-el7-icelake-znver-stanage.rst

This command loads a version of Python, built using a version of GCCcore.
Make sure to adjust the module name to match GCC/GCCcore version of any other modules you need to load.

.. _sys-site-pkgs-note:

.. note:: 

   **System-site packages for:**
    - **Python-bare** *setuptools, pip* 
    - **Python** *flit_core, packaging, pip, setuptools, setuptools-scm, tomli, typing_extensions, wheel*
   
Create a Virtual Environment with venv
--------------------------------------

Now that the optimised Python interpreter is loaded, you can create a virtual environment using ``venv``.

To keep your virtual environments organised, you can set a variable in your ``.bashrc`` to define a consistent location:

.. _venv_home:

.. code-block:: bash

   echo export VENV_HOME=$HOME/.venvs >> ~/.bashrc
   source ~/.bashrc

.. note::

   For larger environments or for heavy users - due to limited home directory storage we recommend building your enviroments in parscratch:

   .. code-block:: bash

      mkdir -p /mnt/parscratch/users/$USER/
      echo export VENV_HOME=/mnt/parscratch/users/$USER/.venvs >> ~/.bashrc
      source ~/.bashrc

Let's create a new ``venv`` which will allow you to install user-specific packages while still inheriting system-wide packages.

.. code-block:: bash

   python -m venv --system-site-packages $VENV_HOME/my_venv

- ``--system-site-packages`` This flag allows the virtual environment to use Python module :ref:`system site packages<sys-site-pkgs-note>`.
- ``my_venv`` This is the name of your virtual environment. You can replace it with any name you prefer.

After creating the virtual environment, activate it to begin using it:

.. code-block:: bash

   source $VENV_HOME/my_venv/bin/activate

Once activated, you will see ``(my_venv)`` appear at the beginning of your command prompt, indicating that you are now working within the virtual environment.

.. code-block:: console

   (my_venv) [te1st@node001 [stanage] ~]$

Within the virtual environment, you can install any additional Python packages you need. For instance, to install ``numpy`` and ``pandas``:

.. code-block:: bash

   pip install numpy pandas

These packages will be installed locally in your virtual environment without affecting the system-wide packages. However, we could instead load a ``SciPy-bundle`` module, which adds its packages to the ``PYTHONPATH``, making ``numpy`` and ``pandas`` importable (refer to :ref:`avail_sys_site_pkgs_stanage`).

To list your virtual environments::

        ls $VENV_HOME

To exit your virtual environment once you’re done, simply run the following command::

        deactivate


Example Using Module System Packages with venv
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your project requires additional software that is available as a module, you can load those modules after creating or activating your virtual environment, such as:

.. code-block:: bash

   module load mpi4py/3.1.4-gompi-2022b

This allows your virtual environment to access additional module system tools or libraries.
The combination of Python packages in ``venv`` and modules like ``mpi4py`` ensures you have the necessary tools configured for optimal performance and compatibility.

Suppose you are working on a project that requires ``mpi4py`` for parallel computing.
You can use the modular system loaded Python interpreter along with a ``venv`` environment that inherits the optimised MPI configuration:

*Since we are demonstrating two approaches — one straightforward and another slightly safer by ensuring Python version consistency — you can choose the option (tab) that best suits your project's needs.*

.. tabs::

    .. group-tab:: Basic

       .. code-block:: bash
       
          module load Python/3.10.8-GCCcore-12.2.0
          python -m venv --system-site-packages $VENV_HOME/my_venv
          source $VENV_HOME/my_venv/bin/activate
          module load mpi4py/3.1.4-gompi-2022b


    .. group-tab:: Safe Option

       Here we append the environment variable ``$EBVERSIONPYTHON`` to the environment name, ensuring that we do not mix python versions.

       .. code-block:: bash
       
          MY_ENV="my_venv"
          module load Python/3.10.8-GCCcore-12.2.0
          python -m venv --system-site-packages ${VENV_HOME}/${MY_ENV}_${EBVERSIONPYTHON}
          source ${VENV_HOME}/${MY_ENV}_${EBVERSIONPYTHON}/bin/activate
          module load mpi4py/3.1.4-gompi-2022b

In this example, ``mpi4py`` is already available from the module loaded system installation, enabling you to ``import mpi4py`` directly in your Python script. The ``OpenMPI`` module is also loaded to ensure proper MPI functionality.

   
**Using this Python venv in a batch job:**

Create a batch job submission script called ``batch.sh`` that is similar to the following:

.. tabs::

    .. group-tab:: Basic

       .. code-block:: bash
       
          #!/bin/bash
          #SBATCH --ntasks-per-node=4
          #SBATCH --time=10:00
          #SBATCH --mem=4000
       
          export SLURM_EXPORT_ENV=ALL
          module load Python/3.10.8-GCCcore-12.2.0
          module load mpi4py/3.1.4-gompi-2022b
       
          # Use the virtual environment 'my_env' you already created
          source $VENV_HOME/my_venv/bin/activate

          # Run the Python script using srun, which can leverage MPI
          srun python my_mpi_work.py 

    .. group-tab:: Safe Option

       Here we append the environment variable ``$EBVERSIONPYTHON`` to the environment name, ensuring that we do not mix python versions.

       .. code-block:: bash
       
          #!/bin/bash
          #SBATCH --ntasks-per-node=4
          #SBATCH --time=10:00
          #SBATCH --mem=4000
       
          export SLURM_EXPORT_ENV=ALL
          module load Python/3.10.8-GCCcore-12.2.0
          module load mpi4py/3.1.4-gompi-2022b
       
          # Use the virtual environment 'my_env' you already created
          MY_ENV="my_venv"
          source ${VENV_HOME}/${MY_ENV}_${EBVERSIONPYTHON}/bin/activate

          # Run the Python script using srun, which can leverage MPI
          srun python my_mpi_work.py 

Where ``import mpi4py`` is present in ``my_mpi_work.py``, and we use ``srun`` so that tasks are properly distributed when using MPI.

Then submit this to Slurm by running:

.. code-block:: bash

   sbatch batch.sh

**Key Points:**
 - Use ``$VENV_HOME`` for consistent environment location.
 - Append ``$EBVERSIONPYTHON`` to environment names to avoid mixing Pyhton versions,
   whilst also setting ``$MY_ENV``.

.. _avail_sys_site_pkgs_stanage:

Available System-site packages
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Examples of some module system packages that have python modules which can be imported are:

.. list-table::
   :widths: 30 30 30 30

   * - AmberTools
     - GitPython
     - mpi4py
     - PyYAML
   * - archspec
     - gmpy2
     - networkx
     - scikit-build
   * - astropy
     - h5py
     - ParaView
     - snakemake
   * - Biopython
     - hypothesis
     - Pillow
     - sympy
   * - CDFlib
     - IPython
     - pkgconfig
     - TensorFlow
   * - cppy
     - jax
     - poetry
     - VTK
   * - dill
     - Mako
     - protobuf-python
     - Xarray
   * - flatbuffers
     - matplotlib
     - pybind11
     - yaff
   * - GDAL
     - Meson
     - pytest-xdist
     - 

With a given module loaded you can list available python modules which we can import with the command: 

.. code-block:: console

        ls $(echo $PYTHONPATH | awk -F: '{print $1}')

For example:

.. code-block:: console

        [te1st@node001 [stanage] ~]$ module load SciPy-bundle/2023.02-gfbf-2022b
        [te1st@node001 [stanage] ~]$ ls $(echo $PYTHONPATH | awk -F: '{print $1}')
        beniget  bottleneck  deap  gast  mpmath  numexpr  numpy  omp  pandas  ply  pythran  scipy 

Here, for clairty we have removed files such as ``*.dist-info`` from the terminal output.

.. attention::

   Always load a package module that matches the toolchain (i.e., GCC version) as that of the Python module you are using.

------------------

Installing a personal copy of miniconda
---------------------------------------


.. dropdown:: We recommend using Anaconda3, however if you would like to have a personal install of miniconda please follow this guide
   

        Miniconda is a free, minimal installer for Conda. It is a small bootstrap version of Anaconda that includes only Conda,
        Python, their dependencies, and a small number of other useful packages such as pip and zlib.

        .. include:: /referenceinfo/imports/scheduler/SLURM/common_commands/srun_start_interactive_session_import_stanage.rst

        The latest miniconda releases will be made available at https://docs.conda.io/en/latest/miniconda.html

        On Stanage you should install miniconda to your home directory, e.g., ``$HOME/software/installs/miniconda``.

        Download the installer: ::

            wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

        Check the sums match with those listed on the `downloads page <https://repo.anaconda.com/miniconda/>`_ : ::

            sha256sum Miniconda3-latest-Linux-x86_64.sh

        Make the file executable: ::

            chmod +x Miniconda3-latest-Linux-x86_64.sh

        Create a directory for the installation::

             mkdir -p $HOME/software/installs

        Run the installer and choose the directory created above (``$HOME/software/installs``) as the installation target. 

        .. code-block:: bash

            ./Miniconda3-latest-Linux-x86_64.sh

        After accepting the licence agreement enter the path for the installation: ::

            $HOME/software/installs/miniconda
            
        .. caution:: 

             Ensure you **do not initialise Miniconda** as this will interfere with loading other Anaconda modules on the cluster.

        Make a modules directory for yourself: ::

            mkdir $HOME/modules
            
        Make a module file for the install you just made: ::

            nano $HOME/modules/miniconda.lua
            
        With content as follows:

        .. code-block:: lua
           :force:

                    ------------------------------------------------------------------------------------------------ 
                    -- ~/modules/miniconda.lua:
                    ------------------------------------------------------------------------------------------------ 
                    
                    help([[
                    Description
                    ===========
                    Makes my personal Miniconda install available.

                    More information
                    ================
                    - Homepage: https://www.anaconda.org
                    ]])

                    -- Describe the module.
                    whatis("Description: Makes my personal Miniconda install available.")
                    whatis("Homepage: https://docs.anaconda.com/miniconda")
                    whatis("URL: https://www.anaconda.org")

                    conflict("Anaconda3")

                    local HOME = os.getenv("HOME")
                    local MINICONDA_DIR = HOME .. "/software/installs/miniconda"
                    setenv("MINICONDA_DIR", MINICONDA_DIR)

                    -- Add directories to environment variables.
                    prepend_path("PATH", pathJoin(MINICONDA_DIR, "bin"))
            
        Note: We have added Anaconda3 as a conflicting module to prevent Miniconda from being loaded concurrently.
        This partially avoids conflicts, as loading Anaconda3 after Miniconda would still result in conflicts.

        .. warning::

            Module files are written in LUA not bash. See :ref:`Making software available via a custom module file<custom-module-files>`.
            
        Make your custom modules available: ::

            module use $HOME/modules

        To skip doing the above everytime you login, you can add this line to the .bashrc file in your home directory::

            echo "module use $HOME/modules" >> ~/.bashrc

        Assuming the module file is named ``miniconda``, you can load it with: ::

            module load miniconda

        .. warning::
        
            Since Miniconda is installed as a module, use the ``source`` command instead of ``conda`` to activate or deactivate environments.

        To activate the base environment:

        .. code-block:: bash
           
            source activate

--------------

Installation notes
------------------

Anaconda3
^^^^^^^^^
Anaconda was installed using Easybuild, build details can be found in folder ``$EBROOTANACONDA3/easybuild`` with a given module loaded.

Python
^^^^^^
Python was installed using Easybuild, build details can be found in folder ``$EBROOTPYTHON/easybuild`` with a given module loaded.

.. include:: /referenceinfo/imports/stanage/packages/python-dpnd-el7-icelake-znver-stanage.rst
