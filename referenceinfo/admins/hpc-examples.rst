.. _hpc-examples-admin:

HPC Example Scripts
===================

The hpc-examples_ repository is a public collection of useful scripts and workflows for HPC systems.

Integrating ``hpc-examples`` into ``sheffield_hpc``
---------------------------------------------------

The hpc-examples_ repository is included as a `Git submodule <https://git-scm.com/book/en/v2/Git-Tools-Submodules>`_ within sheffield_hpc.
Site maintainers can reference or include content from hpc-examples_ within the sheffield_hpc_
repository to enhance documentation and provide users with tested examples.

Steps to Use ``hpc-examples``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Clone or Update the Repository**

Clone the hpc-examples_ repository if it has not already been added, or pull the latest changes:
To ensure the hpc-examples submodule is included:

.. code-block:: bash

   git clone --recurse-submodules git@github.com:rcgsheffield/sheffield_hpc.git

If you've already cloned without submodules:

.. code-block:: bash

   git submodule update --init --recursive


2. **Reference Examples in Documentation**

Use Sphinx's ``literalinclude`` directive to embed scripts directly from the hpc-examples_ submodule. For example:

.. code-block:: bash

  .. literalinclude:: /hpc-examples/path/to/example_script.sh
     :language: bash
     :lines: 3-10

See more :ref:`examples <hpc-examples-example>`.

3. **Add download commands for individual files**

.. code-block:: bash
     
   wget https://raw.githubusercontent.com/rcgsheffield/hpc-examples/refs/heads/main/$PATH_TO_FILE

Where ``$PATH_TO_FILE`` is the relative path within hpc-examples_ , for example ``examples/slurm-ansys-2023R2.sh``

4. **Ensure Consistency**

- Review the examples in ``hpc-examples`` to ensure they align with the configurations and best practices of Sheffield’s HPC systems.
- If customisation is required, document the changes explicitly in the sheffield_hpc_ repository without modifying the original hpc-examples_ content.

-------------------

.. include:: /referenceinfo/imports/hpc-examples-submodule-guide.rst

Repository Structure
--------------------

**Directory Layout** - see hpc-examples_

- ``docs/``: Documentation files.
- ``examples/``: Example SLURM scripts and application workflows (referenced from sheffield_hpc_).
- ``tests/``: Example SLURM test scripts (referenced from sheffield_hpc_).
- ``tools/``: Utilities or helper scripts for HPC users.
- ``scripts/``: General-purpose scripts for users.
- ``configs/``: Example configuration files for applications or tools.
- ``tutorials/``: Scripts required for tutorials within our docmentation.

**Naming Conventions**

- Use descriptive, kebab-case names for files and directories.
- Include version numbers in filenames for version-specific examples.

.. _sheffield_hpc: https://github.com/rcgsheffield/sheffield_hpc
.. _hpc-examples: https://github.com/rcgsheffield/hpc-examples
