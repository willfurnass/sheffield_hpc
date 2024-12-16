.. |softwarename| replace:: EnergyPlus
.. |currentver| replace:: 24.1.0
.. |ebtoolchain| replace:: foss2022b

.. _energyplus_stanage:

|softwarename|
==============

.. sidebar:: EnergyPlus

   :Dependencies: None
   :URL: https://energyplus.net/
   :Latest version: 24.1.0

|softwarename| is a whole building energy simulation program that engineers, architects, and researchers use to model energy and water use in buildings.

--------

Interactive usage
-----------------

The latest version of |softwarename| (currently version |currentver|) is made available with the command:

.. code-block:: console

	$ module load EnergyPlus/24.1.0-foss-2022b


After this any of the |softwarename| commands can be run from the terminal prompt. The available 
commands can be obtained using:

.. code-block:: console

	$ energyplus --help


Batch Submission
----------------


An example batch submission script for this file is :

.. code-block:: console

  #!/bin/bash
  # Request 4 gigabytes of real memory
  #SBATCH --mem=4G
  #SBATCH --mail-user=a.person@sheffield.ac.uk
  #SBATCH --mail-type=ALL
  
  module load EnergyPlus/24.1.0-foss-2022b
  
  # Set paths to input and output directories
  INPUT_FILE=/path/to/your/input_file.idf
  WEATHER_FILE=/path/to/weather_file.epw
  OUTPUT_DIR=/path/to/output_directory
  
  # Run EnergyPlus
  energyplus -w $WEATHER_FILE -r -d $OUTPUT_DIR -i $INPUT_FILE

The above is saved to a file called ``run_job.sh`` and submitted with :

.. code-block:: console

  sbatch run_job.sh

Installation notes
------------------

|softwarename| version 24.1.0 was installed using Easybuild 4.9.2, build details can be found 
in ``$EBROOTENERGYPLUS/easybuild`` with the module loaded.