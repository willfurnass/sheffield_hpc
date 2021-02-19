OpenFOAM
==========

.. sidebar:: OpenFOAM

   :Versions: 8.0
   :URL: https://openfoam.org/
   :Dependencies: Easybuild foss/2020a toolchain, NCurses 6.2, METIS 5.1.0, SCOTCH 6.0.9, CGAL 4.14.3 and Paraview 5.8.0

OpenFOAM is leading software for computational fluid dynamics (CFD). It is licensed free and open source only under the GNU General Public Licence (GPL) by the OpenFOAM Foundation.

Usage
-----
OpenFOAM can be used in an interactive or batch job. OpenFOAM 8.0 can be activated using the module file and sourcing the OpenFOAM environment script::

    module load OpenFOAM/8-foss-2020a
    source $FOAM_BASH

------------

Interactive Usage
--------------------

The following is an example interactive session running the pitzDaily example model.

After connecting to Bessemer (see :ref:`ssh`), you can start an `interactive graphical session <https://docs.hpc.shef.ac.uk/en/latest/hpc/scheduler/submit.html#interactive-sessions>`_. ::

    module load OpenFOAM/8-foss-2020a
    source $FOAM_BASH
    mkdir -p /fastdata/$USER/tests/openfoam/run
    rm -r /fastdata/$USER/tests/openfoam/run/*
    cd /fastdata/$USER/tests/openfoam/run
    cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily .
    chmod 700 -R pitzDaily && cd pitzDaily
    srun blockMesh
    srun simpleFoam
    srun paraFoam #To view the output.

------------

Batch Usage
--------------------

The following is an example batch job running the pitzDaily example model:

.. note::

    You will need to supply a `decomposeParDict <https://cfd.direct/openfoam/user-guide/v8-running-applications-parallel/>`_ in the system subdirectory of the case :

::

    #!/bin/bash
    #SBATCH --nodes=1
    #SBATCH --ntasks-per-node=4
    #SBATCH --mem=8000
    #SBATCH --job-name=name_OpenFOAM_smp_4
    #SBATCH --output=output_OpenFOAM_smp_4
    #SBATCH --time=01:00:00
    #SBATCH --mail-user=joe.bloggs@sheffield.ac.uk
    #SBATCH --mail-type=ALL
    mkdir -p /fastdata/$USER/tests/openfoam/run
    rm -r /fastdata/$USER/tests/openfoam/run/*
    cd /fastdata/$USER/tests/openfoam/run
    module load OpenFOAM/8-foss-2020a
    source $FOAM_BASH
    cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily .
    chmod 700 -R pitzDaily && cd pitzDaily
    srun blockMesh
    srun decomposePar
    srun simpleFoam -parallel

------------

Installation note for Administrators:
-------------------------------------



OpenFOAM 8 has been installed using Easybuild with all third party modules (NCurses 6.2, METIS 5.1.0, SCOTCH 6.0.9, CGAL 4.14.3 and Paraview 5.8.0)

Installation was tested as follows as above with the :download:`example batch script modified </bessemer/software/modulefiles/OpenFOAM/test_OpenFOAM_parallel.sbatch>` (Getting Started example from https://openfoam.org/download/8-source/) with the below decomposeParDict.

Using the provided decomposeParDict from https://openfoamwiki.net/index.php/DecomposePar ::

    /*---------------------------------------------------------------------------*\
    | =========                 |                                                 |
    | \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\    /   O peration     | Version:  1.3                                   |
    |   \\  /    A nd           | Web:      http://www.openfoam.org               |
    |    \\/     M anipulation  |                                                 |
    \*---------------------------------------------------------------------------*/

    FoamFile
    {
        version         2.0;
        format          ascii;

        root            "";
        case            "";
        instance        "";
        local           "";

        class           dictionary;
        object          decomposeParDict;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


    numberOfSubdomains 4;

    method          simple;

    simpleCoeffs
    {
        n               (1 4 1);
        delta           0.001;
    }

    hierarchicalCoeffs
    {
        n               (1 1 1);
        delta           0.001;
        order           xyz;
    }

    metisCoeffs
    {
        processorWeights
        (
            1
            1
            1
        );
    }

    manualCoeffs
    {
        dataFile        "";
    }

    distributed     no;

    roots
    (
    );


    // ************************************************************************* //


Module files are available below:

- :download:`/usr/local/modulefiles/live/eb/all/OpenFOAM/8-foss-2020a </bessemer/software/modulefiles/OpenFOAM/8-foss-2020a>`
