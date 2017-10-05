OpenFOAM
========

.. sidebar:: OpenFOAM
   
   :Version: 4.1
   :Dependencies: GCC compiler and Open MPI; Scotch/PT-Scotch v6.0.3
   :URL: https://openfoam.org 
   :Documentation: https://cfd.direct/openfoam/user-guide/


OpenFOAM is leading software for computational fluid dynamics (CFD). It is licensed free and open source only under the GNU General Public Licence (GPL) by the OpenFOAM Foundation.


Usage
-----

OpenFOAM 4.1 can be activated using the module file::

    module load apps/openfoam/4.1/gcc-6.2-openmpi-2.0.1

Installation notes
------------------

OpenFOAM 4.1 was installed using the
:download:`install_openfoam.sh </sharc/software/install_scripts/apps/openfoam/4.1/install_openfoam.sh>` script, the module
file is
:download:`gcc-6.2-openmpi-2.0.1 </sharc/software/modulefiles/apps/openfoam/4.1/gcc-6.2-openmpi-2.0.1>`. The following optional dependencies were built as part of the installation process: Scotch/PT-Scotch v6.0.3 (located in /usr/local/packages/apps/openfoam/4.1/gcc-6.2-openmpi-2.0.1/ThirdParty-4.1). The following optional dependencies were not installed: ParaView and CGAL.

Installation was tested as follows (Getting Started example from https://openfoam.org/download/4-1-source/)::

    $ mkdir /data/$USER/tests/openfoam/run
    $ cd /data/$USER/tests/openfoam/run
    $ module load apps/openfoam/4.1/gcc-6.2-openmpi-2.0.1
    $ cp -r $FOAM_TUTORIALS/incompressible/simpleFoam/pitzDaily .
    $ cd pitzDaily
    $ blockMesh
    $ simpleFoam
    $ paraFoam
