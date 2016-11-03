.. _iceberg_comsol:

COMSOL Multiphysics
===================

.. sidebar:: COMSOL

   :Version:  5.2
   :Support Level: ?
   :Dependencies: None
   :URL: https://www.comsol.com/comsol-multiphysics

COMSOL Multiphysics provides a range of simulation, finite element analysis and solver functionality. 

Availability, licensing and modules
-----------------------------------

COMSOL Multiphysics is installed on Iceberg but is not presently supported by CiCS.  There is no site license for COMSOL so to run it on Iceberg you will need either your own license or access to licenses managed by others.  For example, a research group may operate a *license server* that can issue COMSOL licenses to (approved) users.  If you need access to a license on a temporary basis then you could try contacting `Prof. Will Zimmerman <https://www.sheffield.ac.uk/cbe/staff/staffprofiles/wzimmerman>`_, who has purchased several licenses for his research group and may be willing to hire out licenses to other researchers.  See `Interactive usage` for information on how to specify a COMSOL license when starting COMSOL on Iceberg.

The following COMSOL Multiphysics 5.2 modules are installed:

* ACDC
* Acoustics
* Batteries and Fuel Cells
* CAD Import
* CFD
* Chemical Reaction Engineering
* COMSOL Multiphysics
* Corrosion
* Design
* ECAD Import
* Electrochemistry
* Electrodeposition
* Fatigue
* Geomechanics
* Heat Transfer
* LiveLink for Excel
* LiveLink for MATLAB
* MEMS
* Microfluidics
* Mixer
* Molecular Flow
* Multibody Dynamics
* Nonlinear Structural Materials
* Optimization
* Particle Tracing
* Pipe Flow
* Plasma
* Ray Optics
* RF
* Semiconductor
* Structural Mechanics
* Subsurface Flow
* Wave Optics

Interactive usage
-----------------

After connecting to iceberg (see :ref:`ssh`),  start an interactive graphical session with the :code:`qrshx` command. 
Alternatively, if you require more memory, for example 16 GB, use the command :code:`qrshx -l rmem=16G` 

Next, run the following to make COMSOL Multiphysics available in your current session: ::

        module load apps/binapps/comsol/5.2

Specify the location of your license information: ::

        export LM_LICENSE_SERVER=/home/myusername/path/to/mylicensefile.dat:$LM_LICENSE_SERVER

Finally, start COMSOL Multiphysics: ::

	comsol

The COMSOL Multiphysics user interface should then appear.  Here ``mylicensefile.dat`` is a file containing either:

* details of your license (which components you can use) *or*
* details of the *license server* that you want COMSOL to request licenses from.

If you are using a license server then your license file needs to contain no more than the following: ::

        SERVER mylicenseserver.sheffield.ac.uk ANY 65321
        USE_SERVER
        
where ``mylicenseserver.sheffield.ac.uk`` is the hostname of your license server and ``654321`` is the *port* to connect to on that machine to request a COMSOL license.

The person responsible for managing the license server may ask for your Iceberg username to allow you to request licenses (whilst preventing others from doing so).

Installation instructions
-------------------------

This section is primarily of interest

Installation notes
------------------

This section is primarily of interest to system administrators.

Version 5.2
^^^^^^^^^^^

No installation notes are available.

:download:`This modulefile </iceberg/software/modulefiles/apps/binapps/comsol/5.2>` was installed as ``/usr/local/modulefiles/apps/binapps/comsol/5.2``.
