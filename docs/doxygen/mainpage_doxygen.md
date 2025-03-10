TOSCA       {#mainpage}
============

\image html tosca.png width=400px
\image latex tosca.png width=400px

Toolbox fOr Stratified Convective Atmospheres

## Recent Highlights/Additions

    - Direct/Indirect profile assimilation techniques to drive LES with observations/mesoscale models
    - New scale-dependent LES model
    - New stability dependent wall model for IBM
    - IBM also for concurrent precursor (to model terrain features)
    - Lateral fringe region.

Now we can simulate real cases, wind rotating 360 degs, with terrain and turbines, all while damping gravity waves in all directions !!


## Executables

    - tosca    : transient solver for stratified incompressible flows. Temperature stratification is accounted via Boussinesq approximation.
    - tosca2PV : post processor for ParaView visualization. Writes data in XMF/HDF format.

## Boundary Conditions

    - noSlip               : on all patches
    - slip                 : on all patches
    - zeroGradient         : on all patches
    - fixedGradient        : on all patches (only T equation)
    - inflowFunction       : only on kLeft patch. Different types of inflow can be prescribed: ABL inflow using teorethical laws, mapped inflow from a precursor
                             simulation that can be also periodized and/or interpolated if meshes are not consistent.
    - velocityWallFunction : shear stress model (only U equation). Can be prescribed on jRight and jLeft patches
    - thetaWallFunction    : potential temperature wall model (only T equation). Can be prescribed on jRight and jLeft patches

## Wall Models

    - Shumann    : applies wall shear stress in the momentum equation according to similarity theory of Paulson and Obuhkov.
                   Velocity BC is dependent (velocityWallFunction/thetaWallFunction, type -3).
    - Cabot      : prescribes velocity BC according to Cabot formulation (velocityWallFunction, type -1).

## Turbulence Models

    - Smagorinsky standard model (les 1).
    - Dynamic Smagorinsky LES turbulence model with cubic (les 3) or lagrangian (les 4) averaging.

## Turbine Models

    - UADM : base actuator disk model with imposed Ct and yaw controller
    - ADM  : advanced actuator disk model with rotor dynamics and yaw/pitch/rotation controllers
    - ALM  : advanced actuator line model with rotor dynamics and yaw/pitch/rotation controllers (anisotropic gaussuan projection avail.)
    - AFM  : anisotropic actuator farm model for coarser meshes, integral/point sampling available

## IBM Method

    TOSCA uses a sharp-interface IBM method with ghost cells applied on the flow side of the body. This allows for more flexibility with
    thin bodies but complicated wall modeling. New wall models have been recently added, validated with both hi-Re terrains (Shumann) or small objects (Cabot).
    Moving IBM is also possible, we are currently validating the latter capability.

## Acquisition System

    - ABL wall-parallel planes-averaged statistics as a function of time
    - field averages
    - domain sections
    - probes with advanced parallel I/O writing
    - turbine data with advanced parallel I/O writing
    - ABL perturbations w.r.t. reference state
    - layer averaged fields (three layers available - 3LM)
    - mechanical energy budgets within user-defined boxes (the code is not energy-conservative so MKE eq. is not fullfilled)

## Future Implementations:

    - Overset with fringe interpolation
    - IBM smooth normals and thin-body augmented search

## Notes:

    - test cases which contain an empty 'inflowDatabase' file require the inflow database. A sample database can be downloaded
      at https://drive.google.com/file/d/17F5wtI5Jc1XGh8crmOVJYVXabC8iQXq1/view?usp=sharing, simply substitute the 'inflowDatabase'
      file with the downloaded folder.

## Installation:

In order to be installed, TOSCA requires a working C/C++ compiler, PETSc (version 3.14.x, 3.15.x), Open MPI (version 4.0.x, 4.1.x), HDF5 and
HYPRE (needed by PETSs in order to build some of the matrix solvers we use). TOSCA has been tested with the above version combinations,
it could work with other combinations or versions but it has not been tested (especially older versions).
We recommend the following versions of the above libraries:
 * gcc      : 9.2.0  (https://gcc.gnu.org/).
 * PETSc    : 3.15.5 (https://ftp.mcs.anl.gov/pub/petsc/).
 * Open MPI : 4.1.2  (https://www.open-mpi.org/software/ompi/v4.1/).
 * HYPRE    : 2.20.0 (https://github.com/hypre-space/hypre/tree/hypre_petsc) (check version in /src/CMakeLists.txt).
 * HDF5     : 1.12.1 (https://www.hdfgroup.org/downloads/hdf5/).

Prior to install TOSCA, we suggest to create a folder named `Software` inside `$HOME`, where the PETSc, HYPRE and TOSCA folders will be located.
In order to compile TOSCA on your system, please follow these steps:

* Check your compiler version with `gcc --version`

* Download PETSc into `$HOME/Software/`

* Download HYPRE `$HOME/Software/`

* Download Open MPI: you can download the binaries or compile from source (the latter is recommended if use `environment-modules`).
  If you have only one version of Open MPI installed on your system in the `/usr` directory (using sudo for example), you can omit the
  `--with-mpi-dir='your--path--to--mpicc'` at point 4: Open MPI will be found by the 'ld' library locator.

* Configure PETSc (will automatically compile HYPRE). We suggest the following configure options:
  `./configure --with-fc=0 --download-f2cblaslapack --with-mpi-dir='your--path--to--mpicc' --download-hypre='your--path--to--hypre \
  --with-64-bit-indices=1 --with-debugging=0`

* Make PETSc with `make all`

* Test PETSc with `make check`

* Save an environment variable that will tell TOSCA where PETSc is installed in your .bashrc:
  `echo "export PETSC_DIR=$HOME/your--path--to--petsc" >> $HOME/.bashrc`

* Save an environment variable that will tell TOSCA which PETSc architecture is required in your .bashrc. Note: this is the folder within $PETSC_DIR with a name beginning with "arch-". In a typical installation, it will be "arch-linux-c-opt":
  `echo "export PETSC_ARCH=arch-linux-c-opt" >> $HOME/.bashrc`

* Add the PETSc shared libraries to your library path environment variable in your .bashrc:
  `echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib" >> $HOME/.bashrc`

* Reload the environment with `source $HOME/.bashrc`

* Go inside TOSCA/src directory and compile the executables with `make tosca` and `make tosca2PV`

* Test the installation by copying `tosca` and `tosca2PV` in one of the example cases and run the simulation
  and the post-processing with `./tosca` and `./tosca2PV` respectively. To run in parallel you have to use
  `mpirun -np 'your-number-of-processors' ./tosca`

## Contribute to the TOSCA Project

  The TOSCA repository is open-source, so anyone can download and use the code. If you want to contribute to the project by adding code to TOSCA repository you need to open a pull-request that has to be approved by our team. In order to do so, please use the following steps:

  * Clone the TOSCA package locally on your machine with `git clone https://github.com/sebastipa/TOSCA.git`

  * Create a new local branch with `git checkout -b your-branch-name`

  * Make the desired changes to the TOSCA code, then check which files have been modified with `git status`

  * Add changes to the git stack with `git add modified-files`

  * Commit the changes using a short but exhaustive comment with `git commit -m "your-commit-description"`

  * Push your local branch online with `git push origin your-branch-name`

  * Go to github, select your branch, click on "Contribute" and open a pull-request describing the motivation of your changes, their effect on the code and the tests you performed.

  * After approval of the pull-request by our team, commits will be added to the main TOSCA version

  * To stay up-to-date, rebase your local master with the your new commits by first checking out in your local master branch with `git checkout master` and then rebase with `git pull --rebase origin master`

  * Delete your local branch as it is not useful anymore with `git branch --delete -d your-branch-name` and operate on the master until you want to make new changes

## Paraview-Catalyst2 OS-Rendering

  TOSCA provides full interface with Paraview-catalyst2 through the USE_CATALYST flag in the makefile. If this is the case, the CATALYST environment variable should point to the catalyst2 installation directory. Paraview-catalyst is optional and can be disabled by setting USE_CATALYST=0

## Usage

  In order to activate off-screen rendering capabilities, -pvCatalyst=1 should be set in the control.dat file. A file called catalystProperties will be required inside the sampling directory. Entries to this file are

  - ioType          can be set to 'script' or 'general'
  - outputType      can be set to 'timeStep' or 'adjustableTime'
  - startTime       model-time at which catalyst actions start
  - timeInterval    acquisition period in seconds if outputType=adjustableTime or iterations if outputType=timeStep
  - scriptName      name of the catalyst actions python script. Only required if ioType=script

## What does it do

  If ioType=general, 3D fields of velocity magnitude, pressure and q-criterion are saved inside the catalyst/ folder.
  If ioType=script, Praview actions defined in the python script are executed and e.g. png images can be saved at runtime.

## Installation:

In order to be installed, Paraview-catalyst2 requires a working C/C++ compiler, Open MPI (version 4.0.x, 4.1.x), Python3 and cmake. In order for Paraview to work, OpenGL must be available at runtime or mesa libraries are required to mimic some hardware components. These are usually available on supercomputers through the 'mesa' module, which should be loaded at runtime. Lastly, paraview and catalyst2 should be manually compiled of the system. As Paraview-5.10 contains a bug in the definition of rectilinear mesh (used by TOSCA), Paraview-5.11 or later is recommended.

Prior to install TOSCA, we suggest to create a folder named `Software` inside `$HOME`, where catalyst2 and paraview-5.11 will be located.
In order to re-compile TOSCA with Paraview-catalyst2 capabilities on your system, please follow these steps:

## 1. install Catalyst2

* `export LOCATION=$HOME/Software`

* `cd $LOCATION`

* `mkdir catalyst2 && cd catalyst2`

* `git clone https://gitlab.kitware.com/paraview/catalyst.git catalyst-src && cd catalyst-src`

* `mkdir -p build && cd build`

* `cmake .. -DCMAKE_INSTALL_PREFIX=$LOCATION/catalyst2/install`

* `make`

* `make install`

* `echo "export CATALYST=$LOCATION" >> $HOME/.bashrc`

* Add the Catalyst2 shared libraries to your library path environment variable in your .bashrc:
`echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/catalyst2/install/lib" >> $HOME/.bashrc` Note: in some cases you may have to replace lib with lib64

## 2. install Paraview-5.11

* `cd $LOCATION`

* `mkdir paraview-5.11.0 && cd paraview-5.11.0`

* `wget https://www.paraview.org/files/v5.11/ParaView-v5.11.0.tar.xz`

* `tar -xvf ParaView-v5.11.0.tar.xz`

* `mv ParaView-v5.11.0 paraview-src && cd paraview-src`

* `mkdir -p build && cd build`

* `export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$LOCATION/catalyst2/install/lib/cmake/catalyst-2.0` Note: in some cases you may have to replace lib with lib64

* `FLAGS=(-DCMAKE_INSTALL_PREFIX=$LOCATION/paraview-5.11.0/install -DVTK_OPENGL_HAS_OSMESA=ON -DPARAVIEW_USE_MPI=ON -DBUILD_TESTING=OFF -DVTK_USE_X=OFF -DPARAVIEW_USE_QT=OFF -DPARAVIEW_USE_PYTHON=ON -DPython3_FIND_STRATEGY=LOCATION -DPython3_ROOT_DIR=$EBROOTPYTHON -DPARAVIEW_BUILD_SHARED_LIBS=ON -DPARAVIEW_ENABLE_RAYTRACING=OFF -DPARAVIEW_ENABLE_CATALYST=ON )`

* `cmake .. "${FLAGS[@]}"`

* `make -j8`

* `make install`

* Add the Paraview shared libraries to your library path environment variable in your .bashrc:
`echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/paraview-5.11.0/install/lib" >> $HOME/.bashrc` Note: in some cases you may have to replace lib with lib64

## 3. Re-compile TOSCA

* Reload the environment with `source $HOME/.bashrc`

* Go inside `TOSCA/src` directory and recompile the solver with `make tosca` ensuring that `-DUSE_CATALYST=1` in the makefile

## 4. Running

In order for catalyst2 to find paraview library, the following envronment variables should be set at runtime:

* `export CATALYST_IMPLEMENTATION_PATHS=$LOCATION/paraview-5.11.0/install/lib/catalyst` Note: in some cases you may have to replace lib with lib64

* `export CATALYST_IMPLEMENTATION_NAME=paraview`


Credits and Copyright: Sebastiano Stipa - Arjun Ajay - Mohammad Hadi - The University of British Columbia
