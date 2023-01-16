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

* `cmake .. -DCMAKE_INSTALL_PREFIX=$LOCATION/catalyst2`

* `make`

* `make install`

* `echo "export CATALYST=$LOCATION >> $HOME/.bashrc`

* Add the Catalyst2 shared libraries to your library path environment variable in your .bashrc:
  `echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/catalyst2/install/lib" >> $HOME/.bashrc` Note: in some cases you may have to replace lib with lib64

## 2. install Paraview-5.11

* `cd $LOCATION`

* `mkdir paraview-5.11.0 && cd paraview-5.11.0`

* `wget https://www.paraview.org/files/v5.11/ParaView-v5.11.0.tar.xz`

* `tar -xvf ParaView-v5.11.0`

* `mv ParaView-v5.11.0 paraview-src && cd paraview-src`

* `mkdir -p build && cd build`

* `export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$LOCATION/catalyst2/install/lib/cmake/catalyst-2.0` Note: in some cases you may have to replace lib with lib64

* `FLAGS=( -DCMAKE_INSTALL_PREFIX=$LOCATION/paraview-5.11.0/install \
        -DVTK_OPENGL_HAS_OSMESA=ON -DPARAVIEW_USE_MPI=ON -DBUILD_TESTING=OFF -DVTK_USE_X=OFF    
        -DPARAVIEW_USE_QT=OFF -DPARAVIEW_USE_PYTHON=ON -DPython3_FIND_STRATEGY=LOCATION \
        -DPython3_ROOT_DIR=$EBROOTPYTHON -DPARAVIEW_BUILD_SHARED_LIBS=ON \
        -DPARAVIEW_ENABLE_RAYTRACING=OFF -DPARAVIEW_ENABLE_CATALYST=ON ) `

* `cmake .. "${FLAGS[@]}"`

* `make -j8`

* `make install`

* Add the Paraview shared libraries to your library path environment variable in your .bashrc:
  `echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/paraview-5.11.0/install/lib" >> $HOME/.bashrc` Note: in some cases you may have to replace lib with lib64

## 3. Re-compile TOSCA

* Reload the environment with `source $HOME/.bashrc`

* Go inside `TOSCA/src` directory and recompile the solver with `make tosca` ensuring that `-DUSE_CATALYST=1` in the makefile

Credits and Copyright: Sebastiano Stipa - Arjun Ajay - Mohammad Hadi - The University of British Columbia
