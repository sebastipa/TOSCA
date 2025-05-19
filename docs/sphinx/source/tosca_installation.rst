.. _installation_section:

Installation
============

In order to be installed, TOSCA requires a working C/C++ compiler, PETSc (version 3.14.x, 3.15.x), Open MPI (version 4.0.x, 4.1.x), HDF5 and
HYPRE (needed by PETSs in order to build some of the matrix solvers we use). TOSCA has been tested with the above version combinations,
it could work with other combinations or versions but it has not been tested (especially older versions).
We recommend the following versions of the above libraries:

- gcc      : 9.2.0  (https://gcc.gnu.org/).
- PETSc    : 3.15.5 (https://ftp.mcs.anl.gov/pub/petsc/).
- Open MPI : 4.1.2  (https://www.open-mpi.org/software/ompi/v4.1/).
- HYPRE    : 2.20.0 (https://github.com/hypre-space/hypre/tree/master) (check version in /src/CMakeLists.txt).
- HDF5     : 1.12.1 (https://www.hdfgroup.org/downloads/hdf5/).

Prior to install TOSCA, we suggest to create a folder named ``software`` inside ``$HOME``, where the PETSc, TOSCA and potentially other libraries directories will be located.
In order to compile TOSCA on your system, please follow these steps:

1. Check your compiler version with ``gcc --version``
2. Download PETSc into ``$HOME/software/``
3. Download Open MPI: you can download the binaries or compile from source (the latter is recommended if using ``environment-modules`` e.g. on an HPC architecture).
   If you have only one version of Open MPI installed on your system in the ``/usr`` directory (installed using sudo for example), you can omit the
   ``--with-mpi-dir='your--path--to--mpicc'`` at point 6: Open MPI will be found by the library locator.
4. Download HYPRE into ``$HOME/software/``. This can be omitted if you use the ``--download-hypre`` option in the PETSc configure step (suggested).
5. Download HDF5 into ``$HOME/software/``. This can be omitted if you use the ``--download-hdf5`` option in the PETSc configure step (suggested).
6. Configure PETSc. We suggest the following configure options, which will automatically compile HYPRE and HDF5. 

   .. code-block:: bash
   
   	./configure --with-fc=0 --download-f2cblaslapack --download-hypre --download-hdf5 --with-64-bit-indices=1 --with-debugging=0
   	
   Note that other options are available. For example, if one wants to specify paths to already installed OpenMPI libraries, 
   the options ``--with-mpi-dir='your--path--to--mpicc'`` can be used. We strongly suggest to use the ``--download-hypre`` and ``--download-hdf5`` options, 
   as they will download and compile the libraries automatically, making hte installation process easier.

7. Make PETSc (PETSc will suggest a command after the configure, we advice to use that):
	
   .. code-block:: bash
   
   	make all

8. Check the PETSc installation (PETSc will suggest a command after compilation, we advice to use that):
	
   .. code-block:: bash
   
   	make check
   
9. Save an environment variable that will tell TOSCA where PETSc is installed in your ``.bashrc``. For HPC installations
   it might be a better practice to use the ``.bash_profile`` instead.

   .. code-block:: bash
   
   	echo "export PETSC_DIR=$HOME/your--path--to--petsc--dir" >> $HOME/.bashrc
   	
   or, for HPCs
   
   .. code-block:: bash
   
   	echo "export PETSC_DIR=$HOME/your--path--to--petsc--dir" >> $HOME/.bash_profile
   
10. Save an environment variable that will tell TOSCA which PETSc architecture is required in your ``.bashrc``. For HPC installations
    it might be a better practice to use the ``.bash_profile`` instead.

   .. code-block:: bash
    
    echo "export PETSC_ARCH=arch-linux-c-opt" >> $HOME/.bashrc
    
   or, for HPCs
    
   .. code-block:: bash
    
    echo "export PETSC_ARCH=arch-linux-c-opt" >> $HOME/.bash_profile
    
.. note::

	This is the folder within ``$PETSC_DIR`` with a name beginning with ``arch-``. In a typical installation, it will be ``arch-linux-c-opt``.

11. Save an environment variable that will tell TOSCA where HDF5 libraries/include files will be located in your ``.bashrc``. For HPC installations 
    it might be a better practice to use the ``.bash_profile`` instead. This variable will be used in the TOSCA ``makefile`` to compile ``tosca2PV``. 
    If you have installed HDF5 using the ``--download-hdf5`` option in the PETSc configure step, the HDF5 `lib` and `include` directories will be located in 
    the ``$PETSC_DIR/$PETSC_ARCH/`` folder, together with HYPRE and PETSc ones.

    .. code-block:: bash
    
     echo "export HDF5_DIR=$PETSC_DIR/$PETSC_ARCH" >> $HOME/.bashrc
     
    or, for HPCs
    
    .. code-block:: bash
    
     echo "export HDF5_DIR=$PETSC_DIR/$PETSC_ARCH" >> $HOME/.bash_profile

12. Add the PETSc shared libraries to your library path environment variable in your ``.bashrc``. For HPC installations
    it might be a better practice to use the ``.bash_profile`` instead.

    .. code-block:: bash
   
   	 echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib" >> $HOME/.bashrc
   	 
    or, for HPCs
    
    .. code-block:: bash
   
   	 echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib" >> $HOME/.bash_profile

13. Reload the environment with 

	 .. code-block:: bash
	
	  source $HOME/.bashrc
	 
	 of, for HPCs
	 
	 .. code-block:: bash
	
	  source $HOME/.bash_profile
    
14. Go inside ``TOSCA/src`` directory and compile the executables. Make sure to have added the path to custom libraries (i.e. compiled
    locally by the user) if any, such as OpenMPI and HDF5, to the ``LD_LIBRARY_PATH`` as done for the PETSc libraries at step 10. 

    .. code-block:: bash 
    
     make tosca
     make tosca2PV

15. Test the installation by copying ``tosca`` and ``tosca2PV`` in one of the example cases and run the simulation
    and the post-processing with ``./tosca`` and ``./tosca2PV`` respectively. To run in parallel you have to use
    ``mpirun -np 'your-number-of-processors' ./tosca``
    
.. tip:: 

	If you run using ``mpirun ./tosca``, MPI will use the maximum number of processors available. 

.. tip::

	You can add ``tosca`` and ``tosca2PV`` to the ``PATH`` so that they will be found as executables without the need to copy them inside the case directory.

.. _contribute_subsection: 

Contribute to the TOSCA Project
-------------------------------

The TOSCA repository is open-source, so anyone can download and use the code. If you want to contribute to the project by adding 
code to TOSCA repository you need to open a pull-request that has to be approved by our team. In order to do so, please use the following steps:

1. Clone the TOSCA package locally on your machine with ``git clone https://github.com/sebastipa/TOSCA.git``
2. Create a new local branch with ``git checkout -b your-branch-name``
3. Make the desired changes to the TOSCA code, then check which files have been modified with ``git status``
4. Add changes to the git stack with ``git add modified-files``
5. Commit the changes using a short but exhaustive comment with ``git commit -m "your-commit-description"``
6. Push your local branch online with ``git push origin your-branch-name``
7. Go to github, select your branch, click on *Contribute* and open a pull-request describing the motivation of your changes, their effect on the code and the tests you performed.
8. After approval of the pull-request by our team, commits will be added to the main TOSCA version
9. To stay up-to-date, rebase your local master with the your new commits by first checking out in your local master branch with ``git checkout master`` and then rebase with ``git pull --rebase origin master``
10. Delete your local branch as it is not useful anymore with ``git branch --delete -d your-branch-name`` and use the master until you want to make new changes

.. _paraview-catalyst-section:

Paraview-Catalyst2 OS-Rendering
-------------------------------

TOSCA provides full interface with Paraview-catalyst2 through the ``USE_CATALYST`` flag in the makefile. If this is the case, the ``CATALYST`` environment variable should point to the catalyst2 installation directory. Paraview-catalyst is optional and can be disabled by setting ``USE_CATALYST=0``

Usage
*****

In order to activate off-screen rendering capabilities, ``-pvCatalyst=1`` should be set in the *control.dat* file. A file called *catalystProperties* will be required inside the sampling directory. Entries to this file are

+------------------+-----------------------------------------------------+
| ``ioType``       | can be set to *script* or *general*                 |
+------------------+-----------------------------------------------------+
| ``outputType``   | can be set to *timeStep* or *adjustableTime*        |
+------------------+-----------------------------------------------------+
| ``startTime``    | model-time at which catalyst actions start          |
+------------------+-----------------------------------------------------+
| ``timeInterval`` | acquisition period in seconds if ``outputType``     |
|                  | is set to *adjustableTime* or iterations if         |
|                  | ``outputType`` is set to *timeStep*                 |
+------------------+-----------------------------------------------------+
| ``scriptName``   | name of the catalyst actions python script.         |
|                  | Only required if ``ioType=script``                  |
+------------------+-----------------------------------------------------+

If ioType=general, 3D fields of velocity magnitude, pressure and q-criterion are saved inside the catalyst/ folder.
If ioType=script, Praview actions defined in the python script are executed and e.g. png images can be saved at runtime.

Installation
************

In order to be installed, Paraview-catalyst2 requires a working C/C++ compiler, Open MPI (version 4.0.x, 4.1.x), Python3 and cmake. In order for Paraview to work, OpenGL must be available at runtime or mesa libraries are required to mimic some hardware components. These are usually available on supercomputers through the 'mesa' module, which should be loaded at runtime. Lastly, paraview and catalyst2 should be manually compiled of the system. As Paraview-5.10 contains a bug in the definition of rectilinear mesh (used by TOSCA), Paraview-5.11 or later is recommended.

Prior to install TOSCA, we suggest to create a folder named ``Software`` inside ``$HOME``, where catalyst2 and paraview-5.11 will be located.
In order to re-compile TOSCA with Paraview-catalyst2 capabilities on your system, please follow these steps:

1. Install catalyst2:

.. code-block:: bash

	export LOCATION=$HOME/Software
	cd $LOCATION
	mkdir catalyst2 && cd catalyst2
	git clone https://gitlab.kitware.com/paraview/catalyst.git catalyst-src && cd catalyst-src
	mkdir -p build && cd build
	cmake .. -DCMAKE_INSTALL_PREFIX=$LOCATION/catalyst2/install
	make
	make install
	echo "export CATALYST=$LOCATION" >> $HOME/.bashrc

Add the Catalyst2 shared libraries to your library path environment variable in your .bashrc with ``echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/catalyst2/install/lib" >> $HOME/.bashrc`` 

.. note::

	In some cases you may have to replace ``lib`` with ``lib64``

2. Install Paraview-5.11:

.. code-block:: bash

	cd $LOCATION
	mkdir paraview-5.11.0 && cd paraview-5.11.0
	wget https://www.paraview.org/files/v5.11/ParaView-v5.11.0.tar.xz
	tar -xvf ParaView-v5.11.0.tar.xz
	mv ParaView-v5.11.0 paraview-src && cd paraview-src
	mkdir -p build && cd build
	export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$LOCATION/catalyst2/install/lib/cmake/catalyst-2.0 
	
.. note::

	In some cases you may have to replace lib with lib64	

.. code-block:: bash

	FLAGS=(-DCMAKE_INSTALL_PREFIX=$LOCATION/paraview-5.11.0/install -DVTK_OPENGL_HAS_OSMESA=ON -DPARAVIEW_USE_MPI=ON -DBUILD_TESTING=OFF -DVTK_USE_X=OFF -DPARAVIEW_USE_QT=OFF -DPARAVIEW_USE_PYTHON=ON -DPython3_FIND_STRATEGY=LOCATION -DPython3_ROOT_DIR=$EBROOTPYTHON -DPARAVIEW_BUILD_SHARED_LIBS=ON -DPARAVIEW_ENABLE_RAYTRACING=OFF -DPARAVIEW_ENABLE_CATALYST=ON )
	cmake .. "${FLAGS[@]}"
	make -j8
	make install
	
Add the Paraview shared libraries to your library path environment variable in your .bashrc with ``echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LOCATION/paraview-5.11.0/install/lib" >> $HOME/.bashrc``

.. note::

	In some cases you may have to replace ``lib`` with ``lib64``

3. Re-compile TOSCA

   * Reload the environment with ``source $HOME/.bashrc``
   * Go inside ``TOSCA/src`` directory and recompile the solver with ``make tosca``, ensuring that ``-DUSE_CATALYST=1`` in the makefile

4. Running

   In order for catalyst2 to find paraview library, the following envronment variables should be set at runtime:

   * ``export CATALYST_IMPLEMENTATION_PATHS=$LOCATION/paraview-5.11.0/install/lib/catalyst`` 
   * ``export CATALYST_IMPLEMENTATION_NAME=paraview``
   
.. note::

	In some cases you may have to replace ``lib`` with ``lib64``
