User Guide
==========

TOSCA uses ASCII input files, organized in files, dictionaries and subdictionaries. The latter two levels of description are always embodied using ’{}’ parenthesis. The code provides some level of input checking, meaning that non-recognized inputs are followed by an error message listing available possibilities. 

TOSCA has a standardized case structure. The minimum-required case structure is depicted on the right of the following figure, while the case structure required to run ABL and wind farm simulations under potential temperature stratification (i.e. with the addition of the ``boundary/T`` and ``ABLProperties.dat`` files) is shown on the right. The principal control file for a TOSCA simulation is the `control.dat` file, located in the case directory. Depending on the type of simulation that one wishes to perform, flags can be activated in the `control.dat` which prompt TOSCA to read additional input files and data. 

.. image:: ../images/structure_1.png
   :align: center
   
Inside the ``boundary`` directory, initial and boundary conditions are set in files having the same name as the field they refer to. 

.. _spatial-mesh-section: 

Spatial Mesh 
------------

The spatial mesh can be provided in two different formats. These can be selected with the ``-meshFileType`` keyword in the `control.dat` file, which can be set to ``cartesian`` or ``curvilinear`` for ``.xyz`` and ``.grid`` formats, respectively. While both formats describe a structured mesh (TOSCA can only operate with this type of mesh), the first one is additionally cartesian. Both types of mesh can be defined such that grid points are stretched in all spatial directions. The ``.xyz`` mesh format is TOSCA's custom Cartesian mesh format, and it is very convenient since the user only has to provide the discretization along each of the three cartesian axes. The ``.xyz`` mesh file format can be generalized as

.. code-block:: bash

	Nx Ny Nz
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=2,j=1,i=1) y(k=2,j=1,i=1) z(k=2,j=1,i=1)
	:
	x(k=Nk ,j=1,i=1) y(k=Nk ,j=1,i=1) z(k=Nk ,j=1,i=1)
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=1,j=1,i=2) y(k=1,j=1,i=2) z(k=1,j=1,i=2)
	:
	x(k=1,j=1,i=Ni) y(k=1,j=1,i=Ni) z(k=1,j=1,i=Ni)
	x(k=1,j=1,i=1) y(k=1,j=1,i=1) z(k=1,j=1,i=1)
	x(k=1,j=2,i=1) y(k=1,j=2,i=1) z(k=1,j=2,i=1)
	:
	x(k=1,j=Nj ,i=1) y(k=1,j=Nj ,i=1) z(k=1,j=Nj ,i=1)

which, for a cube extending 50 m in each direction, with 10 m cells along each direction (hence 6 points), is given by 

.. code-block:: bash

	6 6 6 
	0 0 0
	10 0 0 
	20 0 0
	30 0 0
	40 0 0
	50 0 0
	0 0 0 
	0 10 0 
	0 20 0 
	0 30 0 
	0 40 0
	0 50 0
	0 0 0 
	0 0 10
	0 0 20
	0 0 30
	0 0 40
	0 0 50

As can be noticed, the mesh file is extremely light, even for large meshes, as only the discretization along the three cartesian axes is provided. The coordinates of the remaining mesh cells are inferred iternally within TOSCA, starting from the axes discretization.   

The ``.grid`` format is a curvilinear mesh format (available for example in the meshing software `Pointwise`), where the cartesian coordinates of each mesh point are provided. This becomes useful when the mesh is deformed, and a unique relationship between the cartesian and curvilinear sets of coordinates does not exist anymore. This is a much heavier mesh format, which can be generalized as follows

.. code-block:: bash

	Nj Nk Ni
	x(k=1,j=1,i=1) ..... x(k=1,j=Nj,i=1)
	:
	:
	x(k=Nk,j=1 ,i=1) ....x(k=Nk,j=Nj ,i=1)
	x(k=1,j=1,i=2) ....x(k=1,j=Nj,i=2)
	:
	:
	x(k=Nk ,j=1 ,i=Ni) ...x(k=Nk ,j=Nj ,i=Ni)
	y(k=1,j=1,i=1) ..... y(k=1,j=Nj,i=1)
	:
	:
	y(k=Nk,j=1 ,i=1) ....y(k=Nk,j=Nj ,i=1)
	y(k=1,j=1,i=2) ....y(k=1,j=Nj,i=2)
	:
	:
	y(k=Nk ,j=1 ,i=Ni) ...y(k=Nk ,j=Nj ,i=Ni)
	z(k=1,j=1,i=1) ..... z(k=1,j=Nj,i=1)
	:
	:
	z(k=Nk,j=1 ,i=1) ....z(k=Nk,j=Nj ,i=1)
	z(k=1,j=1,i=2) ....z(k=1,j=Nj,i=2)
	:
	:
	z(k=Nk ,j=1 ,i=Ni) ...z(k=Nk ,j=Nj ,i=Ni)

Taking as an example the same box extending 50 m in each direction, with 10 m cells along each direction, this can be expressed using the ``.grid`` format as

.. code-block:: bash

   6 6 6
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   10 10 10 10 10 10
   20 20 20 20 20 20
   30 30 30 30 30 30
   40 40 40 40 40 40
   50 50 50 50 50 50
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   0  0  0  0  0  0
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   10 10 10 10 10 10
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   20 20 20 20 20 20
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   30 30 30 30 30 30
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   40 40 40 40 40 40
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50
   50 50 50 50 50 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50
   0 10 20 30 40 50 
 
.. note::

   Both in the ``.grid`` and ``.xyz`` format, the user has to provide the *mesh points* rather then the *cell centers* coordinates.
   
Notably, every cartesian mesh can be expressed using the curvilinear format, while the opposite is not true. TOSCA always operates using curvilinear coordinates internally, and since boundary conditions are given on curvilinear patches one has to bear in mind a few concepts, explained below, in order to assign BCs in the way he intends to. Firstly, one can think of curvilinear coordinates as the mesh cell indices in each of the three directions, expressed in TOSCA as `k`, `j` and `i`. With this in mind, it appears clear that each boundary can be identified by keeping one of the indices fixed while the others two vary. In fact, boundaries are identifyed in TOSCA using `kLeft`, `kRight`, `jLeft`, `jRight`, `iLeft` and `iRight`. When the mesh is cartesian (and the axes are aligned to the edges of the domain) boundaries can be also identifyed with the minimum and maximum of each cartesian coordinates. However, this is not general and it would not hold for a curvilinear mesh, or even a cartesian mesh with the domain edges not aligned with the cartesian coordinates. For this reason, boundaries are always identified using the curvilinear definition in TOSCA. To know how to assign boundary conditions, e.g. what boundary corresponds to the `kLeft` patch and so on, one has to keep in mind that the `.grid` file is read by TOSCA using a `i,k,j` nested loop order, namely

.. code-block:: C

   // loop through the mesh cells 
   for (i=0; i<Ni; i++) 
   {
      for (k=0; k<Nk; k++)
      {
         for (j=0; j<Nj; j++)
         {
            // read point k,j,i in the .grid file 
         }
      }
   }
   
Looking at the ``.grid`` file format, it is clear that there is no unique relation between `k,j,i` and `x,y,z`. However, when the mesh is cartesian, this relation depends on how the mesh file is created. It is easy to notice that, in the provided ``.grid`` example, `k,j,i` correspond to `x,z,y`. For the ``.xyz`` file format, where a unique relation exists, TOSCA constructs the mesh such that `k,j,i` correspond to `x,z,y` (this can be also appreciated in the provided example). In an ABL or wind farm simulation, where x is usually aligned with the wind direction, and y and z are the spanwise and vertical coordinates, respectively, the inlet will corrispond to `kLeft`, the outlet to `kRight`, the bottom and top to `jLeft` and `jRight`, respectively, while the spanwise boundaries will be `iLeft` and `iRight`. This relation will always apply as long as a cartesian mesh is used. The relations between curvilinear and cartesian coordinates depending on the mesh and domain type are summarized in the following table. 

+------------------+----------------------------+------------------------------+
|  **format**      | **cartesian domain**       | **curvilinear domain**       |
+------------------+----------------------------+------------------------------+
|   ``.grid``      | case-dependent             | unique relation impossible   |
+------------------+----------------------------+------------------------------+
|   ``.xyz``       | `k,j,i` -> `x,z,y`         | non-representable            |
+------------------+----------------------------+------------------------------+

.. _input-files-section:

Input Files
-----------

This section describes all available entries for each of the TOSCA code's input files. 

`control.dat`
~~~~~~~~~~~~~

As already mentioned, TOSCA expects a specific case structure, where the `control.dat` file contains the main simulation settings, 
some of which prompt the code to read additional input files. The `control.dat` is subdivided in 7 subsections, namely

1. Time Controls
2. I/O Controls 
3. Solution Flags 
4. Solution Controls 
5. Solution Constants 
6. Acquisition Controls 
7. Post Processing Controls 

The last group is only needed when the ``tosca2PV`` executable is launched, which converts the binary TOSCA output in ``.xmf`` 
format, which can be read from e.g. `ParaView`. The following tables summarize all available entries for each of the `control.dat` 
file subsections. 

Time Controls 
*************

.. table:: 
   :widths: 35, 65
   :align: center
                                                                                                       
   ===================== =====================================================================================================
   ``-startFrom``        It can be set to `startTime` (requires ``-startTime`` entry), or `latestTime`. In the last case, 
                         the latest time available in the ``fields`` directory is used as the initial time (the code exits if 
                         the directory is not present, or it doesn't contain any time 
                         directories).                                  
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-startTime``        Specifies the initial time of the simulation if the ``-startFrom`` keyword is set to `startTime`.  
                         Disregarded otherwise.                                                                              
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-endTime``          Defines the time at which the simulation ends.                                                      
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-timeStep``         Initial time step size in seconds.                                                                  
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-adjustTimeStep``   If set to 0, the time step will remain fixed and equal to the specified ``-timeStep``. If set to 
                         1 (requires ``-cfl``), the time step will be adjusted based on the CFL value and I/O settings. This  
                         means that time step size will be varied based on the solution, and in order to land on those    
                         time values designated for data writing by the acquisition 
                         system.                                                            
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-cfl``              Specifies the CFL value to be maintained. Disregarded if ``-adjustTimeStep`` is set to 0.            
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-timePrecision``    Specifies the number of digits after the comma that are used to write files and expected to read    
                         them. Note that, if the time folder storing the initial condition has a different number of digits, 
                         the simulation will throw an error. To solve this one can act both on the ``-timePrecision``, or     
                         rename the folder using a ``-timePrecision`` that matches that specified inside the `control.dat` 
                         file.                                       
   ===================== =====================================================================================================

I/O Controls 
************

.. table:: 
   :widths: 35, 65
   :align: center 
                                                                                                          
   ===================== =====================================================================================================
   ``-intervalType``     It can be set to `adjustableTime`, `timeStep` or `writeNow`. In the last case, a checkpoint 
                         followed by termination of the simulation will be triggered.
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-timeInterval``     Specifies how often a checkpoint file is written. If the ``-timeInterval`` is set to 
                         `adjustableTime`, the time interval between two checkpoints is expressed in seconds. If 
                         ``-timeInterval`` is set to `timeStep`, the time interval refers to the number of iterations.
   --------------------- -----------------------------------------------------------------------------------------------------
   ``-purgeWrite``       If set to 1, eliminates all previous checkpoint files every time that a checkpoint is written 
                         (in order to save disk space).
   ===================== =====================================================================================================

Solution Flags 
**************
 
.. table:: 
   :widths: 35, 65
   :align: center 
                                                                                                           
   ====================== =====================================================================================================  
   ``-les``               Specifies if and which LES model is active: 
                          (1) standard Smagorinsky
                          (2) stability dependent
                          (3) dynamic Smagorinsky with box averaging
                          (4) dynamic Smagorinsky scale invariant with lagrangian averaging (LASI)
                          (5) dynamic smagorinsky scale dependent with lagrangian averaging (LASD)
                          (6) dynamic smagorinsky scale dependent with plane averaging (PASD).
                          TOSCA has been extensively and most used adopting model (4) for ABL and wind plant simulations. 
                          Model (5) provides a good alternative in stably stratified flows.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-potentialT``        Specifies if potential temperature transport equation is solved (set to 1) or not (set to 0).
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-abl``               Specifies if an ABL simulation is run. Requires additional file ``ABLProperties.dat``.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-ibm``               Specifies if immersed bodies are present in the simulation (set to 1) or not (set to 0). Requires 
                          additional input in ``IBM`` directory.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-overset``           Specifies if overset mesh is present in the simulation. Requires additional input in ``Overset`` 
                          directory if activated.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-zDampingLayer``     Specifies if vertical Rayleigh damping layer is present in the simulation. Requires additional 
                          input in ``ABLProperties.dat`` file if activated.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-xDampingLayer``     Specifies if horizontal inlet fringe region is present in the simulation. Requires additional input  
                          in ``ABLProperties.dat`` file if activated. Concurrent precursor can be enabled with this flag.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-yDampingLayer``     Specifies if horizontal lateral fringe region is present in the simulation. Requires additional input  
                          in ``ABLProperties.dat`` file if activated. Concurrent precursor **has to** be enabled with this 
                          flag. 
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-advectionDamping``  Specifies if horizontal inlet advection damping region is present in the simulation. Requires 
                          additional input in ``ABLProperties.dat`` file if activated. 
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-precursorSpinUp``   Concurrent precursor **has to** be enabled with this flag. Specifies the type of inlet boundary 
                          condition and initial condition for the concurrent precursor. If set to 0, streamwise periodic 
                          boundary conditions are applied and initial condition is read from ``fields/precursor`` directory.
                          If set to 1, inlet planes from previous simulation are applied and flow is initialized by 
                          spreading the flow given in the plane corresponding to the simulation start time throughout the 
                          domain. If set to 2 is equivalent to 1, but the initial condition is read from the 
                          ``fields/precursor`` directory. This is used for coarse concurrent precursor, where a good solution
                          has to be continuously feed because the simulation cannot be really self-sustained in the concurrent
                          precursor. 
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-kLeftRayleigh``     Specifies if horizontal Rayleigh damping at ``kLeft`` boundary is present in the simulation. 
                          Requires additional input in ``ABLProperties.dat`` file if activated.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-kRightRayleigh``    Specifies if horizontal Rayleigh damping at ``kRight`` boundary is present in the simulation. 
                          Requires additional input in ``ABLProperties.dat`` file if activated.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-canopy``            Specifies if wind farm canopy model is present in the simulation. Requires additional input in 
                          ``ABLProperties.dat`` file if activated.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-windplant``         Specifies if wind turbines are present in the simulation (set to 1) or not (set to 0). Requires 
                          turbine models definitions in ``turbines`` directory.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-inviscid``          If set to 1, allows disabling viscous terms. Default value is 0.
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-computeContinuity`` Computes the divergence field within the entire domain and writes it to checkpoint files. 
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-pvCatalyst``        Enables `ParaView-Catalyst` off-screen rendering capabilities. Useful to create nice videos of very 
                          large simulations. More details are given in Sec. :ref:`paraview-catalyst-section`.
   ====================== =====================================================================================================

Solution Controls 
*****************

.. table:: 
   :widths: 35, 65
   :align: center 
                                                                                                           
   ========================= ====================================================================================================
   ``-meshFileType``         Defines the format of the mesh input file. It can be set to ``cartesian`` or ``curvilinear``.
                             More details are given in Sec. :ref:`spatial-mesh-section`.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-dUdtScheme``           Time discretization scheme, it can be set to ``forwardEuler`` (explicit first order, usually 
                             unstable), ``rungeKutta4`` (explicit fourth-order Runge-Kutta) or ``backwardEuler``, 
                             which corresponds to the second-order implicit Crank-Nicholson scheme (explicit selection of 
                             the Crank-Nicholson scheme will be made available). For long simulations the 
                             ``backwardEuler`` scheme is preferred, as it can run with CFL greater than 1 and is 
                             unconditionally stable. For simulations affected by constraints other than the CFL (e.g. blade 
                             rotation in actuator line model), ``rungeKutta4`` is a good alternative.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-divScheme``            Determines which divergence scheme is used for the discretization of the advection fluxes. It 
                             can be set to ``central`` (second-order symmetric scheme, dispersive), ``quickDiv`` (third-order 
                             upwind-biased quadratic scheme, diffusive), ``weno3`` (fourth-order weighted essentially 
                             non-oscillatory scheme, diffusive), ``centralUpwind`` (vanLeer blending of central and quadratic 
                             scheme, to balance diffusion and dispersion), ``centralUpwindW`` (weighted version, for 
                             graded/non-uniform meshes).
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-relTolU``              Requires ``-dUdtScheme`` set to ``backwardEuler``, discarded otherwise. Allows to set the relative 
                             exit tolerance for the Newton method used to solve implicit discretized momentum equation, default 
                             value 1e-30.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-absTolU``              Requires ``-dUdtScheme`` set to ``backwardEuler``, discarded otherwise. Allows to set the absolute 
                             exit tolerance for the Newton method used to solve implicit discretized momentum equation, 
                             default value 1e-5.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-poissonSolver``        Allows to specify the library used to solve the pressure equation, it can be set to ``HYPRE`` or 
                             ``PETSc``. ``HYPRE`` is suggested, as it has proved to work better than ``PETSc``.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-hypreSolverType``      Allows to choose the solution method for the linear system if ``-poissonSolver`` is set to 
                             ``HYPRE``, discarded otherwise. Set to 1 to use the Generalized Minimum Residual (GMRES), set 
                             to 2 to use the preconditioned Conjugate-Gradient (PCG) method. Default value is 1.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-poissonTol``           Allows to set the exit tolerance for the pressure solver. Default value is 1e-8.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-poissonIt``            Set the maximum number of iterations for the pressure solver. Default value is 8.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-amgCoarsenType``       Since TOSCA uses the Algebraic Multi-Grid (AMG) preconditioner when the ``-poissonSolver`` is set 
                             to ``HYPRE``, this entry allows to set the coarsening method. Available entries are 0 (CLJP), 
                             6 (Falgout), 8 (PMIS), 10 (HMIS). Default value is 10.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-amgThresh``            Allows to set the AMG threshold. Default value is 0.5. For distorted meshes, a value of 0.6 is 
                             suggested.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-amgAgg``               Allows to set the level of aggressive coarsening. Default value is 0 (not used).
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-pTildeBuoyancy``       If set to 1, buoyancy force is recast into a buoyancy gradient and pressure is defined accordingly. 
                             Default value is 0 (not used).
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-dTdtScheme``           Can be set to ``backwardEuler`` (implicit first-order) or ``rungeKutta4`` (explicit fourth-order). 
                             For ABL simulations, ``backwardEuler`` is suggested. Crank-Nicholson has been removed from TOSCA
                             due to poor velocity-temperature coupling. 
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-relTolT``              Requires ``-dTdtScheme`` set to ``backwardEuler``. Allows to set the relative exit tolerance for 
                             the Newton method used to solve implicit discretized temperature equation, default value 1e-30.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-absTolT``              Requires ``-dTdtScheme`` set to ``backwardEuler``. Allows to set the absolute exit tolerance for 
                             the Newton method used to solve implicit discretized temperature equation, default value 1e-5.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-max_cs``               Maximum value for the LES model :math:`C_s` coefficient, default value is set to 0.5. Only used when
                             ``-les`` is greater than 1. 
   ========================= ====================================================================================================

Solution Constants 
******************

.. table:: 
   :widths: 35, 65
   :align: center 
         
   ========================= ====================================================================================================
   ``-nu``                   Sets the molecular (kinematic) viscosity of the working fluid.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-rho``                  Sets the density of the working fluid (used e.g. to compute forces).
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-Pr``                   Requires ``-potentialT`` to be se to 1. Sets the Prandtl number of the working fluid.
   ------------------------- ----------------------------------------------------------------------------------------------------
   ``-tRef``                 It is a required parameter when ``-potentialT`` is active and ``-abl`` is not. Sets the reference 
                             potential temperature of the flow, otherwise ``-tRef`` is set inside the ``ABLProperties.dat`` file.
   ========================= ====================================================================================================

Acquisition Controls 
********************

.. table:: 
   :widths: 35, 65
   :align: center 
     
   ============================== ========================================================================================================================
   ``-probes``                    Activates probes acquisition. Requires additional input files inside ``sampling/probes`` directory.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-sections``                  Activates acquisition of sections to be visualized in `ParaView`. Requires additional input files in 
                                  ``sampling/surfaces`` directory.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-averageABL``                Activates planar averages at every cell-level in the z-direction. Requires ``-abl`` to be active.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-averageABLPeriod``          Output period of the ABL planar averages. It is a required parameter, even if ``-averageABL`` is set to 
                                  0, for concurrent-precursor simulations, where these averages are always active. 
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-averageABLStartTime``       Time at which ABL planar averages are started. It is a required parameter, even if ``-averageABL`` is set to 
                                  0, for concurrent-precursor simulations.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-average3LM``                Activates vertical averages within layer at user-defined points. Requires additional inputs in 
                                  ``sampling`` directory.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-perturbABL``                Activates acquisition of perturbation fields at the same location as sections to be visualized in ParaView. 
                                  Requires additional inputs in ``sampling`` directory.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-averaging``                 It can be activated by setting to 1, 2, or 3 to get a higher amount of three-dimensional averaged fields.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-avgPeriod``                 Average period of three-dimensional averages. Fields are written at checkpoint times in the correspondent time 
                                  folder.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-avgStartTime``              Start time of three-dimensional averages.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-phaseAveraging``            These averages are a duplicate of the averages, but are useful if one wants to perform both unconditioned-averages 
                                  and phase-averages, e.g. at multiples of some characteristic time, in the same simulation.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-phaseAvgPeriod``            Average period of three-dimensional phase averages. Fields are written at checkpoint times in the correspondent time 
                                  folder.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-phaseAvgStartTime``         Start time of three-dimensional phase averages.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-keBudgets``                 Set to 1 to activate mechanical energy budgets. Requires additional inputs in ``sampling`` directory.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-writePressureForce``        Writes pressure force on the IBM surface.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-computeQ``                  Writes 3D field of Q-criterion at checkpoint times.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-computeL2``                 Writes 3D field of Lambda2-criterion at checkpoint times.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-computeFarmForce``          Writes 3D field of wind farm body force at checkpoint times.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-computeSources``            Compute source fields to be written in checkpoint files, to be converted in ``.xmf`` format by 
                                  ``tosca2PV`` executable. Depending on what is active, it calculates the coriolis force and driving 
                                  pressure gradient (require ``-abl`` set to 1), the inlet and lateral fringe source terms (require 
                                  ``-xDampingLayer`` and ``-yDampingLayer`` set to 1, respectively) and the body force from the canopy
                                  model (requires ``-canopy`` set to 1). 
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-computeBuoyancy``           Writes 3D field of buoyancy term in the momentum equation at checkpoint times.
   ============================== ========================================================================================================================

Post Processing Controls 
************************

.. table:: 
   :widths: 35, 65
   :align: center 
    
   ============================== ========================================================================================================================
   ``-postProcessFields``         Activate to post process 3D fields. It should be deactivated (set to 0) for too big cases to be fit in the memory of a 
                                  single node, as field conversion from binary to ``.xmf`` is not parallelized. Note that there is no plans within the 
                                  TOSCA developers team to parallelize this feature, as too big cases could not be visualized in `ParaView` anyways due 
                                  to RAM limitations. We suggest to use ``-sections`` instead, where parallel writing is enabled. 
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-writeRaster``               Activate to write raster file from jSections.
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-sections``                  Activate to post process binary sections and write ``.xmf`` and ``.hdf5`` files to be visualized in e.g. `ParaView`.
                                  This feature also works if ``tosca2PV`` is launched in parallel. When concurrent precursor is activated, those sections
                                  which fall inside the concurrent precursor domain are also processed, so that the user only has to provide one 
                                  section definition inside the ``sampling/surfaces`` directory. 
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-postProcessPrecursor``      Activate to also post process fields from the concurrent precursor simulation. Similarly to the ``-postProcessFields``
                                  flag, this option is not available in parallel. 
   ============================== ========================================================================================================================

`boundary`
~~~~~~~~~~       
    
This sections describes the input file definition used to prescribe boundary conditions in TOSCA. TOSCA expects a file for each solved field,
where boundary and initial conditions are specified. Pressure boundary conditions are not required, as they always are of the Neumann type, except when 
periodicity is present. Therefore, boundary conditions are expected for *U* (velocity) and *nut* (turbulent viscosity) fields, as well as for *T* (potential temperature), 
when this is activated in the *control.dat* file.
Specifically, TOSCA expects to find ``U``, ``nut`` and ``T`` (if applicable) files, inside a directory named ``boundary``. Within each file, boundary conditions 
have to be specified for each side of the domain. These are identified as ``kLeft``, ``kRight``, ``jLeft``, ``jRight``, ``iLeft`` and ``iRight``. In addition, 
the keyword ``internalField`` has to be specified, which sets the type of initial condition for each of these fields. An example of a ``boundary/U`` file is given below:

.. code-block:: C

   # TOSCA Input file - U boundary conditions 
   # ----------------------------------------

   internalField spreadInflow

   iLeft   periodic
   iRight  periodic
   jLeft   noSlip      
   jRight  slip
   kLeft   fixedValue (5.0 0.0 0.0)
   kRight  zeroGradient   

In addition to common boundary conditions, which are available at all boundaries, TOSCA also features some special boundary conditions, which are only 
available on specific boundaries. These are used for specific simulations and the user is suggested to consider such 
constraint when designing and setting up a TOSCA case. 

Within TOSCA, specific entry types are defined, which are used to define boundary conditions as well as other inputs throughout the code. These are 

1. *bool*       : an integer which is 0 (false) or 1 (true)
2. *integer*    : an integer number
3. *scalar*     : a floating point number (if an integer is provided, this is cast into a floating point number).
4. *vector*     : a group of floating point numbers defined as *(scalar scalar scalar)*, where the scalar **must** be a floating point number. There shouldn't be any space between parantheses and the first and last vector components. 
5. *string*     : a word without spaces. 
6. *dictionary* : a group of additional entries, where the number is dependant on the name of the dictionary, embodied within `{}` parentheses. Dictionaries can containmultiple integer, bool, scalar, vector and string entries. An example is given below. 
           
.. code-block:: C

   dictionary
   {
      boolEntry    1
      integerEntry 10
      scalarEntry  1.0 // or 1
      vectorEntry  (1.0 1.0 1.0)
      stringEntry  someRandomString
   }
   
Boundary Conditions
*******************

The boundary conditions available in TOSCA are summarized in the following table. The table also specifies the fields for which a given boundary condition is available, together with their entry type and description. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **BC type**                    **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``fixedValue``                 scalar or vector    Available for ``U``, ``nut`` and ``T``. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``fixedGradient``              scalar              Available for ``T``. Entry is a scalar value that refers to the gradient
                                                      normal to the patch.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``zeroGradient``               none                Available for ``U``, ``nut`` and ``T``.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``slip``                       none                Available for ``U``.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``noSlip``                     none                Available for ``U``.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``periodic``                   none                Available for ``U``, ``nut`` and ``T``. Requires keywords 
                                                      ``-iPeriodicType``, ``-jPeriodicType``, ``-kPeriodicType`` at the top of 
                                                      the mesh file, depending on the patch to which it is applied. These should 
                                                      be placed in the header, before the number of points in each direction, and
                                                      indicate the type of parallel connectivity between the periodic boundaries. 
                                                      For type 1, the domain is decomposed such that periodic cell-pairs lie 
                                                      within the same processor. Type 2 allows for a more general decomposition 
                                                      (suggested). 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``velocityWallFunction``       dictionary          Available only for ``U``. Applies wall models to velocity and temperature 
                                                      (if active). Only available for ``iLeft``, ``iRight``, ``jLeft`` and 
                                                      ``jRight`` patches. Not really tested on i-patches. It is usually applied 
                                                      on the ``jLeft`` patch in combination with an ``inletFunction`` on the 
                                                      ``kLeft`` patch (see below). 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``inletFunction``              dictionary          Available for ``U``, ``nut`` and ``T`` only at the ``kLeft`` patch. 
   ============================== =================== ============================================================================

The last two boundary conditions listed in the previous table are not available on all patches. On the one hand, the 
``velocityWallFunction``, used to prescribe wall models, is not available on ``kLeft`` and ``kRight`` boundaries. On the 
other hand, the ``inletFunction`` type, used to prescribe special boundary conditions for inlet velocity, temperature and turbulent viscosity, is only available on the ``kLeft`` patch. The user should keep this constraint in mind when setting
up a TOSCA case. 

Wall Models 
***********

The ``velocityWallFunction`` applies the Shumann wall model, used for ABL flows. In particular, the modeled wall shear 
stress is directly enforced in the stress term of the momentum equations, while a constant normal gradient condition 
is applied to the velocity. In order to avoid double counting at the wall when using this condition (the wall shear 
stress should be entirely and only modeled), the turbulent viscosity at the wall should be set to zero in the 
``boundary/nut`` file. The mandatory entries for this boundary condition type are defined as follows:

.. code-block:: C

   velocityWallFunction
   {
      type      -3        // Shumann wall model (only -3 is available)
      kRough    0.003     // equivalent roughness length
      gammaM    4.9       // Shumann model constant
      kappa     0.4       // von Karman constant
      thetaRef  300.0     // reference potential temperature 
      uStarEval averaged  // for laterally homogeneous flows (e.g. ABL), 
                          // otherwise set to 'localized' (e.g. wind farm flows)
   }
   
Inlet Functions 
***************

The ``inletFunction`` boundary condition has several different types, each corresponding to a different inlet boundary 
condition. These are summarized in the following table:

.. table:: 
   :widths: 15, 20, 65
   :align: center
   
   =========== =================================== ==================================================
   **IF type** **name**                            **description**
   ----------- ----------------------------------- --------------------------------------------------
   1           power law profile                   Power law velocity profile :math:`\mathbf{U}=
                                                   \mathbf{U}_\text{ref}\left(z/H_\text{ref}\right)
                                                   ^\alpha`, with :math:`\alpha=0.107027`.
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type      1
                                                         Uref      vector 
                                                         Href      scalar
                                                         uPrimeRMS scalar
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------
   2           logarithmic profile                 Logarithmic velocity profile :math:`\mathbf{U}=
                                                   u*/0.4\ln(z/z_0)\mathbf{e}_U`. 
                                                   :math:`\mathbf{U}` is constant above the 
                                                   inversion height :math:`H` and equal to 
                                                   :math:`\mathbf{U}(H)`.
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type       2
                                                         directionU vector 
                                                         hInversion scalar
                                                         frictionU  scalar
                                                         kRough     scalar
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------
   3           unsteady mapped inflow              Maps 2D planes contained in ``inflowDatabase``
                                                   directory at the inlet. Planes can be tiled 
                                                   laterally or vertically. Data can also be 
                                                   extrapolated in the vertical direction. The 
                                                   extrapolated value can be unsteady or steady. The 
                                                   latter is calculated by averaging the top data 
                                                   contained in the 2D planes and slowly merging with 
                                                   the unsteady data within the last 10 cells of the 
                                                   2D plane. Note: 2D planes and inlet mesh should be 
                                                   identical as data is only mapped.  
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type      3
                                                         n1Inflow  integer 
                                                         n2Inflow  integer 
                                                         n1Periods integer 
                                                         n2Periods integer 
                                                         n1Merge   bool
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------
   4           unsteady interpolated inflow        Same as type 3, but data is interpolated, hence 
                                                   the 2D plane mesh and inlet mesh can be different. 
                                                   Spanwise shift of the inflow data is also 
                                                   available by setting a lateral shift velocity. 
                                                   Note that this velocity is not added to the flow 
                                                   velocity, but rather data are shifted with this
                                                   constant lateral velocity. 
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type       4
                                                         n1Inflow   integer 
                                                         n2Inflow   integer 
                                                         n1Periods  integer 
                                                         n2Periods  integer
                                                         n1Merge    bool
                                                         n2Shift    bool
                                                         shiftSpeed scalar 
                                                         sourceType string
                                                         cellWidth1 integer 
                                                         cellWidth2 integer
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------
   5           Nieuwstadt (1983) model             Applies the Nieuwstadt (1983) model. This is more 
                                                   sophisticated than the logarithmic profile, as it
                                                   also contains wind veer. The wind profile is 
                                                   rotated such that the prescribed direction is 
                                                   imposed at the prescribed reference height.  
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type       5
                                                         directionU vector
                                                         hInversion scalar
                                                         hReference scalar
                                                         frictionU  scalar
                                                         kRough     scalar
                                                         latitude   scalar
                                                      }
                                                            
   ----------- ----------------------------------- --------------------------------------------------
   6           sinusoidally varying i-th component Uniform inflow, where the spanwise component 
                                                   varies sinusoidally with given amplitude and 
                                                   frequency. Useful to test turbine yaw controllers. 
                                                   
                                                   Usage:
                                                    
                                                   .. code-block:: C
                                                    
                                                      inletFunction
                                                      {
                                                         type      6
                                                         Uref      vector
                                                         amplitude scalar
                                                         periods   scalar
                                                      }
                                                      
   =========== =================================== ==================================================

The different entries required in the ``inletFunction`` dictionary for each function type are detailed below:

.. table:: 
   :widths: 25, 20, 55
   :align: center
   
   ======================== ================== =====================================================================================
   **entry**                **entry type**     **description**
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 1 - power law profile*   
   ---------------------------------------------------------------------------------------------------------------------------------
   ``Uref``                 vector             reference velocity in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``Href``                 scalar             height where :math:`\mathbf{U}= \mathbf{U}_\text{ref}` in m
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``uPrimeRMS``            scalar             isotropic fluctuation value in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 2 - logarithmic profile*
   ---------------------------------------------------------------------------------------------------------------------------------
   ``directionU``           vector             velocity direction (will be normalized)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``hInversion``           scalar             inversion height in m
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``frictionU``            scalar             friction velocity in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``kRough``               scalar             equivalent roughness length in m
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 3 - unsteady mapped inflow*
   ---------------------------------------------------------------------------------------------------------------------------------   
   ``n1Inflow``             integer            number of points in direction 1 (j for kLeft patch)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n2Inflow``             integer            number of points in direction 2 (i for kLeft patch)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n1Periods``            integer            number of points in direction 1 for tiling (if target is larger data is extrapolated)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n2Periods``            integer            number of points in direction 2 for tiling (if target is larger data is extrapolated)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n1Merge``              bool               average top cell values from 2D planes and merge within 10 top cells. Useful to make  
                                               the geostrophic region steady, removing e.g. inertial oscillations.
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 4 - unsteady interpolated inflow*
   ---------------------------------------------------------------------------------------------------------------------------------
   ``n1Inflow``             integer            number of points in direction 1 (j for kLeft patch)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n2Inflow``             integer            number of points in direction 2 (i for kLeft patch)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n1Periods``            integer            number of points in direction 1 for tiling (if target is larger data is extrapolated)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n2Periods``            integer            number of points in direction 2 for tiling (if target is larger data is extrapolated)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n1Merge``              bool               average top cell from 2D planes and merge within 10 top cells
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``n2Shift``              bool               spanwise shift of inflow data
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``shiftSpeed``           scalar             spanwise shift velocity in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``sourceType``           string             *uniform* (requires next two entries) or *grading* (requires the mesh file used in 
                                               the simulation that generated the 2D planes, renamed ``inflowMesh.xyz`` to be located
                                               inside the ``inflowDatabase`` directory)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``cellWidth1``           integer            inflow mesh cell width in direction 1 if ``sourceType`` *uniform*
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``cellWidth2``           integer            inflow mesh cell width in direction 2 if ``sourceType`` *uniform*
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 5 - Nieuwstadt (1983) model*
   ---------------------------------------------------------------------------------------------------------------------------------
   ``directionU``           vector             velocity direction (will be normalized)
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``hInversion``           scalar             inversion height in m
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``hReference``           scalar             reference height in m at which the wind has the provided direction
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``frictionU``            scalar             friction velocity in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``hRough``               scalar             equivalent roughness length in m
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``latitude``             scalar             latitude in degrees
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   *type 6 - sinusoidally varying i-th component*
   ---------------------------------------------------------------------------------------------------------------------------------
   ``Uref``                 vector             reference velocity in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``amplitude``            scalar             amplitude of the spanwise velocity oscillation in m/s
   ------------------------ ------------------ -------------------------------------------------------------------------------------
   ``periods``              scalar             number of periods contained in the domain length along the i curvilinear direction
   ======================== ================== =====================================================================================
      
As an example, a logarithmic profile on the ``kLeft`` patch is prescribed as follows

   .. code-block:: C

      kLeft inletFunction
            {
               type        2             // type 2 is logarithmic profile
               directionU  (1.0 0.0 0.0) // velocity is along the x direction
               hInversion  500           // ABL height is 500 m
               frictionU   0.5           // friction velocity is 0.5 m/s
               kRough      0.001         // equivalent roughness lenght is 0.001 m
            }

Initial Condition 
*****************
      
The initial field is finally prescribed after the ``internalField`` keyword within each boundary condition file. Initial conditions 
available within TOSCA are summarized in the following table:

.. table:: 
   :widths: 25, 75
   :align: center
   
   ==================================== ============================================================================
   **initial condition type**           **description**
   ------------------------------------ ----------------------------------------------------------------------------
   ``uniform``                          Available for ``U``, ``nu``, ``T``. The *value* entry can be vector or 
                                        scalar.  If the *perturbations* entry is set to 1, sinusoidal perturbations 
                                        are applied to trigger turbulence (only required for ``U`` field)
                                        
                                        Usage:
                                        
                                        .. code-block:: C
                                        
                                           uniform
                                           {
                                              value          vector or scalar
                                              perturbations  bool
                                           }
                                                                
   ------------------------------------ ----------------------------------------------------------------------------
   ``readField``                        Available for ``U``, ``nu``, ``T``. Reads field data from the ``fields`` 
                                        directory. Used for simulation restart. 
   ------------------------------------ ----------------------------------------------------------------------------
   ``ABLFlow``                          Sets initial log profile for ``U``, while ``T`` is set according to the 
                                        Rampanelli and Zardi (2003) model. Requires ``-abl`` set to 1 in the
                                        ``control.dat`` file, as it picks input parameters from the
                                        ``ABLProperties.dat`` file. Perturbations to trigger turbulence can be 
                                        added using the *perturbations* entry in the ``ABLProperties.dat`` file.
   ------------------------------------ ----------------------------------------------------------------------------
   ``spreadInflow``                     Copies the flow from the ``kLeft`` ghost cells at every ``k``-plane. Useful
                                        for ensuring perfect consistency between the inlet boundary condition and
                                        the internal field when using ``inletFunction``. In practice, sets the 
                                        internal field to whatever the ``inletFunction`` is. 
   ------------------------------------ ----------------------------------------------------------------------------
   ``linear``                           Only available for ``T``. Sets an initial temperature profile characterized
                                        by a linear lapse rate *tLapse* along ``j`` and a ground temperature 
                                        *tRef*. Used to prescribe an initially-linear potential temperature 
                                        stratification.
                                        
                                        Usage:
                                        
                                        .. code-block:: C
                                        
                                           linear
                                           {
                                              tRef   scalar 
                                              tLapse scalar 
                                           }
                                                                
   ==================================== ============================================================================













 
    
    
    
    
    
    
    
    
    
    
    
    
    
     
