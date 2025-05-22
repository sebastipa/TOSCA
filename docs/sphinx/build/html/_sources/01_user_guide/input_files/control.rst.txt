.. _control-subsection:

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
   ``-les``               Set to 1 to activate the LES model. If set to 0, the simulation runs in DNS mode.
   ---------------------- -----------------------------------------------------------------------------------------------------

   ``-lesModel``          Required when ``-les`` is set to 1. It specifies the LES model to be used. The available entries are:
                          
                          - *smagorinsky*: standard Smagorinsky 
                          - *stabilityDependent*: stability dependent
                          - *dynamicSmagorinsky*: dynamic Smagorinsky with box averaging
                          - *dynamicLASI*: dynamic Smagorinsky scale-invariant with lagrangian averaging
                          - *dynamicLASD*: dynamic Smagorinsky scale-dependent with lagrangian averaging
                          - *dynamicPASD*: dynamic Smagorinsky scale-dependent with planar averaging
                          - *amd*: anisotropic minimum dissipation 

                          TOSCA has been extensively and most used adopting the *dynamicLASI* model for ABL and wind plant 
                          simulations. Other models have been recently added to the code and might perform better in some 
                          cases (especially *dynamicLASD* and *amd*).
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
   ``-advectionDampingX`` Specifies if x-direction horizontal inlet advection damping region is present in the simulation. 
                          Requires additional input in ``ABLProperties.dat`` file if activated. 
   ---------------------- -----------------------------------------------------------------------------------------------------
   ``-advectionDampingY`` Specifies if y-direction horizontal inlet advection damping region is present in the simulation. 
                          Requires additional input in ``ABLProperties.dat`` file if activated.
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
   ``-avgABLPeriod``              Output period of the ABL planar averages. It is a required parameter, even if ``-averageABL`` is set to 
                                  0, for concurrent-precursor simulations, where these averages are always active. 
   ------------------------------ ------------------------------------------------------------------------------------------------------------------------
   ``-avgABLStartTime``           Time at which ABL planar averages are started. It is a required parameter, even if ``-averageABL`` is set to 
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

