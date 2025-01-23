.. _boundary-subsection:

`boundary`
~~~~~~~~~~       
    
This section describes the input file definition used to prescribe boundary conditions in TOSCA. TOSCA expects a file for each solved field,
where boundary and initial conditions are specified. Pressure boundary conditions are not required, as they are always of the Neumann type except when 
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
available on specific boundaries. These are used for specific simulations, the user is suggested to consider such 
constraint when designing and setting up a TOSCA case. 
   
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
   ``velocityWallFunction``       dictionary          Available only for ``U``. Applies wall models to velocity. 
                                                      Only available for ``iLeft``, ``iRight``, ``jLeft`` and 
                                                      ``jRight`` patches. Not really tested on i-patches. It is usually applied 
                                                      on the ``jLeft`` patch in combination with an ``inletFunction`` on the 
                                                      ``kLeft`` patch. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``thetaWallFunction``          dictionary          Available only for ``T``. Applies wall models to temperature. 
                                                      Only available for ``iLeft``, ``iRight``, ``jLeft`` and 
                                                      ``jRight`` patches. Not really tested on i-patches. It is usually applied 
                                                      on the ``jLeft`` patch in combination with an ``inletFunction`` on the 
                                                      ``kLeft`` patch. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``inletFunction``              dictionary          Available for ``U``, ``nut`` and ``T`` only at the ``kLeft`` patch. 
   ============================== =================== ============================================================================

The last three boundary conditions listed in the previous table are not available on all patches. On the one hand, the 
``velocityWallFunction`` and ``thetaWallFunction``, used to prescribe wall models, are not available on ``kLeft`` and ``kRight`` boundaries. On the 
other hand, the ``inletFunction`` type, used to prescribe special inlet boundary conditions for inlet velocity, temperature and turbulent 
viscosity, is only available on the ``kLeft`` patch. The user should keep this constraint in mind when setting up a TOSCA case. 

Wall Functions 
**************

The ``velocityWallFunction`` applies the Shumann wall model, used for ABL flows, or the Cabot wall model, used for wall surfaces, to the velocoty field.  
In the Shumann model, the modeled wall shear stress is directly enforced in the stress term of the momentum equations, while a constant normal gradient condition 
is applied to the velocity. In order to avoid double counting at the wall when using this condition (the wall shear 
stress should be entirely and only modeled), the turbulent viscosity at the wall should be set to zero in the 
``boundary/nut`` file. Entries for velocity wall functions are summarized in the folllwing table:

.. table:: 
   :widths: 15, 20, 65
   :align: center
   
   =========== =================================== ================================================================================
   **WF type** **name**                            **description**
   ----------- ----------------------------------- --------------------------------------------------------------------------------
   -1          Cabot wall model                    Usage:
                                                    
                                                   .. code-block:: C

                                                      velocityWallFunction
                                                      {
                                                         type      -1        
                                                         kRough    0.003     // equivalent roughness length
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------------------------------------
   -3          Shumann wall model                  Usage:
                                                    
                                                   .. code-block:: C

                                                      velocityWallFunction
                                                      {
                                                         type      -3        
                                                         kRough    0.003     // equivalent roughness length
                                                         gammaM    4.9       // Shumann model constant
                                                         kappa     0.4       // von Karman constant
                                                         thetaRef  300.0     // reference potential 
                                                                             // temperature in Kelvin
                                                         uStarEval averaged  // for laterally homogeneous 
                                                                             // flows (e.g. ABL), otherwise 
                                                                             // set to 'localized' (e.g. 
                                                                             // wind farm flows)
                                                      }
   =========== =================================== ================================================================================
   
The ``thetaWallFunction`` is used to select different formulations of the Shumann wall model for potential temperature.  
In particular, the modeled heat flux is directly enforced in the stress term of the potential temperature equations, 
while a zero normal gradient condition is applied to the temperature field. In order to avoid double counting at the wall when using this condition, 
the turbulent viscosity at the wall should be set to zero in the ``boundary/nut`` file. 
Entries for temperature wall functions are summarized in the folllwing table:
   
.. table:: 
   :widths: 15, 20, 65
   :align: center
   
   =========== =================================== ================================================================================
   **WF type** **name**                            **description**
   ----------- ----------------------------------- --------------------------------------------------------------------------------
   -1          Shumann wall model,                 Usage:
               standard.                                     
                                                   .. code-block:: C

                                                      thetaWallFunction
                                                      {
                                                         type      -3        
                                                         kRough    0.003      // equivalent roughness length
                                                         gammaM    4.9        // Shumann model constant
                                                         gammaH    7.8        // Shumann model constant
                                                         alphaH    1.0        // Shumann model constant
                                                         kappa     0.4        // von Karman constant
                                                         thetaRef  300.0      // reference potential 
                                                                              // temperature in Kelvin
                                                         uStarEval averaged   // for laterally homogeneous 
                                                                              // flows (e.g. ABL), otherwise 
                                                                              // set to 'localized' (e.g. 
                                                                              // wind farm flows)
                                                         heatingRate -0.00014 // K/s
                                                      }
                                                      
   ----------- ----------------------------------- --------------------------------------------------------------------------------
   -3          Shumann wall model, with            Usage:
               specified heat flux.                                    
                                                   .. code-block:: C

                                                      thetaWallFunction
                                                      {
                                                         type      -2        
                                                         qWall     0.003     // wall heat flux in J/m2 
                                                      }
   ----------- ----------------------------------- --------------------------------------------------------------------------------
   -4          Shumann wall model, with            Only available at ``jLeft`` patch. 
               specified time history of           Usage:
               surface temperature and Obhukhov 
               length.                                    
                                                   .. code-block:: C

                                                      thetaWallFunction
                                                      {
                                                         type      -4        
                                                         kRough    0.003      // equivalent roughness length
                                                         gammaM    4.9        // Shumann model constant
                                                         gammaH    7.8        // Shumann model constant
                                                         alphaH    1.0        // Shumann model constant
                                                         kappa     0.4        // von Karman constant
                                                                              // temperature in Kelvin
                                                         uStarEval averaged   // for laterally homogeneous 
                                                                              // flows (e.g. ABL), otherwise 
                                                                              // set to 'localized' (e.g. 
                                                                              // wind farm flows)
                                                      }
                                                      
                                                   Look up tables of time (s), surface temperature (K), and Obhukhov length (m), 
                                                   stored in ``inflowDatabase/mesoscaleData/time``, 
                                                   ``inflowDatabase/mesoscaleData/surfTemp`` and 
                                                   ``inflowDatabase/mesoscaleData/L``, respectively. All vectors should have same 
                                                   size. 
   =========== =================================== ================================================================================
  
.. _inlet-functions-subsubsection:
   
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


