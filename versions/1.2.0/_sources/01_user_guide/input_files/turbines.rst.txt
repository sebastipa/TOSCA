.. _turbines-section:

`turbines`
~~~~~~~~~~

This section describes how wind turbines can be activated in TOSCA. By setting the ``-windplant``
flag to 1 in the *control.dat* file, TOSCA is prompted to read turbine and farm parameters 
inside the ``turbines`` directory. 

windFarmProperties 
******************

The structure of a wind farm definition in TOSCA works as
follows. The list of all wind turbines present in the simulation is defined in the `windFarmProperties`
file. Here, the user can define how frequently wind turbine variables are written to file, as well
as the type, model, location and activation of global wind farm controller for each wind 
turbine in the list. 

If TOSCA is coupled with OpenFAST, i.e. the compilation flag ``USE_OPENFAST`` has been set to 1 in the `makefile`,
the additional ``fastSettings`` dictionary is required in the `windFarmProperties` file. Notably, the ``nFastSubSteps``
entry defines the number of OpenFAST sub-iterations per TOSCA time step and should be consistent with the values provided in the
OpenFAST input file and TOSCA `control.dat` file. Currently, OpenFAST coupling is only supported with the actuator line model and 
when using a fixed time step.

The `windFarmProperties` file has in general the following syntax:

.. code-block:: C

   # TOSCA Input file - Wind Farm Properties
   # -------------------------------------------

   windFarmName            string
   arraySpecification      string
   debug                   bool

   writeSettings
   {
       timeStart           scalar
       intervalType        string
       timeInterval        scalar
   }

   fastSettings
   {
       nFastSubSteps       integer
       timeInterval        scalar
       restart             bool
   }

   turbineArray
   {
       turbineID_1
       (
           turbineType         string
           turbineModel        string
           baseLocation        vector
           windFarmController  bool
       )
       :
       :
       :
       turbineID_N
       (
           turbineType         string
           turbineModel        string
           baseLocation        vector
           windFarmController  bool
       )
   }
   
The ``windFarmName`` identifies the name of the wind farm, and it is totally up to the user to set it. Regarding ``arraySpecification``,
TOSCA currenty only supports the *onebyone* specification type, which is the most generic and allows to define all types of 
wind farms. In the future, we might add the *aligned* and *staggered* types, which do not require to specify each wind turbine
individually. Notably, these layouts can still be obtained with the current format. The ``debug`` flag prints to screen 
useful wind farm information at each iteration, but requires more processor communication when it is activated, so we suggest 
to set it to zero for production runs. The ``writeSettings`` dictionary allows to set parameters for the wind turbine acquisition
system. The ``timeStart`` entry, defined in seconds, defines the time at which TOSCA starts to write turbine data to file. The 
frequency at which this operation happens can be specified using the next two entries. Specifically, ``intervalType`` can be set 
to *timeStep* or *adjustableTime*. In the first case, turbine data are written to file every ``timeInterval`` iterations (which
has to be an integer), in the second case data are written every ``timeInterval`` seconds (which can be a scalar value with decimals). 
Notably, in this case the ``-adjustTimeStep`` flag has to be set to 1 in the *control.dat* file, and TOSCA will adjust the time 
step to land on multiples of the ``timeInterval`` specified by the user. 

Regarding the definition of the wind turbine list, this can be of arbitrary length. The name of each entry (e.g. *turbineID_1*, *turbineID_N*) 
is a string which must be specified by the user. This will be the name of the file where TOSCA will write each turbine's acquisition data,
and it will be saved in the ``postProcessing/turbines/startTime`` directory, which will be created at the start of the simulation. 
At each simulation restart, the ``startTime`` will differ, so TOSCA will not overwrite the turbine data. In the follwing table,
the entries required for each turbine in the `windFarmProperties` file are detailed. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``turbineType``                string              Name of the file defining the wind turbine model. Different wind turbines
                                                      can be present in the same wind farm. For example, if the wind turbine is
                                                      *ModelA*, a file with this name should be contained inside the ``turbines``
                                                      directory. See below for entries to this file. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``turbineModel``               string              Actuator model used to represent the wind turbine in the flow domain. 
                                                      Inputs to the model are contained in the wind turbine type file. Available 
                                                      models are:
                                                      
                                                      - ``ALM``: actuator line model 
                                                      - ``ADM``: actuator disk model 
                                                      - ``uniformADM``: uniform actuator disk model 
                                                      - ``AFM``: actuator farm model 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``baseLocation``               vector              vector defining the coordinates of the base location.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``windFarmController``         bool                Wether or not the wind farm controller is active. If the controller is 
                                                      active for this turbine, an additional file is required, located in 
                                                      ``turbines/control``. This a look up table that should have the following 
                                                      format: 
                                                      
                                                      .. code-block:: C

                                                         # wind farm control table given as a function of time and pitch/DeltaCT 
                                                         time_1 pitch_1/DeltaCT_1
                                                         :      : 
                                                         time_n pitch_n/DeltaCT_n
                                                      
                                                      Notably, the file should contain only one header line. 
                                                      If ``turbineModel`` is set to ``ALM`` or ``ADM``, pitch in degrees should 
                                                      be provided in the second column, while for ``uniformADM`` and ``AFM``
                                                      the second column corresponds to a delta in thrust coefficient. Whether 
                                                      this is interpreted as disk based or freestream, it depends on the settings
                                                      provided in the wind turbine type file. This is a simple method to impose 
                                                      a wind farm controller (in the context of wind farm wake mixing) with no
                                                      feedback action. 
   ============================== =================== ============================================================================

As TOSCA requires individual definition of each turbine, this comes with some perks. In particular, the entries described in the 
above table can be different for each wind turbine. This means that the user can use different actuator models, controller 
activation and turbine model specification within the same simulation. This becomes useful for studies of multiple wind farm
clusters which might have different turbines and be controlled differently. 

turbineTypeFile 
***************

The turbine type specification is contained in a file named as the ``turbineType`` entry in the `windFarmProperties` file. Its 
syntax is defined as follows
 
.. code-block:: C

   # TOSCA Input file - Turbine Type File
   # -------------------------------------------

   # Global wind turbine parameters
   rTip                scalar 
   rHub                scalar
   hTower              scalar
   overHang            scalar
   precone             scalar
   towerDir            vector
   rotorDir            vector
   upTilt              scalar 
   initialOmega        scalar  // for ALM and ADM only 
   nBlades             integer // for ALM and ADM only 
   rotationDir         string  // for ALM and ADM only 
   includeTower        bool
   includeNacelle      bool
   useOpenFAST         bool   // required if compiled with OpenFAST support
   debug               bool

   # Controllers 
   genControllerType   string  // for ALM and ADM only 
   pitchControllerType string  // for ALM and ADM only
   yawControllerType   string

   # Actuator model parameters
   nRadPts             integer // for ALM, ADM and uniformADM only 
   nAziPts             integer // for ADM and uniformADM only 
   Uref                scalar  
   Ct                  scalar  // for uniformADM and AFM only
   CtType              string  // for uniformADM and AFM only
   projection          string  // for ALM and AFM only 
   sampleType          string  // for ALM, uniformADM and AFM only 
   epsilon             scalar  // for ALM, ADM and unformADM only 
   epsilon_x           scalar  // for anisotropic AFM only 
   epsilon_y           scalar  // for anisotropic AFM only
   epsilon_z           scalar  // for anisotropic AFM only
   gaussexp_x          scalar  // for gaussexp AFM only
   gaussexp_r          scalar  // for gaussexp AFM only
   gaussexp_f          scalar  // for gaussexp AFM only
   epsilonFactor_x     scalar  // for anisotropic ALM only
   epsilonFactor_y     scalar  // for anisotropic ALM only
   epsilonFactor_z     scalar  // for anisotropic ALM only
   
   # Tower properties 
   towerData
   {
       Cd              scalar 
       epsilon         scalar 
       nLinPts         integer 
       rBase           scalar 
       rTop            scalar 
   }

   # Nacelle properties
   nacelleData
   {
       Cd              scalar 
       epsilon         scalar 
   }

   # names of the foils (for ALM and ADM only)
   airfoils
   {
       airfoil_1
       :
       airfoil_n
   }


   # Blade parameters given as:
   # (radius(m) c(m) twist(deg) thickness(m, anisotropic ALM only) airfoilID)  
   # for ALM and ADM only
   bladeData
   { 
       (scalar     scalar   scalar     integer)
       (:          :        :          :      )
       (scalar     scalar   scalar     integer)
   }

   # Ct table against a reference velocity given as:
   # (Uref(m/s) Ct)  
   # for UniformADM and AFM only with CtType set to variable
   CtTable
   { 
       (scalar     scalar)
       (:          :     )
       (scalar     scalar)
   }

Notably, depending on the chosen actuator model type, different parameters are needed. Indication about this are provided in the 
comments of the example file above. Parameters which do not contain a comment line are always required. Each of the above parameters is detailed
in the table below:

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rTip``                       scalar              turbine tip radius in m.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rHub``                       scalar              turbine hub radius in m. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``hTower``                     scalar              turbine tower height in m. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``overHang``                   scalar              distance between rotor center and tower top in m. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``precone``                    scalar              angle in degrees formed between turbine blades and the rotor plane. The 
                                                      rotor plane is defined as the unique plane that is formed when all blades 
                                                      lie in the same plane.
                                                      The more negative the angle, the furthest the blade tip from the tower.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``towerDir``                   vector              vector defining the tower direction from base to top. Normalized by TOSCA. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotorDir``                   vector              vector defining the direction of the rotor, facing the wind. Normalized by 
                                                      TOSCA. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``upTilt``                     scalar              angle in degrees formed between the normal of the rotor plane and the 
                                                      horizontal plane. The more positive the angle, the furthest the blade tip 
                                                      from the tower.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``initialOmega``               scalar              initial blade angular velocity in rpm. This is updated by the angular  
                                                      velocity controller if this is activated, otherwise it is kept constant. 
                                                      Only required when ``turbineModel`` is ``ALM`` or ``ADM``. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nBlades``                    integer             number of turbine blades. Only required when ``turbineModel`` is ``ALM`` 
                                                      or ``ADM``.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotationDir``                string              it can be set to *cw* (clockwise) or *ccw* (counter-clockwise). It is 
                                                      defined by looking at the wind turbine from upstream. Only required when 
                                                      ``turbineModel`` is ``ALM`` or ``ADM``.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``includeTower``               bool                whether or not to model also the tower as an actuator line. Requires the 
                                                      ``towerData`` dictionary. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``includeNacelle``             bool                whether or not to model also the nacelle as an actuator point. Requires the 
                                                      ``nacelleData`` dictionary. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``useOpenFAST``                bool                whether or not this turbine model is coupled with OpenFAST. Notably, a 
                                                      given turbine model can also not be coupled with OpenFAST. 
                                                      Only required when ``USE_OPENFAST`` flag is set to 1 in the `makefile`, i.e. 
                                                      TOSCA is compiled with OpenFAST support. When activated, the entries 
                                                      ``genControllerType``, ``pitchControllerType``, ``yawControllerType`` and 
                                                      tables ``airfoils``, ``bladeData`` and ``CtTable`` are not required, as this 
                                                      data is gathered from OpenFAST.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``debug``                      bool                debug switch. Prints turbine-level information and projection error. 
                                                      Deactivate for performance. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``genControllerType``          string              name of the file that contains the generator control properties of this 
                                                      wind turbine type. Should be located inside the ``turbines/control`` 
                                                      directory and it will be described later. 
                                                      Only required when ``turbineModel`` is ``ALM`` or ``ADM``. Set to 
                                                      *none* to disable.  
                                                      Set to *rpmControlCurve* to read a wind speed vs rpm curve from the file 
                                                      ``turbines/control/rotorRpmCurve``. This curve will be used to set
                                                      the optimal rotor speed as a function of the incoming wind speed, which 
                                                      is sampled 2.5 diameters upstream of the rotor disk. The file containing 
                                                      the curve should have the following format:
                                                      
                                                      .. code-block:: C
                                                      
                                                         rpmTable
                                                         { 
                                                             (scalar   scalar)
                                                             (:        :     )
                                                             (scalar   scalar)
                                                         }
                                                         
                                                      where the first column is the wind speed in m/s and the second column is
                                                      the rotor speed in rpm.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchControllerType``        string              name of the file that contains the collective pitch control properties of  
                                                      this wind turbine type. Should be located inside the ``turbines/control`` 
                                                      directory and it will be described later.
                                                      Only required when ``turbineModel`` is ``ALM`` or ``ADM``. Set to 
                                                      *none* to disable. 
                                                      Set to *bladePitchCurve* to read a wind speed vs pitch angle curve from the  
                                                      file ``turbines/control/bladePitchCurve``. This curve will be used to set
                                                      the optimal pitch angle as a function of the incoming wind speed, which 
                                                      is sampled 2.5 diameters upstream of the rotor disk. The file containing 
                                                      the curve should have the following format:

                                                      .. code-block:: C
                                                       
                                                            pitchTable
                                                            { 
                                                                (scalar   scalar)
                                                                (:        :     )
                                                                (scalar   scalar)
                                                            }
                                                             
                                                      where the first column is the wind speed in m/s and the second column is
                                                      the pitch angle in degrees. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawControllerType``          string              name of the file that contains the yaw control properties of  
                                                      this wind turbine type. Should be located inside the ``turbines/control`` 
                                                      directory and it will be described later. Set to 
                                                      *none* to disable.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``dipcControllerType``         string              name of the file that contains the collective dynamoc induction control 
                                                      properties of this wind turbine type. Should be located inside the 
                                                      ``turbines/control`` directory and it will be described later.
                                                      Only required when ``turbineModel`` is ``ALM`` or ``ADM``. Set to 
                                                      *none* to disable. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nRadPts``                    integer             number of radial points for the actuator model. Only required when 
                                                      ``turbineModel`` is ``ALM``, ``ADM`` or ``uniformADM``. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nAziPts``                    integer             number of azimuthal points for the actuator model. Only required when 
                                                      ``turbineModel`` is ``ADM`` or ``uniformADM``. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``Uref``                       scalar              freestream velocity for this wind turbine in m/s. The meaning of this  
                                                      parameter depends on the chosen actuator model. 
                                                      
                                                      - for ``ALM``, ``ADM`` and ``AFM`` it is used to compute the freestream   
                                                        thrust coefficient when writing the turbine data to file. It does not   
                                                        affect thrust computation.  
                                                      - for ``uniformADM`` it is used to compute thrust when ``sampleType`` is 
                                                        set to *givenVelocity*. For other sample types it is used to compute the 
                                                        freestream thrust coefficient when writing the turbine data to file and 
                                                        it does not affect thrust computation.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``Ct``                         scalar              thrust coefficient for this wind turbine. Only required when 
                                                      ``turbineModel`` is ``uniformADM`` or ``AFM``. The meaning of this parameter 
                                                      depends on the chosen actuator model. 
                                                      
                                                      - for ``uniformADM`` it should be the freestream thrust coefficient when
                                                        ``sampleType`` is set to *givenVelocity* or *rotorUpstream*. Conversely, 
                                                        when using *rotorDisk* it should be the disk-based thrust coefficient. 
                                                      - for ``AFM`` it should be the freestream thrust coefficient when
                                                        ``sampleType`` is set to *momentumTheory*. Conversely, when using 
                                                        *rotorDisk* or *integral* it should be the disk-based thrust coefficient.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``CtType``                     string              Ct type can be set to *constant* or *variable*.  Only required when 
                                                      ``turbineModel`` is ``uniformADM`` or ``AFM``. *constant* type requires a 
                                                      single thrust coefficient value. For *variable* type, a ``CtTable`` needs 
                                                      to be provided with the reference velocity and corresponding freestream 
                                                      thrust coefficient. This is regardless of ``sampleType`` used. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``projection``                 string              body force projection type. Only required for ``ALM`` and ``AFM`` models.
                                                      Entries specific to each model are described below.
   
                                                      ``ALM``
                                                       - *isotropic*: standard actuator line model.
                                                       - *anisotropic*: advanced actuator line model, where the gaussian projection 
                                                         function has different widths along the blade span, chord and thickness. 
                                                         The latter should be provided in the ``bladeData`` table. 
                                                      
                                                      ``AFM``
                                                       - *anisotropic*: anisotropic gaussian projection
                                                       - *gaussexp*: gauss-exponential projection. The latter has been subject to 
                                                         more validation, hence is suggested. 
                                           
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``sampleType``                 string              velocity sampling type. Only required for ``ALM``, ``uniformADM`` and 
                                                      ``AFM`` models. Entries specific to each model are described below.
                                                      
                                                      ``ALM``
                                                       - *rotorDisk*: tri-linearly interpolates velocity at actuator 
                                                         points 
                                                       - *integral*: integrates the velocity around the blade, weighted by 
                                                         the projection function. 
                                                        
                                                      ``uniformADM``
                                                       - *rotorUpstream*: averages velocity on a fictitious 
                                                         actuator disk, located 2.5 diamters upstream of the wind turbine 
                                                       - *givenVelocity*: does not sample any velocity but uses the provided value 
                                                         of ``Uref`` 
                                                       - *rotorDisk*: averages velocity on the rotor disk. 
                                                       
                                                      ``AFM``
                                                       - *rotorDisk*: tri-linearly interpolates velocity at the rotor 
                                                         center
                                                       - *integral*: integrates velocity around the rotor center, weighted  
                                                         by the projection function
                                                       - *momentumTheory*: computes the freestream 
                                                         velocity as :math:`U_{disk}/(1-a)`, where :math:`a` is the induction 
                                                         factor, defined as :math:`\frac{1 - \sqrt{1-C_T}}{2}`.
                                                        
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon``                    scalar              standard deviation of the projection function in m for ``ALM`` (when
                                                      ``projection`` is set to *isotropic*), ``ADM`` and ``uniformADM``. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_x``                  scalar              standard deviation of the projection function in along x for ``AFM`` 
                                                      (when ``projection`` is set to *anisotropic*).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_y``                  scalar              standard deviation of the projection function in along y for ``AFM`` 
                                                      (when ``projection`` is set to *anisotropic*).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_z``                  scalar              standard deviation of the projection function in along z for ``AFM`` 
                                                      (when ``projection`` is set to *anisotropic*).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_x``                 scalar              standard deviation of the projection function in along x for ``AFM`` 
                                                      (when ``projection`` is set to *gaussexp*).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_r``                 scalar              radial helf-decay length in m of the exponential projection function for 
                                                      ``AFM`` model, when ``projection`` is set to *gaussexp*. Usually set equal  
                                                      to the turbine radius. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_f``                 scalar              smoothing parameter of the exponential projection function for ``AFM`` 
                                                      model, when ``projection`` is set to *gaussexp*.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_x``            scalar              chord multiplication factor defining the standard deviation along the blade  
                                                      chord of the gaussian projection function for ``ALM`` when 
                                                      ``projection`` is set to *anisotropic*.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_y``            scalar              thickness multiplication factor defining the standard deviation along the  
                                                      blade thickness of the gaussian projection function for ``ALM`` when 
                                                      ``projection`` is set to *anisotropic*.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_z``            scalar              :math:`dr` multiplication factor defining the standard deviation along the  
                                                      blade radius of the gaussian projection function for ``ALM`` when 
                                                      ``projection`` is set to *anisotropic*.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``towerData``                  dictionary          Dictionary defining the parameters required to represent the tower as an 
                                                      actuator line. Usage:
                                                      
                                                      .. code-block:: C
                                                      
                                                         towerData
                                                         {
                                                             Cd      scalar  // drag coefficient 
                                                             epsilon scalar  // projection length 
                                                                             // in m
                                                             nLinPts integer // num. of actuator 
                                                                             // points
                                                             rBase   scalar  // tower base radius 
                                                             rTop    scalar  // tower top radius 
                                                         }
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nacelleData``                dictionary          Dictionary defining the parameters required to represent the nacelle as an 
                                                      actuator point. Usage:
                                                      
                                                      .. code-block:: C
                                                      
                                                         nacelleData
                                                         {
                                                             Cd      scalar // drag coefficient 
                                                             epsilon scalar // projection length 
                                                                            // in m
                                                         }
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``airfoils``                   table               List of files storing the blade airfoil information. The list has an 
                                                      arbitrary length files named as the element of the list should be located 
                                                      inside ``turbines/airfoils``. Usage:
                                                      
                                                      .. code-block:: C
                                                      
                                                         airfoils
                                                         {
                                                             airfoil_1
                                                             :
                                                             airfoil_n
                                                         }
                                                         
                                                      Each file should have the following format: 
                                                      
                                                      .. code-block:: C
                                                      
                                                         # data provided as | alpha | C_l | C_d

                                                         table
                                                         {
                                                              (scalar    scalar   scalar)
                                                              (:         :        :     )
                                                              (scalar    scalar   scalar)
                                                         }
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bladeData``                  table               List of properties for each blade radial station. The list has an arbitrary 
                                                      length and should store radial position, chord, twist and airfoil ID. 
                                                      The latter is the index of one of the airfoils present in the airfoils list,
                                                      starting from 0. Specifically, that which corresponds to the given radial 
                                                      position. When using ``ALM`` with ``projection`` set to *anisotropic*, 
                                                      blade thickness is added to the list, before the airfoil ID. Usage:
                                                      
                                                      .. code-block:: C
                                                      
                                                         bladeData
                                                         { 
                                                           (scalar scalar scalar integer)
                                                           (:      :      :      :      )
                                                           (scalar scalar scalar integer)
                                                         }
                                                         
                                                      or 
                                                      
                                                      .. code-block:: C
                                                      
                                                         bladeData
                                                         { 
                                                           (scalar scalar scalar scalar integer)
                                                           (:      :      :      :      :      )
                                                           (scalar scalar scalar scalar integer)
                                                         }
                                                         
                                                      when using ``ALM`` with ``projection`` set to *anisotropic*. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``CtTable``                    table               Table of reference velocity and corresponding free stream veocity based 
                                                      thrust coefficient Ct. Usage:
                                                      
                                                      .. code-block:: C
                                                      
                                                         CtTable
                                                         { 
                                                           (scalar scalar)
                                                           (:      :     )
                                                           (scalar scalar)
                                                         }
                                                         
                                                      Required only for ``UniformADM`` and ``AFM`` models when ``CtType`` set to
                                                      variable.                                                   
   ============================== =================== ============================================================================

control 
*******

In order to complete the definition of a given wind turbine, information regarding individual turbine control should be 
provided. While yaw control can be applied to any actuator model, angular velocity and pitch control are only available 
for ``ALM`` and ``ADM``, as these feature a more detailed description of the wind turbine. The name of the turbine control 
files that TOSCA will use for each specific wind turbine are defined with the entries ``genControllerType``, ``pitchControllerType``, 
``yawControllerType`` and ``dipcControllerType``, detailed in the previous table. The general idea is that a given controller is entirely described in a file. 
Hence, in most cases these entries will be equal if all three controllers are active on a given wind turbine. However, a single 
turbine can be also controlled with different controllers, defined in different files. Moreover, as many different controllers 
as desired by the user can be used within a given wind farm. The turbine control file should be contained inside the 
``turbines/control`` directory, and consists of the the entires described in the following table. 

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``genInertia``                 scalar              rotational inertia of the electrical generator in kg/m2.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``hubInertia``                 scalar              rotational inertia of the hub in kg/m2.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bldIntertia``                scalar              rotational inertia of each blade in kg/m2.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gbxRatioG2R``                scalar              gearbox ratio (generator / rotor).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gbxEfficiency``              scalar              gearbox efficiency.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``genEff``                     scalar              generator efficiency.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``initialGenTorque``           scalar              initial torque of the electrical generator.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``genTqControllerParameters``  dictionary          parameters defining the generator torque controller.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchControllerParameters``  dictionary          parameters defining the collective pitch controller. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawControllerParameters``    dictionary          parameters defining the nacelle yaw controller. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``dipcControllerParameters``   dictionary          parameters defining the dynamic induction controller.
   ============================== =================== ============================================================================

Entries for the last four dictionaries are described in the following table

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *genTqControllerParameters* 
   -------------------------------------------------------------------------------------------------------------------------------
   ``genSpeedFilterFreq``         scalar              corner frequency of generator-speed low-pass filter in Hz.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``cutInGenSpeed``              scalar              transitional generator speed between regions 1 and 1.5 in rpm.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``cutInGenTorque``             scalar              generator torque at cut-in wind speed in Nm. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``regTwoStartGenSpeed``        scalar              transitional generator speed between regions 1.5 and 2 in rpm.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``regTwoEndGenSpeed``          scalar              transitional generator speed between regions 2 and 3 in rpm.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``ratedGenTorque``             scalar              generator torque at rated wind speed in Nm. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``controllerPGain``            scalar              proportional gain of torque controller in Nm/rpm2
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``torqueRateLimiter``          bool                activation of limit on rate of torque variation. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rtrSpeedLimiter``            bool                activation of limit on rotor speed. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``torqueMaxRate``              scalar              maximum rate of torque variation in Nm/s.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``ratedRotorSpeed``            scalar              rotor speed at rated wind speed in rpm. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *pitchControllerParameters* 
   -------------------------------------------------------------------------------------------------------------------------------
   ``pitchRateLimiter``           bool                activation of limit on rate of pitch variation. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchAngleLimiter``          bool                activation of limit on pitch angle. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchMaxRate``               scalar              maximum rate of pitch variation in deg/s. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchMin``                   scalar              minimum pitch value in deg. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchMax``                   scalar              maximum pitch value in deg. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``controllerPGain``            scalar              proportional gain of pitch controller in s. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``controllerIGain``            scalar              integral gain of pitch controller. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``controllerDGain``            scalar              derivative gain of pitch controller. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchS2R``                   scalar              blade pitch angle in deg at which the rotor power has doubled. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *yawControllerParameters* 
   -------------------------------------------------------------------------------------------------------------------------------
   ``sampleType``                 string              type of velocity sampling to calculate the incoming wind angle. Currently
                                                      the only possibility is *hubUpDistance*, where the wind is sampled one 
                                                      rotor diameter upstream with respect to the rotor center. We are planning
                                                      to implement also *anemometer* sampling, where the wind is sampled from the 
                                                      nacelle-mounted anemometer (right now this is just a place-holder and it 
                                                      does not do anything).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``avgWindow``                  scalar              averaging window to filter the wind measurement at the sampling location.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawMin``                     scalar              minimum yaw angle in deg (range is -180 to 180).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawMax``                     scalar              maximum yaw angle in deg (range is -180 to 180).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawSpeed``                   scalar              constant nacelle rotational speed in deg/s. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``allowedError``               scalar              allowed misalignment error in deg under which yaw is not changed. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``initialFlowAngle``           scalar              initial angle in deg between the incoming flow and the direction of the 
                                                      rotor plane with zero upTilt. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
    *dipcControllerParameters*
   -------------------------------------------------------------------------------------------------------------------------------
   ``dipcHelixAmp``               scalar              amplitude of the helicoidal perturbation in degrees added to the calculated 
                                                      airfoil angle of attack.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``dipcHelixDir``               string              direction of the helicoidal perturbation as viewed from upstream. 
                                                      Possible entries are *cw* (clockwise) and *ccw* (counter clockwise).
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``dipcHelixFreq``              scalar              frequency of the helicoidal perturbation in Hz.
   ============================== =================== ============================================================================

OpenFAST Coupling
*****************

When coupling TOSCA with OpenFAST, the user must provide an OpenFAST input file for each wind turbine listed in the 
`windFarmProperties` file. The OpenFAST input file is unique for each turbine and should be located inside
the ``turbines`` directory. Moreover, it should be named as ``<turbineType>.ID.fst``, where ``<turbineType>`` is the value of the 
``turbineType`` entry in the `windFarmProperties` file, for the given wind turbine, and ``ID`` is the integer 
identifying each wind turbine in the list, starting from 0. Conversely, OpenFAST module files (such as AeroDyn, ElastoDyn etc.) can be shared among different 
wind turbines, and their path is specified within each of the OpenFAST input files. 

When a given turbine is coupled with OpenFAST, TOSCA's outputs - for example rotor speed, rotor direction etc. - are derived from OpenFAST 
and not calculated by TOSCA. However, in order to integrate with previous TOSCA functionalities and workflows, turbine output data are still written 
to file by TOSCA, as if TOSCA was calculating them. These outputs are coincident with the OpenFAST outputs.

Restart is also possible when coupling with OpenFAST. However, all or no turbines should be restarted by activating the global ``restart`` flag inside 
the ``fastSettings`` dictionary of the `windFarmProperties` file. When restarting, TOSCA will read the OpenFAST restart files named as 
``<turbineType>.ID.iteration.chkp``, located inside the ``turbines`` directory. The ``iteration`` number corresponds to the iteration number at which OpenFAST 
restart is performed. These files are created by OpenFAST according to the settings specified in the OpenFAST input file.

To familiarize with TOSCA-OpenFAST coupling, the ``tests/NREL5MWOpenFastTest`` is provided as an example case in this repository and consisting of the 
``tests/NREL5MWTest`` example case, with the required modifications to feature OpenFAST coupling.



 
