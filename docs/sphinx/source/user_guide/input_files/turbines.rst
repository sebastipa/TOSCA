.. _turbines-section:

`turbines`
~~~~~~~~~~

This section describes how wind turbines can be activated in TOSCA. By setting the ``-windplant``
flag to 1 in the *control.dat* file, TOSCA is prompted to read turbine and farm parameters 
inside the ``turbines`` directory. The structure of a wind farm definition in TOSCA works as
follows. The list of all wind turbines present in the simulation is defined in the ``windFarmProperties``
file. Here, the user can define how frequently wind turbine variables are written to file, as well
as the type, model, location and activation of global wind farm controller for each wind 
turbine in the list. Such file has the following syntax:

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
and it will be saved in the ``postProcessing/turbines/startTime`` folder, which will be created at the start of the simulation. 
At each simulation restart, the ``startTime`` will differ, so TOSCA will not overwrite the turbine data. In the follwing table,
the entries required for each turbine in the ``windFarmProperties`` file are detailed. 

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
                                                      Inputs to the model are contained in the wind turbine model file. Available 
                                                      models are:
                                                      
                                                      - ``ALM``: actuator line model 
                                                      - ``ADM``: actuator disk model 
                                                      - ``uniformADM``: uniform actuator disk model 
                                                      - ``AFM``: actuator farm model 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``baseLocation``               vector              vector defining the coordinates of the base location.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``windFarmController``         bool                Wether or not the wind farm controller is active. If the controller is 
                                                      active for this turbine, an addition file is required, located in 
                                                      ``turbines/control``. This a look up table that should have the following 
                                                      format: 
                                                      
                                                      .. code-block:: C

                                                         # wind farm control table given as a function of time and pitch/DeltaCT 
                                                         time_1 pitch_1/DeltaCT_1
                                                         :      : 
                                                         time_n pitch_n/DeltaCT_n
                                                      
                                                      Notably, the file should contain only one header line. 
                                                      If ``turbineModel`` is set to ``ALM`` or ``ADM``, pitch in degrees should 
                                                      be provided in the second column, while for ``uniformADM`` and ``ADM``
                                                      the second column corresponds to a delta in thrust coefficient. Whether 
                                                      this is interpreted as disk based or freestream, it depends on the settings
                                                      provided in the wind turbine model file. This is a simple method to impose 
                                                      a wind farm controller (in the context of wind farm wake mixing) with no
                                                      feedback action. 
   ============================== =================== ============================================================================

As TOSCA requires individual definition of each turbine, this comes with some perks. In particular, the entires described in the 
above table can be different for each wind turbine. This means that the user can use different actuator models, controller 
activation and turbine model specification within the same simulation. This becomes useful for studies of multiple wind farm
clusters which might have different turbines and be controlled differently. 

The turbine type specification is contained in a file named as the ``turbineType`` entry in the ``windFarmProperties`` file. Its 
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

Notably, depending on the chosen actuator model type, different parameters are needed. Indication about this are provided in the 
above file example. Parameters which do not contain a comment line are always required. Each of the above parameters is detailed
in the table below:

.. table:: 
   :widths: 30, 20, 50
   :align: center
   
   ============================== =================== ============================================================================
   **entry**                      **entry type**      **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rTip``                 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rHub``               
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``hTower``               
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``overHang``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``precone``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``towerDir``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotorDir``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``upTilt``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``initialOmega``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nBlades``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``rotationDir``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``includeTower``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``includeNacelle``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``debug``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``genControllerType``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``pitchControllerType``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``yawControllerType``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nRadPts``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nAziPts``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``Uref``    
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``Ct``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``projection``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``sampleType``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_x``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_y``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilon_z``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_x``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_r``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``gaussexp_f``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_x``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_y``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``epsilonFactor_z``  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``towerData``   
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``nacelleData`` 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``airfoils`` 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``bladeData``                        
   ============================== =================== ============================================================================










 
