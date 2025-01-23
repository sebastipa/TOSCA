.. _examples_neutral_abl_test:

Neutral ABL Example
-------------------

This example case shows how to simulate a conventionally neutral atmospheric boundary 
layer (CNBL) in TOSCA. CNBLs are characterized by a strong capping inversion layer of 
thickness :math:`\Delta h`, across which a potential temperature jump :math:`\Delta \theta`
is observed. Stratification is neutral below the inversion layer, while it is 
characterized by a linear lapse rate :math:`\gamma` aloft. 

These types of flows are simulated using lateral periodic boundary conditions, 
with the flow driven by a large-scale driving pressure gradient imposed using TOSCA's 
velocity controller. In this example, the pressure controller will be used, which controls
the large-scale horizontal pressure gradient in order to obtain a wind speed of 8 m/s at
the height of 90 m. 

The flow is initialized using a logarithmic profile below the inversion height :math:`H`,
and a uniform wind aloft. As the simulation progressess, turbulence is created and starts 
developing inside the flow domain. To accelerate this process, spanwise fluctuations can be 
imposed in the initial flow to anticipate the turbulence startup. 

All the above aspects of the simulations can be controlled from the ``ABLProperties.dat`` file 
(see Sec. :ref:`ablProperties-section` for a detailed description of all entries). 
In particular, entries controlling the initial temperature and velocity profile are defined as 
follows: 

.. code-block:: bash

   uRef    8.0    # wind speed at hRef in m/s
   hRef    90.0   # reference height in m
   tRef    300.0  # ground potential temperature in K
   hRough  0.003  # equivalent roughness length in m 
   hInv    750.0  # height of the inversion center in m
   dInv    100.0  # width of the inversion center in m 
   gInv    5.0    # potential temperature jump in K across inversion 
   gTop    0.003  # gradient of the linear lapse rate above inversion in K/m
   gABL    0.0    # gradient of the linear lapse rate below inversion in K/m (zero for neutral)
   vkConst 0.4    # value of the Von Karman constant 
   smearT  0.33   # defines the smearing of the initial potential temperature profile
   
Then, the following flags are set in order to activate the velocity controller, the temperature controller 
(which maintains constant the horizontally averaged potential temperature profile), the Coriolis force and the 
initial flow perturbations in the velocity field: 

.. code-block:: bash

   controllerActive  1               # activates velocity controller 
   controllerActiveT 1               # activates temperature controller 
   controllerTypeT   initial         # target temperature is the initial one 
   coriolisActive    1               # activates Coriolis force 
   fCoriolis         5.156303966e-5  # 2 times the Coriolis parameter
   perturbations     1               # flow perturbations 
   
Notably, the ``fCoriolis`` should be set to :math:`7.272205217\cdot 10^{-5} \sin(\phi)`, where :math:`\phi` is 
the latitude. Finally, the parameters defining the velocity controller must be specified, namely 

.. code-block:: bash

   controllerProperties
   {
       controllerAction        write        # controller writes the momentumSource inside postProcessing/
       controllerType          pressure     # controller tries to maintain uRef at hRef
       relaxPI                 0.7          # proportional gain 
       alphaPI                 0.9          # proportional over integral ratio 
       timeWindowPI            3600         # time constant of the controller low-pass filter 
       geostrophicDamping      1            # only for pressure controller, eliminates inertial oscillations 
       geoDampingAlpha         1.0          # alpha 1 corresponds to critical damping 
       geoDampingStartTime     5000         # should start damping after turbulence has developed 
       geoDampingTimeWindow    4500         # should be large (ideally one inertial oscillation)
       controllerMaxHeight     100000.0     # the controller is applied in the entire domain 
   }
   
In order to tell TOSCA to initialize the flow using the parameters defined within the ``ABLProperties.dat`` file, the 
``internalField`` keyword inside the boundary condition files should be set to *ABLFlow*. In order to check that 
everything is correct, we run TOSCA one second and check the resulting fields in *ParaView*. In order to 
do so, we set 

.. code-block:: bash

   -startFrom                 startTime
   -startTime                 0
   -endTime                   1
   -adjustTimeStep            0
   -timeStep                  0.5
   
   -timeInterval              1
   
   -sections                  0
   
in the ``control.dat`` file. This will stop the simulation after 2 iterations of 0.5 s each (i.e. 1 s of simulation). 
Notably, the flow sections are also deactivated, as two seconds of simulation will not produce any. 
Once TOSCA executables ``tosca`` and ``tosca2PV`` are copied inside the case directory can be run in serial using 

.. code-block:: bash

   ./tosca 
   
and in parallel (using 4 processors) with 

.. code-block:: bash

   mpirun -np 4 ./tosca 

Once the initial test completes, the user should first check that the entry ``-postProcessFields`` is active in the ``control.dat`` file, then run 

.. code-block:: bash

   ./tosca2PV
   
The last executable will create an ``XMF`` directory in which the solution fields will be stored in *.hdf* format. The 
log file should look as follows:

.. image:: ./images/neutral_abl_test_log_1.png
    :width: 100%

In order to visualize the data, the user should navigate inside the ``XMF`` directory and open the *.xmf* file it e.g. in *ParaView*. 
The following image shows the result of this operation. In particular, velocity and potential temperature fields are depicted on 
the left and right panels, respectively. 

.. image:: ./images/neutral_abl_test_initial_field.png
    :width: 100%

Now the spinup phase of the simulation can be conducted. This is required to reach statistical independence of the 
turbulent fluctuations and, although it does not produce any meaningful data, it is required before one can start to 
gather flow statistics and sampling sections. In particular, the following parameters are set in the 
``control.dat`` file:

.. code-block:: bash

   -startFrom                 startTime
   -startTime                 0
   -endTime                   10000
   -adjustTimeStep            1
   
   -timeInterval              500
   
   -sections                  0
   -probes                    0
   -averageABL                1
   -average3LM                0
   -avgABLPeriod              10
   -avgABLStartTime           0

   -averaging                 0
   -avgPeriod                 50 
   -avgStartTime              1000

The simulation will restart from time 0 and run for 10000 s. The time step will be adjusted based on the CFL value and in 
order to average ABL fields every 10 s. The user can notice that every 500 s a new time directory is created inside the 
``fields`` directory, in which checkpoint files are saved. Moreover, ABL averages are stored every 10 s inside the 
``postProcessing/averaging/startTime`` directory. The figure below shows the velocity (left) and temperature fields (left)
at the end of the spinup phase.

.. image:: ./images/neutral_abl_test_spinup_field.png
    :width: 100%

Now the simulation can be restarted from 10000 s and flow statistics can be gathered. If one wish to use flow sections from 
this simulation as inlet boundary condition for a subsequent simulation (e.g. with wind turbines) sampling planes should be 
activated. These are also useful for large cases when the domain is to big to be visualized in its entirety. 
Regarding sampling planes intended to be later used as inlet boundary conditions (*kSections*), they have to be written with the highest 
frequency possible. As a consequence, the ``sampling/surfaces/kSections`` will have the following syntax:

.. code-block:: bash
   
   surfaceNumber 1         # one k-surface
   timeStart     10000     # acquisition starts at 10000
   intervalType  timeStep  # writes every timeInterval iterations 
   timeInterval  1         # writes every iteration
   coordinates   10.0      # plane located at x = 10 m 
   
while sections used for visualization (usually *iSections* and *jSections*) are be defined as

.. code-block:: bash
   
   surfaceNumber 2               # two j-surface
   timeStart     1000            # acquisition starts at 10000
   intervalType  adjustableTime  # writes every timeInterval seconds
   timeInterval  10              # writes every 10 seconds 
   coordinates   90.0  500.0     # flow section at hub and inversion heights

In order for TOSCA not to restart the initial field using the ABL parameters, but instead reading the latest time 
available in the fields directory and running until 12000 s, one can set in the ``control.dat`` file

.. code-block:: bash
   
   -startFrom    latestTime      # reads the latest time inside the fields directory 
   -endTime      12000           # simulate 2000 more s with data acquisition
   
and change the ``internalField`` entry in the boundary files to ``readField`` (in all ``U``, ``p`` and ``nut``). 
Finally, 3D field averaging and the sampling sections are activated by setting

.. code-block:: bash
   
   -sections                  1
   -averageABL                1
   -avgABLPeriod              10
   -avgABLStartTime           0

   -averaging                 1
   -avgPeriod                 5 
   -avgStartTime              10000

in the ``control.dat`` file. Once TOSCA is submitted, the user can verify from the log file (see below) that *kSections* are stored
every iteration inside the ``postProcessing/kSurfaces/k-ID`` directory, where ``k-ID`` corresponds to the k-index of the
plane of cells that the selected surface intercepts. Notably, the name of each section file corresponds to the time at which 
the section was taken, where the number of decimal places can be controlled with the ``-timePrecision`` keyword in the ``control.dat`` file. 

.. image:: ./images/neutral_abl_test_log_2.png
    :width: 100%
   
   
   
