.. _acquisition-section: 

Data Acquisition System 
-----------------------

TOSCA writes 3D fields at each checkpoint time inside the ``fields`` directory. However, when running 
large simulations, these checkpoint files can be very heavy and can not really be used for post processing. 
In these cases, their main purpose is to be used only for restarting the simulation from one of the save times. 
Moreover, the utility that converts them into a file format which can be visualized (e.g. in *ParaView*),
``tosca2PV``, only works in serial when the flag ``-postProcessFields`` is active in the ``control.dat`` file.
Notably, even if ``tosca2PV`` would be able to post process large 3D fields in parallel, these could 
still not be visualized, as visualization softwares do not work with distributed memory and so the data would not 
fit in the available RAM memory of most architectures. 

Hence, visualizing 3D field is only intended for debugging of small cases and should not be included in the 
workflow of large production runs. 

For these types of cases, TOSCA offers a powerful data acquisition system (DAS) with many different functionalities
to gather data in a format that can be easily post processed by means of user-defined routines in Matlab or 
python. Specifically, the following different types of data acquisition can be included in a production 
TOSCA worflow:

 - sections: user-defined planes that can be converted in a format which can be visualized 
   in e.g. *ParaView*. Available data to be saved are velocity, pressure, potential temperature and effective 
   viscosity. 
 - probes: user-defined points where the time history of velocity, potential temperature and pressure can be saved. 
 - 3D averages and phase averages: 3D fields averages that are saved to file at every checkpoint time. These can 
   be visualized entirely (in serial) or can be sliced by ``tosca2PV`` after the simulation, and visualized in *ParaView*.
 - mechanical energy budgets: these budgets can be calculated in user-defined boxes (e.g. around each wind turbine and in 
   their wake). Notably, TOSCA is not energy conservative so the various terms do not add up to zero. 
 - ABL statistics: only works for cartesian meshes. Writes the time history of horizontally averaged fields to 
   file. This is very useful when running precursor simulations of the ABL. Assumes horizontally-homogeneous flow. 
 - three-layer model variables: computes depth-averaged 3LM fields on a user-defined horizonta grid. This utility was
   implemented to compare LES with the reduced-order three-layer model, but it has not been used much otherwise. 
 - flow perturbations: computes 2D planes, to be visualized in paraview, of gravity waves perturbations. This is very 
   useful for wind farm simulation featuring gravity waves and gravity wave analysis. 
 - turbine acquisition: time series of different turbine variables, for each turbine in the computational domain. 
 - IBM forces: writes pressure and viscous forces on IBM surfaces. 
 - additional fields: serveral additional 3D fields can be written to file for visualization (in serial),
   mainly for debugging.  
   
For standard ABL applications, setions, ABL statistics and 3D averages are suggested (probes can be added for 
velocity correlations in time and space), while for wind farm simulations sections, probes, averages, turbine data
and flow perturbations are suggested. Each individual data acquisition system is detailed in the following sections. 

.. _sections-subsection:

`-sections`
~~~~~~~~~~~ 

Settings related to section acquisition should be contained inside the ``sampling/surfaces`` directory. These are of two types 
 
 - curvilinear sections
 - user-defined sections
 
Curvilinear sections are really meant to be used when the mesh is cartesian in the surface-parallel directions. Otherwise, they become a curvilinear surface. 

Curvilinear sections are defined by specifying the cartesian coordinate in one of the curvilinear directions (*k,j,i*). At this point, the curvilinear
index which defines the surface is selected by looking for the cell with the specified cartesian coordinate using TOSCA convention (*k is x*, *j is z* and *i is y*), 
while the remaining coordinates are set to be in the domain center point. Then, the index in the direction normal to the section is selected, and all cells having the same 
index are selected to form the surface. As can be understood, if the mesh is cartesian in the surface-parallel directions, the surface 
is a plane. Conversely, for curvilinear meshes the surface follows the curvilinear directions, as it is selected through mesh indexing rather 
than coordinates. For these reason, these sections are faster to write.  

Conversely, user-defined sections are defined by providing the cartesian coordinates of each point, hence they can be arbitrary surfaces
within the computational domain and data is tri-linearly interpolated at the surface points. User-defined sections can be saved during the simulation
similar to curvilinear section or processed aposteriori by the utility ``tosca2PV``, in the latter case they are used to slice averaged fields after the simulation 
has completed. Hence, ``-averaging`` should be activated during the run in order to produce the average data to be later sliced. 

The workflow to save sections in TOSCA is as follows: there is the option to save instantaneous sections during the simulation or averaged sections after. 
If the ``-averaging`` flag is activated in the ``control.dat`` file during the simulation, average fields are produced as explained 
in Sec. :ref:`averaging-subsection`. These fields can be sliced after the simulation by ``tosca2PV``, both by curvilinear and user-defined sections.  

TOSCA only saves the slices in binary format, these are then converted to *hdf5* format, which can be visualized in *ParaView*, by running 
``tosca2PV``. If, after converting the slices, the user wishes to save additional slices of the average fields at different locations, the section files 
contained inside the ``sampling/surfaces`` directory can be edited, and ``tosca2PV`` should be run again to produce the new slices. If any slices are already present 
inside the ``postProcessing/iSurfaces``, ``postProcessing/jSurfaces`` or ``postProcessing/kSurfaces`` directories,
``tosca2PV`` will not produce the new slices, so the user should rename the slices saved during the simulation by renaming these directories 
to something else, depending if ``sampling/surfaces/iSections``, ``sampling/surfaces/jSections`` or ``sampling/surfaces/kSections`` have been edited, respectively. 
This procedure is identified in TOSCA as on-the-fly section re-acquisition, and it is printed in the terminal (or the log file) when running ``tosca2PV``.

The classification of each different section type available in TOSCA is provided in the following table: 

.. table:: 
   :widths: 16, 17, 67
   :align: center
   
   ============================== =================== ============================================================================
   **name**                       **type**            **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *kSections*                    k-normal section    Automatically saves velocity, pressure, effective viscosity and potential 
                                                      temperature (if applicable). Defined in ``sampling/surfaces/kSections`` 
                                                      file.  
                                                      Usage: 
                                                      
                                                      .. code-block:: C
                                                      
                                                         surfaceNumber 2        // number of surfaces 
                                                         timeStart     0        // start of acquisition in s
                                                         intervalType  timeStep // or adjustableTime 
                                                         timeInterval  1        // iterations or s based 
                                                                                // on above entry
                                                         coordinates   0 500    // list of x coordinates
                                                         
                                                      Note: no comments should be present in this file. Above comments are 
                                                      only for explanation in the context of the user guide. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *jSections*                    j-normal section    Automatically saves velocity, pressure, effective viscosity and potential 
                                                      temperature (if applicable). Defined in ``sampling/surfaces/jSections`` 
                                                      file. 
                                                      Usage: 
                                                      
                                                      .. code-block:: C
                                                      
                                                         surfaceNumber 1              // number of surfaces 
                                                         timeStart     0              // start of acquisition 
                                                                                      // in s
                                                         intervalType  adjustableTime // or timeStep 
                                                         timeInterval  1              // iterations or s 
                                                         coordinates   90             // list of z coordinates
                                                         
                                                      Note: no comments should be present in this file. Above comments are 
                                                      only for explanation in the context of the user guide. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *iSections*                    i-normal section    Automatically saves velocity, pressure, effective viscosity and potential 
                                                      temperature (if applicable). Defined in ``sampling/surfaces/iSections`` 
                                                      file. 
                                                      Usage: 
                                                      
                                                      .. code-block:: C
                                                      
                                                         surfaceNumber 2        // number of surfaces 
                                                         timeStart     0        // start of acquisition in s
                                                         intervalType  timeStep // or adjustableTime 
                                                         timeInterval  1        // iterations or s based 
                                                                                // on above entry
                                                         coordinates   0 500    // list of y coordinates
                                                         
                                                      Note: no comments should be present in this file. Above comments are 
                                                      only for explanation in the context of the user guide. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   *userSections*                 user-defined        Requires ``-averaging`` set to 1 during the simulation, then sections are 
                                                      produced when running ``tosca2PV``. Automatically saves mean velocity, 
                                                      pressure, effective viscosity and potential 
                                                      temperature (if applicable). Defined in ``sampling/surfaces/userSections``
                                                      directory. This should contain **ONLY** files where sections are defined,  
                                                      the name can be arbitrary and **ALL** files are read by TOSCA. The file 
                                                      syntax is as follows: 
                                                      
                                                      .. code-block:: C
                                                      
                                                         timeStart      0        // start of acquisition in s
                                                         intervalType   timeStep // or adjustableTime 
                                                         timeInterval   1        // iterations or s based 
                                                                                 // on above entry
                                                         flipIndexOrder 0        // if 1 flips ny and nx order
                                                         meshPoints     ny nx
                                                         x_0 y_0 z_0
                                                         :   :   : 
                                                         x_n y_n z_n
                                                         
                                                      Note: no comments should be present in this file. Above comments are 
                                                      only for explanation in the context of the user guide. Total number of 
                                                      points should be ny times nx. For each ny, all nx are read, i.e. 
                                                      coordinates are read as follows:
                                                      
                                                      .. code-block:: C
                                                      
                                                         for(j=0; j<ny; j++)
                                                         {
                                                             for(i=0; i<nx; i++)
                                                             {
                                                                  // read x 
                                                                  // read y
                                                                  // read z
                                                             }
                                                         }
                                                      
                                                      Hence, the file should be created accordingly. 
   ============================== =================== ============================================================================   


.. _probes-subsection:

`-probes`
~~~~~~~~~

Probes acquisition is probably the best acqusition utility in TOSCA and it is fully parallelized. In fact, 
probes in the domain can be defined by an arbitrary number of files, and each file can contain multiple probes. 
All probes contained in each file are controlled by a different group of processors, defined as the ensamle 
of processors that own the cells where the probes in the file are contained. Hence, all probe files can be ideally 
updated simultaneously if they are defined in a clever way. For example, each file may contain an array 
of probes at a given x and y location, which only varies in z. In this manner, only a few processors control 
this given probe rake, and different rakes will be likely controlled by other processors, making the writing 
very efficient. Conversley, if probes are split among files casually, there might be some processor groups 
that feature a large number of processors, while others might only have a few, making the writing more unbalanced. 
Even worse, defining all probes in one file renders the writing un-parallelized, since a single group of processor is 
created for probe acquisition. 

When overset mesh is used (see Sec. :ref:`overset-section`), a given probe may be contained in two domains, i.e. the parent 
and the child. In this case, TOSCA automatically assigns the probe to the child domain, which is ideally that with the 
higher resolution among the two. 

Probe files are contained inside the ``sampling/probes`` directory. Each file contains the definition of a given number of
probe coordinates. Notably, the ``probes`` directory should **only** contain probe files, as all files inside this directory 
are read assuming to be probe files. Probe files should have the following syntax:

.. code-block:: C
                                                      
   probesNumber integer  
   timeStart    scalar 
   intervalType string
   timeInterval scalar or integer 
   fields       string

   locations

   scalar scalar scalar 
   :      :      : 
   scalar scalar scalar 
   
The various entries are detailed in the following table:

.. table:: 
   :widths: 20, 17, 63
   :align: center
   
   ============================== =================== ============================================================================
   **name**                       **type**            **description**
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``probesNumber``               integer             defines the number of entries that are read after the ``locations`` keyword,
                                                      hence the number of probes in the file.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``timeStart``                  scalar              time at which the probe acquisition starts.     
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``intervalType``               string              if it is set to *adjustableTime*, then ``timeInterval`` is in seconds, 
                                                      while if it is set to *timeStep* then ``timeInterval`` indicates the number 
                                                      of interations. In the former case, if ``-adjustTimeStep`` is set to 1 in 
                                                      the ``control.dat`` file, then the time step is adjusted to land on 
                                                      ``timeInterval`` multiples of the simulation start time when writing probes.  
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``timeInterval``               scalar or integer   probes acquisition frequency. 
                                                      It is a scalar in s if ``intervalType`` is set to *adjustableTime*, while it 
                                                      is an integer, indicating the number of iterations, if ``intervalType`` is 
                                                      set to *timeStep*. 
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``fields``                     string              specifies the fields that are sampled at the probe locations. 
                                                      Only velocity, pressure and potential temperature can be sampled. 
                                                      This should be a string **with no spaces** that selects all or some of these
                                                      fields, e.g. *U,T,p* for all, or *T,U* for potential temperature and 
                                                      velocity only.
   ------------------------------ ------------------- ----------------------------------------------------------------------------
   ``locations``                  table               specifies the coordinates of the probe locations as a list of three scalars
                                                      separated by one blank line from the ``locations`` keyword, where the 
                                                      scalars indicate x, y and z coordinates, respectively. 
                                                      There should be as many lines as the value of the ``probesNumber`` 
                                                      keyword. 
   ============================== =================== ============================================================================

Probe data are written inside the ``postProcessing/<probeName>/<startTime>`` directory, where ``<probeName>`` is the name of the 
file where a given set of probes are defined inside ``sampling/probes``, and ``<startTime>`` is the start time of the simulation. If the simulation is 
restarted, a new ``<startTime>`` directory will be created, so that data is not overwritten. 

.. _averaging-subsection:

`-averaging` and `-phaseAveraging`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

The averaging functionalities in TOSCA are activated by setting ``-averaging`` and ``-phaseAveraging`` to a number greater than 0 in the 
``control.dat`` file. The former is used to compute 3D averages, while the latter is used to compute phase averages. Phase averages are a duplicated 
of the average fields, useful when the user wants to co,pute both averages and phase averages in the same simulation. 
Averages are accumulated during the simulation and written to disk at every checkpoint time. They can be visualized in *ParaView*
by running ``tosca2PV`` in serial (entire 3D field), or in parallel (on-the-fly sections described in Sec. :ref:`sections-subsection`).
The various fields that are accumulated by TOSCA based on the values of ``-averaging`` and ``-phaseAveraging`` are detailed in the following table:

.. table:: 
   :widths: 30, 70
   :align: center
   
   ============================== =================================================================================================
   **Output Field**                      **Description**
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgU``                       averaged velocity :math:`\overline{u_i}`. 
                                  Activated with ``-averaging`` greater than 0 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgP``                       averaged pressure :math:`\overline{p}`. 
                                  Activated with ``-averaging`` greater than 0 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgUU``                      averaged resolved Reynolds stress tensor components :math:`\overline{u'_i u'_j}`. 
                                  Activated with ``-averaging`` greater than 0 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgDUU``                     averaged resolved dispersive stress tensor components :math:`\overline{u''_i u''_j}`. 
                                  Activated with ``-averaging`` greater than 0 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgOmega``                   averaged vorticity :math:`\overline{\omega_i}`.
                                  Activated with ``-averaging`` greater than 1 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgPsq``                     averaged pressure squared :math:`\overline{p^2}`.
                                  Activated with ``-averaging`` greater than 2 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgUdotGradP``               averaged dot product between velocity and pressure gradient 
                                  :math:`\overline{u_i \partial p / \partial x_i}`.
                                  Activated with ``-averaging`` greater than 2 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgMagGradU``                averaged square root of two times the strain rate tensor 
                                  :math:`\overline{\sqrt{2 \partial u_i / \partial x_j \partial u_i / \partial x_j}}`.
                                  Activated with ``-averaging`` greater than 2 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgMagUU``                   averaged velocity scaled by its module squared :math:`\overline{(u_i u_i)u_j}`.
                                  Activated with ``-averaging`` greater than 2 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgOmegaOmega``              averaged vorticity fuctuation tensor :math:`\overline{\omega'_i \omega'_j}`.
                                  Activated with ``-averaging`` greater than 2 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgNut``                     average effective viscosity :math:`\overline{\nu_t}`. Activated with ``-averaging`` greater than 
                                  0 and  ``-les`` greater than 1 in the ``control.dat`` file.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``avgCs``                      averaged SGS model coefficient :math:`\overline{C_s}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgU``                     duplicate of ``avgU``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgP``                     duplicate of ``avgP``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgUU``                    duplicate of ``avgUU``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgOmega``                 duplicate of ``avgOmega``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgPsq``                   duplicate of ``avgPsq``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgUdotGradP``             duplicate of ``avgUdotGradP``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgMagGradU``              duplicate of ``avgMagGradU``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgMagUU``                 duplicate of ``avgMagUU``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgOmegaOmega``            duplicate of ``avgOmegaOmega``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgNut``                   duplicate of ``avgNut``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``phAvgCs``                    duplicate of ``avgCs``.
                                  Controlled by ``-phaseAveraging`` in the ``control.dat`` file instead of ``-averaging``.
   ============================== =================================================================================================

Notably, when the simulation has to be restarted and averaging has been performed in the previous run, TOSCA automatically takes into 
account the average cumulation that has been previously performed up to the checkpoint time if average fields have been written in the
``fields/<checkpointTime>`` directory at checkpoint. In fact, weighting info are saved in the ``fields/<checkpointTime>/fieldInfo`` file.    

.. _mke-budgets-subsection:

`-keBudgets`
~~~~~~~~~~~~

The ``-keBudgets`` functionality in TOSCA allows for the computation and storage of mean kinetic energy (MKE) budgets within user-defined boxes. 
These budgets can be calculated in Cartesian components or general curvilinear coordinates (GCC). The results are saved at every checkpoint 
time and can be used for detailed analysis of energy transfer and dissipation within the specified regions. MKE budgets are activated by setting 
``-keBudgets`` to 1 in the ``control.dat`` file, and require the additional ``sampling/keBudgets`` file which is characterized by the following
syntax:

.. code-block:: C

   cartesian    bool         // 1 to perform the budget in Cartesian components, 0 for GCC
   debug        bool         // 1 to activate debug mode, 0 to deactivate
   avgStartTime scalar       // Start time for averaging
   avgPeriod    scalar       // Averaging period

   boxArray
   {
       string                // box 1 name
       {
           center  vector    // Center of the box
           sizeX   scalar    // Size of the box in the x-direction
           sizeY   scalar    // Size of the box in the y-direction
           sizeZ   scalar    // Size of the box in the z-direction
       }
       string                // box 2 name
       {
           center  vector    // Center of the box
           sizeX   scalar    // Size of the box in the x-direction
           sizeY   scalar    // Size of the box in the y-direction
           sizeZ   scalar    // Size of the box in the z-direction
       }
       :
       string                // box n name
       {
           center  vector    // Center of the box
           sizeX   scalar    // Size of the box in the x-direction
           sizeY   scalar    // Size of the box in the y-direction
           sizeZ   scalar    // Size of the box in the z-direction
       }
   }

The MKE budget data is written, one file per box, inside the ``postProcessing/keBoxes/<startTime>`` directory, where ``<startTime>`` is the start time 
of the simulation. If the simulation is restarted, a new ``<startTime>`` directory will be created, so that data is not overwritten.
The various fields that are accumulated by TOSCA, representing the different terms of the mean kinetic energy equation within each box,
either as flux terms on the faces or as volume terms, are detailed in the following table:

.. table:: 
   :widths: 20, 70
   :align: center
   
   ============================== =================================================================================================
   **Output Field**               **Description**
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``name``                       Name of the box.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``center``                     Center coordinates of the box.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``size``                       Size of the box in x, y, and z directions.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dum``                        Sum of mean kinetic energy advection in x, y, and z directions.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dup``                        Sum of turbulent kinetic energy advection in x, y, and z directions.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dpm``                        Sum of mean pressure advection in x, y, and z directions.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dpp``                        Sum of fluctuating pressure advection in x, y, and z directions.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dumx``                       Mean kinetic energy advection in the x direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dupx``                       Turbulent kinetic energy advection in the x direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dpmx``                       Mean pressure advection in the x direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dppx``                       Fluctuating pressure advection in the x direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dumy``                       Mean kinetic energy advection in the y direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dupy``                       Turbulent kinetic energy advection in the y direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dpmy``                       Mean pressure advection in the y direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dppy``                       Fluctuating pressure advection in the y direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dumz``                       Mean kinetic energy advection in the z direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dupz``                       Turbulent kinetic energy advection in the z direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dpmz``                       Mean pressure advection in the z direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Dppz``                       Fluctuating pressure advection in the z direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Fx``                         Turbulent flux in the x direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Fy``                         Turbulent flux in the y direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Fz``                         Turbulent flux in the z direction.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``eps``                        Turbulent dissipation.
   ============================== =================================================================================================

Moreover, the following fields are saved at every checkpoint and can be visualized  with ``tosca2PV`` (this can only be done in serial
and it has been included for debugging purposes):

.. table:: 
   :widths: 20, 70
   :align: center
   
   ============================== =================================================================================================
   **Output Field**               **Description**
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keErr``                      Error on the MKE equation.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keEm``                       Mechanical energy.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keD``                        Dissipation term.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keF``                        Turbulent fluxes.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keEps``                      Turbulent dissipation.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keUpUp``                     Reynolds stress tensor components.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``keFarm``                     Wind turbine power (if wind farm is active).
   ============================== =================================================================================================

.. _averaging-abl-subsection:

`-averageABL`
~~~~~~~~~~~~~

The ``-averageABL`` functionality in TOSCA allows for the computation and storage of horizontally-averaged fields within the atmospheric boundary layer (ABL).
These fields are computed at each curvilinear level in the j direction. Given these restrictions, it only makes sense to activate this functionality when **the mesh is cartesian** 
and **the flow is horizontally homogeneous**. Moreover, the simulation should be set up such that the j direction is the vertical direction. 
The ABL statistics are activated by setting ``-averageABL`` to 1 in the ``control.dat`` file. Moreover, the user is required to set the
``-avgABLPeriod`` and the ``-avgABLStartTime`` parameters in the ``control.dat`` file, which control the frequency at which the ABL statistics are written to disk and 
the time at which they start to be written, respectively. Statistics are written inside the ``postProcessing/averaging/<startTime>`` directory, where ``<startTime>`` is the start time
of the simulation. They consist of different files, one for each field, with the following format: 

.. code-block:: C

   time_0  timeStep_0 stat_level_0 stat_level_1 .. stat_level_n
   time_1  timeStep_1 stat_level_0 stat_level_1 .. stat_level_n
   :       :          :            :            : 
   time_n  timeStep_n stat_level_0 stat_level_1 .. stat_level_n

where *stat* is the field name, and *level* is the curvilinear level in the j direction. In addition, a file called 
``hLevelsCell`` is written in the same directory, which contains the mesh levels in meters. 

The various fields that are accumulated by TOSCA are detailed in the following table:

.. table:: 
   :widths: 20, 70
   :align: center
   
   ============================== =================================================================================================
   **Output Field**               **Description**
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``U_mean``                     averaged x-component velocity :math:`\overline{u}`. 
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``V_mean``                     averaged y-component velocity :math:`\overline{v}`. 
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``W_mean``                     averaged z-component velocity :math:`\overline{w}`.  
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``nu_SGS_mean``                averaged effective viscosity :math:`\overline{\nu}_t`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``uu_mean``                    averaged resolved Reynolds stress tensor 11 component :math:`\overline{u'u'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``vv_mean``                    averaged resolved Reynolds stress tensor 22 component :math:`\overline{v'v'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``ww_mean``                    averaged resolved Reynolds stress tensor 33 component :math:`\overline{w'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``uv_mean``                    averaged resolved Reynolds stress tensor 12 component :math:`\overline{u'v'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``uw_mean``                    averaged resolved Reynolds stress tensor 13 component :math:`\overline{u'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``vw_mean``                    averaged resolved Reynolds stress tensor 23 component :math:`\overline{v'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R11_mean``                   averaged SGS Reynolds stress tensor 11 component :math:`\overline{\tau}_{SGS,11}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R22_mean``                   averaged SGS Reynolds stress tensor 22 component :math:`\overline{\tau}_{SGS,22}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R33_mean``                   averaged SGS Reynolds stress tensor 33 component :math:`\overline{\tau}_{SGS,33}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R12_mean``                   averaged SGS Reynolds stress tensor 12 component :math:`\overline{\tau}_{SGS,12}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R13_mean``                   averaged SGS Reynolds stress tensor 13 component :math:`\overline{\tau}_{SGS,13}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``R23_mean``                   averaged SGS Reynolds stress tensor 23 component :math:`\overline{\tau}_{SGS,23}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``T_mean``                     averaged potential temperature :math:`\overline{\theta}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``q1_mean``                    averaged resolved heat flux 1 component :math:`\overline{q}_1`.  
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``q2_mean``                    averaged resolved heat flux 2 component :math:`\overline{q}_2`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``q3_mean``                    averaged resolved heat flux 3 component :math:`\overline{q}_3`.  
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Tu_mean``                    averaged fluctuation potential temperature times x-velocity velocity fluctuation vector
                                  :math:`\overline{\theta'u'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------    
   ``Tv_mean``                    averaged fluctuation potential temperature times y-velocity velocity fluctuation vector
                                  :math:`\overline{\theta'v'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``Tw_mean``                    averaged fluctuation potential temperature times z-velocity velocity fluctuation vector
                                  :math:`\overline{\theta'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------  
   ``wuu_mean``                   averaged resolved triple velocity fluctiation tensor 311 component :math:`\overline{w'u'u'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``wvv_mean``                   averaged resolved triple velocity fluctiation tensor 322 component :math:`\overline{w'v'v'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``www_mean``                   averaged resolved triple velocity fluctiation tensor 333 component :math:`\overline{w'w'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``wuv_mean``                   averaged resolved triple velocity fluctiation tensor 312 component :math:`\overline{w'u'v'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``wuw_mean``                   averaged resolved triple velocity fluctiation tensor 313 component :math:`\overline{w'u'w'}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``wvw_mean``                   averaged resolved triple velocity fluctiation tensor 323 component :math:`\overline{w'v'w'}`.
   ============================== =================================================================================================

Notably, the averaging process happens horizontally only, hence the above mentioned files actually contain a time series of horizontally-averaged 
quantities. The user can then use these files to compute further statistics, such as boundary layer height history, friction velocity history, 
capping inversion width history, etc. 

.. _averaging-3lm-subsection:

`-average3LM`
~~~~~~~~~~~~~

The ``-average3LM`` functionality in TOSCA allows for the computation and storage of depth-averaged fields within the atmospheric boundary layer (ABL).
Specifically, depth-averaged pressure, velocity and boundary layer height are calculated on a uniform 2D user-defined grid. This utility can be useful 
when comparing LES data against the reduced-order three-layer model (3LM) for the ABL. The 3LM is a simplified model that describes the ABL as a
three-layer system, with two layers inside the ABL and one layer aloft, and it is used to study the interaction between wind farms and atmopsheric gravity waves. 
The 3LM field computation is activated by setting ``-average3LM`` to 1 in the ``control.dat`` file. This prompts TOSCA to read the ``sampling/3LM`` file, 
which contains the 2D mesh definition, the sampling frequency and start time, as well as the definition of the three layers. Specifically, the file
should have the following syntax:

.. code-block:: C

   nstw      scalar        // number of sampling points in the streamwise direction
   nspw      scalar        // number of sampling points in the spanwise direction
   upDir     vector        // upward direction vector
   streamDir vector        // streamwise direction vector

   avgStartTime scalar     // start time for averaging
   avgPeriod    scalar     // averaging period

   level1 
   {
        hStart     scalar  // starting height of the first layer 
        hEnd       scalar  // ending height of the first layer
   }

    level2 
    {
          hStart     scalar  // starting height of the second layer 
          hEnd       scalar  // ending height of the second layer
    }

    level3 
    {
          hStart     scalar  // starting height of the third layer 
          hEnd       scalar  // ending height of the third layer
    }

3LM statistics are written inside the ``postProcessing/3LM/<startTime>`` directory, where ``<startTime>`` is the start time of the simulation. 
They consist of different files, one for each field. Each file has three header lines, one blank line and then the data, which is organized as 
a matric with nstw rows and nspw columns and the value of the field at each element. The various fields that are accumulated by TOSCA are detailed
in the following table:

.. table:: 
   :widths: 20, 70
   :align: center
   
   ============================== =================================================================================================
   **Output Field**               **Description**
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``IBL``                        wind farm boundary layer height :math:`h_{IBL}`. 
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``LBL``                        capping inversion layer bottom height :math:`h_{LBL}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``TBL``                        capping inversion layer top height :math:`h_{TBL}`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``P_L1``                       depth-averaged pressure in the first layer :math:`\overline{p}_1`. 
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``P_L2``                       depth-averaged pressure in the second layer :math:`\overline{p}_2`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``P_L3``                       depth-averaged pressure in the third layer :math:`\overline{p}_3`.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``U_L1``                       depth-averaged velocity in the first layer :math:`\overline{u}_1`. Each element of the matrix
                                  data is a TOSCA vector.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``U_L2``                       depth-averaged velocity in the second layer :math:`\overline{u}_2`. Each element of the matrix
                                  data is a TOSCA vector.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``U_L3``                       depth-averaged velocity in the third layer :math:`\overline{u}_3`. Each element of the matrix
                                  data is a TOSCA vector.
   ------------------------------ -------------------------------------------------------------------------------------------------
   ``points``                     2D mesh points. Each element of the matrix data is a TOSCA vector.
   ============================== =================================================================================================


.. _perturb-abl-subsection:

`-perturbABL`
~~~~~~~~~~~~~

The ``-perturbABL`` functionality in TOSCA allows for the computation and storage of gravity wave perturbations (or in general any other perturbation field) within the atmospheric boundary layer (ABL).
ABL perturbations of a generic variable :math:`\phi` are calculated as follows:

.. math:: 

   \phi_p(x,y,z) = \overline{\phi}(x,y,z) - \overline{\phi}(x_s,y,z)

where :math:`\overline{\phi}` is the temporal mean and :math:`\overline{\phi}(x_s,y,z)` is the vertical plane of :math:`\overline{\phi}`, normal to the streamwise direction and at a given streamwise location :math:`x_s`. 
Notably, the location at which the latter is evaluated should be representative of non-perturbed conditions. This utility can be useful when studying the interaction between wind farms and atmopsheric gravity waves, 
as it allows to quantify and visualize the perturbations to which the ABL is subjected to. This functionality is activated by setting ``-perturbABL`` to 1 in the ``control.dat`` file. This prompts TOSCA to 
read the ``sampling/perturbABL`` file, which has the following syntax:

.. code-block:: C

   avgStartTime scalar     // start time for averaging
   avgPeriod    scalar     // averaging period
   xSample      scalar     // streamwise location at which the non-perturbed profile is evaluated

Perturbation 3D fields ``UpABL`` (velocity perturbation), ``PpABL`` (pressure perturbation) and ``TpABL`` (potential temperature perturbation) are written inside the ``fields/<checkpointTime>`` directory. 
these fields can be either later visualized in serial with ``tosca2PV``. For large cases this is not possible, and ``-perturbABL`` should be used in conjunction with ``-sections`` to save slices of the perturbation fields
at the same locations where instantaneous sections are saved. This is again done by the ``tosca2PV`` utility, so the user does not have to decide where to slice the perturbation fields at simulation runtime. 

.. _turbine-data-subsection:

`turbine data`
~~~~~~~~~~~~~~

Wind turbine data acquisition is always active in TOSCA once the wind farm is activated. The user can specify the acquisition settings in the ``writeSettings`` dictionary of the ``windFarmProperties`` file (see Sec. :ref:`turbines-section`).
A file is written for each turbine inside the  ``postProcessing/turbines/<startTime>`` directory, where ``<startTime>`` is the start time of the simulation. `` directory. The file name matches the wind turbine ID defined in the ``windFarmProperties`` file. 
The amount of data that is written in each file depends on the actuator model adopted for the wind turbine. Each wind turbine output file has as many rows as the number of time samples, and as many columns as the number of fields that are written.
All available fields that can be accumulated are detailed in the following table, where it is also shown on which actuator model these fields are available:

.. table:: 
   :widths: 20, 20, 20, 40
   :align: center

   ============================== ====================== =================== ============================================================================
   **Output Field**               **Actuator Model**     **Units**           **Description**
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``time``                       All                    s                   simulation time.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``rtrAvgMagU``                 ADM, uniformADM, ALM   m/s                 rotor average wind speed magnitude.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``rtrAvgUpMagU``               ADM, uniformADM, ALM   m/s                 rotor average upstream wind speed magnitude.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``rtrThrust``                  ADM, uniformADM, AFM,  kN                  rotor thrust force.
                                  ALM
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``aeroPwr``                    ADM, uniformADM, AFM,  MW                  aerodynamic power.
                                  ALM
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``CtInf``                      ADM, uniformADM, ALM   n/a                 thrust coefficient based on inflow wind speed.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``CtLoc``                      ADM, uniformADM, ALM   n/a                 thrust coefficient based on local wind speed.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``CtUp``                       ADM, uniformADM, ALM   n/a                 thrust coefficient based on upstream wind speed.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``rtrTorque``                  ADM, ALM               kNm                 rotor torque.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``rtrOmega``                   ADM, ALM               rpm                 rotor rotational speed.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``genTorque``                  ADM, ALM               kNm                 generator torque.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``genPwr``                     ADM, ALM               MW                  generator power.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``genOmega``                   ADM, ALM               rpm                 generator rotational speed.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``collPitch``                  ADM, ALM               deg                 collective blade pitch angle.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``flowAngle``                  All                    deg                 flow angle. Only when yaw control is active.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``yawAngle``                   All                    deg                 yaw angle. Only when yaw control is active.
   ------------------------------ ---------------------- ------------------- ----------------------------------------------------------------------------
   ``azimuth``                    ALM                    deg                 azimuth angle. 
   ============================== ====================== =================== ============================================================================

For the ALM model, additional files are created, which contain the following fields as a function of the radial coordinate (columns) for each time sample (rows). 
These files are written inside the ``postProcessing/turbines/<startTime>`` directory, where ``<startTime>`` is the start time of the simulation. The name of 
the files is ``<turbineID>_<var>``, where *var* is any of the fields listed in the following table:

.. table::
    :widths: 20, 20, 20, 40
    :align: center

    ============================== =================== =================== ============================================================================
    **Output File**                **Actuator Model**  **Units**           **Description**
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``aoa``                        ALM                 degs                angle of attack at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``relVelMag``                  ALM                 m?s                 relative velocity magnitude at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``phi``                        ALM                 degs                relative velocity angle at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``uAxl``                       ALM                 m/s                 inflow velocity in axial direction at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``uRad``                       ALM                 m/s                 inflow velocity in radial direction at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``uTan``                       ALM                 m/s                 inflow velocity in tangential direction at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``Cl``                         ALM                 n/a                 lift coefficient at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``Cd``                         ALM                 n/a                 drag coefficient at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``axialF``                     ALM                 N                   axial force at radial points.
    ------------------------------ ------------------- ------------------- ----------------------------------------------------------------------------
    ``tangF``                      ALM                 N                   tangential force at radial points.
    ============================== =================== =================== ============================================================================

.. _ibm-force-subsection:

`-writePressureForce`
~~~~~~~~~~~~~~~~~~~~~

This functionality allows to write the pressure and viscous forces exerted by the flow onto the IBM body. 

.. _add-fields-subsection:

`additional fields`
~~~~~~~~~~~~~~~~~~~

For debugging purposes, TOSCA can output several additional 3D fields at every checkpoint file. These can be then visualized in serial when running ``tosca2PV``.
The following table summarizes the available fields and how to activate them in the ``control.dat`` file:

.. table:: 
   :widths: 20, 80
   :align: center

   ============================== ============================================================================
   **Output Field**               **Description**
   ------------------------------ ----------------------------------------------------------------------------
   ``Q``                          Q-criterion, activated with ``-computeQ`` set to 1 in the ``control.dat`` 
                                  file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Qg``                         second invariant of velocity gradient tensor, activated with ``-computeQg`` 
                                  set to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Rg``                         third invariant of velocity gradient tensor, activated with ``-computeRg`` 
                                  set to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Qs``                         second invariant of strain rate tensor, activated with ``-computeQs`` 
                                  set to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Rs``                         third invariant of strain rate tensor, activated with ``-computeRs`` 
                                  set to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Qr``                         second invariant of rotational rate tensor, activated with ``-computeQr`` 
                                  set to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``bf``                         wind farm body force field, activated with ``-computeFarmForce`` set to 1 
                                  in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Coriolis``                   Coriolis force field, activated with ``-computeSources`` set to 1 in the 
                                  ``control.dat`` file and ``coriolisActive`` set to 1 in the 
                                  ``ABLProperties.dat`` file. 
   ------------------------------ ----------------------------------------------------------------------------
   ``Driving``                    driving pressure gradient field, activated with ``-computeSources`` set to 1 
                                  in the ``control.dat`` file and ``controllerActive`` set to 1 in the 
                                  ``ABLProperties.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``Damping``                    body force in the fringe and damping layers, activated with 
                                  ``-computeSources`` set to 1 and one or more 
                                  between ``-xDampingLayer``, ``-zDampingLayer``, 
                                  ``-kLeftRayleigh`` and ``-kRightRayleigh`` set to 1 in the ``control.dat``
                                  file. 
   ------------------------------ ----------------------------------------------------------------------------
   ``CanopyForce``                canopy force field, activated with ``-computeSources`` and ``-canopy`` set 
                                  to 1 in the ``control.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``yDampU``                     mapped velocity from the x fringe region to the y fringe region to check the 
                                  result of the tiled mapping. Activated when ``-computeSources`` is set to 1
                                  in the ``control.dat`` file and both ``-xDampingLayer`` and 
                                  ``-yDampingLayer`` are set to 1 in the ``control.dat`` file. Notably, 
                                  ``uBarSelectionType`` should be set to 3 in the ``xDampingProperties`` 
                                  inside the ``ABLProperties.dat`` file.
   ------------------------------ ----------------------------------------------------------------------------
   ``divU``                       divergence of velocity, activated with ``-computeContinuity`` set to 1 
                                  in the ``control.dat`` file. 
   ------------------------------ ----------------------------------------------------------------------------
   ``buoyancy``                   buoyancy term of the momentum equation, activated with ``-computeBuoyancy`` 
                                  and ``-potentialT`` set to 1 in the ``control.dat`` file. 
   ============================== ============================================================================
   














