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
within the computational domain and data is tri-linearly interpolated at the surface points. The classification of each different section 
type available in TOSCA is given in the following table: 

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
   *userSections*                 user-defined        Automatically saves velocity, pressure, effective viscosity and potential 
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

.. _mke-budgets-subsection:

`-keBudgets`
~~~~~~~~~~~~

.. _averaging-abl-subsection:

`-averageABL`
~~~~~~~~~~~~~

.. _averaging-3lm-subsection:

`-average3LM`
~~~~~~~~~~~~~

.. _perturb-abl-subsection:

`-perturbABL`
~~~~~~~~~~~~~

.. _turbine-data-subsection:

`turbine data`
~~~~~~~~~~~~~~

.. _ibm-force-subsection:

`-writePressureForce`
~~~~~~~~~~~~~~~~~~~~~

.. _add-fields-subsection:

`additional fields`
~~~~~~~~~~~~~~~~~~~

















