.. _applications_section: 

TOSCA Applications
====================

This section describes the main capabilities of the TOSCA package, highlighting steps to follow and considerations to make when 
setting up different types of simulations. The idea behind this section of the user guide is not that of explaining how to run 
a given case, but rather to show and explain how TOSCA can be leveraged to adopt best practices in the diverse engineering and 
scientific applications that can be studied with the code. 

The main focus of the TOSCA code has been, since the start of its development, atmospheric boundary layer flows and wind farm 
flows. The idea was to develop a tool that was able to tackle large-scale LES problems, such as atmospheric gravity waves (either
triggered by terrain features or wind farms), global blockage effects of wind farms, as well as the evolution of wakes shed by
large arrays of wind turbines. TOSCA is currently the only code that implements characteristic features of both spectral and 
finite volume codes. For example, the fringe region and concurrent precursor method are the only viable solution when one wants to  
accurately capture the wind farm interaction with the atmospheric gravity waves, as they ensure that these waves are not reflected at the 
domain boundaries. However, it does not exist a finite volume code that allows to run two simulations concurrently. Conversely, 
the parallel solution of a concurrent simulation is a well-known method, adopted in spectral codes, to nudge to an arbitrary inflow 
condition the flow reintroduced in the domain by periodic boundary conditions, which are implied in spectral methods. Luckily for 
those who possessed a spectral code, the concurrent precursor method resulted extremely efficient in damping gravity waves 
perturbations close to the boundary. As a result, prior to the introduction of TOSCA, the only LESs of wind farm induced gravity waves 
available in literature have been performed with spectral codes. 

On the other hand, spectral codes impose big limitation and poor flexibility. The mesh has to be cartesian and uniform 
in the majority of the methods, and it is very difficult to implement sharp-interface IBM methods due to the difficulty of 
reconstructing the boundary conditions. Conversely, this is easier in finite volume codes, which are more suited to simulate the 
flow over complex terrains both using IBM or terrain-following coordinates (if the ground is not too complex). 

Regarding the type of mesh, an unstructured mesh increases the flexibility of the code at the expense of its computational 
efficiency, as a connectivity matrix is required. This in turn requires more operations to perform cell access or to compute gradients. 
On the contrary, a structured mesh is much faster to solve and techniques such as grid nesting (overset) can be used to increase the 
resolution around specific regions of interest. 

TOSCA has been developed with all these considerations in mind. To make it a parallel-efficient code by construction, TOSCA is 
heavily based on PETSc, a state of the art library for the parallel solution of large numerical problems developed by the Argonne 
National Laboratory. The choiche was made to solve the Navier-Stokes equations in curvilinear coordinates, which offer some degree 
of flexibility (grid stretching and terrain-following meshes can be used) but still ensure simpler numerics with respect to an 
unstructured solver. Grid nesting has also been introduced after a major refactoring of the code in January 2022, which allows to 
introduce refined regions within the main background domain, and in the latest release the coupling has been made two-way. TOSCA features a dynamic 
sharp-interface immersed boundary method (IBM), meaning that any type of 3D object can be introduced and moved arbitrarily in the domain. The actual motion 
of the object is currently implemented through motion functions, but the introduction of additional types of motions, or even 
arbitrary motion time series, is straightforward as the IBM method is general by construction and does not assume any specific 
type of motion. TOSCA also features the state-of-the art turbine actuator models, such as
actuator line, actuator disk (uniform and non-uniform), actuator farm, as well as simpler canopy models. 
Depending on the selected model, turbines can be controlled in pitch, angular frequency and nacelle yaw. 

We are currently introducing in TOSCA a handful of data assimilation techniques, which allow to nudge the atmospheric boundary layer 
variables, making them follow observational data or data from weather prediction codes such as the Weather Research and Forecasting (WRF)
model. New turbulence models have been recently introduced, as well as integration with python in order 
to potentially conduct Machine Learning applications within the code. 

In the future, we plan on coupling TOSCA with the OpenFAST library, hence introducing blade flexibility, as well as introducing 
boundary coupling with WRF (this is already underway). There are ongoing efforts to port parts of the TOSCA code to GPUs and to also include multiphase 
capabilities to model ABL-wave interaction and floating wind turbines.  

Currently, TOSCA can be used to simulate atmopsheric boundary layers and wind farms over sea or complex terrain, including thermal 
stratification. The code can also be used to study urban environments and wind engineering problems (such as flow around ships, 
bridges and buildings). Furthermore, TOSCA can be (and it has been) used to simulate the flow around vertical axis wind turbines 
in a fully resolved manner, by leveraging the dynamic IBM method. 

Best practices and workflows on how to use the code are given for the following main applications. We highlight that this list is 
not complete and that the code could be used for additional problems that even the authors did not foresee. 

.. toctree::
   :maxdepth: 1

   03_applications/abl.rst
   03_applications/ibm.rst
   03_applications/precursor.rst
   03_applications/concurrent_precursor.rst
   03_applications/terrain.rst
   03_applications/assimilation.rst
   03_applications/turbulence.rst
