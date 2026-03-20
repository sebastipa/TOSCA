.. _theory-section: 

Theory Guide
============

TOSCA is a finite-volume code, formulated in generalized curvilinear coordinates, allowing it to take as input also non-Cartesian structured meshes. 
The present section is organized as follows. The governing equations in Cartesian coordinates are reported in Sec. :ref:`gov-equations-section`, TOSCA's 
numerical method — including the generalized curvilinear coordinate framework, the fractional-step projection method, velocity–temperature 
buoyancy coupling — are described in Sec. :ref:`numerics-section`. The available time-integration schemes for both the momentum and temperature equations are 
detailed in Sec. :ref:`time-integration-section`, while advection schemes for the divergence term are summarized in Sec. :ref:`advection-schemes-section`. 
The sub-grid scale turbulence models in curvilinear coordinates available in TOSCA are detailed in Sec. :ref:`sgs-model-section`. Velocity and temperature controllers
either used to avoid inertial oscillations produced by the Coriolis force in the free atmosphere when the initial condition is not consistent with the geostrophic wind,
to fix a given velocity at a reference height or to couple the TOSCA model with mesoscale model outputs (micro-scale to meso-scale coupling) are described in 
Sec. :ref:`controllers-section`. Sec. :ref:`precursor-section` details TOSCA's hybrid off-line/concurrent precursor methodology, which saves 
computational resources when performing the turbulence initialization in the precursor phase. The sharp-interface immersed boundary method (IBM), which 
enables the simulation of flows over complex terrain, objects and moving bodies, is described in Sec. :ref:`ibm-theory-section`. Actuator models used to 
represent wind turbines in the domain are described in Sec. :ref:`turbine-models-section`. Finally, an overview of TOSCA's parallel efficiency is given in 
Sec. :ref:`parallel_eff-section`, where the time per iteration is analyzed with increasing number of nodes and mesh elements on the Niagara 
high-performance computer at the SciNet HPC Consortium in Canada. Notably, TOSCA has been used to run finite wind farm simulations on the entire 
Niagara cluster (2024 nodes, 40 cores per node) and on all Cascade nodes of the University of British Columbia Sockeye cluster, at the  Advanced 
Research and Computing Lab, demonstrating its capability to handle massively-parallel computations.


.. toctree::
   :maxdepth: 1

   04_theory_guide/governing_equations.rst
   04_theory_guide/numerics.rst
   04_theory_guide/time_integration.rst
   04_theory_guide/advection_schemes.rst
   04_theory_guide/parallel_efficiency.rst
   04_theory_guide/turbulence_models.rst
   04_theory_guide/flow_controllers.rst
   04_theory_guide/precursor.rst
   04_theory_guide/ibm.rst
   04_theory_guide/turbine_models.rst

