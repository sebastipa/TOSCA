Theory Guide
============

TOSCA is a finite-volume code, formulated in generalized curvilinear coordinates, allowing it to take as input also non-Cartesian structured meshes. The present section is organized as follows. The governing equations in Cartesian coordinates are reported in Sec. :ref:`gov-equations-section`, while actuator models used to represent wind turbines in the domain are described in Sec. :ref:`turbine-models-section`.

TOSCA's numerical method, the governing equations in curvilinear coordinates that we actually solve, and a brief overview of generalized curvilinear coordinates are reported in Sec. :ref:`numerics-section`, while the LES turbulence model in the curvilinear frame is detailed in Sec. :ref:`sgs-model-section`. 

An overview of TOSCA's parallel efficiency is given in Sec. :ref:`parallel_eff-section`, where we analyze the time per iteration with increasing number of nodes and mesh elements on the Niagara high-performance computer at the SciNet HPC Consortium. In addition, TOSCA has been used to run finite wind farm simulations on the whole Niagara cluster (2024 nodes, 40 cores per node) and on all Cascade nodes of the UBC-ARC Sockeye cluster, demonstrating its capability to handle massively-parallel computations.

In order to run ABL simulations, we developed a novel methodology, described in Sec. :ref:`controllers-section`, that enforces a desired hub-height wind speed while simultaneously avoiding inertial oscillations produced by the Coriolis force above the boundary layer. In addition, the disagreement exists between different CFD codes in predicting the final mean potential temperature profile inside the boundary layer, we developed a mean temperature controller which maintains a prescribed average potential temperature profile, harmonizing the comparison of simulation results produced with different codes in future studies. Finally, Sec. :ref:`precursor-section` details TOSCA's hybrid off-line/concurrent precursor methodology, which saves computational resources when performing the turbulence initialization in the precursor phase. 

.. toctree::
   :maxdepth: 1

   04_theory_guide/governing_equations.rst
   04_theory_guide/numerics.rst
   04_theory_guide/parallel_efficiency.rst
   04_theory_guide/turbulence_models.rst
   04_theory_guide/flow_controllers.rst
   04_theory_guide/precursor.rst
   04_theory_guide/turbine_models.rst

