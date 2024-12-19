Theory Guide
============

TOSCA is a finite-volume code, formulated in generalized curvilinear coordinates, allowing it to take as input also non-Cartesian structured meshes. The present section is organized as follows. The governing equations in Cartesian coordinates are reported in Sec. :ref:`gov-equations-section`, while actuator models used to represent wind turbines in the domain are described in Sec. :ref:`turbine-models-section`.

TOSCA's numerical method, the governing equations in curvilinear coordinates that we actually solve, and a brief overview of generalized curvilinear coordinates are reported in Sec. :ref:`numerics-section`, while the LES turbulence model in the curvilinear frame is detailed in Sec. :ref:`sgs-model-section`. 

An overview of TOSCA's parallel efficiency is given in Sec. :ref:`parallel_eff-section`, where we analyze the time per iteration with increasing number of nodes and mesh elements on the Niagara high-performance computer at the SciNet HPC Consortium. In addition, TOSCA has been used to run finite wind farm simulations on the whole Niagara cluster (2024 nodes, 40 cores per node) and on all Cascade nodes of the UBC-ARC Sockeye cluster, demonstrating its capability to handle massively-parallel computations.

In order to run ABL simulations, we developed a novel methodology, described in Sec. :ref:`controllers-section`, that enforces a desired hub-height wind speed while simultaneously avoiding inertial oscillations produced by the Coriolis force above the boundary layer. In addition, the disagreement exists between different CFD codes in predicting the final mean potential temperature profile inside the boundary layer, we developed a mean temperature controller which maintains a prescribed average potential temperature profile, harmonizing the comparison of simulation results produced with different codes in future studies. Finally, Sec. :ref:`precursor-section` details TOSCA's hybrid off-line/concurrent precursor methodology, which saves computational resources when performing the turbulence initialization in the precursor phase. 

.. _gov-equations-section:

Governing Equations
-------------------

Governing equations correspond to mass and momentum conservation for an incompressible flow with Coriolis forces and Boussinesq approximation for the buoyancy term. The latter is calculated using the modified density :math:`\rho_k`, evaluated by solving a transport equation for the potential temperature. These equations, expressed in Cartesian coordinates using tensor notation read
   
.. math::
    :label: eq:massCartesian

    \frac{\partial u_i}{\partial x_i} = 0

.. math::
    :label: eq:momentumCartesian

    \frac{\partial u_i}{\partial t} + \frac{\partial}{\partial x_j}\left(u_j u_i\right) &=
    -\frac{1}{\rho_0}\frac{\partial p}{\partial x_i} 
    + \frac{\partial}{\partial x_j}\left[\nu_\text{eff}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)\right] \\
    & - \frac{1}{\rho_0}\frac{\partial p_{\infty}}{\partial x_i} 
    + \frac{\rho_k}{\rho_0}g_i - 2\epsilon_{ijk}\Omega_j u_k 
    + f_i + s^v_i + s^h_i

.. math::
    :label: eq:temperatureCartesian

    \frac{\partial\theta}{\partial t} + \frac{\partial}{\partial x_j}\left(u_j\theta\right) =
    \frac{\partial}{\partial x_j}\left(k_\text{eff}\frac{\partial \theta}{\partial x_j}\right)

where :math:`u_i` is the Cartesian velocity, :math:`p/\rho_0` is the kinematic pressure, :math:`\theta` is the potential temperature, defined as :math:`\theta = T(p_0/p)^{R/c_p}` (`T` is the absolute temperature, :math:`R` is the gas specific constant, :math:`c_p` is the specific heat at constant pressure and :math:`p_0` is the reference pressure), :math:`g_i` is the gravitational acceleration vector, :math:`\Omega_j` is the rotation rate vector at an arbitrary location on the planetary surface (defined as :math:`\omega\cos\phi\widehat{y}+\omega\sin\phi\widehat{z}`, where :math:`\phi` is the latitude, in a local reference frame having :math:`\widehat{z}` aligned and opposite to the gravitational acceleration vector, :math:`\widehat{x}` tangent to Earth's parallels and :math:`\widehat{y}` such that the frame is right-handed). Source terms :math:`f_i`, :math:`s^v_i`, and :math:`s^h_i` are body forces introduced by turbines, and by vertical and horizontal damping regions, respectively. Moreover, the modified density :math:`\rho_k` is defined as

.. math::
    :label: eq:modifiedDensity
    
	\frac{\rho_k}{\rho_0} = 1 - \left(\frac{\theta-\theta_0}{\theta_0}\right) 


where :math:`\theta_0` is a reference potential temperature, chosen as the ground temperature. Parameters :math:`\nu_\text{eff}` and :math:`\kappa_\text{eff}` are the effective viscosity and thermal diffusivity respectively. The former is the sum of the kinematic viscosity :math:`\nu` and the sub-grid scale viscosity :math:`\nu_t`, while the latter is sum between the thermal diffusivity :math:`\kappa = \nu/Pr` and the turbulent thermal diffusivity :math:`\kappa_t`. The computation of both :math:`\nu_t` and :math:`\kappa_t` is detailed in Sec. :ref:`sgs-model-section`, while the Prandtl number `Pr` is set to 0.7 by default, although it can be changed in the *control.dat* file. The third term on the right-hand side of Eq. :eq:`eq:momentumCartesian` is a uniform horizontal pressure gradient that balances turbulent stresses and the Coriolis force, allowing the boundary layer to reach a statistically steady state. This term is commonly referred to as velocity controller, and it is explained in Sec. :ref:`controllers-section`.

.. _turbine-models-section:

Turbine Models
--------------

.. _sgs-model-section:

Sub-grid Scale Model
--------------------

.. _controllers-section:

Controllers
-----------

.. _precursor-section:

Precursor, Fringe Regions & Damping Layers 
------------------------------------------

.. _numerics-section:

Numerical Method & Curvilinear Coordinates
------------------------------------------

.. _parallel_eff-section:

Parallel Efficiency 
-------------------
