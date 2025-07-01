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