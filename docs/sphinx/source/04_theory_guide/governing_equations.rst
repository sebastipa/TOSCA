.. _gov-equations-section:

Governing Equations
-------------------

TOSCA solves the incompressible Navier–Stokes equations for a flow with Coriolis forces augmented by the Boussinesq approximation for buoyancy. The buoyancy term is evaluated through the modified density :math:`\rho_k`, derived from a transport equation for the potential temperature. In Cartesian coordinates the equations read

.. math::
    :label: eq:massCartesian

    \frac{\partial u_i}{\partial x_i} = 0

.. math::
    :label: eq:momentumCartesian

    \frac{\partial u_i}{\partial t} + \frac{\partial}{\partial x_j}\left(u_j u_i\right)
    &= -\frac{1}{\rho_0}\frac{\partial p}{\partial x_i}
    + \frac{\partial}{\partial x_j}\left[\nu_\text{eff}\left(
        \frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}
    \right)\right] \\
    & - \frac{1}{\rho_0}\frac{\partial p_{\infty}}{\partial x_i}
    + \frac{\rho_k}{\rho_0}g_i
    - 2\epsilon_{ijk}\Omega_j u_k
    + f_i + s^v_i + s^h_i

.. math::
    :label: eq:temperatureCartesian

    \frac{\partial\theta}{\partial t} + \frac{\partial}{\partial x_j}\left(u_j\theta\right)
    = \frac{\partial}{\partial x_j}\left(\kappa_\text{eff}\frac{\partial \theta}{\partial x_j}\right)

where :math:`u_i` is the Cartesian velocity, :math:`p/\rho_0` is the kinematic pressure, and :math:`\theta` is the potential temperature defined as :math:`\theta = T(p_0/p)^{R/c_p}` (with :math:`T` the absolute temperature, :math:`R` the gas specific constant, :math:`c_p` the specific heat at constant pressure, and :math:`p_0` a reference pressure). The vector :math:`g_i` is the gravitational acceleration, and :math:`\Omega_j` is the planetary rotation rate vector, defined as :math:`\omega\cos\phi\,\hat{y}+\omega\sin\phi\,\hat{z}` at latitude :math:`\phi` in a local frame with :math:`\hat{z}` pointing against gravity, :math:`\hat{x}` tangent to Earth's parallels, and :math:`\hat{y}` completing the right-handed triad. Source terms :math:`f_i`, :math:`s^v_i`, and :math:`s^h_i` are body forces due to wind turbines, and vertical and horizontal damping regions, respectively.

The modified density :math:`\rho_k` is

.. math::
    :label: eq:modifiedDensity

    \frac{\rho_k}{\rho_0} = 1 - \left(\frac{\theta-\theta_0}{\theta_0}\right)

where :math:`\theta_0` is a reference potential temperature chosen as the ground temperature. The effective viscosity :math:`\nu_\text{eff} = \nu + \nu_t` and effective thermal diffusivity :math:`\kappa_\text{eff} = \kappa + \kappa_t` include both molecular and sub-grid scale (SGS) contributions; the SGS quantities :math:`\nu_t` and :math:`\kappa_t` are described in :ref:`sgs-model-section`. The Prandtl number :math:`Pr = \nu/\kappa` is set to 0.7 by default. The term :math:`-\rho_0^{-1}\partial p_\infty/\partial x_i` is a prescribed horizontal pressure gradient that drives the flow and represents the velocity controller described in :ref:`controllers-section`.

The discretized form of the governing equations solved by TOSCA — including the curvilinear coordinate transformation, contravariant flux variables, staggered grid arrangement, and fractional-step projection — is described in :ref:`numerics-section`.
