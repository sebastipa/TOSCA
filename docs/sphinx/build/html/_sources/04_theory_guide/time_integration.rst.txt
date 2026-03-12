.. _time-integration-section:

Time Integration Schemes
------------------------

Four time-stepping schemes are available for both the momentum and temperature equations.  The momentum scheme is selected via the keyword ``-dUdtScheme`` and
the temperature scheme via ``-dTdtScheme`` in ``control.dat``.

Momentum Equation
~~~~~~~~~~~~~~~~~~

Let :math:`\mathcal{H}^q_U = H^q_U + L^q` denote the full right-hand side of the momentum equation (advection, viscous diffusion, and all source terms, excluding
the pressure gradient).  The pressure gradient :math:`G^q(p^n)` is always treated explicitly and removed by the projection step (see :ref:`numerics-section`).

**Forward Euler**

First-order explicit scheme:

.. math::

    V^{q,n+1}_* = V^{q,n} + \Delta t\,\left[\mathcal{H}^q_U(V^n, u^n) - G^q(p^n)\right].

Conditionally stable; requires :math:`\mathrm{CFL} < 1`.  Recommended only for coarse (fast) exploratory runs.

**Runge–Kutta 4**

Classical four-stage explicit scheme:

.. math::

    V^{q,n+1}_* = V^{q,n} + \Delta t\sum_{s=1}^{4} b_s K_s,
    \qquad
    K_s = \mathcal{H}^q_U\bigl(V^n + a_s\,\Delta t\,K_{s-1},\,u^n + a_s\,\Delta t\,K_{s-1}\bigr),

with Butcher coefficients :math:`\mathbf{b} = (\tfrac{1}{6},\tfrac{1}{3},\tfrac{1}{3},\tfrac{1}{6})^\top`
and :math:`\mathbf{a} = (0,\tfrac{1}{2},\tfrac{1}{2},1)^\top`. Fourth-order accurate in time; conditionally stable. Requires four full
right-hand-side evaluations per step.

**Crank–Nicolson**

Fully implicit second-order scheme solved with PETSc's Newton trust-region SNES solver.  The update reads:

.. math::

    \frac{V^{q,n+1}_* - V^{q,n}}{\Delta t} = \frac{1}{2}\,\mathcal{H}^q_U(V^{n+1}_*, u^{n+1}_*)
    + \frac{1}{2}\,\mathcal{H}^q_U(V^n, u^n) - G^q(p^n).

The nonlinear system is solved using Newton–Krylov iterations with a matrix-free Jacobian (MATMFFD) and a BiCGStab or GMRES (more robust slower) inner solver.  
To reduce the number of Newton iterations, a source buffer :math:`\mathbf{b}_U` precomputes all terms that are independent of the current iterate (pressure
gradient, :math:`\tfrac{1}{2}\mathcal{H}^q_U(V^n,u^n)`, source terms), and an explicit Forward Euler prediction provides the initial guess.  The scheme is
unconditionally stable and second-order accurate in time. Recommended for production runs where the cost of nonlinear solves is offset by the ability to use 
large time steps (CFL :math:`\to 1`).

**IMEX-CNAB (Implicit–Explicit Crank–Nicolson Adams–Bashforth)**

A split scheme where convection is treated **explicitly** with Adams–Bashforth 2 (AB2) and viscous diffusion is treated **implicitly** with Crank–Nicolson.
The update reads:

.. math::

    \frac{V^{q,n+1}_* - V^{q,n}}{\Delta t} = \mathbf{b}^q_U + \tfrac{1}{2}\,\mathcal{L}^q(V^{n+1}_*),

where the **explicit buffer** :math:`\mathbf{b}^q_U` is assembled once per step from fields at time levels :math:`n` and :math:`n-1`:

.. math::

    \mathbf{b}^q_U
    = -G^q(p^n)
    + L^{q,\text{buoy}}
    + \tfrac{3}{2}\,\mathcal{N}^q_U(V^n, u^n)
    - \tfrac{1}{2}\,\mathcal{N}^q_U(V^{n-1}, u^{n-1})
    + \tfrac{1}{2}\,\mathcal{L}^q(V^n)
    + \mathcal{S}^q,

with :math:`\mathcal{N}^q_U` the convective operator, :math:`\mathcal{L}^q` the viscous-diffusion operator, and :math:`\mathcal{S}^q` the remaining source terms.
The first step uses Forward Euler for the AB2 convection stencil.  The resulting **linear** system

.. math::

    A\,V^{q,n+1}_* = V^{q,n} + \Delta t\,\mathbf{b}^q_U,
    \qquad A\,v \equiv v - \tfrac{\Delta t}{2}\,\mathcal{L}^q(v),

is solved directly with a standalone PETSc KSP (GMRES or BiCGStab, with :math:`L=2`) solver. The scheme is second-order accurate in time and overall faster 
than Crank-Nicholson in terms of time advancement of the simulation. However, because convection is treated explicitly, the scheme is only conditionally stable and 
often requires stabilization through the ``-central4`` + ``-hyperVisc`` combination (see :ref:`advection-schemes-section`) or the addition of biharmonic hyperviscosity 
(see :ref:`biharmonic-hyperviscosity` below) to suppress Nyquist-mode noise at high CFL numbers. The recommended CFL limit is :math:`\mathrm{CFL} \lesssim 0.5`, if no explicit
filtering is used, or :math:`\mathrm{CFL} \lesssim 0.6` when biharmonic hyperviscosity is active (either as a standalone correction or as a mild filter on top of the ``-central4`` scheme).
Running this scheme with CFL close to or above 1.0 will lead to numerical instability.

.. _biharmonic-hyperviscosity:

IMEX-CNAB Biharmonic Hyperviscosity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the IMEX-CNAB scheme, an optional biharmonic (fourth-order) hyperviscosity correction is available to suppress odd–even checkerboard noise and instability
at the Nyquist wavenumber.  The per-step flux correction is

.. math::

    \Delta V^q_{\mathrm{hv}} = -\varepsilon_4
    \bigl(\delta^4_\xi + \delta^4_\eta + \delta^4_\zeta\bigr)\,V^{q,n},

where the one-dimensional fourth-order difference operator is

.. math::

    \delta^4_d\,v_n = v_{n+2} - 4v_{n+1} + 6v_n - 4v_{n-1} + v_{n-2},

with Fourier symbol :math:`\widehat{\delta^4}(\kappa) = 4(\cos\kappa h - 1)^2`. The symbol vanishes at :math:`\kappa h = 0` (no effect on resolved scales) and
reaches 16 at the Nyquist wavenumber :math:`\kappa h = \pi`.  For stable, non-amplifying operation the coefficient must satisfy

.. math::

    \varepsilon_4 \leq \frac{1}{48} \approx 0.021.

The coefficient is directional and it is set by the keywords ``-imexHyperVisc4U_i``, ``-imexHyperVisc4U_j`` and ``-imexHyperVisc4U_k`` in ``control.dat``. 
When hyperviscosity is active, the effective CFL limit for IMEX-CNAB increases from 0.5 to approximately 0.6-0.7. Notably, hyperviscosity addition can be 
directional, allowing for selective damping of the Nyquist mode in the horizontal or vertical directios. 

Temperature Equation
~~~~~~~~~~~~~~~~~~~~~

Let :math:`\mathcal{F}(T) = \mathcal{N}(T) + \mathcal{D}(T)` be the full right-hand side of the temperature equation (convection :math:`\mathcal{N}`
plus diffusion :math:`\mathcal{D}`), including Rayleigh-damping and source terms.

**Runge–Kutta 4**

Four-stage explicit scheme, identical in structure to the momentum counterpart:

.. math::

    T^{n+1} = T^n + \Delta t\sum_{s=1}^{4} b_s K_s,
    \qquad
    K_s = \mathcal{F}\!\bigl(T^n + a_s\,\Delta t\,K_{s-1}\bigr),

with the same Butcher coefficients as above.  Conditionally stable; requires :math:`\mathrm{CFL} < 1`.

**Backward Euler**

Fully implicit first-order scheme.  The nonlinear residual

.. math::

    \mathcal{G}(T^{n+1}) = -\frac{T^{n+1} - T^n}{\Delta t} + \mathcal{F}(T^{n+1}) = 0

is solved using PETSc's SNES with the same Newton–Krylov configuration used for the momentum equation.  Unconditionally stable; first-order accurate in time.

**BDF2 (Second-order Backward Differentiation Formula)**

Second-order fully implicit scheme that uses a three-level time-derivative approximation.  The nonlinear residual is

.. math::

    \mathcal{G}(T^{n+1}) = -\frac{\tfrac{3}{2}T^{n+1} - 2T^n + \tfrac{1}{2}T^{n-1}}{\Delta t}
    + \mathcal{F}(T^{n+1}) = 0,

which is solved with the same Newton–Krylov SNES configuration as Backward Euler.  On the first time step BDF2 falls back to Backward Euler (only :math:`T^n` is
available).  BDF2 is A-stable and mildly dissipative at high wavenumbers (:math:`|G| < 1`), unlike Crank–Nicolson which has :math:`|G| = 1` exactly and
can accumulate high-frequency errors in stratified LES.  An extra storage vector :math:`T^{n-1}` (``Tmprt_oo``) is allocated only when this scheme is selected.
Recommended for stable stratification runs where second-order accuracy in time is required.

**IMEX-BEAB (Implicit–Explicit Backward Euler Adams–Bashforth 2)**

Implicit–explicit scheme: convection treated explicitly with Adams–Bashforth 2 and diffusion treated implicitly with Backward Euler.  The scheme avoids the
need for :math:`\mathbf{u}^{n+1}` and yields a **linear** system in :math:`T^{n+1}`:

.. math::

    \underbrace{\bigl(I - \Delta t\,\mathcal{D}\bigr)}_{A}\,T^{n+1}
    = T^n + \Delta t\,b_T,

where the explicit buffer is

.. math::

    b_T = c_1\,\mathcal{N}(T^n) + c_2\,\mathcal{N}(T^{n-1})
        + \mathcal{D}_\text{damp}(T^n) + \frac{s_T}{\Delta t},

with :math:`c_1 = 3/2`, :math:`c_2 = -1/2` (AB2; on the first step: :math:`c_1 = 1`, :math:`c_2 = 0`), :math:`\mathcal{D}` the diffusion operator,
:math:`\mathcal{D}_\text{damp}` the Rayleigh damping operator, and :math:`s_T` a prescribed source term.  The operator :math:`A` is symmetric positive definite;
the system is solved with the same standalone KSP (GMRES or BiCGStab) used for the IMEX momentum equation.  

Practical Combination of Schemes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following table, the recommended combinations of momentum and temperature time-stepping schemes are listed for different types of workflow.

.. list-table::
   :header-rows: 1
   :widths: 17 16 19 48

   * - U Scheme
     - T Scheme
     - Linear Solver
     - Notes
   * - Forward Euler
     - Runge–Kutta 4
     - N/A
     - :math:`\mathrm{CFL} < 1`. Fast exploratory runs only, first-order accurate in time. Not recommended for production.
   * - Runge–Kutta 4
     - Runge–Kutta 4
     - N/A
     - :math:`\mathrm{CFL} < 1`. Conditionally stable, fourth-order accurate in time. Suitable for non-turbulent or idealized runs.
   * - Crank–Nicolson
     - BDF2 / Backward Euler
     - BiCGStab inside SNES, GMRES as fallback.
     - *robust choice*, unconditionally stable (:math:`\mathrm{CFL} \gtrsim 0.9`). Production runs where the nonlinear-solve cost is offset by large time steps. BDF2 gives second-order accuracy in time with mild A-stable dissipation; Backward Euler is first-order but simpler.
   * - IMEX-CNAB
     - BDF2 / Backward Euler
     - *fastest choice*, standalone GMRES KSP; BiCGStab as a fallback.
     - :math:`\mathrm{CFL} \leq 0.5` without any hyperviscosity, :math:`\mathrm{CFL} \leq 0.6`–:math:`0.7` with explicit filtering (recommended horizontal coefficient around 0.003; set via ``-imexHyperVisc4U_i/j/k``) or added advection scheme viscosity (only for ``-central4``). Production runs balancing stability and efficiency; second-order accurate in time.
