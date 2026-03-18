.. _numerics-section:

Numerical Method
----------------

Generalized Curvilinear Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TOSCA solves the governing equations in **generalized curvilinear coordinates**
:math:`(l_1, l_2, l_3) \equiv (\xi, \eta, \zeta)`, related to the Cartesian
coordinates :math:`(x_1, x_2, x_3)` through a differentiable, invertible mapping
:math:`x_i = x_i(\xi,\eta,\zeta)`.  The **Jacobian** of the transformation is

.. math::

    J = \det\!\left(\frac{\partial x_i}{\partial l_q}\right).

The **face area vectors** are defined as

.. math::

    S^q_i = \frac{1}{J}\,\frac{\partial l_q}{\partial x_i},
    \qquad q,i \in \{1,2,3\}.

Geometrically, :math:`\mathbf{S}^q \equiv (S^q_1, S^q_2, S^q_3)` is the outward
unit normal of the :math:`q`-th cell face scaled by the ratio of face area to cell
volume.  The **contravariant flux** through the :math:`q`-th face is

.. math::

    V^q = J S^q_i\, u_i \equiv J\,\mathbf{S}^q \cdot \mathbf{u}.

The **contravariant metric tensor** (absorbing :math:`J` for numerical convenience)
is the symmetric positive-definite tensor

.. math::

    S^{qr} = J S^q_i S^r_i.

Differential Operators in Curvilinear Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The four fundamental operators used throughout the discretization are:

**Gradient** of a scalar :math:`\phi`:

.. math::

    \frac{\partial \phi}{\partial x_i} = S^q_i\,\frac{\partial \phi}{\partial l_q}
    \qquad (\text{sum over }q).

**Divergence** of a vector field :math:`\mathbf{a}`:

.. math::

    \frac{\partial a_i}{\partial x_i}
    = \frac{1}{J}\,\frac{\partial}{\partial l_q}\!\left(J S^q_i a_i\right)
    = \frac{1}{J}\,\frac{\partial V^q_a}{\partial l_q},

where :math:`V^q_a = J S^q_i a_i` is the contravariant flux of :math:`\mathbf{a}`.
For an incompressible velocity field this reduces to :math:`\partial V^q/\partial l_q = 0`.

**Divergence** of a rank-2 tensor :math:`T_{ij}`:

.. math::

    \frac{\partial T_{ij}}{\partial x_j}
    = \frac{1}{J}\,\frac{\partial \tau^r_i}{\partial l_r},
    \qquad \tau^r_i = J S^r_j\, T_{ij}.

The quantity :math:`\tau^r_i` is the :math:`r`-th contravariant flux of the
:math:`i`-th row of the tensor.  For the viscous stress
:math:`T_{ij} = \nu_\text{eff}(\partial u_i/\partial x_j + \partial u_j/\partial x_i)`,
each component is evaluated using the gradient formula above, yielding the viscous
flux terms of the momentum equation.

**Laplacian** of a scalar :math:`\phi`:

.. math::

    \nabla^2\phi
    = \frac{1}{J}\,\frac{\partial}{\partial l_q}
      \!\left(S^{qr}\,\frac{\partial \phi}{\partial l_r}\right).

The diagonal entries of :math:`S^{qr}` drive face-normal diffusion, while the
off-diagonal entries couple adjacent curvilinear directions on non-orthogonal
grids.

Governing Equations in Curvilinear Form
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Applying the operators above, the governing equations solved by TOSCA take the
**partially-transformed** form (Ge & Sotiropoulos, 2007).  The momentum
equation is written for the contravariant flux :math:`V^q`, which eliminates
Christoffel-symbol correction terms and maintains second-order accuracy on any
smoothly varying mesh.  The system reads:

.. math::
    :nowrap:

    \begin{align}
        \frac{\partial V^q}{\partial l_q} &= 0, \\[6pt]
        \frac{\partial V^q}{\partial t}
        + J S^q_i \frac{\partial}{\partial l_r}\!\left(V^r u_i\right)
        &= -\frac{J^2}{\rho_\text{ref}}\frac{\partial p}{\partial l_r}S^{rq}
        + J S^q_i \frac{\partial}{\partial l_r}
          \!\left[J\nu_\text{eff}\!\left(S^{rn}\frac{\partial u_i}{\partial l_n}
          + S^k_j R_{ij}\right)\right] \\
        &\quad
        - \frac{J^2}{\rho_\text{ref}}\frac{\partial p_\infty}{\partial l_r}S^{rq}
        + S^q_i\frac{\rho_k}{\rho_\text{ref}}g_i
        - 2 S^q_i\epsilon_{ijk}\Omega_j u_k \\
        &\quad
        + S^q_i\!\left(f_i + s^v_i + s^h_i\right), \\[6pt]
        \frac{\partial\theta}{\partial t}
        + \frac{\partial}{\partial l_r}\!\left(V^r\theta\right)
        &= J\frac{\partial}{\partial l_r}
          \!\left(J\kappa_\text{eff}S^{rk}\frac{\partial\theta}{\partial l_k}\right)
        + s_\theta,
    \end{align}

where :math:`R_{ij} = S^k_j\,\partial u_i/\partial l_k` is the transpose-gradient
term that ensures symmetry of the viscous flux, :math:`p_\infty` is the large-scale
background driving pressure, :math:`\rho_k/\rho_\text{ref}` is the Boussinesq
density ratio, and :math:`f_i`, :math:`s^v_i`, :math:`s^h_i`, :math:`s_\theta`
are body forces from turbines, vertical damping, horizontal damping, and the
temperature controller, respectively (see :ref:`gov-equations-section`).

Staggered Grid Arrangement
^^^^^^^^^^^^^^^^^^^^^^^^^^^

TOSCA employs a **hybrid collocated–staggered** discretization
(Ge & Sotiropoulos, 2007).  Contravariant fluxes :math:`V^q` are stored at the
center of their respective cell face (staggered), ensuring exact discrete
mass conservation.  Pressure :math:`p`, potential temperature :math:`\theta`, and
the sub-grid viscosity :math:`\nu_t` are stored at cell centers (collocated).

Cartesian velocity :math:`\mathbf{u}` is required at cell faces for the advection
and viscous operators.  Rather than solving three Cartesian momentum equations at
each face (which would triple the cost), TOSCA first assembles the momentum
right-hand side at cell centers and then interpolates it to the faces using
second-order linear interpolation along each grid direction.  After each step,
Cartesian velocity is reconstructed at cell centers from the updated contravariant
fluxes by inverting :math:`V^q = J S^q_i u_i`.

Fractional-Step Projection Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pressure–velocity coupling uses a second-order **fractional-step projection**
method (van Kan, 1986).  Given the state
:math:`(\mathbf{V}^n, p^n, \theta^n)` at time level :math:`n`, one time step
proceeds through the four stages below.

**Stage 1 — Momentum predict.**
Advance the contravariant fluxes to a predicted state
:math:`\mathbf{V}^{n+1}_*` that satisfies the momentum equation but not yet the
divergence-free constraint.  Schematically (Crank–Nicolson form):

.. math::

    \frac{V^{q,n+1}_* - V^{q,n}}{\Delta t}
    = \frac{1}{2}\!\left[H^q_U(\mathbf{V}^{n+1}_*,\mathbf{u}^{n+1}_*)
                       + H^q_U(\mathbf{V}^n,\mathbf{u}^n)\right]
    + L^q(\mathbf{u}^n) + G^q(p^n).

Here :math:`H^q_U` contains the advection and viscous fluxes (see
:ref:`advection-schemes-section`); :math:`L^q` contains buoyancy, Coriolis,
turbine, and damping source terms; and :math:`G^q(p)` is the discrete
pressure-gradient operator evaluated at the :math:`q`-th faces:

.. math::

    G^q(p)\big|_\text{face} = J\!\left(
    S^{q1}\frac{\partial p}{\partial l_1} +
    S^{q2}\frac{\partial p}{\partial l_2} +
    S^{q3}\frac{\partial p}{\partial l_3}
    \right)_\text{face}.

A full description of the available time-stepping choices is given in
:ref:`time-integration-section`.

**Stage 2 — Pressure correction (Poisson solve).**
Compute the scalar correction :math:`\phi^{n+1}` from:

.. math::

    D\!\left(G^q(\phi^{n+1})\right) = -\frac{1}{\Delta t}\,D\!\left(V^{q,n+1}_*\right),

where the **discrete divergence operator** at cell :math:`(i,j,k)` is

.. math::

    D_{ijk}(\gamma) =
      \gamma_{i+\frac{1}{2},j,k}
    + \gamma_{i,j+\frac{1}{2},k}
    + \gamma_{i,j,k+\frac{1}{2}}
    - \gamma_{i-\frac{1}{2},j,k}
    - \gamma_{i,j-\frac{1}{2},k}
    - \gamma_{i,j,k-\frac{1}{2}}.

Substituting the pressure-gradient operator into the divergence gives a sparse
19-point linear system per cell (1 diagonal, 6 face-centred, 12 edge-centred
cross-derivative contributions from :math:`S^{qr}`).  The system is solved with
HYPRE's BoomerAMG algebraic multigrid preconditioner (GMRES or PCG outer Krylov
method), providing :math:`\mathcal{O}(N)` scalable performance.
The constant mode of :math:`\phi` is removed after the solve.

**Stage 3 — Flux projection.**
Correct the predicted fluxes to enforce mass conservation:

.. math::

    V^{q,n+1} = V^{q,n+1}_* + \Delta t\, G^q(\phi^{n+1}).

The pressure is updated as :math:`p^{n+1} = p^n + \phi^{n+1}`.

**Stage 4 — Temperature solve.**
With the divergence-free :math:`\mathbf{V}^{n+1}` available, advance
:math:`\theta^n \to \theta^{n+1}` using the selected temperature
time-integration scheme (see :ref:`time-integration-section`).

Velocity–Temperature Buoyancy Coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The buoyancy source :math:`L^{q,\text{buoy}}` in the momentum equation depends on
the face-interpolated potential temperature through the Boussinesq density ratio.
Because :math:`\theta` is solved only in Stage 4, it cannot be used implicitly in
Stage 1.

Before the momentum solve, the buoyancy body-force vector
:math:`\mathbf{b}(\theta^n)` is computed from the current temperature field
:math:`\theta^n` and stored.  The previous step's buoyancy
:math:`\mathbf{b}(\theta^{n-1})` is also retained.  The strategy depends on
the momentum time-integration scheme selected:

- **Forward Euler and Runge–Kutta 4** — buoyancy is an explicit source at the
  current time level:

  .. math::

      L^{q,\text{buoy}} = \mathbf{b}(\theta^n).

  This introduces a first-order :math:`\mathcal{O}(\Delta t)` lag between the
  temperature field and its effect on the momentum.

- **Crank–Nicolson/SNES and IMEX-CNAB** — buoyancy is extrapolated to the
  mid-step level :math:`n+\tfrac{1}{2}` using Adams–Bashforth 2:

  .. math::

      L^{q,\text{buoy}} = \tfrac{3}{2}\,\mathbf{b}(\theta^n)
                        - \tfrac{1}{2}\,\mathbf{b}(\theta^{n-1}).

  For CN/SNES this is consistent with the Crank–Nicolson time-averaging of the
  advection and diffusion operators; for IMEX-CNAB it matches the AB2 treatment
  of the convective operator.  Both reduce the temporal coupling error to
  :math:`\mathcal{O}(\Delta t^2)`.  On the first time step both formulae revert
  to Forward Euler.

No temperature predictor step is performed; the only temperature field available
at the start of Stage 1 is :math:`\theta^n`.
