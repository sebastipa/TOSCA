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

    V^q = S^q_i\, u_i \equiv \,\mathbf{S}^q \cdot \mathbf{u}.

The **contravariant metric tensor**, multiplied by a factor :math:`1/J^2` for numerical convenience, is the symmetric positive-definite tensor

.. math::

    S^{qr} = S^q_i S^r_i.

Differential Operators in Curvilinear Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The four fundamental operators used throughout the discretization are:

**Gradient** of a scalar :math:`\phi`:

.. math::

    \frac{\partial \phi}{\partial x_i} = J S^k_i\,\frac{\partial \phi}{\partial l_k}

**Gradient** of a vector field :math:`\mathbf{a}`:

.. math::

    \frac{\partial a_j}{\partial x_i} = J S^k_i\,\frac{\partial a_j}{\partial l_k}

**Divergence** of a vector field :math:`\mathbf{a}`:

.. math::

    \frac{\partial a_i}{\partial x_i}
    = J\,\frac{\partial}{\partial l_k}\!\left(J S^k_i a_i\right)
    = J\,\frac{\partial A^k}{\partial l_k},

where :math:`A^k = S^k_i a_i` is the contravariant flux of :math:`\mathbf{a}`.
For an incompressible velocity field this is equal to :math:`\partial V^q/\partial l_q = 0`.

**Divergence** of a rank-2 tensor :math:`\tau_{ij}`:

.. math::

    \frac{\partial}{\partial x_i}\left(\tau_{ij}\right)
    = J\,\frac{\partial}{\partial l_k} \left(S^k_i \tau_{ij}\right)
    = J\,\frac{\partial}{\partial l_k} T^k_j,

The quantity :math:`T^k_j` is the :math:`k`-th contravariant flux of the
:math:`j`-th row of the tensor.  Using the above formula, the divergence of the advection 
term in the momentum equation can be written as:

.. math::

    \frac{\partial}{\partial x_i}\!\left(u_i u_j\right)
    = J\,\frac{\partial}{\partial l_k} \left(S^k_i u_i u_j\right)
    = J\,\frac{\partial}{\partial l_k} \left(V^k u_j\right).

Using the above rules, the full viscous term in the momentum equation can be rewritten as:

.. math::

    \frac{\partial}{\partial x_j}\left[\nu_\text{eff}\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)\right] 
    = J\frac{\partial}{\partial l_k}\left[\nu_\text{eff}\left(JS^{kn}\frac{\partial u_i}{\partial l_n} + JS^k_jR_{ij}\right)\right].

where :math:`R_{ij}`  is a non-symmetric tensor defined as :math:`R_{ij}=S^m_i\partial u_j/\partial l_m`, and :math:`S^{kn}` is the contravariant metric tensor evaluated 
with the face area vectors, namely :math:`S^{kn} = S^k_i S^n_i = 1/J^2 \partial l_k/\partial x_i \partial l_n / \partial x_i`.

Governing Equations in Curvilinear Form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~

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
fluxes by inverting :math:`V^q = S^q_i u_i`.

Fractional-Step Projection Method
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pressure–velocity coupling uses a second-order **fractional-step projection**
method.  Given the state
:math:`(\mathbf{V}^n, p^n, \theta^n)` at time level :math:`n`, one time step
proceeds through the four stages below.

**Stage 1 — Momentum Predictor**: 

Advance the contravariant fluxes to a predicted state
:math:`\mathbf{V}^{*}` that satisfies the momentum equation but not yet the
divergence-free constraint.  Schematically, assuming to choose the Crank–Nicolson scheme (other time-stepping schemes, described in
:ref:`time-integration-section`, follow in a similar manner), the predictor step reads:

.. math::

    \frac{V^{q,*} - V^{q,n}}{\Delta t}
    = \frac{1}{2}\!\left[H^q_U(\mathbf{V}^{*},\mathbf{u}^{*})
                       + H^q_U(\mathbf{V}^n,\mathbf{u}^n)\right]
    + S^{q,n} + G^{q,n} + B^{q,*}.

Here :math:`H^q_U` contains the advection and viscous fluxes (see
:ref:`advection-schemes-section`); :math:`S^q` contains buoyancy, Coriolis,
turbine, and damping source terms; and :math:`G^q` is the discrete
pressure-gradient operator. The buoyancy source :math:`B^{q}_U` in the momentum equation depends on
the face-interpolated potential temperature through the Boussinesq density ratio.
Because :math:`\theta` is solved only in Stage 4, it cannot be used implicitly in
Stage 1. Hence, before the momentum predictor step, the buoyancy body-force vector
:math:`B^{q,n}` is computed from the current temperature field
:math:`\theta^n` and stored.  The previous step's buoyancy
:math:`B^{q,n-1}` is also retained, so that buoyancy can be extrapolated to the
mid-step level :math:`n+\tfrac{1}{2}` using an Adams–Bashforth 2 formula, which reads:

  .. math::

      B^{q,*} = \tfrac{3}{2}\,B^{q}(\theta^n)
                        - \tfrac{1}{2}\,B^{q}(\theta^{n-1}).

This approach reduces the temporal coupling error to
:math:`\mathcal{O}(\Delta t^2)` without the need to perform a temperature predictor step. 
On the first time step, a simple Forward Euler formula is used.

**Stage 2 — Pressure Correction**: 

The scalar correction :math:`\phi^{n+1}` is calculated as:

.. math::

    D\!\left(G^q(\phi^{n+1})\right) = -\frac{1}{\Delta t}\,D\!\left(V^{q,*}\right),

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
The constant mode of :math:`\phi` is removed after the solve by subtracting the domain-averaged value, ensuring that the pressure is defined up to a constant.

**Stage 3 — Contravariant Flux Projection**: 

The predicted contravariant fluxes are corrected to enforce mass conservation:

.. math::

    V^{q,n+1} = V^{q,*} + \Delta t\, G^q(\phi^{n+1}).

The pressure is updated as :math:`p^{n+1} = p^n + \phi^{n+1}`.

**Stage 4 — Temperature Solution**: 

With the divergence-free :math:`\mathbf{V}^{n+1}` available,
:math:`\theta^n` is advanced to  :math:`\theta^{n+1}` using the selected temperature
time-integration scheme (see :ref:`time-integration-section`).
