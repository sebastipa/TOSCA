.. _advection-schemes-section:

Advection Schemes
-----------------

The advection (divergence) term in both the momentum and temperature equations is
discretized using one of six selectable schemes, chosen via the keyword
``-divScheme`` in ``control.dat``.  All schemes operate on a one-dimensional,
four-point stencil.  Assuming a positive advecting contravariant
flux (the stencil mirrors for the opposite sign), for a cell face located between cells :math:`i` and 
:math:`i+1`, the fluxes at the upwind and downind nodes are :math:`f_U = f_i` and :math:`f_D = f_{i+1}`, respectively, 
while those at the far-upwind and far-downwind nodes are :math:`f_{UU} = f_{i-1}` and
:math:`f_{DD} = f_{i+2}`, respectively. Using this notation, the six available schemes are:

**central** — Second-order linear interpolation (non-dissipative, dispersive):

.. math::

    f_{i+1/2} = \frac{f_U + f_D}{2}.

**central4** — Fourth-order central interpolation (4-point stencil):

.. math::

    f_{i+1/2} = \frac{-f_{UU} + 9f_U + 9f_D - f_{DD}}{16}.

A mild explicit biharmonic hyperviscosity correction :math:`-\nu_{h}(f_{UU} - 4f_U + 6f_{i+1/2} - 4f_D + f_{DD})`
is applied to damp Nyquist-wavenumber ringing (coefficient :math:`\nu_h` set by ``-hyperVisc``).

**quickDiv** — Third-order QUICK (Quadratic Upstream Interpolation for ConvectiveKinematics, Leonard 1979) for uniform mesh:

.. math::

    f_{i+1/2} = \frac{f_U + f_D}{2} + \frac{2f_U - f_D - f_{UU}}{8}.

The second term is a third-order correction to the linear interpolant, which biases the stencil upstream.

**centralUpwind** — Van Leer-limited blend of central and QUICK for uniform mesh.  The face value is

.. math::

    f_{i+1/2}
    = \underbrace{\frac{f_U + f_D}{2}}_{\text{central}}
    + \bigl(1 - \psi(r)\bigr)
      \underbrace{\frac{2f_U - f_D - f_{UU}}{8}}_{\text{QUICK correction}},

where the **Van Leer limiter** and the upwind-ratio :math:`r` are

.. math::

    \psi(r) = \frac{r + |r|}{2(1 + |r|)},
    \qquad
    r = \frac{f_U - f_{UU}}{f_D - f_U}.

When :math:`r < 0` (the solution has a local extremum) :math:`\psi = 0` and the full QUICK upwind correction is applied, giving TVD behavior.  When
:math:`r \to \infty` (the solution is smooth and gently varying) :math:`\psi \to 1` and the scheme reduces to pure central interpolation, recovering low-dissipation
properties on well-resolved scales.  The Van Leer limiter uses contravariant flux ratios rather than scalar ratios for the momentum equation, providing
flow-direction awareness on non-Cartesian grids.

This scheme is the default and recommended choice for most LES simulations, as it combines low numerical dissipation in smooth regions with robust behavior near
sharp gradients or steep velocity layers.

**centralUpwindW** — Identical to ``centralUpwind`` but with cell-distance-weighted interpolation coefficients to preserve second-order
accuracy on non-uniform (stretched) meshes.

**weno3** — Third-order WENO scheme with two candidate stencils (Liu et al., 1994).  The two candidate interpolants are

.. math::

    u_0 = \frac{f_U + f_D}{2}, \qquad u_1 = \frac{3f_U - f_{UU}}{2},

with ideal convex weights :math:`d_0 = 2/3`, :math:`d_1 = 1/3`. The smoothness indicators are

.. math::

    \beta_0 = (f_U - f_D)^2, \qquad \beta_1 = (f_{UU} - f_U)^2,

and the nonlinear WENO weights are

.. math::

    \omega_k = \frac{\alpha_k}{\alpha_0 + \alpha_1},
    \qquad
    \alpha_k = \frac{d_k}{(\varepsilon + \beta_k)^2},
    \qquad \varepsilon = 10^{-6}.

The final face value is :math:`f_{i+1/2} = \omega_0 u_0 + \omega_1 u_1`. Near discontinuities :math:`\beta_0 \gg \beta_1`, driving :math:`\omega_0 \to 0`
so that the upwind stencil dominates.  In smooth regions the weights recover the optimal third-order interpolant.

.. note::

    Viscous-stress and pressure-gradient terms always use pure second-order central differences, regardless of the advection scheme selected.
    The Laplacian on the curvilinear mesh couples all three coordinate directions through the off-diagonal metric tensor components :math:`S^{qr}` (see
    :ref:`numerics-section`).
