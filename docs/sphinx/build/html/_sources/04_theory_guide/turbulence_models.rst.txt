
.. _sgs-model-section:

Sub-grid Scale Models
---------------------

TOSCA provides ten SGS models that fall into three categories: *functional* (eddy-viscosity) models, which reproduce the forward energy cascade but
discard structural information, purely *structural* models, and *mixed* models, which combine both approaches. All models operate in generalized curvilinear
coordinates following the partial-transformation framework of Ge and Sotiropoulos (2007) and the GCC extension of Armenio and Piomelli (2000). SGS turbulence models
are selected with the keyword ``-lesModel`` keyword in the ``control.dat`` file.

Under the linear eddy-viscosity hypothesis, the deviatoric part of the SGS stress tensor is

.. math::

    \tau_{ij} - \tfrac{1}{3}\delta_{ij}\tau_{kk}
    = -2\nu_t S_{ij}
    = -2\nu_t\,\tfrac{1}{2}\!\left(\frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i}\right),

and its contravariant projection (used inside the momentum equation) is

.. math::

    T^r_i = \frac{1}{J}\frac{\partial l^r}{\partial x_j}\tau_{ij}.

The resolved strain-rate magnitude, required by most models, is :math:`|\overline{S}| = \sqrt{2S_{ij}S_{ij}}`, and the LES filter width is the local cell
volume cubic root :math:`\Delta = J^{-1/3}`.

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Model name
     - keyword
     - Category
   * - Smagorinsky
     - ``smagorinsky``
     - Functional — fixed-coefficient eddy viscosity
   * - Stability-Dependent
     - ``stabilityDependent``
     - Functional — stability-corrected eddy viscosity
   * - Dynamic Smagorinsky (box)
     - ``dynamicSmagorinsky``
     - Functional — dynamically computed :math:`C_s`, local box averaging
   * - Dynamic Lagrangian Scale-Invariant
     - ``dynamicLASI``
     - Functional — dynamically computed :math:`C_s`, Lagrangian averaging, scale-invariant
   * - Dynamic Lagrangian Scale-Dependent
     - ``dynamicLASD``
     - Functional — dynamically computed :math:`C_s`, Lagrangian averaging, scale-dependent
   * - Dynamic Plane-Averaged Scale-Dependent
     - ``dynamicPASD``
     - Functional — dynamically computed :math:`C_s`, plane averaging, scale-dependent
   * - Anisotropic Minimum Dissipation (AMD)
     - ``amd``
     - Functional — minimum-dissipation eddy viscosity
   * - Vreman
     - ``vreman``
     - Functional — frame-invariant eddy viscosity
   * - Bardina–Vreman (BV)
     - ``bardinaVreman``
     - Mixed — Bardina scale-similarity (structural) + Vreman eddy viscosity
   * - Bardina–AMD (BAMD)
     - ``bardinaAMD``
     - Mixed — Bardina scale-similarity (structural) + AMD eddy viscosity

Functional Models
~~~~~~~~~~~~~~~~~

Smagorinsky Model
^^^^^^^^^^^^^^^^^

The classical Smagorinsky (1963) model prescribes a fixed coefficient :math:`C_s = 0.17` (corresponding to the TOSCA default ``std_cs`` = 0.0289 when
expressed as :math:`C_s^2`):

.. math::

    \nu_t = \left(C_s\,\Delta\right)^2\,|\overline{S}|.

The coefficient is assigned globally for the entire domain.  This model generally over-predicts SGS dissipation near solid walls.  

Stability-Dependent Smagorinsky Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In stratified ABL flows the SGS length scale must shrink under stable conditions.  The stability-dependent model adjusts the effective mixing length using the
local gradient Richardson number :math:`Ri = N^2/|\overline{S}|^2`, where :math:`N = \sqrt{(g/\theta_0)\,\partial\theta/\partial z}` is the Brunt–Väisälä
frequency.  Three stability regimes are distinguished:

- :math:`Ri \geq Ri_c = 0.23` (strongly stable): turbulence is suppressed; the length scale is set effectively to zero.
- :math:`Ri \leq Ri_s = 0.059` (neutral or weakly stable): the mixing length equals the grid filter width, :math:`l = \Delta`.
- :math:`Ri_s < Ri < Ri_c` (intermediate stability): the length scale :math:`l` decreases smoothly according to a TKE-based parametrization, and is
  additionally bounded by the von Kármán length scale :math:`0.41\,z` near the surface.

The effective constant is then

.. math::

    C_s^2 = \frac{c_m^3}{c_\varepsilon}\frac{l^2}{\Delta^2}\sqrt{1 - \frac{c_h}{c_m}Ri},

where :math:`c_m, c_h, c_\varepsilon` are model constants.  The temperature diffusivity follows from the stability-dependent turbulent Prandtl number
described in the :ref:`sgs-temp-flux` section below.

Dynamic Smagorinsky Models
^^^^^^^^^^^^^^^^^^^^^^^^^^

All dynamic variants determine :math:`C_s(\mathbf{x},t)` from the resolved flow using the Germano identity (Germano et al. 1991, Lilly 1992).  A test filter
:math:`\widetilde{(\cdot)}` at scale :math:`\widetilde{\Delta} = 3\Delta` separates the resolved Leonard stress tensor

.. math::

    L^r_i = \widetilde{\overline{V}^r\,\overline{u}_i} - \widetilde{\overline{V}^r}\,\widetilde{\overline{u}}_i,

computable directly from resolved quantities from the model prediction.  The forcing tensor is

.. math::

    M^r_i = 2\left[
        \widetilde{\overline{\Delta}}^2\,|\widetilde{\overline{S}}|\,\widetilde{S}^r_i
        - \Delta^2\,\widetilde{|\overline{S}|\,S^r_i}
    \right].

Because the tensors implicitly contain the face-area vectors :math:`\partial l^r/\partial x_j`, Galilean invariance requires contracting them with the
covariant metric tensor :math:`g_{ij}` before the minimisation.  The scale-invariant dynamic coefficient is

.. math::

    C_s^2(\mathbf{x},t) =
    \frac{\langle L^r_i M^m_i\,g_{rm}\rangle}{\langle M^q_j M^n_j\,g_{qn}\rangle}.

Four averaging strategies for the angle brackets are available:

**DSM** (``dynamicSmagorinsky``): local box averaging over the 3×3×3 stencil with a Simpson integration kernel (Ghosal et al. 1995), best for homogeneous
or weakly inhomogeneous flows.

**DLASI** (``dynamicLASI``): Lagrangian averaging along fluid path lines (Meneveau et al. 1996), assumes *scale invariance* (:math:`C_s` is the same at
the grid and test-filter scales).  The averaged numerator :math:`I_{LM}` and denominator :math:`I_{MM}` evolve in time via the relaxation

.. math::

    I^n_{LM} = \epsilon\,\langle L^r_i M^m_i g_{rm}\rangle
    + (1-\epsilon)\,I^{n-1}_{LM}(\mathbf{x}-\mathbf{u}\Delta t),
    \quad
    \epsilon = \frac{\Delta t/T_s}{1+\Delta t/T_s},

with Lagrangian time scale :math:`T_s = 1.5\Delta\,[|I_{LM}\,I_{MM}|]^{-1/8} + \varepsilon`.  The lookup uses trilinear backward-in-time interpolation.
This is one of the recommended combination for production LES runs in TOSCA.

**DLASD** (``dynamicLASD``) — Lagrangian averaging with *scale dependence* (Porté-Agel et al. 2000, Bou-Zeid et al. 2005).  A second, coarser test filter at
scale :math:`4\Delta` provides the scale-dependence parameter

.. math::

    \beta = \frac{C^2_{s,4\Delta}}{C^2_{s,2\Delta}},
    \quad
    C^2_{s,\Delta} = \frac{C^2_{s,2\Delta}}{\max(\beta,\,0.125)}.

:math:`\beta` is clipped at 0.125 to prevent numerical instability in regions where the scale similarity assumption breaks down.
This is one of the recommended combination for production LES runs in TOSCA.

**DPASD** (``dynamicPASD``) — same scale-dependent formulation as DLASD, but uses horizontal plane averaging instead of Lagrangian path-line averaging.
Suitable for flows with a clear direction of statistical homogeneity (e.g., flat-terrain atmospheric boundary layers with no turbines).

Anisotropic Minimum Dissipation (AMD) Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The AMD model (Rozema et al. 2015, Abkar and Moin 2016) enforces that the sub-filter energy must not grow in time anywhere in the domain, using a Poincaré
inequality for the box-filter energy.  The resulting eddy viscosity is

.. math::

    \nu_t = \max\!\left(
    \frac{-\,(\hat\partial_k\overline{u}_i)(\hat\partial_k\overline{u}_j)\,\overline{S}_{ij}}
         {(\partial_l\overline{u}_m)(\partial_l\overline{u}_m)},\;0
    \right),

where :math:`\hat\partial_j = \sqrt{C_p}\,\delta_j\,\partial_j` is the scaled gradient operator (:math:`\delta_j` the cell width in direction :math:`j`,
:math:`C_p` a modified Poincaré constant set by ``-amdCs``, with default value 0.1).  The model is purely dissipative (no backscatter), low-overhead (no filtering),
and well-suited to wall-bounded flows because :math:`\nu_t \to 0` in the absence of velocity gradients.  In curvilinear coordinates the gradients are
computed in curvilinear directions and then projected to Cartesian components using the metric terms.

Vreman Model
^^^^^^^^^^^^

The Vreman model (Vreman 2004) offers an alternative frame-invariant eddy viscosity that scales correctly near walls and requires no dynamic procedure.  The
eddy viscosity is

.. math::

    \nu_t = 2.5\,C_s^2\sqrt{\frac{B_{\beta_V}}{A_{ij}A_{ij}}},

where :math:`A_{ij} = \partial\overline{u}_i/\partial x_j` is the velocity-gradient tensor, :math:`\beta_{V,ij} = \Delta^2 A_{mi}A_{mj}`, and

.. math::

    B_{\beta_V} = \beta_{V,11}\beta_{V,22} - \beta_{V,12}^2
               + \beta_{V,11}\beta_{V,33} - \beta_{V,13}^2
               + \beta_{V,22}\beta_{V,33} - \beta_{V,23}^2.

The constant used in TOSCA is ``vreman_cs`` = 0.07225.  The model vanishes identically in laminar shear and in two-component flows, making it suitable for
transitional or weakly turbulent regions.

Mixed Models
~~~~~~~~~~~~

Mixed models combine a *structural* scale-similarity term (Bardina et al. 1980, 1983) with a functional eddy-viscosity term.  The Bardina term captures the
anisotropic and directional features of the SGS stress tensor but lacks sufficient dissipation on its own.  Adding an eddy-viscosity closure restores the
correct net energy transfer to the sub-grid scales.

The SGS stress tensor in curvilinear form is split as

.. math::

    T^r_i = \underbrace{{T^r_i}^{\alpha,\text{dev}}}_{\text{eddy viscosity}}
          + \underbrace{{T^r_i}^{\beta}}_{\text{Bardina}},

where the structural component is

.. math::

    {T^r_i}^{\beta} = C_b\!\left(
        \widetilde{\overline{V}^r\,\overline{u}_i}
        - \widetilde{\overline{V}^r}\,\widetilde{\overline{u}}_i
    \right),

with :math:`C_b = 1.0` (Galilean invariance, Speziale 1985).  The primary filter corresponds to the grid scale :math:`\Delta` and the secondary test filter
has width :math:`2\Delta`.

Bardina–Vreman (BV) Model
^^^^^^^^^^^^^^^^^^^^^^^^^

``bardinaVreman`` — the Bardina structural term is combined with the Vreman eddy viscosity.  A fixed Smagorinsky-like coefficient ``vreman_cs`` = 0.07225 is
used for the functional part.  This model provides a baseline mixed formulation that is computationally efficient (avoiding a dynamic procedure) and suitable
for flows where a scale-invariant coefficient is a reasonable approximation.

Bardina–AMD (BAMD) Model
^^^^^^^^^^^^^^^^^^^^^^^^

``bardinaAMD`` — the Bardina structural term is combined with the AMD eddy viscosity.  The AMD coefficient is set via ``-amdCs`` (with default value 0.1).  Among all
mixed models tested in the validation study (published in JCP, Ajay et al. 2025), BAMD combined with the ``-central4`` advection scheme (with 0.8 ``-hyperVisc``) shows the
best balance between accuracy and computational efficiency across Taylor–Green vortex flow, turbulent channel flow, and ABL flow over heterogeneous terrain test cases.
This is one of the recommended combination for production LES runs in TOSCA.

.. _sgs-temp-flux:

SGS Temperature Closure
~~~~~~~~~~~~~~~~~~~~~~~

Sub-grid scale temperature fluxes are modeled following Moeng (1984) through a turbulent thermal diffusivity :math:`\kappa_t = \nu_t / Pr_t`, where the
turbulent Prandtl number depends on the local stability:

.. math::

    Pr_t = \frac{1}{1 + 2l/\Delta}

.. math::

    l = \begin{cases}
        \min\!\left(\dfrac{7.6\,\nu_t}{\Delta}\sqrt{\theta_0/|s|},\;\Delta\right) & s < 0 \\
        \Delta & s \geq 0
    \end{cases}

where :math:`s = g_i\,\partial\theta/\partial x_i` is the stratification parameter. Under neutral or unstable conditions (:math:`s \geq 0`) the mixing
length equals the filter width and :math:`Pr_t = 1/3`, reflecting increased scalar mixing. Under stable conditions (:math:`s < 0`) the mixing length
decreases and :math:`Pr_t \to 1`, consistent with the suppression of vertical mixing by stratification.