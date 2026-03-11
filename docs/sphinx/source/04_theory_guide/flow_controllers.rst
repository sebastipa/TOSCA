.. _controllers-section:

Flow Controllers and Mesoscale Coupling
---------------------------------------

TOSCA provides two classes of flow-driving strategies: *precursor controllers*, which sustain a statistically stationary ABL for spin-up or recycling purposes,
and *mesoscale-to-microscale coupling (MMC) methods*, which introduce time- and height-varying mesoscale conditions into the LES through internal forcing source
terms.  Both classes are configured via the ``controllerType`` (momentum) and ``controllerTypeT`` (temperature) keywords in ``ABLProperties.dat``.

.. list-table::
   :header-rows: 1
   :widths: 35 25 40

   * - Method
     - ``controllerType``
     - Purpose
   * - Geostrophic controller
     - ``geostrophic``
     - Fixed pressure gradient from geostrophic balance
   * - Pressure (PI) controller
     - ``pressure``
     - PI feedback to a target hub-height velocity
   * - Direct profile assimilation
     - ``directProfileAssimilation``
     - MMC: pointwise error feedback
   * - Indirect profile assimilation
     - ``indirectProfileAssimilation``
     - MMC: polynomial-smoothed error feedback
   * - Wavelet profile assimilation
     - ``waveletProfileAssimilation``
     - MMC: wavelet-filtered error feedback
   * - Geostrophic profile assimilation
     - ``geostrophicProfileAssimilation``
     - MMC: physics-based pressure-gradient forcing + WPA temperature
   * - Bulk velocity controller
     - ``-meanGradPForce 1`` (CLI flag)
     - Drive volume-averaged velocity to a target bulk value

Precursor Controllers
~~~~~~~~~~~~~~~~~~~~~~

Geostrophic Controller (``geostrophic``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The driving pressure gradient is derived from geostrophic balance:

.. math::
    :label: eq:geoBalance

    \frac{1}{\rho_0}\frac{\partial p_\infty}{\partial x} = f_c V_G,
    \qquad
    \frac{1}{\rho_0}\frac{\partial p_\infty}{\partial y} = -f_c U_G

where :math:`f_c = 2\Omega_z` is the Coriolis parameter and :math:`(U_G, V_G)` are the specified geostrophic wind components.  The geostrophic controller does
not directly control the hub-height wind speed, which is set by the turbulent stress balance inside the boundary layer.  An optional wind-angle controller
(Allaerts and Meyers 2015) can be added to prescribe the wind direction at a reference height by slowly rotating the domain.

Pressure (PI) Controller (``pressure``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pressure controller drives the horizontally averaged velocity at a reference height :math:`h_\text{ref}` to a target value :math:`u_{\text{ref},i}` using
a proportional–integral (PI) algorithm:

.. math::
    :label: eq:pControl

    \frac{1}{\rho_0}\frac{\partial p_\infty}{\partial x_i}
    = r\left(\alpha\, e_{P,i} + (1-\alpha)\,e^n_{I,i}\right),

with proportional and integral errors

.. math::

    e_{P,i} = \frac{u_{\text{ref},i} - \langle u_i(h_{\text{ref}})\rangle_{xy}}{\Delta t},
    \qquad
    e^n_{I,i} = \left(1 - \frac{\Delta t}{T}\right)e^{n-1}_{I,i} + \frac{\Delta t}{T}\,e_{P,i},

where :math:`r` is a relaxation factor, :math:`\alpha` the proportional fraction, and :math:`T` an integral time scale.

Because the geostrophic wind components :math:`(U_G, V_G)` are not prescribed a priori, an initial velocity inconsistent with geostrophic balance drives
inertial oscillations above the boundary layer at frequency :math:`f_c`.  These are suppressed by the **geostrophic damping** term described below.

Geostrophic Damping
"""""""""""""""""""

The momentum equations above the boundary layer are augmented with:

.. math::

    \begin{cases}
        \dfrac{\partial u}{\partial t} + 2\alpha_d f_c(u - U_G) + f_c(V_G - v) = 0 \\[4pt]
        \dfrac{\partial v}{\partial t} + 2\alpha_d f_c(v - V_G) - f_c(U_G - u) = 0
    \end{cases}

where :math:`\alpha_d > 0` determines the damping strength (critically damped at :math:`\alpha_d = 1`).  The e-folding decay time is
:math:`1/(2\alpha_d f_c)`.  For the amplitude to fall below 3% of its initial value the required damping time is
:math:`T_{3\%} = \ln(100/3)/(2\alpha_d f_c)`.

The geostrophic wind components are inferred from the PI pressure gradient via :eq:`eq:geoBalance`, filtered with a time constant :math:`0.2\pi/f_c`.  The
damping is blended to zero below the capping inversion using

.. math::

    f_d = \frac{1}{2}\left[1 + \tanh\!\left(\frac{7(h - H_d)}{\Delta_d}\right)\right],

where :math:`H_d` is the inversion center height and :math:`\Delta_d` the inversion width.

Temperature Controller (``initial``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To keep the mean potential temperature profile — and therefore the capping inversion — constant throughout a precursor simulation, a height-dependent source
term is added to the temperature equation:

.. math::
    :label: eq:temperatureControl

    s_\theta(h) = r\,\frac{\bar\theta(h) - \langle\theta(h)\rangle_{xy}}{\Delta t},

where :math:`\bar\theta(h)` is the target profile (taken from the initial condition) and :math:`r \approx 0.7` is a relaxation coefficient.  This controller
eliminates drift of the inversion layer caused by differences in numerical dissipation between codes or simulation durations.

Bulk Velocity Controller (``-meanGradPForce``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bulk velocity controller adds a uniform body force to the momentum equations every time step to drive the *volume-averaged* velocity toward a prescribed
target bulk velocity :math:`\mathbf{U}_\text{bulk}`.  Unlike the pressure PI controller, which acts on the planar average at a single reference height, this
controller acts on the full three-dimensional volume average:

.. math::

    \langle \mathbf{u} \rangle_V(t)
    = \frac{\displaystyle\sum_\Omega \mathbf{u}\,\delta V}
           {\displaystyle\sum_\Omega \delta V},

where the sums run over all fluid (non-IBM, non-overset) cells.  The per-step forcing applied uniformly to every cell is

.. math::

    \mathbf{F}_\text{bulk}(t) = \mathbf{U}_\text{bulk} - \langle \mathbf{u} \rangle_V(t),

which is a purely proportional correction with unit gain (gain = 1 s\ :sup:`-1`\ ).  The forcing is projected onto the curvilinear contravariant face fluxes so
that it is consistent with the collocated finite-volume discretisation of the non-orthogonal grid.

The target bulk velocity :math:`\mathbf{U}_\text{bulk}` is a three-component vector read from the ``uBulk`` keyword in the initial-condition input file (under
whichever initialisation block is active, e.g. ``ABLFlow``, ``uniform``, ``spreadInflow``).

A running diagnostic accumulator ``meanGradP`` tracks the cumulative correction applied over the simulation.

Mesoscale-to-Microscale Coupling (MMC)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Mesoscale-to-microscale coupling extends LES beyond statistically stationary, idealized conditions by incorporating time- and height-varying mesoscale
information through internal source terms :math:`F_{u_i}` and :math:`F_\theta` added to the momentum and temperature equations, respectively:

.. math::

    \frac{\partial \overline{u}_i}{\partial t} + \cdots = \cdots + F_{u_i}(z,t),
    \qquad
    \frac{\partial \overline{\theta}}{\partial t} + \cdots = \cdots + F_\theta(z,t).

All MMC methods read mesoscale data (wind velocity and virtual potential temperature as functions of height and time) from `inflowDatabase/mesoscaleData`.
For momentum, the region below a transition height can be treated differently from the region above, controlled by the ``lowerLayerForcingType`` keyword:

- ``constantHeight`` — transition height is a fixed value.
- ``ablHeight`` — transition height tracks the dynamically computed ABL height.
- ``mesoDataHeight`` — transition is set to the lowest height at which mesoscale data are available.

Profile Assimilation Methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All profile assimilation techniques compute the forcing from an error between the horizontally averaged LES field
:math:`\langle\bar\varphi_\text{LES}\rangle(z,t)` and the target mesoscale data :math:`\varphi_\text{m}(z,t)`:

.. math::

    e_\varphi(z,t) = \varphi_\text{m}(z,t) - \langle\bar\varphi_\text{LES}\rangle(z,t), \qquad \varphi \in \{u,\,v,\,\theta\}.

The three methods differ in how they filter or smooth :math:`e_\varphi` before applying it as a forcing.

Direct Profile Assimilation — DPA (``directProfileAssimilation``)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The error at each vertical level is used directly as a proportional forcing:

.. math::

    F_\varphi(z,t) = \mathcal{K}_p\,e_\varphi(z,t),

where :math:`\mathcal{K}_p = 0.2\,\mathrm{s^{-1}}` is the proportional gain.  Because the corrections at each height level are independent, DPA can generate
an artificial buildup of turbulent kinetic energy in regions where the LES has limited ability to adjust the mean flow through turbulent mixing.

Indirect Profile Assimilation — IPA (``indirectProfileAssimilation``)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

IPA smooths :math:`e_\varphi` by fitting a polynomial of order :math:`m` (set by ``polynomialOrder``) to the vertical error profile at each time step.
The polynomial coefficients :math:`\hat{\boldsymbol\beta}` are found by weighted least-squares minimisation:

.. math::

    \hat{\boldsymbol\beta}(t) = \left(Z^\top W Z\right)^{-1} Z^\top W\,e_\varphi(t),

where :math:`Z` is the Vandermonde matrix evaluated at the vertical grid points and :math:`W` is a diagonal weight matrix.  The forcing is then

.. math::

    F_\varphi(z,t) = \mathcal{K}_p\,Z\hat{\boldsymbol\beta}(t).

The polynomial fit acts as a vertical low-pass filter that enforces inter-level coherence.  A typical choice is a cubic polynomial (:math:`m=3`) with uniform
weights, but higher orders can be used to capture stronger vertical gradients.  IPA reduces the TKE build-up of DPA but can over-smooth sharp features such as
low-level jets.

Wavelet Profile Assimilation — WPA (``waveletProfileAssimilation``)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

WPA applies a multi-resolution analysis (MRA) to the error profile using the discrete wavelet transform (DWT), providing a scale-aware filtering that
preserves vertical structure while avoiding the over-smoothing of IPA.  The error is decomposed across :math:`m_0` levels (set by ``waveletDecompLevel``) into
approximation and detail coefficients using a chosen mother wavelet (``waveletName``, e.g. ``db7`` — Daubechies of order 7):

.. math::

    e_\varphi(k,t) = \underbrace{\sum_n \mathbb{W}^A_{m_0,n}(t)\,\phi_{m_0,n}(k)}_{\text{low-freq. approx.}}
                   + \sum_{m=0}^{m_0-1}\sum_n \mathbb{W}^D_{m,n}(t)\,\psi_{m,n}(k).

By retaining only the low-frequency approximation coefficients :math:`\mathbb{W}^A_{m_0,n}` (controlled by ``waveletTMethod``), the reconstructed forcing is

.. math::

    F_\varphi(z,t) = \mathcal{K}_p \sum_n \mathbb{W}^A_{m_0,n}(t)\,\phi_{m_0,n}(z).

WPA requires that TOSCA is compiled with Python support (``USE_PYTHON=1`` in the Makefile).  Symmetric signal extension is applied at the domain boundaries to
minimise edge artefacts, which were observed in IPA near the surface.

Geostrophic Profile Assimilation — HGWPA (``geostrophicProfileAssimilation``)
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The hybrid geostrophic-wavelet profile assimilation (HGWPA) method computes the momentum source term from the large-scale pressure gradient force rather than
from an error-based correction.  This is more physically consistent than DPA/IPA/WPA because it avoids momentum error forcing and directly encodes the
relevant physical processes (pressure gradient, Coriolis force, and turbulent friction) into the LES source term.  The temperature source term is computed with
WPA.

Under baroclinic, unsteady conditions the large-scale pressure gradient is

.. math::
    :label: eq:geoBalanceUnsteadyMMC

    -\frac{\partial p_\infty(z,t)}{\partial x}
    = \frac{dG_x(z,t)}{dt} - f_c\,G_y(z,t) - \frac{\partial\langle\mathcal{T}_{xz}(z,t)\rangle}{\partial z},

and analogously for the :math:`y` component, where :math:`G_x(z,t), G_y(z,t)` are the mesoscale (possibly height-varying) geostrophic wind components,
and :math:`\mathcal{T}_{i3}` is the total shear stress (resolved Reynolds stress + SGS stress):

.. math::

    \mathcal{T}_{i3}(z) = -\langle\overline{u'_i u'_3}\rangle - \tau_{i3}^\text{sgs}(z).

The frictional drag term (last term in :eq:`eq:geoBalanceUnsteadyMMC`) is computed directly from the running LES field using horizontal plane averaging and
temporal filtering.  HGWPA requires ``lowerLayerForcingType = ablHeight``.

The frictional drag contribution is active only inside the boundary layer and is blended smoothly to zero above the ABL height using a hyperbolic-tangent
taper:

.. math::

    f_s(z) = 1 - \frac{1}{2}\!\left[1 + \tanh\!\left(\frac{7(z - H_t)}{\Delta_s}\right)\right],

where :math:`H_t = H + 0.5\,\Delta_s`, :math:`H` is the instantaneous ABL height, and :math:`\Delta_s = 0.2H` is the transition width.  The ABL height is
computed dynamically at every time step from the running LES statistics.  The stability regime is detected from the surface heat flux:

- **Convective** (surface heat flux :math:`> 0.02\,\text{K\,m\,s}^{-1}`): :math:`H` is the height of the minimum of the horizontally averaged total heat flux
  profile above its near-surface maximum.  If no minimum is found, the height of the maximum vertical temperature gradient is used as a fallback.
- **Stable/neutral**: :math:`H` is the height at which the total (resolved + SGS) shear stress drops below 10 % of its surface-layer maximum,
  cross-checked against the bulk Richardson number criterion :math:`\text{Ri}_B \geq 0.25`.

The raw estimate is smoothed with an exponential running average

.. math::

    H^{n+1} = \frac{T_H}{T_H + \Delta t}\,H^n + \frac{\Delta t}{T_H + \Delta t}\,H_\text{raw}^{n+1},

where the averaging time scale :math:`T_H` is set by the ``hAverageTime`` keyword in the ``controllerProperties`` sub-dictionary of ``ABLProperties.dat``.