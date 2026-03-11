.. _turbine-models-section:

Turbine Models
--------------

TOSCA implements four actuator-based wind turbine models at increasing levels of fidelity and cost:

.. list-table::
   :header-rows: 1
   :widths: 18 14 14 14 14 14 12

   * - Model
     - BEM aero
     - Blade rotation
     - Ct-based
     - Drive-train ODE
     - OpenFAST coupling
     - Tower/nacelle
   * - ADM
     - yes
     - no (axisymmetric)
     - no
     - yes
     - yes
     - optional
   * - ALM
     - yes
     - yes
     - no
     - yes
     - yes
     - optional
   * - uniformADM
     - no
     - no
     - yes
     - no
     - partial
     - optional
   * - AFM
     - no
     - no
     - yes
     - no
     - partial
     - optional

All models represent the rotor as a distribution of Lagrangian force points whose body forces are projected onto the surrounding CFD mesh using an isotropic Gaussian kernel. Optionally, every turbine in a farm can delegate its structural dynamics and detailed aerodynamics to OpenFAST while TOSCA handles the background LES flow.


Gaussian Force Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~

The Lagrangian body force at each actuator point is smeared onto the surrounding mesh cells using the regularised delta function

.. math::
    :label: eq:gaussianProjection

    g(\mathbf{x}) = \frac{1}{\varepsilon^3 \pi^{3/2}}
    \exp\!\left(-\frac{|\mathbf{x}-\mathbf{x}_0|^2}{\varepsilon^2}\right),

where :math:`\mathbf{x}_0` is the actuator point location and :math:`\varepsilon` is the regularisation width. The integration domain is truncated once 99 % of the Gaussian volume has been accounted for. To avoid systematic under-projection and numerical instabilities, :math:`\varepsilon/\Delta > 2` is recommended, where :math:`\Delta` is the local mesh size.


Aerodynamic Force Computation (ADM and ALM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both the ADM and the ALM share the same blade-element aerodynamic routine. This section describes the common formulation; model-specific aspects are covered below.

**Blade reference frame.** At each actuator point with position :math:`\mathbf{r}_p` relative to the rotor centre, a local reference frame :math:`(\hat{\mathbf{x}}_b, \hat{\mathbf{y}}_b, \hat{\mathbf{z}}_b)` is constructed:

.. math::

   \hat{\mathbf{z}}_b = \pm\,\hat{\mathbf{r}}_p, \qquad
   \hat{\mathbf{x}}_b = -\hat{\boldsymbol{\Omega}}, \qquad
   \hat{\mathbf{y}}_b = \hat{\mathbf{z}}_b \times \hat{\mathbf{x}}_b,

where :math:`\hat{\boldsymbol{\Omega}}` is the unit rotor-axis vector (positive from nacelle to wind), and the sign of :math:`\hat{\mathbf{z}}_b` follows the rotation convention: for a clockwise rotor :math:`\hat{\mathbf{z}}_b = +\hat{\mathbf{r}}_p` (root to tip), for counter-clockwise :math:`\hat{\mathbf{z}}_b = -\hat{\mathbf{r}}_p`. This convention keeps the local velocity vector in the first quadrant of the :math:`(x_b, y_b)` plane in both cases. The CFD velocity :math:`\mathbf{u}` is projected onto this frame to give the axial component :math:`u_{b,x}` (aligned with wind direction) and tangential component :math:`u_{b,y}` (aligned with blade motion).

**Angle of attack.** The geometric angle of attack at radial station :math:`r` is

.. math::
    :label: eq:aoa

    \alpha = \arctan\!\left(\frac{u_{b,x}}{u_{b,y}}\right) - \theta_t - \beta_c - \Delta\beta_\text{wf},

where :math:`\theta_t` is the local pre-twist, :math:`\beta_c` is the collective blade pitch, and :math:`\Delta\beta_\text{wf}` is an optional wind-farm pitch offset provided by the external scheduler.

**Airfoil table lookup.** The lift and drag coefficients are obtained by 2-D interpolation in pre-supplied tabular data. Each radial station :math:`r` lies between two airfoil tables at radial positions :math:`r_1` and :math:`r_2`, and the coefficients are linearly interpolated in both the radial direction and in angle of attack:

.. math::

   C_{l/d}(r, \alpha) = w_1\,C_{l/d}^{(1)}(\alpha) + w_2\,C_{l/d}^{(2)}(\alpha),
   \qquad w_1 = \frac{r_2 - r}{r_2 - r_1},\; w_2 = 1 - w_1.

**Prandtl tip/root-loss correction.** To account for the finite number of blades, the Prandtl correction factor :math:`F = F_\text{tip} \cdot F_\text{root}` is applied, where

.. math::
    :label: eq:prandtl

    F_\text{tip}  &= \frac{2}{\pi}\arccos\!\exp\!\left(-\frac{N_b}{2}\,
                     \frac{R-r}{r\sin\varphi}\right), \\
    F_\text{root} &= \frac{2}{\pi}\arccos\!\exp\!\left(-\frac{N_b}{2}\,
                     \frac{r-r_\text{hub}}{r\sin\varphi}\right),

with :math:`N_b` the number of blades, :math:`R` the tip radius, :math:`r_\text{hub}` the hub radius, and :math:`\varphi` the local flow angle. The corrected coefficients are :math:`\tilde{C}_{l/d} = F\,C_{l/d}`.

**Lift and drag forces.** The elemental lift and drag at a blade segment of chord :math:`c` and span :math:`dr` are

.. math::
    :label: eq:liftdrag

    L = \tfrac{1}{2}|\mathbf{u}|^2\,c\,dr\,\tilde{C}_l, \qquad
    D = \tfrac{1}{2}|\mathbf{u}|^2\,c\,dr\,\tilde{C}_d.

The aerodynamic reference frame :math:`(\hat{\mathbf{x}}_a, \hat{\mathbf{y}}_a, \hat{\mathbf{z}}_a)` places :math:`\hat{\mathbf{x}}_a` along the local velocity vector, :math:`\hat{\mathbf{y}}_a` in the radial direction, and :math:`\hat{\mathbf{z}}_a = \hat{\mathbf{x}}_a \times \hat{\mathbf{y}}_a` as the lift direction. The body force imposed on the flow (reaction principle) at point :math:`p` is

.. math::
    :label: eq:bodyforce

    \mathbf{B}_p = -(L\,\hat{\mathbf{z}}_a + D\,\hat{\mathbf{x}}_a).

**Rotor loads.** Thrust and torque are accumulated across all actuator points weighted by their solidity :math:`\sigma_p`:

.. math::
    :label: eq:thrust_torque

    T &= \sum_p f_{a,p}\,\sigma_p\,\rho, \\
    Q &= \sum_p f_{t,p}\,\sigma_p\,\rho\,r_p\cos\psi_\text{cone},

where :math:`f_{a,p}` is the axial force component, :math:`f_{t,p}` the tangential component, :math:`\rho` the air density, and :math:`\psi_\text{cone}` the pre-cone angle. The aerodynamic power is :math:`P_\text{aero} = Q\,\Omega`, and the electrical power is :math:`P_\text{gen} = Q_\text{gen}\,\Omega_\text{gen}\,\eta_\text{gen}`.


Actuator Disk Model (ADM)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ADM represents each rotor as a fixed axisymmetric disk of actuator points distributed in the radial and azimuthal directions. Because the disk is time-steady, the forcing represents the azimuthal average of the blade loading and is suitable for wake studies that do not require blade-passage unsteadiness.

Each point on the ADM disk is assigned a solidity weight that distributes the discrete blade loading over the full revolution. The BEM aerodynamic routine described above is applied at every disk point with the local pitch, chord, and twist from the blade geometry tables.

**DIPC helix pitch control.** The ADM supports an individual-blade-pitch strategy that creates a helical wake structure (Dynamic Induction and Phase Control, DIPC). An oscillating pitch offset is superimposed on the collective pitch at each azimuthal position :math:`\psi`:

.. math::
    :label: eq:dipc

    \Delta\beta_\text{DIPC}(\psi, t) = A\sin\!\left(\psi \pm 2\pi f t\right),

where :math:`A` is the amplitude (degrees), :math:`f` is the helix frequency (Hz), and the sign selects the helix handedness. This option is not available for the ALM.


Actuator Line Model (ALM)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ALM represents each rotor blade as a radial line of actuator points that physically rotates with the rotor at every time step. This makes it an inherently unsteady three-dimensional model capable of resolving blade-passage effects and tip vortices.

At each time step the blade azimuth :math:`\psi` is advanced by :math:`\Omega\,\Delta t` (from the drive-train solver described below), and each actuator point is rotated accordingly. The BEM aerodynamic routine is then applied at every point using the instantaneous sampled CFD velocity. The azimuth is synchronised to the master writer rank at checkpoint write times via a single global reduction.

The ALM uses the same aerodynamic and control routines as the ADM; the DIPC helix option is not implemented for rotating blades.


Uniform Actuator Disk Model (uniformADM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The uniformADM does not require detailed blade geometry. The rotor thrust is calculated from a prescribed thrust coefficient :math:`C_t`, which may be constant or variable (tabulated as a function of a reference velocity). The body force at the :math:`p`-th disk element is

.. math::
    :label: eq:uniformADM

    \mathbf{f}_p = \tfrac{1}{2}\,U_\text{ref}^2\,dA_p\,C_t\,(-\hat{\boldsymbol{\Omega}}),

where :math:`dA_p` is the rotor area fraction of point :math:`p` and :math:`U_\text{ref}` is a reference velocity whose computation depends on the chosen sampling strategy (``sampleType``):

- ``rotorUpstream``: :math:`U_\text{ref}` is sampled at the upstream face of the rotor disk.
- ``givenVelocity``: :math:`U_\text{ref}` is set directly by the user.
- ``rotorDisk``: the disk-averaged velocity :math:`\bar{U}_\text{disk}` is measured on the rotor plane and :math:`U_\text{ref}` is recovered via momentum theory, :math:`U_\text{ref} = \bar{U}_\text{disk}/(1-a)`, where :math:`a` satisfies :math:`C_t = 4a(1-a)`. When :math:`C_t` is variable this requires an iterative solve.

When a variable-:math:`C_p` table is provided, the aerodynamic power is computed as :math:`P_\text{aero} = \tfrac{1}{2}\,U_\text{ref}^3\,A_\text{rotor}\,C_p(U_\text{ref})`; otherwise it is derived from thrust and induction: :math:`P_\text{aero} = T\,(1-a)\,U_\text{ref}`.


Actuator Farm Model (AFM)
~~~~~~~~~~~~~~~~~~~~~~~~~~

The AFM is a single-point coarse-LES model intended for very large wind farms where resolving the rotor disk is unnecessary. The entire rotor is represented by a single actuator point at the hub location, and the same Ct-based thrust formulation as the uniformADM is applied. Three sampling strategies are available:

- ``momentumTheory``: the hub-height velocity is used to back-compute :math:`U_\text{ref}` from one-dimensional momentum theory.
- ``rotorDisk``: disk-averaged velocity is measured on the rotor plane and corrected for induction.
- ``integral``: :math:`U_\text{ref}` is computed from the area-integrated velocity over the upstream rotor face.


Tower and Nacelle Models
~~~~~~~~~~~~~~~~~~~~~~~~~

Optional tower and nacelle drag models can be activated independently for any turbine model (``includeTwr``, ``includeNacelle``). The tower is represented as a series of actuator points along its axis; the nacelle as a single point at the hub. The drag body force at each point follows a standard drag formulation and is projected onto the mesh via the same Gaussian kernel as the blade forces.


Rotor Dynamics and Drive Train
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the ADM and ALM (without OpenFAST), TOSCA solves a lumped-parameter drive-train model. The equation of motion for the low-speed shaft is

.. math::
    :label: eq:drivetrain

    J\,\frac{d\Omega}{dt} = \eta_\text{gbx}\,Q_\text{aero} - i_\text{gbx}\,Q_\text{gen},

where :math:`J` is the total drive-train inertia referenced to the low-speed shaft (``driveTrainInertia``), :math:`\eta_\text{gbx}` is the gearbox efficiency (``gbxEff``), :math:`i_\text{gbx}` is the gearbox ratio (``gbxRatioG2R``, rotor-to-generator), :math:`Q_\text{aero}` is the aerodynamic torque accumulated from the BEM loop, and :math:`Q_\text{gen}` is the generator torque. The generator speed is :math:`\Omega_\text{gen} = \Omega\,i_\text{gbx}`.

A first-order low-pass filter smooths the rotor speed signal used by the generator controller:

.. math::
    :label: eq:lowpass

    a &= e^{-\Delta t\,f_\text{filt}}, \\
    \Omega_\text{filt}(t) &= (1-a)\,\Omega(t) + a\,\Omega_\text{filt}(t-\Delta t),

where :math:`f_\text{filt}` is the filter cut-off frequency (``rtrSpdFilterFreq``, rad/s).


Control Systems
~~~~~~~~~~~~~~~

**Generator torque controller.**
TOSCA implements a classic four-region generator torque controller. Below the cut-in generator speed :math:`\Omega_\text{ci}` (region 1) the torque is zero. Between cut-in and the rated region-2 speed :math:`\Omega_{2}`, the torque rises linearly (region 1½). In region 2 the torque follows the optimal power curve

.. math::
    :label: eq:genCtrl

    Q_\text{gen} = k_P\,\Omega_\text{gen}^2, \qquad
    k_P = \frac{\pi\,\rho\,R^5\,C_{p,\text{opt}}}{2\,\lambda_\text{opt}^3\,i_\text{gbx}^3},

where :math:`\lambda_\text{opt}` is the optimal tip-speed ratio. Between the top of region 2 and rated speed (region 2½) the torque is interpolated linearly, and above rated (region 3) it is held at :math:`Q_\text{rated}`. A torque rate limiter constrains :math:`|dQ_\text{gen}/dt| \leq \dot{Q}_\text{max}` when enabled. Alternatively, the rotor speed can be assigned directly from a lookup table of :math:`\Omega(U_\text{ref})` by setting ``genControllerType = rpmControlCurve``.

**Blade pitch controller.**
In region 3 the blade pitch is regulated by a gain-scheduled PID controller. The scheduling gain is

.. math::
    :label: eq:pitchGain

    G(\beta) = \frac{1}{1 + \beta / \beta_{s2r}},

where :math:`\beta` is the current pitch angle and :math:`\beta_{s2r}` is the pitch angle at which the sensitivity halves (``pitchGainSchedule``). The error signal is :math:`e = \Omega - \Omega_\text{rated}`, and the commanded pitch is

.. math::
    :label: eq:pitchPID

    \beta_\text{cmd} = G\!\left(K_P\,e + K_I\int_0^t e\,d\tau + K_D\,\dot{e}\right),

saturated to :math:`[\beta_\text{min}, \beta_\text{max}]`. The integral is clamped to prevent wind-up. Setting ``pitchControllerType = bladePitchCurve`` replaces the PID with a static :math:`\beta(U_\text{ref})` lookup table.

**Nacelle yaw controller.**
The nacelle yaw controller estimates the mean wind direction from a velocity sample taken 2R upstream of the rotor centre (``yawSamplingType = hubUpDistance``). An exponential moving average filters the instantaneous flow angle :math:`\phi_0 = \text{atan2}(u_y, u_x)`:

.. math::
    :label: eq:yawFilter

    \phi(t) = \phi(t-\Delta t) + \frac{\Delta t}{\tau_\text{yaw}}\!\left[\phi_0(t) - \phi(t-\Delta t)\right],

where :math:`\tau_\text{yaw}` is the yaw time constant (``yawRelaxFactor``). When the yaw error :math:`|\phi - \psi_\text{yaw}|` exceeds the deadband ``yawAllowedError``, the nacelle rotates at the prescribed yaw rate ``yawSpeed`` (degrees/s) within the angular bounds ``yawMin``/``yawMax``. All rotor geometry—rotor axis, rotor direction, rotor centre, and all actuator points—is rotated accordingly. Yaw data are synchronised to the master rank only at output write times, requiring a single global MPI reduction per time series write.

**Wind farm controller.**
An external wind farm controller applies time-varying pitch or thrust offsets to all turbines simultaneously. A plain-text time-series file supplies pitch offsets :math:`\Delta\beta_\text{wf}(t)` for ADM/ALM turbines, or thrust-coefficient offsets :math:`\Delta C_t(t)` for uniformADM/AFM turbines. The values are linearly interpolated in time and added to the individual turbine setpoints at each step.


Upstream Velocity Sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~

All turbine models that require a reference velocity (uniformADM rotorDisk sampling, yaw control, Uref monitoring) use a rotor-shaped upsampling grid at 2.5 diameters upstream of the rotor centre. The grid consists of a :math:`12 \times 12` azimuth–radial array of points spanning from hub to tip radius. The area weight of each point accounts for the annular area of its radial band, including a redistribution of the hub area over the swept annulus. The sampled velocity components are area-averaged to yield the scalar reference velocity :math:`U_\text{ref}`.


Parallel Implementation
~~~~~~~~~~~~~~~~~~~~~~~~

To minimise global communication, each turbine is managed by a dedicated sub-communicator ``TRB_COMM`` of at most 8 MPI ranks. Only the mesh cells within the turbine's influence sphere—the smallest bounding box that contains the rotor at any yaw angle—are mapped to a turbine. Force projection, velocity sampling, and checkpoint I/O are all confined to these rank subsets.

The upsampling grid has its own communicator ``UPW_COMM``. When the nacelle yaws, ``trbMoved`` is set and the controlled-cell search is repeated. A single global ``MPI_Reduce`` is issued only at output write times to synchronise the yaw angle across the full domain. Each turbine designates a writer rank (rank 0 of ``TRB_COMM``) that handles all checkpoint files and time-series logs.


OpenFAST Coupling
~~~~~~~~~~~~~~~~~

When compiled with ``-DUSE_OPENFAST``, any turbine in the farm can be delegated to `OpenFAST <https://openfast.readthedocs.io>`_ for its structural dynamics, detailed aerodynamics (AeroDyn), and generator/pitch control (ServoDyn). The remaining turbines continue to use TOSCA-native models. The two solvers exchange data every CFD time step through a partitioned coupling scheme.

**Initialisation.**
TOSCA creates one ``fast::OpenFAST`` instance per coupled turbine and populates the coupling metadata structure ``fi``:

- ``dtFAST = dtCFD / nFastSubSteps`` — OpenFAST internal time step;
- FST input file: ``turbines/<type>.<id>.fst``;
- restart checkpoint: ``turbines/<type>.<id>.chkpt``.

OpenFAST exposes two independent grids of coupling points per turbine: *velocity points* (AeroDyn inflow nodes, classified as hub=0, blade=1, or tower=2) and *force points* (AeroDyn/ElastoDyn structural nodes). All coordinates are maintained in the TOSCA global Cartesian frame. After ``FAST->init()`` the number of force and velocity points is confirmed and storage for ``velPtsBlade``, ``velValsBlade``, ``forcePtsBlade``, ``forceValsBlade``, and the tower equivalents is allocated.

**Per-step protocol (one CFD time step).**
The coupling follows a sequential predictor–corrector protocol:

1. **Point discovery.** ``findControlledPointsRotorOpenFAST`` locates the nearest CFD cell to every velocity and force point. Ties between processors are broken by a small rank-dependent perturbation so that each point is owned by exactly one rank.

2. **Velocity transfer (TOSCA → OpenFAST).** ``computeWindVectorsRotorOpenFAST`` and ``computeWindVectorsTowerOpenFAST`` trilinearly interpolate the TOSCA velocity field at every OpenFAST velocity point, then call ``FAST->setVelocity`` on the writer rank.

3. **OpenFAST advance.** ``stepOpenFAST`` calls ``FAST->step()`` *nFastSubSteps* times per CFD time step, advancing the OpenFAST internal state at the finer resolution ``dtFAST``.

4. **Geometry update.** ``getForcePtsBladeOpenFAST`` and ``getForcePtsTwrOpenFAST`` retrieve the updated structural positions from OpenFAST and broadcast them to all ranks. ``mapOFDisplToActPts`` then copies these positions to the TOSCA actuator-point arrays (e.g. ``alm.points``).

5. **Force transfer (OpenFAST → TOSCA).** ``mapOFForcesToActPts`` reads the aerodynamic body forces from the OpenFAST force points and writes them to the TOSCA actuator-point force arrays, scaled by :math:`1/\rho`. These forces are subsequently projected onto the CFD mesh by the standard Gaussian projector.

**Global rotor state.**
``getGlobParamsOpenFAST`` extracts all rotor-level scalars needed by the TOSCA I/O and upsampling infrastructure:

- *Hub position* and *shaft direction*: obtained via ``FAST->getHubPos`` / ``FAST->getHubShftDir`` and broadcast to all ranks.
- *Rotor speed* :math:`\Omega`: derived from the angular displacement of the last blade force point (blade-3 tip) between successive steps, :math:`\Omega = \Delta\theta / \Delta t_\text{CFD}`, where :math:`\Delta\theta = \arccos(|\hat{\mathbf{r}}_\text{new}\!\cdot\!\hat{\mathbf{r}}_\text{old}|)`.
- *Azimuth*: computed from the projection of the tip vector onto the rotor plane; 0° when blade 0 points vertically up, increasing in the clockwise direction when viewed from the front.
- *Thrust and torque*: obtained from ``FAST->computeTorqueThrust`` and projected onto the rotor axis; also independently accumulated from the force-point array for cross-validation.
- *Rotation direction*: inferred from :math:`\operatorname{sign}\!\left((\hat{\mathbf{r}}_\text{old}\times\hat{\mathbf{r}}_\text{new})\cdot\hat{\boldsymbol{\Omega}}\right)`.

**Yaw handling.**
Because OpenFAST (ServoDyn) controls nacelle yaw, TOSCA reads the current yaw angle from the updated shaft direction rather than running its own yaw controller. Whenever the yaw angle changes by more than 1°, the upsampling rig is fully rebuilt: the 12 × 12 azimuth–radial grid is recomputed at 2.5D upstream in the new rotor-normal direction, and ``trbMoved`` is set to trigger a new controlled-cell search.
