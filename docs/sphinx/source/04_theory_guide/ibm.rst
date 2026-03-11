.. _ibm-theory-section:

Immersed Boundary Method
========================

TOSCA implements a **ghost-cell, direct-forcing** immersed boundary method (IBM)
that allows complex geometries — terrain features, wind turbine towers and
nacelles, and moving objects — to be represented on a structured background grid
without mesh deformation.  The immersed surface enters the momentum equation
through a body-force term :math:`\mathbf{f}^{IB}`, which is applied implicitly:
at each time step the velocity (and temperature) field is prescribed in *IBM
interface cells* (ghost cells) to values consistent with the desired wall
boundary condition, so the forcing is never evaluated as an explicit source
term.

.. figure:: ../../images/ibm-sketch.jpg
   :width: 70%
   :align: center

   Schematic of the ghost-cell IBM.  The immersed surface (triangular mesh,
   black) divides the background grid into fluid cells (white), solid cells
   (dark), and IBM interface cells (grey).  Points **C** and **D** are the two
   wall-normal sampling locations used during boundary-condition
   reconstruction when a wall model is active.


Surface Representation
----------------------

The immersed-body surface is described as an **unstructured triangular mesh**
provided in STL (``.stl``) or Abaqus (``.inp``) format and placed in the
``IBM/`` directory.  Multiple bodies can be defined simultaneously in
``IBM/IBMProperties.dat`` as ``object0``, ``object1``, etc.

To accelerate element searches the code builds a two-level spatial data
structure during initialisation:

1. **Axis-aligned bounding box (AABB)** — a tight box around each body
   (``findBodyBoundingBox``), used to exclude bodies entirely when a query
   point lies far outside.

2. **Search-cell list** — the bounding box is subdivided into uniform
   search cells whose size is proportional to the mean element edge length
   (``findSearchCellDim``); every triangular element is inserted into all
   search cells it overlaps (``createSearchCellList``).  Only elements in
   the relevant search cells are tested during subsequent point-in-triangle
   queries, reducing the search cost from :math:`O(N_e)` to approximately
   :math:`O(1)`.


Cell Classification (Topology Check)
-------------------------------------

The function ``ibmSearch`` classifies every background-grid cell into one of
three categories by **ray-casting**: a ray is cast from each cell centre and the
number of intersections with the triangular IBM surface is counted.

* **Fluid cell** — even number of intersections: the cell lies outside the body.
* **Solid cell** — odd number of intersections: the cell lies inside the body.
* **IBM interface cell (ghost cell)** — any solid cell immediately adjacent to
  at least one fluid cell.  These cells are not advanced by the flow solver;
  instead their field values are overwritten each time step by the
  IBM boundary-condition reconstruction.

For each IBM interface cell the code identifies the closest triangular element
and stores its outward unit normal :math:`\hat{\mathbf{n}}` and the projection
of the cell centre onto the surface (the *wall point* :math:`\mathbf{p}_\text{w}`).
The closest-element search uses a bounding-sphere acceleration: each element's
circumsphere centre :math:`\mathbf{q}` and radius :math:`r` are pre-computed,
and the element with the smallest value of
:math:`\|\mathbf{x} - \mathbf{q}\| - r` is selected as the candidate before
the precise face/edge/vertex projection is performed.


Boundary-Condition Reconstruction
----------------------------------

Two reconstruction modes are available, selected via ``wallShearOn`` (and the
corresponding ``velocityBC``) in ``IBMProperties.dat``.

Without Wall Model
~~~~~~~~~~~~~~~~~~

When ``wallShearOn = false`` and ``velocityBC`` is ``slip`` or ``noSlip``,
the ghost-cell velocity is set by **one-point interpolation**.  A single
background point **B** is placed a distance :math:`\delta` along
:math:`\hat{\mathbf{n}}` from the ghost-cell centre; :math:`\delta` equals the
minimum local cell size (``CurvibTrilinear``) or the distance to an interception
point on an adjacent grid plane (``CurvibTriangular``).  Denoting the
wall-normal signed distances from the IBM surface as :math:`s_b` (ghost cell)
and :math:`s_c` (background point **B**), the reconstruction gives

.. math::
   \mathbf{u}^{(b)} = \frac{s_b}{s_c}\,\mathbf{u}_B +
                       \left(1 - \frac{s_b}{s_c}\right)\mathbf{u}_\text{wall},

where :math:`\mathbf{u}_\text{wall}` is the IBM-surface velocity interpolated
from the triangular mesh nodes using barycentric weights.  For a ``slip``
boundary the normal component of the ghost-cell velocity is set to zero and the
tangential component is reflected from **B**.

With Wall Model (Ghost-Cell Hybrid Approach)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``wallShearOn = true`` and ``velocityBC = velocityWallFunction``, the
reconstruction uses **two background sampling points** placed at distances
:math:`s_c` and :math:`s_d` from the IBM wall.  Setting
:math:`\delta = \texttt{interpDist}`:

.. math::
   s_b = -\varepsilon, \qquad
   s_c = \delta, \qquad
   s_d = 3\,\delta,

where :math:`\varepsilon` is the distance from the ghost cell centre to the
projected wall point.  The arrangement along the wall normal is therefore

.. math::
   \underbrace{b}_{\text{ghost}} \;\cdots\; \text{wall} \;\cdots\;
   \underbrace{c}_{\mathbf{C}} \;\cdots\; \underbrace{d}_{\mathbf{D}}.

Velocities at **C** (:math:`\mathbf{u}_C`) and **D** (:math:`\mathbf{u}_D`) are
obtained by trilinear interpolation from the background grid.  The ghost-cell
velocity is then set by the selected wall-shear model
(``CurvibInterpolationInternalCell``) using all three points
:math:`(b, c, d)` and the IBM-surface velocity :math:`\mathbf{u}_\text{wall}`.

Pressure in the ghost cell is set by linear wall-normal extrapolation from
**C**, corrected for the IBM-surface acceleration:

.. math::
   p^{(b)} = p^{(c)} -
              \left(\frac{\partial \mathbf{u}_\text{wall}}{\partial t}
              \cdot \hat{\mathbf{n}}\right)(s_c - s_b).


Wall Models
-----------

The wall boundary condition for velocity and temperature is configured per-body
in ``IBMProperties.dat`` (or matched from a domain boundary face via
``velocityBCSetType = matchU*``).

Velocity Wall Models
~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 12 22 66

   * - Type
     - Name
     - Description
   * - ``-1``
     - Cabot
     - Neutral log-law with optional aerodynamic roughness ``kRough``.
       The friction velocity :math:`u_*` is solved iteratively from
       :math:`u_t = (u_*/\kappa)\ln[(s_c + z_0)/z_0]`, where :math:`\kappa`
       is the von Kármán constant and :math:`z_0 = \texttt{kRough}`.  The
       smooth-wall limit is recovered when ``kRough = 0``.
   * - ``-3``
     - Schumann
     - Stability-corrected log-law.  The relation
       :math:`u_t = (u_*/\kappa)[\ln(s_c/z_0) - \psi_m(s_c/L)]`
       is used, where :math:`L` is the Obukhov length and :math:`\psi_m` is
       the Monin–Obukhov stability function.  Must be paired with the
       Schumann temperature wall model (``-2`` or ``-4``) when temperature
       is active.
   * - ``-4``
     - Power-law APG
     - Power-law wall profile for adverse-pressure-gradient flows.
   * - ``-5``
     - Log-law APG
     - Log-law wall profile for adverse-pressure-gradient flows.

Temperature Wall Models
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 12 28 60

   * - Type
     - Name
     - Description
   * - ``-2``
     - Shumann — constant flux
     - Prescribed constant wall heat flux ``qWall``
       :math:`[\mathrm{K\,m\,s^{-1}}]`.
   * - ``-3``
     - Shumann — heating rate
     - Prescribed constant surface heating rate.
       *Not yet implemented for IBM.*
   * - ``-4``
     - Shumann — time history
     - Surface potential temperature and Obukhov length read from a
       time-series file (``readSurfaceTempData``).  The ghost-cell temperature
       is interpolated linearly between the IBM surface temperature (looked
       up at the current time) and the value at point **C**.

In addition to ``thetaWallFunction``, the temperature boundary condition at an
IBM surface can be set to ``zeroGradient`` (ghost-cell temperature equals the
value at **C**) or ``fixedValue`` (ghost-cell temperature equals the prescribed
``fixedValueT``).


Interpolation Schemes
---------------------

The interpolation model is selected via ``IBInterpolationModel`` in
``IBMProperties.dat``.  The available options and their governing functions are
listed below.

.. list-table::
   :header-rows: 1
   :widths: 22 22 16 40

   * - ``IBInterpolationModel``
     - ``curvibType``
     - ``curvibOrder``
     - Function / Notes
   * - ``CURVIB``
     - ``CurvibTrilinear``
     - ``linear``
     - ``CurvibInterpolation`` — one background point at one minimum cell
       size along :math:`\hat{\mathbf{n}}`; supports Cabot, Schumann, slip,
       noSlip boundary conditions.
   * - ``CURVIB``
     - ``CurvibTrilinear``
     - ``quadratic``
     - ``CurvibInterpolationQuadratic`` — as above but uses a quadratic
       wall-normal profile.
   * - ``CURVIB``
     - ``CurvibTriangular``
     - —
     - ``CurvibInterpolationTriangular`` — the image point is the
       intersection of the wall-normal ray with an adjacent curvilinear grid
       plane; barycentric weights of the enclosing triangle are used for
       interpolation.  Supports the same boundary conditions as trilinear.
   * - ``CURVIB`` + ``wallShearOn``
     - —
     - —
     - ``CurvibInterpolationInternalCell`` — two background points; full
       hybrid wall-shear reconstruction with quadratic pressure extrapolation
       (see :ref:`ibm-theory-section`).
   * - ``MLS``
     - —
     - —
     - ``MLSInterpolation`` — Moving Least Squares reconstruction over a
       cloud of fluid support nodes identified by ``findFluidSupportNodes``
       and ``findIBMMeshSupportNodes``.

.. note::
   The ``wallShearOn`` flag requires ``IBInterpolationModel = CURVIB`` and
   a ``velocityBC = velocityWallFunction`` entry; it overrides the
   ``curvibType`` / ``curvibOrder`` settings and always uses
   ``CurvibInterpolationInternalCell``.


Dynamic IBM
-----------

Setting ``dynamic = true`` in ``IBMProperties.dat`` enables a moving immersed
body.  At every time step ``UpdateIBM`` rebuilds the full cell-classification
data structures:

1. Previous lists are destroyed (``destroyLists``).
2. The mesh is advanced by ``UpdateIBMesh`` (see motion types below).
3. The bounding box, search-cell list, outward normals, ray-casting
   classification, IBM fluid-cell list, and closest-element search are all
   recomputed from scratch.

This per-step re-initialisation supports arbitrarily large motions with no
restriction on the displacement per time step.

The mesh motion is performed by ``UpdateIBMesh``, which dispatches to the
following ``bodyMotion`` types:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - ``bodyMotion``
     - Description
   * - ``static``
     - No motion; lists are built once at initialisation and not rebuilt.
   * - ``rotation``
     - Rigid rotation about ``rotCenter`` and ``rotAxis`` at angular speed
       ``angSpeed`` :math:`[\mathrm{rad\,s^{-1}}]`.  The speed may ramp at
       a constant rate ``angAcc`` :math:`[\mathrm{rad\,s^{-2}}]`.  Node
       coordinates and velocities are updated by ``rotateIBMesh``.  The
       instantaneous angular position is written to
       ``fields/<mesh>/ibm/<time>/<bodyName>``.
   * - ``sinusoidal``
     - Sinusoidal translation along ``motionDir`` with ``amplitude``
       :math:`A\,[\mathrm{m}]` and ``frequency`` :math:`f\,[\mathrm{Hz}]`.
       The displacement applied each step is
       :math:`A\bigl(\cos 2\pi f t_\text{prev} - \cos 2\pi f t\bigr)`
       (function ``sineMotion``).
   * - ``pitchingOscillation``
     - Sinusoidal pitching rotation about ``pitchCenter`` / ``pitchAxis``
       with amplitude :math:`A\,[\mathrm{rad}]` and frequency :math:`f`
       (function ``pitchingMotion``).


Force and Moment Integration
-----------------------------

When ``computeForce = true``, the aerodynamic load on each IBM object is
evaluated every ``timeInterval`` (or physical-time interval, depending on
``intervalType``) by ``ComputeForceMoment``.  For each triangular element of
the IBM surface mesh the procedure is:

1. A pressure-interpolation point is placed at a distance
   :math:`\ell = \sqrt{2\,A_e}` (where :math:`A_e` is the element area)
   outward from the element centre along :math:`\hat{\mathbf{n}}`.  If that
   point lies inside solid cells it is shifted outward in increments of
   :math:`0.2\,\ell` until a valid trilinear interpolation stencil is found.

2. The static pressure :math:`p` and velocity :math:`\mathbf{u}` are
   trilinearly interpolated at the interpolation point.

3. The **pressure force** on the element is

   .. math::
      \mathbf{F}_P^{(e)} = -\rho\!\left(p - \frac{\partial u_{\text{wall},n}}
      {\partial t}\,\ell\right)A_e\,\hat{\mathbf{n}},

   where :math:`u_{\text{wall},n} = \partial\mathbf{u}_\text{wall}/\partial t
   \cdot\hat{\mathbf{n}}` is the IBM-surface normal acceleration, which
   provides a correction for the unsteady pressure gradient next to a moving
   wall.

4. The **viscous shear force** is computed from the friction velocity
   :math:`u_* = u_\tau(\mathbf{u}_B - \mathbf{u}_\text{wall})` using the
   Cabot model and applied in the direction of the relative tangential velocity.

Net pressure force, viscous force, aerodynamic torque, and shaft power (for
``rotation`` bodies) are reduced across all MPI ranks and appended to

.. code-block:: none

   postProcessing/<meshName>/IBM/<startTime>/<bodyName>/netForce/<bodyName>_netForce

Per-element pressure loads (in Abaqus ``.dlo`` format) are written to the
``elementForce/`` subfolder when ``writePForce = true`` in the ``io`` sub-dict.
