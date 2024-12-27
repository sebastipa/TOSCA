.. _ablProperties-section:

`ABLProperties.dat`
~~~~~~~~~~~~~~~~~~~

When the ``-abl`` flag is activated in the ``control.dat`` file, TOSCA is prompted to read the ``ABLProperties.dat`` file.
This contains important settings that can be activated for ABL or wind farm simulations with potential temperature
stratification. For example, this file contains settings for the potential temperature profile, it allows to activate
velocity controllers, used to drive the flow when laterally-periodic simulations are applied, the potential temperature
controller, used to fix the potential temperature profile over time, fringe and damping regions, as well as the geostrophic damper,
which avoids inertial oscillations produced when the flow is not initialized in geostrophic balance and the Coriolis force is
active. It is also in this file that one can activate TOSCA's canopy model or the concurrent precursor method.

Primarily, this file contains entries that are required by TOSCA when the ``-abl`` flag is activated. In
addition, it contains dictionaries that are required when additional flags are activated in the ``control.dat`` file. For
example, the ``xDampingProperties`` dictionary is only read if ``-xDampingLayer`` is set to 1 in the ``control.dat`` file.

The folliwing table summarizes all entries and dictionaries available in the ``ABLProperties.dat``. A description of
entries specific to each dictionaries is given in the subsequent table.

.. table::
   :widths: 32, 15, 53
   :align: center

   =================================== ============== ==================================================
   **entry name**                      **entry type**    **description**
   ----------------------------------- -------------- --------------------------------------------------
   ``hRough``                          scalar         equivalent roughness length in m used to set
                                                      set the initial velocity condition when
                                                      ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/U`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``uRef``                            scalar         reference velocity in m/s used to set
                                                      set the initial condition when
                                                      ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/U`` file. The initial wind will be
                                                      a logarthmic profile with this value of velocity
                                                      at ``hRef``. It is also used by the
                                                      velocity controller when ``controllerType`` is set
                                                      to *pressure* inside the ``controllerProperties``
                                                      dictionary.
   ----------------------------------- -------------- --------------------------------------------------
   ``hRef``                            scalar         reference height in m used to set
                                                      set the initial condition when
                                                      ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/U`` file. The initial wind will be
                                                      a logarthmic profile with velocity ``uRef``
                                                      at this height. It is also used by the
                                                      velocity controller when ``controllerType`` is set
                                                      to *pressure* inside the ``controllerProperties``
                                                      dictionary.
   ----------------------------------- -------------- --------------------------------------------------
   ``hInv``                            scalar         inversion layer height in m used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``dInv``                            scalar         inversion layer wifth in m used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``gInv``                            scalar         inversion layer strength in K used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``gTop``                            scalar         lapse rate above inversion in K/m used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``gABL``                            scalar         lapse rate below inversion in K/m used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``tRef``                            scalar         reference potential temperaure used to set
                                                      set the initial potential temperature condition
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file. This is also assumed as a
                                                      reference temperature throughout the simulation.
   ----------------------------------- -------------- --------------------------------------------------
   ``vkConst``                         scalar         von Karman constant.
   ----------------------------------- -------------- --------------------------------------------------
   ``smearT``                          scalar         smearing parameter used to define how sharp the
                                                      initial potential temperaure profile is, when
                                                      when ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/T`` file. Suggested value is 0.33.
   ----------------------------------- -------------- --------------------------------------------------
   ``coriolisActive``                  bool           Coriolis force activation flag.
   ----------------------------------- -------------- --------------------------------------------------
   ``fCoriolis``                       scalar         Coriolis parameter. Should be computed as
                                                      :math:`10^{-5} \cdot 7.272205217sin(\phi)`,
                                                      where :math:`\phi` is latitude.
   ----------------------------------- -------------- --------------------------------------------------
   ``controllerActive``                bool           Velocity controller activation flag. Requires
                                                      ``controllerProperties`` dictionary.
   ----------------------------------- -------------- --------------------------------------------------
   ``controllerActiveT``               bool           Temperature controller activation flag. The target
                                                      laterally-averaged potential temperature profile
                                                      that the controller aims at maintaining is the
                                                      initial one. Requires
                                                      ``controllerProperties`` dictionary.
   ----------------------------------- -------------- --------------------------------------------------
   ``controllerActivePrecursorT``      bool           Temperature controller activation flag. Same as
                                                      ``controllerActiveT``, but for the concurrent
                                                      precursor simulation. Requires
                                                      ``controllerProperties`` dictionary.
   ----------------------------------- -------------- --------------------------------------------------
   ``perturbations``                   bool           Adds divergence-free sinusoidal perturbations to
                                                      triger turbulence in the initial condition when
                                                      ``internalField`` is set to *ABLFlow* in the
                                                      ``boundary/U`` file.
   ----------------------------------- -------------- --------------------------------------------------
   ``controllerProperties``            dictionary     Contains inputs for velocity and temperature
                                                      controllers. Required when
                                                      ``controllerActive``, ``controllerActiveT``
                                                      or ``controllerActivePrecursorT`` are set to 1.

                                                      Usage:

                                                      .. code-block:: C

                                                         controllerProperties
                                                         {
                                                            relaxPI                scalar
                                                            controllerMaxHeight    scalar
                                                            controllerType         string
                                                            alphaPI                scalar
                                                            timeWindowPI           scalar
                                                            geostrophicDamping     bool
                                                            geoDampingAlpha        scalar
                                                            geoDampingStartTime    scalar
                                                            geoDampingTimeWindow   scalar
                                                            hGeo                   scalar
                                                            alphaGeo               scalar
                                                            uGeoMag                scalar
                                                            controllerAvgStartTime scalar
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``xDampingProperties``              dictionary     Defines fringe region parameters, activated with
                                                      ``-xDampingLayer`` 1 in ``control.dat``.

                                                      Usage:

                                                      .. code-block:: C

                                                         xDampingProperties
                                                         {
                                                            xDampingStart             scalar
                                                            xDampingEnd               scalar
                                                            xDampingDelta             scalar
                                                            xDampingAlpha             scalar
                                                            xDampingAlphaControlType  scalar
                                                            xDampingLineSamplingYmin  scalar
                                                            xDampingLineSamplingYmax  scalar
                                                            xDampingTimeWindow        scalar
                                                            uBarSelectionType         integer
                                                            // additional parameters depending
                                                            // on uBarSelectionType (see next
                                                            // table)
                                                         }

                                                      The *uBarSelectionType* entry defines how the
                                                      reference wind field is calculated inside the
                                                      fringe region, and it requires additional
                                                      parameters depending on the type. The **concurrent
                                                      precursor** (i.e. when this reference field is
                                                      solved concurrently with the main simulation) is
                                                      activated by setting the *uBarSelectionType* to
                                                      3. TOSCA creates the second simulation instance
                                                      automatically, without requiring additional user
                                                      parameters.
   ----------------------------------- -------------- --------------------------------------------------
   ``yDampingProperties``              dictionary     Defines lateral fringe region parameters,
                                                      activated with  ``-yDampingLayer`` 1 in
                                                      ``control.dat``.

                                                      Usage:

                                                      .. code-block:: C

                                                         yDampingProperties
                                                         {
                                                            yDampingStart  scalar
                                                            yDampingEnd    scalar
                                                            yDampingDelta  scalar
                                                            yDampingAlpha  scalar
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``zDampingProperties``              dictionary     Defines Reyleigh damping layer parameters,
                                                      activated with  ``-zDampingLayer`` 1 in
                                                      ``control.dat``.

                                                      Usage:

                                                      .. code-block:: C

                                                         zDampingProperties
                                                         {
                                                            zDampingStart   scalar
                                                            zDampingEnd     scalar
                                                            zDampingAlpha   scalar
                                                            zDampingAlsoXY  bool
                                                            zDampingXYType  integer
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``advectionDampingProperties``      dictionary     Defines advection damping regions parameters.
                                                      This corresponds to the technique developed by
                                                      Lanzilao and Meyers (2022a). It is activated with
                                                      ``-advectionDamping`` 1 in ``control.dat``.

                                                      Usage:

                                                      .. code-block:: C

                                                         advectionDampingProperties
                                                         {
                                                            advDampingStart       scalar
                                                            advDampingEnd         scalar
                                                            advDampingDeltaStart  scalar
                                                            advDampingDeltaEnd    scalar
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``kLeftDampingProperties``          dictionary     Defines Rayleigh damping layer at the *kLeft*
                                                      patch. Requires ``-kLeftRayleigh`` set to 1 in
                                                      ``control.dat``. Damping transitions from zero
                                                      to max across a layer of width
                                                      *kLeftFilterWidth* centered at
                                                      *kLeftFilterHeight* and is applied between the
                                                      *kLeft* patch and a plane at a distance
                                                      *kLeftPatchDist* from the *kLeft* patch to
                                                      obtain the desired velocity *kLeftDampingUBar*.

                                                      Usage:

                                                      .. code-block:: C

                                                         kLeftDampingProperties
                                                         {
                                                            kLeftPatchDist     scalar
                                                            kLeftDampingAlpha  scalar
                                                            kLeftDampingUBar   vector
                                                            kLeftFilterHeight  scalar
                                                            kLeftFilterWidth   scalar
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``kRightDampingProperties``         dictionary     Defines Rayleigh damping layer at the *kRight*
                                                      patch. Requires ``-kRightRayleigh`` set to 1 in
                                                      ``control.dat``. Damping transitions from zero
                                                      to max across a layer of width
                                                      *kRightFilterWidth* centered at
                                                      *kRightFilterHeight* and is applied between the
                                                      *kRight* patch and a plane at a distance
                                                      *kRightPatchDist* from the *kRight* patch to
                                                      obtain the desired velocity *kRightDampingUBar*.

                                                      Usage:

                                                      .. code-block:: C

                                                         kRightDampingProperties
                                                         {
                                                            kRightPatchDist     scalar
                                                            kRightDampingAlpha  scalar
                                                            kRightDampingUBar   vector
                                                            kRightFilterHeight  scalar
                                                            kRightFilterWidth   scalar
                                                         }

   ----------------------------------- -------------- --------------------------------------------------
   ``canopyProperties``                dictionary     Defines input parameters for the canopy model.
                                                      Requires ``-canopy`` set to 1 in ``control.dat``
                                                      file.

                                                      Usage:

                                                      .. code-block:: C

                                                         canopyProperties
                                                         {
                                                            xStartCanopy     scalar
                                                            xEndCanopy       scalar
                                                            yStartCanopy     scalar
                                                            yEndCanopy       scalar
                                                            zStartCanopy     scalar
                                                            zEndCanopy       scalar
                                                            cftCanopy        scalar
                                                            diskDirCanopy    vector
                                                         }

   =================================== ============== ==================================================

The meaning of the entires required in the dictionaries listed in the above table are described below.

.. table::
   :widths: 35, 20, 45
   :align: center

   ============================= ================== =====================================================================================
   **entry**                     **entry type**     **description**
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *controllerProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``relaxPI``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``controllerMaxHeight``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``controllerType``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``alphaPI``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``timeWindowPI``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``geostrophicDamping``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``geoDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``geoDampingStartTime``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``geoDampingTimeWindow``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``hGeo``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``alphaGeo``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``uGeoMag``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``controllerAvgStartTime``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *xDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``xDampingStart``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingEnd``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingDelta``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingAlphaControlType``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingLineSamplingYmin``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingLineSamplingYmax``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xDampingTimeWindow``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``uBarSelectionType``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *yDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``yDampingStart``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``yDampingEnd``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``yDampingDelta``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``yDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *zDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``zDampingStart``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zDampingEnd``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zDampingAlsoXY``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zDampingXYType``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *advectionDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``advDampingStart``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``advDampingEnd``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``advDampingDeltaStart``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``advDampingDeltaEnd``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *kLeftDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``kLeftPatchDist``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kLeftDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kLeftDampingUBar``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kLeftFilterHeight``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kLeftFilterWidth``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *kRightDampingProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``kRightPatchDist``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kRightDampingAlpha``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kRightDampingUBar``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kRightFilterHeight``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``kRightFilterWidth``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   *canopyProperties*
   --------------------------------------------------------------------------------------------------------------------------------------
   ``xStartCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``xEndCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``yStartCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``yEndCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zStartCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``zEndCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``cftCanopy``
   ----------------------------- ------------------ -------------------------------------------------------------------------------------
   ``diskDirCanopy``
   ============================= ================== =====================================================================================
