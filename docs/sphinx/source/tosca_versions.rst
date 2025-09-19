Versions 
========

Development
-----------

The **master** branch is the latest stable development branch of TOSCA. New working features that are tested to a certain extent are merged into this branch.
For ongoing development and testing, please refer to the respective branches based on the specific feature of interest. Cloning the TOSCA repository will 
automatically set the **master** as the default branch. If you want to use previous releases, please refer to the `respective tags <https://github.com/sebastipa/TOSCA/releases>`_.

The documentation for the **master** branch, which is updated alongside the introduction of new features and/or code refactoring, is available at 
`TOSCA Stable Development Documentation <https://sebastipa.github.io/TOSCA/index.html>`_. We save snapshots of the documentation for each release, which can 
be found below.     

New working features since release (V1.1.0):  
 - Added integral damping for `geostrophic` ABL controller. 
 - Added geostrophic profile assimilation method for micro- to meso-scale coupling.
 - Code refactoring: added lesScalar.c to handle eddy diffusivity SGS models.
 - Added computation of 5 invariants of the Reynolds stress tensor.
 - Generalized velocity wall models for non-cartesian grids.
 - Added Arnold–Beltrami–Childress (ABC) flow and homogeneous isotropic turbulence (HIT) initial conditions. 
 - Added Bardina LES model. 
 - Added 4th-order gradient scheme. 
 - Removed IBM half-edge structure when using IBM bodies that are not fully closed. 
 - Added variable turbine Ct based on a Ct vs Uref lookup table. (only for AFM and Uniform ADM). 
 - Added additional fields as output when using geostrophic profile assimilation. 
 - Added possibility to read temperature sources when using `timeHeightSeries` in temperature controller. 
 - Updated concurrent-precusror method to include latest meso- to micro-scale coupling implementation. 
 - Changed NREL5MW test case to use variable Ct based on lookup table. 
 - Added additional test case ChannelFlowRe395Test (channel flow at Re = 395).
 - Added control action `read` or `write` also for temperature controller. 
 - Updated and extended documentation to incorporate new changes & features.

Bug fixes: 
- Fields are only written within the `writeFields()` function now in order to avoid problem with temporary file names while writing is in progress.


V1.1.0
------

This is release 1.1.0 of TOSCA. An image of the original documentation is 
available at `TOSCA v1.1.0 Documentation <https://sebastipa.github.io/TOSCA/versions/1.1.0/index.html>`_

New working features since release (V1.0.0):  
 - Two-way coupled overset mesh method: this feature enables complete block-structed mesh refinement in TOSCA. 
 - Additional scale-dependent (``dynamicLASD``) and anisotropic minimum dissipation (``amd``) LES models.  
 - Additional high-order numerical scheme (``-central4`` with ``-hyperVisc`` parameter) to better capture the high-frequency tail of turbulence spectra. 
 - Improved immersed boundary method: faster algorithms and time-varying wall modeling. 
 - Lateral fringe region (``-yDampingLayer``): allows for gravity wave problems with wind direction changes and assimilated atmospheric data. 
 - Data assimilation methods: direct/indirect profile assimilation techniques, as well as wavelet assimilation method (selectable through ``controllerType`` and ``controllerTypeT`` flags in *ABLProperties.dat*). 
 - Additional example cases: *IndirectProfileAssimilationTest*, *WaveletAssimilationTest*, *MultiDomainOverset*, *PitchingAerofoilOversetTest*, *SuccessorOversetTest*, *SuccessorPeriodizedLateralFringeTest* and *TaylorGreenVortexTest*.
 - Updated and extended documentation to incorporate new changes & features. 

Bug fixes: 
 - Corrected bug in the fixedGradient BC for temperature. 
 - Corrected IBM shear stress-based wall model search algorithm.
 - Corrected bug when writing fields using the `-perturbABL` flag.
 - Corrected bug in IBM interpolation that was cousing oscillations in the pressure field. 
 - Corrected bug when creating the timelist for `inletFunctions` 3 and 4, now any possible duplicate time is removed. 

Improvements:
 - Improved `-purgeWrite` so that more than one time step can be kept. 
 - Improved overset initialization algorithm using octree search method.
 - Improved theory guide. 


V1.0.0
------

This is the first official release of TOSCA. An image of the original documentation is 
available at `TOSCA v1.0.0 Documentation <https://sebastipa.github.io/TOSCA/versions/1.0.0/index.html>`_