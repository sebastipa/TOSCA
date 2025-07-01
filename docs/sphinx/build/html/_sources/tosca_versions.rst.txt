Versions 
========

Development
-----------

The **master** branch is the current stable development branch of TOSCA. Dodumentation is available at 
`TOSCA Stable Development Documentation <https://sebastipa.github.io/TOSCA/index.html>`_

New working features since last release (V1.0.0):  

- Two-way coupled overset mesh method: this feature enables complete block-structed mesh refinement in TOSCA. 
- Additional scale-dependent (``dynamicLASD``) and anisotropic minimum dissipation (``amd``) LES models.  
- Additional high-order numerical scheme (``-central4`` with ``-hyperVisc`` parameter) to better capture the high-frequency tail of turbulence spectra. 
- Improved immersed boundary method: faster algorithms and time-varying wall modeling. 
- Lateral fringe region (``-yDampingLayer``): allows for gravity wave problems with wind direction changes and assimilated atmospheric data. 
- Data assimilation methods: direct/indirect profile assimilation techniques, as well as wavelet assimilation method (selectable through ``controllerType`` and ``controllerTypeT`` flags in *ABLProperties.dat*). 
- Additional example cases: *IndirectProfileAssimilationTest*, *WaveletAssimilationTest*, *MultiDomainOverset*, *PitchingAerofoilOversetTest*, *SuccessorOversetTest*, *SuccessorPeriodizedLateralFringeTest* and *TaylorGreenVortexTest*.
- Updated and extended documentation to incorporate new changes & features. 


V1.0.0
------

This is the first official release of TOSCA. An image of the original documentation is 
available at `TOSCA v1.0.0 Documentation <https://sebastipa.github.io/TOSCA/versions/1.0.0/index.html>`_