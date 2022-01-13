# TOSCA
Toolbox fOr Stratified Convective Atmospheres

## Executables

    - tosca    : transient solver for stratified incompressible flows. Temperature stratification is accounted via Boussinesq approximation. 
    - windToPW : post processor for ParaView visualization. Writes data in XMF/HDF format. 
    
## Boundary Conditions

    - noSlip         : on all patches 
    - slip           : on all patches 
    - zeroGradient   : on all patches 
    - fixedGradient  : on all patches (only T equation)
    - inflowFunction : only on kLeft patch. Different types of inflow can be prescribed: ABL inflow using teorethical laws, mapped inflow from a precursor                                simulation that can be also periodized and/or interpolated if meshes are not consistent

## Wall Models

    - Shumann    : applies wall shear stress in the momentum equation according to similarity theory of Paulson and Obuhkov. Velocity BC is dependent. 
    - Cabot      : prescribes velocity BC according to Cabot formulation.
      
## Turbulence Models

    - Smagorinsky standard model.
    - Dynamic Smagorinsky LES turbulence model with cubic or lagrangian (Meneveau) averaging.
    
## Turbine Models

    - UADM : base actuator disk model with imposed Ct and yaw controller
    - ADM  : advanced actuator disk model with rotor dynamics and yaw/pitch/rotation controllers
    
## Acquisition System 

    - ABL wall-parallel planes-averaged statistics as a function of time 
    - field averages 
    - domain sections
    - probes with advanced parallel I/O writing 
    - turbine data with advanced parallel I/O writing

## Future Implementations:

    - Boundary condition for wall heat flux with similarity theories
    
## Notes:
  
    - test cases which contain an empty 'inflowDatabase' folder require the inflow database. A sample database can be downloaded at https://drive.google.com/file/d/17F5wtI5Jc1XGh8crmOVJYVXabC8iQXq1/view?usp=sharing, simply substitute the empty 'inflowDatabase' folder with the downloaded one. 
    
## Installation:

We suggest to create a folder named 'software', where the PETSc, HYPRE and TOSCA folders will be located. 

    1. Download PETSc (https://petsc.org/release/download/). Version 3.14.6 is recommended, available from (https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.14.6.tar.gz)
    
    2. Download HYPRE (https://github.com/hypre-space/hypre)
    
    3. TOSCA requires PETSc to be compiled with OpenMPI version 4.0.3 or earlier. To determine your version of OpenMPI, type
    
    $ ompi_info
    
    4. Configure PETSc (will automatically compile HYPRE). We suggest the following configure options: 
       './configure --with-cc=gcc --with-fc=0 --download-f2cblaslapack --with-mpi-dir='your--path--to--mpicc' --download-hypre='your--path--to--hypre' --with-64-bit-indices=1 --with-debugging=0'
       Be sure to have the correct path to your OpenMPI installation in your environment variables. 
    
    5. Make PETSc with 'make all'
    
    6. Test PETSc with 'make check'
    
    7. Compile TOSCA executables with 'make tosca' and 'make windToPW'

Credits and Copyright: Sebastiano Stipa - Arjun Ajay - Mohammad Hadi - The University of British Columbia
