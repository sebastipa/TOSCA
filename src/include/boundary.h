//! \file  boundary.h
//! \brief Boundary conditions header file.

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

//! \brief Structure dafining the type of boundary conditions for a scalar
struct scalarBC
{
    word   iLeft, iRight,  //!< type of boundary condition
           jLeft, jRight,
           kLeft, kRight;
    PetscReal iLval, iRval ,  //!< value (defined if BC prefix is 'fixed' only )
           jLval, jRval ,
           kLval, kRval ;
    PetscInt    iLWF , iRWF  ,  //!< wall function type (defined for velocity BC only)
           jLWF , jRWF  ,
           kLWF , kRWF  ;
};

//! \brief Structure dafining the type of boundary conditions for a vector
struct  vectorBC
{
    word   iLeft, iRight,  //!< type of boundary condition
           jLeft, jRight,
           kLeft, kRight;
    Cmpnts iLval, iRval ,  //!< value (defined if BC prefix is 'fixed' only )
           jLval, jRval ,
           kLval, kRval ;
    PetscInt    iLWF , iRWF  ,  //!< wall function type (defined for velocity BC only)
           jLWF , jRWF  ,
           kLWF , kRWF  ;
};

//! \brief Struct defining inflow information
struct inflowData
{
    PetscReal   *inflowTimes;                    //!< array with the available times (stored at the beginning)
    PetscInt       nInflowTimes;                   //!< number of times in the inflow database
    PetscInt       currentCloseIdx;                //!< save the current closest index at each iteration to speed up the interpolation search
};

//! \brief Struct storing the inlet function for a patch
struct inletFunctionTypes
{
    // type of inlet function
    PetscInt           typeU;
    PetscInt           typeT;
    PetscInt           typeNut;

    // type 1: power law (velocity only)
    Cmpnts        Uref;
    PetscReal        Href;
    PetscReal        uPrimeRMS;

    // type 2: log law + neutral stratification
    Cmpnts        Udir;                       //!< velocity direction vector
    PetscReal        roughness;                  //!< equivalent roughness size
    PetscReal        hInv;                       //!< height from bottom patch of inversion layer
    PetscReal        uTau;                       //!< friction velocity used to compute log profile

    // type 3/4: unsteady mapped/interpolated inflow (mapped assumes precursor and successor mesh are equal at the inflow location, interpolated assumes only periodicity)
    PetscInt           n1;                         //!< number of inflow cells along the 1st index (ordered as accessed: k,j,i. Eg. for k-normal patch this is n cells along j)
    PetscInt           n2;                         //!< number of inflow cells along the 2nd index (ordered as accessed: k,j,i. Eg. for k-normal patch this is n cells along i)
    PetscInt           n1wg;                       //!< n1 with ghosts (n1 + 2)
    PetscInt           n2wg;                       //!< n2 with ghosts (n2 + 2)
    PetscInt           prds1;                      //!< indicates how many times inflow data should be periodized along 1st index
    PetscInt           prds2;                      //!< indicates how many times inflow data should be periodized along 2nd index

    PetscInt           mapT;                       //!< flag telling if also T is mapped (if temperatureTransport is active is mandatory)
    PetscInt           mapNut;                     //!< flag telling if also nut is mapped (optional)

    PetscReal        ***inflowWeights;           //!< array of 4 weights for each cell center (ordered as a 2D plane array)
    cellIds       ***closestCells;            //!< array of 4 closest cells for each cell center
    PetscReal        width1;                     //!< inflow cell width in the 1st direction
    PetscReal        width2;                     //!< inflow cell width in the 2nd direction

    Cmpnts        **ucat_plane;               //!< cartesian velocity inflow data
    PetscReal        **t_plane;                  //!< temperature inflow data
    PetscReal        **nut_plane;                //!< nut inflow data
    inflowData    inflowU;                    //!< velocity inflow database inflormation for unsteadyMappedInflow BC
    inflowData    inflowT;                    //!< temperature inflow database inflormation for unsteadyMappedInflow BC
    inflowData    inflowNut;                  //!< temperature inflow database inflormation for unsteadyMappedInflow BC
};

//! \brief Struct storing inlet functions data
struct inletFunctions
{
    inletFunctionTypes *iLeft;
    inletFunctionTypes *iRight;
    inletFunctionTypes *jLeft;
    inletFunctionTypes *jRight;
    inletFunctionTypes *kLeft;
    inletFunctionTypes *kRight;
};

#endif

//! \brief Reads boundary conditions for a scalar field
PetscErrorCode readScalarBC(const word &location, const word &field, scalarBC *bc);

//! \brief Reads bundary conditions for a vector field
PetscErrorCode readVectorBC(const word &location, const word &field, vectorBC *bc);

//! \brief Set boundary conditions
PetscErrorCode SetBoundaryConditions(mesh_ *mesh);

//! \brief Checks available boundary conditions
PetscErrorCode checkBoundaryConditions(mesh_ *mesh);

//! \brief Set periodicconnectivity type
PetscErrorCode SetPeriodicConnectivity(mesh_ *mesh, word &meshFileName);

//! \brief Set wall models
PetscErrorCode SetWallModels(ueqn_ *ueqn);

//! \brief Update contravariant fluxes boundary conditions
PetscErrorCode UpdateContravariantBCs(ueqn_ *ueqn);

//! \brief Update cartesian boundary conditions
PetscErrorCode UpdateCartesianBCs(ueqn_ *ueqn);

//! \brief Update temperature boundary conditions
PetscErrorCode UpdateTemperatureBCs(teqn_ *teqn);

//! \brief Update effective viscosity boundary conditions
PetscErrorCode UpdateNutBCs(les_ *les);

//! \brief Update pressure boundary conditions
PetscErrorCode UpdatePressureBCs(peqn_ *peqn);

//! \brief Update pressure correction boundary conditions
PetscErrorCode UpdatePhiBCs(peqn_ *peqn);

//! \brief Update ibm boundary conditions
PetscErrorCode UpdateImmersedBCs(ibm_ *ibm);

//! \brief Update wall model for specified wall shear stress
PetscErrorCode UpdateWallModels(ueqn_ *ueqn);
