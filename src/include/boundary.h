//! \file  boundary.h
//! \brief Boundary conditions header file.

#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_

//! \brief Structure dafining the type of boundary conditions for a scalar
struct scalarBC
{
    word      iLeft, iRight,  //!< type of boundary condition
              jLeft, jRight,
              kLeft, kRight;
    PetscReal iLval, iRval ,  //!< value (defined if BC prefix is 'fixed' only )
              jLval, jRval ,
              kLval, kRval ;
    PetscInt  iLWF , iRWF  ,  //!< wall function type (defined for velocity BC only)
              jLWF , jRWF  ,
              kLWF , kRWF  ;
};

//! \brief Structure dafining the type of boundary conditions for a vector

struct  vectorBC
{
    word     iLeft, iRight,  //!< type of boundary condition
             jLeft, jRight,
             kLeft, kRight;
    Cmpnts   iLval, iRval ,  //!< value (defined if BC prefix is 'fixed' only )
             jLval, jRval ,
             kLval, kRval ;
    PetscInt iLWF , iRWF  ,  //!< wall function type (defined for velocity BC only)
             jLWF , jRWF  ,
             kLWF , kRWF  ;
};

//! \brief Struct defining inflow information
struct inflowData
{
    PetscReal   *inflowTimes;                   //!< array with the available times (stored at the beginning)
    PetscInt    nInflowTimes;                   //!< number of times in the inflow database
    PetscInt    currentCloseIdx;                //!< save the current closest index at each iteration to speed up the interpolation search
};

//! \brief Struct storing the inlet function for a patch
struct inletFunctionTypes
{
    // type of inlet function
    PetscInt      typeU;
    PetscInt      typeT;
    PetscInt      typeNut;

    // type 1: power law (velocity only)
    Cmpnts        Uref;
    PetscReal     Href;
    PetscReal     uPrimeRMS;

    // type 2/5: log law/with veer + neutral stratification
    Cmpnts        Udir;                       //!< velocity direction vector
    PetscReal     roughness;                  //!< equivalent roughness size
    PetscReal     hInv;                       //!< height from bottom patch of inversion layer
    PetscReal     uTau;                       //!< friction velocity used to compute log profile
    PetscReal     dInv;                       //!< inversion width
    PetscReal     gInv;                       //!< delta T across inversion layer
    PetscReal     tRef;                       //!< reference potential temperature
    PetscReal     gTop;                       //!< temperature gradient above the inversion layer
	PetscReal     gABL;                       //!< temperature gradient below the inversion layer
    PetscReal     smear;                      //!< Rampanelli Zardi model parameter
    PetscReal     latitude;                   //!< for Nieuwstadt model
    PetscReal     fc;                         //!< Coriolis parameter

    // type 3/4: unsteady mapped/interpolated inflow (mapped assumes precursor and successor mesh are equal at the inflow location, interpolated assumes only periodicity)
    PetscInt      n1;                         //!< number of inflow cells along the 1st index (ordered as accessed: k,j,i. Eg. for k-normal patch this is n cells along j)
    PetscInt      n2;                         //!< number of inflow cells along the 2nd index (ordered as accessed: k,j,i. Eg. for k-normal patch this is n cells along i)
    PetscInt      n1wg;                       //!< n1 with ghosts (n1 + 2)
    PetscInt      n2wg;                       //!< n2 with ghosts (n2 + 2)
    PetscInt      prds1;                      //!< indicates how many times inflow data should be periodized along 1st index
    PetscInt      prds2;                      //!< indicates how many times inflow data should be periodized along 2nd index
    PetscInt      merge1;                     //!< average data at 5 top cells
	PetscInt      shift2;                     //!< apply shift on inflow slice data along direction 2

    //type 6: Synthetic turbulence inlet using random fourier series model.
    PetscInt      FSumNum;                    //!< number of fourier series to include in summations
    PetscInt      iterTKE;                    //!< number of iteration for Energy spectrum adjustment. should be <= 10.
    PetscReal     Urms;                       //!< user defined rms velocity
    PetscReal     kolLScale;                  //!< kolomagrov length scale in m
    PetscReal     intLScale;                //!< length scale in m of largest eddy
    PetscReal     dkn;                        //!< calculated wave vector interval, based on log scale
    PetscReal     kMax;                       //!< largerst wave number to include in driving energy spectrum
    PetscReal     kMin;                       //!< smallest wave number to include in driving energy spectrum
    PetscReal     kKol;                       //!< wave number assoicated with smallest eddies
    PetscReal     kEng;                       //!< wave number assoicated with most energetic eddies (peak of driving energy spectrum)
    Cmpnts        meanU;                      //!< mean velocity vector
    Cmpnts           *kn;                         //!< randomized wave vector
    Cmpnts           *Gn;                         //!< randomized direction unit vector
    PetscReal           *phaseN;                     //!< randomized phase angle
    PetscReal           *uMagN;                      //!< randomized velocity magnitudes
    PetscReal           *knMag;                      //!< wave number magnitudes
    PetscReal           *Ek;                         //!< Driving energy spectra values
    word                genType;                //!< which type of random generator is being used?


    PetscInt      mapT;                       //!< flag telling if also T is mapped (if temperatureTransport is active is mandatory)
    PetscInt      mapNut;                     //!< flag telling if also nut is mapped (optional)

    PetscReal     ***inflowWeights;           //!< array of 4 weights for each cell center (ordered as a 2D plane array)
    cellIds       ***closestCells;            //!< array of 4 closest cells for each cell center
    PetscReal     ***inflowWeights_1;         //!< array of 6 weights for each cell center (ordered as a 2D plane array) - j spline interp
    cellIds       ***closestCells_1;          //!< array of 6 closest cells for each cell center - j spline interp
    PetscReal     ***inflowWeights_2;         //!< array of 6 weights for each cell center (ordered as a 2D plane array) - i spline interp
    cellIds       ***closestCells_2;          //!< array of 6 closest cells for each cell center - i spline interp
    word          sourceType;                 //!< source mesh type (uniform or grading)
    word          interpMethod;               //!< interpolation method (linear or nullDiv)
    PetscReal     width1;                     //!< inflow cell width in the 1st direction
    PetscReal     width2;                     //!< inflow cell width in the 2nd direction
    PetscReal     inflowHeigth;               //!< inflow height for gradient extrapolation aloft

    Cmpnts        **ucat_plane;               //!< cartesian velocity inflow data
    PetscReal     **t_plane;                  //!< temperature inflow data
    PetscReal     **nut_plane;                //!< nut inflow data
    inflowData    inflowU;                    //!< velocity inflow database inflormation for unsteadyMappedInflow BC
    inflowData    inflowT;                    //!< temperature inflow database inflormation for unsteadyMappedInflow BC
    inflowData    inflowNut;                  //!< temperature inflow database inflormation for unsteadyMappedInflow BC

    Cmpnts       *uBarAvgTopX;                //!< velocity average from inflow database at top 10 points
    PetscReal    *tBarAvgTopX;                //!< temperature average from inflow database at top 10 points
    PetscReal    *avgTopPointCoords;          //!< z coordinates of the 10 average top points
    PetscReal    avgTopDelta;                 //!< length of the five top cells
    PetscReal    avgTopLength;                //!< height of the inflow slices

	PetscReal    shiftSpeed;                  //!< speed of the i-direction shift
    PetscReal    *ycent;                      //!< vector storing y cell center coordinates (assuming they don't vary vertically)
    PetscInt     *yIDs;                       //!< IDs for the right interpolation point (the left would be IDs[i]-1)
    PetscReal    *yWeights;                   //!< vector storing the right weight for the interpolation (left weight is 1.0 - right weight)
    MPI_Comm     IFFCN_COMM;                  //!< communicator involving all processors touching k-boundaries

    // type 6: i-dir sinusoidal inflow
    PetscReal    amplitude;                   //!< oscillation amplitude w.r.t. reference velocity magnitude
    PetscReal    periods;                     //!< number of periods in the spanwise direction
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

//! \brief Reads boundary conditions for a scalarMoments field multiple scalar fields
PetscErrorCode readScalarMomentsBC(const word &location, const word &field, SMObj_ *smObject);

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
PetscErrorCode UpdateWallModelsU(ueqn_ *ueqn);

//! \brief Update wall model for specified wall heat flux
PetscErrorCode UpdateWallModelsT(teqn_ *teqn);

//! \brief read surface temperature and obhukhov length data
PetscErrorCode readSurfaceTempData(Shumann *wm);

//! \brief Update Scalar Moment boundary conditions for scalar moment object ii
PetscErrorCode UpdateScalarMomentBCs(sm_ *sm, PetscInt ii);
