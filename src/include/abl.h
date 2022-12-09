//! \file  abl.h
//! \brief ABL-related object definition

#ifndef _ABL_H_
#define _ABL_H_

struct abl_
{
    // flags
    PetscInt     controllerActive;               //!< activate velocity controller
    PetscInt     controllerActiveT;              //!< activate temperature controller
    PetscInt     coriolisActive;                 //!< activate coriolis force

    // physical quantities
    PetscReal    uTau;                           //!< friction Velocity
    PetscReal    hRough;                         //!< equivalent roughness length
    PetscReal    uRef;                           //!< reference velocity
    PetscReal    hRef;                           //!< reference height
    PetscReal    hGeo;                           //!< geostrophic height (required for geostrophic controller)
    PetscReal    hInv;                           //!< inversion height
    PetscReal    dInv;                           //!< inversion width
    PetscReal    gInv;                           //!< delta T across inversion layer
    PetscReal    tRef;                           //!< reference potential temperature
    PetscReal    gTop;                           //!< temperature gradient above the inversion layer
    PetscReal    vkConst;                        //!< von Karman constant
    PetscReal    smear;                          //!< Rampanelli Zardi model parameter
    PetscReal    fc;                             //!< Coriolis parameter (omegaEarth * sin(latitude) = 7.292115e-5 * sin(latitude))

    PetscReal   *cellLevels;                     //!< heights of the averaging planes

    // temperature controller
    PetscReal    *tDes;                          //!< initial temperature to be maintained

    // velocity controller (common)
    word         controllerType;                 //!< velocity controller type: write/read (writes in postProcessing/momentumSource, reads from momentumSource)
    PetscReal    relax;                          //!< source term relaxation factor
    PetscReal    alpha;                          //!< proportional over integral controller action ratio
    PetscReal    timeWindow;                     //!< time window of the integral part
    PetscReal    *totVolPerLevel;                //!< total volume at each cell level
    PetscInt     *totCelPerLevel;                //!< total number of cells per level
    PetscInt     *closestLabels;                 //!< closest heights w.r.t. controller height
    PetscReal    *levelWeights;                  //!< weights for variables interpolated at closest heights w.r.t. controller height
    PetscReal    controllerMaxHeight;            //!< max height of influence of the velocity controller

    // geostrophic damping for pressure controller
    PetscInt     geostrophicDampingActive;       //!< geosptrophic oscillation damping
	PetscReal    geoDampAvgDT;                   //!< average time step from simulation start
	Cmpnts       geoDampAvgS;                    //!< expected geostrophic velocity
	Cmpnts       geoDampUBar;                    //!< expected geostrophic velocity
	Cmpnts       *geoDampU;                      //!< average horizontal velocity at current iteration
    PetscReal    geoDampAlpha;                   //!< alpha = 1 critical damping (>1 over-damped, <1 under-damped)
    PetscReal    geoDampH;                       //!< H at which damping begins to be applied
    PetscReal    geoDampDelta;                   //!< rising distance for the application function
    PetscReal    geoDampC;                       //!< critical damping coefficient
    PetscReal    geoDampStart;                   //!< starting time of geostrophic damping
    PetscReal    geoDampWindow;                  //!< averaging window for the geostrophic velocity

    // geostrophic controller
    PetscInt     *closestLabelsGeo;              //!< closest heights w.r.t. controller height
    PetscReal    *levelWeightsGeo;               //!< weights for variables interpolated at closest heights w.r.t. controller height
    Cmpnts       uGeoBar;                        //!< desired geostrophic wind speed magnitude
    PetscReal    omegaBar;                       //!< rotation velocity
    PetscReal    hubAngle;                       //!< filtered angle at hub height
    PetscReal    geoAngle;                       //!< cumulated wind angle given by the sum of all rotations and filtered
    Cmpnts       a,b;                            //!< the two constant parts of the controller (a = geo forcing, b = wind angle controller)

    // read and average
    PetscInt     currentCloseIdx;                //!< save the current closest index at each iteration to speed up the interpolation search
    PetscReal    sourceAvgStartTime;             //!< if controllerType is 'average', average sources from this time value
    PetscReal    **preCompSources;               //!< table of given sources [ntimesteps][time|sourceX|sourceY|sourceZ] for velocity controller type = read.
    PetscInt     nSourceTimes;                   //!< number of times in the pre-computed sources
    Cmpnts       cumulatedSource;                //!< cumulated error of the velocity controller (equalt to gradP at steady state)
    PetscReal    avgTimeStep;                    //!< average time step from the momentum source file

    // z damping layer (Rayleigh damping)
    PetscReal    zDampingStart;                  //!< starting height of the Rayleigh damping layer
    PetscReal    zDampingEnd;                    //!< ending height of the Rayleigh damping layer
    PetscReal    zDampingAlpha;                  //!< damping paramter (it is equal to the Rayleigh viscosity at hEnd)
    Cmpnts       *uBarMeanZ;                     //!< array storing the reference velocity (averages along i at the k = 1 faces, size = my)
    PetscInt     avgWeight;
    PetscInt     zDampingAlsoXY;                 //!< damp also x and y velocity components (if 0 only z component is damped)
    PetscInt     zDampingXYType;                 //!< type 1 (default) averages at inlet, type (2) requires concurrent precursor and xDamping and uses planar averages

    // x damping layer (recycling fringe region)
    PetscReal    xDampingStart;                  //!< starting x of the fringe layer
    PetscReal    xDampingEnd;                    //!< ending x of the fringe layer
    PetscReal    xDampingDelta;                  //!< damping raise/decay distance (must be less than 0.5*(xDampingEnd - xDampingStart))
    PetscReal    xDampingAlpha;                  //!< damping paramter (see Inoue, Matheou, Teixeira 2014)

    // x damping layer controller parameters
    word         xDampingControlType;            //!< type of controller: alphaFixed or alphaOptimized
    PetscInt     advectionDampingType;           //!< type of advection damping (0: none, 1: LanzilaoMeyers2022)
    PetscReal    advDampingStart;                //!< starting x of the adv damping layer
    PetscReal    advDampingEnd;                  //!< ending x of the adv damping layer
    PetscReal    advDampingDeltaStart;           //!< damping raise/decay distance
    PetscReal    advDampingDeltaEnd;             //!< damping raise/decay distance
    PetscReal    xDampingTimeWindow;             //!< time constant for velocity filtering
    PetscReal    xDampingVBar;                   //!< desired y-velocity sampled from precursor
    PetscReal    vEnd, vStart;                   //!< line-averaged start and end y-velocity at the fringe region extrema
    PetscReal    xDampingError;                  //!< error on y-velocity at fringe end w.r.t precursor
    PetscReal    xDampingLineSamplingYmin;       //!< starting y of the sampling lines
    PetscReal    xDampingLineSamplingYmax;       //!< ending y of the sampling lines
    PetscReal    xDampingTimeStart;              //!< last time that alpha was modified
    PetscReal    xDampingDeltaV;                 //!< y-velocity jump across the fringe
    PetscReal    xDampingCoeff;                  //!< coeff = |Vend - Vstart| / alpha
    PetscInt    *closestLabelsFringe;            //!< closest height w.r.t. fringe controller height
    PetscReal   *levelWeightsFringe;             //!< weights for variables interpolated at closest heights w.r.t. fringe controller height

    // type of uBar computation
    PetscInt     xFringeUBarSelectionType;       //!< read type of fringe region in uBarSelectionType
    Cmpnts       **uBarInstX;                    //!< array storing the instantaneous velocity field for x damping layer
    PetscReal    **tBarInstX;                    //!< array storing the instantaneous temperature field for x damping layer
    PetscInt     **nProcsKLine;                  //!< number of processors in each k-line, used to average after MPI_Allreduce

    // side force for fringe region testing
    PetscReal    xStartSideF, xEndSideF;         //!< x start and ending coordinates of the region where the side force is applied
    PetscReal    zStartSideF, zEndSideF;         //!< z start and ending coordinates of the region where the side force is applied

    // turbulent flow initialization
    PetscReal    zPeak;
    PetscReal    deltaV;
    PetscReal    deltaU;
    PetscReal    Uperiods;
    PetscReal    Vperiods;

    // concurrent precursor
    precursor_    *precursor;                    //!< concurrent precursor data structure

    // access database
    access_       *access;
};

#endif

//! \brief Read from ABLProperties.dat and initialize the ABL parameters
PetscErrorCode InitializeABL(abl_ *abl);

//! \brief Evaluate geostrophic speed using Nieuwstadt model
PetscReal NieuwstadtGeostrophicWind(abl_ *abl);
