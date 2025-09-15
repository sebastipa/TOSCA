//! \file  abl.h
//! \brief ABL-related object definition

#ifndef _ABL_H_
#define _ABL_H_

struct abl_
{
    // flags
    PetscInt     controllerActive;               //!< activate velocity controller
    PetscInt     controllerActiveT;              //!< activate temperature controller
    PetscInt     controllerTypeMismatch;         //!< activate if the controllers used in successor and precursor differ
    PetscInt     coriolisActive;                 //!< activate coriolis force

    // physical quantities
    PetscReal    uTau;                           //!< friction Velocity
    PetscReal    hRough;                         //!< equivalent roughness length
    Cpt2D        uRef;                           //!< reference velocity
    PetscReal    hRef;                           //!< reference height
    PetscReal    hGeo;                           //!< geostrophic height (required for geostrophic controller)
    PetscReal    hInv;                           //!< inversion height
    PetscReal    dInv;                           //!< inversion width
    PetscReal    gInv;                           //!< delta T across inversion layer
    PetscReal    gABL;                           //!< temperature gradient below the inversion layer
    PetscReal    tRef;                           //!< reference potential temperature
    PetscReal    gTop;                           //!< temperature gradient above the inversion layer
    PetscReal    vkConst;                        //!< von Karman constant
    PetscReal    smear;                          //!< Rampanelli Zardi model parameter
    PetscReal    fc;                             //!< Coriolis parameter (omegaEarth * sin(latitude) = 7.292115e-5 * sin(latitude))
    PetscInt     perturbations;                  //!< add turbulent perturbations (only if initialization is set to ABLFlow)

    PetscReal   *cellLevels;                     //!< heights of the averaging planes

    // temperature controller
    PetscReal    *tDes;                          //!< initial temperature to be maintained
    word         controllerTypeT;                //!< initial or directProfileAssimilation
    
    // velocity controller (common)
    word         controllerType;                 //!< velocity controller type: write/read (writes in postProcessing/momentumSource, reads from momentumSource)
    word         controllerAction;
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
    PetscInt     mesoScaleInputActive;           //!< use mesoscale data for pressure controller uDes
    
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
    PetscReal    hubAnglePrev;                   //!< filtered angle at hub height at the previous timeStep
    PetscReal    geoAngle;                       //!< cumulated wind angle given by the sum of all rotations and filtered
    PetscReal    refHubAngle;                    //!< reference hub angle (typically 0 for canonical ABL cases)
    PetscReal    cumulatedAngle;
    Cmpnts       a,b;                            //!< the two constant parts of the controller (a = geo forcing, b = wind angle controller)

    // read and average
    PetscInt     currentCloseIdx;                //!< save the current closest index at each iteration to speed up the interpolation search
    PetscInt     currentCloseIdxT;    
    PetscReal    sourceAvgStartTime;             //!< if controllerType is 'average', average sources from this time value
    PetscReal    **preCompSources;               //!< table of given sources [ntimesteps][time|sourceX|sourceY|sourceZ] for velocity controller type = timeSeries/timeAverageSeries
    PetscReal    ***timeHtSources;               //!< table of given timeheight sources [time|sourceX|sourceY|sourceZ] at each cell level, controller type = timeHeightSeries
    PetscReal    ***timeHtSourcesT;               //!< table of given timeheight sources [time|source|] at each cell level, controller type = timeHeightSeries
    PetscInt     nSourceTimes;                   //!< number of times in the pre-computed sources
    Cmpnts       cumulatedSource;                //!< cumulated error of the velocity controller (equalt to gradP at steady state)
    Cmpnts       *cumulatedSourceHt;              //!< cumulated error of the velocity controller for every mesh level
    PetscReal    avgTimeStep;                    //!< average time step from the momentum source file

    // z damping layer (Rayleigh damping)
    PetscReal    zDampingStart;                  //!< starting height of the Rayleigh damping layer
    PetscReal    zDampingEnd;                    //!< ending height of the Rayleigh damping layer
    PetscReal    zDampingAlpha;                  //!< damping paramter (it is equal to the Rayleigh viscosity at hEnd)
    Cmpnts       *uBarMeanZ;                     //!< array storing the reference velocity (averages along i at the k = 1 faces, size = my)
    PetscInt     avgWeight;
    PetscInt     zDampingAlsoXY;                 //!< damp also x and y velocity components (if 0 only z component is damped)
    PetscInt     zDampingXYType;                 //!< type 1 (default) averages at inlet, type (2) requires concurrent precursor and xDamping and uses planar averages

    // kLeft damping layer
    PetscReal    kLeftPatchDist;                 //!< width of the kLeft Rayleigh damping layer
    PetscReal    kLeftDampingAlpha;              //!< kLeft Rayleigh damping coefficient
    Cmpnts       kLeftDampingUBar;               //!< kLeft bar velocity with respect to which the flow is damped
    PetscReal    kLeftDampingFilterHeight;       //!< above this height damping is unity (transitions to zero at kLeftDampingFilterHeight - kLeftDampingFilterWidth)
    PetscReal    kLeftDampingFilterWidth;        //!< width of transition region from no damping to damping

    // kRight damping layer
    PetscReal    kRightPatchDist;                 //!< width of the kRight Rayleigh damping layer
    PetscReal    kRightDampingAlpha;              //!< kRight Rayleigh damping coefficient
    Cmpnts       kRightDampingUBar;               //!< kRight bar velocity with respect to which the flow is damped
    PetscReal    kRightDampingFilterHeight;       //!< above this height damping is unity (transitions to zero at kRightDampingFilterHeight - kRightDampingFilterWidth)
    PetscReal    kRightDampingFilterWidth;        //!< width of transition region from no damping to damping

    // x damping layer (recycling fringe region)
    PetscReal    xDampingStart;                  //!< starting x of the fringe layer
    PetscReal    xDampingEnd;                    //!< ending x of the fringe layer
    PetscReal    xDampingDelta;                  //!< damping raise/decay distance (must be less than 0.5*(xDampingEnd - xDampingStart))
    PetscReal    xDampingAlpha;                  //!< damping paramter (see Inoue, Matheou, Teixeira 2014)

    // x damping layer controller parameters
    word         xDampingControlType;            //!< type of controller: alphaFixed or alphaOptimized
    PetscReal    advDampingXStart;                //!< starting x of the adv damping layer
    PetscReal    advDampingXEnd;                  //!< ending x of the adv damping layer
    PetscReal    advDampingXDeltaStart;           //!< damping raise/decay distance
    PetscReal    advDampingXDeltaEnd;             //!< damping raise/decay distance
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
    PetscReal    advDampingYStart;               //!< starting y of the adv damping y layer
    PetscReal    advDampingYEnd;                 //!< ending y of the adv damping y layer
    PetscReal    advDampingYDeltaStart;          //!< damping raise/decay distance
    PetscReal    advDampingYDeltaEnd;            //!< damping raise/decay distance

    // y damping layer (can only be applied together with the x damping layer)
    PetscReal    yDampingStart;                  //!< starting y of the fringe layer
    PetscReal    yDampingEnd;                    //!< ending y of the fringe layer
    PetscReal    yDampingDelta;                  //!< damping raise/decay distance (must be less than 0.5*(yDampingEnd - yDampingStart))
    PetscReal    yDampingAlpha;                  //!< damping paramter
    PetscInt     yDampingNumPeriods;             //!< number of periodizations in the streamwise direction of the x fringe region
    Vec          uBarInstY;                      //!< instantaneous bar velocity for y-fringe region (only used for concurent precursor)
    Vec          tBarInstY;                      //!< instantaneous bar temperature for y-fringe region (only used for concurent precursor)
   
    PetscInt     numSourceProc;                  //!< global to all procs - number of processors within the source (part of x fringe region within the lateral fringe region)
    PetscInt     *sourceProcList;                //!< global to all procs - list of processors in the source
    MPI_Comm     *yDamp_comm;                    //!< communicator that links each source processor to its corresponding periodization processors(destination) in the lateral fringe region - each source processor has a separate communicator for its set of source-destination processors
    PetscMPIInt  *srcCommLocalRank;              //!< global to communicator procs - local rank of the source processor within the source-destination communicator
    cellIds      *srcMinInd;                     //!< global to communicator procs - minimum k,j,i index of each processor within the source domain
    cellIds      *srcMaxInd;                     //!< global to communicator procs - maximum k,j,i index of each processor within the source domain
    PetscInt     *isdestProc;                    //!< flag indicating if a given processor is within the destination region of a source processor
    cellIds      **destMinInd;                   //!< local to each proc - minimum k,j,i index of a processor in the destination domain
    cellIds      **destMaxInd;                   //!< local to each proc - maximum k,j,i index of a processor in the destination domain
    PetscInt      *srcNumI;                      //!< global to communicator procs - number of i index for each processor in source
    PetscInt      *srcNumJ;                      //!< global to communicator procs - number of j index for each processor in source
    PetscInt      *srcNumK;                      //!< global to communicator procs - number of k index for each processor in source
    MPI_Request  *mapRequest;                    //!< MPI variable to perform non blocking broadcast operation
    Cmpnts       **velMapped;                    //!< one d array of the mapped velocity of each source processor
    PetscReal    **tMapped;
    PetscInt     **closestKCell;                 //!< closest 2 k index of the fictitious mesh(mapped from source) to the k indexes of the sucessor domain 
    PetscReal    **wtsKCell;                     //!< weights of the 2 closest k cells based on their distance 

    // type of uBar computation
    PetscInt     xFringeUBarSelectionType;       //!< read type of fringe region in uBarSelectionType
    Cmpnts       **uBarInstX;                    //!< array storing the instantaneous velocity field for x damping layer
    PetscReal    **tBarInstX;                    //!< array storing the instantaneous temperature field for x damping layer
    PetscInt     **nProcsKLine;                  //!< number of processors in each k-line, used to average after MPI_Allreduce

    // side force for fringe region testing
    PetscReal    xStartCanopy, xEndCanopy;         //!< x start and ending coordinates of the region where the side force is applied
    PetscReal    yStartCanopy, yEndCanopy;         //!< x start and ending coordinates of the region where the side force is applied
    PetscReal    zStartCanopy, zEndCanopy;         //!< z start and ending coordinates of the region where the side force is applied
    PetscReal    cftCanopy;                        //!< thrust coefficient of the entire wind canopy
    Cmpnts       diskDirCanopy;                    //!< disk direction of turbines inside the canopy

    // turbulent flow initialization
    PetscReal    zPeak;
    PetscReal    deltaV;
    PetscReal    deltaU;
    PetscReal    Uperiods;
    PetscReal    Vperiods;

    // mesoscale input parameters
    PetscReal    *timeV;
    PetscReal    *hV;
    PetscReal    *timeT;
    PetscReal    *hT;
    Cmpnts       **uMeso;
    Cmpnts       *uGeoPrev;
    PetscReal    **tMeso;

    PetscInt     numhV;
    PetscInt     numhT;
    PetscInt     numtV;
    PetscInt     numtT;

    PetscInt     **velInterpIdx;
    PetscReal    **velInterpWts;

    PetscInt     **tempInterpIdx;
    PetscReal    **tempInterpWts;

    PetscInt     lowestIndV;
    PetscInt     lowestIndT;
    PetscInt     highestIndV;
    PetscInt     highestIndT;
    PetscInt     lMesoIndV;
    PetscInt     hMesoIndV;
    PetscReal    lowestSrcHt;

    Cmpnts       *luMean;
    Cmpnts       *guMean;
    Cmpnts       *srcPA;
    
    PetscInt     closestTimeIndV;                        //!< closest index in time for the mesoscale timevarying data
    PetscReal    closestTimeWtV;
    PetscInt     closestTimeIndT;                        //!< closest index in time for the mesoscale timevarying data
    PetscReal    closestTimeWtT;
    
    PetscInt     polyOrder;
    PetscInt     polyOrderT;
    word         wtDist;
    PetscReal    **polyCoeffM;
    PetscReal    **polyCoeffT;

    //averaging 
    PetscInt     averageSource;
    PetscReal    currAvgtime;
    PetscReal    tAvgWindow;
    Cmpnts       *avgsrc;

    word         flType;
    word         flTypeT;
    PetscReal    *avgTotalStress;
    Cmpnts       *avgStress;
    PetscReal    *avgHeatFlux;
    PetscReal    hAvgTime;
    PetscReal    bottomSrcHtV;
    PetscReal    bottomSrcHtT;

    PetscReal    ablHt;
    PetscReal    ablHtStartTime;
    PetscReal    heatFluxSwitch;
    
    //continuous wavelet transform parameters 
    PetscInt     kernelRadius;
    PetscReal    sigma;  
    PetscReal    omega;

    //discrete wavelet transform parameters 
    word        waveName;
    word        waveTMethod;
    word        waveExtn;
    word        waveConv;
    PetscInt    waveLevel;
    PetscInt    waveletBlend;

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

//! \brief read the mesoscale driving velocity and potential temperature profile
PetscErrorCode readMesoScaleTemperatureData(abl_ *abl);

PetscErrorCode readMesoScaleVelocityData(abl_ *abl);

PetscErrorCode findVelocityInterpolationWeights(abl_ *abl);

PetscErrorCode findVelocityInterpolationWeightsOnePt(abl_ *abl);

PetscErrorCode findTemperatureInterpolationWeights(abl_ *abl);

PetscErrorCode initializeYDampingMapping(abl_ *abl);

PetscErrorCode setWeightsYDamping(abl_ *abl);

PetscErrorCode computeLSqPolynomialCoefficientMatrix(abl_ *abl);

PetscErrorCode computeLSqPolynomialCoefficientMatrixT(abl_ *abl);

PetscErrorCode findTimeHeightSeriesInterpolationWts(abl_ *abl);

PetscErrorCode findTimeHeightSeriesInterpolationWtsT(abl_ *abl);

PetscErrorCode findABLHeight(abl_ *abl);

PetscErrorCode waveletTransformContinuousVector(abl_ *abl, Cmpnts *srcPAIn, Cmpnts *srcPAOut, PetscInt sigLength);

PetscErrorCode waveletTransformContinuousScalar(abl_ *abl, PetscReal *srcPAIn, PetscReal *srcPAOut, PetscInt sigLength);

#if USE_PYTHON
    #ifdef __cplusplus
    extern "C" {
    #endif
        PetscErrorCode pywavedecVector(abl_ *abl, Cmpnts *src, Cmpnts *dest, PetscInt nlevels);
        PetscErrorCode pywavedecScalar(abl_ *abl, PetscReal *srcPAIn, PetscReal *srcPAOut, PetscInt sigLength);
    #ifdef __cplusplus
    }
    #endif
#endif
