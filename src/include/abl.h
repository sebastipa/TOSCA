//! \file  abl.h
//! \brief ABL-related object definition

#ifndef _ABL_H_
#define _ABL_H_

struct abl_
{
    // physical quantities
    PetscReal    uTau;                           //!< friction Velocity
    PetscReal    hRough;                         //!< equivalent roughness length
    PetscReal    uRef;                           //!< reference velocity
    PetscReal    hRef;                           //!< reference height
    PetscReal    hInv;                           //!< inversion height
    PetscReal    dInv;                           //!< inversion width
    PetscReal    gInv;                           //!< delta T across inversion layer
    PetscReal    tRef;                           //!< reference potential temperature
    PetscReal    gTop;                           //!< temperature gradient above the inversion layer
    PetscReal    vkConst;                        //!< von Karman constant
    PetscReal    smear;                          //!< Rampanelli Zardi model parameter
    PetscReal    fc;                             //!< Coriolis parameter (omegaEarth * sin(latitude) = 7.292115e-5 * sin(latitude))

    // velocity controller
    PetscReal    relax;                          //!< source term relaxation factor
    PetscReal    alpha;                          //!< proportional over integral controller action ratio
    PetscReal    timeWindow;                     //!< time window of the integral part
    word         controllerType;                 //!< velocity controller type: write/read (writes in postProcessing/momentumSource, reads from momentumSource)
    PetscReal    controllerHeight;               //!< max height of influence of the velocity controller
    PetscReal    sourceAvgStartTime;             //!< if controllerType is 'average', average sources from this time value
    PetscReal  **preCompSources;                 //!< table of given sources [ntimesteps][time|sourceX|sourceY|sourceZ] for velocity controller type = read.
    PetscInt     nSourceTimes;                   //!< number of times in the pre-computed sources
    PetscInt     currentCloseIdx;                //!< save the current closest index at each iteration to speed up the interpolation search
    PetscReal   *cellLevels;                     //!< heights of the averaging planes
    PetscReal   *totVolPerLevel;                 //!< total volume at each cell level
    PetscInt    *totCelPerLevel;                 //!< total number of cells per level
    Cmpnts      cumulatedSource;                 //!< cumulated error of the velocity controller (equalt to gradP at steady state)

    std::vector<std::vector<PetscReal>>
                 sourceHystory;                  //!< momentum source hystory (only if read)

    // z damping layer (Rayleigh damping)
    PetscReal    zDampingStart;                  //!< starting height of the Rayleigh damping layer
    PetscReal    zDampingEnd;                    //!< ending height of the Rayleigh damping layer
    PetscReal    zDampingAlpha;                  //!< damping paramter (it is equal to the Rayleigh viscosity at hEnd)
    Cmpnts       *uBarMeanZ;                     //!< array storing the reference velocity (averages along i at the k = 1 faces, size = my)
    PetscInt     avgWeight;
    PetscInt     zDampingAlsoXY;                 //!< damp also x and y velocity components (if 0 only z component is damped)

    // x damping layer (recycling fringe region)
    PetscReal    xDampingStart;                  //!< starting x of the fringe layer
    PetscReal    xDampingEnd;                    //!< ending x of the fringe layer
    PetscReal    xDampingDelta;                  //!< damping raise/decay distance (must be less than 0.5*(xDampingEnd - xDampingStart))
    PetscReal    xDampingAlpha;                  //!< damping paramter (see Inoue, Matheou, Teixeira 2014)
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
