//! \file  ablAcquisition.h
//! \brief Three layer model validation test. Depth averaged fields are gathered
//! on a gridded point mesh for subsequent post processing, namely velocity,
//! pressure and wind farm body force. The mesh must be cartesian and the vertical
//! direction is assumed to be the j-direction. The same hypothesis holds for the
//! ABL data acquisition, which performes spatial averages as a function of time
//! on the planes parallel to the wall defined by each cell's center.

#ifndef _ABLACQUISITION_H_
#define _ABLACQUISITION_H_

//! \brief Struct defining a gridded point mesh level for three layer model (used for its validation)
typedef struct
{
    word     levelName;                      //!< name of the mesh level
    Cmpnts   hStart;                         //!< starting height for depth average
    Cmpnts   hEnd;                           //!< ending height for depth average

    Cmpnts   **U;                            //!< depth average velocity averaged in time at the 3LM mesh points
    PetscReal   **P;                            //!< depth average pressure averaged in time at the 3LM mesh points
} level3LM;

//! \brief Struct defining a gridded point mesh for three layer model data (used for its validation)
struct data3LM
{
    PetscInt      nstw;                           //!< number of points in streamwise direction
    PetscInt      nspw;                           //!< number of points in the spanwise direction
    Cmpnts   upDir;                          //!< vertical direction unit vector
    Cmpnts   streamDir;                      //!< streamwise direction unit vector
    Cmpnts   spanDir;                        //!< spanwise direction unit vector
    Cmpnts   **points;                       //!< array of [nstw, nspw] where point coordinates are stored
    cellIds  **closestCells;                 //!< array of [nstw, nspw] where i and k ids of the averaging line are stored (j is meaningless)

    PetscReal   avgStartTime;                   //!< start time of acquisition system
    PetscReal   avgPrd;                         //!< acquisition time interval (overrides simulation time step if smaller)

    PetscInt      avgWeight;                      //!< number of averages taken up to this time step (snapshot weighting)

    level3LM **levels;                       //!< array of pointers to the three levels of the 3LM model
};

//! \brief Struct defining ABL simulation acquisition data
struct dataABL
{
    word     timeName;                       //!< name of the time directory where ABL averages are written
    PetscReal   avgStartTime;                   //!< start time of acquisition system
    PetscReal   avgPrd;                         //!< acquisition time interval (overrides simulation time step if smaller)                     //!< number of averages taken up to this time step (snapshot weighting)

    PetscReal   *cellLevels;                    //!< heights of the averaging planes
    PetscReal   *totVolPerLevel;                //!< total volume at each cell level
    PetscInt      *totCelPerLevel;                //!< total number of cells per level
    PetscReal   *UMean;                         //!< mean x velocity at each level
    PetscReal   *VMean;                         //!< mean y velocity at each level
    PetscReal   *WMean;                         //!< mean z velocity at each level
    PetscReal   *TMean;                         //!< mean temperature at each level
    PetscReal   *nutMean;                       //!< mean turbulent viscosity at each level
    Vec      UPrime;                         //!< fluctuating velocity field
    Vec      TPrime;                         //!< fluctuating temperature field

    PetscReal   *uuMean;
    PetscReal   *uvMean;
    PetscReal   *uwMean;
    PetscReal   *vvMean;
    PetscReal   *vwMean;
    PetscReal   *wwMean;

    PetscReal   *wuuMean;
    PetscReal   *wuvMean;
    PetscReal   *wuwMean;
    PetscReal   *wvvMean;
    PetscReal   *wvwMean;
    PetscReal   *wwwMean;

    PetscReal   *R11Mean;
    PetscReal   *R12Mean;
    PetscReal   *R13Mean;
    PetscReal   *R22Mean;
    PetscReal   *R23Mean;
    PetscReal   *R33Mean;

    PetscReal   *TuMean;
    PetscReal   *TvMean;
    PetscReal   *TwMean;

    PetscReal   *q1Mean;
    PetscReal   *q2Mean;
    PetscReal   *q3Mean;
};

#endif
