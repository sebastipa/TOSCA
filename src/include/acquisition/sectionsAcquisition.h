//! \file  sectionsAcquisition.h
//! \brief Sections along curvilinear coordinates object declaration

#ifndef _SECTIONSACQUISITION_H_
#define _SECTIONSACQUISITION_H_

//! \brief Struct defining a generic section to save
struct sections
{
    PetscInt           available;                  //!< does this struct contain data?
    PetscInt           nSections;                  //!< number of sections in k,j,i-direction
    PetscReal          timeStart;                  //!< start time of acquisition system
    word               intervalType;               //!< timeStep if sample at every time step, adjustableTime to enforce acquisition time interval
    PetscReal          timeInterval;               //!< acquisition time interval (overrides simulation time step if smaller), not used if intervalType = timeStep
    PetscInt           *indices;
    PetscReal          *coordinates;               //!< cutting coordinate of the section
    symmTensor         **symmTensorSec;            //!< temporary saves symmetric tensor section data
    Cmpnts             **vectorSec;                //!< temporarly saves vector section data
    PetscReal          **scalarSec;                //!< temporarly saves scalar section data
};

//! \brief Struct defining a user defined section not necessarily along curvilinear co-ordinates
struct uSections
{
    PetscInt           nodes;
    word               sectionName;
    Cmpnts             **coor;
    cellIds            **closestId;
    PetscInt           **hasCoor;                  //!< flag which sets if a given processor has this coordinate
    PetscReal          timeStart;                  //!< start time of acquisition system
    word               intervalType;               //!< timeStep if sample at every time step, adjustableTime to enforce acquisition time interval
    PetscReal          timeInterval;               //!< acquisition time interval (overrides simulation time step if smaller), not used if intervalType = timeStep
    PetscInt           nx;
    PetscInt           ny;
    Cmpnts             **vectorSec;                //!< temporarly saves vector section data
    PetscReal          **scalarSec;                //!< temporarly saves scalar section data
};


struct udSections
{
    PetscInt           available;                  //!< does this struct contain data?
    PetscInt           nSections;                  //!< number of sections in k,j,i-direction
    uSections          **uSection;
};

#endif
