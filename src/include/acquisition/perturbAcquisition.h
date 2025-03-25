//! \file  perturbAcquisition.h
//! \brief perturbation fields w.r.t. an x-normal sample plane

#ifndef _PRTACQUISITION_H
#define _PRTACQUISITION_H

//! \brief Struct defining the ke fields
struct perturbFields
{
    // time controls
    PetscReal     avgStartTime;
    PetscReal     avgPrd;
    PetscReal     avgWindowTime;

    Cmpnts        **bkgU;
    PetscReal     **bkgP;
    PetscReal     **bkgT;
    PetscReal     xSample;
    PetscInt      kSample;
    PetscInt      avgWeight;

    // fields
    Vec           pertU;
    Vec           pertP;
    Vec           pertT;

    // flags
    PetscInt      read;
    PetscInt      initialized;
    PetscInt      movAvgActive;
};

#endif
