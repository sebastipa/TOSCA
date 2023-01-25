//! \brief solution access flags

#ifndef _FLAGS_H_
#define _FLAGS_H_

struct flags_
{
    // simulation flags
    PetscInt isOversetActive;
    PetscInt isAdjustableTime;
    PetscInt isLesActive;
    PetscInt isTeqnActive;
    PetscInt isAquisitionActive;
    PetscInt isWindFarmActive;
    PetscInt isAblActive;
    PetscInt isIBMActive;
    PetscInt isZDampingActive;
    PetscInt isXDampingActive;
    PetscInt isYDampingActive;
    PetscInt isSideForceActive;
    PetscInt isPrecursorSpinUp;
    PetscInt isConcurrentPrecursorActive;
    PetscInt isPvCatalystActive;
};

#endif
