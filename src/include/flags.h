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
    PetscInt isCanopyActive;
    PetscInt isPrecursorSpinUp;
    PetscInt isConcurrentPrecursorActive;
    PetscInt isPvCatalystActive;
    PetscInt isKLeftRayleighDampingActive;
    PetscInt isKRightRayleighDampingActive;
    PetscInt isAdvectionDampingActive;
    PetscInt isAdvectionDampingYActive;
    PetscInt isSideForceActive;
    PetscInt isNonInertialFrameActive;
    PetscInt isGravityWaveModelingActive;
};

#endif
