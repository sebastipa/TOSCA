#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_parallel.h
//! \brief Turbine parallel computation routines

#ifndef TRB_PAR_H
#define TRB_PAR_H

//! \brief Initialize sphere cells and see what proc controls which turbine
PetscErrorCode initControlledCells(farm_ *farm);

//! \brief Initialize sphere cells and see what proc controls which sample point
PetscErrorCode initSampleControlledCells(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each rotor are controlled by this processor
PetscErrorCode findControlledPointsRotor(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each tower are controlled by this processor
PetscErrorCode findControlledPointsTower(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each nacelle are controlled by this processor
PetscErrorCode findControlledPointsNacelle(farm_ *farm);

//! \brief Discrimination algorithm: find out which points of each sample rig are controlled by this processor
PetscErrorCode findControlledPointsSample(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on the rotor
PetscErrorCode checkPointDiscriminationRotor(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on the tower
PetscErrorCode checkPointDiscriminationTower(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on the nacelle
PetscErrorCode checkPointDiscriminationNacelle(farm_ *farm);

//! \brief Debug check for the discrimination algorithm on sample points
PetscErrorCode checkPointDiscriminationSample(farm_ *farm);

#endif