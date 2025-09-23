#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_sampling.h
//! \brief Turbine velocity sampling routines

#ifndef TRB_SAMPLING_H
#define TRB_SAMPLING_H

//! \brief Compute wind velocity at the rotor mesh points
PetscErrorCode computeWindVectorsRotor(farm_ *farm);

//! \brief Compute wind velocity at the tower mesh points
PetscErrorCode computeWindVectorsTower(farm_ *farm);

//! \brief Compute wind velocity at the nacelle mesh points
PetscErrorCode computeWindVectorsNacelle(farm_ *farm);

//! \brief Compute wind velocity at the sample mesh points
PetscErrorCode computeWindVectorsSample(farm_ *farm);

#endif