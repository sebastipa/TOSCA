#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_force.h
//! \brief Turbine force computation and projection routines

#ifndef TRB_FORCE_H
#define TRB_FORCE_H

//! \brief Compute aerodynamic forces at the turbine mesh points
PetscErrorCode computeBladeForce(farm_ *farm);

//! \brief Project the wind turbine forces on the background mesh
PetscErrorCode projectBladeForce(farm_ *farm);

//! \brief Compute the tower forces on the background mesh
PetscErrorCode computeTowerForce(farm_ *farm);

//! \brief Compute project the tower forces on the background mesh
PetscErrorCode projectTowerForce(farm_ *farm);

//! \brief Compute and project the nacelle forces on the background mesh
PetscErrorCode projectNacelleForce(farm_ *farm);

//! \brief Transform the cartesian body force to contravariant
PetscErrorCode bodyForceCartesian2Contravariant(farm_ *farm);

#endif