#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"
#include "../include/inflow.h"

//! \file  ueqn_sources.h
//! \brief Private momentum-equation source-term function declarations.

#ifndef UEQN_SOURCES_H
#define UEQN_SOURCES_H

//! \brief Update driving source terms
PetscErrorCode CorrectSourceTerms(ueqn_ *ueqn, PetscInt print);

//! \brief Map y-direction damping reference state
PetscErrorCode mapYDamping(ueqn_ *ueqn);

//! \brief Correct damping source terms
PetscErrorCode correctDampingSources(ueqn_ *ueqn);

//! \brief Apply buoyancy force
PetscErrorCode Buoyancy(ueqn_ *ueqn, PetscReal scale);

//! \brief Apply U control source terms
PetscErrorCode sourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply fringe region damping to momentum
PetscErrorCode dampingSourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply Coriolis force
PetscErrorCode Coriolis(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply canopy drag force
PetscErrorCode CanopyForce(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply forcing based on bulk velocity difference 
PetscErrorCode bulkGradPForcing(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply forcing based on prescribed mean pressure gradient
PetscErrorCode meanGradPForcing(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Apply hyper-viscosity to momentum equation
PetscErrorCode hyperViscosityU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);


#endif
