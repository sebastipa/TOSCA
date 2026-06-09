#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"

//! \file  teqn_sources.h
//! \brief Private temperature-equation source-term function declarations.

#ifndef TEQN_SOURCES_H
#define TEQN_SOURCES_H

//! \brief Compute temperature control source term
PetscErrorCode CorrectSourceTermsT(teqn_ *teqn, PetscInt print);

//! \brief Compute tBar state for lateral damping region
PetscErrorCode correctDampingSourcesT(teqn_ *teqn);

//! \brief Apply fringe region damping
PetscErrorCode dampingSourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

//! \brief Apply temperature control
PetscErrorCode sourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

//! \brief Apply biharmonic hyperviscosity to the T equation RHS 
PetscErrorCode hyperViscosityT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

#endif
