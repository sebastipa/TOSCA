#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"

//! \file  aeqn_solvers.h
//! \brief Private Alpha-Water-equation solver function declarations.

#ifndef AEQN_SOLVERS_H
#define AEQN_SOLVERS_H

//! \brief Viscous and divergence terms
PetscErrorCode FormA(aeqn_ *aeqn, Vec &Rhs, PetscReal scale, PetscInt formMode = 0);

//! \brief Implicit SNES solver for Alpha-Water equation
PetscErrorCode AeqnSNES(aeqn_ *aeqn);

//! \brief SNES evaluation function
PetscErrorCode SNESFuncEvalA(SNES snes, Vec Alpha, Vec Rhs, void *ptr);

//! \brief Compute RHS of Alpha-Water equation using current lAlpha
PetscErrorCode FormExplicitRhsA(aeqn_ *aeqn);

//! \brief Solve Alpha-Water equation using RungeKutta 4
PetscErrorCode AeqnRK4(aeqn_ *aeqn);

#endif