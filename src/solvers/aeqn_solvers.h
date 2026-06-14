#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"

//! \file  aeqn_solvers.h
//! \brief Private Alpha-Water-equation solver function declarations.

#ifndef AEQN_SOLVERS_H
#define AEQN_SOLVERS_H

//! \brief low-order (upwind) face fluxes into lDivA
PetscErrorCode FormALowOrder(aeqn_ *aeqn);

//! \brief high-order (central4) face fluxes into lDivAHO
PetscErrorCode FormAHighOrder(aeqn_ *aeqn);

//! \brief compute low-order provisional cell update alpha^{LO} into lAlphaLO
PetscErrorCode ComputeLowOrderUpdate(aeqn_ *aeqn);

//! \brief apply compression flux for interface sharpening
PetscErrorCode AddCompressionFlux(aeqn_ *aeqn);

//! \brief apply mules limiter: compute lambda_f and accumulate limited rhs into Rhs
PetscErrorCode ApplyMULESLimiter(aeqn_ *aeqn);

//! \brief implicit snes solver for alpha-water equation
PetscErrorCode AeqnSNES(aeqn_ *aeqn);

//! \brief snes evaluation function
PetscErrorCode SNESFuncEvalA(SNES snes, Vec Alpha, Vec Rhs, void *ptr);

//! \brief compute rhs of alpha-water equation using current lAlpha
PetscErrorCode FormExplicitRhsA(aeqn_ *aeqn);

//! \brief solve alpha-water equation using runge-kutta 4
PetscErrorCode AeqnRK4(aeqn_ *aeqn);

#endif