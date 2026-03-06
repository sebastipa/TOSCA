#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"
#include "../sources/teqn_sources.h"

//! \file  teqn_solvers.h
//! \brief Private temperature-equation solver function declarations.

#ifndef TEQN_SOLVERS_H
#define TEQN_SOLVERS_H

//! \brief Solve Teqn using RungeKutta 4
PetscErrorCode TeqnRK4(teqn_ *teqn);

//! \brief IMEX-CNAB time step for T equation (AB2 convection + backward-Euler diffusion)
PetscErrorCode TeqnIMEX(teqn_ *teqn);

//! \brief SNES evaluation function
PetscErrorCode TeqnSNES(SNES snes, Vec T, Vec Rhs, void *ptr);

//! \brief Compute RHS of temperature equation using current lTmprt
PetscErrorCode FormExplicitRhsT(teqn_ *teqn);

//! \brief RHS of the potential temperature transport equation
//!        formMode: 0 = full (conv+diff), 1 = conv-only, 2 = diff-only
PetscErrorCode FormT(teqn_ *teqn, Vec &Rhs, PetscReal scale, PetscInt formMode = 0);

//! \brief MatShell operator for IMEX system: A*v = v - dt*D(v)
PetscErrorCode IMEXTMatVec(Mat A, Vec v, Vec Av);

#endif
