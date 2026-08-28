#include "../include/base.h"
#include "../include/domain.h"
#include "../include/io.h"
#include "../include/inline.h"  
#include "../include/inflow.h"
#include "../sources/ueqn_sources.h"

//! \file  ueqn_solvers.h
//! \brief Private momentum-equation solver function declarations.

#ifndef UEQN_SOLVERS_H
#define UEQN_SOLVERS_H

//! \brief Viscous and divergence terms
PetscErrorCode FormU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale, PetscInt formMode = 0);

//! \brief SNES implicit time step for momentum equation
PetscErrorCode UeqnSNES(ueqn_ *ueqn);

//! \brief SNES function evaluation for momentum equation
PetscErrorCode SNESFuncEval(SNES snes, Vec Ucont, Vec Rhs, void *ptr);

//! \brief Explicit RHS for momentum equation
PetscErrorCode ExplicitRhsU(ueqn_ *ueqn, PetscInt formMode = 0);

//! \brief Forward-Euler time step for momentum equation
PetscErrorCode UeqnFE(ueqn_ *ueqn);

//! \brief RungeKutta 4 time step for momentum equation
PetscErrorCode UeqnRK4(ueqn_ *ueqn);

//! \brief ABCN time step for momentum equation (AB2 convection + CN diffusion)
PetscErrorCode UeqnABCN(ueqn_ *ueqn);

//! \brief CNRK3 time step for momentum equation (RK3 convection + CN diffusion)
PetscErrorCode UeqnRK3CN_SO(ueqn_ *ueqn);
PetscErrorCode UeqnRK3CN_W(ueqn_ *ueqn);

//! \brief MatShell operator for ABCN linear system
PetscErrorCode IMEXMatVec(Mat A, Vec v, Vec Av);

#endif
