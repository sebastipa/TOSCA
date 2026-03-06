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

//! \brief MatShell operator for IMEX system: A*v = v - dt*scale*Visc(v)
PetscErrorCode IMEXMatVec(Mat A, Vec v, Vec Av);

//! \brief SNES evaluation function
PetscErrorCode UeqnSNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr);

//! \brief Compute explicit RHS for momentum equation
PetscErrorCode FormExplicitRhsU(ueqn_ *ueqn);

//! \brief Forward-Euler time step for momentum equation
PetscErrorCode UeqnEuler(ueqn_ *ueqn);

//! \brief Solve Ueqn using RungeKutta 4
PetscErrorCode UeqnRK4(ueqn_ *ueqn);

//! \brief Apply hyper-viscosity to momentum equation
PetscErrorCode hyperViscosityU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief IMEX-CNAB time step for momentum equation (AB2 convection + backward-Euler diffusion)
PetscErrorCode UeqnIMEX(ueqn_ *ueqn);

#endif
