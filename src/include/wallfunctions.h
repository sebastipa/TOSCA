//! \file  wallfunctions.h
//! \brief wallfunction functions
#include "io.h"

// WALL FUNCTIONS
// ============================================================================================================= //

PetscReal uTauCabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn);

PetscReal utau_wf(PetscReal nu, PetscReal ks, PetscReal sb, PetscReal Ut_mag);

void      wallFunctionCabot(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallFunctionPowerlaw(double nu, double sc, double sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);
