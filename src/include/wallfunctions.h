//! \file  wallfunctions.h
//! \brief wallfunction functions
#include "io.h"

// WALL FUNCTIONS
// ============================================================================================================= //

PetscReal uTauCabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn);

PetscReal utau_wf(PetscReal nu, PetscReal ks, PetscReal sb, PetscReal Ut_mag);

PetscReal uTauCabotRoughness(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn, PetscReal ks);

PetscReal uTauCabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn);

void      wallFunctionCabot(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallFunctionCabotRoughness(PetscReal nu, PetscReal ks, PetscReal sc, PetscReal sb,
            Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallFunctionPowerlaw(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallFunctionPowerlawAPG(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness, PetscReal kappa, Cmpnts Ua,
                                  Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf, PetscReal dpdx, PetscReal dpdy, PetscReal dpdz);

void      wallFunctionLogLawAPG(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness, PetscReal kappa, Cmpnts Ua,
                                  Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf, PetscReal dpdx, PetscReal dpdy, PetscReal dpdz);

void      slipBC(PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, Cmpnts nf);

void      wallFunctionSchumann(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness,
                                    PetscReal kappa, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallShearVelocityBC(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness,
                                    PetscReal kappa, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallShearVelocityBCQuadratic(PetscReal nu,  PetscReal sd, PetscReal sc, PetscReal sb, PetscReal roughness,
                                    PetscReal kappa, Cmpnts Ua, Cmpnts Ud, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void      wallShearGhostVelocityBC(PetscReal nu,  PetscReal sd, PetscReal sc, PetscReal sb, PetscReal roughness,
                                    PetscReal kappa, Cmpnts Ua, Cmpnts Ud, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf);

void uStarShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &roughness,
    PetscReal &gammaM, PetscReal &kappa, PetscReal &qwall, PetscReal &thetaRef,
    PetscReal &uStar, PetscReal &phiM, PetscReal &L
);

void qWallShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &z0,
    PetscReal &gammaM, PetscReal &gammaH, PetscReal &alphaH,
    PetscReal &thetaRef, PetscReal &deltaTheta, PetscReal &kappa,
    PetscReal &qWall, PetscReal &uStar, PetscReal &phiM, PetscReal &phiH, PetscReal &L
);
