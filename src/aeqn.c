//! \file  aeqn.c
//! \brief Contains Alpha equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

#include "solvers/aeqn_solvers.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorA(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeAEqn(aeqn_ *aeqn)
{
    if(aeqn != NULL)
    {
        // set pointer to mesh
        mesh_  *mesh  = aeqn->access->mesh;
        flags_ *flags = aeqn->access->flags;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        // alpha fields
        VecDuplicate(mesh->Nvert, &(aeqn->AlphaTmp));    VecSet(aeqn->AlphaTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Alpha));       VecSet(aeqn->Alpha,    0.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lAlpha));      VecSet(aeqn->lAlpha,   0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Alpha_o));     VecSet(aeqn->Alpha_o,  0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Rhs));         VecSet(aeqn->Rhs,      0.0);
        
        // density fields
        VecDuplicate(mesh->lAj,   &(aeqn->lRho));        VecSet(aeqn->lRho,   0.0);
        VecDuplicate(mesh->Cent,  &(aeqn->dRho));        VecSet(aeqn->dRho,      0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lRhoFace));    VecSet(aeqn->lRhoFace,0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lRhoFace_o));  VecSet(aeqn->lRhoFace_o,0.0);

        // MULES fields
        VecDuplicate(mesh->lCent, &(aeqn->lDivA));       VecSet(aeqn->lDivA,    0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lDivALO));     VecSet(aeqn->lDivALO,  0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lDivAHO));     VecSet(aeqn->lDivAHO,  0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lDivACor));    VecSet(aeqn->lDivACor, 0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lLambdaA));    VecSet(aeqn->lLambdaA, 1.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lAlphaLO));    VecSet(aeqn->lAlphaLO, 0.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lRplusA));     VecSet(aeqn->lRplusA,  1.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lRminusA));    VecSet(aeqn->lRminusA, 1.0);

        // alpha sub-cycling steps
        aeqn->nAlphaSubCycles = 3;
        PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-nAlphaSubCycles", &(aeqn->nAlphaSubCycles), PETSC_NULL);

        // numerical compression coefficient
        aeqn->compCoeff = 0.2;
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-alphaCompCoeff",  &(aeqn->compCoeff), PETSC_NULL);

        // read time discretization scheme
        readDictWord("control.dat", "-dAdtScheme", &(aeqn->ddtScheme));

        // create the SNES solver
        if(aeqn->ddtScheme=="BE" || aeqn->ddtScheme=="BDF2")
        {
            // default parameters
            aeqn->solverType       = "SNES";
            aeqn->absExitTol        = 1e-5;
            aeqn->relExitTol        = 1e-30;

            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolA",  &(aeqn->absExitTol), PETSC_NULL);
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolA",  &(aeqn->relExitTol), PETSC_NULL);

            SNESCreate    (mesh->MESH_COMM,&(aeqn->snesA));
            SNESMonitorSet(aeqn->snesA, SNESMonitorA, (void*)&(mesh->MESH_COMM), PETSC_NULL);

            // set the SNES evaluating function
            SNESSetFunction(aeqn->snesA, aeqn->Rhs, SNESFuncEvalA, (void *)aeqn);

            // create jacobian matrix
            MatCreateSNESMF(aeqn->snesA, &(aeqn->JA));
            MatMFFDSetType (aeqn->JA, MATMFFD_WP);
            SNESSetJacobian(aeqn->snesA, aeqn->JA, aeqn->JA, MatMFFDComputeJacobian, (void *)aeqn);

            // set SNES solver type
            SNESSetType(aeqn->snesA, SNESNEWTONTR);          // trust region: more robust when linear solve is imperfect
            //SNESSetType(aeqn->snesA, SNESNEWTONLS);

            // set SNES solve and step failures
            SNESSetMaxLinearSolveFailures  (aeqn->snesA,10000);
            SNESSetMaxNonlinearStepFailures(aeqn->snesA,10000);
            SNESKSPSetUseEW(aeqn->snesA, PETSC_TRUE);

            // set SNES Krylov Sub-Space parameters
            SNESKSPSetParametersEW(aeqn->snesA,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

            // SNES tolerances
            // 2nd arg: absolute tolerance
            // 3rd arg: relative tolerance
            // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
            // 5th arg: maximum number of iterations
            // 6th arg: maximum function evaluations
            SNESSetTolerances(aeqn->snesA, aeqn->absExitTol, aeqn->relExitTol, 1e-30, 20, 1000);

            // KSP solver
            SNESGetKSP(aeqn->snesA, &(aeqn->ksp));
            KSPGetPC  (aeqn->ksp,&(aeqn->pc));
            KSPSetType(aeqn->ksp, KSPGMRES);

            // maybe one day add a preconditioner (at least Jacobi)
            PCSetType(aeqn->pc, PCNONE);
            PetscReal rtol=aeqn->relExitTol, atol=aeqn->absExitTol, dtol=PETSC_DEFAULT;
            KSPSetTolerances(aeqn->ksp, rtol, atol, dtol, 1000);

            // previous solution for BDF2
            if(aeqn->ddtScheme == "BDF2")
            {
                VecDuplicate(mesh->Nvert, &(aeqn->Alpha_oo));  VecSet(aeqn->Alpha_oo, 0.0);
            }

            // print info
            PetscPrintf(mesh->MESH_COMM, "selected %s time scheme for Alpha-Water equation\n", aeqn->ddtScheme.c_str());
            PetscPrintf(mesh->MESH_COMM, " > relTolA = %e\n", aeqn->relExitTol);
            PetscPrintf(mesh->MESH_COMM, " > absTolA = %e\n", aeqn->absExitTol);
        }
        else if (aeqn->ddtScheme=="RK3")
        {
            aeqn->solverType       = "EXP";
        }
        else
        {
            char error[512];
            sprintf(error,  "unknown ddtScheme %s for Alpha-Water equation, available schemes are\n"
                            "    - RK3 (Runge Kutta 3)\n"
                            "    - BE (Backward Euler - BDF1)\n"
                            "    - BDF2 (2nd order Backward Differentiation Formula)\n",
                            aeqn->ddtScheme.c_str());
            fatalErrorInFunction("InitializeAEqn", error);
        }

        // read right damping region parameters 
        if(flags->isKRightAlphaDampingActive)
        {
            readSubDictDouble("waveProperties", "kRightAlphaDampingProperties", "patchDistance",   &(aeqn->kRightAlphaDampingDelta));
            readSubDictDouble("waveProperties", "kRightAlphaDampingProperties", "dampingCoeff",    &(aeqn->kRightAlphaDampingCoeff));
            readSubDictDouble("waveProperties", "kRightAlphaDampingProperties", "stillWaterLevel", &(aeqn->kRightAlphaDampingWaterLevel));
        }

        // read left damping region parameters
        if(flags->isKLeftAlphaDampingActive)
        {
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "dampingStart",     &(aeqn->kLeftAlphaDampingStart));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "dampingEnd",       &(aeqn->kLeftAlphaDampingEnd));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "dampingDelta",     &(aeqn->kLeftAlphaDampingDelta));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "dampingCoeff",     &(aeqn->kLeftAlphaDampingCoeff));

            readSubDictWord  ("waveProperties", "kLeftAlphaDampingProperties", "waveType",         &(aeqn->waveType));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "waveHeight",       &(aeqn->waveHeight));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "wavePeriod",       &(aeqn->wavePeriod));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "waveLevel",        &(aeqn->waveLevel));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "waveDirection",    &(aeqn->waveDirection));
            readSubDictDouble("waveProperties", "kLeftAlphaDampingProperties", "wavePhase",        &(aeqn->wavePhase));
            
            // validate wave type
            if(aeqn->waveType != "linear" && aeqn->waveType != "stokes2")
            {
                char error[512];
                sprintf(error, "unknown waveType %s, available types are:\n        linear  : Linear (Airy) wave theory\n        stokes2 : Stokes 2nd order theory\n", aeqn->waveType.c_str());
                fatalErrorInFunction("InitializeAEqn", error);
            }
            
            // compute derived wave quantities
            PetscReal g = 9.81;  
            PetscReal T = aeqn->wavePeriod;
            PetscReal d = aeqn->waveLevel;
            
            // angular frequency
            aeqn->waveOmega = 2.0 * M_PI / T;
            PetscReal omega = aeqn->waveOmega;
            
            // solve dispersion relation: omega^2 = g*k*tanh(k*d)
            // iterative solution for wave number k

            // deep water initial guess
            PetscReal k = omega * omega / g;  

            for(PetscInt iter=0; iter<200; iter++)
            {
                PetscReal f  = omega * omega - g * k * std::tanh(k * d);
                PetscReal df = -g * std::tanh(k * d) - g * k * d * (1.0 - std::tanh(k * d) * std::tanh(k * d));
                
                // Newton-Raphson
                k = k - f / df;  
            }

            // wave number and wave length
            aeqn->waveNumber = k;
            aeqn->waveLambda = 2.0 * M_PI / k;
            
            PetscPrintf(mesh->MESH_COMM, "\nkLeft alpha-water damping region parameters:\n");
            PetscPrintf(mesh->MESH_COMM, "   %s waves:\n", aeqn->waveType.c_str());
            PetscPrintf(mesh->MESH_COMM, "   -> H = %.4f m, T = %.4f s\n", aeqn->waveHeight, aeqn->wavePeriod);
            PetscPrintf(mesh->MESH_COMM, "   -> waveLevel = %.4f m (MWL position in domain)\n", aeqn->waveLevel);
            PetscPrintf(mesh->MESH_COMM, "   -> k = %.6f rad/m, lambda = %.4f m, omega = %.6f rad/s\n", aeqn->waveNumber, aeqn->waveLambda, aeqn->waveOmega);
            PetscPrintf(mesh->MESH_COMM, "   -> kd = %.3f (%s water)\n\n", aeqn->waveNumber * aeqn->waveLevel, aeqn->waveNumber * aeqn->waveLevel > 3.14 ? "deep" : (aeqn->waveNumber * aeqn->waveLevel < 0.5 ? "shallow" : "intermediate"));
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveAEqn(aeqn_ *aeqn)
{
    mesh_          *mesh  = aeqn->access->mesh;
    clock_         *clock = aeqn->access->clock;

    if(aeqn->ddtScheme=="BE" || aeqn->ddtScheme=="BDF2")
        AeqnSNES(aeqn);
    else if (aeqn->ddtScheme=="RK3")
        AeqnRK3(aeqn);

    // reset periodicity (strict enforcement is crucial)
    enforceInteriorCellPeriodicity(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar");
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode boundAlpha(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    PetscReal     ***alpha;

    DMDAVecGetArray(da,  aeqn->Alpha, &alpha);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                alpha[k][j][i] = PetscMax(0.0, PetscMin(1.0, alpha[k][j][i]));
            }
        }
    }

    // scatter to global to local 
    DMGlobalToLocalBegin(da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    DMDAVecRestoreArray(da,  aeqn->Alpha, &alpha);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode GradRho(aeqn_ *aeqn)
{
    mesh_         *mesh = aeqn->access->mesh;
    ueqn_         *ueqn = aeqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***icsi, ***ieta, ***izet,
                  ***jcsi, ***jeta, ***jzet,
                  ***kcsi, ***keta, ***kzet,
                  ***drho, ***cent;

    PetscReal     ***rho, ***nvert, ***meshTag;
    PetscReal     ***iaj, ***jaj, ***kaj;

    Cmpnts        gravity  = nSet(aeqn->access->constants->gravity);
    Cmpnts        vertical = nUnit(nScale(-1.0, gravity));
    PetscReal     Href     = aeqn->access->constants->Href;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     drdc, drde, drdz;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscInt periodic_i = mesh->i_periodic + mesh->ii_periodic,
             periodic_j = mesh->j_periodic + mesh->jj_periodic,
             periodic_k = mesh->k_periodic + mesh->kk_periodic;

    VecSet(aeqn->dRho, 0.0);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);
    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);
    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);
    DMDAVecGetArray(da,  mesh->lNvert,&nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,&meshTag);
    DMDAVecGetArray(da,  mesh->lIAj,  &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,  &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,  &kaj);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, aeqn->dRho, &drho);
    DMDAVecGetArray(da,  aeqn->lRho, &rho);

    PetscReal coeff_i, coeff_j, coeff_k;
    Cmpnts    cent_f;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal g11_i = (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y + icsi[k][j][i].z * icsi[k][j][i].z);
                PetscReal g12_i = (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y + ieta[k][j][i].z * icsi[k][j][i].z);
                PetscReal g13_i = (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y + izet[k][j][i].z * icsi[k][j][i].z);
                PetscReal g21_j = (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z);
                PetscReal g22_j = (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z);
                PetscReal g23_j = (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z);
                PetscReal g31_k = (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z);
                PetscReal g32_k = (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z);
                PetscReal g33_k = (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);

                // density gradient in the i-direction
                drdc = rho[k][j][i+1] - rho[k][j][i];

                if
                (
                    // j-right boundary -> use upwind only at the corner faces
                    (j==my-2 && i==mx-2) || isOversetIFace(k, j+1, i, i+1, meshTag)
                )
                {
                    drde = (rho[k][j  ][i  ] - rho[k][j-1][i  ] + rho[k][j  ][i+1] - rho[k][j-1][i+1]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (j == 1 && i == mx-2) || isOversetIFace(k, j-1, i, i+1, meshTag)
                )
                {
                    drde = (rho[k][j+1][i  ] - rho[k][j  ][i  ] + rho[k][j+1][i+1] - rho[k][j  ][i+1]) * 0.5;
                }
                else
                {
                    drde = (rho[k][j+1][i] - rho[k][j-1][i] + rho[k][j+1][i+1] - rho[k][j-1][i+1]) * 0.25;
                }

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (k == mz-2 && i==mx-2) || isOversetIFace(k+1, j, i, i+1, meshTag)
                )
                {
                    drdz = (rho[k][j][i  ] - rho[k-1][j][i  ] + rho[k][j][i+1] - rho[k-1][j][i+1]) * 0.5;
                }
                else if
                (
                    // k-left boundary  -> use upwind  only at the corner faces
                    (k == 1 && i==mx-2) || isOversetIFace(k-1, j, i, i+1, meshTag)
                )
                {
                    drdz = (rho[k+1][j][i  ] - rho[k][j][i  ] + rho[k+1][j][i+1] - rho[k][j][i+1]) * 0.5;
                }
                else
                {
                    drdz = (rho[k+1][j][i] - rho[k-1][j][i] + rho[k+1][j][i+1] - rho[k-1][j][i+1]) * 0.25;
                }  

                cent_f  = nScale(0.5, nSum(cent[k][j][i], cent[k][j][i+1]));
                coeff_i = nMag(gravity)*(nDot(vertical, cent_f) - Href);

                drho[k][j][i].x = (drdc * g11_i + drde *  g12_i + drdz * g13_i) * iaj[k][j][i] * coeff_i;

                // pressure gradient in the j-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (i == mx-2 && j==my-2) || isOversetJFace(k, j, i+1, j+1, meshTag)
                )
                {
                    drdc = (rho[k][j  ][i] - rho[k][j  ][i-1] + rho[k][j+1][i] - rho[k][j+1][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (i == 1 && j==my-2) || isOversetJFace(k, j, i-1, j+1, meshTag)
                )
                {
                    drdc = (rho[k][j  ][i+1] - rho[k][j  ][i] + rho[k][j+1][i+1] - rho[k][j+1][i]) * 0.5;
                }
                else
                {
                    drdc = (rho[k][j  ][i+1] - rho[k][j  ][i-1] + rho[k][j+1][i+1] - rho[k][j+1][i-1]) * 0.25;
                }

                drde = rho[k][j+1][i] - rho[k][j][i];

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (k == mz-2 && j== my-2) || isOversetJFace(k+1, j, i, j+1, meshTag)
                )
                {
                    drdz = (rho[k][j  ][i] - rho[k-1][j  ][i] + rho[k][j+1][i] - rho[k-1][j+1][i]) * 0.5;
                }
                else if
                (
                    // k-left boundary -> use upwind  only at the corner faces
                    (k == 1 && j== my-2) || isOversetJFace(k-1, j, i, j+1, meshTag)
                )
                {
                    drdz = (rho[k+1][j  ][i] - rho[k][j  ][i] + rho[k+1][j+1][i] - rho[k][j+1][i]) * 0.5;
                }
                else
                {
                    drdz = (rho[k+1][j  ][i] - rho[k-1][j  ][i] + rho[k+1][j+1][i] - rho[k-1][j+1][i]) * 0.25;
                }
                
                cent_f  = nScale(0.5, nSum(cent[k][j][i], cent[k][j+1][i]));
                coeff_j = nMag(gravity)*(nDot(vertical, cent_f) - Href);

                drho[k][j][i].y = (drdc * g21_j + drde * g22_j + drdz * g23_j) * jaj[k][j][i] * coeff_j;

                // pressure gradient in the k-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (i == mx - 2 && k==mz-2) || isOversetKFace(k, j, i+1, k+1, meshTag)
                )
                {
                    drdc = (rho[k  ][j][i] - rho[k  ][j][i-1] + rho[k+1][j][i] - rho[k+1][j][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (i == 1 && k == mz - 2) || isOversetKFace(k, j, i-1, k+1, meshTag)
                )
                {
                    drdc = (rho[k  ][j][i+1] - rho[k  ][j][i] + rho[k+1][j][i+1] - rho[k+1][j][i]) * 0.5;
                }
                else
                {
                    drdc = (rho[k  ][j][i+1] - rho[k  ][j][i-1] + rho[k+1][j][i+1] - rho[k+1][j][i-1]) * 0.25;
                }

                if
                (
                    // j-right boundary -> use upwind  only at the corner faces
                    (j == my - 2 && k==mz-2) || isOversetKFace(k, j+1, i, k+1, meshTag)
                )
                {
                    drde = (rho[k  ][j][i] - rho[k  ][j-1][i] + rho[k+1][j][i] - rho[k+1][j-1][i]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (j == 1 && k==mz-2)  || isOversetKFace(k, j-1, i, k+1, meshTag)
                )
                {
                    drde = (rho[k  ][j+1][i] - rho[k  ][j][i] + rho[k+1][j+1][i] - rho[k+1][j][i]) * 0.5;
                }
                else
                {
                    drde = (rho[k  ][j+1][i] - rho[k  ][j-1][i] + rho[k+1][j+1][i] - rho[k+1][j-1][i]) * 0.25;
                }

                drdz = (rho[k+1][j][i] - rho[k][j][i]);
                   
                cent_f  = nScale(0.5, nSum(cent[k][j][i], cent[k+1][j][i]));
                coeff_k = nMag(gravity)*(nDot(vertical, cent_f) - Href);

                drho[k][j][i].z = (drdc * g31_k + drde * g32_k + drdz * g33_k) * kaj[k][j][i] * coeff_k;

                // periodic: set to zero only on left boundaries since the contrav. velocity is not solved there
                // non-periodic: set to zero also on right boundaries since the contrav. velocity is not solved there
                if
                (
                    i==0 || (!mesh->i_periodic && !mesh->ii_periodic && i==mx-2) || isOversetIFace(k, j, i, i+1, meshTag)
                )
                {
                    drho[k][j][i].x = 0;
                }
                if
                (
                    j==0 || (!mesh->j_periodic && !mesh->jj_periodic && j==my-2) || isOversetJFace(k, j, i, j+1, meshTag)
                )
                {
                    drho[k][j][i].y = 0;
                }
                if
                (
                    k==0 || (!mesh->k_periodic && !mesh->kk_periodic && k==mz-2) || isOversetKFace(k, j, i, k+1, meshTag)
                )
                {
                    drho[k][j][i].z = 0;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);
    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);
    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert,&nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,&meshTag);
    DMDAVecRestoreArray(da,  mesh->lIAj,  &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,  &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,  &kaj);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(fda, aeqn->dRho, &drho);
    DMDAVecRestoreArray(da,  aeqn->lRho, &rho);

    return(0);
}