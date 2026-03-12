//! \file  teqn.c
//! \brief Contains T equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

#include "sources/teqn_sources.h"
#include "solvers/teqn_solvers.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorT(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeTEqn(teqn_ *teqn)
{
    if(teqn != NULL)
    {
        // set pointer to mesh
        mesh_  *mesh  = teqn->access->mesh;
        flags_ *flags = teqn->access->flags;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        // initialize p-tilde formulation of momentum equation
        teqn->pTildeFormulation = 0;
        PetscOptionsGetInt(PETSC_NULL, PETSC_NULL,  "-pTildeBuoyancy", &(teqn->pTildeFormulation),   PETSC_NULL);

        VecDuplicate(mesh->Nvert, &(teqn->TmprtTmp));    VecSet(teqn->TmprtTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt));       VecSet(teqn->Tmprt,    0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt_o));     VecSet(teqn->Tmprt_o,  0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Rhs));         VecSet(teqn->Rhs,      0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt));      VecSet(teqn->lTmprt,   0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt_o));    VecSet(teqn->lTmprt_o, 0.0);

        VecDuplicate(mesh->lCent, &(teqn->lDivT));       VecSet(teqn->lDivT,    0.0);
        VecDuplicate(mesh->lCent, &(teqn->lViscT));      VecSet(teqn->lViscT,   0.0);
        VecDuplicate(mesh->Nvert, &(teqn->sourceT));     VecSet(teqn->sourceT,  0.0);

        if(teqn->pTildeFormulation)
        {
            VecDuplicate(mesh->lAj,   &(teqn->lRhoK));       VecSet(teqn->lRhoK,   0.0);
            VecDuplicate(mesh->Cent,  &(teqn->ghGradRhok));   VecSet(teqn->ghGradRhok,   0.0);
            VecDuplicate(mesh->Cent,  &(teqn->ghGradRhok_o)); VecSet(teqn->ghGradRhok_o, 0.0);
        }

        if(flags->isIBMActive)
        {
            VecDuplicate(mesh->lCent, &(teqn->lViscIBMT));   VecSet(teqn->lViscIBMT,  0.0);
        }

        // read time discretization scheme
        readDictWord("control.dat", "-dTdtScheme", &(teqn->ddtScheme));

        // create the SNES solver
        if(teqn->ddtScheme=="BE" || teqn->ddtScheme=="BDF2")
        {
            // default parameters
            teqn->absExitTol        = 1e-5;
            teqn->relExitTol        = 1e-30;

            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolT",  &(teqn->absExitTol), PETSC_NULL);
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolT",  &(teqn->relExitTol), PETSC_NULL);

            SNESCreate(mesh->MESH_COMM,&(teqn->snesT));
            SNESMonitorSet(teqn->snesT, SNESMonitorT, (void*)&(mesh->MESH_COMM), PETSC_NULL);

            // set the SNES evaluating function
            SNESSetFunction(teqn->snesT, teqn->Rhs, TeqnSNES, (void *)teqn);

            // create jacobian matrix
            MatCreateSNESMF(teqn->snesT, &(teqn->JT));
            SNESSetJacobian(teqn->snesT, teqn->JT, teqn->JT, MatMFFDComputeJacobian, (void *)teqn);

            // set SNES solver type
            //SNESSetType(teqn->snesT, SNESNEWTONTR);           //SNESTR
            SNESSetType(teqn->snesT, SNESNEWTONLS);        //SNESLS is better for stiff PDEs such as the one including IB but slower

            // set SNES solve and step failures
            SNESSetMaxLinearSolveFailures(teqn->snesT,10000);
            SNESSetMaxNonlinearStepFailures(teqn->snesT,10000);
            SNESKSPSetUseEW(teqn->snesT, PETSC_TRUE);

            // set SNES Krylov Sub-Space parameters
            SNESKSPSetParametersEW(teqn->snesT,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

            // SNES tolerances
            // 2nd arg: absolute tolerance
            // 3rd arg: relative tolerance
            // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
            // 5th arg: maximum number of iterations
            // 6th arg: maximum function evaluations
            SNESSetTolerances(teqn->snesT, teqn->absExitTol, teqn->relExitTol, 1e-30, 20, 1000);

            SNESGetKSP(teqn->snesT, &(teqn->ksp));
            KSPGetPC(teqn->ksp,&(teqn->pc));

            // set KSP solver type
            KSPSetType(teqn->ksp, KSPGMRES);

            //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);    //2009.09.22 poor performance
            //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);      //2009.09.22

            //KSPFischerGuess itg;
            //KSPFischerGuessCreate(ksp,1,100,&itg);
            //KSPSetFischerGuess(ksp, itg);                  //2009.09.22

            //KSPGMRESSetPreAllocateVectors(ksp);            --> crazy thing consumes memory

            PCSetType(teqn->pc, PCNONE);
            PetscReal rtol=teqn->relExitTol, atol=teqn->absExitTol, dtol=PETSC_DEFAULT;
            KSPSetTolerances(teqn->ksp, rtol, atol, dtol, 1000);

            // BDF2 needs T^{n-1} storage
            if(teqn->ddtScheme == "BDF2")
            {
                VecDuplicate(mesh->Nvert, &(teqn->Tmprt_oo));  VecSet(teqn->Tmprt_oo, 0.0);
            }

            // print info
            if(teqn->ddtScheme == "BE")
                PetscPrintf(mesh->MESH_COMM, "selected BE time scheme for T equation\n");
            else
                PetscPrintf(mesh->MESH_COMM, "selected BDF2 time scheme for T equation\n");
            PetscPrintf(mesh->MESH_COMM, " > relTolT = %e\n", teqn->relExitTol);
            PetscPrintf(mesh->MESH_COMM, " > absTolT = %e\n", teqn->absExitTol);
        }
        else if (teqn->ddtScheme=="RK4")
        {

        }
        else if (teqn->ddtScheme=="BEAB")
        {
            // IMEX-CNAB: Adams-Bashforth 2 for convection (explicit) + backward-Euler for diffusion (implicit).
            // Solved using a standalone KSP for the linear system A*T^{n+1} = b,
            // where A*v = v - dt*D(v) is assembled via MatShell, and
            //       b = T^n + dt*bT is prebuilt once per time step.

            teqn->gmresRestart = 30;
            PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-kspGMRESRestartT", &(teqn->gmresRestart), PETSC_NULL);

            // KSP type: BiCGStab or GMRES set via -kspTypeT in control.dat
            readDictWord("control.dat", "-kspTypeT", &(teqn->kspType));

            if(teqn->kspType != "BiCGStab" && teqn->kspType != "GMRES")
            {
                char error[512];
                sprintf(error,
                    "unknown kspTypeT %s for T equation (IMEX), available types are\n"
                    "    BiCGStab\n"
                    "    GMRES",
                    teqn->kspType.c_str());
                fatalErrorInFunction("InitializeTEqn", error);
            }

            // read KSP exit tolerances
            teqn->relExitTol = 1e-30;
            teqn->absExitTol = 1e-5;
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolT", &(teqn->relExitTol), PETSC_NULL);
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolT", &(teqn->absExitTol), PETSC_NULL);

            // Create standalone KSP for the IMEX linear system A*T^{n+1} = b.
            {
                PetscInt n, N;
                VecGetLocalSize(teqn->Tmprt, &n);
                VecGetSize(teqn->Tmprt, &N);
                MatCreateShell(mesh->MESH_COMM, n, n, N, N, (void*)teqn, &(teqn->JvIMEX));
                MatShellSetOperation(teqn->JvIMEX, MATOP_MULT, (void(*)(void))IMEXTMatVec);
            }

            KSPCreate(mesh->MESH_COMM, &(teqn->kspIMEX));
            KSPSetOperators(teqn->kspIMEX, teqn->JvIMEX, teqn->JvIMEX);

            {
                PC pcIMEX;
                KSPGetPC(teqn->kspIMEX, &pcIMEX);
                PCSetType(pcIMEX, PCNONE);
            }

            if(teqn->kspType == "BiCGStab")
            {
                KSPSetType(teqn->kspIMEX, KSPBCGSL);
                KSPBCGSLSetEll(teqn->kspIMEX, 2);
            }
            else
            {
                KSPSetType(teqn->kspIMEX, KSPGMRES);
                KSPGMRESSetRestart(teqn->kspIMEX, teqn->gmresRestart);
            }

            // use explicit-Euler estimate as non-zero initial guess
            KSPSetInitialGuessNonzero(teqn->kspIMEX, PETSC_TRUE);
            KSPSetTolerances(teqn->kspIMEX, teqn->relExitTol, teqn->absExitTol, 1e30, 1000);

            // allocate IMEX buffers
            VecDuplicate(mesh->Nvert, &(teqn->bT));        VecSet(teqn->bT,        0.0);
            VecDuplicate(mesh->Nvert, &(teqn->RhsConv));   VecSet(teqn->RhsConv,   0.0);
            VecDuplicate(mesh->Nvert, &(teqn->RhsConv_o)); VecSet(teqn->RhsConv_o, 0.0);

            // print info
            PetscPrintf(mesh->MESH_COMM, "selected IMEX time scheme for T equation\n");
            PetscPrintf(mesh->MESH_COMM, " > scheme: IMEX-CNAB (AB2 conv + backward-Euler diff)\n");
            PetscPrintf(mesh->MESH_COMM, " > relTolT = %e\n", teqn->relExitTol);
            PetscPrintf(mesh->MESH_COMM, " > absTolT = %e\n", teqn->absExitTol);
            PetscPrintf(mesh->MESH_COMM, " > kspTypeT = %s\n", teqn->kspType.c_str());
            if(teqn->kspType == "GMRES")
            {
                PetscPrintf(mesh->MESH_COMM, " > kspGMRESRestartT = %d\n", teqn->gmresRestart);
            }
        }
        else
        {
            char error[512];
            sprintf(error, "unknown ddtScheme %s for T equation, available schemes are\n    1. BE\n    2. RK4\n    3. BEAB\n    4. BDF2", teqn->ddtScheme.c_str());
            fatalErrorInFunction("InitializeTEqn", error);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveTEqn(teqn_ *teqn)
{
    mesh_          *mesh  = teqn->access->mesh;
    clock_         *clock = teqn->access->clock;

    // set the right hand side to zero
    VecSet (teqn->Rhs, 0.0);

    if(teqn->ddtScheme=="BE" || teqn->ddtScheme=="BDF2")
    {
        PetscReal     norm;
        PetscInt      iter;
        PetscReal     ts, te;

        PetscTime(&ts);
        PetscPrintf(mesh->MESH_COMM, "TRSNES: Solving for T, Initial residual = ");
   
        // compute initial guess: conv-only forward Euler (avoids explicit diffusion instability)
        // VecSet(teqn->Rhs, 0.0);
        // FormT(teqn, teqn->Rhs, 1.0, 1);                    // conv-only, scale=1
        VecCopy(teqn->Tmprt_o, teqn->TmprtTmp);            // TmprtTmp = T^n
        //VecAXPY(teqn->TmprtTmp, clock->dt, teqn->Rhs);     // TmprtTmp = T^n + dt * N(T^n)

        SNESSolve(teqn->snesT, PETSC_NULL, teqn->TmprtTmp);
        SNESGetFunctionNorm(teqn->snesT, &norm);
        SNESGetIterationNumber(teqn->snesT, &iter);

        // report total inner linear iterations (= total MatMFFD FormT calls for J*v)
        PetscInt linIter = 0;
        SNESGetLinearSolveIterations(teqn->snesT, &linIter);

        VecCopy(teqn->TmprtTmp, teqn->Tmprt);

        DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

        PetscTime(&te);
        PetscPrintf(mesh->MESH_COMM, "Final residual = %e, Iterations = %ld, Linear iterations = %ld, Elapsed Time = %lf\n", norm, iter, linIter, te-ts);
    }
    else if (teqn->ddtScheme=="RK4")
    {
        TeqnRK4(teqn);
    }
    else if (teqn->ddtScheme=="BEAB")
    {
        TeqnIMEX(teqn);
    }

    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ghGradRhoK(teqn_ *teqn)
{
    mesh_         *mesh = teqn->access->mesh;
    ueqn_         *ueqn = teqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***icsi, ***ieta, ***izet,
                  ***jcsi, ***jeta, ***jzet,
                  ***kcsi, ***keta, ***kzet,
                  ***coor, ***db;
    Vec           Coor;

    PetscReal     ***tmprt, ***rhok, ***nvert;

    PetscReal     ***iaj, ***jaj, ***kaj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    const PetscReal  g = -9.81, tRef = ueqn->access->abl->tRef;;

    PetscReal     dbdc, dbde, dbdz;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

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
    DMDAVecGetArray(da,  mesh->lIAj,  &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,  &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,  &kaj);

    DMDAVecGetArray(da,  teqn->lTmprt, &tmprt);
    DMDAVecGetArray(fda, teqn->ghGradRhok,  &db);

    VecSet(teqn->lRhoK, 0.0);
    DMDAVecGetArray(da, teqn->lRhoK,  &rhok);

    // create the field ghrhok = (rhok/rho)
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                rhok[k][j][i]
                =
                (2* tRef - tmprt[k][j][i]) / tRef;
            }
        }
    }

    DMDAVecRestoreArray(da, teqn->lRhoK,  &rhok);

    // scatter Phi from global to local
    DMLocalToLocalBegin(da, teqn->lRhoK, INSERT_VALUES, teqn->lRhoK);
    DMLocalToLocalEnd(da, teqn->lRhoK, INSERT_VALUES, teqn->lRhoK);

    DMDAVecGetArray(da, teqn->lRhoK,  &rhok);

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

                // pressure gradient in the i-direction
                if( i==mx-2 && mesh->ii_periodic)
                {
                    dbdc = rhok[k][j][mx+1] - rhok[k][j][i];
                }
                else
                {
                    dbdc = rhok[k][j][i+1] - rhok[k][j][i];
                }

                if
                (
                    // j-right boundary -> use upwind only at the corner faces
                    (
                        j==my-2 &&
                        (
                            i==mx-2
                        )
                    )
                )
                {
                    dbde = (rhok[k][j  ][i  ] - rhok[k][j-1][i  ] + rhok[k][j  ][i+1] - rhok[k][j-1][i+1]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (
                        j == 1 &&
                        (
                            i == mx - 2
                        )
                     )
                )
                {
                    dbde = (rhok[k][j+1][i  ] - rhok[k][j  ][i  ] + rhok[k][j+1][i+1] - rhok[k][j  ][i+1]) * 0.5;
                }
                else
                {
                    dbde = (rhok[k][j+1][i] - rhok[k][j-1][i] + rhok[k][j+1][i+1] - rhok[k][j-1][i+1]) * 0.25;
                }

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (
                        k == mz - 2 &&
                        (
                            i==mx-2
                        )
                    )
                )
                {
                    dbdz = (rhok[k][j][i  ] - rhok[k-1][j][i  ] + rhok[k][j][i+1] - rhok[k-1][j][i+1]) * 0.5;
                }
                else if
                (
                    // k-left boundary  -> use upwind  only at the corner faces
                    (
                        k == 1 &&
                        (
                            i==mx-2
                        )
                    )
                )
                {
                    dbdz = (rhok[k+1][j][i  ] - rhok[k][j][i  ] + rhok[k+1][j][i+1] - rhok[k][j][i+1]) * 0.5;
                }
                else
                {
                    dbdz = (rhok[k+1][j][i] - rhok[k-1][j][i] + rhok[k+1][j][i+1] - rhok[k-1][j][i+1]) * 0.25;
                }

                db[k][j][i].x = g * coor[k][j][i].z * (dbdc * g11_i + dbde *  g12_i + dbdz * g13_i ) * iaj[k][j][i];

                // pressure gradient in the j-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (
                        i == mx-2 &&
                        (
                            j==my-2
                        )
                    )
                )
                {
                    dbdc = (rhok[k][j  ][i] - rhok[k][j  ][i-1] + rhok[k][j+1][i] - rhok[k][j+1][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (
                        i == 1 &&
                        (
                            j==my-2
                        )
                    )
                )
                {
                    dbdc = (rhok[k][j  ][i+1] - rhok[k][j  ][i] + rhok[k][j+1][i+1] - rhok[k][j+1][i]) * 0.5;
                }
                else
                {
                    dbdc = (rhok[k][j  ][i+1] - rhok[k][j  ][i-1] + rhok[k][j+1][i+1] - rhok[k][j+1][i-1]) * 0.25;
                }

                if( j==my-2 && mesh->jj_periodic)
                {
                    dbde = rhok[k][my+1][i] - rhok[k][j][i];
                }
                else
                {
                    dbde = rhok[k][j+1][i] - rhok[k][j][i];
                }

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (
                        k == mz-2 &&
                        (
                            j== my-2
                        )
                    )
                )
                {
                    dbdz = (rhok[k][j  ][i] - rhok[k-1][j  ][i] + rhok[k][j+1][i] - rhok[k-1][j+1][i]) * 0.5;
                }
                else if
                (
                    // k-left boundary -> use upwind  only at the corner faces
                    (
                        k == 1 &&
                        (
                            j== my-2
                        )
                    )
                )
                {
                    dbdz = (rhok[k+1][j  ][i] - rhok[k][j  ][i] + rhok[k+1][j+1][i] - rhok[k][j+1][i]) * 0.5;
                }
                else
                {
                    dbdz = (rhok[k+1][j  ][i] - rhok[k-1][j  ][i] + rhok[k+1][j+1][i] - rhok[k-1][j+1][i]) * 0.25;
                }

                db[k][j][i].y = g * coor[k][j][i].z * (dbdc * g21_j + dbde * g22_j + dbdz * g23_j ) * jaj[k][j][i];

                // pressure gradient in the k-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (
                        i == mx - 2 &&
                        (
                            k==mz-2
                        )
                    )
                )
                {
                    dbdc = (rhok[k  ][j][i] - rhok[k  ][j][i-1] + rhok[k+1][j][i] - rhok[k+1][j][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (
                        i == 1 &&
                        (
                            k == mz - 2
                        )
                    )
                )
                {
                    dbdc = (rhok[k  ][j][i+1] - rhok[k  ][j][i] + rhok[k+1][j][i+1] - rhok[k+1][j][i]) * 0.5;
                }
                else
                {
                    dbdc = (rhok[k  ][j][i+1] - rhok[k  ][j][i-1] + rhok[k+1][j][i+1] - rhok[k+1][j][i-1]) * 0.25;
                }

                if
                (
                    // j-right boundary -> use upwind  only at the corner faces
                    (
                        j == my - 2 &&
                        (
                            k==mz-2
                        )
                    )
                )
                {
                    dbde = (rhok[k  ][j][i] - rhok[k  ][j-1][i] + rhok[k+1][j][i] - rhok[k+1][j-1][i]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (
                        j == 1 &&
                        (
                            k==mz-2
                        )
                    )
                )
                {
                    dbde = (rhok[k  ][j+1][i] - rhok[k  ][j][i] + rhok[k+1][j+1][i] - rhok[k+1][j][i]) * 0.5;
                }
                else
                {
                    dbde = (rhok[k  ][j+1][i] - rhok[k  ][j-1][i] + rhok[k+1][j+1][i] - rhok[k+1][j-1][i]) * 0.25;
                }

                if( k==mz-2 && mesh->kk_periodic)
                {
                    dbdz = rhok[mz+1][j][i] - rhok[k][j][i];
                }
                else
                {
                    dbdz = (rhok[k+1][j][i] - rhok[k][j][i]);
                }

                db[k][j][i].z = g * coor[k][j][i].z * (dbdc * g31_k + dbde * g32_k + dbdz * g33_k ) * kaj[k][j][i];

                // periodic: set to zero only on left boundaries since the contrav. velocity is not solved there
                // non-periodic: set to zero also on right boundaries since the contrav. velocity is not solved there
                if
                (
                    i==0 || (!mesh->i_periodic && !mesh->ii_periodic && i==mx-2)
                )
                {
                    db[k][j][i].x = 0;
                }
                if
                (
                    j==0 || (!mesh->j_periodic && !mesh->jj_periodic && j==my-2)
                )
                {
                    db[k][j][i].y = 0;
                }
                if
                (
                    k==0 || (!mesh->k_periodic && !mesh->kk_periodic && k==mz-2)
                )
                {
                    db[k][j][i].z = 0;
                }
            }
        }
    }

    DMDAVecRestoreArray(da, teqn->lRhoK,  &rhok);

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
    DMDAVecRestoreArray(da,  mesh->lIAj,  &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,  &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,  &kaj);

    DMDAVecRestoreArray(da,  teqn->lTmprt,      &tmprt);
    DMDAVecRestoreArray(fda, teqn->ghGradRhok,  &db);

    DMDAVecRestoreArray(fda, Coor, &coor);

    return(0);

}

//***************************************************************************************************************//
PetscErrorCode TmprtPredictor(teqn_ *teqn)
{
    // Explicit forward-Euler conv-only predictor: T* = T^n + dt * N(T^n)
    // Advances Tmprt/lTmprt from T^n to T* so that Buoyancy() in SolveUEqn
    // sees a forward-extrapolated temperature field.
    // Tmprt_o retains T^n so TmprtRestoreFromOld() can undo this before SolveTEqn.

    mesh_  *mesh  = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;

    PetscReal ts, te;
    PetscTime(&ts);
    //PetscPrintf(mesh->MESH_COMM, "Predicted potential temperature in ");

    VecSet(teqn->Rhs, 0.0);

    // conv-only explicit RHS: avoids explicit diffusion stability limit (dt << dx^2/kappa_eff)
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");
    FormT(teqn, teqn->Rhs, 1.0, 1);
    resetNonResolvedCellCentersScalar(mesh, teqn->Rhs);

    // T* = T^n + 0.5*dt * N(T^n)  — half-step: midpoint estimate consistent with CN U integration
    VecAXPY(teqn->Tmprt, 0.5 * clock->dt, teqn->Rhs);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // fix periodic ghost cells: ghost slots at i=0/mx-1 were left at T^n, not T*
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    PetscTime(&te);
    //PetscPrintf(mesh->MESH_COMM, "%f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TmprtRestoreFromOld(teqn_ *teqn)
{
    // Restore Tmprt/lTmprt to T^n (saved in Tmprt_o) after TmprtPredictor has
    // temporarily advanced them to T* for buoyancy in SolveUEqn.
    // Must be called before SolveTEqn so the T solve starts from T^n.

    mesh_ *mesh = teqn->access->mesh;

    VecCopy(teqn->Tmprt_o, teqn->Tmprt);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    return(0);
}
