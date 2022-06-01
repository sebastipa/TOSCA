//! \file  ueqn.c
//! \brief Contains U equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorU(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeUEqn(ueqn_ *ueqn)
{
    // set pointer to mesh
    mesh_ *mesh = ueqn->access->mesh;

    // input file
    PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

    // default parameters
    ueqn->inviscid          = 0;
    ueqn->centralDiv        = 0;
    ueqn->centralUpwindDiv  = 0;
    ueqn->centralUpwindWDiv = 0;
    ueqn->weno3Div          = 0;
    ueqn->quickDiv          = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL,  "-inviscid", &(ueqn->inviscid),   PETSC_NULL);

    // read time discretization scheme
    readDictWord("control.dat", "-dUdtScheme", &(ueqn->ddtScheme));

    // read divergence scheme
    readDictWord("control.dat", "-divScheme", &(ueqn->divScheme));

    if(     ueqn->divScheme == "centralUpwind")  ueqn->centralUpwindDiv  = 1;
    else if(ueqn->divScheme == "centralUpwindW") ueqn->centralUpwindWDiv = 1;
    else if(ueqn->divScheme == "central")        ueqn->centralDiv        = 1;
    else if(ueqn->divScheme == "weno3")          ueqn->weno3Div          = 1;
    else if(ueqn->divScheme == "quickDiv")       ueqn->quickDiv          = 1;

    VecDuplicate(mesh->Cent, &(ueqn->Utmp));      VecSet(ueqn->Utmp,    0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Rhs));       VecSet(ueqn->Rhs,     0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Rhs_o));     VecSet(ueqn->Rhs_o,   0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucont));     VecSet(ueqn->Ucont,   0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucont_o));   VecSet(ueqn->Ucont_o, 0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucat));      VecSet(ueqn->Ucat,    0.0);
    VecDuplicate(mesh->Cent, &(ueqn->dP));        VecSet(ueqn->dP,      0.0);
    VecDuplicate(mesh->Cent, &(ueqn->sourceU));   VecSet(ueqn->sourceU, 0.0);
    VecDuplicate(mesh->Cent, &(ueqn->gCont));     VecSet(ueqn->gCont,   0.0);

    VecDuplicate(mesh->lCent, &(ueqn->lUcont));   VecSet(ueqn->lUcont,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lUcat));    VecSet(ueqn->lUcat,   0.0);

    VecDuplicate(mesh->lCent, &(ueqn->lFp));      VecSet(ueqn->lFp,     0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv1));    VecSet(ueqn->lDiv1,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv2));    VecSet(ueqn->lDiv2,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv3));    VecSet(ueqn->lDiv3,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc1));   VecSet(ueqn->lVisc1,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc2));   VecSet(ueqn->lVisc2,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc3));   VecSet(ueqn->lVisc3,  0.0);

    VecDuplicate(mesh->lAj, &(ueqn->lUstar));     VecSet(ueqn->lUstar,  0.0);

    if(ueqn->ddtScheme == "backwardEuler")
    {
        ueqn->relExitTol        = 1e-30;
        ueqn->absExitTol        = 1e-5;

        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolU",  &(ueqn->relExitTol), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolU",  &(ueqn->absExitTol), PETSC_NULL);

        // create the SNES solver
        SNESCreate(mesh->MESH_COMM, &(ueqn->snesU));
        SNESMonitorSet(ueqn->snesU, SNESMonitorU, (void*)&(mesh->MESH_COMM), PETSC_NULL);

        // set the SNES evaluating function
        SNESSetFunction(ueqn->snesU, ueqn->Rhs, UeqnSNES, (void *)ueqn);

        // create jacobian matrix
        MatCreateSNESMF(ueqn->snesU, &(ueqn->JU));
        SNESSetJacobian(ueqn->snesU, ueqn->JU, ueqn->JU, MatMFFDComputeJacobian, (void *)ueqn);

        // set SNES solver type
        SNESSetType(ueqn->snesU, SNESNEWTONTR);          //SNESTR
        //SNESSetType(ueqn->snesU, SNESNEWTONLS);        //SNESLS is better for stiff PDEs such as the one including IB but slower

        // set SNES solve and step failures
        SNESSetMaxLinearSolveFailures(ueqn->snesU,10000);
        SNESSetMaxNonlinearStepFailures(ueqn->snesU,10000);
        SNESKSPSetUseEW(ueqn->snesU, PETSC_TRUE);

        // set SNES Krylov Sub-Space parameters
        SNESKSPSetParametersEW(ueqn->snesU,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

        // SNES tolerances
        // 2nd arg: absolute tolerance
        // 3rd arg: relative tolerance
        // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
        // 5th arg: maximum number of iterations
        // 6th arg: maximum function evaluations
        SNESSetTolerances(ueqn->snesU, ueqn->absExitTol, 1e-30, 1e-30, 20, 1000);

        SNESGetKSP(ueqn->snesU, &(ueqn->ksp));
        KSPGetPC(ueqn->ksp,&(ueqn->pc));

        // set KSP solver type
        KSPSetType(ueqn->ksp, KSPGMRES);

        //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);    //2009.09.22 poor performance
        //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);      //2009.09.22

        //KSPFischerGuess itg;
        //KSPFischerGuessCreate(ksp,1,100,&itg);
        //KSPSetFischerGuess(ksp, itg);                  //2009.09.22

        //KSPGMRESSetPreAllocateVectors(ksp);            --> crazy thing consumes memory

        PCSetType(ueqn->pc, PCNONE);
        PetscReal rtol=ueqn->relExitTol, atol=ueqn->absExitTol, dtol=PETSC_DEFAULT;
        KSPSetTolerances(ueqn->ksp, rtol, atol, dtol, 1000);
    }
    else if(ueqn->ddtScheme == "forwardEuler" || ueqn->ddtScheme == "rungeKutta4")
    {
        // explicit schemes
    }
    else
    {
         char error[512];
         sprintf(error, "unknown ddtScheme %s for U equation, available schemes are\n    1. backwardEuler\n    2. forwardEuler\n    3. rungeKutta4", ueqn->ddtScheme.c_str());
         fatalErrorInFunction("InitializeUEqn", error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateFluxLimiter(ueqn_ *ueqn)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucont, ***limiter;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(j!=0 && k!=0)
                {
                    PetscInt    iL, iR;
                    getFace2Face3StencilCsi(mesh, i, mx, &iL, &iR);

                    limiter[k][j][i].x
                    =
                    vanLeer
                    (
                        ucont[k][j][iL].x,
                        ucont[k][j][i].x,
                        ucont[k][j][iR].x
                    );
                }

                if(i!=0 && k!=0)
                {
                    PetscInt    jL, jR;
                    getFace2Face3StencilEta(mesh, j, my, &jL, &jR);

                    limiter[k][j][i].y
                    =
                    vanLeer
                    (
                        ucont[k][jL][i].y,
                        ucont[k][j][i].y,
                        ucont[k][jR][i].y
                    );
                }

                if(j!=0 && i!=0)
                {
                    PetscInt    kL, kR;
                    getFace2Face3StencilZet(mesh, k, mz, &kL, &kR);

                    limiter[k][j][i].z
                    =
                    vanLeer
                    (
                        ucont[kL][j][i].z,
                        ucont[k][j][i].z,
                        ucont[kR][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CorrectSourceTerms(ueqn_ *ueqn, PetscInt print)
{
    mesh_         *mesh = ueqn->access->mesh;
    clock_        *clock= ueqn->access->clock;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucat,  ***source;
    PetscReal     ***nvert, ***aj;
    Cmpnts        ***cent;

    Cmpnts        uDes;
    Cmpnts        s;

    PetscReal     meanCellVolume = 0, sumVolume = 0;
    PetscReal     rASum = 0, rAMean = 0;

    abl_          *abl  = ueqn->access->abl;

    PetscReal     relax = abl->relax;
    PetscReal     alpha = abl->alpha;
    PetscReal     T     = abl->timeWindow;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->sourceU, &source);

    // compute source terms using the velocity controller
    if(abl->controllerType=="write")
    {
        // set the wanted velocity
        uDes.x = abl->uRef;
        uDes.y = 0.0;
        uDes.z = 0.0;

        DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
        DMDAVecGetArray(da, mesh->lAj, &aj);

        // find the first two closest levels
        PetscReal nLevels = my-2;

        Cmpnts luMean1; luMean1.x = 0.0; luMean1.y = 0.0; luMean1.z = 0.0;
        Cmpnts luMean2; luMean2.x = 0.0; luMean2.y = 0.0; luMean2.z = 0.0;
        Cmpnts guMean1; guMean1.x = 0.0; guMean1.y = 0.0; guMean1.z = 0.0;
        Cmpnts guMean2; guMean2.x = 0.0; guMean2.y = 0.0; guMean2.z = 0.0;

        for (k=zs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    if(j==abl->closestLabels[0])
                    {
                        luMean1.x += ucat[k][j][i].x / aj[k][j][i];
                        luMean1.y += ucat[k][j][i].y / aj[k][j][i];
                        luMean1.z += ucat[k][j][i].z / aj[k][j][i];
                    }
                    else if(j==abl->closestLabels[1])
                    {
                        luMean2.x += ucat[k][j][i].x / aj[k][j][i];
                        luMean2.y += ucat[k][j][i].y / aj[k][j][i];
                        luMean2.z += ucat[k][j][i].z / aj[k][j][i];
                    }
                }
            }
        }

        DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
        DMDAVecRestoreArray(da, mesh->lAj, &aj);

        MPI_Allreduce(&luMean1, &guMean1, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&luMean2, &guMean2, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        guMean1.x = guMean1.x / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean1.y = guMean1.y / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean1.z = guMean1.z / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean2.x = guMean2.x / abl->totVolPerLevel[abl->closestLabels[1]-1];
        guMean2.y = guMean2.y / abl->totVolPerLevel[abl->closestLabels[1]-1];
        guMean2.z = guMean2.z / abl->totVolPerLevel[abl->closestLabels[1]-1];

        Cmpnts uMean;

        uMean.x = guMean1.x * abl->levelWeights[0] + guMean2.x * abl->levelWeights[1];
        uMean.y = guMean1.y * abl->levelWeights[0] + guMean2.y * abl->levelWeights[1];
        uMean.z = guMean1.z * abl->levelWeights[0] + guMean2.z * abl->levelWeights[1];

        if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: wind height is %lf m, h1 = %lf m, h2 = %lf m\n", abl->hRef, abl->cellLevels[abl->closestLabels[0]-1], abl->cellLevels[abl->closestLabels[1]-1]);

        PetscReal magUMean = std::sqrt(uMean.x*uMean.x + uMean.y*uMean.y + uMean.z*uMean.z);
        PetscReal magUDes  = std::sqrt(uDes.x*uDes.x + uDes.y*uDes.y + uDes.z*uDes.z);

        // compute the error w.r.t reference
        Cmpnts error;
        error.x = (uDes.x - uMean.x);
        error.y = (uDes.y - uMean.y);
        error.z = (uDes.z - uMean.z);

        // compute relative error for print statement
        Cmpnts relError;
        relError.x = error.x / magUDes;
        relError.y = error.y / magUDes;
        relError.z = error.z / magUDes;

        if(print) PetscPrintf(mesh->MESH_COMM, "                         avg mag U at hRef = %lf m/s, U desired = %lf m/s\n", magUMean, magUDes);
        if(print) PetscPrintf(mesh->MESH_COMM, "                         U error in perc. of desired: (%.3f %.3f %.3f)\n", relError.x*100,relError.y*100, relError.z*100);

        // cumulate the error (integral part of the controller)
        abl->cumulatedSource.x = std::exp(-clock->dt / T) * abl->cumulatedSource.x + (clock->dt / T) * error.x;
        abl->cumulatedSource.y = std::exp(-clock->dt / T) * abl->cumulatedSource.y + (clock->dt / T) * error.y;
        abl->cumulatedSource.z = 0.0;

        // compute the uniform source terms (PI controller with adjustable gains)
        s.x = relax * (alpha * error.x + (1.0-alpha) * abl->cumulatedSource.x) ;
        s.y = relax * (alpha * error.y + (1.0-alpha) * abl->cumulatedSource.y) ;
        s.z = 0.0;

        // write the uniform source term
        if(!rank)
        {
            if(clock->it == clock->itStart)
            {
                errno = 0;
                PetscInt dirRes = mkdir("./postProcessing", 0777);
                if(dirRes != 0 && errno != EEXIST)
                {
                   char error[512];
                    sprintf(error, "could not create postProcessing directory\n");
                    fatalErrorInFunction("correctSourceTerm",  error);
                }
                else
                {
                    unlink("postProcessing/momentumSource");
                }
            }

            FILE *fp = fopen("postProcessing/momentumSource", "a");
            if(!fp)
            {
               char error[512];
                sprintf(error, "cannot open file postProcessing/momentumSource\n");
                fatalErrorInFunction("correctSourceTerm",  error);
            }
            else
            {
                PetscInt width = -15;

                if(clock->it == clock->itStart)
                {
                    word w1 = "time";
                    word w2 = "source x [m/s2]";
                    word w3 = "source y [m/s2]";
                    word w4 = "source z [m/s2]";

                    PetscFPrintf(mesh->MESH_COMM, fp, "%*s\t%*s\t%*s\t%*s\n", width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str());
                }

                PetscFPrintf(mesh->MESH_COMM, fp, "%*.2f\t%*.5e\t%*.5e%*.5e\t\n", width, clock->time, width, s.x, width, s.y, width, s.z);

                fclose(fp);
            }
        }
    }
    else if(abl->controllerType=="read")
    {
        // find the two closest times in the pre-computed sources
        // get the 2 time values closest to the current time
        PetscInt idx_1 = abl->currentCloseIdx;
        PetscInt idx_2 = abl->currentCloseIdx + 1;

        PetscInt lwrBound = 0;
        PetscInt uprBound = abl->nSourceTimes;

        // if past first iteration do the search on a subset to speed up the process
        if(clock->it > clock->itStart)
        {
            lwrBound = PetscMax(0, (abl->currentCloseIdx - 50));
            uprBound = PetscMin(abl->nSourceTimes, (abl->currentCloseIdx + 50));
        }

        // build error vector for the time search
        PetscReal  diff[abl->nSourceTimes];

        for(PetscInt i=lwrBound; i<uprBound; i++)
        {
            diff[i] = fabs(abl->preCompSources[i][0] - clock->time);
        }

        // find the two closest times
        for(PetscInt i=lwrBound; i<uprBound; i++)
        {
            if(diff[i] < diff[idx_1])
            {
                idx_2 = idx_1;
                idx_1 = i;
            }
            if(diff[i] < diff[idx_2] && i != idx_1)
            {
                idx_2 = i;
            }
        }

        // always put the lower time at idx_1 and higher at idx_2
        if(abl->preCompSources[idx_2][0] < abl->preCompSources[idx_1][0])
        {
            PetscInt idx_tmp = idx_2;
            idx_2 = idx_1;
            idx_1 = idx_tmp;
        }

        // find interpolation weights
        PetscReal idx = (idx_2 - idx_1) / (abl->preCompSources[idx_2][0] - abl->preCompSources[idx_1][0]) * (clock->time - abl->preCompSources[idx_1][0]) + idx_1;
        PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
        PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

        if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading sources\n", w1 * abl->preCompSources[idx_1][0] + w2 * abl->preCompSources[idx_2][0]);
        if(print) PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
        if(print) PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->preCompSources[idx_1][0], abl->preCompSources[idx_2][0]);

        // reset the closest index for nex iteration
        abl->currentCloseIdx = idx_1;

        // get also the dt at the time the source was calculated
        double dtSource1 = std::max((abl->preCompSources[idx_1+1][0] - abl->preCompSources[idx_1][0]), 1e-5);
        double dtSource2 = std::max((abl->preCompSources[idx_2+1][0] - abl->preCompSources[idx_2][0]), 1e-5);

        // do not scale with time step
        // s.x = w1 * abl->preCompSources[idx_1][1] + w2 * abl->preCompSources[idx_2][1];
        // s.y = w1 * abl->preCompSources[idx_1][2] + w2 * abl->preCompSources[idx_2][2];
        // s.z = w1 * abl->preCompSources[idx_1][3] + w2 * abl->preCompSources[idx_2][3];

        // scale with time step
        s.x = (w1 * abl->preCompSources[idx_1][1] / dtSource1 + w2 * abl->preCompSources[idx_2][1] / dtSource2) * clock->dt;
        s.y = (w1 * abl->preCompSources[idx_1][2] / dtSource1 + w2 * abl->preCompSources[idx_2][2] / dtSource2) * clock->dt;
        s.z = (w1 * abl->preCompSources[idx_1][3] / dtSource1 + w2 * abl->preCompSources[idx_2][3] / dtSource2) * clock->dt;


    }
    else if(abl->controllerType=="average")
    {
        // PetscPrintf(mesh->MESH_COMM, "Correcting source terms using source average from time %lf\n", abl->sourceAvgStartTime);

        // do not scale with dt
        // s.x = abl->cumulatedSource.x;
        // s.y = abl->cumulatedSource.y;
        // s.z = abl->cumulatedSource.z;

        // scale with dt
        PetscReal dtScale = clock->dt / abl->avgTimeStep;

        s.x = abl->cumulatedSource.x * dtScale;
        s.y = abl->cumulatedSource.y * dtScale;
        s.z = abl->cumulatedSource.z * dtScale;


    }

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // update the source term
    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(cent[k][j][i].z <= abl->controllerHeight)
                {
                    source[k][j][i].x = s.x;
                    source[k][j][i].y = s.y;
                    source[k][j][i].z = s.z;
                }
                else
                {
                    source[k][j][i].x = 0.0;
                    source[k][j][i].y = 0.0;
                    source[k][j][i].z = 0.0;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->sourceU, &source);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode sourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***rhs, ***source, ***cent;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // damping viscosity for fringe region exclusion
    double nudI;
    double nudJ;
    double nudK;

    // fringe region parameters (set only if active)
    double xS;
    double xE;
    double xD;

    if(ueqn->access->flags->isXDampingActive)
    {
        xS     = ueqn->access->abl->xDampingStart;
        xE     = ueqn->access->abl->xDampingEnd;
        xD     = ueqn->access->abl->xDampingDelta;
    }
    else
    {
        nudI = nudJ = nudK = 1.0;
    }

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(fda, Rhs,  &rhs);

    DMDAVecGetArray(fda, ueqn->sourceU, &source);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if
                (
                    isFluidIFace(k, j, i, i+1, nvert)
                )
                {
                    if(ueqn->access->flags->isXDampingActive)
                    {
                        // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                        double x     = cent[k][j][i].x;
                        double xi    = cent[k][j][i+1].x;
                        double xj    = cent[k][j+1][i].x;
                        double xk    = cent[k+1][j][i].x;

                        // compute Stipa viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                        double nud_x   = viscStipa(xS, xE, xD, x);
                        double nudi_x  = viscStipa(xS, xE, xD, xi);
                        double nudj_x  = viscStipa(xS, xE, xD, xj);
                        double nudk_x  = viscStipa(xS, xE, xD, xk);

                        nudI           = central(nud_x, nudi_x);
                        nudJ           = central(nud_x, nudj_x);
                        nudK           = central(nud_x, nudk_x);
                    }

                    rhs[k][j][i].x
                    +=
                    scale * nudI *
                    (
                          source[k][j][i].x * icsi[k][j][i].x +
                          source[k][j][i].y * icsi[k][j][i].y +
                          source[k][j][i].z * icsi[k][j][i].z
                    );
                }

                if
                (
                        isFluidJFace(k, j, i, j+1, nvert)
                )
                {
                    rhs[k][j][i].y
                    +=
                    scale * nudJ *
                    (
                          source[k][j][i].x * jeta[k][j][i].x +
                          source[k][j][i].y * jeta[k][j][i].y +
                          source[k][j][i].z * jeta[k][j][i].z
                    );
                }

                if
                (
                        isFluidKFace(k, j, i, k+1, nvert)
                )
                {
                    rhs[k][j][i].z
                    +=
                    scale * nudK *
                    (
                          source[k][j][i].x * kzet[k][j][i].x +
                          source[k][j][i].y * kzet[k][j][i].y +
                          source[k][j][i].z * kzet[k][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->sourceU, &source);

    DMDAVecRestoreArray(fda, Rhs,  &rhs);

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode correctDampingSources(ueqn_ *ueqn)
{
    mesh_         *mesh = ueqn->access->mesh;
    abl_          *abl  = ueqn->access->abl;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***cent, ***ucat, ***ucatP;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // update uBar for xDampingLayer (read from inflow data)
    if(ueqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2)
        {
            // get inflow database info pointer
            inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

            // variables to recover lapse rate above data end
            PetscReal ldataHeight = 0, gdataHeight = 0;

            // define local uBar and tBar vectors
            std::vector<std::vector<Cmpnts>> luBarInstX(my);
            std::vector<std::vector<PetscReal>> ltBarInstX(my);

            // set it to zero
            for(j=0; j<my; j++)
            {
                luBarInstX[j].resize(mx);
                ltBarInstX[j].resize(mx);

                for(i=0; i<mx; i++)
                {
                    luBarInstX[j][i].x = 0.0;
                    luBarInstX[j][i].y = 0.0;
                    luBarInstX[j][i].z = 0.0;

                    ltBarInstX[j][i]   = 0.0;
                }
            }

            DMDAVecGetArray(fda, mesh->lCent, &cent);

            // read inflow if necessary
            {
                readInflowU(ifPtr, ueqn->access->clock);
            }

            // read T data from database
            if(ueqn->access->flags->isTeqnActive)
            {
                readInflowT(ifPtr, ueqn->access->clock);

                // compute hight at which temperature inflow data ends

                // type 1: inflow and actual meshes have same cell dimensions
                if (ifPtr->typeT == 1)
                {
                    PetscInt lcount = 0, gcount = 0;

                    // make sure this processor can access data
                    if(lys <= ifPtr->n1*ifPtr->prds1 && ifPtr->n1*ifPtr->prds1 <= lye)
                    {
                        // compute end of data height (j is vertical direction)
                        i = std::floor(0.5*(lxe-lxs) + lxs);
                        k = std::floor(0.5*(lze-lzs) + lzs);
                        ldataHeight = cent[k][ifPtr->n1*ifPtr->prds1][i].z;
                        lcount      = 1;
                    }

                    MPI_Allreduce(&ldataHeight, &gdataHeight, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                    MPI_Allreduce(&lcount, &gcount, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                    gdataHeight = gdataHeight / gcount;
                }
                // type 2: inflow and actual meshes are different, use inflow width
                else if (ifPtr->typeT == 2)
                {
                    gdataHeight = ifPtr->n1 * ifPtr->prds1 * ifPtr->width1;
                }
            }

            // update uBarInstX at this time step for this processor. These fields are defined at
            // k-face centers, so the indexing is equal to cell centers.

            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // unsteady mapped
                    if (ifPtr->typeU == 1)
                    {
                        // periodize inflow according to input

                        // compute period fraction (handle index = n case)
                        PetscInt jif = j % ifPtr->n1 == 0 ? ifPtr->n1 : j % ifPtr->n1;
                        PetscInt iif = i % ifPtr->n2 == 0 ? ifPtr->n2 : i % ifPtr->n2;

                        // index is less than nPrds times inflow points: have data
                        if
                        (
                            j<=ifPtr->n1*ifPtr->prds1 &&
                            i<=ifPtr->n2*ifPtr->prds2
                        )
                        {
                            luBarInstX[j][i].x = ifPtr->ucat_plane[jif][iif].x;
                            luBarInstX[j][i].y = ifPtr->ucat_plane[jif][iif].y;
                            luBarInstX[j][i].z = ifPtr->ucat_plane[jif][iif].z;
                        }
                        // index is more than nPrds times inflow points: extrapolate
                        else
                        {
                            // extrapolate along j
                            if(j>ifPtr->n1*ifPtr->prds1) jif = ifPtr->n1;

                            // extrapolate along i
                            if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                            luBarInstX[j][i].x = ifPtr->ucat_plane[jif][iif].x;
                            luBarInstX[j][i].y = ifPtr->ucat_plane[jif][iif].y;
                            luBarInstX[j][i].z = ifPtr->ucat_plane[jif][iif].z;
                        }
                    }
                    // unsteady mapped interpolated
                    else if (ifPtr->typeU == 2)
                    {
                        luBarInstX[j][i].x
                        =
                        ifPtr->inflowWeights[j][i][0] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].x +
                        ifPtr->inflowWeights[j][i][1] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].x +
                        ifPtr->inflowWeights[j][i][2] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].x +
                        ifPtr->inflowWeights[j][i][3] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].x;

                        luBarInstX[j][i].y
                        =
                        ifPtr->inflowWeights[j][i][0] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].y +
                        ifPtr->inflowWeights[j][i][1] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].y +
                        ifPtr->inflowWeights[j][i][2] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].y +
                        ifPtr->inflowWeights[j][i][3] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].y;

                        luBarInstX[j][i].z
                        =
                        ifPtr->inflowWeights[j][i][0] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].z +
                        ifPtr->inflowWeights[j][i][1] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].z +
                        ifPtr->inflowWeights[j][i][2] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].z +
                        ifPtr->inflowWeights[j][i][3] *
                        ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].z;
                    }

                    // set k to the starting value of this processor. It is needed to
                    // get the z coordinate. This is only valid for cartesian meshes
                    PetscInt k_idx = lzs;

                    if(ueqn->access->flags->isTeqnActive)
                    {
                        // periodized mapped inflow
                        if (ifPtr->typeT == 1)
                        {
                            // periodize inflow according to input

                            // compute period fraction (handle index = n case)
                            PetscInt jif = j % ifPtr->n1 == 0 ? ifPtr->n1 : j % ifPtr->n1;
                            PetscInt iif = i % ifPtr->n2 == 0 ? ifPtr->n2 : i % ifPtr->n2;

                            // index is less than nPrds times inflow points: have data
                            if
                            (
                                j<=ifPtr->n1*ifPtr->prds1 &&
                                i<=ifPtr->n2*ifPtr->prds2
                            )
                            {
                                ltBarInstX[j][i] = ifPtr->t_plane[jif][iif];
                            }
                            // index is more than nPrds times inflow points: apply lapse rate
                            else
                            {
                                PetscReal delta = 0;

                                // extrapolate along j
                                if(j>ifPtr->n1*ifPtr->prds1)
                                {
                                    jif   = ifPtr->n1;
                                    delta = cent[k_idx][j][i].z - gdataHeight;
                                }

                                // extrapolate along i
                                if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                                ltBarInstX[j][i] = ifPtr->t_plane[jif][iif] + delta * abl->gTop;
                            }
                        }

                        // interpolated periodized mapped inflow
                        else if (ifPtr->typeT == 2)
                        {
                            PetscReal delta = PetscMax(0.0, cent[k_idx][j][i].z - gdataHeight);

                            ltBarInstX[j][i]
                            =
                            ifPtr->inflowWeights[j][i][0] *
                            ifPtr->t_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i] +
                            ifPtr->inflowWeights[j][i][1] *
                            ifPtr->t_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i] +
                            ifPtr->inflowWeights[j][i][2] *
                            ifPtr->t_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i] +
                            ifPtr->inflowWeights[j][i][3] *
                            ifPtr->t_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i] +
                            delta * abl->gTop;
                        }
                    }
                }
            }

            // scatter to all processors
            for(j=0; j<my; j++)
            {
                MPI_Allreduce(&luBarInstX[j][0], &abl->uBarInstX[j][0], 3*mx, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

                if(ueqn->access->flags->isTeqnActive)
                {
                    MPI_Allreduce(&ltBarInstX[j][0], &abl->tBarInstX[j][0],   mx, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                }
            }

            // divide by number of processors in the line
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    abl->uBarInstX[j][i].x /= (PetscReal)abl->nProcsKLine[j][i];
                    abl->uBarInstX[j][i].y /= (PetscReal)abl->nProcsKLine[j][i];
                    abl->uBarInstX[j][i].z /= (PetscReal)abl->nProcsKLine[j][i];

                    if(ueqn->access->flags->isTeqnActive)
                    {
                        abl->tBarInstX[j][i]   /= (PetscReal)abl->nProcsKLine[j][i];
                    }
                }
            }

            // clear local vectors
            for(j=0; j<my; j++)
            {
                std::vector<Cmpnts> ().swap(luBarInstX[j]);
                std::vector<PetscReal> ().swap(ltBarInstX[j]);
            }

            // restore cell centers array
            DMDAVecRestoreArray(fda, mesh->lCent, &cent);

            // periodize in i and j directions (we follow cell indexing)
            for (j=0; j<my; j++)
            {
                for (i=0; i<mx; i++)
                {
                    if(i==0)
                    {
                        abl->uBarInstX[j][i].x = abl->uBarInstX[j][mx-2].x;
                        abl->uBarInstX[j][i].y = abl->uBarInstX[j][mx-2].y;
                        abl->uBarInstX[j][i].z = abl->uBarInstX[j][mx-2].z;
                        if(ueqn->access->flags->isTeqnActive) abl->tBarInstX[j][i]   = abl->tBarInstX[j][mx-2]  ;
                    }
                    if(i==mx-1)
                    {
                        abl->uBarInstX[j][i].x = abl->uBarInstX[j][1].x;
                        abl->uBarInstX[j][i].y = abl->uBarInstX[j][1].y;
                        abl->uBarInstX[j][i].z = abl->uBarInstX[j][1].z;
                        if(ueqn->access->flags->isTeqnActive) abl->tBarInstX[j][i]   = abl->tBarInstX[j][1]  ;
                    }
                    if(j==0)
                    {
                        abl->uBarInstX[j][i].x = abl->uBarInstX[my-2][i].x;
                        abl->uBarInstX[j][i].y = abl->uBarInstX[my-2][i].y;
                        abl->uBarInstX[j][i].z = abl->uBarInstX[my-2][i].z;
                        if(ueqn->access->flags->isTeqnActive) abl->tBarInstX[j][i]   = abl->tBarInstX[my-2][i]  ;
                    }
                    if(j==my-1)
                    {
                        abl->uBarInstX[j][i].x = abl->uBarInstX[1][i].x;
                        abl->uBarInstX[j][i].y = abl->uBarInstX[1][i].y;
                        abl->uBarInstX[j][i].z = abl->uBarInstX[1][i].z;
                        if(ueqn->access->flags->isTeqnActive) abl->tBarInstX[j][i]   = abl->tBarInstX[1][i]  ;
                    }
                }
            }
        }
        else if(abl->xFringeUBarSelectionType == 3)
        {
            if(abl->xDampingControlType == "alphaOptimized")
            {
                clock_ *clock = ueqn->access->clock;

                // fringe region starting cell lines
                PetscReal lsumStart1   = 0.0, gsumStart1   = 0.0,
                          lsumStart2   = 0.0, gsumStart2   = 0.0;
                PetscInt  lcountStart1 = 0,   gcountStart1 = 0,
                          lcountStart2 = 0,   gcountStart2 = 0;

                // fringe region ending cell lines
                PetscReal lsumEnd1   = 0.0, gsumEnd1   = 0.0,
                          lsumEnd2   = 0.0, gsumEnd2   = 0.0;
                PetscInt  lcountEnd1 = 0,   gcountEnd1 = 0,
                          lcountEnd2 = 0,   gcountEnd2 = 0;

                // precursor average desired velocity
                PetscReal lsum1   = 0.0, gsum1   = 0.0,
                          lsum2   = 0.0, gsum2   = 0.0;
                PetscInt  lcount1 = 0,   gcount1 = 0,
                          lcount2 = 0,   gcount2 = 0;

                precursor_ *precursor = ueqn->access->abl->precursor;

                if(precursor->thisProcessorInFringe)
                {
                    DMDAVecGetArray(fda, precursor->domain->ueqn->lUcat, &ucatP);
                }

                DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
                DMDAVecGetArray(fda, mesh->lCent, &cent);

                for (k=zs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            // velocity values at fringe start in successor domain (two levels to interpolate)
                            if(j == abl->closestLabels[0] && k == precursor->map.kStart + 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumStart1 += ucat[k][j][i].y;
                                lcountStart1 ++;
                            }
                            else if(j==abl->closestLabels[1] && k == precursor->map.kStart + 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumStart2 += ucat[k][j][i].y;
                                lcountStart2 ++;
                            }

                            // velocity values at fringe end in successor domain (two levels to interpolate)
                            else if(j == abl->closestLabels[0] && k == precursor->map.kEnd - 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumEnd1 += ucat[k][j][i].y;
                                lcountEnd1 ++;
                            }
                            else if(j==abl->closestLabels[1] && k == precursor->map.kEnd - 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumEnd2 += ucat[k][j][i].y;
                                lcountEnd2 ++;
                            }

                            // velocity values at reference height in precursor domain (two levels to interpolate)
                            if(precursor->thisProcessorInFringe)
                            {
                                if(j == precursor->domain->abl->closestLabels[0])
                                {
                                    lsum1 += ucatP[k][j][i].y;
                                    lcount1 ++;
                                }
                                else if(j == precursor->domain->abl->closestLabels[1])
                                {
                                    lsum2 += ucatP[k][j][i].y;
                                    lcount2 ++;
                                }
                            }
                        }
                    }
                }

                // start line: reduce the values on all processors
                MPI_Allreduce(&lsumStart1,   &gsumStart1,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lsumStart2,   &gsumStart2,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lcountStart1, &gcountStart1, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);
                MPI_Allreduce(&lcountStart2, &gcountStart2, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);
                // end line: reduce the values on all processors
                MPI_Allreduce(&lsumEnd1,   &gsumEnd1,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lsumEnd2,   &gsumEnd2,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lcountEnd1, &gcountEnd1, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);
                MPI_Allreduce(&lcountEnd2, &gcountEnd2, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);
                // precursor: reduce the values on all processors
                MPI_Allreduce(&lsum1,   &gsum1,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lsum2,   &gsum2,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lcount1, &gcount1, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);
                MPI_Allreduce(&lcount2, &gcount2, 1, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);

                if(precursor->thisProcessorInFringe)
                {
                    DMDAVecRestoreArray(fda, precursor->domain->ueqn->lUcat, &ucatP);
                }

                DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
                DMDAVecRestoreArray(fda, mesh->lCent, &cent);

                gsumStart1 = gsumStart1 / gcountStart1;
                gsumStart2 = gsumStart2 / gcountStart2;
                gsumEnd1   = gsumEnd1   / gcountEnd1;
                gsumEnd2   = gsumEnd2   / gcountEnd2;
                gsum1      = gsum1      / gcount1;
                gsum2      = gsum2      / gcount2;

                PetscReal vStart   = gsumStart1 * abl->levelWeights[0] + gsumStart2 * abl->levelWeights[1];
                PetscReal vEnd     = gsumEnd1   * abl->levelWeights[0] + gsumEnd2   * abl->levelWeights[1];
                PetscReal vBarPrec = gsum1      * abl->levelWeights[0] + gsum2      * abl->levelWeights[1];

                // if first iteration initialize abl->vStart and abl->vEnd after computing instantaneous values for the first time
                if(clock->it == 0)
                {
                    abl->vStart             = vStart;
                    abl->vEnd               = vEnd;
                    abl->xDampingVBar       = vBarPrec;
                    abl->xDampingTimeStart  = clock->startTime;
                }

                // window-filtering starting and ending y-velocities
                abl->vStart       = std::exp(-clock->dt / abl->xDampingTimeWindow) * abl->vStart + (clock->dt / abl->xDampingTimeWindow) * vStart;
                abl->vEnd         = std::exp(-clock->dt / abl->xDampingTimeWindow) * abl->vEnd   + (clock->dt / abl->xDampingTimeWindow) * vEnd;
                abl->xDampingVBar = std::exp(-clock->dt / abl->xDampingTimeWindow) * abl->xDampingVBar + (clock->dt / abl->xDampingTimeWindow) * vBarPrec;

                // error at fringe exit
                abl->xDampingError    = fabs(abl->xDampingVBar - abl->vEnd);

                // velocity jump across the fringe
                abl->xDampingDeltaV   = fabs(abl->vEnd - abl->vStart);

                // predicted angular coefficient for the fringe alpha-relation
                abl->xDampingCoeff    = abl->xDampingDeltaV / abl->xDampingAlpha;

                PetscReal fringeTime  = clock->time - abl->xDampingTimeStart;
                PetscReal waitTime    = abl->xDampingTimeWindow; // (abl->xDampingEnd - abl->xDampingStart)/abl->uRef;


                // see if must correct alpha
                if
                (
                    (abl->xDampingError/abl->uRef) > 0.01 && // check if error is greater than 1%
                    fringeTime > waitTime                    // check if allowed to override alpha
                )
                {
                    abl->xDampingAlpha = fabs(abl->xDampingVBar - abl->vStart) / abl->xDampingCoeff;
                    abl->xDampingTimeStart = clock->time;
                }

                // bound alpha
                abl->xDampingAlpha = std::max(std::min(abl->xDampingAlpha, 1.0), 0.0);

                // compute useful print quantities
                PetscReal percErrVStart  = fabs(abl->xDampingVBar - abl->vStart) / abl->uRef * 100.0;
                PetscReal percErrVEnd    = abl->xDampingError / abl->uRef * 100.0;
                PetscReal percTimeFringe = fringeTime / waitTime * 100.0;

                PetscPrintf(mesh->MESH_COMM, "Correcting fringe region: errStart = %.3lf %%, errEnd = %.3lf %%, vBar = %.3lf m/s, mCoeff = %.3lf, tFringe = %.1lf %%, alpha = %.5lf\n", percErrVStart, percErrVEnd, abl->xDampingVBar, abl->xDampingCoeff, percTimeFringe, abl->xDampingAlpha);

                // get current process
                PetscMPIInt   rank;
                MPI_Comm_rank(mesh->MESH_COMM, &rank);

                // write file
                if (!rank)
                {
                    FILE *f;
                    char filen[80];
                    PetscInt width = -20;
                    sprintf(filen, "fringeRegionData");

                    if(clock->it == clock->itStart)
                    {
                        // eliminate previous file
                        unlink(filen);

                        // open a new file
                        f = fopen(filen, "a");

                        // write header line
                        word w1 = "time";
                        word w2 = "errStartPercent";
                        word w3 = "errEndPercent";
                        word w4 = "vBar";
                        word w5 = "mCoeff";
                        word w6 = "tFringePercent";
                        word w7 = "alpha";
                        PetscFPrintf(PETSC_COMM_WORLD, f, "%*s\t%*s\t%*s\t%*s\t%*s\t%*s\t%*s\n", width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str());
                    }

                    f = fopen(filen, "a");

                    PetscFPrintf(PETSC_COMM_WORLD, f, "%*.3f\t%*.3f\t%*.3f\t%*.5f\t%*.3f\t%*.3f\t%*.3f\n", width, clock->time, width, percErrVStart, width, percErrVEnd, width, abl->xDampingVBar, width, abl->xDampingCoeff, width, percTimeFringe, width, abl->xDampingAlpha);

                    fclose(f);
                }
            }

            PetscPrintf(mesh->MESH_COMM, "Solving concurrent precursor:\n");

            // solve concurrent precursor
            concurrentPrecursorSolve(abl);

            PetscPrintf(mesh->MESH_COMM, "Solving successor:\n");
        }
    }

    // update uBar for zDampingLayer (average at kLeft patch)
    if(ueqn->access->flags->isZDampingActive)
    {
        // do it only if also x and y velocity component have to be damped, otherwise
        // if only z component damping is active the uBar is zero.
        if(abl->zDampingAlsoXY)
        {
            std::vector<Cmpnts>   luBar(my);
            std::vector<Cmpnts>   guBar(my);

            if(abl->zDampingXYType == 1)
            {
                std::vector<PetscInt> ln(my);
                std::vector<PetscInt> gn(my);

                DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

                // compute uBar: the i-line averaged contravariant fluxes at the first internal
                // face. They only depend on the j-index.

                // test if this processor is on k-left boundary
                if(zs == 0)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            luBar[j].x += ucat[lzs][j][i].x;
                            luBar[j].y += ucat[lzs][j][i].y;
                            luBar[j].z += ucat[lzs][j][i].z;

                            ln[j]++;
                        }
                    }
                }

                // reduce the value
                MPI_Allreduce(&(luBar[0]), &(guBar[0]), 3*my, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&(ln[0]),    &(gn[0]),      my, MPIU_INT,  MPI_SUM,  mesh->MESH_COMM);

                // compute global mean
                for(j=1; j<my-1; j++)
                {
                    mScale(1.0/(PetscReal)gn[j], guBar[j]);
                }

                // set time averaging weights
                PetscReal mN = (PetscReal)abl->avgWeight;
                PetscReal m1 = mN  / (mN + 1.0);
                PetscReal m2 = 1.0 / (mN + 1.0);

                // cumulate uBarMean
                for(j=1; j<my-1; j++)
                {
                    abl->uBarMeanZ[j].x = m1 * abl->uBarMeanZ[j].x + m2 * guBar[j].x;
                    abl->uBarMeanZ[j].y = m1 * abl->uBarMeanZ[j].y + m2 * guBar[j].y;
                    abl->uBarMeanZ[j].z = m1 * abl->uBarMeanZ[j].z + m2 * guBar[j].z;
                }

                // increase snapshot weighting
                abl->avgWeight++;

                DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

                std::vector<PetscInt>    ().swap(ln);
                std::vector<PetscInt>    ().swap(gn);
            }
            else if(abl->zDampingXYType == 2)
            {
                precursor_ *precursor = abl->precursor;
                dataABL    *ablStat;

                if(precursor->thisProcessorInFringe)
                {
                    ablStat = precursor->domain->acquisition->statisticsABL;

                    PetscMPIInt rank; MPI_Comm_rank(precursor->domain->mesh->MESH_COMM, &rank);

                    // get the averages from the master rank of the precursor mesh comm
                    if(!rank)
                    {
                        for(j=1; j<my-1; j++)
                        {
                            luBar[j].x = ablStat->UMean[j-1];
                            luBar[j].y = ablStat->VMean[j-1];
                            luBar[j].z = 0.0;
                        }
                    }
                }

                // reduce the value
                MPI_Allreduce(&(luBar[0]), &(guBar[0]), 3*my, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

                // assign to this processor
                for(j=1; j<my-1; j++)
                {
                    abl->uBarMeanZ[j].x = guBar[j].x;
                    abl->uBarMeanZ[j].y = guBar[j].y;
                    abl->uBarMeanZ[j].z = guBar[j].z;
                }
            }

            // clean local vectors
            std::vector<Cmpnts>      ().swap(luBar);
            std::vector<Cmpnts>      ().swap(guBar);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode dampingSourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    abl_          *abl  = ueqn->access->abl;
    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***rhs, ***ucat, ***cent;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    Cmpnts        ***ucont, ***ucontP;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    precursor_    *precursor;
    domain_       *pdomain;
    PetscInt      kStart;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, Rhs,  &rhs);

    if(ueqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->fda, pdomain->ueqn->lUcont,  &ucontP);
                kStart = precursor->map.kStart;
            }
        }
    }

    // z damping layer
    PetscReal alphaZ = abl->zDampingAlpha;
    PetscReal zS     = abl->zDampingStart;
    PetscReal zE     = abl->zDampingEnd;

    // x damping layer
    PetscReal alphaX = abl->xDampingAlpha;
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;
    PetscReal xD     = abl->xDampingDelta;

    // cell center coordinates
    PetscReal x, xi, xj, xk;
    PetscReal z, zi, zj, zk;

    // Nordstrom viscosities
    PetscReal nud_x, nudi_x, nudj_x, nudk_x;

    // Rayleigh viscosities
    PetscReal nud_z, nudi_z, nudj_z, nudk_z;

    // Stipa viscosities
    PetscReal nud_x_s, nudi_x_s, nudj_x_s, nudk_x_s;

    // damping viscosity for zDamping exlusion in fringe region
    PetscReal nudI = 1.0;
    PetscReal nudJ = 1.0;
    PetscReal nudK = 1.0;

    // loop over internal cell faces - include right boundary faces which will be periodic
    // at the beginning. Then they will be zeroed if applicable when building the SNES rhs.
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(ueqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    x     = cent[k][j][i].x;
                    xi    = cent[k][j][i+1].x;
                    xj    = cent[k][j+1][i].x;
                    xk    = cent[k+1][j][i].x;

                    // compute Nordstrom viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_x   = viscNordstrom(alphaX, xS, xE, xD, x);
                    nudi_x  = viscNordstrom(alphaX, xS, xE, xD, xi);
                    nudj_x  = viscNordstrom(alphaX, xS, xE, xD, xj);
                    nudk_x  = viscNordstrom(alphaX, xS, xE, xD, xk);

                    // compute Stipa viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_x_s    = viscStipa(xS, xE, xD, x);
                    nudi_x_s   = viscStipa(xS, xE, xD, xi);
                    nudj_x_s   = viscStipa(xS, xE, xD, xj);
                    nudk_x_s   = viscStipa(xS, xE, xD, xk);

                    // interpolate Stipa viscosity at cell faces
                    nudI       = central(nud_x_s, nudi_x_s);
                    nudJ       = central(nud_x_s, nudj_x_s);
                    nudK       = central(nud_x_s, nudk_x_s);

                    // X DAMPING LAYER
                    // ---------------

                    if(abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2)
                    {
                        Cmpnts uBarInstX  = nSet(abl->uBarInstX[j][i]);
                        Cmpnts uBarInstXi = nSet(abl->uBarInstX[j][i+1]);
                        Cmpnts uBarInstXj = nSet(abl->uBarInstX[j+1][i]);

                        // i-fluxes
                        rhs[k][j][i].x
                        +=
                        scale * central(nud_x, nudi_x) *
                        (
                            (
                                (central(uBarInstX.x, uBarInstXi.x) - central(ucat[k][j][i].x, ucat[k][j][i+1].x)) * icsi[k][j][i].x +
                                (central(uBarInstX.y, uBarInstXi.y) - central(ucat[k][j][i].y, ucat[k][j][i+1].y)) * icsi[k][j][i].y +
                                (central(uBarInstX.z, uBarInstXi.z) - central(ucat[k][j][i].z, ucat[k][j][i+1].z)) * icsi[k][j][i].z
                            )
                        );

                        // j-fluxes
                        rhs[k][j][i].y
                        +=
                        scale * central(nud_x, nudj_x) *
                        (
                            (
                                (central(uBarInstX.x, uBarInstXj.x) - central(ucat[k][j][i].x, ucat[k][j+1][i].x)) * jeta[k][j][i].x +
                                (central(uBarInstX.y, uBarInstXj.y) - central(ucat[k][j][i].y, ucat[k][j+1][i].y)) * jeta[k][j][i].y +
                                (central(uBarInstX.z, uBarInstXj.z) - central(ucat[k][j][i].z, ucat[k][j+1][i].z)) * jeta[k][j][i].z
                            )
                        );

                        // k-fluxes
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_x, nudk_x) *
                        (
                            (
                                (uBarInstX.x - central(ucat[k][j][i].x, ucat[k+1][j][i].x)) * kzet[k][j][i].x +
                                (uBarInstX.y - central(ucat[k][j][i].y, ucat[k+1][j][i].y)) * kzet[k][j][i].y +
                                (uBarInstX.z - central(ucat[k][j][i].z, ucat[k+1][j][i].z)) * kzet[k][j][i].z
                            )
                        );
                    }
                    else if(abl->xFringeUBarSelectionType == 3)
                    {
                        PetscReal uBarContK;
                        PetscReal uBarContJ;
                        PetscReal uBarContI;

                        // note: here we can use contravariant fluxes since uBar is defined at every point

                        if(precursor->thisProcessorInFringe)
                        {
                            uBarContI = ucontP[k+kStart][j][i].x;
                            uBarContJ = ucontP[k+kStart][j][i].y;
                            uBarContK = ucontP[k+kStart][j][i].z;
                        }
                        else
                        {
                            uBarContI = ucont[k][j][i].x;
                            uBarContJ = ucont[k][j][i].y;
                            uBarContK = ucont[k][j][i].z;
                        }

                        // i-fluxes
                        rhs[k][j][i].x
                        +=
                        scale * central(nud_x, nudi_x) *
                        (
                            uBarContI - ucont[k][j][i].x
                        );

                        // j-fluxes
                        rhs[k][j][i].y
                        +=
                        scale * central(nud_x, nudj_x) *
                        (
                            uBarContJ - ucont[k][j][i].y
                        );

                        // k-fluxes
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_x, nudk_x) *
                        (
                            uBarContK - ucont[k][j][i].z
                        );
                    }
                }

                // Z DAMPING LAYER
                // ---------------
                if(ueqn->access->flags->isZDampingActive)
                {
                    // compute cell center z at i,j,k and i,j+1,k points
                    z     = (cent[k][j][i].z   - mesh->bounds.zmin);
                    zj    = (cent[k][j+1][i].z - mesh->bounds.zmin);

                    // compute Rayleigh viscosity at i,j,k and i,j+1,k points
                    nud_z   = viscRayleigh(alphaZ, zS, zE, z);
                    nudj_z  = viscRayleigh(alphaZ, zS, zE, zj);

                    // damp also x and y components (exclude in xFringe if present)
                    if(abl->zDampingAlsoXY)
                    {
                        // compute cell center z at i+1,j,k and i,j,k+1 points
                        zi    = (cent[k][j][i+1].z - mesh->bounds.zmin);
                        zk    = (cent[k+1][j][i].z - mesh->bounds.zmin);

                        // compute Rayleigh viscosity at i+1,j,k and i,j,k+1 points
                        nudi_z  = viscRayleigh(alphaZ, zS, zE, zi);
                        nudk_z  = viscRayleigh(alphaZ, zS, zE, zk);

                        // i-fluxes: dampen w.r.t. uBarMean
                        rhs[k][j][i].x
                        +=
                        scale * central(nud_z, nudi_z) * nudI *
                        (
                            (
                                (abl->uBarMeanZ[j].x - central(ucat[k][j][i].x, ucat[k][j][i+1].x)) * icsi[k][j][i].x +
                                (abl->uBarMeanZ[j].y - central(ucat[k][j][i].y, ucat[k][j][i+1].y)) * icsi[k][j][i].y +
                                (abl->uBarMeanZ[j].z - central(ucat[k][j][i].z, ucat[k][j][i+1].z)) * icsi[k][j][i].z
                            )
                        );

                        // k-fluxes: dampen w.r.t. uBarMean
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_z, nudk_z) * nudK *
                        (
                            (
                                (abl->uBarMeanZ[j].x - central(ucat[k][j][i].x, ucat[k+1][j][i].x)) * kzet[k][j][i].x +
                                (abl->uBarMeanZ[j].y - central(ucat[k][j][i].y, ucat[k+1][j][i].y)) * kzet[k][j][i].y +
                                (abl->uBarMeanZ[j].z - central(ucat[k][j][i].z, ucat[k+1][j][i].z)) * kzet[k][j][i].z
                            )
                        );
                    }

                    // j-fluxes: total damping to reach no penetration at jRight (damp also in xFringe if present)
                    rhs[k][j][i].y
                    +=
                    -1.0 * scale * central(nud_z, nudj_z) *
                    (
                        central(ucat[k][j][i].x, ucat[k][j+1][i].x) * jeta[k][j][i].x +
                        central(ucat[k][j][i].y, ucat[k][j+1][i].y) * jeta[k][j][i].y +
                        central(ucat[k][j][i].z, ucat[k][j+1][i].z) * jeta[k][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, Rhs,  &rhs);

    if(ueqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->fda, pdomain->ueqn->lUcont,  &ucontP);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Coriolis(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***rhs, ***ucat, ***cent;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    PetscReal     fc      = ueqn->access->abl->fc; // coriolis parameter

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // The coriolis force is assumed to be present only in the x and y (wall parallel)
    // directions. Since the equations are projected along the generalized curvilinear
    // coordinates we have to dot this vector term with the face area vectors. Note that
    // if the eta axis is aligned with the z cartesian direction, the dotting will output
    // zero since eta face area vectors have zero component in x and y cartesian directions.

    // damping viscosity for fringe region exclusion
    double nudI;
    double nudJ;
    double nudK;

    // fringe region parameters (set only if active)
    double xS;
    double xE;
    double xD;

    if(ueqn->access->flags->isXDampingActive)
    {
        xS     = ueqn->access->abl->xDampingStart;
        xE     = ueqn->access->abl->xDampingEnd;
        xD     = ueqn->access->abl->xDampingDelta;
    }
    else
    {
        nudI = nudJ = nudK = 1.0;
    }

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(fda, Rhs,  &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(ueqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    double x     = cent[k][j][i].x;
                    double xi    = cent[k][j][i+1].x;
                    double xj    = cent[k][j+1][i].x;
                    double xk    = cent[k+1][j][i].x;

                    // compute Stipa viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    double nud_x   = viscStipa(xS, xE, xD, x);
                    double nudi_x  = viscStipa(xS, xE, xD, xi);
                    double nudj_x  = viscStipa(xS, xE, xD, xj);
                    double nudk_x  = viscStipa(xS, xE, xD, xk);

                    nudI           = central(nud_x, nudi_x);
                    nudJ           = central(nud_x, nudj_x);
                    nudK           = central(nud_x, nudk_x);
                }

                if
                (
                    isFluidIFace(k, j, i, i+1, nvert)
                )
                {
                    rhs[k][j][i].x
                    +=
                    scale * nudI *
                    (
                        -2.0 *
                        (
                            - fc * central(ucat[k][j][i].y, ucat[k][j][i+1].y) * icsi[k][j][i].x +
                              fc * central(ucat[k][j][i].x, ucat[k][j][i+1].x) * icsi[k][j][i].y
                        )
                    );
                }

                if
                (
                    isFluidJFace(k, j, i, j+1, nvert)
                )
                {
                    rhs[k][j][i].y
                    +=
                    scale * nudJ *
                    (
                        -2.0 *
                        (
                          - fc * central(ucat[k][j][i].y, ucat[k][j+1][i].y) * jeta[k][j][i].x +
                            fc * central(ucat[k][j][i].x, ucat[k][j+1][i].x) * jeta[k][j][i].y
                        )
                    );
                }

                if
                (
                    isFluidKFace(k, j, i, k+1, nvert)
                )
                {
                    rhs[k][j][i].z
                    +=
                    scale * nudK *
                    (
                        -2.0 *
                        (
                          - fc * central(ucat[k][j][i].y, ucat[k+1][j][i].y) * kzet[k][j][i].x +
                            fc * central(ucat[k][j][i].x, ucat[k+1][j][i].x) * kzet[k][j][i].y
                        )
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    DMDAVecRestoreArray(fda, Rhs,  &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SideForce(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    // fictitious force evaluated as the Coriolis force with opposite direction:
    // SideForce = -K * Coriolis

    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***rhs, ***ucat, ***cent;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    PetscReal     fc      = ueqn->access->abl->fc; // coriolis parameter
    PetscReal     xStart  = ueqn->access->abl->xStartSideF,
                  zStart  = ueqn->access->abl->zStartSideF,
                  xEnd    = ueqn->access->abl->xEndSideF,
                  zEnd    = ueqn->access->abl->zEndSideF;
    double        K       = 5.0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(fda, Rhs,  &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                double coeff_i = 0.0, coeff_j = 0.0, coeff_k = 0.0;

                // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                double xi    = 0.5 * (cent[k][j][i].x + cent[k][j][i+1].x) - mesh->bounds.xmin;
                double xj    = 0.5 * (cent[k][j][i].x + cent[k][j+1][i].x) - mesh->bounds.xmin;
                double xk    = 0.5 * (cent[k][j][i].x + cent[k+1][j][i].x) - mesh->bounds.xmin;

                // compute cell center z at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                double zi    = 0.5 * (cent[k][j][i].z + cent[k][j][i+1].z) - mesh->bounds.zmin;
                double zj    = 0.5 * (cent[k][j][i].z + cent[k][j+1][i].z) - mesh->bounds.zmin;
                double zk    = 0.5 * (cent[k][j][i].z + cent[k+1][j][i].z) - mesh->bounds.zmin;

                if(xi < xEnd && xi > xStart && zi < zEnd && zi > zStart) coeff_i = 1;
                if(xj < xEnd && xj > xStart && zj < zEnd && zj > zStart) coeff_j = 1;
                if(xk < xEnd && xk > xStart && zk < zEnd && zk > zStart) coeff_k = 1;

                rhs[k][j][i].x
                +=
                scale * coeff_i *
                (
                    2.0 * K *
                    (
                        fc * 10.0 * icsi[k][j][i].y
                    )
                );

                if
                (
                    isFluidJFace(k, j, i, j+1, nvert)
                )
                {
                    rhs[k][j][i].y
                    +=
                    scale * coeff_j *
                    0.0 ;
                }

                if
                (
                    isFluidKFace(k, j, i, k+1, nvert)
                )
                {
                    rhs[k][j][i].z
                    +=
                    scale * coeff_k *
                    0.0 ;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    DMDAVecRestoreArray(fda, Rhs,  &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Buoyancy(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh  = ueqn->access->mesh;
    teqn_         *teqn = ueqn->access->teqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***rhs, ***gcont;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert, ***tmprt;

    PetscInt           i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        gravity;
                  gravity.x = 0;
                  gravity.y = 0;
                  gravity.z = -9.81;
    PetscReal     tRef      = ueqn->access->abl->tRef;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->gCont,  &gcont);

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);

    // loop over i-face centers
    for(k=zs; k<lze; k++)
    {
        for(j=ys; j<lye; j++)
        {
            for(i=xs; i<lxe; i++)
            {
                if(j > 0 && k > 0)
                {
                    gcont[k][j][i].x
                    =
                    (
                        gravity.x * icsi[k][j][i].x +
                        gravity.y * icsi[k][j][i].y +
                        gravity.z * icsi[k][j][i].z
                    );
                }

                if(i > 0 && k > 0)
                {
                    gcont[k][j][i].y
                    =
                    (
                        gravity.x * jeta[k][j][i].x +
                        gravity.y * jeta[k][j][i].y +
                        gravity.z * jeta[k][j][i].z
                    );
                }

                if(i > 0 && j > 0)
                {
                    gcont[k][j][i].z
                    =
                    (
                        gravity.x * kzet[k][j][i].x +
                        gravity.y * kzet[k][j][i].y +
                        gravity.z * kzet[k][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, teqn->lTmprt, &tmprt);
    DMDAVecGetArray(fda, Rhs,  &rhs);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // interpolate at cell faces
                PetscReal tempI = 0.5*(tmprt[k][j][i] + tmprt[k][j][i+1]);

                if (isFluidIFace(k, j, i, i+1, nvert))
                {
                    rhs[k][j][i].x
                    +=
                    scale *
                    (
                        gcont[k][j][i].x *
                        (tRef - tempI) / tRef
                    );
                }


                // interpolate at cell faces
                PetscReal tempJ =  0.5*(tmprt[k][j][i] + tmprt[k][j+1][i]);

                if (isFluidJFace(k, j, i, j+1, nvert))
                {
                    rhs[k][j][i].y
                    +=
                    scale *
                    (
                        gcont[k][j][i].y *
                        (tRef - tempJ) / tRef
                    );
                }


                // interpolate at cell faces
                PetscReal tempK =  0.5*(tmprt[k][j][i] + tmprt[k+1][j][i]);

                if(isFluidKFace(k, j, i, k+1, nvert))
                {
                    rhs[k][j][i].z
                    +=
                    scale *
                    (
                        gcont[k][j][i].z *
                        (tRef - tempK) / tRef
                    );
                }

            }
        }
    }

    DMDAVecRestoreArray(fda, Rhs,  &rhs);
    DMDAVecRestoreArray(da,  teqn->lTmprt, &tmprt);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, ueqn->gCont, &gcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode contravariantToCartesian(ueqn_ *ueqn)
{
    mesh_           *mesh = ueqn->access->mesh;
    DM               da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    PetscReal        mat[3][3], det, det0, det1, det2;

    PetscReal        ***aj, ***nvert;
    Cmpnts           ***ucat, ***lucont;

    PetscReal q[3];  // local working array

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    Cmpnts    ***csi,  ***eta,  ***zet;

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);
    DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);

    // interpolate ucont from faces to cell centers
    // transform from contravariant to cartesian

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if ( isFluidCell(k, j, i, nvert) )
                {
                    mat[0][0] = (csi[k][j][i].x);
                    mat[0][1] = (csi[k][j][i].y);
                    mat[0][2] = (csi[k][j][i].z);

                    mat[1][0] = (eta[k][j][i].x);
                    mat[1][1] = (eta[k][j][i].y);
                    mat[1][2] = (eta[k][j][i].z);

                    mat[2][0] = (zet[k][j][i].x);
                    mat[2][1] = (zet[k][j][i].y);
                    mat[2][2] = (zet[k][j][i].z);


                    PetscInt iL = i-1,  iR = i;
                    PetscInt jL = j-1,  jR = j;
                    PetscInt kL = k-1,  kR = k;

                    // interpolate ucont from opposite faces to cell centers
                    q[0] = central( lucont[k][j][iL].x, lucont[k][j][iR].x);
                    q[1] = central( lucont[k][jL][i].y, lucont[k][jR][i].y);
                    q[2] = central( lucont[kL][j][i].z, lucont[kR][j][i].z);

                    // transform to cartesian
                    det =   mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

                    det0 =  q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                            q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                    det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                            q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                    det2 =  q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                            q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                            q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                    ucat[k][j][i].x = det0 / det;
                    ucat[k][j][i].y = det1 / det;
                    ucat[k][j][i].z = det2 / det;
                }
            }
        }
    }

    DMDAVecRestoreArray (fda, ueqn->Ucat,  &ucat);
    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd  (fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMDAVecRestoreArray (fda, ueqn->lUcont, &lucont);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode contravariantToCartesianGeneric(mesh_ *mesh, Vec &lCont, Vec &lCat)
{
    DM               da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    PetscReal        mat[3][3], det, det0, det1, det2;

    PetscReal        ***aj, ***nvert;
    Cmpnts           ***lcat, ***lcont;
    Cmpnts           ***csi,  ***eta,  ***zet;

    PetscReal q[3];  // local working array

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, lCont, &lcont);
    DMDAVecGetArray(fda, lCat,  &lcat);

    // interpolate ucont from faces to cell centers
    // transform from contravariant to cartesian

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if ( isFluidCell(k, j, i, nvert) )
                {
                    mat[0][0] = (csi[k][j][i].x);
                    mat[0][1] = (csi[k][j][i].y);
                    mat[0][2] = (csi[k][j][i].z);

                    mat[1][0] = (eta[k][j][i].x);
                    mat[1][1] = (eta[k][j][i].y);
                    mat[1][2] = (eta[k][j][i].z);

                    mat[2][0] = (zet[k][j][i].x);
                    mat[2][1] = (zet[k][j][i].y);
                    mat[2][2] = (zet[k][j][i].z);


                    PetscInt iL = i-1,  iR = i;
                    PetscInt jL = j-1,  jR = j;
                    PetscInt kL = k-1,  kR = k;

                    // interpolate ucont from opposite faces to cell centers
                    q[0] = central( lcont[k][j][iL].x, lcont[k][j][iR].x);
                    q[1] = central( lcont[k][jL][i].y, lcont[k][jR][i].y);
                    q[2] = central( lcont[kL][j][i].z, lcont[kR][j][i].z);

                    // transform to cartesian
                    det =   mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

                    det0 =  q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                            q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                    det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                            q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                    det2 =  q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                            q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                            q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                    lcat[k][j][i].x = det0 / det;
                    lcat[k][j][i].y = det1 / det;
                    lcat[k][j][i].z = det2 / det;
                }
            }
        }
    }

    DMDAVecRestoreArray (fda, lCat,  &lcat);
    DMLocalToLocalBegin(fda, lCat, INSERT_VALUES, lCat);
    DMLocalToLocalEnd  (fda, lCat, INSERT_VALUES, lCat);
    DMDAVecRestoreArray (fda, lCont, &lcont);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode adjustFluxes(ueqn_ *ueqn)
{
    mesh_           *mesh = ueqn->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont, ***lucont,
                     ***coor;
    PetscReal        epsilon = 1.e-10;

    PetscScalar      ratio;
    PetscScalar      lFluxIn = 0.0, lFluxOut = 0.0,
                     FluxIn  = 0.0, FluxOut  = 0.0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // compute inflow flux at k-left boundary
        if (zs==0)
        {
            // k-right boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }

        // compute outflow flux at k-right boundary
        if (ze==mz)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0)
        {
            // k-right boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my)
        {
            // k-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0)
        {
            // i-right boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon)
                    {
                        lFluxIn   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx)
        {
            // k-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }
    }

    // cumulate the net influx and net outflux
    MPI_Allreduce(&lFluxIn, &FluxIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //PetscPrintf(mesh->MESH_COMM, "\n fluxin = %lf, fluxout = %lf\n", FluxIn, FluxOut);

    if(FluxOut > 1.e-15)
    {
        ratio = (FluxIn) / FluxOut;
    }
    else
    {
        ratio = 1.0;
    }

    // correct only outflow patches - outflow patches are assumed to be the right side patches
    // this is not generic but for inflow - outflow problems with outflow on the right patches this should work
    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {

        // correct outflow velocity on the k-right patch

            if (ze==mz  && (mesh->boundaryU.kRight!="inletFunction") && (mesh->boundaryU.kRight!="fixedValue"))
            {
                k = mz-2;

                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].z) > epsilon)
                        {
                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].z *= ratio;

                        }
                    }
                }
            }

            if (zs==0  && (mesh->boundaryU.kLeft!="inletFunction") && (mesh->boundaryU.kLeft!="fixedValue"))
            {
                k = 0;

                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].z) > epsilon)
                        {
                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].z /= ratio;

                        }
                    }
                }
            }

    }

    if
    (
       !mesh->j_periodic && !mesh->jj_periodic
    )
    {

        // correct outflow velocity on the j-right patch

            if (ye==my  && (mesh->boundaryU.jRight!="inletFunction") && (mesh->boundaryU.jRight!="fixedValue"))
            {
                j = my-2;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].y) > epsilon)
                        {

                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].y *= ratio;

                        }
                    }
                }
            }

            if (ys==0  && (mesh->boundaryU.jLeft!="inletFunction") && (mesh->boundaryU.jLeft!="fixedValue"))
            {
                j = 0;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].y) > epsilon)
                        {

                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].y /= ratio;

                        }
                    }
                }
            }

    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {

        // correct outflow velocity on the i-right patch

        if (xe==mx  && (mesh->boundaryU.iRight!="inletFunction") && (mesh->boundaryU.iRight!="fixedValue"))
        {
            i = mx-2;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    if (fabs(ucont[k][j][i].x) > epsilon)
                    {

                        // correct contrav. vel. at boundary faces
                        ucont[k][j][i].x *= ratio;
                    }
                }
            }
        }

        if (xs==0  && (mesh->boundaryU.iLeft!="inletFunction") && (mesh->boundaryU.iLeft!="fixedValue"))
        {
            i = 0;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    if (fabs(ucont[k][j][i].x) > epsilon)
                    {

                        // correct contrav. vel. at boundary faces
                        ucont[k][j][i].x /= ratio;
                    }
                }
            }
        }
    }

    // Recalculate flux to check if corrected
    lFluxIn = 0.0, lFluxOut = 0.0,
    FluxIn  = 0.0, FluxOut  = 0.0;

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // compute inflow flux at k-left boundary
        if (zs==0)
        {
            // k-right boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }

        // compute outflow flux at k-right boundary
        if (ze==mz)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0)
        {
            // k-right boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my)
        {
            // k-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0)
        {
            // i-right boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon)
                    {
                        lFluxIn   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx)
        {
            // k-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }
    }

    // cumulate the net influx and net outflux
    MPI_Allreduce(&lFluxIn, &FluxIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //PetscPrintf(mesh->MESH_COMM, "After correction fluxin = %lf, fluxout = %lf\n", FluxIn, FluxOut);

    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &lucont);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode adjustFluxesOverset(ueqn_ *ueqn)
{
    mesh_           *mesh = ueqn->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont, ***lucont,
                     ***coor;
    PetscReal        epsilon = 1.e-10;

    PetscScalar      ratio;
    PetscScalar      lFluxIn = 0.0, lFluxOut = 0.0,
                     FluxIn  = 0.0, FluxOut  = 0.0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // compute inflow flux at k-left boundary
        if (zs==0)
        {
            // k-right boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }

        // compute outflow flux at k-right boundary
        if (ze==mz)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0)
        {
            // k-right boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my)
        {
            // k-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0)
        {
            // i-right boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon)
                    {
                        lFluxIn   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx)
        {
            // k-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }
    }

    // cumulate the net influx and net outflux
    MPI_Allreduce(&lFluxIn, &FluxIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //PetscPrintf(mesh->MESH_COMM, "\n fluxin = %lf, fluxout = %lf\n", FluxIn, FluxOut);

    if(FluxOut > 1.e-15)
    {
        ratio = (FluxIn) / FluxOut;
    }
    else
    {
        ratio = 1.0;
    }

    // correct only outflow patches - outflow patches are assumed to be the right side patches
    // this is not generic but for inflow - outflow problems with outflow on the right patches this should work
    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {

        // correct outflow velocity on the k-right patch

            if (ze==mz)
            {
                k = mz-2;

                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].z) > epsilon)
                        {
                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].z *= sqrt(ratio);

                        }
                    }
                }
            }

            if (zs==0)
            {
                k = 0;

                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].z) > epsilon)
                        {
                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].z /= sqrt(ratio);

                        }
                    }
                }
            }

    }

    if
    (
       !mesh->j_periodic && !mesh->jj_periodic
    )
    {

        // correct outflow velocity on the j-right patch

            if (ye==my)
            {
                j = my-2;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].y) > epsilon)
                        {

                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].y *= sqrt(ratio);

                        }
                    }
                }
            }

            if (ys==0)
            {
                j = 0;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if (fabs(ucont[k][j][i].y) > epsilon)
                        {

                            // correct contrav. vel. at boundary faces
                            ucont[k][j][i].y /= sqrt(ratio);

                        }
                    }
                }
            }

    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {

        // correct outflow velocity on the i-right patch

        if (xe==mx)
        {
            i = mx-2;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    if (fabs(ucont[k][j][i].x) > epsilon)
                    {

                        // correct contrav. vel. at boundary faces
                        ucont[k][j][i].x *= sqrt(ratio);
                    }
                }
            }
        }

        if (xs==0)
        {
            i = 0;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    if (fabs(ucont[k][j][i].x) > epsilon)
                    {

                        // correct contrav. vel. at boundary faces
                        ucont[k][j][i].x /= sqrt(ratio);
                    }
                }
            }
        }
    }

    // Recalculate flux to check if corrected
    lFluxIn = 0.0, lFluxOut = 0.0,
    FluxIn  = 0.0, FluxOut  = 0.0;

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // compute inflow flux at k-left boundary
        if (zs==0)
        {
            // k-right boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }

        // compute outflow flux at k-right boundary
        if (ze==mz)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].z) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].z;

                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0)
        {
            // k-right boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxIn   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my)
        {
            // k-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].y) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].y;
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }

                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0)
        {
            // i-right boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon)
                    {
                        lFluxIn   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx)
        {
            // k-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux and area

                    if(fabs(ucont[k][j][i].x) > epsilon){
                        lFluxOut   +=  ucont[k][j][i].x;
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }

                }
            }
        }
    }

    // cumulate the net influx and net outflux
    MPI_Allreduce(&lFluxIn, &FluxIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //PetscPrintf(mesh->MESH_COMM, "After correction fluxin = %lf, fluxout = %lf\n", FluxIn, FluxOut);

    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &lucont);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    // In this function the viscous + divergence term of the momentum equation are
    // discretized at cell faces. The disposition is staggered, meaning that the
    // component of the face-centered quantities are not located at the same point
    // but rather at the corresponding faces.
    // First at every cell face (in the 1st 2nd and 3rd curvilinear direction) the
    // fluxes are evaluated in cartesian coordinates. Secondly the net flux is evaluated
    // for every cell in cartesian coordinates, summing up the previously calculated fluxes.
    // Third, this cell centered flux is interpolated at the faces and dotted with that face's
    // area vector, meaning that it will point in that face's curiviliear coordinate direction.
    // Note: if the flow is non-periodic in a give direction, only the contravariant velocity at
    // the internal faces is solved for. If the flow is periodic in that direction also the
    // contravariant velocity at the right boundary faces is solved for and fluxes at the left
    // boundary faces are calculated as if that face was the right boundary face (as if left and
    // right boundaries were the same boundary).

    mesh_            *mesh   = ueqn->access->mesh;
    les_             *les    = ueqn->access->les;
    constants_       *cst    = ueqn->access->constants;
    DM               da      = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info    = mesh->info;
    PetscInt         xs      = info.xs, xe = info.xs + info.xm;
    PetscInt         ys      = info.ys, ye = info.ys + info.ym;
    PetscInt         zs      = info.zs, ze = info.zs + info.zm;
    PetscInt         mx      = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont, ***ucat;

    Cmpnts           ***csi,  ***eta,  ***zet;
    Cmpnts           ***icsi, ***ieta, ***izet;
    Cmpnts           ***jcsi, ***jeta, ***jzet;
    Cmpnts           ***kcsi, ***keta, ***kzet;

    PetscReal        ***nvert, ***lnu_t;

    Cmpnts           ***div1,  ***div2,  ***div3;                               // divergence & cumulative fluxes
    Cmpnts           ***visc1, ***visc2, ***visc3;                              // viscous terms
    Cmpnts           ***rhs,   ***fp,    ***limiter;                            // right hand side
    PetscReal        ***aj,    ***iaj,   ***jaj, ***kaj;                        // cell and face jacobians

    PetscReal        dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords
    PetscReal        csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;      // surface area vectors components
    PetscReal        g11, g21, g31;                                             // metric tensor components
    PetscReal        r11, r21, r31,
                     r12, r22, r32,
                     r13, r23, r33;

    PetscScalar      solid = 0.5;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, Rhs,  &rhs);

    // get fundamental distributed arrays
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);

    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);

    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lIAj, &iaj);
    DMDAVecGetArray(da,  mesh->lJAj, &jaj);
    DMDAVecGetArray(da,  mesh->lKAj, &kaj);

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    // ---------------------------------------------------------------------- //
    // FORM DIVERGENCE AND VISCOUS CONTRIBUTIONS                              //
    // ---------------------------------------------------------------------- //

    VecSet(ueqn->lFp,     0.0);
    VecSet(ueqn->lDiv1,   0.0);
    VecSet(ueqn->lDiv2,   0.0);
    VecSet(ueqn->lDiv3,   0.0);
    VecSet(ueqn->lVisc1,  0.0);
    VecSet(ueqn->lVisc2,  0.0);
    VecSet(ueqn->lVisc3,  0.0);

    // get distributed arrays
    DMDAVecGetArray(fda, ueqn->lDiv1, &div1);
    DMDAVecGetArray(fda, ueqn->lDiv2, &div2);
    DMDAVecGetArray(fda, ueqn->lDiv3, &div3);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
    }

    // i direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(j==0 || k==0) continue;

                // get 1/V at the i-face
                PetscReal ajc = iaj[k][j][i];

                // get face normals
                csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_i (mesh, i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                // compute metric tensor - WARNING: there is a factor of 1/J^2 if using face area vectors
                //                                  must multiply for ajc in viscous term!!!
                g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                // compute cartesian velocity derivatives w.r.t. cartesian coords
                r11 = dudc * csi0 + dude * eta0 + dudz * zet0;    //du_dx / J -> another factor of / J is added in the viscous term
                r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;    //dv_dx / J -> another factor of / J is added in the viscous term
                r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;    //dw_dx / J -> another factor of / J is added in the viscous term

                r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

                Compute_du_dxyz
                (
                    mesh,
                    csi0,   csi1,   csi2,   eta0,   eta1,   eta2,   zet0,   zet1,   zet2,   ajc,
                    dudc,   dvdc,   dwdc,   dude,   dvde,   dwde,   dudz,   dvdz,   dwdz,
                    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                PetscInt    iL, iR;
                PetscReal denom;
                getFace2Cell4StencilCsi(mesh, i, mx, &iL, &iR, &denom);

                if(isIBMIFace(k, j, iL, iR, nvert) )
                {
                    iL = i;
                    iR = i+1;
                }

                if(ueqn->centralDiv)
                {
                    iL = i, iR=i+1;
                    denom=1.;
                }

                PetscReal ucon = ucont[k][j][i].x;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                // test: IBM active and not on the i-right boundary face
                if (i!=mx-2 && isIBMCell(k, j, i, nvert))
                {
                  div1[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k][j][i+2].x - 2. * ucat[k][j][i+1].x + 3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
                  -up * (0.125 * (-    ucat[k][j][i  ].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) +  ucat[k][j][i  ].x);

                  div1[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k][j][i+2].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) +  ucat[k][j][i+1].y) +
                  -up * (0.125 * (-    ucat[k][j][i  ].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) +  ucat[k][j][i  ].y);

                  div1[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k][j][i+2].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
                  -up * (0.125 * (-    ucat[k][j][i  ].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) +  ucat[k][j][i  ].z);
                }
                // test: IBM active and not on the i-left boundary face
                else if (i!=0 && isIBMCell(k, j, i+1, nvert) )
                {
                  div1[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k][j][i+1].x -  2. * ucat[k][j][i+1].x +  3. * ucat[k][j][i  ].x) + ucat[k][j][i+1].x) +
                  -up * (0.125 * (-    ucat[k][j][i-1].x -  2. * ucat[k][j][i  ].x +  3. * ucat[k][j][i+1].x) + ucat[k][j][i  ].x);

                  div1[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k][j][i+1].y -  2. * ucat[k][j][i+1].y +  3. * ucat[k][j][i  ].y) + ucat[k][j][i+1].y) +
                  -up * (0.125 * (-    ucat[k][j][i-1].y -  2. * ucat[k][j][i  ].y +  3. * ucat[k][j][i+1].y) + ucat[k][j][i  ].y);

                  div1[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k][j][i+1].z -  2. * ucat[k][j][i+1].z +  3. * ucat[k][j][i  ].z) + ucat[k][j][i+1].z) +
                  -up * (0.125 * (-    ucat[k][j][i-1].z -  2. * ucat[k][j][i  ].z +  3. * ucat[k][j][i+1].z) + ucat[k][j][i  ].z);
                }
                // not IBM
                else
                {
                    // test: inviscid flow or weno3Div
                    if(ueqn->inviscid || ueqn->weno3Div)
                    {
                        div1[k][j][i].x = - ucont[k][j][i].x *
                        weno3
                        (
                            ucat[k][j][iL].x, ucat[k][j][i].x, ucat[k][j][i+1].x, ucat[k][j][iR].x,
                            ucont[k][j][i].x
                        );

                        div1[k][j][i].y = - ucont[k][j][i].x *
                        weno3
                        (
                            ucat[k][j][iL].y, ucat[k][j][i].y, ucat[k][j][i+1].y, ucat[k][j][iR].y,
                            ucont[k][j][i].x
                        );

                        div1[k][j][i].z = - ucont[k][j][i].x *
                        weno3
                        (
                            ucat[k][j][iL].z, ucat[k][j][i].z, ucat[k][j][i+1].z, ucat[k][j][iR].z,
                            ucont[k][j][i].x
                        );
                    }
                    else
                    {
                        // central scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
                            div1[k][j][i].x = - ucont[k][j][i].x *
                            central
                            (
                                ucat[k][j][i].x, ucat[k][j][i+1].x
                            );

                            div1[k][j][i].y = - ucont[k][j][i].x *
                            central
                            (
                                ucat[k][j][i].y, ucat[k][j][i+1].y
                            );

                            div1[k][j][i].z = - ucont[k][j][i].x *
                            central
                            (
                                ucat[k][j][i].z, ucat[k][j][i+1].z
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div1[k][j][i].x = - ucont[k][j][i].x *
                            centralUpwind
                            (
                                ucat[k][j][iL].x, ucat[k][j][i].x, ucat[k][j][i+1].x, ucat[k][j][iR].x,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );

                            div1[k][j][i].y = - ucont[k][j][i].x *
                            centralUpwind
                            (
                                ucat[k][j][iL].y, ucat[k][j][i].y, ucat[k][j][i+1].y, ucat[k][j][iR].y,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );

                            div1[k][j][i].z = - ucont[k][j][i].x *
                            centralUpwind
                            (
                                ucat[k][j][iL].z, ucat[k][j][i].z, ucat[k][j][i+1].z, ucat[k][j][iR].z,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][j][iL ] * nMag(csi[k][j][iL ]));
                            PetscReal d1 = 1.0 / (aj[k][j][i  ] * nMag(csi[k][j][i  ]));
                            PetscReal d2 = 1.0 / (aj[k][j][i+1] * nMag(csi[k][j][i+1]));
                            PetscReal d3 = 1.0 / (aj[k][j][iR ] * nMag(csi[k][j][iR ]));

                            div1[k][j][i].x = - ucont[k][j][i].x *
                            wCentralUpwind
                            (
                                ucat[k][j][iL].x, ucat[k][j][i].x, ucat[k][j][i+1].x, ucat[k][j][iR].x,
                                d0, d1, d2, d3,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );

                            div1[k][j][i].y = - ucont[k][j][i].x *
                            wCentralUpwind
                            (
                                ucat[k][j][iL].y, ucat[k][j][i].y, ucat[k][j][i+1].y, ucat[k][j][iR].y,
                                d0, d1, d2, d3,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );

                            div1[k][j][i].z = - ucont[k][j][i].x *
                            wCentralUpwind
                            (
                                ucat[k][j][iL].z, ucat[k][j][i].z, ucat[k][j][i+1].z, ucat[k][j][iR].z,
                                d0, d1, d2, d3,
                                ucont[k][j][i].x, limiter[k][j][i].x
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div1[k][j][i].x = - ucont[k][j][i].x *
                            quadraticUpwind
                            (
                                ucat[k][j][iL].x, ucat[k][j][i].x, ucat[k][j][i+1].x, ucat[k][j][iR].x,
                                ucont[k][j][i].x
                            );

                            div1[k][j][i].y = - ucont[k][j][i].x *
                            quadraticUpwind
                            (
                                ucat[k][j][iL].y, ucat[k][j][i].y, ucat[k][j][i+1].y, ucat[k][j][iR].y,
                                ucont[k][j][i].x
                            );

                            div1[k][j][i].z = - ucont[k][j][i].x *
                            quadraticUpwind
                            (
                                ucat[k][j][iL].z, ucat[k][j][i].z, ucat[k][j][i+1].z, ucat[k][j][iR].z,
                                ucont[k][j][i].x
                            );
                        }
                    }
                }

                PetscReal nuEff, nu = cst->nu, nut;

                // viscous terms
                if
                (
                    ueqn->access->flags->isLesActive
                )
                {
                    if(isIBMCell(k, j, i, nvert))
                    {
                        nut = lnu_t[k][j][i+1];
                    }
                    else if (isIBMCell(k, j, i+1, nvert))
                    {
                        nut = lnu_t[k][j][i];
                    }
                    else
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
                    }

                    // wall model i-left patch
                    if
                    (
                            (mesh->boundaryU.iLeft=="velocityWallFunction"  && i==0) ||
                            (mesh->boundaryU.iRight=="velocityWallFunction" && i==mx-2)
                    )
                    {
                        PetscReal signTau =  1.0;
                        if(i==0)  signTau = -1.0;

                        visc1[k][j][i].x = signTau * ueqn->iLWM->tauWall.x[k][j];
                        visc1[k][j][i].y = signTau * ueqn->iLWM->tauWall.y[k][j];
                        visc1[k][j][i].z = signTau * ueqn->iLWM->tauWall.z[k][j];

                        nut = 0.0;
                    }

                    // slip boundary condition on U (set nuEff to 0)
                    if
                    (
                        (mesh->boundaryU.iLeft =="slip" && i==0   ) ||
                        (mesh->boundaryU.iRight=="slip" && i==mx-2)
                    )
                    {
                        nuEff = 0.0;
                    }
                    else
                    {
                        nuEff = nu + nut;
                    }
                }
                else
                {
                    nuEff = nu;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nuEff);
                visc1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nuEff);
                visc1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nuEff);
            }
        }
    }

    // j direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || k==0) continue;

                // get 1/V at the j-face
                PetscReal ajc = jaj[k][j][i];

                // get face normals
                csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_j (mesh, i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                // compute metric tensor
                g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
                g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
                g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

                // compute cartesian velocity derivatives w.r.t. cartesian coords
                r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

                Compute_du_dxyz
                (
                    mesh,
                    csi0,   csi1,   csi2,   eta0,   eta1,   eta2,   zet0,   zet1,   zet2,   ajc,
                    dudc,   dvdc,   dwdc,   dude,   dvde,   dwde,   dudz,   dvdz,   dwdz,
                    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                PetscInt    jL, jR;
                PetscReal denom;
                getFace2Cell4StencilEta(mesh, j, my, &jL, &jR, &denom);

                if( isIBMJFace(k, jL, i, jR, nvert) )
                {
                    jL = j;
                    jR = j+1;
                }

                if(ueqn->centralDiv)
                {
                    jL = j, jR=j+1;
                    denom=1.;
                }

                PetscReal ucon = ucont[k][j][i].y;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                // test: IBM active and not on the j-right boundary face
                if (j!=my-2 && isIBMCell(k, j, i, nvert) )
                {
                  div2[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k][j+2][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) +  ucat[k][j+1][i].x) +
                  -up * (0.125 * (-    ucat[k][j  ][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) +  ucat[k][j][i  ].x);

                  div2[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k][j+2][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +
                  -up * (0.125 * (-    ucat[k][j  ][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) +  ucat[k][j][i  ].y);

                  div2[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k][j+2][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +
                  -up * (0.125 * (-    ucat[k][j  ][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) + ucat[k][j][i  ].z);
                }
                // test: IBM active and not on the j-left boundary face
                else if (j!=0 && isIBMCell(k, j+1, i, nvert) )
                {
                  div2[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k][j+1][i].x -  2. * ucat[k][j+1][i].x +  3. * ucat[k][j  ][i].x) + ucat[k][j+1][i].x) +
                  -up * (0.125 * (-    ucat[k][j-1][i].x -  2. * ucat[k][j  ][i].x +  3. * ucat[k][j+1][i].x) + ucat[k][j][i  ].x);

                  div2[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k][j+1][i].y -  2. * ucat[k][j+1][i].y +  3. * ucat[k][j  ][i].y) + ucat[k][j+1][i].y) +
                  -up * (0.125 * (-    ucat[k][j-1][i].y -  2. * ucat[k][j  ][i].y +  3. * ucat[k][j+1][i].y) + ucat[k][j][i  ].y);

                  div2[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k][j+1][i].z -  2. * ucat[k][j+1][i].z +  3. * ucat[k][j  ][i].z) +  ucat[k][j+1][i].z) +
                  -up * (0.125 * (-    ucat[k][j-1][i].z -  2. * ucat[k][j  ][i].z +  3. * ucat[k][j+1][i].z) +  ucat[k][j][i  ].z);
                }
                // not IBM
                else
                {
                    // test: inviscid flow or weno3Div
                    if( ueqn->inviscid || ueqn->weno3Div)
                    {
                        div2[k][j][i].x = - ucont[k][j][i].y *
                        weno3
                        (
                            ucat[k][jL][i].x, ucat[k][j][i].x, ucat[k][j+1][i].x, ucat[k][jR][i].x,
                            ucont[k][j][i].y
                        );

                        div2[k][j][i].y = - ucont[k][j][i].y *
                        weno3
                        (
                            ucat[k][jL][i].y, ucat[k][j][i].y, ucat[k][j+1][i].y, ucat[k][jR][i].y,
                            ucont[k][j][i].y
                        );

                        div2[k][j][i].z = - ucont[k][j][i].y *
                        weno3
                        (
                            ucat[k][jL][i].z, ucat[k][j][i].z, ucat[k][j+1][i].z, ucat[k][jR][i].z,
                            ucont[k][j][i].y
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            div2[k][j][i].x = - ucont[k][j][i].y *
                            central
                            (
                                ucat[k][j][i].x, ucat[k][j+1][i].x
                            );

                            div2[k][j][i].y = - ucont[k][j][i].y *

                            central
                            (
                                ucat[k][j][i].y, ucat[k][j+1][i].y
                            );

                            div2[k][j][i].z = - ucont[k][j][i].y *
                            central
                            (
                                ucat[k][j][i].z, ucat[k][j+1][i].z
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div2[k][j][i].x = - ucont[k][j][i].y *
                            centralUpwind
                            (
                                ucat[k][jL][i].x, ucat[k][j][i].x, ucat[k][j+1][i].x, ucat[k][jR][i].x,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );

                            div2[k][j][i].y = - ucont[k][j][i].y *
                            centralUpwind
                            (
                                ucat[k][jL][i].y, ucat[k][j][i].y, ucat[k][j+1][i].y, ucat[k][jR][i].y,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );

                            div2[k][j][i].z = - ucont[k][j][i].y *
                            centralUpwind
                            (
                                ucat[k][jL][i].z, ucat[k][j][i].z, ucat[k][j+1][i].z, ucat[k][jR][i].z,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][jL ][i] * nMag(eta[k][jL ][i]));
                            PetscReal d1 = 1.0 / (aj[k][j  ][i] * nMag(eta[k][j  ][i]));
                            PetscReal d2 = 1.0 / (aj[k][j+1][i] * nMag(eta[k][j+1][i]));
                            PetscReal d3 = 1.0 / (aj[k][jR ][i] * nMag(eta[k][jR ][i]));

                            div2[k][j][i].x = - ucont[k][j][i].y *
                            wCentralUpwind
                            (
                                ucat[k][jL][i].x, ucat[k][j][i].x, ucat[k][j+1][i].x, ucat[k][jR][i].x,
                                d0, d1, d2, d3,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );

                            div2[k][j][i].y = - ucont[k][j][i].y *
                            wCentralUpwind
                            (
                                ucat[k][jL][i].y, ucat[k][j][i].y, ucat[k][j+1][i].y, ucat[k][jR][i].y,
                                d0, d1, d2, d3,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );

                            div2[k][j][i].z = - ucont[k][j][i].y *
                            wCentralUpwind
                            (
                                ucat[k][jL][i].z, ucat[k][j][i].z, ucat[k][j+1][i].z, ucat[k][jR][i].z,
                                d0, d1, d2, d3,
                                ucont[k][j][i].y, limiter[k][j][i].y
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div2[k][j][i].x = - ucont[k][j][i].y *
                            quadraticUpwind
                            (
                                ucat[k][jL][i].x, ucat[k][j][i].x, ucat[k][j+1][i].x, ucat[k][jR][i].x,
                                ucont[k][j][i].y
                            );

                            div2[k][j][i].y = - ucont[k][j][i].y *
                            quadraticUpwind
                            (
                                ucat[k][jL][i].y, ucat[k][j][i].y, ucat[k][j+1][i].y, ucat[k][jR][i].y,
                                ucont[k][j][i].y
                            );

                            div2[k][j][i].z = - ucont[k][j][i].y *
                            quadraticUpwind
                            (
                                ucat[k][jL][i].z, ucat[k][j][i].z, ucat[k][j+1][i].z, ucat[k][jR][i].z,
                                ucont[k][j][i].y
                            );
                        }
                    }
                }

                PetscReal nuEff, nu = cst->nu, nut;

                if
                (
                    ueqn->access->flags->isLesActive
                )
                {
                    if(isIBMCell(k, j, i, nvert))
                    {
                        nut = lnu_t[k][j+1][i];
                    }
                    else if (isIBMCell(k, j+1, i, nvert))
                    {
                        nut = lnu_t[k][j][i];
                    }
                    else
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
                    }

                    // wall model j-left patch
                    if
                    (
                            (mesh->boundaryU.jLeft=="velocityWallFunction" && j==0) ||
                            (mesh->boundaryU.jRight=="velocityWallFunction" && j==my-2)
                    )
                    {
                        PetscReal signTau =  1.0;
                        if(j==0)  signTau = -1.0;

                        visc2[k][j][i].x = signTau * ueqn->jLWM->tauWall.x[k][i];
                        visc2[k][j][i].y = signTau * ueqn->jLWM->tauWall.y[k][i];
                        visc2[k][j][i].z = signTau * ueqn->jLWM->tauWall.z[k][i];

                        nut = 0.0;
                    }

                    // slip boundary condition on U (set nuEff to 0)
                    if
                    (
                        (mesh->boundaryU.jLeft =="slip" && j==0   ) ||
                        (mesh->boundaryU.jRight=="slip" && j==my-2)
                    )
                    {
                        nuEff = 0.0;
                    }
                    else
                    {
                        nuEff = nu + nut;
                    }
                }
                else
                {
                    nuEff = nu;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nuEff);
                visc2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nuEff);
                visc2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nuEff);
            }
        }
    }

    // k direction faces are from 0 to mz-2
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || j==0) continue;

                // get 1/V at the j-face
                PetscReal ajc = kaj[k][j][i];

                // get face normals
                csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_k (mesh, i, j, k, mx, my, mz, ucat, nvert, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                // compute metric tensor
                g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                // compute cartesian velocity derivatives w.r.t. cartesian coords
                r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

                Compute_du_dxyz
                (
                    mesh,
                    csi0,   csi1,   csi2,   eta0,   eta1,   eta2,   zet0,   zet1,   zet2,   ajc,
                    dudc,   dvdc,   dwdc,   dude,   dvde,   dwde,   dudz,   dvdz,   dwdz,
                    &du_dx, &dv_dx, &dw_dx, &du_dy, &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                PetscInt    kL, kR;
                PetscReal denom;
                getFace2Cell4StencilZet(mesh, k, mz, &kL, &kR, &denom);

                if( isIBMKFace(kL, j, i, kR, nvert) )
                {
                    kL = k;
                    kR=k+1;
                    denom=1.;
                }

                if(ueqn->centralDiv)
                {
                    kL = k, kR=k+1;
                    denom=1.;
                }

                PetscReal ucon = ucont[k][j][i].z;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                // test: IBM active and not on the k-right boundary face
                if (k!=mz-2 && isIBMCell(k, j, i, nvert))
                {
                  div3[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k+2][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +
                  -up * (0.125 * (-    ucat[k  ][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);

                  div3[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k+2][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) + ucat[k+1][j][i].y)   +
                  -up * (0.125 * (-    ucat[k  ][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) +  ucat[k][j][i  ].y);

                  div3[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k+2][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +
                  -up * (0.125 * (-    ucat[k  ][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
                }
                // test: IBM active and not on the k-left boundary face
                else if (k!=0 && isIBMCell(k+1, j, i, nvert) )
                {
                  div3[k][j][i].x
                  =
                  -um * (0.125 * (-    ucat[k+1][j][i].x -  2. * ucat[k+1][j][i].x +  3. * ucat[k  ][j][i].x) + ucat[k+1][j][i].x)   +
                  -up * (0.125 * (-    ucat[k-1][j][i].x -  2. * ucat[k  ][j][i].x +  3. * ucat[k+1][j][i].x) +  ucat[k][j][i  ].x);

                  div3[k][j][i].y
                  =
                  -um * (0.125 * (-    ucat[k+1][j][i].y -  2. * ucat[k+1][j][i].y +  3. * ucat[k  ][j][i].y) +  ucat[k+1][j][i].y)   +
                  -up * (0.125 * (-    ucat[k-1][j][i].y -  2. * ucat[k  ][j][i].y +  3. * ucat[k+1][j][i].y) + ucat[k][j][i  ].y);

                  div3[k][j][i].z
                  =
                  -um * (0.125 * (-    ucat[k+1][j][i].z -  2. * ucat[k+1][j][i].z +  3. * ucat[k  ][j][i].z) + ucat[k+1][j][i].z)   +
                  -up * (0.125 * (-    ucat[k-1][j][i].z -  2. * ucat[k  ][j][i].z +  3. * ucat[k+1][j][i].z) + ucat[k][j][i  ].z);
                }
                // not IBM
                else
                {
                    // test: inviscid flow or weno3Div
                    if(ueqn->inviscid || ueqn->weno3Div)
                    {
                        div3[k][j][i].x = - ucont[k][j][i].z *
                        weno3
                        (
                            ucat[kL][j][i].x, ucat[k][j][i].x, ucat[k+1][j][i].x, ucat[kR][j][i].x,
                            ucont[k][j][i].z
                        );

                        div3[k][j][i].y = - ucont[k][j][i].z *
                        weno3
                        (
                            ucat[kL][j][i].y, ucat[k][j][i].y, ucat[k+1][j][i].y, ucat[kR][j][i].y,
                            ucont[k][j][i].z
                        );

                        div3[k][j][i].z = - ucont[k][j][i].z *
                        weno3
                        (
                            ucat[kL][j][i].z, ucat[k][j][i].z, ucat[k+1][j][i].z, ucat[kR][j][i].z,
                            ucont[k][j][i].z
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
<<<<<<< HEAD
                            div3[k][j][i].x = - ucont[k][j][i].z  *
=======
                            div3[k][j][i].x = - ucont[k][j][i].z *
>>>>>>> fixed bug in rayleigh damping, added pTilde buoyancy formulation, added hydrostatic assumption
                            central
                            (
                                ucat[k][j][i].x, ucat[k+1][j][i].x
                            );

                            div3[k][j][i].y = - ucont[k][j][i].z *
                            central
                            (
                                ucat[k][j][i].y, ucat[k+1][j][i].y
                            );

                            div3[k][j][i].z = - ucont[k][j][i].z *
                            central
                            (
                                ucat[k][j][i].z, ucat[k+1][j][i].z
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div3[k][j][i].x = - ucont[k][j][i].z *
                            centralUpwind
                            (
                                ucat[kL][j][i].x, ucat[k][j][i].x, ucat[k+1][j][i].x, ucat[kR][j][i].x,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );

                            div3[k][j][i].y = - ucont[k][j][i].z *
                            centralUpwind
                            (
                                ucat[kL][j][i].y, ucat[k][j][i].y, ucat[k+1][j][i].y, ucat[kR][j][i].y,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );

                            div3[k][j][i].z = - ucont[k][j][i].z *
                            centralUpwind
                            (
                                ucat[kL][j][i].z, ucat[k][j][i].z, ucat[k+1][j][i].z, ucat[kR][j][i].z,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[kL ][j][i] * nMag(zet[kL ][j][i]));
                            PetscReal d1 = 1.0 / (aj[k  ][j][i] * nMag(zet[k  ][j][i]));
                            PetscReal d2 = 1.0 / (aj[k+1][j][i] * nMag(zet[k+1][j][i]));
                            PetscReal d3 = 1.0 / (aj[kR ][j][i] * nMag(zet[kR ][j][i]));

                            div3[k][j][i].x = - ucont[k][j][i].z *
                            wCentralUpwind
                            (
                                ucat[kL][j][i].x, ucat[k][j][i].x, ucat[k+1][j][i].x, ucat[kR][j][i].x,
                                d0, d1, d2, d3,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );

                            div3[k][j][i].y = - ucont[k][j][i].z *
                            wCentralUpwind
                            (
                                ucat[kL][j][i].y, ucat[k][j][i].y, ucat[k+1][j][i].y, ucat[kR][j][i].y,
                                d0, d1, d2, d3,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );

                            div3[k][j][i].z = - ucont[k][j][i].z *
                            wCentralUpwind
                            (
                                ucat[kL][j][i].z, ucat[k][j][i].z, ucat[k+1][j][i].z, ucat[kR][j][i].z,
                                d0, d1, d2, d3,
                                ucont[k][j][i].z, limiter[k][j][i].z
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div3[k][j][i].x = - ucont[k][j][i].z *
                            quadraticUpwind
                            (
                                ucat[kL][j][i].x, ucat[k][j][i].x, ucat[k+1][j][i].x, ucat[kR][j][i].x,
                                ucont[k][j][i].z
                            );

                            div3[k][j][i].y = - ucont[k][j][i].z *
                            quadraticUpwind
                            (
                                ucat[kL][j][i].y, ucat[k][j][i].y, ucat[k+1][j][i].y, ucat[kR][j][i].y,
                                ucont[k][j][i].z
                            );

                            div3[k][j][i].z = - ucont[k][j][i].z *
                            quadraticUpwind
                            (
                                ucat[kL][j][i].z, ucat[k][j][i].z, ucat[k+1][j][i].z, ucat[kR][j][i].z,
                                ucont[k][j][i].z
                            );
                        }
                    }
                }

                PetscReal nuEff, nu = cst->nu, nut;

                if
                (
                    ueqn->access->flags->isLesActive
                )
                {
                    if(isIBMCell(k, j, i, nvert))
                    {
                        nut = lnu_t[k+1][j][i];
                    }
                    else if (isIBMCell(k+1, j, i, nvert))
                    {
                        nut = lnu_t[k][j][i];
                    }
                    else
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
                    }

                    // wall models not present in k direction

                    // slip boundary condition on U (set nuEff to 0)
                    if
                    (
                        (mesh->boundaryU.kLeft =="slip" && k==0   ) ||
                        (mesh->boundaryU.kRight=="slip" && k==mz-2)
                    )
                    {
                        nuEff = 0.0;
                    }
                    else
                    {
                        nuEff = nu + nut;
                    }
                }
                else
                {
                    nuEff = nu;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nuEff);
                visc3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nuEff);
                visc3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nuEff);
            }
        }
    }

    // restore distributed arrays and scatter local to local
    DMDAVecRestoreArray(fda, ueqn->lDiv1, &div1);
    DMDAVecRestoreArray(fda, ueqn->lDiv2, &div2);
    DMDAVecRestoreArray(fda, ueqn->lDiv3, &div3);

    DMLocalToLocalBegin(fda, ueqn->lDiv1, INSERT_VALUES, ueqn->lDiv1);
    DMLocalToLocalEnd  (fda, ueqn->lDiv1, INSERT_VALUES, ueqn->lDiv1);
    DMLocalToLocalBegin(fda, ueqn->lDiv2, INSERT_VALUES, ueqn->lDiv2);
    DMLocalToLocalEnd  (fda, ueqn->lDiv2, INSERT_VALUES, ueqn->lDiv2);
    DMLocalToLocalBegin(fda, ueqn->lDiv3, INSERT_VALUES, ueqn->lDiv3);
    DMLocalToLocalEnd  (fda, ueqn->lDiv3, INSERT_VALUES, ueqn->lDiv3);

    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);

    DMLocalToLocalBegin(fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalEnd  (fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalBegin(fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalEnd  (fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalBegin(fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);
    DMLocalToLocalEnd  (fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);

    // ---------------------------------------------------------------------- //
    // FORM THE RIGHT HAND SIDE
    // ---------------------------------------------------------------------- //

    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv1,   ueqn->lDiv1,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv2,   ueqn->lDiv1,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv3,   ueqn->lDiv1,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc1,  ueqn->lVisc1,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc2,  ueqn->lVisc1,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc3,  ueqn->lVisc1,  "localToLocal");

    // get the arrays
    DMDAVecGetArray(fda, ueqn->lDiv1, &div1);
    DMDAVecGetArray(fda, ueqn->lDiv2, &div2);
    DMDAVecGetArray(fda, ueqn->lDiv3, &div3);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);

    DMDAVecGetArray(fda, ueqn->lFp, &fp);

    // loop over cells and compute fp (viscous + divergence contributions)
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // divergence contribution
                fp[k][j][i].x
                +=
                (
                    div1[k][j][i].x - div1[k][j][i-1].x +
                    div2[k][j][i].x - div2[k][j-1][i].x +
                    div3[k][j][i].x - div3[k-1][j][i].x
                );

                fp[k][j][i].y
                +=
                (
                    div1[k][j][i].y - div1[k][j][i-1].y +
                    div2[k][j][i].y - div2[k][j-1][i].y +
                    div3[k][j][i].y - div3[k-1][j][i].y
                );

                fp[k][j][i].z
                +=
                (
                    div1[k][j][i].z - div1[k][j][i-1].z +
                    div2[k][j][i].z - div2[k][j-1][i].z +
                    div3[k][j][i].z - div3[k-1][j][i].z
                );

                // viscous contribution
                if(ueqn->inviscid)
                {

                }
                else
                {
                    fp[k][j][i].x
                    +=
                    (
                        visc1[k][j][i].x - visc1[k][j][i-1].x +
                        visc2[k][j][i].x - visc2[k][j-1][i].x +
                        visc3[k][j][i].x - visc3[k-1][j][i].x
                    );

                    fp[k][j][i].y
                    +=
                    (
                        visc1[k][j][i].y - visc1[k][j][i-1].y +
                        visc2[k][j][i].y - visc2[k][j-1][i].y +
                        visc3[k][j][i].y - visc3[k-1][j][i].y
                    );

                    fp[k][j][i].z
                    +=
                    (
                        visc1[k][j][i].z - visc1[k][j][i-1].z +
                        visc2[k][j][i].z - visc2[k][j-1][i].z +
                        visc3[k][j][i].z - visc3[k-1][j][i].z
                    );

                    if( (i == 1) && !(mesh->i_periodic) && !(mesh->ii_periodic))
                    {
                        fp[k][j][i-1].x = fp[k][j][i].x;
                        fp[k][j][i-1].y = fp[k][j][i].y;
                        fp[k][j][i-1].z = fp[k][j][i].z;
                    }

                    if( (j == 1) && !(mesh->j_periodic) && !(mesh->jj_periodic))
                    {
                        fp[k][j-1][i].x = fp[k][j][i].x;
                        fp[k][j-1][i].y = fp[k][j][i].y;
                        fp[k][j-1][i].z = fp[k][j][i].z;
                    }

                    if( (k == 1) && !(mesh->k_periodic) && !(mesh->kk_periodic))
                    {
                        fp[k-1][j][i].x = fp[k][j][i].x;
                        fp[k-1][j][i].y = fp[k][j][i].y;
                        fp[k-1][j][i].z = fp[k][j][i].z;
                    }

                    if( (i == mx-2) && !(mesh->i_periodic) && !(mesh->ii_periodic))
                    {
                        fp[k][j][i+1].x = fp[k][j][i].x;
                        fp[k][j][i+1].y = fp[k][j][i].y;
                        fp[k][j][i+1].z = fp[k][j][i].z;
                    }

                    if( (j == my-2) && !(mesh->j_periodic) && !(mesh->jj_periodic))
                    {
                        fp[k][j+1][i].x = fp[k][j][i].x;
                        fp[k][j+1][i].y = fp[k][j][i].y;
                        fp[k][j+1][i].z = fp[k][j][i].z;
                    }

                    if((k == mz-2) && !(mesh->k_periodic) && !(mesh->kk_periodic))
                    {
                        fp[k+1][j][i].x = fp[k][j][i].x;
                        fp[k+1][j][i].y = fp[k][j][i].y;
                        fp[k+1][j][i].z = fp[k][j][i].z;
                    }

                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lFp, &fp);

    DMLocalToLocalBegin(fda, ueqn->lFp, INSERT_VALUES, ueqn->lFp);
    DMLocalToLocalEnd  (fda, ueqn->lFp, INSERT_VALUES, ueqn->lFp);

    resetCellPeriodicFluxes(mesh, ueqn->lFp, ueqn->lFp, "vector", "localToLocal");

    DMDAVecGetArray(fda, ueqn->lFp, &fp);

    // interpolate the convective and viscous contributions at the surface
    // centers. This is done by projecting the flux fp along the curvilinear
    // coordinates and averaging between two adhiacent cells center values to
    // get the face value
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // get the stencils
                PetscReal denom;
                PetscInt    iL, iR;
                PetscInt    jL, jR;
                PetscInt    kL, kR;
                getFace2Cell4StencilCsi(mesh, i, mx, &iL, &iR, &denom);
                getFace2Cell4StencilEta(mesh, j, my, &jL, &jR, &denom);
                getFace2Cell4StencilZet(mesh, k, mz, &kL, &kR, &denom);

                PetscReal v0, v1, v2, v3;
                PetscReal d0, d1, d2, d3;

                v0 = csi[k][j][iL ].x * fp[k][j][iL ].x +
                     csi[k][j][iL ].y * fp[k][j][iL ].y +
                     csi[k][j][iL ].z * fp[k][j][iL ].z;
                v1 = csi[k][j][i  ].x * fp[k][j][i  ].x +
                     csi[k][j][i  ].y * fp[k][j][i  ].y +
                     csi[k][j][i  ].z * fp[k][j][i  ].z;
                v2 = csi[k][j][i+1].x * fp[k][j][i+1].x +
                     csi[k][j][i+1].y * fp[k][j][i+1].y +
                     csi[k][j][i+1].z * fp[k][j][i+1].z;
                v3 = csi[k][j][iR ].x * fp[k][j][iR ].x +
                     csi[k][j][iR ].y * fp[k][j][iR ].y +
                     csi[k][j][iR ].z * fp[k][j][iR ].z;

                // compute cell widths
                d0 = 1.0 / (aj[k][j][iL ] * nMag(csi[k][j][iL  ]));
                d1 = 1.0 / (aj[k][j][i  ] * nMag(csi[k][j][i  ]));
                d2 = 1.0 / (aj[k][j][i+1] * nMag(csi[k][j][i+1]));
                d3 = 1.0 / (aj[k][j][iR ] * nMag(csi[k][j][iR ]));

                rhs[k][j][i].x += scale * iaj[k][j][i] * central(v1, v2);

                v0 = eta[k][jL ][i].x * fp[k][jL ][i].x +
                     eta[k][jL ][i].y * fp[k][jL ][i].y +
                     eta[k][jL ][i].z * fp[k][jL ][i].z;
                v1 = eta[k][j  ][i].x * fp[k][j  ][i].x +
                     eta[k][j  ][i].y * fp[k][j  ][i].y +
                     eta[k][j  ][i].z * fp[k][j  ][i].z;
                v2 = eta[k][j+1][i].x * fp[k][j+1][i].x +
                     eta[k][j+1][i].y * fp[k][j+1][i].y +
                     eta[k][j+1][i].z * fp[k][j+1][i].z;
                v3 = eta[k][jR ][i].x * fp[k][jR ][i].x +
                     eta[k][jR ][i].y * fp[k][jR ][i].y +
                     eta[k][jR ][i].z * fp[k][jR ][i].z;

                // compute cell widths
                d0 = 1.0 / (aj[k][jL ][i] * nMag(eta[k][jL ][i]));
                d1 = 1.0 / (aj[k][j  ][i] * nMag(eta[k][j  ][i]));
                d2 = 1.0 / (aj[k][j+1][i] * nMag(eta[k][j+1][i]));
                d3 = 1.0 / (aj[k][jR ][i] * nMag(eta[k][jR ][i]));

                rhs[k][j][i].y += scale * jaj[k][j][i] * central(v1, v2);

                v0 = zet[kL ][j][i].x * fp[kL ][j][i].x +
                     zet[kL ][j][i].y * fp[kL ][j][i].y +
                     zet[kL ][j][i].z * fp[kL ][j][i].z;
                v1 = zet[k  ][j][i].x * fp[k  ][j][i].x +
                     zet[k  ][j][i].y * fp[k  ][j][i].y +
                     zet[k  ][j][i].z * fp[k  ][j][i].z;
                v2 = zet[k+1][j][i].x * fp[k+1][j][i].x +
                     zet[k+1][j][i].y * fp[k+1][j][i].y +
                     zet[k+1][j][i].z * fp[k+1][j][i].z;
                v3 = zet[kR ][j][i].x * fp[kR ][j][i].x +
                     zet[kR ][j][i].y * fp[kR ][j][i].y +
                     zet[kR ][j][i].z * fp[kR ][j][i].z;

                // compute cell widths
                d0 = 1.0 / (aj[kL ][j][i] * nMag(zet[kL ][j][i]));
                d1 = 1.0 / (aj[k  ][j][i] * nMag(zet[k  ][j][i]));
                d2 = 1.0 / (aj[k+1][j][i] * nMag(zet[k+1][j][i]));
                d3 = 1.0 / (aj[kR ][j][i] * nMag(zet[kR ][j][i]));

                rhs[k][j][i].z += scale * kaj[k][j][i] * central(v1, v2);

                /*
                rhs[k][j][i].z
                +=
                scale *
                (
                      0.5 *
                    (
                        zet[k][j][i].x * fp[k][j][i].x +
                        zet[k][j][i].y * fp[k][j][i].y +
                        zet[k][j][i].z * fp[k][j][i].z
                    )
                    + 0.5 *
                    (
                        zet[k+1][j][i].x * fp[k+1][j][i].x +
                        zet[k+1][j][i].y * fp[k+1][j][i].y +
                        zet[k+1][j][i].z * fp[k+1][j][i].z
                    )
                ) * kaj[k][j][i];
                */

                if
                (
                    isIBMIFace(k, j, i, i+1, nvert)
                )
                {
                    rhs[k][j][i].x = 0;
                }

                if
                (
                    isIBMJFace(k, j, i, j+1, nvert)
                )
                {
                    rhs[k][j][i].y = 0;
                }

                if
                (
                    isIBMKFace(k, j, i, k+1, nvert)
                )
                {
                    rhs[k][j][i].z = 0;
                }
            }
        }
    }

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    DMDAVecRestoreArray(fda, ueqn->lFp, &fp);
    DMDAVecRestoreArray(fda, ueqn->lDiv1, &div1);
    DMDAVecRestoreArray(fda, ueqn->lDiv2, &div2);
    DMDAVecRestoreArray(fda, ueqn->lDiv3, &div3);

    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, Rhs,  &rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnSNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
{
    ueqn_  *ueqn  = (ueqn_*)ptr;
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    VecCopy(Ucont, ueqn->Ucont);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset contravariant periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // transform to cartesian
    contravariantToCartesian(ueqn);

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // update wall model (optional)
    UpdateWallModelsU(ueqn);

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    PetscReal dt = clock->dt;

    // add pressure gradient term
    VecAXPY(Rhs, -1.0, ueqn->dP);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        Coriolis(ueqn, Rhs, 1.0);
    }

    // add side force term
    if(ueqn->access->flags->isSideForceActive)
    {
        SideForce(ueqn, Rhs, 1.0);
    }

    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;

        // add buoyancy gradient term
        if(teqn->pTildeFormulation)
        {
            VecAXPY(Rhs, -1.0, teqn->ghGradRhok);
        }
        else
        {
            Buoyancy(ueqn, Rhs, 1.0);
        }
    }

    // add viscous and transport terms
    if(clock->it > clock->itStart)
    {
        VecAXPY(Rhs, 0.5, ueqn->Rhs_o);
        FormU(ueqn, Rhs, 0.5);
    }
    else
    {
        FormU(ueqn, Rhs, 1.0);
    }

    // add wind farm model
    if(ueqn->access->flags->isWindFarmActive)
    {
        VecAXPY(Rhs, 1.0, ueqn->access->farm->sourceFarmCont);
    }

    if
    (
        ueqn->access->flags->isXDampingActive ||
        ueqn->access->flags->isZDampingActive
    )
    {
        dampingSourceU(ueqn, Rhs, 1.0);
    }


    // multiply for dt
    VecScale(Rhs, dt);

    // add driving source terms after as it is not scaled by 1/dt
    if(ueqn->access->flags->isAblActive)
    {
        sourceU(ueqn, Rhs, 1.0);
    }

    resetNonResolvedCellFaces(mesh, Rhs);

    // add time derivative term

    VecAXPY(Rhs, -1.0, Ucont);
    VecAXPY(Rhs,  1.0, ueqn->Ucont_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsU(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset contravariant periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // transform to cartesian
    contravariantToCartesian(ueqn);

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // update wall model (optional)
    UpdateWallModelsU(ueqn);

    // initialize the rhs vector
    VecSet(ueqn->Rhs, 0.0);

    // add pressure gradient term
    VecAXPY(ueqn->Rhs, -1.0, ueqn->dP);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        Coriolis(ueqn, ueqn->Rhs, 1.0);
    }

    // add side force term
    if(ueqn->access->flags->isSideForceActive)
    {
        SideForce(ueqn, ueqn->Rhs, 1.0);
    }

    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;

        // add buoyancy gradient term
        if(teqn->pTildeFormulation)
        {
            VecAXPY(ueqn->Rhs, -1.0, teqn->ghGradRhok);
        }
        else
        {
            Buoyancy(ueqn, ueqn->Rhs, 1.0);
        }
    }

    // add viscous and transport terms
    FormU(ueqn, ueqn->Rhs, 1.0);

    // add wind farm model
    if(ueqn->access->flags->isWindFarmActive)
    {
        VecAXPY(ueqn->Rhs, 1.0, ueqn->access->farm->sourceFarmCont);
    }

    if
    (
        ueqn->access->flags->isXDampingActive ||
        ueqn->access->flags->isZDampingActive
    )
    {
        dampingSourceU(ueqn, ueqn->Rhs, 1.0);
    }

    // source term is pre-scaled by dt, so here we have to divide it again since
    // the rhs is not multiplied by dt in this function
    if(ueqn->access->flags->isAblActive)
    {
        sourceU(ueqn, ueqn->Rhs, 1.0 / clock->dt);
    }

    resetNonResolvedCellFaces(mesh, ueqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnEuler(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    // compute the right hand side
    FormExplicitRhsU(ueqn);

    // multiply for dt
    VecScale(ueqn->Rhs, clock->dt);

    // add time derivative term
    VecAXPY(ueqn->Rhs, 1, ueqn->Ucont);

    // store the solution
    VecCopy(ueqn->Rhs,  ueqn->Ucont);

    // scatter to local values
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnRK4(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for Ucont, Stage ");

    PetscInt  s = 4;
    PetscReal b[4];
    PetscReal a[4];

    b[0] = 1.0 / 6.0;
    b[1] = 1.0 / 3.0;
    b[2] = 1.0 / 3.0;
    b[3] = 1.0 / 6.0;

    a[0] = 0.0;
    a[1] = 0.5;
    a[2] = 0.5;
    a[3] = 1.0;

    PetscReal dt = clock->dt;

    // Ucont_o contribution
    VecCopy(ueqn->Ucont_o, ueqn->Utmp);

    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate U guess and evaluate RHS
        if(i!=0)
        {
            VecWAXPY(ueqn->Ucont, a[i] * dt, ueqn->Rhs, ueqn->Ucont_o);
            DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
            DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
        }

        // compute function guess
        FormExplicitRhsU(ueqn);

        // add contribution from K1, K2, K3, K4
        VecAXPY(ueqn->Utmp, dt * b[i], ueqn->Rhs);
    }

    VecCopy(ueqn->Utmp, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveUEqn(ueqn_ *ueqn)
{
    mesh_          *mesh = ueqn->access->mesh;
    clock_         *clock = ueqn->access->clock;

    // set the right hand side to zero
    VecSet (ueqn->Rhs, 0.0);

    if(ueqn->ddtScheme == "backwardEuler")
    {
        PetscReal     norm;
        PetscInt      iter;
        PetscReal     ts,te;

        PetscTime(&ts);
        PetscPrintf(mesh->MESH_COMM, "TRSNES: Solving for Ucont, Initial residual = ");

        // copy the solution in the tmp variable
        VecCopy(ueqn->Ucont, ueqn->Utmp);

        // solve momentum equation and compute solution norm and iteration number
        SNESSolve(ueqn->snesU, PETSC_NULL, ueqn->Utmp);
        SNESGetFunctionNorm(ueqn->snesU, &norm);
        SNESGetIterationNumber(ueqn->snesU, &iter);

        // store the solution and global to local scatter
        VecCopy(ueqn->Utmp, ueqn->Ucont);

        // scatter
        DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
        DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

        PetscTime(&te);
        PetscPrintf(mesh->MESH_COMM,"Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n", norm, iter, te-ts);
    }
    else if (ueqn->ddtScheme == "forwardEuler")
    {
        UeqnEuler(ueqn);
    }
    else if (ueqn->ddtScheme == "rungeKutta4")
    {
        UeqnRK4(ueqn);
    }

    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // adjust inflow/outflow fluxes to ensure mass conservation
    if(ueqn->access->flags->isOversetActive && *(ueqn->access->domainID) != 0)
    {
        adjustFluxesOverset(ueqn);
    }
    else
    {
        adjustFluxes(ueqn);
    }

    return(0);
}
