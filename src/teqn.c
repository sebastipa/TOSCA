//! \file  teqn.c
//! \brief Contains T equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

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

        // initialize p-tilde formulation of temperature equation
        teqn->pTildeFormulation = 0;
        PetscOptionsGetInt(PETSC_NULL, PETSC_NULL,  "-pTildeBuoyancy", &(teqn->pTildeFormulation),   PETSC_NULL);

        VecDuplicate(mesh->Nvert, &(teqn->TmprtTmp));    VecSet(teqn->TmprtTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt));       VecSet(teqn->Tmprt,    0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt_o));     VecSet(teqn->Tmprt_o,  0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Rhs));         VecSet(teqn->Rhs,      0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Rhs_o));       VecSet(teqn->Rhs_o,    0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt));      VecSet(teqn->lTmprt,   0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt_o));    VecSet(teqn->lTmprt_o, 0.0);

        VecDuplicate(mesh->lCent, &(teqn->lDivT));       VecSet(teqn->lDivT,    0.0);
        VecDuplicate(mesh->lCent, &(teqn->lViscT));      VecSet(teqn->lViscT,   0.0);
        VecDuplicate(mesh->Nvert, &(teqn->sourceT));     VecSet(teqn->sourceT,  0.0);

        if(teqn->pTildeFormulation)
        {
            VecDuplicate(mesh->lAj,   &(teqn->lRhoK));       VecSet(teqn->lRhoK,   0.0);
            VecDuplicate(mesh->Cent,  &(teqn->ghGradRhok));  VecSet(teqn->ghGradRhok,   0.0);
        }

        if(flags->isIBMActive)
        {
            VecDuplicate(mesh->lCent, &(teqn->lViscIBMT));   VecSet(teqn->lViscIBMT,  0.0);
        }

        // read time discretization scheme
        readDictWord("control.dat", "-dTdtScheme", &(teqn->ddtScheme));

        // create the SNES solver
        if(teqn->ddtScheme=="backwardEuler")
        {
            // readDictDouble("control.dat", "-absTolT", &(teqn->absExitTol));
            // readDictDouble("control.dat", "-relTolT", &(teqn->relExitTol));

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
        }
        else if (teqn->ddtScheme=="rungeKutta4")
        {

        }
        else
        {
            char error[512];
            sprintf(error, "unknown ddtScheme %s for T equation, available schemes are\n    1. backwardEuler\n    2. rungeKutta4", teqn->ddtScheme.c_str());
            fatalErrorInFunction("InitializeTEqn", error);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CorrectSourceTermsT(teqn_ *teqn, PetscInt print)
{
    mesh_         *mesh = teqn->access->mesh;
    clock_        *clock= teqn->access->clock;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscReal     ***t,  ***source;
    PetscReal     ***nvert, ***aj;
    Cmpnts        ***cent;

    abl_          *abl  = teqn->access->abl;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, teqn->sourceT, &source);
    DMDAVecGetArray(da, teqn->lTmprt,  &t);

    PetscReal nLevels = my-2;

    std::vector<PetscReal> ltMean(nLevels);
    std::vector<PetscReal> gtMean(nLevels);
    std::vector<PetscReal> tD(nLevels);
    std::vector<PetscReal> tM(nLevels);

    for(j=0; j<nLevels; j++)
    {
        ltMean[j] = 0.0;
        gtMean[j] = 0.0;
    }

    DMDAVecGetArray(da, mesh->lAj, &aj);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                ltMean[j-1] += t[k][j][i] / aj[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    MPI_Allreduce(&ltMean[0], &gtMean[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    for(j=0; j<nLevels; j++)
    {
        gtMean[j] = gtMean[j] / abl->totVolPerLevel[j];
        tM[j]     = gtMean[j];
    }

    if(abl->controllerTypeT == "initial")
    {
        // save initial temperature
        if(clock->it == clock->itStart)
        {
            for(j=0; j<nLevels; j++)
            {
                abl->tDes[j] = gtMean[j];
            }
        }
    }
    else if((abl->controllerTypeT == "directProfileAssimilation") || (abl->controllerTypeT == "indirectProfileAssimilation"))
    {
        PetscInt  idxh1, idxh2, idxt1;
        PetscReal wth1, wth2, wtt1;
        PetscReal tH1, tH2;

        //find the two closest available mesoscale data in time
        PetscInt idx_1 = abl->closestTimeIndT;
        PetscInt idx_2 = abl->closestTimeIndT + 1;

        PetscInt lwrBound = 0;
        PetscInt uprBound = abl->numtT;

        if(clock->it > clock->itStart)
        {
            lwrBound = PetscMax(0, (abl->closestTimeIndT - 50));
            uprBound = PetscMin(abl->numtT, (abl->closestTimeIndT + 50));
        }

        // build error vector for the time search
        PetscReal  diff[abl->numtT];

        for(PetscInt i=lwrBound; i<uprBound; i++)
        {
            diff[i] = fabs(abl->timeT[i] - clock->time);
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
        if(abl->timeT[idx_2] < abl->timeT[idx_1])
        {
            PetscInt idx_tmp = idx_2;
            idx_2 = idx_1;
            idx_1 = idx_tmp;
        }

        // find interpolation weights
        PetscReal idx = (idx_2 - idx_1) / (abl->timeT[idx_2] - abl->timeT[idx_1]) * (clock->time - abl->timeT[idx_1]) + idx_1;
        PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
        PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

        // reset the closest index for next iteration
        abl->closestTimeIndT = idx_1;

        for(j=0; j<nLevels; j++)
        {
            idxh1 = abl->tempInterpIdx[j][0];
            idxh2 = abl->tempInterpIdx[j][1];

            wth1  = abl->tempInterpWts[j][0];
            wth2  = abl->tempInterpWts[j][1];

            tH1 = w1 * abl->tMeso[idxh1][idx_1] + (w2) * abl->tMeso[idxh1][idx_2];
            tH2 = w1 * abl->tMeso[idxh2][idx_1] + (w2) * abl->tMeso[idxh2][idx_2];
            tD[j] = wth1*tH1 + wth2*tH2;
            abl->tDes[j] = tD[j];
            // PetscPrintf(PETSC_COMM_WORLD, "tdes = %lf, gtmean = %lf, interpolated from time:%lf %lf, wts %lf %lf \n", abl->tDes[j], gtMean[j], abl->timeT[idx_1], abl->timeT[idx_2], w1, w2);
        }

        // for(j=0; j<nLevels; j++)
        // {
        //     // if cellLevels is below the lowest mesoscale data point use the source at the last available height
        //     if(abl->cellLevels[j] < abl->lowestSrcHt)
        //     {
        //         abl->tDes[j] = abl->tDes[abl->lowestIndT];
        //         gtMean[j] = gtMean[abl->lowestIndT];
        //         tD[j]     = tD[abl->lowestIndT];
        //         tM[j]     = tM[abl->lowestIndT];
        //     }

        //     if(abl->cellLevels[i] > abl->hT[abl->numhT - 1])
        //     {
        //         abl->tDes[j] = abl->tDes[abl->highestIndT];
        //         gtMean[j] = gtMean[abl->highestIndT];
        //         tD[j]     = tD[abl->highestIndT];
        //         tM[j]     = tM[abl->highestIndT];
        //     }
        // }

        if((abl->controllerTypeT == "indirectProfileAssimilation"))
        {
            //smooth the source based on the polynomial interpolation
            for (PetscInt iCtr = 0; iCtr < nLevels; iCtr++)
            {
                PetscReal dotProduct1 = 0.0, dotProduct2 = 0.0;

                for (PetscInt kCtr = 0; kCtr < nLevels; kCtr++)
                {
                    dotProduct1 += abl->polyCoeffM[iCtr][kCtr] * tD[kCtr];
                    dotProduct2 += abl->polyCoeffM[iCtr][kCtr] * tM[kCtr];
                }
                abl->tDes[iCtr] = dotProduct1;
                gtMean[iCtr]    = dotProduct2;
            }

            for(j=0; j<nLevels; j++)
            {
                if(abl->cellLevels[j] < abl->lowestSrcHt)
                {
                    abl->tDes[j] = abl->tDes[abl->lowestIndT];
                    gtMean[j] = gtMean[abl->lowestIndT];
                }

                if(abl->cellLevels[i] > abl->hT[abl->numhT - 1])
                {
                    abl->tDes[j] = abl->tDes[abl->highestIndT];
                    gtMean[j] = gtMean[abl->highestIndT];
                }
            }
        }
    }

    // compute level-wise source
    std::vector<PetscReal> lsource(lye-lys);
    for(j=lys; j<lye; j++)
    {
        lsource[j-lys] = abl->relax*(abl->tDes[j-1] - gtMean[j-1]);
    }

    // apply source and compute error
    PetscReal lrelError = 0.0, grelError = 0.0;
    PetscInt lcells = 0, gcells = 0;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                source[k][j][i]  = lsource[j-lys];
                lrelError += lsource[j-lys] / abl->tDes[j-1];
                lcells++;
            }
        }
    }

    MPI_Allreduce(&lrelError, &grelError, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lcells, &gcells, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

    DMDAVecRestoreArray(da, teqn->sourceT, &source);
    DMDAVecRestoreArray(da, teqn->lTmprt,  &t);

    if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: global avg error on theta = %.5f percent\n", grelError*100/gcells);

    // write source terms to file (radiative fluxes)
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
                fatalErrorInFunction("correctSourceTermT",  error);
            }
        }

        word fileName = "postProcessing/temperatureSource_" + getStartTimeName(clock);
        FILE *fp = fopen(fileName.c_str(), "a");

        if(!fp)
        {
            char error[512];
            sprintf(error, "cannot open file postProcessing/temperatureSource\n");
            fatalErrorInFunction("correctSourceTermT",  error);
        }
        else
        {
            PetscInt width = -15;

            if(clock->it == clock->itStart)
            {
                word w1 = "levels";
                PetscFPrintf(mesh->MESH_COMM, fp, "%*s\n", width, w1.c_str());

                for(j=0; j<nLevels; j++)
                {
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, abl->cellLevels[j]);
                }

                PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                word w2 = "time";
                PetscFPrintf(mesh->MESH_COMM, fp, "%*s\n", width, w2.c_str());
            }

            PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, clock->time);

            for(j=0; j<nLevels; j++)
            {
                PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e\t", width, abl->tDes[j] - gtMean[j]);
            }

            PetscFPrintf(mesh->MESH_COMM, fp, "\n");

            fclose(fp);
        }
    }

    std::vector<PetscReal> ().swap(ltMean);
    std::vector<PetscReal> ().swap(gtMean);
    std::vector<PetscReal> ().swap(lsource);
    std::vector<PetscReal> ().swap(tD);
    std::vector<PetscReal> ().swap(tM);
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

    PetscReal     g, tRef;

    if(teqn->access->flags->isAblActive)
    {
         g = -9.81, tRef = ueqn->access->abl->tRef;
    }
    else
    {
         g = -9.81, tRef = ueqn->access->constants->tRef;
    }


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
PetscErrorCode correctDampingSourcesT(teqn_ *teqn)
{
    mesh_         *mesh = teqn->access->mesh;
    abl_          *abl  = teqn->access->abl;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***cent;

    PetscReal     ***tP;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    if(teqn->access->flags->isYDampingActive)
    {
        cellIds     *srcMinInd = abl->srcMinInd, *srcMaxInd = abl->srcMaxInd;
        PetscReal   **tMapped = abl->tMapped;
        PetscInt    numJ, numK, numI;
        PetscInt    imin, imax, jmin, jmax, kmin, kmax;
        precursor_  *precursor;
        domain_     *pdomain;

        // update the unsteady uBar state
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
            }

            for(PetscInt p=0; p < abl->numSourceProc; p++)
            {
                if(rank == abl->sourceProcList[p])
                {

                    imin = abl->srcMinInd[p].i;
                    imax = abl->srcMaxInd[p].i;
                    jmin = abl->srcMinInd[p].j;
                    jmax = abl->srcMaxInd[p].j;
                    kmin = abl->srcMinInd[p].k;
                    kmax = abl->srcMaxInd[p].k;

                    numI = abl->srcNumI[p];
                    numJ = abl->srcNumJ[p];
                    numK = abl->srcNumK[p];

                    //loop through the sub domain of this processor
                    for (k=kmin; k<=kmax; k++)
                    {
                        for (j=jmin; j<=jmax; j++)
                        {
                            for (i=imin; i<=imax; i++)
                            {
                                tMapped[p][(k-kmin) * numJ * numI + (j-jmin) * numI + (i-imin)] = tP[k][j][i];
                            }
                        }
                    }
                }

                if(abl->isdestProc[p] == 1)
                {

                    numI = abl->srcNumI[p];
                    numJ = abl->srcNumJ[p];
                    numK = abl->srcNumK[p];

                    MPI_Ibcast(tMapped[p], numI * numJ * numK, MPIU_REAL, abl->srcCommLocalRank[p], abl->yDamp_comm[p], &(abl->mapRequest[p]));
                }

            }

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
            }
        }
    }
    return (0);
}
//***************************************************************************************************************//

PetscErrorCode dampingSourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
{
    abl_          *abl  = teqn->access->abl;
    mesh_         *mesh = teqn->access->mesh;
    ueqn_         *ueqn = teqn->access->ueqn;
    les_          *les  = teqn->access->les;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***cent;
    PetscReal     ***rhs, ***t, ***tP, ***tBarY;

    precursor_    *precursor;
    domain_       *pdomain;
    PetscInt      kStart, kEnd;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(da,  teqn->lTmprt, &t);
    DMDAVecGetArray(da,  Rhs,  &rhs);

    if(teqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
                kStart = precursor->map.kStart;
                kEnd   = precursor->map.kEnd;
            }
        }
    }

    if(teqn->access->flags->isYDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            DMDAVecGetArray(da, abl->tBarInstY, &tBarY);
        }
    }

    // x damping layer
    PetscReal alphaX = abl->xDampingAlpha;
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;
    PetscReal xD     = abl->xDampingDelta;

    // y damping layer
    PetscReal alphaY = abl->yDampingAlpha;
    PetscReal yS     = abl->yDampingStart;
    PetscReal yE     = abl->yDampingEnd;
    PetscReal yD     = abl->yDampingDelta;


    // loop over internal cell faces
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(teqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k
                    PetscReal x     = cent[k][j][i].x;

                    // compute Nordstrom viscosity at i,j,k
                    PetscReal nud_x   = viscNordstrom(alphaX, xS, xE, xD, x);
                    PetscReal tBarInstX;

                    if
                    (
                        abl->xFringeUBarSelectionType == 0 ||
                        abl->xFringeUBarSelectionType == 1 ||
                        abl->xFringeUBarSelectionType == 2 ||
                        abl->xFringeUBarSelectionType == 4
                    )
                    {
                        // set desired temperature
                        tBarInstX  = abl->tBarInstX[j][i];
                    }
                    else if(abl->xFringeUBarSelectionType == 3)
                    {
                        if
                        (
                            precursor->thisProcessorInFringe && // is this processor in the fringe?
                            k >= kStart && k <= kEnd            // is this face in the fringe?
                        )
                        {
                            // set desired temperature
                            tBarInstX  = tP[k-kStart][j][i];
                        }
                        else
                        {
                            tBarInstX  = t[k][j][i];
                        }
                    }

                    rhs[k][j][i] += scale * nud_x * (tBarInstX - t[k][j][i]);
                }

                if(teqn->access->flags->isYDampingActive)
                {
                    PetscReal y = cent[k][j][i].y;
                    PetscReal x = cent[k][j][i].x;

                    PetscReal nud_y  = viscNordstrom(alphaY, yS, yE, yD, y);
                    PetscReal fAsc_x = viscNordstromNoVertFilter(xS, xE, xD, x);

                    if
                    (
                        abl->xFringeUBarSelectionType == 0 ||
                        abl->xFringeUBarSelectionType == 1 ||
                        abl->xFringeUBarSelectionType == 2 ||
                        abl->xFringeUBarSelectionType == 4
                    )
                    {

                    }
                    else if(abl->xFringeUBarSelectionType == 3)
                    {
                        rhs[k][j][i] += scale * nud_y * fAsc_x * (tBarY[k][j][i] - t[k][j][i]);
                    }

                }

            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  teqn->lTmprt, &t);
    DMDAVecRestoreArray(da,  Rhs,  &rhs);

    if(teqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
            }
        }
    }

    if(teqn->access->flags->isYDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            DMDAVecRestoreArray(da, abl->tBarInstY, &tBarY);
        }
    }

    return(0);
}

//***************************************************************************************************************//


PetscErrorCode sourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = teqn->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***source, ***rhs;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  teqn->sourceT, &source);
    DMDAVecGetArray(da,  Rhs,  &rhs);

    // loop over internal cell faces
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhs[k][j][i] += scale * source[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(da,  teqn->sourceT, &source);
    DMDAVecRestoreArray(da,  Rhs,  &rhs);

    return(0);
}

//***************************************************************************************************************//


PetscErrorCode FormT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
{
    // In this function the viscous + divergence term of the temperature equation are
    // discretized at cell centers.
    // First the divergence and viscous fluxes are evaluated at cell faces, then
    // their budget is evaluated at the internal cells, forming the Rhs.

    mesh_         *mesh  = teqn->access->mesh;
    ueqn_         *ueqn  = teqn->access->ueqn;
    les_          *les   = teqn->access->les;
    constants_    *cst   = teqn->access->constants;
    flags_        *flags = teqn->access->flags;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***ucont;

    Cmpnts        ***csi, ***eta, ***zet, ***cent;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;

    PetscReal     ***tmprt, ***rhs, ***nvert;

    Cmpnts        ***div, ***visc, ***viscIBM;                                                // divergence and viscous terms
    Cmpnts        ***limiter;                                                     // flux limiter
    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;                                  // cell and face jacobians
    PetscReal     ***lnu_t, ***Sabs, ***lch, ***lcs;

    PetscReal     dtdc, dtde, dtdz;                                              // velocity der. w.r.t. curvil. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components

    PetscReal     tRef;
    if(ueqn->access->flags->isAblActive) tRef = teqn->access->abl->tRef;
    else                                 tRef = teqn->access->constants->tRef;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  Rhs,               &rhs);
    DMDAVecGetArray(da,  teqn->lTmprt,      &tmprt);
    DMDAVecGetArray(fda, ueqn->lUcont,      &ucont);

    DMDAVecGetArray(fda, mesh->lCsi,        &csi);
    DMDAVecGetArray(fda, mesh->lEta,        &eta);
    DMDAVecGetArray(fda, mesh->lZet,        &zet);
    DMDAVecGetArray(fda, mesh->lICsi,       &icsi);
    DMDAVecGetArray(fda, mesh->lIEta,       &ieta);
    DMDAVecGetArray(fda, mesh->lIZet,       &izet);
    DMDAVecGetArray(fda, mesh->lJCsi,       &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta,       &jeta);
    DMDAVecGetArray(fda, mesh->lJZet,       &jzet);
    DMDAVecGetArray(fda, mesh->lKCsi,       &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta,       &keta);
    DMDAVecGetArray(fda, mesh->lKZet,       &kzet);
    DMDAVecGetArray(da,  mesh->lAj,         &aj);
    DMDAVecGetArray(da,  mesh->lIAj,        &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,        &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,        &kaj);
    DMDAVecGetArray(fda, mesh->lCent,       &cent);
    DMDAVecGetArray(da,  mesh->lNvert,      &nvert);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    if(flags->isIBMActive)
    {
        DMDAVecGetArray(fda, teqn->lViscIBMT, &viscIBM);
    }

    VecSet(teqn->lDivT,  0.0);
    VecSet(teqn->lViscT, 0.0);;

    DMDAVecGetArray(fda, teqn->lDivT, &div);
    DMDAVecGetArray(fda, teqn->lViscT, &visc);

    if (teqn->access->flags->isLesActive)
    {
        if(teqn->access->flags->isLesActive != 2)
        {
            DMDAVecGetArray(da, les->lNu_t, &lnu_t);
        }
        else
        {
            DMDAVecGetArray(da,  les->lS,  &Sabs);
            DMDAVecGetArray(da,  les->lCh, &lch);
            DMDAVecGetArray(da,  les->lCs, &lcs);
        }
    }

    // ---------------------------------------------------------------------- //
    //                FORM DIVERGENCE AND VISCOUS CONTRIBUTIONS               //
    // ---------------------------------------------------------------------- //

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

                // compute metric tensor - WARNING: there is a factor of 1/J^2 if using face area vectors
                //                                  must multiply for ajc in viscous term!!!
                g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                Compute_dscalar_i(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    iL, iR;
                PetscReal denom;
                getFace2Cell4StencilCsi(mesh, k, j, i, mx, &iL, &iR, &denom, nvert);

                div[k][j][i].x =
                - ucont[k][j][i].x
                * centralUpwind
                (
                    tmprt[k][j][iL],
                    tmprt[k][j][i],
                    tmprt[k][j][i+1],
                    tmprt[k][j][iR],
                    ucont[k][j][i].x
                    ,limiter[k][j][i].x
                );

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                // viscous terms
                if
                (
                    teqn->access->flags->isLesActive
                )
                {

                    if(teqn->access->flags->isLesActive != 2)
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                        // compute stability dependent turbulent Prandtl number
                        PetscReal gradTdotG = dtde*(-9.81);
                        PetscReal l, delta = pow( 1./ajc, 1./3. );
                        if(gradTdotG < 0.)
                        {
                            l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                        }
                        else
                        {
                            l = delta;
                        }

                        PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));


                        kappaEff = (nu / cst->Pr) + (nut / Prt);
                    }
                    else
                    {
                        kappaEff = (nu / cst->Pr) + 0.5 *
                        (
                            lch[k][j][i+1] / 0.1 * Sabs[k][j][i+1] * lcs[k][j][i+1] +
                            lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                        ) * pow(1.0 / ajc, 1.0/3.0);
                    }

                    // wall model i-left/right patch
                    if
                    (
                        (mesh->boundaryT.iLeft=="thetaWallFunction" && i==0) ||
                        (mesh->boundaryT.iRight=="thetaWallFunction" && i==mx-2)
                    )
                    {
                        PetscReal signQ =  1.0;
                        if(i==0)  signQ = -1.0;

                        visc[k][j][i].x
                        =
                        signQ *
                        (
                            teqn->iRWM->qWall.x[k-zs][j-ys] * icsi[k][j][i].x +
                            teqn->iRWM->qWall.y[k-zs][j-ys] * icsi[k][j][i].y +
                            teqn->iRWM->qWall.z[k-zs][j-ys] * icsi[k][j][i].z
                        );

                        kappaEff = 0.0;
                    }

                    //IBM wall model
                    if(isIBMFluidIFace(k, j, i, i+1, nvert))
                    {
                        if(teqn->access->ibm->wallShearOn && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                        {
                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc[k][j][i].x = viscIBM[k][j][i].x;
                            }
                            else if(isIBMFluidCell(k, j, i+1, nvert))
                            {
                                visc[k][j][i].x = viscIBM[k][j][i+1].x;
                            }

                            kappaEff = 0;
                        }
                    }


                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].x += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
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

                // compute metric tensor
                g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
                g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
                g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

                Compute_dscalar_j(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    jL, jR;
                PetscReal denom;
                getFace2Cell4StencilEta(mesh, k, j, i, my, &jL, &jR, &denom, nvert);

                div[k][j][i].y =
                - ucont[k][j][i].y
                * centralUpwind
                (
                    tmprt[k][jL][i],
                    tmprt[k][j][i],
                    tmprt[k][j+1][i],
                    tmprt[k][jR][i],
                    ucont[k][j][i].y
                    ,limiter[k][j][i].y
                );

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                if
                (
                    teqn->access->flags->isLesActive
                )
                {

                    if(teqn->access->flags->isLesActive != 2)
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                        // compute stability depentend turbulent Prandtl number
                        PetscReal gradTdotG = dtde*(-9.81);
                        PetscReal l, delta = pow( 1./ajc, 1./3. );
                        if(gradTdotG < 0.)
                        {
                            l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                        }
                        else
                        {
                            l = delta;
                        }

                        PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));

                        kappaEff = (nu / cst->Pr) + (nut / Prt);
                    }
                    else
                    {
                        kappaEff = (nu / cst->Pr) + 0.5 *
                        (
                            lch[k][j+1][i] / 0.1 * Sabs[k][j+1][i] * lcs[k][j+1][i] +
                            lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                        ) * pow(1.0 / ajc, 1.0/3.0);
                    }

                    // wall model j-left patch
                    if
                    (
                        (mesh->boundaryT.jLeft=="thetaWallFunction"  && j==0) ||
                        (mesh->boundaryT.jRight=="thetaWallFunction" && j==my-2)
                    )
                    {
                        PetscReal signQ =  1.0;
                        if(j==0)  signQ = -1.0;

                        visc[k][j][i].y
                        =
                        signQ *
                        (
                            teqn->jLWM->qWall.x[k-zs][i-xs] * jeta[k][j][i].x +
                            teqn->jLWM->qWall.y[k-zs][i-xs] * jeta[k][j][i].y +
                            teqn->jLWM->qWall.z[k-zs][i-xs] * jeta[k][j][i].z
                        );

                        kappaEff = 0.0;
                    }

                    if(isIBMFluidJFace(k, j, i, j+1, nvert) && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                    {
                        if(teqn->access->ibm->wallShearOn)
                        {

                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc[k][j][i].y = viscIBM[k][j][i].y;
                            }
                            else if(isIBMFluidCell(k, j+1, i, nvert))
                            {
                                visc[k][j][i].y = viscIBM[k][j+1][i].y;
                            }

                            kappaEff = 0;
                        }

                    }
                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].y += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
            }
        }
    }

    // k direction
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

                // compute metric tensor
                g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                Compute_dscalar_k(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    kL, kR;
                PetscReal denom;
                getFace2Cell4StencilZet(mesh, k, j, i, mz, &kL, &kR, &denom, nvert);

                div[k][j][i].z =
                - ucont[k][j][i].z
                * centralUpwind
                (
                    tmprt[kL][j][i],
                    tmprt[k][j][i],
                    tmprt[k+1][j][i],
                    tmprt[kR][j][i],
                    ucont[k][j][i].z
                    ,limiter[k][j][i].z
                );

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                if
                (
                    teqn->access->flags->isLesActive
                )
                {

                    if(teqn->access->flags->isLesActive != 2)
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

                        // compute stability depentend turbulent Prandtl number
                        PetscReal gradTdotG = dtde*(-9.81);
                        PetscReal l, delta = pow( 1./ajc, 1./3. );
                        if(gradTdotG < 0.)
                        {
                            l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                        }
                        else
                        {
                            l = delta;
                        }

                        PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));

                        kappaEff = (nu / cst->Pr) + (nut / Prt);
                    }
                    else
                    {
                        kappaEff = (nu / cst->Pr) + 0.5 *
                        (
                            lch[k+1][j][i] / 0.1 * Sabs[k+1][j][i] * lcs[k+1][j][i] +
                            lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                        ) * pow(1.0 / ajc, 1.0/3.0);
                    }

                    if(isIBMFluidKFace(k, j, i, k+1, nvert) && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                    {
                        if(teqn->access->ibm->wallShearOn)
                        {
                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc[k][j][i].z = viscIBM[k][j][i].z;
                            }
                            else if(isIBMFluidCell(k+1, j, i, nvert))
                            {
                                visc[k][j][i].z = viscIBM[k+1][j][i].z;
                            }

                            kappaEff = 0;
                        }
                    }
                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].z += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
            }
        }
    }

    DMDAVecRestoreArray(fda, teqn->lDivT, &div);
    DMDAVecRestoreArray(fda, teqn->lViscT, &visc);

    if(flags->isIBMActive)
    {
        DMDAVecRestoreArray(fda, teqn->lViscIBMT, &viscIBM);
    }

    DMLocalToLocalBegin(fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalEnd  (fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalBegin(fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);
    DMLocalToLocalEnd  (fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);

    resetFacePeriodicFluxesVector(mesh, teqn->lDivT,   teqn->lDivT, "localToLocal");
    resetFacePeriodicFluxesVector(mesh, teqn->lViscT, teqn->lViscT, "localToLocal");

    DMDAVecGetArray(fda, teqn->lDivT, &div);
    DMDAVecGetArray(fda, teqn->lViscT, &visc);

    // ---------------------------------------------------------------------- //
    //                      FORM THE CUMULATIVE FLUXES                        //
    // ---------------------------------------------------------------------- //

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhs[k][j][i]
                =
                scale *
                (
                        div[k][j][i].x  - div[k  ][j  ][i-1].x       +
                        div[k][j][i].y  - div[k  ][j-1][i  ].y       +
                        div[k][j][i].z  - div[k-1][j  ][i  ].z       +
                        visc[k][j][i].x - visc[k  ][j  ][i-1].x      +
                        visc[k][j][i].y - visc[k  ][j-1][i  ].y      +
                        visc[k][j][i].z - visc[k-1][j  ][i  ].z
                ) * aj[k][j][i];
            }
        }
    }

    if(teqn->access->flags->isLesActive)
    {
        if(teqn->access->flags->isLesActive != 2)
        {
            DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
        }
        else
        {
            DMDAVecRestoreArray(da,  les->lS,  &Sabs);
            DMDAVecRestoreArray(da,  les->lCh, &lch);
            DMDAVecRestoreArray(da,  les->lCs, &lcs);
        }
    }

    DMDAVecRestoreArray(fda, teqn->lDivT,       &div);
    DMDAVecRestoreArray(fda, teqn->lViscT,      &visc);

    DMDAVecRestoreArray(da,  Rhs,              &rhs);
    DMDAVecRestoreArray(da,  teqn->lTmprt,      &tmprt);
    DMDAVecRestoreArray(fda, ueqn->lUcont,     &ucont);

    DMDAVecRestoreArray(fda, mesh->lCsi,        &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,        &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,        &zet);
    DMDAVecRestoreArray(fda, mesh->lICsi,       &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta,       &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet,       &izet);
    DMDAVecRestoreArray(fda, mesh->lJCsi,       &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,       &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet,       &jzet);
    DMDAVecRestoreArray(fda, mesh->lKCsi,       &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta,       &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet,       &kzet);
    DMDAVecRestoreArray(da,  mesh->lAj,         &aj);
    DMDAVecRestoreArray(da,  mesh->lIAj,        &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,        &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,        &kaj);
    DMDAVecRestoreArray(fda, mesh->lCent,       &cent);
    DMDAVecRestoreArray(da,  mesh->lNvert,      &nvert);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TeqnSNES(SNES snes, Vec T, Vec Rhs, void *ptr)
{
    teqn_ *teqn   = (teqn_*)ptr;
    mesh_ *mesh   = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;
    VecCopy(T, teqn->Tmprt);

    // scatter temperature from global to local
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // reset temperature periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // update wall model (optional)
    UpdateWallModelsT(teqn);

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    const PetscReal dt = clock->dt;

    // add viscous and transport terms
    FormT(teqn, Rhs, 1.0);

    if(teqn->access->flags->isXDampingActive)
    {
        // this is causing spurious oscillation: deactivate
        dampingSourceT(teqn, Rhs, 1.0);
    }

    // multiply for dt
    VecScale(Rhs, dt);

    // add driving source terms after as it is not scaled by 1/dt
    if(teqn->access->flags->isAblActive)
    {
        if(teqn->access->abl->controllerActiveT)
        {
            sourceT(teqn, Rhs, 1.0);
        }
    }

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, Rhs);

    VecAXPY(Rhs, -1.0, T);
    VecAXPY(Rhs,  1.0, teqn->Tmprt_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsT(teqn_ *teqn)
{
    mesh_ *mesh   = teqn->access->mesh;

    // reset temperature periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // update wall model (optional)
    UpdateWallModelsT(teqn);

    // initialize the rhs vector
    VecSet(teqn->Rhs, 0.0);

    // add viscous and transport terms
    FormT(teqn, teqn->Rhs, 1.0);

    if(teqn->access->flags->isXDampingActive)
    {
        // this is causing spurious oscillation: deactivate
        dampingSourceT(teqn, teqn->Rhs, 1.0);
    }

    // add driving source terms after as it is not scaled by 1/dt
    if(teqn->access->flags->isAblActive)
    {
        if(teqn->access->abl->controllerActiveT)
        {
            sourceT(teqn, teqn->Rhs, 1.0 / teqn->access->clock->dt);
        }
    }

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, teqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TeqnRK4(teqn_ *teqn)
{
    mesh_  *mesh  = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for T, Stage ");

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

    // Tmprt_o contribution
    VecCopy(teqn->Tmprt_o, teqn->TmprtTmp);

    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate U guess and evaluate RHS
        if(i!=0)
        {
            VecWAXPY(teqn->Tmprt, a[i] * dt, teqn->Rhs, teqn->Tmprt_o);
            DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
            DMGlobalToLocalEnd(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        }

        // compute function guess
        FormExplicitRhsT(teqn);

        // add contribution from K1, K2, K3, K4
        VecAXPY(teqn->TmprtTmp, dt * b[i], teqn->Rhs);
    }

    VecCopy(teqn->TmprtTmp, teqn->Tmprt);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveTEqn(teqn_ *teqn)
{
    mesh_          *mesh = teqn->access->mesh;

    // set the right hand side to zero
    VecSet (teqn->Rhs, 0.0);

    if(teqn->ddtScheme=="backwardEuler")
    {
        PetscReal     norm;
        PetscInt      iter;
        PetscReal     ts,te;

        PetscTime(&ts);
        PetscPrintf(mesh->MESH_COMM, "TRSNES: Solving for T, Initial residual = ");

        // copy the solution in the tmp variable
        VecCopy(teqn->Tmprt, teqn->TmprtTmp);

        // solve temperature equation and compute solution norm and iteration number
        SNESSolve(teqn->snesT, PETSC_NULL, teqn->TmprtTmp);
        SNESGetFunctionNorm(teqn->snesT, &norm);
        SNESGetIterationNumber(teqn->snesT, &iter);

        // store the solution and global to local scatter
        VecCopy(teqn->TmprtTmp, teqn->Tmprt);

        // scatter
        DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

        // compute elapsed time
        PetscTime(&te);
        PetscPrintf(mesh->MESH_COMM,"Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n", norm, iter, te-ts);
    }
    else if (teqn->ddtScheme=="rungeKutta4")
    {
        TeqnRK4(teqn);
    }

    return(0);
}
