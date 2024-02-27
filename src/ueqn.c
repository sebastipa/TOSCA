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
    mesh_  *mesh  = ueqn->access->mesh;
    flags_ *flags = ueqn->access->flags;

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
	VecDuplicate(mesh->lCent, &(ueqn->sourceU));  VecSet(ueqn->sourceU, 0.0);

    VecDuplicate(mesh->lAj, &(ueqn->lUstar));     VecSet(ueqn->lUstar,  0.0);

    if(ueqn->access->flags->isTeqnActive)
    {
        VecDuplicate(mesh->Cent, &(ueqn->bTheta)); VecSet(ueqn->bTheta,    0.0);
    }

    if(flags->isIBMActive)
    {
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM1));   VecSet(ueqn->lViscIBM1,  0.0);
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM2));   VecSet(ueqn->lViscIBM2,  0.0);
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM3));   VecSet(ueqn->lViscIBM3,  0.0);
    }

	if(ueqn->access->flags->isGravityWaveModelingActive)
    {
		VecDuplicate(mesh->Cent, &(ueqn->dPAGW));        VecSet(ueqn->dPAGW,      0.0);
	}

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
        SNESSetTolerances(ueqn->snesU, ueqn->absExitTol, ueqn->relExitTol, 1e-30, 20, 1000);

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
        PetscReal rtol=ueqn->relExitTol, atol=ueqn->absExitTol, dtol=1e30;
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
    Cmpnts        *src;

    PetscInt      applyGeoDamping = 0;

    abl_          *abl  = ueqn->access->abl;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->sourceU, &source);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    // compute source terms using the velocity controller
    if(abl->controllerType=="pressure")
    {
        PetscReal     relax = abl->relax;
        PetscReal     alpha = abl->alpha;
        PetscReal     T     = abl->timeWindow;

        // set the wanted velocity
        uDes.x = abl->uRef;
        uDes.y = 0.0;
        uDes.z = 0.0;

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
        abl->cumulatedSource.x = (1.0 - clock->dt / T) * abl->cumulatedSource.x + (clock->dt / T) * error.x;
        abl->cumulatedSource.y = (1.0 - clock->dt / T) * abl->cumulatedSource.y + (clock->dt / T) * error.y;
        abl->cumulatedSource.z = 0.0;

        // compute the uniform source terms (PI controller with adjustable gains)
        s.x = relax * (alpha * error.x + (1.0-alpha) * abl->cumulatedSource.x) ;
        s.y = relax * (alpha * error.y + (1.0-alpha) * abl->cumulatedSource.y) ;
        s.z = 0.0;

		if(abl->geostrophicDampingActive)
        {
			// spatially average velocity
			std::vector<Cmpnts>    lgDes(nLevels);

			// set to zero
			for(j=0; j<nLevels; j++) lgDes[j] = nSetZero();

			DMDAVecGetArray(da, mesh->lAj, &aj);

			for(j=lys; j<lye; j++)
			{
				for(k=lzs; k<lze; k++)
				{
					for(i=lxs; i<lxe; i++)
					{
						mSum
						(
							lgDes[j-1],
							nScale(1.0/aj[k][j][i], ucat[k][j][i])
						);
					}
				}
			}

			DMDAVecRestoreArray(da, mesh->lAj, &aj);

			MPI_Allreduce(&lgDes[0], &(abl->geoDampU[0]), 3*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

			// divide by total volume per level
			for(j=0; j<nLevels; j++)
			{
				mScale(1.0/abl->totVolPerLevel[j], abl->geoDampU[j]);
			}

            // filter geostrophic wind
            Cmpnts    Sc      = nSetFromComponents(s.x, s.y, 0.0);
            PetscReal m1      = (1.0 - clock->dt/abl->geoDampWindow),
                      m2      = clock->dt/abl->geoDampWindow;
	        abl->geoDampAvgDT = m1 * abl->geoDampAvgDT + m2 * clock->dt;
			abl->geoDampAvgS  = nSum(nScale(m1, abl->geoDampAvgS), nScale(m2, Sc));
            abl->geoDampUBar  = nSetFromComponents(abl->geoDampAvgS.y / (2.0 * abl->fc * abl->geoDampAvgDT), -1.0 * abl->geoDampAvgS.x / (2.0 * abl->fc * abl->geoDampAvgDT), 0.0);

            if(print) PetscPrintf(mesh->MESH_COMM, "                         Filtered geostrophic wind Ug = (%.3lf, %.3lf, 0.00)\n", abl->geoDampUBar.x, abl->geoDampUBar.y);

            // set damping flag
            if(clock->time > abl->geoDampStart)
            {
                applyGeoDamping = 1;
            }

            if(mesh->access->io->runTimeWrite)
            {
                if(!rank)
                {
                    word location = "./fields/" + mesh->meshName + "/" + getTimeName(clock);
                    word fileName = location + "/geostrophicDampingInfo";

                    FILE *fp=fopen(fileName.c_str(), "w");

                    if(fp==NULL)
                    {
                        char error[512];
                        sprintf(error, "cannot open file %s\n", fileName.c_str());
                        fatalErrorInFunction("ABLInitialize",  error);
                    }
                    else
                    {
                        fprintf(fp, "filteredS \t\t(%.6e %.6e %.6e)\n", abl->geoDampAvgS.x, abl->geoDampAvgS.y, 0.0);
						fprintf(fp, "filteredDT\t\t %.5lf\n", abl->geoDampAvgDT);
                        fclose(fp);
                    }
                }
            }
		}

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
            }

            word fileName = "postProcessing/momentumSource_" + getStartTimeName(clock);
            FILE *fp = fopen(fileName.c_str(), "a");

            if(!fp)
            {
                char error[512];
                sprintf(error, "cannot open file %s\n", fileName.c_str());
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

                PetscFPrintf(mesh->MESH_COMM, fp, "%*.6f\t%*.5e\t%*.5e%*.5e\t\n", width, clock->time, width, s.x, width, s.y, width, s.z);

                fclose(fp);
            }
        }
    }
    // compute source terms using the velocity controller
    else if(abl->controllerType=="geostrophic")
    {
        PetscReal     relax = abl->relax;
        PetscReal     alpha = abl->alpha;
        PetscReal     T     = abl->timeWindow;

        DMDAVecGetArray(da, mesh->lAj, &aj);

        // find the first two closest levels
        PetscReal nLevels = my-2;

        Cmpnts luMean1; luMean1.x = 0.0; luMean1.y = 0.0; luMean1.z = 0.0;
        Cmpnts luMean2; luMean2.x = 0.0; luMean2.y = 0.0; luMean2.z = 0.0;
        Cmpnts guMean1; guMean1.x = 0.0; guMean1.y = 0.0; guMean1.z = 0.0;
        Cmpnts guMean2; guMean2.x = 0.0; guMean2.y = 0.0; guMean2.z = 0.0;
        Cmpnts luMeanGeo1; luMeanGeo1.x = 0.0; luMeanGeo1.y = 0.0; luMeanGeo1.z = 0.0;
        Cmpnts luMeanGeo2; luMeanGeo2.x = 0.0; luMeanGeo2.y = 0.0; luMeanGeo2.z = 0.0;
        Cmpnts guMeanGeo1; guMeanGeo1.x = 0.0; guMeanGeo1.y = 0.0; guMeanGeo1.z = 0.0;
        Cmpnts guMeanGeo2; guMeanGeo2.x = 0.0; guMeanGeo2.y = 0.0; guMeanGeo2.z = 0.0;

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

                    if(j==abl->closestLabelsGeo[0])
                    {
                        luMeanGeo1.x += ucat[k][j][i].x / aj[k][j][i];
                        luMeanGeo1.y += ucat[k][j][i].y / aj[k][j][i];
                        luMeanGeo1.z += ucat[k][j][i].z / aj[k][j][i];
                    }
                    else if(j==abl->closestLabelsGeo[1])
                    {
                        luMeanGeo2.x += ucat[k][j][i].x / aj[k][j][i];
                        luMeanGeo2.y += ucat[k][j][i].y / aj[k][j][i];
                        luMeanGeo2.z += ucat[k][j][i].z / aj[k][j][i];
                    }
                }
            }
        }

        DMDAVecRestoreArray(da, mesh->lAj, &aj);

        MPI_Allreduce(&luMean1, &guMean1, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&luMean2, &guMean2, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&luMeanGeo1, &guMeanGeo1, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&luMeanGeo2, &guMeanGeo2, 3, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        guMean1.x = guMean1.x / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean1.y = guMean1.y / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean1.z = guMean1.z / abl->totVolPerLevel[abl->closestLabels[0]-1];
        guMean2.x = guMean2.x / abl->totVolPerLevel[abl->closestLabels[1]-1];
        guMean2.y = guMean2.y / abl->totVolPerLevel[abl->closestLabels[1]-1];
        guMean2.z = guMean2.z / abl->totVolPerLevel[abl->closestLabels[1]-1];

        guMeanGeo1.x = guMeanGeo1.x / abl->totVolPerLevel[abl->closestLabelsGeo[0]-1];
        guMeanGeo1.y = guMeanGeo1.y / abl->totVolPerLevel[abl->closestLabelsGeo[0]-1];
        guMeanGeo1.z = guMeanGeo1.z / abl->totVolPerLevel[abl->closestLabelsGeo[0]-1];
        guMeanGeo2.x = guMeanGeo2.x / abl->totVolPerLevel[abl->closestLabelsGeo[1]-1];
        guMeanGeo2.y = guMeanGeo2.y / abl->totVolPerLevel[abl->closestLabelsGeo[1]-1];
        guMeanGeo2.z = guMeanGeo2.z / abl->totVolPerLevel[abl->closestLabelsGeo[1]-1];

        Cmpnts uMean, uMeanGeo;

        uMean.x = guMean1.x * abl->levelWeights[0] + guMean2.x * abl->levelWeights[1];
        uMean.y = guMean1.y * abl->levelWeights[0] + guMean2.y * abl->levelWeights[1];
        uMean.z = guMean1.z * abl->levelWeights[0] + guMean2.z * abl->levelWeights[1];

        uMeanGeo.x = guMeanGeo1.x * abl->levelWeightsGeo[0] + guMeanGeo2.x * abl->levelWeightsGeo[1];
        uMeanGeo.y = guMeanGeo1.y * abl->levelWeightsGeo[0] + guMeanGeo2.y * abl->levelWeightsGeo[1];
        uMeanGeo.z = guMeanGeo1.z * abl->levelWeightsGeo[0] + guMeanGeo2.z * abl->levelWeightsGeo[1];

        if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: wind height is %lf m, h1 = %lf m, h2 = %lf m\n", abl->hRef, abl->cellLevels[abl->closestLabels[0]-1], abl->cellLevels[abl->closestLabels[1]-1]);
		if(print) PetscPrintf(mesh->MESH_COMM, "                         U at hRef: (%.3f %.3f %.3f) m/s, U at hGeo: (%.3f %.3f %.3f) m/s\n", uMean.x, uMean.y, uMean.z, uMeanGeo.x, uMeanGeo.y, uMeanGeo.z);

        // compute actual angles
        PetscReal hubAngle    = std::atan(uMean.y/uMean.x);
        PetscReal geoAngleNew = std::atan(uMeanGeo.y/uMeanGeo.x);

        // compute filtered hub angle
        abl->hubAngle = (1.0 - clock->dt / T) * abl->hubAngle + (clock->dt / T) * hubAngle;

        // compute filtered geostrophic angle
        abl->geoAngle = (1.0 - clock->dt / T) * abl->geoAngle + (clock->dt / T) * geoAngleNew;

		if(print) PetscPrintf(mesh->MESH_COMM, "                         Alpha at hRef: %.3f deg, Alpha at hGeo: %.3f deg\n",abl->hubAngle/M_PI*180, abl->geoAngle/M_PI*180);

        // compute the uniform source terms
        abl->a.x = -2.0*abl->fc*nMag(abl->uGeoBar)*sin(abl->geoAngle);
        abl->a.y =  2.0*abl->fc*nMag(abl->uGeoBar)*cos(abl->geoAngle);
        abl->a.z =  0.0;

        if(!relax)
        {
            abl->b.x = 0.0;
            abl->b.y = 0.0;
            abl->b.z = 0.0;
        }
        else
        {
            // compute previous geostrophic angle and delta angle w.r.t actual
            PetscReal geoAngleOld = std::atan(abl->uGeoBar.y/abl->uGeoBar.x);
            PetscReal geoDelta    = geoAngleNew - geoAngleOld;

            // compute rotation angle at hub height
            PetscReal hubDelta = hubAngle;
            PetscReal omega    = hubDelta / clock->dt;

            // time constant
            PetscReal sigma    = clock->dt / 200; // filter hardcoded to 200s

            // compute omega bar
            abl->omegaBar    = sigma * omega + (1.0 - sigma) * abl->omegaBar;

            // rotate geostrophic speed
            Cmpnts uGeoBarTmp = nSetZero();
            uGeoBarTmp.x = std::cos(geoDelta) * abl->uGeoBar.x - std::sin(geoDelta) * abl->uGeoBar.y;
            uGeoBarTmp.y = std::sin(geoDelta) * abl->uGeoBar.x + std::cos(geoDelta) * abl->uGeoBar.y;
            mSet(abl->uGeoBar, uGeoBarTmp);

            abl->b.x =   relax*(abl->omegaBar + 1.0/T*(abl->hubAngle));
            abl->b.y = - relax*(abl->omegaBar + 1.0/T*(abl->hubAngle));
            abl->b.z =   0.0;
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
    else if(abl->controllerType=="directProfileAssimilation")
    {
        PetscReal relax   = abl->relax;
        PetscInt  nlevels = my-2;

        Cmpnts    *luMean = abl->luMean;
        Cmpnts    *guMean = abl->guMean;
        Cmpnts    **uMeso = abl->uMeso;
        Cmpnts    uDes, uH1, uH2;
        PetscInt  idxh1, idxh2, idxt1;
        PetscReal wth1, wth2, wtt1;

        src    = abl->source;

        for(j=0; j<nlevels; j++)
        {
            luMean[j] = nSetZero();
            guMean[j] = nSetZero();
            src[j]    = nSetZero();
        }

        DMDAVecGetArray(da, mesh->lAj, &aj);

        //loop through the cells and find the mean at each vertical cell level
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    luMean[j-1].x += ucat[k][j][i].x / aj[k][j][i];
                    luMean[j-1].y += ucat[k][j][i].y / aj[k][j][i];
                    luMean[j-1].z += ucat[k][j][i].z / aj[k][j][i];
                }
            }
        }

        DMDAVecRestoreArray(da, mesh->lAj, &aj);

        MPI_Allreduce(&(luMean[0]), &(guMean[0]), 3*nlevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        if(abl->assimilationType == "constantTime")
        {
            for(j=0; j<nlevels; j++)
            {

                mScale(1.0/abl->totVolPerLevel[j], guMean[j]);
                
                idxh1 = abl->velInterpIdx[j][0];
                idxh2 = abl->velInterpIdx[j][1];

                wth1  = abl->velInterpWts[j][0];
                wth2  = abl->velInterpWts[j][1];

                idxt1 = abl->closestTimeIndV;
                wtt1 = abl->closestTimeWtV;

                //interpolating in time
                uH1 = nSum(nScale(wtt1, abl->uMeso[idxh1][idxt1]), nScale(1 - wtt1, abl->uMeso[idxh1][idxt1 + 1]));
                uH2 = nSum(nScale(wtt1, abl->uMeso[idxh2][idxt1]), nScale(1 - wtt1, abl->uMeso[idxh2][idxt1 + 1]));

                //interpolating in vertical direction 
                uDes = nSum(nScale(wth1, uH1), nScale(wth2, uH2));

                // PetscPrintf(PETSC_COMM_WORLD, "uDes = %lf %lf %lf, guMean = %lf %lf %lf\n", uDes.x, uDes.y, uDes.z, guMean[j].x, guMean[j].y, guMean[j].z);
                
                src[j].x = (uDes.x - guMean[j].x);
                src[j].y = (uDes.y - guMean[j].y);
                src[j].z = (uDes.z - guMean[j].z);

                mScale(relax, src[j]);
            }
        }
        else if(abl->assimilationType == "variableTime")
        {
            //find the two closest available mesoscale data in time
            PetscInt idx_1 = abl->closestTimeIndV;
            PetscInt idx_2 = abl->closestTimeIndV + 1;

            PetscInt lwrBound = 0;
            PetscInt uprBound = abl->numtV;

            if(clock->it > clock->itStart)
            {
                lwrBound = PetscMax(0, (abl->closestTimeIndV - 50));
                uprBound = PetscMin(abl->numtV, (abl->closestTimeIndV + 50));
            }

            // build error vector for the time search
            PetscReal  diff[abl->numtV];

            for(PetscInt i=lwrBound; i<uprBound; i++)
            {
                diff[i] = fabs(abl->timeV[i] - clock->time);
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
            if(abl->timeV[idx_2] < abl->timeV[idx_1])
            {
                PetscInt idx_tmp = idx_2;
                idx_2 = idx_1;
                idx_1 = idx_tmp;
            }

            // find interpolation weights
            PetscReal idx = (idx_2 - idx_1) / (abl->timeV[idx_2] - abl->timeV[idx_1]) * (clock->time - abl->timeV[idx_1]) + idx_1;
            PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
            PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

            PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading mesoscale data\n", w1 * abl->timeV[idx_1] + w2 * abl->timeV[idx_2]);
            PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
            PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeV[idx_1], abl->timeV[idx_2]);

            // reset the closest index for nex iteration
            abl->closestTimeIndV = idx_1;

            for(j=0; j<nlevels; j++)
            {

                mScale(1.0/abl->totVolPerLevel[j], guMean[j]);
                
                idxh1 = abl->velInterpIdx[j][0];
                idxh2 = abl->velInterpIdx[j][1];

                wth1  = abl->velInterpWts[j][0];
                wth2  = abl->velInterpWts[j][1];

                //interpolating in time
                uH1 = nSum(nScale(w1, abl->uMeso[idxh1][idx_1]), nScale(w2, abl->uMeso[idxh1][idx_2]));
                uH2 = nSum(nScale(w1, abl->uMeso[idxh2][idx_1]), nScale(w2, abl->uMeso[idxh2][idx_2]));

                //interpolating in vertical direction 
                uDes = nSum(nScale(wth1, uH1), nScale(wth2, uH2));

                // PetscPrintf(PETSC_COMM_WORLD, "uDes = %lf %lf %lf, guMean = %lf %lf %lf, interpolated from %lf and %lf, wts = %lf %lf\n", uDes.x, uDes.y, uDes.z, guMean[j].x, guMean[j].y, guMean[j].z, abl->timeV[idx_1], abl->timeV[idx_2], w1, w2);
                
                src[j].x = (uDes.x - guMean[j].x);
                src[j].y = (uDes.y - guMean[j].y);
                src[j].z = (uDes.z - guMean[j].z);

                mScale(relax, src[j]);
            }

        }
        else 
        {
            char error[512];
            sprintf(error, "wrong assimilation method chosen. Available options are constantTime or variableTime\n");
            fatalErrorInFunction("CorrectSourceTerms",  error);
        }

        // if cellLevels is below the lowest mesoscale data point use the source at the last available height
        for(j=0; j<nlevels; j++)
        {
            if(abl->cellLevels[j] < abl->hV[0])
            {
                src[j].x = src[abl->lowestIndV].x;
                src[j].y = src[abl->lowestIndV].y;
                src[j].z = src[abl->lowestIndV].z;
            }
        }

        //write the source terms 
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
            }

            word fileName = "postProcessing/momentumSource_" + getStartTimeName(clock);
            FILE *fp = fopen(fileName.c_str(), "a");

            if(!fp)
            {
                char error[512];
                sprintf(error, "cannot open file postProcessing/momentumSource\n");
                fatalErrorInFunction("correctSourceTermT",  error);
            }
            else
            {

                PetscInt width = -15;

                if(clock->it == clock->itStart)
                {
                    word w1 = "levels";
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*s\n", width, w1.c_str());

                    for(j=0; j<nlevels; j++)
                    {
                        PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, abl->cellLevels[j]);
                    }

                    PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                    word w2 = "time";
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*s\n", width, w2.c_str());
                } 

                PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, clock->time);

                for(j=0; j<nlevels; j++)
                {
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e  %*.5e  %*.5e\t", width, src[j].x,  width, src[j].y,  width, src[j].z);
                }

                PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                fclose(fp);                  
            
            }
        }
    }

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // update the source term
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(cent[k][j][i].z <= abl->controllerMaxHeight)
                {
                    if(abl->controllerType=="geostrophic")
                    {
						// multiply by dt for how TOSCA builds the source RHS
                        source[k][j][i].x = (abl->a.x + abl->b.x*ucat[k][j][i].y) * clock->dt;
                        source[k][j][i].y = (abl->a.y + abl->b.y*ucat[k][j][i].x) * clock->dt;
                        source[k][j][i].z = 0.0;
                    }
                    else if(abl->controllerType=="pressure")
                    {
                        if(applyGeoDamping)
                        {
							// multiply by dt for how TOSCA builds the source RHS (only geo damping)

                            PetscReal  height = cent[k][j][i].z - mesh->bounds.zmin;
                            source[k][j][i].x = s.x +
                                                scaleHyperTangTop   (height, abl->geoDampH, abl->geoDampDelta) *
                                                abl->geoDampC*abl->geoDampAlpha*(abl->geoDampUBar.x - abl->geoDampU[j-1].x) * clock->dt;
                            source[k][j][i].y = s.y +
                                                scaleHyperTangTop   (height, abl->geoDampH, abl->geoDampDelta) *
                                                abl->geoDampC*abl->geoDampAlpha*(abl->geoDampUBar.y - abl->geoDampU[j-1].y) * clock->dt;
                            source[k][j][i].z = 0.0;
                        }
                        else
                        {
                            source[k][j][i].x = s.x;
                            source[k][j][i].y = s.y;
                            source[k][j][i].z = 0.0;
                        }
                    }
                    else if(abl->controllerType=="directProfileAssimilation")
                    {
                        source[k][j][i].x = src[j-1].x;
                        source[k][j][i].y = src[j-1].y;
                        source[k][j][i].z = 0.0;
                    } 
                    else
                    {
                        source[k][j][i].x = s.x;
                        source[k][j][i].y = s.y;
                        source[k][j][i].z = 0.0;
                    }
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

	// apply zero-gradient boundary conditions
	for (k=zs; k<ze; k++)
	{
		for (j=ys; j<ye; j++)
		{
			for (i=xs; i<xe; i++)
			{
				PetscInt flag=0, a=i, b=j, c=k;

				if(i==0)         a=1,    flag=1;
				else if(i==mx-1) a=mx-2, flag=1;

				if(j==0)         b=1,    flag=1;
				else if(j==my-1) b=my-2, flag=1;

				if(k==0)         c=1,    flag=1;
				else if(k==mz-1) c=mz-2, flag=1;

				if(flag)
				{
					source[k][j][i] = nSet(source[c][b][a]);
				}
			}
		}
	}

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->sourceU, &source);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

	DMLocalToLocalBegin(fda, ueqn->sourceU, INSERT_VALUES, ueqn->sourceU);
    DMLocalToLocalEnd  (fda, ueqn->sourceU, INSERT_VALUES, ueqn->sourceU);

	resetCellPeriodicFluxes(mesh, ueqn->sourceU, ueqn->sourceU, "vector", "localToLocal");

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
                    rhs[k][j][i].x
                    +=
                    scale *
                    (
                          central(source[k][j][i].x, source[k][j][i+1].x) * icsi[k][j][i].x +
                          central(source[k][j][i].y, source[k][j][i+1].y) * icsi[k][j][i].y +
                          central(source[k][j][i].z, source[k][j][i+1].z) * icsi[k][j][i].z
                    );
                }

                if
                (
                        isFluidJFace(k, j, i, j+1, nvert)
                )
                {
                    rhs[k][j][i].y
                    +=
                    scale *
                    (
                          central(source[k][j][i].x, source[k][j+1][i].x) * jeta[k][j][i].x +
                          central(source[k][j][i].y, source[k][j+1][i].y) * jeta[k][j][i].y +
                          central(source[k][j][i].z, source[k][j+1][i].z) * jeta[k][j][i].z
                    );
                }

                if
                (
                        isFluidKFace(k, j, i, k+1, nvert)
                )
                {
                    rhs[k][j][i].z
                    +=
                    scale *
                    (
                          central(source[k][j][i].x, source[k+1][j][i].x) * kzet[k][j][i].x +
                          central(source[k][j][i].y, source[k+1][j][i].y) * kzet[k][j][i].y +
                          central(source[k][j][i].z, source[k+1][j][i].z) * kzet[k][j][i].z
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
PetscErrorCode mapYDamping(ueqn_ *ueqn)
{
    mesh_         *mesh = ueqn->access->mesh;
    abl_          *abl  = ueqn->access->abl;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ubar;
    Cmpnts        **velMapped = abl->velMapped;
    PetscInt      numJ, numK, numI;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, p, l;
    PetscInt      imin, imax, jmin, jmax, kmin, kmax;
    PetscInt      iminSrc, jminSrc, kminSrc;
    PetscInt      numDestBounds = abl->yDampingNumPeriods;
    PetscInt      k_src_left, k_src_right;
    PetscReal     w_left, w_right;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;
    
    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, abl->uBarInstY, &ubar);

    if(ueqn->access->flags->isYDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            PetscPrintf(mesh->MESH_COMM, "Mapping lateral fringe region from sources..\n");
            for(p=0; p < abl->numSourceProc; p++)
            {
                if(abl->isdestProc[p] == 1)
                {
                    // Wait for the completion of the non-blocking broadcast performed in CorrectDampingSources function
                    MPI_Status status;
                    MPI_Wait(&(abl->mapRequest[p]), &status);

                    numI = abl->srcNumI[p];
                    numJ = abl->srcNumJ[p];
                    numK = abl->srcNumK[p];  

                    iminSrc = abl->srcMinInd[p].i;
                    jminSrc = abl->srcMinInd[p].j;
                    kminSrc = abl->srcMinInd[p].k;

                    for(l=0; l < numDestBounds; l++)
                    {
                        imin = abl->destMinInd[p][l].i;
                        imax = abl->destMaxInd[p][l].i;

                        jmin = abl->destMinInd[p][l].j;
                        jmax = abl->destMaxInd[p][l].j;

                        kmin = abl->destMinInd[p][l].k;
                        kmax = abl->destMaxInd[p][l].k;

                        if(imin == -1 || jmin == -1 || kmin == -1
                        || imax == -1 || jmax == -1 || kmax == -1)
                        {
                            continue;
                        }

                        //check that the min and max are within the processor boundaries 
                        if(kmin >= lzs && kmax <= lze && jmin >= lys && jmax <= lye && imin >= lxs && imax <= lxe)
                        {
                            for (k=kmin; k<=kmax; k++)
                            {
                                for (j=jmin; j<=jmax; j++)
                                {
                                    for (i=imin; i<=imax; i++)
                                    {
                                        k_src_left    = abl->closestKCell[k][0];
                                        k_src_right   = abl->closestKCell[k][1];
                                        w_left        = abl->wtsKCell[k][0];
                                        w_right       = abl->wtsKCell[k][1];

                                        Cmpnts uLeft  = nSet(velMapped[p][(k_src_left-kminSrc) * numJ * numI + (j-jminSrc) * numI + (i-iminSrc)]);
                                        Cmpnts uRight = nSet(velMapped[p][(k_src_right-kminSrc) * numJ * numI + (j-jminSrc) * numI + (i-iminSrc)]);

                                        ubar[k][j][i] = nSum(nScale(w_left, uLeft), nScale(w_right, uRight));

                                        // PetscPrintf(PETSC_COMM_SELF, "ubar[%ld][%ld][%ld] = %lf %lf %lf\n", k, j, i, ubar[k][j][i].x, ubar[k][j][i].y, ubar[k][j][i].z);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, abl->uBarInstY, &ubar);

    MPI_Barrier(mesh->MESH_COMM);

    DMLocalToLocalBegin (fda,  abl->uBarInstY, INSERT_VALUES, abl->uBarInstY);
    DMLocalToLocalEnd   (fda,  abl->uBarInstY, INSERT_VALUES, abl->uBarInstY);

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
    
    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // update uBar for xDampingLayer (read from inflow data)
    if(ueqn->access->flags->isXDampingActive)
    {
		DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

        if
        (
            abl->xFringeUBarSelectionType == 0 ||
            abl->xFringeUBarSelectionType == 1 ||
            abl->xFringeUBarSelectionType == 2 ||
            abl->xFringeUBarSelectionType == 4
        )
        {
            // get inflow database info pointer
            inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

            // variables to recover lapse rate above data end
            PetscReal ldataHeight = 0, gdataHeight = 0;

            // define local uBar and tBar vectors
            std::vector<std::vector<Cmpnts>>    luBarInstX(my);
            std::vector<std::vector<PetscReal>> ltBarInstX(my);

            // set it to zero
            for(j=0; j<my; j++)
            {
                luBarInstX[j].resize(mx);
                ltBarInstX[j].resize(mx);

                for(i=0; i<mx; i++)
                {
                    luBarInstX[j][i] = nSetZero();
                    ltBarInstX[j][i] = 0.0;
                }
            }

            DMDAVecGetArray(fda, mesh->lCent, &cent);

            if
            (
                abl->xFringeUBarSelectionType == 1 ||
                abl->xFringeUBarSelectionType == 2
            )
            {
                // read inflow if necessary
                readInflowU(ifPtr, ueqn->access->clock);

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
                        gdataHeight = ifPtr->avgTopLength;
                    }
                }
            }

            // update uBarInstX at this time step for this processor. These fields are defined at
            // k-face centers, so the indexing is equal to cell centers.
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // set k to the starting value of this processor. It is needed to
                    // get the z coordinate. This is only valid for cartesian meshes
                    PetscInt k_idx = lzs;

                    // steady prescribed ubar
                    if (ifPtr->typeU == 0)
                    {
                        PetscReal h = cent[k_idx][j][i].z - mesh->bounds.zmin;
                        PetscReal uMag;

                        if(h <= ifPtr->hInv)
                        {
                            uMag
                            =
                            PetscMax
                            (
                                (ifPtr->uTau/0.4)*std::log(h/ifPtr->roughness),
                                1e-5
                            );
                        }
                        else
                        {
                            uMag
                            =
                            (ifPtr->uTau/0.4)*std::log(ifPtr->hInv/ifPtr->roughness);
                        }

						luBarInstX[j][i] =  nScale(uMag, ifPtr->Udir);
                    }
                    // unsteady mapped
                    else if (ifPtr->typeU == 1)
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
                            PetscReal height = cent[k_idx][j][i].z - mesh->bounds.zmin;
                            PetscInt  IDs[2];
                            PetscReal Wg [2];

                            findInterpolationWeigthsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                            luBarInstX[j][i].x = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 ifPtr->ucat_plane[jif][iif].x +
                                                 scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 (
                                                     ifPtr->uBarAvgTopX[IDs[0]].x * Wg[0] +
                                                     ifPtr->uBarAvgTopX[IDs[1]].x * Wg[1]
                                                 );

                            luBarInstX[j][i].y = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 ifPtr->ucat_plane[jif][iif].y +
                                                 scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 (
                                                     ifPtr->uBarAvgTopX[IDs[0]].y * Wg[0] +
                                                     ifPtr->uBarAvgTopX[IDs[1]].y * Wg[1]
                                                 );

                            luBarInstX[j][i].z = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 ifPtr->ucat_plane[jif][iif].z +
                                                 scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                 (
                                                     ifPtr->uBarAvgTopX[IDs[0]].z * Wg[0] +
                                                     ifPtr->uBarAvgTopX[IDs[1]].z * Wg[1]
                                                 );
                        }
                        // index is more than nPrds times inflow points: extrapolate
                        else
                        {
                            // extrapolate along j
                            //if(j>ifPtr->n1*ifPtr->prds1) jif = ifPtr->n1;

                            // extrapolate along i
                            //if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                            luBarInstX[j][i].x = ifPtr->uBarAvgTopX[9].x;
                            luBarInstX[j][i].y = ifPtr->uBarAvgTopX[9].y;
                            luBarInstX[j][i].z = ifPtr->uBarAvgTopX[9].z;
                        }
                    }
                    // unsteady mapped interpolated
                    else if (ifPtr->typeU == 2)
                    {
                        PetscReal height = cent[k_idx][j][i].z - mesh->bounds.zmin;
                        PetscInt  IDs[2];
                        PetscReal Wg [2];

                        findInterpolationWeigthsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                        luBarInstX[j][i].x
                        =
                        scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->inflowWeights[j][i][0] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].x +
                            ifPtr->inflowWeights[j][i][1] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].x +
                            ifPtr->inflowWeights[j][i][2] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].x +
                            ifPtr->inflowWeights[j][i][3] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].x
                        ) +
                        scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->uBarAvgTopX[IDs[0]].x * Wg[0] +
                            ifPtr->uBarAvgTopX[IDs[1]].x * Wg[1]
                        );



                        luBarInstX[j][i].y
                        =
                        scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->inflowWeights[j][i][0] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].y +
                            ifPtr->inflowWeights[j][i][1] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].y +
                            ifPtr->inflowWeights[j][i][2] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].y +
                            ifPtr->inflowWeights[j][i][3] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].y
                        ) +
                        scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->uBarAvgTopX[IDs[0]].y * Wg[0] +
                            ifPtr->uBarAvgTopX[IDs[1]].y * Wg[1]
                        );

                        luBarInstX[j][i].z
                        =
                        scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->inflowWeights[j][i][0] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].z +
                            ifPtr->inflowWeights[j][i][1] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].z +
                            ifPtr->inflowWeights[j][i][2] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].z +
                            ifPtr->inflowWeights[j][i][3] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].z
                        ) +
                        scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            ifPtr->uBarAvgTopX[IDs[0]].z * Wg[0] +
                            ifPtr->uBarAvgTopX[IDs[1]].z * Wg[1]
                        );
                    }
                    // Nieuwstadt model
                    else if (ifPtr->typeU == 4)
                    {
						// get cell height
                        PetscReal h      = cent[k_idx][j][i].z - mesh->bounds.zmin;

                        luBarInstX[j][i] = nSet(NieuwstadtInflowEvaluate(ifPtr, h));
                    }

                    if(ueqn->access->flags->isTeqnActive)
                    {
                        // steady prescribed ubar
                        if (ifPtr->typeT == 0 || ifPtr->typeT == 4)
                        {
                            PetscReal b      = abl->smear * abl->gTop * abl->dInv;
                            PetscReal a      = abl->gInv - b;
                            PetscReal h      = cent[k_idx][j][i].z - mesh->bounds.zmin;
                            PetscReal etaLim = abl->hInv / abl->smear / abl->dInv;

                            // non dimensional height eta
                            PetscReal eta = (h - abl->hInv) / abl->smear / abl->dInv;

                            // below BL and capping
                            if(eta < etaLim)
                            {
                                // non dimensional functions
                                PetscReal f_eta = (std::tanh(eta) + 1.0) / 2.0;
                                PetscReal g_eta = (std::log(2.0 * std::cosh(eta)) + eta) / 2.0;

                                // potential temperature
                                ltBarInstX[j][i] = abl->tRef + a * f_eta + b * g_eta;
                            }
                            // asymptotic behavior
                            else
                            {
                                // potential temperature
                                ltBarInstX[j][i] = abl->tRef + a + b * eta;
                            }
                        }
                        // periodized mapped inflow
                        else if (ifPtr->typeT == 1)
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
                                PetscReal height = cent[k_idx][j][i].z - mesh->bounds.zmin;
                                PetscInt  IDs[2];
                                PetscReal Wg [2];

                                findInterpolationWeigthsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);


                                ltBarInstX[j][i] = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                   ifPtr->t_plane[jif][iif] +
                                                   scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                   (
                                                       ifPtr->tBarAvgTopX[IDs[0]] * Wg[0] +
                                                       ifPtr->tBarAvgTopX[IDs[1]] * Wg[1]
                                                   );
                            }
                            // index is more than nPrds times inflow points: apply lapse rate
                            else
                            {
                                PetscReal delta;

                                // compute distance for gradient addition due to height
                                if(j>ifPtr->n1*ifPtr->prds1)
                                {
                                    delta = cent[k_idx][j][i].z - gdataHeight;
                                }

                                ltBarInstX[j][i] = ifPtr->tBarAvgTopX[9] + delta * abl->gTop;
                            }
                        }
                        // interpolated periodized mapped inflow
                        else if (ifPtr->typeT == 2)
                        {
                            PetscReal delta  = PetscMax(0.0, cent[k_idx][j][i].z - gdataHeight);
                            PetscReal height = cent[k_idx][j][i].z - mesh->bounds.zmin;

                            PetscInt  IDs[2];
                            PetscReal Wg [2];
                            findInterpolationWeigthsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                            ltBarInstX[j][i]
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i] +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i] +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i] +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i]
                            ) +
                            scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->tBarAvgTopX[IDs[0]] * Wg[0] +
                                ifPtr->tBarAvgTopX[IDs[1]] * Wg[1]
                            ) +
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
            for(j=0; j<my; j++)
            {
                for(i=0; i<mx; i++)
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
                std::vector<Cmpnts>    ().swap(luBarInstX[j]);
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
                clock_     *clock     = ueqn->access->clock;
                precursor_ *precursor = ueqn->access->abl->precursor;

                // get precursor mesh info (used to stay within processor bounds)
                DM         da_p       = precursor->domain->mesh->da,
                           fda_p      = precursor->domain->mesh->fda;

                PetscInt      kStart, kEnd;

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

                if(precursor->thisProcessorInFringe)
                {
                    DMDAVecGetArray(fda_p, precursor->domain->ueqn->lUcat, &ucatP);
                    kStart = precursor->map.kStart;
                    kEnd   = precursor->map.kEnd;
                }

                DMDAVecGetArray(fda, mesh->lCent, &cent);

                for (k=zs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            // velocity values at fringe start in successor domain (two levels to interpolate)
                            if(j == abl->closestLabelsFringe[0] && k == precursor->map.kStart + 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumStart1 += ucat[k][j][i].y;
                                lcountStart1 ++;
                            }
                            else if(j==abl->closestLabelsFringe[1] && k == precursor->map.kStart + 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumStart2 += ucat[k][j][i].y;
                                lcountStart2 ++;
                            }

                            // velocity values at fringe end in successor domain (two levels to interpolate)
                            else if(j == abl->closestLabelsFringe[0] && k == precursor->map.kEnd - 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumEnd1 += ucat[k][j][i].y;
                                lcountEnd1 ++;
                            }
                            else if(j==abl->closestLabelsFringe[1] && k == precursor->map.kEnd - 1 && cent[k][j][i].y >= abl->xDampingLineSamplingYmin && cent[k][j][i].y <= abl->xDampingLineSamplingYmax)
                            {
                                lsumEnd2 += ucat[k][j][i].y;
                                lcountEnd2 ++;
                            }

                            // velocity values at reference height in precursor domain (two levels to interpolate)
                            if(precursor->thisProcessorInFringe)
                            {
                                // stay within the precursor processor bounds (boundary procs could have less cells)
                                if(j == abl->closestLabelsFringe[0] && k >= kStart && k <= kEnd)
                                {
                                    lsum1 += ucatP[k-kStart][j][i].y;
                                    lcount1 ++;
                                }
                                // stay within the precursor processor bounds (boundary procs could have less cells)
                                else if(j == abl->closestLabelsFringe[1] && k >= kStart && k <= kEnd)
                                {
                                    lsum2 += ucatP[k-kStart][j][i].y;
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
                    DMDAVecRestoreArray(fda_p, precursor->domain->ueqn->lUcat, &ucatP);
                }

                DMDAVecRestoreArray(fda, mesh->lCent, &cent);

                gsumStart1 = gsumStart1 / gcountStart1;
                gsumStart2 = gsumStart2 / gcountStart2;
                gsumEnd1   = gsumEnd1   / gcountEnd1;
                gsumEnd2   = gsumEnd2   / gcountEnd2;
                gsum1      = gsum1      / gcount1;
                gsum2      = gsum2      / gcount2;

                PetscReal vStart   = gsumStart1 * abl->levelWeightsFringe[0] + gsumStart2 * abl->levelWeightsFringe[1];
                PetscReal vEnd     = gsumEnd1   * abl->levelWeightsFringe[0] + gsumEnd2   * abl->levelWeightsFringe[1];
                PetscReal vBarPrec = gsum1      * abl->levelWeightsFringe[0] + gsum2      * abl->levelWeightsFringe[1];

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

		DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    }

    if(ueqn->access->flags->isYDampingActive)
    {
        cellIds *srcMinInd = abl->srcMinInd, *srcMaxInd = abl->srcMaxInd;
        Cmpnts  **velMapped = abl->velMapped;
        PetscInt numJ, numK, numI;
        PetscInt imin, imax, jmin, jmax, kmin, kmax;

        DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

        // update the unsteady uBar state
        if(abl->xFringeUBarSelectionType == 3)
        {

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
                                velMapped[p][(k-kmin) * numJ * numI + (j-jmin) * numI + (i-imin)] = nSet(ucat[k][j][i]);
                            }
                        }
                    }
                }

                if(abl->isdestProc[p] == 1)
                {

                    numI = abl->srcNumI[p];
                    numJ = abl->srcNumJ[p];
                    numK = abl->srcNumK[p];
                    
                    MPI_Ibcast(velMapped[p], numI * numJ * numK * 3, MPIU_REAL, abl->srcCommLocalRank[p], abl->yDamp_comm[p], &(abl->mapRequest[p]));
                }
                
            }
        }

        DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
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
    Cmpnts        ***ucont, ***ucontP, ***uBarY;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    precursor_    *precursor;
    domain_       *pdomain;
    PetscInt      kStart, kEnd, iStart, iEnd;
	PetscReal     advDampH = ueqn->access->abl->hInv - 0.5*ueqn->access->abl->dInv;

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
                kEnd   = precursor->map.kEnd;
            }
        }
    }

    if(ueqn->access->flags->isYDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            DMDAVecGetArray(fda, abl->uBarInstY, &uBarY);
            iStart = abl->iStart;
            iEnd = abl->iEnd;
        }
    }

    // z damping layer
    PetscReal alphaZ = abl->zDampingAlpha;
    PetscReal zS     = abl->zDampingStart;
    PetscReal zE     = abl->zDampingEnd;

    // y damping layer
    PetscReal alphaY = abl->yDampingAlpha;
    PetscReal yS     = abl->yDampingStart;
    PetscReal yE     = abl->yDampingEnd;
    PetscReal yD     = abl->yDampingDelta;

    // x damping layer
    PetscReal alphaX = abl->xDampingAlpha;
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;
    PetscReal xD     = abl->xDampingDelta;

    // cell center coordinates
    PetscReal x, xi, xj, xk;
    PetscReal y, yi, yj, yk;
    PetscReal z, zi, zj, zk;

    // Nordstrom viscosities
    PetscReal nud_x, nudi_x, nudj_x, nudk_x;
    PetscReal nud_y, nudi_y, nudj_y, nudk_y;

    // Rayleigh viscosities
    PetscReal nud_z, nudi_z, nudj_z, nudk_z;

    // Stipa viscosities
    PetscReal nud_x_s, nudi_x_s, nudj_x_s, nudk_x_s;

    // damping viscosity for zDamping exlusion in fringe region
    PetscReal nudI = 1.0;
    PetscReal nudJ = 1.0;
    PetscReal nudK = 1.0;

	// see if the inflow is of spread type
	PetscInt isInflowSpreadType = 0;

	if
	(
		abl->xFringeUBarSelectionType == 0 ||
		abl->xFringeUBarSelectionType == 1 ||
		abl->xFringeUBarSelectionType == 2 ||
		abl->xFringeUBarSelectionType == 4
	)
	{
		isInflowSpreadType = 1;
	}

    // loop over internal cell faces - include right boundary faces which will be periodic
    // at the beginning. Then they will be zeroed if applicable when building the SNES rhs.
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // x damping layer
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
                    nud_x_s    = viscStipa(xS, xE, xD, x,  z , advDampH);
                    nudi_x_s   = viscStipa(xS, xE, xD, xi, zi, advDampH);
                    nudj_x_s   = viscStipa(xS, xE, xD, xj, zj, advDampH);
                    nudk_x_s   = viscStipa(xS, xE, xD, xk, zk, advDampH);

                    // interpolate Stipa viscosity at cell faces
                    nudI       = central(nud_x_s, nudi_x_s);
                    nudJ       = central(nud_x_s, nudj_x_s);
                    nudK       = central(nud_x_s, nudk_x_s);

                    // X DAMPING LAYER
                    // ---------------

                    if(isInflowSpreadType)
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
					// concurrent precursor
                    else
                    {
                        PetscReal uBarContK;
                        PetscReal uBarContJ;
                        PetscReal uBarContI;

                        // note: here we can use contravariant fluxes since uBar is defined at every point

                        if
                        (
                            precursor->thisProcessorInFringe && // is this processor in the fringe?
                            k >= kStart && k <= kEnd            // is this face in the fringe?

                        )
                        {
                            uBarContI = ucontP[k-kStart][j][i].x;
                            uBarContJ = ucontP[k-kStart][j][i].y;
                            uBarContK = ucontP[k-kStart][j][i].z;
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

                // y damping layer
                if(ueqn->access->flags->isYDampingActive)
                {
                    // compute cell center y at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    y     = cent[k][j][i].y;
                    yi    = cent[k][j][i+1].y;
                    yj    = cent[k][j+1][i].y;
                    yk    = cent[k+1][j][i].y;

                    // compute Nordstrom viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_y   = viscNordstrom(alphaY, yS, yE, yD, y);
                    nudi_y  = viscNordstrom(alphaY, yS, yE, yD, yi);
                    nudj_y  = viscNordstrom(alphaY, yS, yE, yD, yj);
                    nudk_y  = viscNordstrom(alphaY, yS, yE, yD, yk);

                    if(isInflowSpreadType)
                    {

                        // retrieve from fringe region inflow slice
                        Cmpnts uBarInstY  = nSet(abl->uBarInstX[j][i]);
                        Cmpnts uBarInstYi = nSet(abl->uBarInstX[j][i+1]);
                        Cmpnts uBarInstYj = nSet(abl->uBarInstX[j+1][i]);

                        // i-fluxes
                        rhs[k][j][i].x
                        +=
                        scale * central(nud_y, nudi_y) *
                        (
                            (
                                (central(uBarInstY.x, uBarInstYi.x) - central(ucat[k][j][i].x, ucat[k][j][i+1].x)) * icsi[k][j][i].x +
                                (central(uBarInstY.y, uBarInstYi.y) - central(ucat[k][j][i].y, ucat[k][j][i+1].y)) * icsi[k][j][i].y +
                                (central(uBarInstY.z, uBarInstYi.z) - central(ucat[k][j][i].z, ucat[k][j][i+1].z)) * icsi[k][j][i].z
                            )
                        );

                        // j-fluxes
                        rhs[k][j][i].y
                        +=
                        scale * central(nud_y, nudj_y) *
                        (
                            (
                                (central(uBarInstY.x, uBarInstYj.x) - central(ucat[k][j][i].x, ucat[k][j+1][i].x)) * jeta[k][j][i].x +
                                (central(uBarInstY.y, uBarInstYj.y) - central(ucat[k][j][i].y, ucat[k][j+1][i].y)) * jeta[k][j][i].y +
                                (central(uBarInstY.z, uBarInstYj.z) - central(ucat[k][j][i].z, ucat[k][j+1][i].z)) * jeta[k][j][i].z
                            )
                        );

                        // k-fluxes
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_y, nudk_y) *
                        (
                            (
                                (uBarInstY.x - central(ucat[k][j][i].x, ucat[k+1][j][i].x)) * kzet[k][j][i].x +
                                (uBarInstY.y - central(ucat[k][j][i].y, ucat[k+1][j][i].y)) * kzet[k][j][i].y +
                                (uBarInstY.z - central(ucat[k][j][i].z, ucat[k+1][j][i].z)) * kzet[k][j][i].z
                            )
                        );
                    }
                    else
                    {
                        PetscReal uBarContK;
                        PetscReal uBarContJ;
                        PetscReal uBarContI;

                        if
                        (
                            abl->inYFringeRegionOnly && 
                            i >= iStart && i < iEnd            
                        )
                        {
                            uBarContI = central(uBarY[k][j][i].x, uBarY[k][j][i+1].x) * icsi[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k][j][i+1].y) * icsi[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k][j][i+1].z) * icsi[k][j][i].z;
                            uBarContJ = central(uBarY[k][j][i].x, uBarY[k][j+1][i].x) * jeta[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k][j+1][i].y) * jeta[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k][j+1][i].z) * jeta[k][j][i].z;
                            uBarContK = central(uBarY[k][j][i].x, uBarY[k+1][j][i].x) * kzet[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k+1][j][i].y) * kzet[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k+1][j][i].z) * kzet[k][j][i].z;
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
                        scale * central(nud_y, nudi_y) *
                        (
                            uBarContI - ucont[k][j][i].x
                        );

                        // j-fluxes
                        rhs[k][j][i].y
                        +=
                        scale * central(nud_y, nudj_y) *
                        (
                            uBarContJ - ucont[k][j][i].y
                        );

                        // k-fluxes
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_y, nudk_y) *
                        (
                            uBarContK - ucont[k][j][i].z
                        );
                    }
                }

                // z damping layer
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

                // k Left damping layer
                if(ueqn->access->flags->isKLeftRayleighDampingActive)
                {
                    // compute cell center z at i,j,k and i,j+1,k points
                    z     = (cent[k][j][i].z   - mesh->bounds.zmin);
                    zi    = (cent[k][j][i+1].z - mesh->bounds.zmin);
                    zj    = (cent[k][j+1][i].z - mesh->bounds.zmin);
                    zk    = (cent[k+1][j][i].z - mesh->bounds.zmin);

                    // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    x     = cent[k][j][i].x;
                    xi    = cent[k][j][i+1].x;
                    xj    = cent[k][j+1][i].x;
                    xk    = cent[k+1][j][i].x;

                    PetscReal hs = mesh->bounds.xmin,
                              he = mesh->bounds.xmin+abl->kLeftPatchDist;

                    // compute Cosine viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_x   = viscCosDescending(abl->kLeftDampingAlpha, hs, he, x);
                    nudi_x  = viscCosDescending(abl->kLeftDampingAlpha, hs, he, xi);
                    nudj_x  = viscCosDescending(abl->kLeftDampingAlpha, hs, he, xj);
                    nudk_x  = viscCosDescending(abl->kLeftDampingAlpha, hs, he, xk);

                    // i-fluxes
                    rhs[k][j][i].x
                    +=
                    scale * central(nud_x, nudi_x) * scaleHyperTangTop(central(z,zi), abl->kLeftDampingFilterHeight, abl->kLeftDampingFilterWidth) *
                    (
                        (
                            (abl->kLeftDampingUBar.x - central(ucat[k][j][i].x, ucat[k][j][i+1].x)) * icsi[k][j][i].x +
                            (abl->kLeftDampingUBar.y - central(ucat[k][j][i].y, ucat[k][j][i+1].y)) * icsi[k][j][i].y +
                            (abl->kLeftDampingUBar.z - central(ucat[k][j][i].z, ucat[k][j][i+1].z)) * icsi[k][j][i].z
                        )
                    );

                    // j-fluxes
                    rhs[k][j][i].y
                    +=
                    scale * central(nud_x, nudj_x) * scaleHyperTangTop(central(z,zj), abl->kLeftDampingFilterHeight, abl->kLeftDampingFilterWidth) *
                    (
                        (
                            (abl->kLeftDampingUBar.x - central(ucat[k][j][i].x, ucat[k][j+1][i].x)) * jeta[k][j][i].x +
                            (abl->kLeftDampingUBar.y - central(ucat[k][j][i].y, ucat[k][j+1][i].y)) * jeta[k][j][i].y +
                            (abl->kLeftDampingUBar.z - central(ucat[k][j][i].z, ucat[k][j+1][i].z)) * jeta[k][j][i].z
                        )
                    );

                    // k-fluxes
                    rhs[k][j][i].z
                    +=
                    scale * central(nud_x, nudk_x) * scaleHyperTangTop(central(z,zk), abl->kLeftDampingFilterHeight, abl->kLeftDampingFilterWidth) *
                    (
                        (
                            (abl->kLeftDampingUBar.x - central(ucat[k][j][i].x, ucat[k+1][j][i].x)) * kzet[k][j][i].x +
                            (abl->kLeftDampingUBar.y - central(ucat[k][j][i].y, ucat[k+1][j][i].y)) * kzet[k][j][i].y +
                            (abl->kLeftDampingUBar.z - central(ucat[k][j][i].z, ucat[k+1][j][i].z)) * kzet[k][j][i].z
                        )
                    );
                }

                // k Right damping layer
                if(ueqn->access->flags->isKRightRayleighDampingActive)
                {
                    // compute cell center z at i,j,k and i,j+1,k points
                    z     = (cent[k][j][i].z   - mesh->bounds.zmin);
                    zi    = (cent[k][j][i+1].z - mesh->bounds.zmin);
                    zj    = (cent[k][j+1][i].z - mesh->bounds.zmin);
                    zk    = (cent[k+1][j][i].z - mesh->bounds.zmin);

                    // compute cell center x at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    x     = cent[k][j][i].x;
                    xi    = cent[k][j][i+1].x;
                    xj    = cent[k][j+1][i].x;
                    xk    = cent[k+1][j][i].x;

                    PetscReal hs = mesh->bounds.xmax-abl->kRightPatchDist,
                              he = mesh->bounds.xmax;

                    // compute Cosine viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_x   = viscCosAscending(abl->kRightDampingAlpha, hs, he, x);
                    nudi_x  = viscCosAscending(abl->kRightDampingAlpha, hs, he, xi);
                    nudj_x  = viscCosAscending(abl->kRightDampingAlpha, hs, he, xj);
                    nudk_x  = viscCosAscending(abl->kRightDampingAlpha, hs, he, xk);

                    // i-fluxes
                    rhs[k][j][i].x
                    +=
                    scale * central(nud_x, nudi_x) * scaleHyperTangTop(central(z,zi), abl->kRightDampingFilterHeight, abl->kRightDampingFilterWidth) *
                    (
                        (
                            (abl->kRightDampingUBar.x - central(ucat[k][j][i].x, ucat[k][j][i+1].x)) * icsi[k][j][i].x +
                            (abl->kRightDampingUBar.y - central(ucat[k][j][i].y, ucat[k][j][i+1].y)) * icsi[k][j][i].y +
                            (abl->kRightDampingUBar.z - central(ucat[k][j][i].z, ucat[k][j][i+1].z)) * icsi[k][j][i].z
                        )
                    );

                    // j-fluxes
                    rhs[k][j][i].y
                    +=
                    scale * central(nud_x, nudj_x) * scaleHyperTangTop(central(z,zj), abl->kRightDampingFilterHeight, abl->kRightDampingFilterWidth) *
                    (
                        (
                            (abl->kRightDampingUBar.x - central(ucat[k][j][i].x, ucat[k][j+1][i].x)) * jeta[k][j][i].x +
                            (abl->kRightDampingUBar.y - central(ucat[k][j][i].y, ucat[k][j+1][i].y)) * jeta[k][j][i].y +
                            (abl->kRightDampingUBar.z - central(ucat[k][j][i].z, ucat[k][j+1][i].z)) * jeta[k][j][i].z
                        )
                    );

                    // k-fluxes
                    rhs[k][j][i].z
                    +=
                    scale * central(nud_x, nudk_x) * scaleHyperTangTop(central(z,zk), abl->kRightDampingFilterHeight, abl->kRightDampingFilterWidth) * 
                    (
                        (
                            (abl->kRightDampingUBar.x - central(ucat[k][j][i].x, ucat[k+1][j][i].x)) * kzet[k][j][i].x +
                            (abl->kRightDampingUBar.y - central(ucat[k][j][i].y, ucat[k+1][j][i].y)) * kzet[k][j][i].y +
                            (abl->kRightDampingUBar.z - central(ucat[k][j][i].z, ucat[k+1][j][i].z)) * kzet[k][j][i].z
                        )
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

    if(ueqn->access->flags->isYDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            DMDAVecRestoreArray(fda, abl->uBarInstY, &uBarY);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Coriolis(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = ueqn->access->mesh;
    abl_          *abl  = ueqn->access->abl;
    clock_        *clock = ueqn->access->clock;
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

    PetscReal     fc = abl->fc; // coriolis parameter / 2

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // The coriolis force is assumed to be present only in the x and y (wall parallel)
    // directions. Since the equations are projected along the generalized curvilinear
    // coordinates we have to dot this vector term with the face area vectors. Note that
    // if the eta axis is aligned with the z cartesian direction, the dotting will output
    // zero since eta face area vectors have zero component in x and y cartesian directions.


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
                if(isFluidIFace(k, j, i, i+1, nvert))
                {
                    rhs[k][j][i].x
                    +=
                    scale *
                    (
                        -2.0 *
                        (
                            - fc * central(ucat[k][j][i].y, ucat[k][j][i+1].y) * icsi[k][j][i].x +
                              fc * central(ucat[k][j][i].x, ucat[k][j][i+1].x) * icsi[k][j][i].y
                        )
                    );
                }

                if(isFluidJFace(k, j, i, j+1, nvert))
                {
                    rhs[k][j][i].y
                    +=
                    scale *
                    (
                        -2.0 *
                        (
                          - fc * central(ucat[k][j][i].y, ucat[k][j+1][i].y) * jeta[k][j][i].x +
                            fc * central(ucat[k][j][i].x, ucat[k][j+1][i].x) * jeta[k][j][i].y
                        )
                    );
                }

                if(isFluidKFace(k, j, i, k+1, nvert))
                {
                    rhs[k][j][i].z
                    +=
                    scale *
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

PetscErrorCode CanopyForce(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
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
    PetscReal     ***nvert, ***aj;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    // adjust canopy bounds if they are defined outside of the domain
    PetscReal     xStart  = std::max(ueqn->access->abl->xStartCanopy,mesh->bounds.xmin),
                  yStart  = std::max(ueqn->access->abl->yStartCanopy,mesh->bounds.ymin),
                  zStart  = std::max(ueqn->access->abl->zStartCanopy,mesh->bounds.zmin),
                  xEnd    = std::min(ueqn->access->abl->xEndCanopy, mesh->bounds.xmax),
                  yEnd    = std::min(ueqn->access->abl->yEndCanopy, mesh->bounds.ymax),
                  zEnd    = std::min(ueqn->access->abl->zEndCanopy, mesh->bounds.zmax);

    Cmpnts        diskNormal = ueqn->access->abl->diskDirCanopy;

    PetscReal     Hc      = (zEnd-zStart);
    PetscReal     V       = (xEnd-xStart)*(yEnd-yStart)*(zEnd-zStart);
    PetscReal     cft     = ueqn->access->abl->cftCanopy;

    PetscReal     ltotalIntU      = 0, gtotalIntU      = 0;
    PetscReal     ltotalIntThrust = 0, gtotalIntThrust = 0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(da,  mesh->lAj,  &aj);

    DMDAVecGetArray(fda, Rhs,  &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal uSource;
                PetscReal coeff = 0.0, coeff_i = 0.0, coeff_j = 0.0, coeff_k = 0.0;

                // x coord of k-face center
                PetscReal xk    = central(cent[k][j][i].x, cent[k+1][j][i].x);

                // y coord of i-face center
                PetscReal yi    = central(cent[k][j][i].y, cent[k][j][i+1].y);

                // z coord of j-face center
                PetscReal zj    = central(cent[k][j][i].z, cent[k][j][i+1].z);

                // coords of cell center
                PetscReal x     = cent[k][j][i].x;
                PetscReal y     = cent[k][j][i].y;
                PetscReal z     = cent[k][j][i].z;

                // k-face center is located in the canopy box?
                if(xk < xEnd && xk > xStart && y  < yEnd && y  > yStart && z  < zEnd && z  > zStart) coeff_k = 1;
                // i-face center is located in the canopy box?
                if(x  < xEnd && x  > xStart && yi < yEnd && yi > yStart && z  < zEnd && z  > zStart) coeff_i = 1;
                // j-face center is located in the canopy box?
                if(x  < xEnd && x  > xStart && y  < yEnd && y  > yStart && zj < zEnd && zj > zStart) coeff_j = 1;
                // cell center is located in the canopy box?
                if(x  < xEnd && x  > xStart && y  < yEnd && y  > yStart && z  < zEnd && z  > zStart) coeff   = 1;

                // source-velocity at this i-face (use central interpolation scheme)
                uSource   = nMag(centralVec(ucat[k][j][i], ucat[k][j][i+1]));

                // force per unit volume at this i-face
                PetscReal forceMagI = 0.5*cft*uSource*uSource/Hc;

                // body force vector at this k-face
                Cmpnts    forceI    = nScale(forceMagI, diskNormal);

                // body force i-flux
                if(isFluidIFace(k, j, i, i+1, nvert))
                {
                    rhs[k][j][i].x
                    +=
                    coeff_i *
                    (
                        forceI.x * icsi[k][j][i].x +
                        forceI.y * icsi[k][j][i].y +
                        forceI.z * icsi[k][j][i].z
                    );
                }

                // source-velocity at this j-face (use central interpolation scheme)
                uSource   = nMag(centralVec(ucat[k][j][i], ucat[k][j+1][i]));

                // force per unit volume at this j-face
                PetscReal forceMagJ = 0.5*cft*uSource*uSource/Hc;

                // body force vector at this j-face
                Cmpnts    forceJ    = nScale(forceMagJ, diskNormal);

                // body force j-flux
                if(isFluidJFace(k, j, i, j+1, nvert))
                {
                    rhs[k][j][i].y
                    +=
                    coeff_j *
                    (
                        forceJ.x * jeta[k][j][i].x +
                        forceJ.y * jeta[k][j][i].y +
                        forceJ.z * jeta[k][j][i].z
                    );
                }

                // source-velocity at this k-face (use central interpolation scheme)
                uSource   = nMag(centralVec(ucat[k][j][i], ucat[k+1][j][i]));

                // force per unit volume at this k-face
                PetscReal forceMagK = 0.5*cft*uSource*uSource/Hc;

                // body force vector at this k-face
                Cmpnts    forceK    = nScale(forceMagK, diskNormal);

                // body force k-flux
                if(isFluidKFace(k, j, i, k+1, nvert))
                {
                    rhs[k][j][i].z
                    +=
                    coeff_k *
                    (
                        forceK.x * kzet[k][j][i].x +
                        forceK.y * kzet[k][j][i].y +
                        forceK.z * kzet[k][j][i].z
                    );
                }

                // now do the computation on cell centers just to check
                /*
                uSource             = nMag(ucat[k][j][i]);
                PetscReal vCell     = 1.0/aj[k][j][i];
                PetscReal forceMag  = 0.5*cft*uSource*uSource/Hc;

                // integrate velocity on the canopy weghting with volume
                ltotalIntU      += coeff*uSource*vCell/V;

                // integrate the cell-thrust in the canopy
                ltotalIntThrust += coeff*forceMag*vCell;
                */
            }
        }
    }

    /*
    MPI_Allreduce(&ltotalIntU,      &gtotalIntU,      1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&ltotalIntThrust, &gtotalIntThrust, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    PetscReal totalThrust = 0.5*cft*gtotalIntU*gtotalIntU*V/Hc;

    PetscPrintf(mesh->MESH_COMM, "Canopy: actual thrust = %.2f, integrated thrust = %.2f, error % = %.2f\n", totalThrust, gtotalIntThrust, fabs(totalThrust-gtotalIntThrust)/totalThrust*100);
    */

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);

    DMDAVecRestoreArray(fda, Rhs,  &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Buoyancy(ueqn_ *ueqn, PetscReal scale)
{
    mesh_         *mesh  = ueqn->access->mesh;
    teqn_         *teqn = ueqn->access->teqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***btheta, ***gcont;
    Cmpnts        ***icsi, ***jeta, ***kzet, ***cent;
    PetscReal     ***nvert, ***tmprt, ***aj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        gravity;
                  gravity.x = 0;
                  gravity.y = 0;
                  gravity.z = -9.81;

    PetscReal     tRef;
    if(ueqn->access->flags->isAblActive) tRef = ueqn->access->abl->tRef;
    else                                 tRef = ueqn->access->constants->tRef;


    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->gCont,  &gcont);
    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  teqn->lTmprt, &tmprt);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    VecSet(ueqn->bTheta, 0.0);
    DMDAVecGetArray(fda, ueqn->bTheta,  &btheta);

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
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // interpolate temperature at cell faces
                PetscReal tempI    =  0.5*(tmprt[k][j][i] + tmprt[k][j][i+1]);
                PetscReal tempJ    =  0.5*(tmprt[k][j][i] + tmprt[k][j+1][i]);
                PetscReal tempK    =  0.5*(tmprt[k][j][i] + tmprt[k+1][j][i]);

                if(isFluidIFace(k, j, i, i+1, nvert))
                {
                    btheta[k][j][i].x
                    +=
                    scale *
                    (
                        gcont[k][j][i].x *
                        (tRef - tempI) / tRef
                        //(tempRefI - tempI) / tRef
                        //(2* tRef - tempI) / tRef
                    );
                }

                if(isFluidJFace(k, j, i, j+1, nvert))
                {
                    btheta[k][j][i].y
                    +=
                    scale *
                    (
                        gcont[k][j][i].y *
                        (tRef - tempJ) / tRef
                        //(tempRefJ - tempJ) / tRef
                        //(2* tRef - tempJ) / tRef
                    );
                }

                if(isFluidKFace(k, j, i, k+1, nvert))
                {
                    btheta[k][j][i].z
                    +=
                    scale *
                    (
                        gcont[k][j][i].z *
                        (tRef - tempK) / tRef
                        //(tempRefK - tempK) / tRef
                        //(2* tRef - tempK) / tRef
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->bTheta,  &btheta);
    DMDAVecRestoreArray(da,  teqn->lTmprt,  &tmprt);
    DMDAVecRestoreArray(da,  mesh->lNvert,  &nvert);
    DMDAVecRestoreArray(fda, ueqn->gCont,   &gcont);
    DMDAVecRestoreArray(fda, mesh->lCent,   &cent);

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
            // k-left boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux

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
                    // cumulate flux

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
            // j-left boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux

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
            // j-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux

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
            // i-left boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux

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
            // i-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux

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

    // PetscPrintf(mesh->MESH_COMM, "After correction fluxin = %lf, fluxout = %lf\n", FluxIn, FluxOut);

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
    flags_           *flags  = ueqn->access->flags;
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
    Cmpnts           ***cent;

    PetscReal        ***nvert, ***lnu_t;

    Cmpnts           ***div1,  ***div2,  ***div3;                               // divergence & cumulative fluxes
    Cmpnts           ***visc1, ***visc2, ***visc3;
    Cmpnts           ***viscIBM1, ***viscIBM2, ***viscIBM3;                              // viscous terms                             // viscous terms
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

    DMDAVecGetArray(fda, mesh->lCent,  &cent);

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

    if(flags->isIBMActive)
    {
        DMDAVecGetArray(fda, ueqn->lViscIBM1, &viscIBM1);
        DMDAVecGetArray(fda, ueqn->lViscIBM2, &viscIBM2);
        DMDAVecGetArray(fda, ueqn->lViscIBM3, &viscIBM3);
    }

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
    }

    // damping viscosity for fringe region advection damping
    PetscReal nuD;

    // fringe region parameters (set only if active)
    PetscReal xS;
    PetscReal xE;
    PetscReal xDS;
    PetscReal xDE;

    PetscInt  advectionDamping = 0;
	PetscReal advDampH = 0;

    if(ueqn->access->flags->isAdvectionDampingActive)
    {
        xS     = ueqn->access->abl->advDampingStart;
        xE     = ueqn->access->abl->advDampingEnd;
        xDE    = ueqn->access->abl->advDampingDeltaEnd;
        xDS    = ueqn->access->abl->advDampingDeltaStart;

        advectionDamping = 1;

        advDampH = ueqn->access->abl->hInv - 0.5*ueqn->access->abl->dInv;
    }
    else
    {
        nuD = 1.0;
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
                        div1[k][j][i] = nScale
                        (
                            - ucont[k][j][i].x,
                            weno3Vec
                            (
                                ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                ucont[k][j][i].x
                            )
                        );
                    }
                    else
                    {
                        // central scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k][j][i+1]
                                )
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                centralUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    ucont[k][j][i].x, limiter[k][j][i].x
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][j][iL ] * nMag(csi[k][j][iL ]));
                            PetscReal d1 = 1.0 / (aj[k][j][i  ] * nMag(csi[k][j][i  ]));
                            PetscReal d2 = 1.0 / (aj[k][j][i+1] * nMag(csi[k][j][i+1]));
                            PetscReal d3 = 1.0 / (aj[k][j][iR ] * nMag(csi[k][j][iR ]));

                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                wCentralUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].x, limiter[k][j][i].x
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                quadraticUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    ucont[k][j][i].x
                                )
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

                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                    // wall model i-left patch
                    if
                    (
                            (mesh->boundaryU.iLeft=="velocityWallFunction"  && i==0) ||
                            (mesh->boundaryU.iRight=="velocityWallFunction" && i==mx-2)
                    )
                    {
                        visc1[k][j][i].x = - ueqn->iLWM->tauWall.x[k-zs][j-ys];
                        visc1[k][j][i].y = - ueqn->iLWM->tauWall.y[k-zs][j-ys];
                        visc1[k][j][i].z = - ueqn->iLWM->tauWall.z[k-zs][j-ys];

                        nuEff = nu;
                    }
                    // slip boundary condition on U (set nuEff to 0)
                    else if
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

                    if(isIBMFluidIFace(k, j, i, i+1, nvert))
                    {
                        if(ueqn->access->ibm->wallShearOn)
                        {
                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc1[k][j][i] = nSet(viscIBM1[k][j][i]);
                            }
                            else if(isIBMFluidCell(k, j, i+1, nvert))
                            {
                                visc1[k][j][i] = nSet(viscIBM1[k][j][i+1]);
                            }

                            nuEff = 0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

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
                        div2[k][j][i] = nScale
                        (
                            - ucont[k][j][i].y,
                            weno3Vec
                            (
                                ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                ucont[k][j][i].y
                            )
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k][j+1][i]
                                )
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                centralUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    ucont[k][j][i].y, limiter[k][j][i].y
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][jL ][i] * nMag(eta[k][jL ][i]));
                            PetscReal d1 = 1.0 / (aj[k][j  ][i] * nMag(eta[k][j  ][i]));
                            PetscReal d2 = 1.0 / (aj[k][j+1][i] * nMag(eta[k][j+1][i]));
                            PetscReal d3 = 1.0 / (aj[k][jR ][i] * nMag(eta[k][jR ][i]));

                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                wCentralUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].y, limiter[k][j][i].y
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                quadraticUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    ucont[k][j][i].y
                                )
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

                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                    // wall model j-left patch
                    if
                    (
                            (mesh->boundaryU.jLeft=="velocityWallFunction" && j==0) ||
                            (mesh->boundaryU.jRight=="velocityWallFunction" && j==my-2)
                    )
                    {
                        visc2[k][j][i].x = - ueqn->jLWM->tauWall.x[k-zs][i-xs];
                        visc2[k][j][i].y = - ueqn->jLWM->tauWall.y[k-zs][i-xs];
                        visc2[k][j][i].z = - ueqn->jLWM->tauWall.z[k-zs][i-xs];

                        nuEff = nu;
                    }
                    // slip boundary condition on U (set nuEff to 0)
                    else if
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

                    if(isIBMFluidJFace(k, j, i, j+1, nvert))
                    {
                        if(ueqn->access->ibm->wallShearOn)
                        {

                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc2[k][j][i] = nSet(viscIBM2[k][j][i]);
                            }
                            else if(isIBMFluidCell(k, j+1, i, nvert))
                            {
                                visc2[k][j][i] = nSet(viscIBM2[k][j+1][i]);
                            }
                           
                            nuEff = 0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

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
                    // inviscid flow or weno3Div
                    if(ueqn->inviscid || ueqn->weno3Div)
                    {
                        div3[k][j][i] = nScale
                        (
                            - ucont[k][j][i].z,
                            weno3Vec
                            (
                                ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                ucont[k][j][i].z
                            )
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k+1][j][i]
                                )
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                centralUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    ucont[k][j][i].z, limiter[k][j][i].z
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[kL ][j][i] * nMag(zet[kL ][j][i]));
                            PetscReal d1 = 1.0 / (aj[k  ][j][i] * nMag(zet[k  ][j][i]));
                            PetscReal d2 = 1.0 / (aj[k+1][j][i] * nMag(zet[k+1][j][i]));
                            PetscReal d3 = 1.0 / (aj[kR ][j][i] * nMag(zet[kR ][j][i]));

                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                wCentralUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].z, limiter[k][j][i].z
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                quadraticUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    ucont[k][j][i].z
                                )
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
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

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

                    if(isIBMFluidKFace(k, j, i, k+1, nvert))
                    {
                        if(ueqn->access->ibm->wallShearOn)
                        {
                            visc3[k][j][i] = nSet(viscIBM3[k][j][i]);

                            if(isIBMFluidCell(k, j, i, nvert))
                            {
                                visc3[k][j][i] = nSet(viscIBM3[k][j][i]);
                            }
                            else if(isIBMFluidCell(k+1, j, i, nvert))
                            {
                                visc3[k][j][i] = nSet(viscIBM3[k+1][j][i]);
                            }

                            nuEff = 0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

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

    if(flags->isIBMActive)
    {
        DMDAVecRestoreArray(fda, ueqn->lViscIBM1, &viscIBM1);
        DMDAVecRestoreArray(fda, ueqn->lViscIBM2, &viscIBM2);
        DMDAVecRestoreArray(fda, ueqn->lViscIBM3, &viscIBM3);
    }

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
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv2,   ueqn->lDiv2,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv3,   ueqn->lDiv3,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc1,  ueqn->lVisc1,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc2,  ueqn->lVisc2,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc3,  ueqn->lVisc3,  "localToLocal");

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

				// compute Stipa viscosity at i,j,k,  point
				if(advectionDamping)
				{
					PetscReal height = cent[k][j][i].z - mesh->bounds.zmin;
					nuD              = viscStipaDelta(xS, xE, xDS, xDE, cent[k][j][i].x, height, advDampH);
					fp[k][j][i].z    = nuD * fp[k][j][i].z;
				}

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

            }
        }
    }

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

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

    if(ueqn->access->flags->isIBMActive)
    {
        UpdateImmersedBCs(ueqn->access->ibm);
    }

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    PetscReal dt = clock->dt;

    // add pressure gradient term
    VecAXPY(Rhs, -1.0, ueqn->dP);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->coriolisActive)
        {
            Coriolis(ueqn, Rhs, 1.0);
        }
    }

    // add side force term
    if(ueqn->access->flags->isCanopyActive)
    {
        CanopyForce(ueqn, Rhs, 1.0);
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
            VecAXPY(Rhs,  1.0, ueqn->bTheta);
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
        ueqn->access->flags->isYDampingActive ||
        ueqn->access->flags->isZDampingActive ||
        ueqn->access->flags->isKLeftRayleighDampingActive ||
        ueqn->access->flags->isKRightRayleighDampingActive
    )
    {
        dampingSourceU(ueqn, Rhs, 1.0);
    }


    // multiply for dt
    VecScale(Rhs, dt);

    // add driving source terms after as it is not scaled by 1/dt
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->controllerActive)
        {
            sourceU(ueqn, Rhs, 1.0);
        }
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

    if(ueqn->access->flags->isIBMActive)
    {
        UpdateImmersedBCs(ueqn->access->ibm);
    }

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // initialize the rhs vector
    VecSet(ueqn->Rhs, 0.0);

    // add pressure gradient term
    VecAXPY(ueqn->Rhs, -1.0, ueqn->dP);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->coriolisActive)
        {
            Coriolis(ueqn, ueqn->Rhs, 1.0);
        }
    }

    // add side force term
    if(ueqn->access->flags->isCanopyActive)
    {
        CanopyForce(ueqn, ueqn->Rhs, 1.0);
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
            VecAXPY(ueqn->Rhs,  1.0, ueqn->bTheta);
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
        ueqn->access->flags->isZDampingActive ||
        ueqn->access->flags->isKLeftRayleighDampingActive ||
        ueqn->access->flags->isKRightRayleighDampingActive
    )
    {
        dampingSourceU(ueqn, ueqn->Rhs, 1.0);
    }

    // source term is pre-scaled by dt, so here we have to divide it again since
    // the rhs is not multiplied by dt in this function
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->controllerActive)
        {
            sourceU(ueqn, ueqn->Rhs, 1.0 / clock->dt);
        }
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
