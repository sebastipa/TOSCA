//! \file  clock.c
//! \brief Contains simulation time function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

PetscErrorCode adjustTimeStep (domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;
    flags_   *flags   = domain[0].access.flags;

    PetscInt  flag       = 0;
    PetscReal cfl        = 1e10;
    PetscReal maxU       = 0.0;
    PetscReal dxByU_min  = 1e10;
    cellIds   maxUCell;

    clock_        *clock = domain[0].clock;

    // save old time step
    clock->dtOld = clock->dt;

    for(PetscInt d=0; d<nDomains; d++)
    {
        acquisition_  *acquisition = domain->acquisition;

        // set cfl
        cfl = PetscMin(cfl, clock->cfl);

        // time step imposed by the flow on this domain
        timeStepInfo(&domain[d], clock, dxByU_min, maxU, maxUCell);

        // time step imposed by the precursor on this domain (must be corrected to syncronize time steps)
        if(domain[d].flags.isConcurrentPrecursorActive)
        {
            //timeStepInfo(domain[d].abl->precursor->domain, clock, dxByU_min, maxU, maxUCell);
        }

        // try to guess a uniform predicted time step for the acquisition
        PetscReal predictedDt;

        // output fields
        if(flags->isAdjustableTime)
        {
            // 2. takes the local ratio
            clock->dt   = clock->cfl * dxByU_min;

            PetscReal timeStart;
            PetscReal timeInterval;

            timeStart    = clock->startTime;
            timeInterval = domain[d].io->timeInterval;

            MPI_Barrier(domain[d].mesh->MESH_COMM);
            PetscPrintf(PETSC_COMM_WORLD, "Domain %ld, barrier 1\n", d);

            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
            predictedDt  = currentDistanceToWriteTime(clock, timeStart, timeInterval);

            MPI_Barrier(domain[d].mesh->MESH_COMM);
            PetscPrintf(PETSC_COMM_WORLD, "Domain %ld, barrier 2\n", d);

            // averaged fields
            if(domain[d].io->averaging)
            {
                timeStart    = domain[d].io->avgStartTime;
                timeInterval = domain[d].io->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            // phase averaged fields
            if(domain[d].io->phaseAveraging)
            {
                timeStart    = domain[d].io->phAvgStartTime;
                timeInterval = domain[d].io->phAvgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            // ke budgets
            if(domain[d].io->keBudgets)
            {
                timeStart    = acquisition->keBudFields->avgStartTime;
                timeInterval = acquisition->keBudFields->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
            }

            if(flags->isIBMActive)
            {
                ibm_ *ibm = domain[d].ibm;

                if(domain[d].ibm->intervalType == "adjustableTime")
                {
                    timeStart    = ibm->timeStart;
                    timeInterval = ibm->timeInterval;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
            }

            if(flags->isWindFarmActive)
            {
                if(domain[d].farm->intervalType == "adjustableTime")
                {
                    timeStart    = domain[d].farm->timeStart;
                    timeInterval = domain[d].farm->timeInterval;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
            }
            MPI_Barrier(domain[d].mesh->MESH_COMM);
            PetscPrintf(PETSC_COMM_WORLD, "Domain %ld, barrier 3\n", d);

            // catalyst
            if(flags->isPvCatalystActive)
            {
                #ifdef USE_CATALYST
                if(domain[d].io->outputTypeCatalyst == "adjustableTime")
                {
                    timeStart    = domain[d].io->startTimeCatalyst;
                    timeInterval = domain[d].io->timeIntervalCatalyst;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }
                #endif
            }

            if(flags->isAquisitionActive)
            {
                // 3LM averaging for background domain only
                if(acquisition->isAverage3LMActive)
                {
                    timeStart    = domain[0].acquisition->LM3->avgStartTime;
                    timeInterval = domain[0].acquisition->LM3->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                // ABL perturbations
                if(acquisition->isPerturbABLActive)
                {
                    timeStart    = domain[0].acquisition->perturbABL->avgStartTime;
                    timeInterval = domain[0].acquisition->perturbABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                // ABL averaging for background domain only
                if(acquisition->isAverageABLActive)
                {
                    timeStart    = domain[0].acquisition->statisticsABL->avgStartTime;
                    timeInterval = domain[0].acquisition->statisticsABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                }

                if(acquisition->isProbesActive)
                {
                    for(PetscInt r=0; r<acquisition->probes->nRakes; r++)
                    {
                        if(acquisition->probes->rakes[r].intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->probes->rakes[r].timeStart;
                            timeInterval = acquisition->probes->rakes[r].timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));

                            if(acquisition->probes->allSameIO) break;
                        }
                    }
                }

                if(acquisition->isSectionsActive)
                {
                    // iSections
                    if(acquisition->iSections->available)
                    {
                        if(acquisition->iSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->iSections->timeStart;
                            timeInterval = acquisition->iSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }

                    // iSections
                    if(acquisition->jSections->available)
                    {
                        if(acquisition->jSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->jSections->timeStart;
                            timeInterval = acquisition->jSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }

                    // kSections
                    if(acquisition->kSections->available)
                    {
                        if(acquisition->kSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->kSections->timeStart;
                            timeInterval = acquisition->kSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                            predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                        }
                    }
                }

                MPI_Barrier(domain[d].mesh->MESH_COMM);
                PetscPrintf(PETSC_COMM_WORLD, "Domain %ld, barrier 4\n", d);
            }

            // concurrent precursor
            if(domain[d].flags.isConcurrentPrecursorActive)
            {
                /*
                if(domain[d].abl->precursor->domain->flags.isAquisitionActive)
                {
                    if(domain[d].abl->precursor->domain->acquisition->isAverageABLActive)
                    {
                        timeStart    = domain[d].abl->precursor->domain->acquisition->statisticsABL->avgStartTime;
                        timeInterval = domain[d].abl->precursor->domain->acquisition->statisticsABL->avgPrd;

                        timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                        predictedDt  = gcd(predictedDt, currentDistanceToWriteTime(clock, timeStart, timeInterval));
                    }
                }
                */
            }

            // set time step as the gcd of all constraints
            if(clock->it == 0)
            {
                // scale max uniform dt due to acquisition so that it complies the CFL
                while(predictedDt / dxByU_min > clock->cfl)
                {
                    predictedDt /= 2.0;
                }

                // save
                clock->acquisitionDt = predictedDt;
            }

            // check if must limit dt due to ALM rotation
            if(flags->isWindFarmActive)
            {
                farm_ *farm = domain[d].farm;

                // set checkCFL flag
                if(clock->it == 0)
                {
                    // loop over each wind turbine
                    for(PetscInt t=0; t<farm->size; t++)
                    {
                        if((*farm->turbineModels[t]) == "ALM")
                        {
                            farm->checkCFL = 1;
                        }
                    }
                }

                if(farm->checkCFL)
                {
                    computeMaxTipSpeed(domain[d].farm);

                    PetscReal dtFarm = clock->dxMin / farm->maxTipSpeed;

                    clock->dt    = std::min(clock->dt, dtFarm);
                }
            }

            // check cfl limit for rotating ibm
            if(flags->isIBMActive)
            {
                ibm_ *ibm = domain[d].ibm;

                for (PetscInt i=0; i < ibm->numBodies; i++)
                {
                    ibmObject   *ibmBody = ibm->ibmBody[i];

                    if(ibm->dynamic)
                    {
                        if(ibmBody->bodyMotion == "rotation")
                        {
                            ibmRotation *ibmRot  = ibmBody->ibmRot;

                            PetscReal maxSpeed   = ibmRot->maxR * ibmRot->angSpeed;

                            PetscReal dtIBM      = dxByU_min*maxU / maxSpeed;

                            clock->dt    = std::min(clock->dt, dtIBM);
                        }
                    }
                }
            }
        }
        else
        {
            if(clock->it>clock->itStart) clock->cfl =  clock->dt / dxByU_min;

            // added by Arjun for fixed time step as a last control
            /*
			if(clock->cfl > 1.0)
            {
                clock->dt = clock->dt/2.0;
            }
			*/
        }
    }

    // prevent time step from increasing too fast
    if(flags->isAdjustableTime)
    {
        if(clock->dt > 1.5 * clock->dtOld)
        {
            clock->dt    = 1.5 * clock->dtOld;
        }
    }

    // avoid too small time step (will not write, have to fix this)
    if (clock->dt < 1e-10)
    {
        clock->dt = clock->startDt;
    }

    // discard all previous changes if time step is fixed as acquistion gcd
    if(flags->isAdjustableTime == 2)
    {
        if(clock->acquisitionDt / dxByU_min > 1.0)
        {
            clock->acquisitionDt /= 2.0;
        }

        clock->dt = clock->acquisitionDt;
    }

    // make sure to hit last time
    if(clock->time + clock->dt > clock->endTime)
    {
        clock->dt = clock->endTime - clock->time;
    }

    if(flags->isAdjustableTime)
    {
        clock->time = clock->time + clock->dt;
    }
    else 
    {
        //to exit when reaching end time 
        if(clock->time + clock->dt > clock->endTime)
        {
            clock->time = clock->time + clock->dt;
        }
        else 
        {
            // this ensures there is no floating point addition error
            clock->time = clock->startTime + (clock->it + 1) * clock->dt;
        }
    }
    
    cfl = clock->dt / dxByU_min;

    PetscPrintf(PETSC_COMM_WORLD, "\n\nTime: %lf\n\n", clock->time);

    if(clock->it==clock->itStart)
    {
        PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL = %lf, uMax = %.6f, dt = %.6f, dtMaxCFL = %.6f, adjust due to write flag: %ld\n", clock->it, cfl, maxU, clock->dt, dxByU_min, flag);
    }
    else
    {
        PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL = %lf, uMax = %.6f (i,j,k = %ld, %ld, %ld), dt = %.6f, dtMaxCFL = %.6f, adjust due to write flag: %ld\n", clock->it, cfl, maxU, maxUCell.i, maxUCell.j, maxUCell.k, clock->dt, dxByU_min, flag);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode timeStepInfo(domain_ *domain, clock_ *clock, PetscReal &dxByU_min, PetscReal &maxU, cellIds &maxUCell)
{
    flags_        *flags = domain->access.flags;
    mesh_         *mesh  = domain->mesh;
    ueqn_         *ueqn  = domain->ueqn;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucat, ***ucont;
    Cmpnts        ***csi, ***eta, ***zet;
    PetscReal     ***nvert, ***aj;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    PetscReal lmaxU      = 0.0,  gmaxU      = 0.0;
    PetscReal ldx        = 1e10;
    PetscReal ldi_min    = 1e10, ldj_min    = 1e10, ldk_min = 1e10;
    PetscReal ldxByU_min = 1e10, gdxByU_min = 1e10;
    cellIds   lmaxUCell; lmaxUCell.i = 0; lmaxUCell.j = 0; lmaxUCell.k = 0;
    cellIds   gmaxUCell; gmaxUCell.i = 0; gmaxUCell.j = 0; gmaxUCell.k = 0;

    DMDAVecGetArray(fda, ueqn->Ucat,   &ucat);
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    // time step due to flow restrictions
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // cartesian velocity magnitude
                PetscReal Umag = nMag(ucat[k][j][i]);

                // maximum overall velocity in this proc
                if(Umag > lmaxU)
                {
                    lmaxU = Umag;
                    lmaxUCell.i = i;
                    lmaxUCell.j = j;
                    lmaxUCell.k = k;
                }

                // compute cell sizes
                PetscReal ldi = 1./aj[k][j][i]/nMag(csi[k][j][i]);
                PetscReal ldj = 1./aj[k][j][i]/nMag(eta[k][j][i]);
                PetscReal ldk = 1./aj[k][j][i]/nMag(zet[k][j][i]);

                // minimum overall GCC sizes in this proc
                if(clock->it == clock->itStart)
                {
                    ldi_min = PetscMin(ldi_min, ldi);
                    ldj_min = PetscMin(ldj_min, ldj);
                    ldk_min = PetscMin(ldk_min, ldk);
                }

                // compute directional velocity in curvilinear coordinates
                PetscReal Vi = fabs(0.5 * (ucont[k][j][i].x + ucont[k][j][i-1].x) / nMag(csi[k][j][i])),
                          Vj = fabs(0.5 * (ucont[k][j][i].y + ucont[k][j-1][i].y) / nMag(eta[k][j][i])),
                          Vk = fabs(0.5 * (ucont[k][j][i].z + ucont[k-1][j][i].z) / nMag(zet[k][j][i]));

                // compute dxByU
                PetscReal ldxByU = PetscMin(ldi/Vi, PetscMin(ldj/Vj, ldk/Vk));

                if(isFluidCell(k, j, i, nvert))
                {
                    // min overall cell size
                    if(clock->it == clock->itStart)
                    {
                        ldx = PetscMin(ldi_min, PetscMin(ldj_min, ldk_min));
                    }

                    // min overall dxByU
                    ldxByU_min = PetscMin(ldxByU_min, ldxByU);
                }
            }
        }
    }

    if(clock->it == clock->itStart)
    {
        ldx = PetscMin(ldx, clock->dxMin);
        MPI_Allreduce(&ldx, &clock->dxMin, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);
    }

    // max in this domain
    MPI_Allreduce(&lmaxU,        &gmaxU,      1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);
    MPI_Allreduce(&ldxByU_min,   &gdxByU_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

    // check if maxU is in this processor and print where it is
    if(lmaxU != gmaxU)
    {
        lmaxUCell.i = 0; lmaxUCell.j = 0; lmaxUCell.k = 0;
    }

    MPI_Allreduce(&lmaxUCell, &gmaxUCell, 3, MPIU_INT, MPIU_MAX, mesh->MESH_COMM);

    // max in all domains and save ids
    if(gmaxU > maxU)
    {
        maxU     = gmaxU;
        maxUCell = gmaxUCell;
    }
    if(gdxByU_min < dxByU_min)
    {
        dxByU_min = gdxByU_min;
    }

    DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    return(0);
}
