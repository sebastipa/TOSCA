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

    clock_        *clock = domain[0].clock;

    // save old time step
    clock->dtOld = clock->dt;

    for(PetscInt d=0; d<nDomains; d++)
    {
        flags_        *flags = domain[d].access.flags;
        acquisition_  *acquisition = domain[d].acquisition;
        mesh_         *mesh  = domain[d].mesh;
        ueqn_         *ueqn  = domain[d].ueqn;
        DM            da = mesh->da, fda = mesh->fda;
        DMDALocalInfo info = mesh->info;
        PetscInt      xs = info.xs, xe = info.xs + info.xm;
        PetscInt      ys = info.ys, ye = info.ys + info.ym;
        PetscInt      zs = info.zs, ze = info.zs + info.zm;
        PetscInt      mx = info.mx, my = info.my, mz = info.mz;

        Cmpnts        ***ucat;
        Cmpnts        ***csi, ***eta, ***zet;
        PetscReal     ***nvert, ***aj;

        PetscInt           i, j, k;
        PetscInt           lxs, lxe, lys, lye, lzs, lze;

        lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
        lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
        lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

        PetscReal lmaxU      = 0.0;
        PetscReal ldx        = 1e10;
        PetscReal ldi_min    = 1e10, ldj_min = 1e10, ldk_min = 1e10;
        PetscReal ldxByU_min = 1e10;

        cfl = PetscMin(cfl, clock->cfl);

        DMDAVecGetArray(fda, ueqn->Ucat,   &ucat);
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
                    // contravariant velocity magnitude
                    PetscReal Umag = nMag(ucat[k][j][i]);

                    // minimum overall velocity in this proc
                    lmaxU = PetscMax(Umag, lmaxU);

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

                    // compute dxByU
                    PetscReal ldxByU = PetscMin(ldi, PetscMin(ldj, ldk)) / Umag;

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

        lmaxU      = PetscMax(lmaxU, maxU);
        ldxByU_min = PetscMin(ldxByU_min, dxByU_min);

        MPI_Allreduce(&lmaxU,        &maxU,      1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);
        MPI_Allreduce(&ldxByU_min,   &dxByU_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);

        DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);
        DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
        DMDAVecRestoreArray(fda, mesh->lEta, &eta);
        DMDAVecRestoreArray(fda, mesh->lZet, &zet);
        DMDAVecRestoreArray(da,  mesh->lAj, &aj);
        DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

        // output fields
        if(flags->isAdjustableTime)
        {
            // 2. takes the local ratio
            clock->dt = clock->cfl * dxByU_min;

            PetscReal timeStart;
            PetscReal timeInterval;

            timeStart    = clock->startTime;
            timeInterval = domain[d].io->timeInterval;

            timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);

            // averaged fields
            if(domain[d].io->averaging)
            {
                timeStart    = domain[d].io->avgStartTime;
                timeInterval = domain[d].io->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
            }

            // phase averaged fields
            if(domain[d].io->phaseAveraging)
            {
                timeStart    = domain[d].io->phAvgStartTime;
                timeInterval = domain[d].io->phAvgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
            }

            // ke budgets
            if(domain[d].io->keBudgets)
            {
                timeStart    = acquisition->keBudFields->avgStartTime;
                timeInterval = acquisition->keBudFields->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
            }

            // ibm write pressure force
            if(flags->isIBMActive)
            {
                ibm_ *ibm = domain[d].ibm;

                if(domain[d].io->writePForce)
                {
                    timeStart    = ibm->startTime;
                    timeInterval = ibm->writePrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                }
            }

            if(flags->isAquisitionActive)
            {
                // 3LM averaging for background domain only
                if(acquisition->isAverage3LMActive)
                {
                    timeStart    = domain[0].acquisition->LM3->avgStartTime;
                    timeInterval = domain[0].acquisition->LM3->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                }

                // ABL perturbations
                if(acquisition->isPerturbABLActive)
                {
                    timeStart    = domain[0].acquisition->perturbABL->avgStartTime;
                    timeInterval = domain[0].acquisition->perturbABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                }

                // ABL averaging for background domain only
                if(acquisition->isAverageABLActive)
                {
                    timeStart    = domain[0].acquisition->statisticsABL->avgStartTime;
                    timeInterval = domain[0].acquisition->statisticsABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                }

                if(acquisition->isProbesActive)
                {
                    for(PetscInt r=0; r<acquisition->probes->nRakes; r++)
                    {
                        timeStart    = acquisition->probes->rakes[r].timeStart;
                        timeInterval = acquisition->probes->rakes[r].timeInterval;

                        timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
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
                        }
                    }
                }

                if(flags->isWindFarmActive)
                {
                    if(domain[d].farm->intervalType == "adjustableTime")
                    {
                        timeStart    = domain[d].farm->timeStart;
                        timeInterval = domain[d].farm->timeInterval;

                        timeStepSet(clock, timeStart, timeInterval, dxByU_min, flag, cfl);
                    }
                }
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

                    clock->dt = std::min(clock->dt, dtFarm);
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

                            PetscReal maxSpeed = ibmRot->maxR * ibmRot->angSpeed;

                            PetscReal dtIBM = dxByU_min*maxU / maxSpeed;

                            clock->dt = std::min(clock->dt, dtIBM);
                        }
                    }
                }
            }

        }
        else
        {
            if(clock->it>clock->itStart) clock->cfl =  clock->dt / dxByU_min;

            // added by Arjun for fixed time step as a last control
            if(clock->cfl > 1.0)
            {
                clock->dt = clock->dt/2.0;
            }
        }
    }

    // prevent time step from increasing too fast
    if(flags->isAdjustableTime)
    {
        if(clock->dt > 1.5 * clock->dtOld)
        {
            clock->dt = 1.5 * clock->dtOld;
        }
    }

    // avoid too small time step (will not write, have to fix this)
    if (clock->dt < 1e-10)
    {
        clock->dt = clock->startDt;
    }

    // make sure to hit last time
    if(clock->time + clock->dt > clock->endTime)
    {
        clock->dt = clock->endTime - clock->time;
    }

    clock->time = clock->time + clock->dt;

    cfl = clock->dt / dxByU_min;

    PetscPrintf(PETSC_COMM_WORLD, "\n\nTime: %lf\n\n", clock->time);
    PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL = %lf, uMax = %.6f, dt = %.6f, dtMaxCFL = %.6f, adjust due to write flag: %ld\n", clock->it, cfl, maxU, clock->dt, dxByU_min, flag);

    return(0);
}
