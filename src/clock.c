//! \file  clock.c
//! \brief Contains simulation time function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

PetscErrorCode adjustTimeStep (domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;

    PetscInt  flag   = 0;
    PetscReal cfl    = 1e10;
    PetscReal dx_min = 1e10;
    PetscReal maxU   = 0.0;

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

        PetscReal lmaxU   = 0.0;
        PetscReal ldx     = 1e10;
        PetscReal ldi_min = 1e10, ldj_min = 1e10, ldk_min = 1e10;

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
                    PetscReal Umag =
                    sqrt
                    (
                        ucat[k][j][i].x*ucat[k][j][i].x +
                        ucat[k][j][i].y*ucat[k][j][i].y +
                        ucat[k][j][i].z*ucat[k][j][i].z
                    );

                    lmaxU = PetscMax(Umag, lmaxU);

                    PetscReal ldi = 1./aj[k][j][i]/sqrt( csi[k][j][i].x*csi[k][j][i].x + csi[k][j][i].y*csi[k][j][i].y + csi[k][j][i].z*csi[k][j][i].z );
                    PetscReal ldj = 1./aj[k][j][i]/sqrt( eta[k][j][i].x*eta[k][j][i].x + eta[k][j][i].y*eta[k][j][i].y + eta[k][j][i].z*eta[k][j][i].z );
                    PetscReal ldk = 1./aj[k][j][i]/sqrt( zet[k][j][i].x*zet[k][j][i].x + zet[k][j][i].y*zet[k][j][i].y + zet[k][j][i].z*zet[k][j][i].z );

                    ldi_min = PetscMin ( ldi_min, ldi );
                    ldj_min = PetscMin ( ldj_min, ldj );
                    ldk_min = PetscMin ( ldk_min, ldk );

                    if(isFluidCell(k, j, i, nvert))
                    {
                        ldx = PetscMin(ldi_min, PetscMin(ldj_min, ldk_min));
                    }
                }
            }
        }

        ldx   = PetscMin(ldx, dx_min);
        lmaxU = PetscMax(lmaxU, maxU);

        MPI_Allreduce(&ldx,   &dx_min, 1, MPIU_REAL, MPIU_MIN, mesh->MESH_COMM);
        MPI_Allreduce(&lmaxU, &maxU, 1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);

        DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);
        DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
        DMDAVecRestoreArray(fda, mesh->lEta, &eta);
        DMDAVecRestoreArray(fda, mesh->lZet, &zet);
        DMDAVecRestoreArray(da, mesh->lAj, &aj);
        DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

        // output fields
        if(flags->isAdjustableTime)
        {
            clock->dt = clock->cfl * dx_min / maxU;

            PetscReal timeStart;
            PetscReal timeInterval;

            timeStart    = clock->startTime;
            timeInterval = domain[d].io->timeInterval;

            timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);

            // averaged fields
            if(domain[d].io->averaging)
            {
                timeStart    = domain[d].io->avgStartTime;
                timeInterval = domain[d].io->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
            }

            // phase averaged fields
            if(domain[d].io->phaseAveraging)
            {
                timeStart    = domain[d].io->phAvgStartTime;
                timeInterval = domain[d].io->phAvgPrd;

                timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
            }

            // ke budgets
            if(domain[d].io->keBudgets)
            {
                timeStart    = acquisition->keBudFields->avgStartTime;
                timeInterval = acquisition->keBudFields->avgPrd;

                timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
            }

            if(flags->isAquisitionActive)
            {
                // 3LM averaging for background domain only
                if(acquisition->isAverage3LMActive)
                {
                    timeStart    = domain[0].acquisition->LM3->avgStartTime;
                    timeInterval = domain[0].acquisition->LM3->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
                }

                // ABL averaging for background domain only
                if(acquisition->isAverageABLActive)
                {
                    timeStart    = domain[0].acquisition->statisticsABL->avgStartTime;
                    timeInterval = domain[0].acquisition->statisticsABL->avgPrd;

                    timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
                }

                if(acquisition->isProbesActive)
                {
                    for(PetscInt r=0; r<acquisition->probes->nRakes; r++)
                    {
                        timeStart    = acquisition->probes->rakes[r].timeStart;
                        timeInterval = acquisition->probes->rakes[r].timeInterval;

                        timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
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

                            timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
                        }
                    }

                    // iSections
                    if(acquisition->jSections->available)
                    {
                        if(acquisition->jSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->jSections->timeStart;
                            timeInterval = acquisition->jSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
                        }
                    }

                    // kSections
                    if(acquisition->kSections->available)
                    {
                        if(acquisition->kSections->intervalType == "adjustableTime")
                        {
                            timeStart    = acquisition->kSections->timeStart;
                            timeInterval = acquisition->kSections->timeInterval;

                            timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
                        }
                    }
                }

                if(flags->isWindFarmActive)
                {
                    if(domain[d].farm->intervalType == "adjustableTime")
                    {
                        timeStart    = domain[d].farm->timeStart;
                        timeInterval = domain[d].farm->timeInterval;

                        timeStepSet(clock, timeStart, timeInterval, dx_min, maxU, flag, cfl);
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

                    PetscReal dtFarm = dx_min / farm->maxTipSpeed;

                    clock->dt = std::min(clock->dt, dtFarm);
                }
            }

        }
        else
        {
            if(clock->it>clock->itStart) clock->cfl =  maxU * clock->dt / dx_min;
        }
    }

    // prevent time step from increasing too fast
    if(clock->dt > 1.5 * clock->dtOld)
    {
        clock->dt = 1.5 * clock->dtOld;
    }

    // make sure to hit last time
    if(clock->time + clock->dt > clock->endTime)
    {
        clock->dt = clock->endTime - clock->time;
    }

    clock->time = clock->time + clock->dt;

    cfl = clock->dt * maxU / dx_min;

    PetscPrintf(PETSC_COMM_WORLD, "\n\nTime: %lf\n\n", clock->time);
    PetscPrintf(PETSC_COMM_WORLD, "Iteration = %ld, CFL = %lf, uMax = %.6f, dt = %.6f, deltaMin = %.6f, adjust due to write flag: %ld\n", clock->it, cfl, maxU, clock->dt, dx_min, flag);

    return(0);
}
