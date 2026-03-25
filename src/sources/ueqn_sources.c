//! \file  ueqn_sources.c
//! \brief Momentum-equation source-term and auxiliary-flux function implementations.

#include "ueqn_sources.h"

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

    PetscReal     lrelError = 0.0, grelError = 0.0;
    PetscInt      lcells = 0, gcells = 0;
    std::vector<Cmpnts> uMeso(my-2);


    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->sourceU, &source);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    if(abl->controllerAction == "write")
    {
        if(abl->controllerType=="pressure")
        {
            PetscReal     relax = abl->relax;
            PetscReal     alpha = abl->alpha;
            PetscReal     T     = abl->timeWindow;

            // set the wanted velocity
            if(abl->mesoScaleInputActive)
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

                // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
                if(abl->timeV[idx_2] < clock->time)
                {
                    idx_1 = idx_2;
                    idx_2 = idx_1 + 1;
                }

                if(abl->timeV[idx_1] > clock->time)
                {
                    idx_2 = idx_1;
                    idx_1 = idx_2 - 1;
                }

                // find interpolation weights
                PetscReal idx = (idx_2 - idx_1) / (abl->timeV[idx_2] - abl->timeV[idx_1]) * (clock->time - abl->timeV[idx_1]) + idx_1;
                PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
                PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

                PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading mesoscale data\n", w1 * abl->timeV[idx_1] + w2 * abl->timeV[idx_2]);
                PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
                PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeV[idx_1], abl->timeV[idx_2]);

                // reset the closest index for next iteration
                abl->closestTimeIndV = idx_1;

                PetscInt idxh1 = abl->velInterpIdx[0][0], idxh2 = abl->velInterpIdx[0][1];
                PetscReal wth1 = abl->velInterpWts[0][0], wth2  = abl->velInterpWts[0][1];

                //interpolating in time
                Cmpnts uH1 = nSum(nScale(w1, abl->uMeso[idxh1][idx_1]), nScale(w2, abl->uMeso[idxh1][idx_2]));
                Cmpnts uH2 = nSum(nScale(w1, abl->uMeso[idxh2][idx_1]), nScale(w2, abl->uMeso[idxh2][idx_2]));

                //interpolating in vertical direction 
                uDes = nSum(nScale(wth1, uH1), nScale(wth2, uH2));
            }
            else 
            {
                uDes.x = abl->uRef.x;
                uDes.y = abl->uRef.y;
                uDes.z = 0.0;
            }

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
                Cmpnts    Sc      = nSetFromComponents(s.x*clock->dt, s.y*clock->dt, 0.0);
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
        else if(abl->controllerType=="geostrophic")
        {
            PetscReal     relax = abl->relax;

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

            for (k=lzs; k<lze; k++)
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
            if(print) PetscPrintf(mesh->MESH_COMM, "                         U at hRef: (%.3f %.3f %.3f) m/s, U at hGeo: (%.3f %.3f %.3f)  = %0.3f m/s\n", uMean.x, uMean.y, uMean.z, uMeanGeo.x, uMeanGeo.y, uMeanGeo.z, nMag(uMeanGeo));

            // compute actual angles
            PetscReal hubAngle    = std::atan2(uMean.y, uMean.x);
            PetscReal geoAngleNew = std::atan2(uMeanGeo.y, uMeanGeo.x);
            
            if(abl->windAngleController)
            {
                // compute the uniform source terms for geostrophic wind controller
                abl->a.x = -2.0 * abl->fc * nMag(abl->uGeoBar) * sin(geoAngleNew);
                abl->a.y = 2.0 * abl->fc * nMag(abl->uGeoBar) * cos(geoAngleNew);
                abl->a.z = 0.0;

                // compute rotation speed at hub height (omega)
                PetscReal omega = relax * angleDiff(hubAngle, abl->refHubAngle);

                abl->b.x =   omega;
                abl->b.y = - omega;
                abl->b.z =   0.0;
            }
            else 
            {
                // compute the uniform source terms for geostrophic wind controller
                abl->a.x = -2.0 * abl->fc * abl->uGeoBar.y;
                abl->a.y = 2.0 * abl->fc * abl->uGeoBar.x;
                abl->a.z = 0.0;

                abl->b = nSetZero();
            }

            PetscPrintf(mesh->MESH_COMM, "                         Alpha at hRef: %.5f deg, Alpha at hGeo: %.3f deg\n", hubAngle/M_PI*180, geoAngleNew/M_PI*180);

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

                    PetscFPrintf(mesh->MESH_COMM, fp, "%*.6f\t%*.5e\t%*.5e%*.5e\t\n", width, clock->time, width, abl->a.x, width, abl->a.y, width, abl->a.z);

                    fclose(fp);
                }
            }
        }
        else if (abl->controllerType=="geostrophicProfileAssimilation")
        {
            PetscInt    nlevels = my-2;
            Cmpnts      uH1, uH2;
            PetscInt    idxh1, idxh2, idxt1;
            PetscReal   wth1, wth2, wtt1;
            PetscReal   alpha = 1.5; 

            Cmpnts    *srcPA  = abl->srcPA;
            Cmpnts    *luMean = abl->luMean;
            Cmpnts    *guMean = abl->guMean;

            //reset all fields to 0
            for(j=0; j<nlevels; j++)
            {
                luMean[j] = nSetZero();
                guMean[j] = nSetZero();
                srcPA[j]    = nSetZero();
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

            for(j=0; j<nlevels; j++)
            {
                mScale(1.0/abl->totVolPerLevel[j], guMean[j]);
            }

            //for cells above the abl height apply geostrophic forcing based on the mesodata 
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

            // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
            if(abl->timeV[idx_2] < clock->time)
            {
                idx_1 = idx_2;
                idx_2 = idx_1 + 1;
            }

            if(abl->timeV[idx_1] > clock->time)
            {
                idx_2 = idx_1;
                idx_1 = idx_2 - 1;
            }

            // find interpolation weights
            PetscReal idx = (idx_2 - idx_1) / (abl->timeV[idx_2] - abl->timeV[idx_1]) * (clock->time - abl->timeV[idx_1]) + idx_1;
            PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
            PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

            PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading mesoscale data\n", w1 * abl->timeV[idx_1] + w2 * abl->timeV[idx_2]);
            PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
            PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeV[idx_1], abl->timeV[idx_2]);

            // reset the closest index for next iteration
            abl->closestTimeIndV = idx_1;

            // find the atmospheric boundary layer height 
            findABLHeight(abl);

            PetscReal ablHeight = abl->ablHt, scaleFactor;

            // Step 1: Add the coriolis force + unsteady term into the pressure gradient source term
            for(j=0; j<nlevels; j++)
            {                
                idxh1 = abl->velInterpIdx[j][0];
                idxh2 = abl->velInterpIdx[j][1];

                wth1  = abl->velInterpWts[j][0];
                wth2  = abl->velInterpWts[j][1];

                //interpolating in time
                uH1 = nSum(nScale(w1, abl->uMeso[idxh1][idx_1]), nScale(w2, abl->uMeso[idxh1][idx_2]));
                uH2 = nSum(nScale(w1, abl->uMeso[idxh2][idx_1]), nScale(w2, abl->uMeso[idxh2][idx_2]));

                //interpolating in vertical direction 
                uMeso[j] = nSum(nScale(wth1, uH1), nScale(wth2, uH2));

                if(clock->it == clock->itStart)
                {
                    srcPA[j].x  =  - 2.0*abl->fc*uMeso[j].y;
                    srcPA[j].y  =    2.0*abl->fc*uMeso[j].x;
                    srcPA[j].z  =    0.0;
                }
                else
                {
                    srcPA[j].x  =  - 2.0*abl->fc*uMeso[j].y + (uMeso[j].x - abl->uGeoPrev[j].x)/PetscMax(clock->dt, 1e-8);
                    srcPA[j].y  =    2.0*abl->fc*uMeso[j].x + (uMeso[j].y - abl->uGeoPrev[j].y)/PetscMax(clock->dt, 1e-8);
                    srcPA[j].z  =    0.0;
                }
            }

            //save the previous mesoscale velocities 
            for(j=0; j<nlevels; j++)
            {
                abl->uGeoPrev[j] = nSet(uMeso[j]);
            }

            //Step 2: use constant source term in the region where no mesoscale data is available 
            i = 0;    
            PetscInt lMesoIndv, hMesoIndV;
            while(abl->cellLevels[i] < abl->hV[0])
            {
                i++;
            }
            lMesoIndv = i;

            i = nlevels-1;    
            while(abl->cellLevels[i] > abl->hV[abl->numhV - 1])
            {
                i--;
            }
            hMesoIndV = i;
            
            for(j=0; j<nlevels; j++)
            {
                if(abl->cellLevels[j] < abl->hV[0])
                {
                    srcPA[j].x = srcPA[lMesoIndv].x;
                    srcPA[j].y = srcPA[lMesoIndv].y;
                    srcPA[j].z = srcPA[lMesoIndv].z;
                } 

                if(abl->cellLevels[j] > abl->hV[abl->numhV - 1])
                {
                    srcPA[j].x = srcPA[hMesoIndV].x;
                    srcPA[j].y = srcPA[hMesoIndV].y;
                    srcPA[j].z = srcPA[hMesoIndV].z;
                }
            }

            // Step 3: add the frictional drag term 
            PetscReal ablgeoDelta = 0.2 * ablHeight;
            PetscReal ablgeoDampTop = ablHeight + 0.5*ablgeoDelta;

            PetscInt surfaceCell = 0;
            for(l = 0; l < nlevels; l++)
            {
                if(abl->cellLevels[l] > 0)
                {
                    surfaceCell = l;
                    break;
                }
            }

            for (int j = 0; j < nlevels; ++j)
            {
                //apply drag term only if below the ABL height
                scaleFactor = scaleHyperTangTop(abl->cellLevels[j], ablgeoDampTop, ablgeoDelta);

                PetscReal dz;
                if(j < surfaceCell)
                {
                    srcPA[j].x += 0;
                    srcPA[j].y += 0;
                    srcPA[j].z += 0;
                }
                else if(j == surfaceCell)
                {
                    dz = abl->cellLevels[j+1] - abl->cellLevels[j];
                    srcPA[j].x += (1 - scaleFactor) * (abl->avgStress[j+1].x - abl->avgStress[j].x) / dz;
                    srcPA[j].y += (1 - scaleFactor) * (abl->avgStress[j+1].y - abl->avgStress[j].y) / dz;
                    srcPA[j].z =   0.0;
                }
                else if (j == nlevels - 1)
                {
                    dz = abl->cellLevels[j] - abl->cellLevels[j-1];
                    srcPA[j].x += (1 - scaleFactor) * (abl->avgStress[j].x - abl->avgStress[j-1].x) / dz;
                    srcPA[j].y += (1 - scaleFactor) * (abl->avgStress[j].y - abl->avgStress[j-1].y) / dz;
                    srcPA[j].z =   0.0;
                }
                else
                {
                    dz = abl->cellLevels[j+1] - abl->cellLevels[j-1];
                    srcPA[j].x += (1 - scaleFactor) * (abl->avgStress[j+1].x - abl->avgStress[j-1].x) / dz;
                    srcPA[j].y += (1 - scaleFactor) * (abl->avgStress[j+1].y - abl->avgStress[j-1].y) / dz;
                    srcPA[j].z =   0.0;
                }
            }

            // //now add the damping term 
            // PetscReal ablgeoDelta = 200;
            // PetscReal ablgeoDampTop = ablHeight + 0.5*ablgeoDelta;

            // for(j=0; j<nlevels; j++)
            // {                
            //     idxh1 = abl->velInterpIdx[j][0];
            //     idxh2 = abl->velInterpIdx[j][1];

            //     wth1  = abl->velInterpWts[j][0];
            //     wth2  = abl->velInterpWts[j][1];

            //     //interpolating in time
            //     uH1 = nSum(nScale(w1, abl->uMeso[idxh1][idx_1]), nScale(w2, abl->uMeso[idxh1][idx_2]));
            //     uH2 = nSum(nScale(w1, abl->uMeso[idxh2][idx_1]), nScale(w2, abl->uMeso[idxh2][idx_2]));

            //     //interpolating in vertical direction 
            //     uMeso[j] = nSum(nScale(wth1, uH1), nScale(wth2, uH2));

            //     //ensure no damping within the boundary layer
            //     scaleFactor = scaleHyperTangTop(abl->cellLevels[j], ablgeoDampTop, ablgeoDelta);

            //     srcPA[j].x +=  -scaleFactor * 2.0*(2.0*abl->fc)*alpha*(guMean[j].x - uMeso[j].x);
            //     srcPA[j].y +=  -scaleFactor * 2.0*(2.0*abl->fc)*alpha*(guMean[j].y - uMeso[j].y);
            //     srcPA[j].z =   0.0;
            // }

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
                        PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e  %*.5e  %*.5e\t", width, srcPA[j].x,  width, srcPA[j].y,  width, srcPA[j].z);
                    }

                    PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                    fclose(fp);                  
                
                }
            }
        }
        else if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") || (abl->controllerType=="waveletProfileAssimilation"))
        {
            PetscReal relax = abl->relax;
            PetscReal alpha = abl->alpha;
            PetscReal T     = abl->timeWindow;
            PetscInt  nlevels = my-2;

            Cmpnts    *luMean = abl->luMean;
            Cmpnts    *guMean = abl->guMean;
            Cmpnts    *srcPA  = abl->srcPA;

            Cmpnts    uH1, uH2;
            PetscInt  idxh1, idxh2, idxt1;
            PetscReal wth1, wth2, wtt1;

            PetscReal aN, m1, m2;

            PetscMalloc(sizeof(Cmpnts) * (nlevels), &(src));

            for(j=0; j<nlevels; j++)
            {
                luMean[j] = nSetZero();
                guMean[j] = nSetZero();
                src[j]    = nSetZero();
                srcPA[j]    = nSetZero();
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

            for(j=0; j<nlevels; j++)
            {
                mScale(1.0/abl->totVolPerLevel[j], guMean[j]);
            }

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

            // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
            if(abl->timeV[idx_2] < clock->time)
            {
                idx_1 = idx_2;
                idx_2 = idx_1 + 1;
            }

            if(abl->timeV[idx_1] > clock->time)
            {
                idx_2 = idx_1;
                idx_1 = idx_2 - 1;
            }

            // find interpolation weights
            PetscReal idx = (idx_2 - idx_1) / (abl->timeV[idx_2] - abl->timeV[idx_1]) * (clock->time - abl->timeV[idx_1]) + idx_1;
            PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
            PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

            PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading mesoscale data\n", w1 * abl->timeV[idx_1] + w2 * abl->timeV[idx_2]);
            PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
            PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeV[idx_1], abl->timeV[idx_2]);

            // reset the closest index for next iteration
            abl->closestTimeIndV = idx_1;

            for(j=0; j<nlevels; j++)
            {                
                idxh1 = abl->velInterpIdx[j][0];
                idxh2 = abl->velInterpIdx[j][1];

                wth1  = abl->velInterpWts[j][0];
                wth2  = abl->velInterpWts[j][1];

                //interpolating in time
                uH1 = nSum(nScale(w1, abl->uMeso[idxh1][idx_1]), nScale(w2, abl->uMeso[idxh1][idx_2]));
                uH2 = nSum(nScale(w1, abl->uMeso[idxh2][idx_1]), nScale(w2, abl->uMeso[idxh2][idx_2]));

                //interpolating in vertical direction 
                uMeso[j] = nSum(nScale(wth1, uH1), nScale(wth2, uH2));

                // PetscPrintf(PETSC_COMM_WORLD, "uDes = %lf %lf %lf, guMean = %lf %lf %lf, interpolated from %lf and %lf, wts = %lf %lf\n", uMeso[j].x, uMeso[j].y, uMeso[j].z, guMean[j].x, guMean[j].y, guMean[j].z, abl->timeV[idx_1], abl->timeV[idx_2], w1, w2);

                // cumulate the error (integral part of the controller)
                abl->cumulatedSourceHt[j].x = (1.0 - clock->dt / T) * abl->cumulatedSourceHt[j].x + (clock->dt / T) * (uMeso[j].x - guMean[j].x);
                abl->cumulatedSourceHt[j].y = (1.0 - clock->dt / T) * abl->cumulatedSourceHt[j].y + (clock->dt / T) * (uMeso[j].y - guMean[j].y);
                abl->cumulatedSourceHt[j].z = 0.0;
                
                src[j].x = relax * (alpha * (uMeso[j].x - guMean[j].x) + (1.0-alpha) * abl->cumulatedSourceHt[j].x);
                src[j].y = relax * (alpha * (uMeso[j].y - guMean[j].y) + (1.0-alpha) * abl->cumulatedSourceHt[j].y) ;
                src[j].z = 0.0;
                
            }

            // set the source terms outside the provided range. this is set incorrectly above so corrected here
            if(abl->flType == "ablHeight")
            {
                findABLHeight(abl);
                PetscReal ablHeight = abl->ablHt;

                if(ablHeight > abl->hV[abl->numhV - 1])
                {
                    char error[512];
                    sprintf(error, "the abl height is greater than available mesoscale data height %lf. Use different lower layer\n", abl->hV[abl->numhV - 1]);
                    fatalErrorInFunction("CorrectSourceTerm",  error); 
                }

                //if abl height is below the lowest mesoscale data point
                if(ablHeight < abl->hV[0])
                {   
                    i = 0;    
                    while(abl->cellLevels[i] < abl->hV[0])
                    {
                        i++;
                    }

                    abl->lowestIndV = i;
                    abl->bottomSrcHtV = abl->hV[0];
                }
                // if abl height is in between the mesoscale data points
                else if(ablHeight >= abl->hV[0])
                {
                    i = 0;    

                    while(abl->cellLevels[i] < ablHeight)
                    {
                        i++;
                    }   

                    abl->lowestIndV = i;
                    abl->bottomSrcHtV = ablHeight;
                }
            }

            for(j=0; j<nlevels; j++)
            {
                if(abl->cellLevels[j] < abl->bottomSrcHtV)
                {
                    src[j].x = src[abl->lowestIndV].x;
                    src[j].y = src[abl->lowestIndV].y;
                    src[j].z = src[abl->lowestIndV].z;
                }
                
                if(abl->cellLevels[j] > abl->hV[abl->numhV - 1])
                {
                    src[j].x = src[abl->highestIndV].x;
                    src[j].y = src[abl->highestIndV].y;
                    src[j].z = src[abl->highestIndV].z;
                }

                srcPA[j] = nSet(src[j]);

                if(abl->averageSource && abl->controllerType=="directProfileAssimilation")
                {
                    if((abl->currAvgtime < abl->tAvgWindow) && (j == 0))
                    {
                        abl->currAvgtime = abl->currAvgtime + clock->dt;
                    }

                    aN = (PetscReal)abl->currAvgtime;
                    m1 = aN  / (aN + clock->dt);
                    m2 = clock->dt / (aN + clock->dt);

                    abl->avgsrc[j].x = m1*abl->avgsrc[j].x + m2*srcPA[j].x;
                    abl->avgsrc[j].y = m1*abl->avgsrc[j].y + m2*srcPA[j].y;
                    abl->avgsrc[j].z = m1*abl->avgsrc[j].z + m2*srcPA[j].z;
                }
            }

            if (abl->controllerType=="indirectProfileAssimilation")
            {
                for (j=0; j<nlevels; j++) 
                {
                    PetscReal dotProductx = 0.0, dotProducty = 0.0, dotProductz = 0.0;

                    for (k = 0; k < nlevels; k++) 
                    {
                        dotProductx += abl->polyCoeffM[j][k] * src[k].x;
                        dotProducty += abl->polyCoeffM[j][k] * src[k].y;
                        dotProductz += abl->polyCoeffM[j][k] * src[k].z;
                    }

                    srcPA[j].x = dotProductx;
                    srcPA[j].y = dotProducty;
                    srcPA[j].z = dotProductz;
                }

                for (j=0; j<nlevels; j++) 
                {            
                    if(abl->flType == "constantHeight" || abl->flType == "ablHeight")
                    {
                        if(abl->cellLevels[j] < abl->bottomSrcHtV)
                        {
                            srcPA[j].x = srcPA[abl->lowestIndV].x;
                            srcPA[j].y = srcPA[abl->lowestIndV].y;
                            srcPA[j].z = srcPA[abl->lowestIndV].z;
                        }
                    }

                    if(abl->cellLevels[j] > abl->hV[abl->numhV - 1])
                    {
                        srcPA[j].x = srcPA[abl->highestIndV].x;
                        srcPA[j].y = srcPA[abl->highestIndV].y;
                        srcPA[j].z = srcPA[abl->highestIndV].z;
                    }

                    if(abl->averageSource)
                    {
                        if((abl->currAvgtime < abl->tAvgWindow) && (j == 0))
                        {
                            abl->currAvgtime = abl->currAvgtime + clock->dt;
                        }

                        aN = (PetscReal)abl->currAvgtime;
                        m1 = aN  / (aN + clock->dt);
                        m2 = clock->dt / (aN + clock->dt);

                        abl->avgsrc[j].x = m1*abl->avgsrc[j].x + m2*srcPA[j].x;
                        abl->avgsrc[j].y = m1*abl->avgsrc[j].y + m2*srcPA[j].y;
                        abl->avgsrc[j].z = m1*abl->avgsrc[j].z + m2*srcPA[j].z;
                        
                    }

                    // PetscPrintf(PETSC_COMM_WORLD, "src[%ld] = %lf %lf %lf, avgsrc[%ld] = %lf %lf %lf\n", j, srcPA[j].x, srcPA[j].y, srcPA[j].z);
                }
            }

            if (abl->controllerType=="waveletProfileAssimilation")
            {                
                if(abl->waveletBlend)
                {
                    Cmpnts *srcAbove, *srcBelow;                                                                                                                
                    PetscMalloc(sizeof(Cmpnts) * (nlevels), &(srcAbove));
                    PetscMalloc(sizeof(Cmpnts) * (nlevels), &(srcBelow));

                    #if USE_PYTHON
                        pywavedecVector(abl, src, srcBelow, nlevels);
                    #endif
    
                    PetscInt waveLevel = abl->waveLevel;
                    
                    //hardcoded for now
                    abl->waveLevel = abl->waveLevel - 4;
                    if(abl->waveLevel <= 0) abl->waveLevel = 1;
    
                    #if USE_PYTHON
                        pywavedecVector(abl, src, srcAbove, nlevels);
                    #endif
    
                    abl->waveLevel = waveLevel;

                    // blend between the below and above profiles 
                    PetscReal blendHeight = mesh->bounds.zmax-500.0, scaleFactor;
                    PetscReal ablgeoDelta = 300;
                    PetscReal blendTop = blendHeight + 0.5*ablgeoDelta;
    
                    for(j=0; j<nlevels; j++)
                    {
                        scaleFactor = scaleHyperTangTop(abl->cellLevels[j], blendTop, ablgeoDelta);
    
                        srcPA[j].x = (1-scaleFactor)*srcBelow[j].x + scaleFactor*srcAbove[j].x;
                        srcPA[j].y = (1-scaleFactor)*srcBelow[j].y + scaleFactor*srcAbove[j].y;
                        srcPA[j].z = (1-scaleFactor)*srcBelow[j].z + scaleFactor*srcAbove[j].z;
                        // PetscPrintf(PETSC_COMM_WORLD, "%lf %lf %lf %lf\n", src[j].x, srcBelow[j].x, srcAbove[j].x, srcPA[j].x);
    
                    }
    
                    PetscFree(srcAbove);
                    PetscFree(srcBelow);

                }
                else
                {
                    #if USE_PYTHON
                    pywavedecVector(abl, src, srcPA, nlevels);
                    #endif 
                } 

                
                for (j=0; j<nlevels; j++) 
                {             
                    // if(abl->flType == "constantHeight" || abl->flType == "ablHeight")
                    // {
                        if(abl->cellLevels[j] < abl->bottomSrcHtV)
                        {
                            srcPA[j].x = srcPA[abl->lowestIndV].x;
                            srcPA[j].y = srcPA[abl->lowestIndV].y;
                            srcPA[j].z = srcPA[abl->lowestIndV].z;
                        }
                    
                    // }

                    if(abl->cellLevels[j] > abl->hV[abl->numhV - 1])
                    {
                        srcPA[j].x = srcPA[abl->highestIndV].x;
                        srcPA[j].y = srcPA[abl->highestIndV].y;
                        srcPA[j].z = srcPA[abl->highestIndV].z;
                    } 

                    if(abl->averageSource)
                    {
                        if((abl->currAvgtime < abl->tAvgWindow) && (j == 0))
                        {
                            abl->currAvgtime = abl->currAvgtime + clock->dt;
                        }

                        aN = (PetscReal)abl->currAvgtime;
                        m1 = aN  / (aN + clock->dt);
                        m2 = clock->dt / (aN + clock->dt);

                        abl->avgsrc[j].x = m1*abl->avgsrc[j].x + m2*srcPA[j].x;
                        abl->avgsrc[j].y = m1*abl->avgsrc[j].y + m2*srcPA[j].y;
                        abl->avgsrc[j].z = m1*abl->avgsrc[j].z + m2*srcPA[j].z;
                        
                    }

                    // PetscPrintf(PETSC_COMM_WORLD, "src[%ld] = %lf %lf %lf, avgsrc[%ld] = %lf %lf %lf\n", j, srcPA[j].x, srcPA[j].y, srcPA[j].z);
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

                    if(abl->averageSource)
                    {
                        for(j=0; j<nlevels; j++)
                        {
                            PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e  %*.5e  %*.5e\t", width, abl->avgsrc[j].x,  width, abl->avgsrc[j].y,  width, abl->avgsrc[j].z);
                        }                       
                    }
                    else 
                    {
                        for(j=0; j<nlevels; j++)
                        {
                            PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e  %*.5e  %*.5e\t", width, srcPA[j].x,  width, srcPA[j].y,  width, srcPA[j].z);
                        }
                    }

                    PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                    fclose(fp);                  
                
                }
            }
        }
    }
    else if(abl->controllerAction == "read")
    {

        if(abl->controllerType=="timeSeries")
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
            s.x = w1 * abl->preCompSources[idx_1][1] + w2 * abl->preCompSources[idx_2][1];
            s.y = w1 * abl->preCompSources[idx_1][2] + w2 * abl->preCompSources[idx_2][2];
            s.z = w1 * abl->preCompSources[idx_1][3] + w2 * abl->preCompSources[idx_2][3];

            // scale with time step
            // s.x = (w1 * abl->preCompSources[idx_1][1] / dtSource1 + w2 * abl->preCompSources[idx_2][1] / dtSource2) * clock->dt;
            // s.y = (w1 * abl->preCompSources[idx_1][2] / dtSource1 + w2 * abl->preCompSources[idx_2][2] / dtSource2) * clock->dt;
            // s.z = (w1 * abl->preCompSources[idx_1][3] / dtSource1 + w2 * abl->preCompSources[idx_2][3] / dtSource2) * clock->dt;


        }
        else if(abl->controllerType=="timeAverageSeries")
        {
            // PetscPrintf(mesh->MESH_COMM, "Correcting source terms using source average from time %lf\n", abl->sourceAvgStartTime);

            // do not scale with dt
            s.x = abl->cumulatedSource.x;
            s.y = abl->cumulatedSource.y;
            s.z = abl->cumulatedSource.z;

            // scale with dt
            //PetscReal dtScale = clock->dt / abl->avgTimeStep;

            // s.x = abl->cumulatedSource.x * dtScale;
            // s.y = abl->cumulatedSource.y * dtScale;
            // s.z = abl->cumulatedSource.z * dtScale;
        }
        else if(abl->controllerType=="timeHeightSeries")
        {
            PetscInt   nlevels = my-2;
            PetscReal  uT1, uT2;
            PetscReal  wth1, wth2;
            PetscInt   idxh1, idxh2;

            PetscMalloc(sizeof(Cmpnts) * (nlevels), &(src));

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
                diff[i] = fabs(abl->timeHtSources[i][0][0] - clock->time);
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
            if(abl->timeHtSources[idx_2][0][0] < abl->timeHtSources[idx_1][0][0])
            {
                PetscInt idx_tmp = idx_2;
                idx_2 = idx_1;
                idx_1 = idx_tmp;
            }

            // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
            if(abl->timeHtSources[idx_2][0][0] < clock->time)
            {
                idx_1 = idx_2;
                idx_2 = idx_1 + 1;
            }

            if(abl->timeHtSources[idx_1][0][0] > clock->time)
            {
                idx_2 = idx_1;
                idx_1 = idx_2 - 1;
            }
            
            // find interpolation weights
            PetscReal idx = (idx_2 - idx_1) / (abl->timeHtSources[idx_2][0][0] - abl->timeHtSources[idx_1][0][0]) * (clock->time - abl->timeHtSources[idx_1][0][0]) + idx_1;
            PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
            PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);
            
            if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: selected time %lf for reading sources\n", w1 * abl->timeHtSources[idx_1][0][0] + w2 * abl->timeHtSources[idx_2][0][0]);
            if(print) PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
            if(print) PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeHtSources[idx_1][0][0], abl->timeHtSources[idx_2][0][0]);

            // reset the closest index for nex iteration
            abl->currentCloseIdx = idx_1;

            // get also the dt at the time the source was calculated
            double dtSource1 = std::max((abl->timeHtSources[idx_1+1][0][0] - abl->timeHtSources[idx_1][0][0]), 1e-5);
            double dtSource2 = std::max((abl->timeHtSources[idx_2+1][0][0] - abl->timeHtSources[idx_2][0][0]), 1e-5);

            // scale with time step for each height
            for(j=0; j<nlevels; j++)
            {
                if(abl->numhV != nlevels)
                {
                    //interpolating in height for each time
                    idxh1 = abl->velInterpIdx[j][0];
                    idxh2 = abl->velInterpIdx[j][1];

                    wth1  = abl->velInterpWts[j][0];
                    wth2  = abl->velInterpWts[j][1];
                    
                    uT1   = wth1 * abl->timeHtSources[idx_1][idxh1][1] + wth2 * abl->timeHtSources[idx_1][idxh2][1];
                    uT2   = wth1 * abl->timeHtSources[idx_2][idxh1][1] + wth2 * abl->timeHtSources[idx_2][idxh2][1];

                    // src[j].x = (w1 * uT1 / dtSource1 + w2 * uT2 / dtSource2) * clock->dt;

                    src[j].x = w1 * uT1 + w2 * uT2;

                    uT1   = wth1 * abl->timeHtSources[idx_1][idxh1][2] + wth2 * abl->timeHtSources[idx_1][idxh2][2];
                    uT2   = wth1 * abl->timeHtSources[idx_2][idxh1][2] + wth2 * abl->timeHtSources[idx_2][idxh2][2];

                    // src[j].y = (w1 * uT1 / dtSource1 + w2 * uT2 / dtSource2) * clock->dt;
                    src[j].y = w1 * uT1 + w2 * uT2;
                    
                    uT1   = wth1 * abl->timeHtSources[idx_1][idxh1][3] + wth2 * abl->timeHtSources[idx_1][idxh2][3];
                    uT2   = wth1 * abl->timeHtSources[idx_2][idxh1][3] + wth2 * abl->timeHtSources[idx_2][idxh2][3];

                    // src[j].z = (w1 * uT1 / dtSource1 + w2 * uT2 / dtSource2) * clock->dt; 
                    src[j].z = w1 * uT1 + w2 * uT2;
                }
                else 
                {
                    // src[j].x = (w1 * abl->timeHtSources[idx_1][j][1] / dtSource1 + w2 * abl->timeHtSources[idx_2][j][1] / dtSource2) * clock->dt;
                    // src[j].y = (w1 * abl->timeHtSources[idx_1][j][2] / dtSource1 + w2 * abl->timeHtSources[idx_2][j][2] / dtSource2) * clock->dt;
                    // src[j].z = (w1 * abl->timeHtSources[idx_1][j][3] / dtSource1 + w2 * abl->timeHtSources[idx_2][j][3] / dtSource2) * clock->dt; 

                    src[j].x = (w1 * abl->timeHtSources[idx_1][j][1] + w2 * abl->timeHtSources[idx_2][j][1] ) ;
                    src[j].y = (w1 * abl->timeHtSources[idx_1][j][2] + w2 * abl->timeHtSources[idx_2][j][2] ) ;
                    src[j].z = (w1 * abl->timeHtSources[idx_1][j][3] + w2 * abl->timeHtSources[idx_2][j][3] ) ; 

                    // if(j%40 == 0)
                    // PetscPrintf(mesh->MESH_COMM, "src at level %d, height %lf is %lf %lf dt = %lf %lf %lf, orig src = %lf %lf\n", j, abl->cellLevels[j], src[j].x, w1 * abl->timeHtSources[idx_1][j][1] + w2 * abl->timeHtSources[idx_2][j][1], dtSource1, dtSource2, clock->dt, abl->timeHtSources[idx_1][j][1], abl->timeHtSources[idx_2][j][1]);
                    
                } 
            }           
 
        }
        else if(abl->controllerType=="timeSeriesFromPrecursor")
        {
            s.x = abl->preCompSources[0][1];
            s.y = abl->preCompSources[0][2];
            s.z = abl->preCompSources[0][3];
        }
    }

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // update the source term
    // source term made independent of simulation time step dt to prevent variation with time step changes 
    // previously source term = deltaU / dt, now source term = alpha deltaU, where alpha is a constant with units of 1/s
    // this prevents need for source term scaling in read mode based on the time step used when the source term was computed
    // source term is scaled when multiplied by dt in the momentum equation assembly instead
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(cent[k][j][i].z <= abl->controllerMaxHeight)
                {
                    if(abl->controllerAction == "write")
                    {
                        if(abl->controllerType=="geostrophic")
                        {
                            source[k][j][i].x = abl->a.x + abl->b.x * ucat[k][j][i].y;
                            
                            source[k][j][i].y = abl->a.y + abl->b.y * ucat[k][j][i].x;
                            
                            source[k][j][i].z = 0.0;
                        }
                        else if(abl->controllerType=="pressure")
                        {
                            if(applyGeoDamping)
                            {
                                PetscReal  height = cent[k][j][i].z - mesh->grndLevel;
                                source[k][j][i].x = (s.x +
                                                    scaleHyperTangTop   (height, abl->geoDampH, abl->geoDampDelta) *
                                                    abl->geoDampC*abl->geoDampAlpha*(abl->geoDampUBar.x - abl->geoDampU[j-1].x));
                                source[k][j][i].y = (s.y +
                                                    scaleHyperTangTop   (height, abl->geoDampH, abl->geoDampDelta) *
                                                    abl->geoDampC*abl->geoDampAlpha*(abl->geoDampUBar.y - abl->geoDampU[j-1].y));
                                source[k][j][i].z = 0.0;
                            }
                            else
                            {
                                source[k][j][i].x = s.x;
                                source[k][j][i].y = s.y;
                                source[k][j][i].z = 0.0;
                            }
                        }
                        else if(abl->controllerType=="geostrophicProfileAssimilation")
                        {
                            source[k][j][i].x = abl->srcPA[j-1].x;
                            source[k][j][i].y = abl->srcPA[j-1].y;
                            source[k][j][i].z = 0.0;
                        }
                        else if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") || (abl->controllerType=="waveletProfileAssimilation"))
                        {
                            if(abl->averageSource)
                            {
                                source[k][j][i].x = abl->avgsrc[j-1].x;
                                source[k][j][i].y = abl->avgsrc[j-1].y;
                                source[k][j][i].z = 0.0;  

                                if((abl->cellLevels[j-1] > abl->bottomSrcHtV) && (abl->cellLevels[j-1] < abl->hV[abl->numhV - 1]))
                                {
                                    lrelError += std::pow(abl->avgsrc[j-1].x*abl->avgsrc[j-1].x + abl->avgsrc[j-1].y*abl->avgsrc[j-1].y, 0.5)/std::pow(uMeso[j-1].x*uMeso[j-1].x + uMeso[j-1].y*uMeso[j-1].y, 0.5);
                                    lcells ++;
                                }                              
                            }
                            else
                            {
                                source[k][j][i].x = abl->srcPA[j-1].x;
                                source[k][j][i].y = abl->srcPA[j-1].y;
                                source[k][j][i].z = 0.0;
                                                                
                                if((abl->cellLevels[j-1] > abl->bottomSrcHtV) && (abl->cellLevels[j-1] < abl->hV[abl->numhV - 1]))
                                {
                                    lrelError += std::pow(abl->srcPA[j-1].x*abl->srcPA[j-1].x + abl->srcPA[j-1].y*abl->srcPA[j-1].y, 0.5)/std::pow(uMeso[j-1].x*uMeso[j-1].x + uMeso[j-1].y*uMeso[j-1].y, 0.5);
                                    lcells ++;
                                } 
                            }
                        } 
                    }
                    else if(abl->controllerAction == "read")
                    {
                        if(abl->controllerType == "timeSeries" || abl->controllerType == "timeAverageSeries" || abl->controllerType == "timeSeriesFromPrecursor")
                        {
                            source[k][j][i].x = s.x;
                            source[k][j][i].y = s.y;
                            source[k][j][i].z = 0.0;
                        }
                        else if(abl->controllerType == "timeHeightSeries")
                        {
                            source[k][j][i].x = src[j-1].x;
                            source[k][j][i].y = src[j-1].y;
                            source[k][j][i].z = 0.0;
                        }
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

    MPI_Allreduce(&lrelError, &grelError, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lcells, &gcells, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

    if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") || (abl->controllerType=="waveletProfileAssimilation") )
    {
        if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms: global avg error on velocity = %.5f percent\n", grelError*100/gcells);
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

    if(abl->controllerAction == "write")
    {
        if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") || (abl->controllerType=="waveletProfileAssimilation"))
        {
            free(src);
        }
    }
    else if(abl->controllerAction == "read")
    {
        if (abl->controllerType=="timeHeightSeries")
        {
            free(src);
        }       
    }

    std::vector<Cmpnts> ().swap(uMeso);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->sourceU, &source);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

	DMLocalToLocalBegin(fda, ueqn->sourceU, INSERT_VALUES, ueqn->sourceU);
    DMLocalToLocalEnd  (fda, ueqn->sourceU, INSERT_VALUES, ueqn->sourceU);

	resetCellPeriodicFluxes(mesh, ueqn->sourceU, ueqn->sourceU, "vector", "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode mapYDamping(ueqn_ *ueqn)
{
    mesh_         *mesh = ueqn->access->mesh;
    abl_          *abl  = ueqn->access->abl;
    flags_        *flags= ueqn->access->flags;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ubar;
    PetscReal     ***tbar;
    Cmpnts        **velMapped = abl->velMapped;
    PetscReal     **tMapped   = abl->tMapped;
    PetscInt      numJ, numK, numI;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, p, l;
    PetscInt      imin, imax, jmin, jmax, kmin, kmax;
    PetscInt      iminSrc, jminSrc, kminSrc;
    PetscInt      imaxSrc, jmaxSrc, kmaxSrc;

    PetscInt      numDestBounds = abl->yDampingNumPeriods;
    PetscInt      k_src_left, k_src_right;
    PetscReal     w_left, w_right;

    VecCopy(ueqn->lUcat, abl->uBarInstY);
    
    if(flags->isTeqnActive)
    {
        VecCopy(ueqn->access->teqn->lTmprt, abl->tBarInstY);
    }

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;
    
    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, abl->uBarInstY, &ubar);

    if(flags->isTeqnActive)
    {
        DMDAVecGetArray(da, abl->tBarInstY, &tbar);
    }

    if(flags->isYDampingActive)
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

                    imaxSrc = abl->srcMaxInd[p].i;
                    jmaxSrc = abl->srcMaxInd[p].j;
                    kmaxSrc = abl->srcMaxInd[p].k;

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

                                        //if the k closest index go out of processor bounds
                                        if(k_src_left  < kminSrc)  k_src_left  = kminSrc;
                                        if(k_src_right > kmaxSrc)  k_src_right = kmaxSrc;

                                        w_left        = abl->wtsKCell[k][0];
                                        w_right       = abl->wtsKCell[k][1];

                                        //Note: k index has been converted to indices of source processors, j and i are still based on destination processors
                                        Cmpnts uLeft  = nSet(velMapped[p][(k_src_left-kminSrc) * numJ * numI + (j-jmin) * numI + (i-imin)]);
                                        Cmpnts uRight = nSet(velMapped[p][(k_src_right-kminSrc) * numJ * numI + (j-jmin) * numI + (i-imin)]);

                                        ubar[k][j][i] = nSum(nScale(w_left, uLeft), nScale(w_right, uRight));
                                        
                                        if(flags->isTeqnActive)
                                        {
                                            PetscReal tLeft  = tMapped[p][(k_src_left-kminSrc) * numJ * numI + (j-jmin) * numI + (i-imin)];
                                            PetscReal tRight = tMapped[p][(k_src_right-kminSrc) * numJ * numI + (j-jmin) * numI + (i-imin)];

                                            tbar[k][j][i] = w_left * tLeft + w_right * tRight;                                                
                                        }
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

    if(flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, abl->tBarInstY, &tbar);
    }

    MPI_Barrier(mesh->MESH_COMM);

    DMLocalToLocalBegin (fda,  abl->uBarInstY, INSERT_VALUES, abl->uBarInstY);
    DMLocalToLocalEnd   (fda,  abl->uBarInstY, INSERT_VALUES, abl->uBarInstY);

    if(flags->isTeqnActive)
    {
        DMLocalToLocalBegin (da,  abl->tBarInstY, INSERT_VALUES, abl->tBarInstY);
        DMLocalToLocalEnd   (da,  abl->tBarInstY, INSERT_VALUES, abl->tBarInstY);
    }
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
                    // type 2: inflow and actual meshes are different, use avgTopLength
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
                        PetscReal h = cent[k_idx][j][i].z - mesh->grndLevel;
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
                            PetscReal height = cent[k_idx][j][i].z - mesh->grndLevel;
                            PetscInt  IDs[2];
                            PetscReal Wg [2];

                            findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

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
                        Cmpnts uGhost;

                        uGhost.x =
                            ifPtr->inflowWeights[j][i][0] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].x +
                            ifPtr->inflowWeights[j][i][1] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].x +
                            ifPtr->inflowWeights[j][i][2] *
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].x +
                            ifPtr->inflowWeights[j][i][3] ;
                            ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].x;
                        
                        if(ifPtr->interpMethod == "spline")
                        {
                            uGhost.y =
                                ifPtr->inflowWeights_2[j][i][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][0].j][ifPtr->closestCells_2[j][i][0].i].y +
                                ifPtr->inflowWeights_2[j][i][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][1].j][ifPtr->closestCells_2[j][i][1].i].y +
                                ifPtr->inflowWeights_2[j][i][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][2].j][ifPtr->closestCells_2[j][i][2].i].y +
                                ifPtr->inflowWeights_2[j][i][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][3].j][ifPtr->closestCells_2[j][i][3].i].y +
                                ifPtr->inflowWeights_2[j][i][4] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][4].j][ifPtr->closestCells_2[j][i][4].i].y +
                                ifPtr->inflowWeights_2[j][i][5] *
                                ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][5].j][ifPtr->closestCells_2[j][i][5].i].y;

                            uGhost.z =
                                ifPtr->inflowWeights_1[j][i][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][0].j][ifPtr->closestCells_1[j][i][0].i].z +
                                ifPtr->inflowWeights_1[j][i][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][1].j][ifPtr->closestCells_1[j][i][1].i].z +
                                ifPtr->inflowWeights_1[j][i][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][2].j][ifPtr->closestCells_1[j][i][2].i].z +
                                ifPtr->inflowWeights_1[j][i][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][3].j][ifPtr->closestCells_1[j][i][3].i].z +
                                ifPtr->inflowWeights_1[j][i][4] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][4].j][ifPtr->closestCells_1[j][i][4].i].z +
                                ifPtr->inflowWeights_1[j][i][5] *
                                ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][5].j][ifPtr->closestCells_1[j][i][5].i].z;
                        }
                        else
                        {
                            uGhost.y =
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].y +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].y +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].y +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].y;

                            uGhost.z =
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].z +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].z +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].z +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].z;
                        }

                        PetscReal height = cent[k_idx][j][i].z - mesh->grndLevel;
                        PetscInt  IDs[2];
                        PetscReal Wg [2];

                        findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                        luBarInstX[j][i].x
                        =
                        scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                        (
                            uGhost.x
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
                            uGhost.y
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
                            uGhost.z
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
                        PetscReal h      = cent[k_idx][j][i].z - mesh->grndLevel;

                        luBarInstX[j][i] = nSet(NieuwstadtInflowEvaluate(ifPtr, h));
                    }

                    if(ueqn->access->flags->isTeqnActive)
                    {
                        // steady prescribed ubar
                        if (ifPtr->typeT == 0 || ifPtr->typeT == 4)
                        {
                            PetscReal b      = abl->smear * abl->gTop * abl->dInv;
                            PetscReal a      = abl->gInv - b;
                            PetscReal h      = cent[k_idx][j][i].z - mesh->grndLevel;
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
                                PetscReal height = cent[k_idx][j][i].z - mesh->grndLevel;
                                PetscInt  IDs[2];
                                PetscReal Wg [2];

                                findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);


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
						    PetscReal tGhost;

                            tGhost =
    							ifPtr->inflowWeights[j][i][0] *
    							ifPtr->t_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i] +
    							ifPtr->inflowWeights[j][i][1] *
    							ifPtr->t_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i] +
    							ifPtr->inflowWeights[j][i][2] *
    							ifPtr->t_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i] +
    							ifPtr->inflowWeights[j][i][3] *
    							ifPtr->t_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i];

                            PetscReal height = cent[k_idx][j][i].z - mesh->grndLevel;

                            PetscInt  IDs[2];
                            PetscReal Wg [2];
                            findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                            ltBarInstX[j][i]
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                tGhost
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

                PetscReal uRefMag     = PetscSqrtReal(abl->uRef.x*abl->uRef.x + abl->uRef.y*abl->uRef.y);

                PetscReal fringeTime  = clock->time - abl->xDampingTimeStart;
                PetscReal waitTime    = abl->xDampingTimeWindow; // (abl->xDampingEnd - abl->xDampingStart)/uRefMag;


                // see if must correct alpha
                if
                (
                    (abl->xDampingError/uRefMag) > 0.01 && // check if error is greater than 1%
                    fringeTime > waitTime                    // check if allowed to override alpha
                )
                {
                    abl->xDampingAlpha = fabs(abl->xDampingVBar - abl->vStart) / abl->xDampingCoeff;
                    abl->xDampingTimeStart = clock->time;
                }

                // bound alpha
                abl->xDampingAlpha = std::max(std::min(abl->xDampingAlpha, 1.0), 0.0);

                // compute useful print quantities
                PetscReal percErrVStart  = fabs(abl->xDampingVBar - abl->vStart) / uRefMag * 100.0;
                PetscReal percErrVEnd    = abl->xDampingError / uRefMag * 100.0;
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

    // update uBar for yDamping layer (start mapping from the streamwise fringe region)
    if(ueqn->access->flags->isYDampingActive)
    {
        cellIds *srcMinInd = abl->srcMinInd, *srcMaxInd = abl->srcMaxInd;
        Cmpnts  **velMapped = abl->velMapped;
        PetscInt numJ, numK, numI;
        PetscInt imin, imax, jmin, jmax, kmin, kmax;
        precursor_    *precursor;
        domain_       *pdomain;

        // update the unsteady uBar state
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->fda, pdomain->ueqn->lUcat,  &ucatP);
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

                    // PetscPrintf(PETSC_COMM_SELF, "p = %ld,rank = %d, src kmin = %ld %ld, jmin = %ld %ld , imin = %ld %ld, num k,j, i = %ld %ld %ld\n", p,  rank, kmin, kmax, jmin, jmax, imin, imax, numK, numJ, numI);

                    //loop through the sub domain of this processor
                    for (k=kmin; k<=kmax; k++)
                    {
                        for (j=jmin; j<=jmax; j++)
                        {
                            for (i=imin; i<=imax; i++)
                            {
                                velMapped[p][(k-kmin) * numJ * numI + (j-jmin) * numI + (i-imin)] = nSet(ucatP[k][j][i]);
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

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->fda, pdomain->ueqn->lUcat,  &ucatP);
            }
        }
    }

    //update tBar for yDamping layer
    if(ueqn->access->flags->isTeqnActive)
    {
        correctDampingSourcesT(ueqn->access->teqn);
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
    PetscReal     ***nvert, ***tmprt, ***aj, ***meshTag;

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
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

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

                if(isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag))
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

                if(isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag))
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

                if(isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag))
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
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->gCont,   &gcont);
    DMDAVecRestoreArray(fda, mesh->lCent,   &cent);

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
    PetscReal     ***nvert, ***meshTag;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
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
                    isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag)
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
                        isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag)
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
                        isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag)
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
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

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
    PetscInt      kStart, kEnd;
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
    
    // ascend function factor
    PetscReal fAsc_x, fAsc_xi, fAsc_xj, fAsc_xk;

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

                    x     = cent[k][j][i].x;
                    xi    = cent[k][j][i+1].x;
                    xj    = cent[k][j+1][i].x;
                    xk    = cent[k+1][j][i].x;

                    // compute Nordstrom viscosity at i,j,k, i+1,j,k, i,j+1,k and i,j,k+1 points
                    nud_y   = viscNordstrom(alphaY, yS, yE, yD, y);
                    nudi_y  = viscNordstrom(alphaY, yS, yE, yD, yi);
                    nudj_y  = viscNordstrom(alphaY, yS, yE, yD, yj);
                    nudk_y  = viscNordstrom(alphaY, yS, yE, yD, yk);

                    //scaling factor to ensure that there is no lateral fringe region in the streamfringe region domain.
                    fAsc_x  = viscNordstromNoVertFilter(xS, xE, xD, x);
                    fAsc_xi = viscNordstromNoVertFilter(xS, xE, xD, xi);
                    fAsc_xj = viscNordstromNoVertFilter(xS, xE, xD, xj);
                    fAsc_xk = viscNordstromNoVertFilter(xS, xE, xD, xk);

                    if(isInflowSpreadType)
                    {
                        //lateral damping not defined
                    }
                    else
                    {
                        PetscReal uBarContK;
                        PetscReal uBarContJ;
                        PetscReal uBarContI;

                        uBarContI = central(uBarY[k][j][i].x, uBarY[k][j][i+1].x) * icsi[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k][j][i+1].y) * icsi[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k][j][i+1].z) * icsi[k][j][i].z;
                        uBarContJ = central(uBarY[k][j][i].x, uBarY[k][j+1][i].x) * jeta[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k][j+1][i].y) * jeta[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k][j+1][i].z) * jeta[k][j][i].z;
                        uBarContK = central(uBarY[k][j][i].x, uBarY[k+1][j][i].x) * kzet[k][j][i].x + central(uBarY[k][j][i].y, uBarY[k+1][j][i].y) * kzet[k][j][i].y + central(uBarY[k][j][i].z, uBarY[k+1][j][i].z) * kzet[k][j][i].z;

                        // i-fluxes
                        rhs[k][j][i].x
                        +=
                        scale * central(nud_y, nudi_y) * central(fAsc_x, fAsc_xi) * 
                        (
                            uBarContI - ucont[k][j][i].x
                        );

                        // j-fluxes
                        rhs[k][j][i].y
                        +=
                        scale * central(nud_y, nudj_y) * central(fAsc_x, fAsc_xj) * 
                        (
                            uBarContJ - ucont[k][j][i].y
                        );

                        // k-fluxes
                        rhs[k][j][i].z
                        +=
                        scale * central(nud_y, nudk_y) * central(fAsc_x, fAsc_xk) *
                        (
                            uBarContK - ucont[k][j][i].z
                        );
                    }
                }

                // z damping layer
                if(ueqn->access->flags->isZDampingActive)
                {
                    // compute cell center z at i,j,k and i,j+1,k points
                    z     = (cent[k][j][i].z   - mesh->grndLevel);
                    zj    = (cent[k][j+1][i].z - mesh->grndLevel);

                    // compute Rayleigh viscosity at i,j,k and i,j+1,k points
                    nud_z   = viscRayleigh(alphaZ, zS, zE, z);
                    nudj_z  = viscRayleigh(alphaZ, zS, zE, zj);

                    // damp also x and y components (exclude in xFringe if present)
                    if(abl->zDampingAlsoXY)
                    {
                        // compute cell center z at i+1,j,k and i,j,k+1 points
                        zi    = (cent[k][j][i+1].z - mesh->grndLevel);
                        zk    = (cent[k+1][j][i].z - mesh->grndLevel);

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
                    z     = (cent[k][j][i].z   - mesh->grndLevel);
                    zi    = (cent[k][j][i+1].z - mesh->grndLevel);
                    zj    = (cent[k][j+1][i].z - mesh->grndLevel);
                    zk    = (cent[k+1][j][i].z - mesh->grndLevel);

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
                    z     = (cent[k][j][i].z   - mesh->grndLevel);
                    zi    = (cent[k][j][i+1].z - mesh->grndLevel);
                    zj    = (cent[k][j+1][i].z - mesh->grndLevel);
                    zk    = (cent[k+1][j][i].z - mesh->grndLevel);

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
    PetscReal     ***nvert, ***meshTag;

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
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(fda, Rhs,  &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag))
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

                if(isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag))
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

                if(isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag))
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
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
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
    PetscReal     ***nvert, ***meshTag;

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
    PetscReal     cft     = ueqn->access->abl->cftCanopy;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(fda, Rhs,  &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal uSource;
                PetscReal coeff_i = 0.0, coeff_j = 0.0, coeff_k = 0.0;

                // coords of cell center and adjacent face centers
                PetscReal x  = cent[k][j][i].x, y  = cent[k][j][i].y, z  = cent[k][j][i].z;
                PetscReal yi = central(cent[k][j][i].y, cent[k][j][i+1].y);
                PetscReal zj = central(cent[k][j][i].z, cent[k][j][i+1].z);
                PetscReal xk = central(cent[k][j][i].x, cent[k+1][j][i].x);

                // i-face center in canopy box?
                if(x  < xEnd && x  > xStart && yi < yEnd && yi > yStart && z  < zEnd && z  > zStart) coeff_i = 1;
                // j-face center in canopy box?
                if(x  < xEnd && x  > xStart && y  < yEnd && y  > yStart && zj < zEnd && zj > zStart) coeff_j = 1;
                // k-face center in canopy box?
                if(xk < xEnd && xk > xStart && y  < yEnd && y  > yStart && z  < zEnd && z  > zStart) coeff_k = 1;

                if(coeff_i > 0 && isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag))
                {
                    uSource = nMag(centralVec(ucat[k][j][i], ucat[k][j][i+1]));
                    Cmpnts forceI = nScale(0.5*cft*uSource*uSource/Hc, diskNormal);
                    rhs[k][j][i].x += coeff_i * (forceI.x*icsi[k][j][i].x + forceI.y*icsi[k][j][i].y + forceI.z*icsi[k][j][i].z);
                }

                if(coeff_j > 0 && isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag))
                {
                    uSource = nMag(centralVec(ucat[k][j][i], ucat[k][j+1][i]));
                    Cmpnts forceJ = nScale(0.5*cft*uSource*uSource/Hc, diskNormal);
                    rhs[k][j][i].y += coeff_j * (forceJ.x*jeta[k][j][i].x + forceJ.y*jeta[k][j][i].y + forceJ.z*jeta[k][j][i].z);
                }

                if(coeff_k > 0 && isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag))
                {
                    uSource = nMag(centralVec(ucat[k][j][i], ucat[k+1][j][i]));
                    Cmpnts forceK = nScale(0.5*cft*uSource*uSource/Hc, diskNormal);
                    rhs[k][j][i].z += coeff_k * (forceK.x*kzet[k][j][i].x + forceK.y*kzet[k][j][i].y + forceK.z*kzet[k][j][i].z);
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    DMDAVecRestoreArray(fda, Rhs,  &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode meanGradPForcing(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM             da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo  info = mesh->info;
    PetscInt       xs   = info.xs, xe = info.xs + info.xm;
    PetscInt       ys   = info.ys, ye = info.ys + info.ym;
    PetscInt       zs   = info.zs, ze = info.zs + info.zm;
    PetscInt       mx   = info.mx, my = info.my, mz = info.mz;
    
    Cmpnts        ***rhs;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert, ***meshTag, cellVol;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;
    
    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);

    DMDAVecGetArray(fda, Rhs,            &rhs);
    DMDAVecGetArray(da,  mesh->lNvert,   &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);

    Cmpnts meanGradP = ueqn->meanGradP;

    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                if(isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag))
                {
                    rhs[k][j][i].x 
                    +=
                    scale *
                    (
                        meanGradP.x * icsi[k][j][i].x +
                        meanGradP.y * icsi[k][j][i].y +
                        meanGradP.z * icsi[k][j][i].z
                    );
                }

                if(isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag))
                {

                    rhs[k][j][i].y
                    +=
                    scale *
                    (
                        meanGradP.x * jeta[k][j][i].x +
                        meanGradP.y * jeta[k][j][i].y +
                        meanGradP.z * jeta[k][j][i].z
                    );
                }

                if(isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag))
                {

                    rhs[k][j][i].z
                    +=
                    scale *
                    (
                        meanGradP.x * kzet[k][j][i].x +
                        meanGradP.y * kzet[k][j][i].y +
                        meanGradP.z * kzet[k][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

    DMDAVecRestoreArray(fda, Rhs,            &rhs);
    DMDAVecRestoreArray(da,  mesh->lNvert,   &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode bulkGradPForcing(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM             da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo  info = mesh->info;
    PetscInt       xs   = info.xs, xe = info.xs + info.xm;
    PetscInt       ys   = info.ys, ye = info.ys + info.ym;
    PetscInt       zs   = info.zs, ze = info.zs + info.zm;
    PetscInt       mx   = info.mx, my = info.my, mz = info.mz;
    
    Cmpnts        ***rhs, ***ucat;
    Cmpnts        ***icsi, ***jeta, ***kzet;
    PetscReal     ***aj, ***nvert, ***meshTag, cellVol;
    PetscReal      lSumUx = 0.0, lSumUy = 0.0, lSumUz = 0.0;
    PetscReal      lVol   = 0.0;
    PetscReal      gSumUx = 0.0, gSumUy = 0.0, gSumUz = 0.0, gVol = 0.0;
    PetscReal      meanUx, meanUy, meanUz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l;
    
    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lAj,    &aj); 

    DMDAVecGetArray(fda, Rhs,          &rhs);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
   
    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                if(isIBMCell(k,j,i,nvert) || isOversetCell(k,j,i,meshTag)) continue;

                cellVol = 1.0 / aj[k][j][i];
                lSumUx += ucat[k][j][i].x * cellVol;
                lSumUy += ucat[k][j][i].y * cellVol;
                lSumUz += ucat[k][j][i].z * cellVol;
                lVol   += cellVol;
            }
        }
    }
    
    MPI_Allreduce(&lSumUx, &gSumUx, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lSumUy, &gSumUy, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lSumUz, &gSumUz, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lVol,   &gVol,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //volume averaged
    meanUx = gSumUx / gVol;
    meanUy = gSumUy / gVol;
    meanUz = gSumUz / gVol;

    Cmpnts meanU = nSetFromComponents(meanUx, meanUy, meanUz);
    Cmpnts diffU = nSub(ueqn->uBulk, meanU);
        
    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                if(isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag))
                {
                    rhs[k][j][i].x 
                    +=
                    scale *
                    (
                        diffU.x * icsi[k][j][i].x +
                        diffU.y * icsi[k][j][i].y +
                        diffU.z * icsi[k][j][i].z
                    );
                }

                if(isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag))
                {

                    rhs[k][j][i].y
                    +=
                    scale *
                    (
                        diffU.x * jeta[k][j][i].x +
                        diffU.y * jeta[k][j][i].y +
                        diffU.z * jeta[k][j][i].z
                    );
                }

                if(isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag))
                {

                    rhs[k][j][i].z
                    +=
                    scale *
                    (
                        diffU.x * kzet[k][j][i].x +
                        diffU.y * kzet[k][j][i].y +
                        diffU.z * kzet[k][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    DMDAVecRestoreArray(fda, Rhs,          &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode hyperViscosityU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    // based on Knievel et al 2007 and adapted to TOSCA's staggered-hybrid contravaiant-cartesian formulation 
    // 1. transform cell-centered cartesian velocity to contravariant 
    // 2. compute 6th order diffusion fluxes and apply them at cell centers 
    // 3. interpolate the cell centered rhs to the faces
    
    if(ueqn->hyperVisc4i == 0.0 && ueqn->hyperVisc4j == 0.0 && ueqn->hyperVisc4k == 0.0) return(0);

    mesh_        *mesh  = ueqn->access->mesh;
    clock_       *clock = ueqn->access->clock;
    DM            fda   = mesh->fda;
    DM            da    = mesh->da;
    DMDALocalInfo info  = mesh->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs = xs; if(lxs == 0) lxs=lxs+1; PetscInt lxe = xe; if(lxe == mx) lxe=lxe-1;
    PetscInt lys = ys; if(lys == 0) lys=lys+1; PetscInt lye = ye; if(lye == my) lye=lye-1;
    PetscInt lzs = zs; if(lzs == 0) lzs=lzs+1; PetscInt lze = ze; if(lze == mz) lze=lze-1;

    Cmpnts           ***csi,  ***eta,  ***zet;
    Cmpnts           ***icsi, ***ieta, ***izet;
    Cmpnts           ***jcsi, ***jeta, ***jzet;
    Cmpnts           ***kcsi, ***keta, ***kzet;

    PetscInt    i, j, k;
    Vec         Ucont_c, Rhs_c;
    Cmpnts      ***ucont, ***ucat, ***ucont_c, ***rhs, ***rhs_c;
    PetscReal   ***nvert;
    PetscScalar solid = 0.5;

    PetscReal eps_i = scale * 0.015625 * ueqn->hyperVisc4i;
    PetscReal eps_j = scale * 0.015625 * ueqn->hyperVisc4j;
    PetscReal eps_k = scale * 0.015625 * ueqn->hyperVisc4k;

    PetscReal Fp_x;
    PetscReal Fm_x;
    PetscReal Fp_y;
    PetscReal Fm_y;
    PetscReal Fp_z;
    PetscReal Fm_z;

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    VecDuplicate( ueqn->lUcont, &Ucont_c); VecSet(Ucont_c, 0.0);
    VecDuplicate( ueqn->lUcont, &Rhs_c); VecSet(Rhs_c, 0.0);
    DMDAVecGetArray(fda, Ucont_c, &ucont_c);

    // 1. transform cartesian velocity to contravariant at internal face centers 
    for(k = lzs; k < lze; k++)
    for(j = lys; j < lye; j++)
    for(i = lxs; i < lxe; i++)
    {
        ucont_c[k][j][i].x  = csi[k][j][i].x * ucat[k][j][i].x + csi[k][j][i].y * ucat[k][j][i].y + csi[k][j][i].z * ucat[k][j][i].z;
        ucont_c[k][j][i].y  = eta[k][j][i].x * ucat[k][j][i].x + eta[k][j][i].y * ucat[k][j][i].y + eta[k][j][i].z * ucat[k][j][i].z;
        ucont_c[k][j][i].z  = zet[k][j][i].x * ucat[k][j][i].x + zet[k][j][i].y * ucat[k][j][i].y + zet[k][j][i].z * ucat[k][j][i].z;

        // handle zerogradient and no penetration boundaries 
        if(!(mesh->i_periodic) && !(mesh->ii_periodic))
        {
            if (i == 1)
            {
                // generic zero gradient 
                ucont_c[k][j][i-1].x = ucont_c[k][j][i].x;
                ucont_c[k][j][i-1].y = ucont_c[k][j][i].y;
                ucont_c[k][j][i-1].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.iLeftPatchType == 1) // no penetration
                {
                    ucont_c[k][j][i-1].x = -ucont_c[k][j][i].x;
                }
            }
            else if(i == mx-2)
            {
                // generic zero gradient 
                ucont_c[k][j][i+1].x = ucont_c[k][j][i].x;
                ucont_c[k][j][i+1].y = ucont_c[k][j][i].y;
                ucont_c[k][j][i+1].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.iRightPatchType == 1) // no penetration
                {
                    ucont_c[k][j][i+1].x = -ucont_c[k][j][i].x;
                }
            }
        }

        if(!(mesh->j_periodic) && !(mesh->jj_periodic))
        {
            if (j == 1)
            {
                // generic zero gradient 
                ucont_c[k][j-1][i].x = ucont_c[k][j][i].x;
                ucont_c[k][j-1][i].y = ucont_c[k][j][i].y;
                ucont_c[k][j-1][i].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.jLeftPatchType == 1) // no penetration
                {
                    ucont_c[k][j-1][i].y = -ucont_c[k][j][i].y;
                }
            }
            else if(j == my-2)
            {
                // generic zero gradient 
                ucont_c[k][j+1][i].x = ucont_c[k][j][i].x;
                ucont_c[k][j+1][i].y = ucont_c[k][j][i].y;
                ucont_c[k][j+1][i].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.jRightPatchType == 1) // no penetration
                {
                    ucont_c[k][j+1][i].y = -ucont_c[k][j][i].y;
                }
            }
        }

        if(!(mesh->k_periodic) && !(mesh->kk_periodic))
        {
            if (k == 1)
            {
                // generic zero gradient 
                ucont_c[k-1][j][i].x = ucont_c[k][j][i].x;
                ucont_c[k-1][j][i].y = ucont_c[k][j][i].y;
                ucont_c[k-1][j][i].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.kLeftPatchType == 1) // no penetration
                {
                    ucont_c[k-1][j][i].z = -ucont_c[k][j][i].z;
                }
            }
            else if(k == mz-2)
            {
                // generic zero gradient 
                ucont_c[k+1][j][i].x = ucont_c[k][j][i].x;
                ucont_c[k+1][j][i].y = ucont_c[k][j][i].y;
                ucont_c[k+1][j][i].z = ucont_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.kRightPatchType == 1) // no penetration
                {
                    ucont_c[k+1][j][i].z = -ucont_c[k][j][i].z;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, Ucont_c, &ucont_c);
    DMLocalToLocalBegin(fda, Ucont_c, INSERT_VALUES, Ucont_c);
    DMLocalToLocalEnd  (fda, Ucont_c, INSERT_VALUES, Ucont_c);

    // handle periodic values and populate stencils
    resetCellPeriodicFluxes(mesh, Ucont_c, Ucont_c, "vector", "localToLocal");

    // 2. compute diffusion term at cell centers (stay on internal faces)
    DMDAVecGetArray(fda, Rhs_c, &rhs_c);
    DMDAVecGetArray(fda, Ucont_c, &ucont_c);

    lxs = xs; if(lxs == 0) lxs=lxs+1;  lxe = xe; if(lxe == mx) lxe=lxe-1;
    lys = ys; if(lys == 0) lys=lys+1;  lye = ye; if(lye == my) lye=lye-1;
    lzs = zs; if(lzs == 0) lzs=lzs+1;  lze = ze; if(lze == mz) lze=lze-1;

    for(k = lzs; k < lze; k++)
    for(j = lys; j < lye; j++)
    for(i = lxs; i < lxe; i++)
    {

        // damping i direction
        PetscInt im3, im2, im1, ip1, ip2, ip3;
        getCell2Cell7StencilCsi(mesh, i, mx, &im3, &im2, &im1, &ip1, &ip2, &ip3);

        PetscReal ux_ip1 = ucont_c[k][j][ip1].x, uy_ip1 = ucont_c[k][j][ip1].y, uz_ip1 = ucont_c[k][j][ip1].z,
                  ux_im1 = ucont_c[k][j][im1].x, uy_im1 = ucont_c[k][j][im1].y, uz_im1 = ucont_c[k][j][im1].z,
                  ux_ip2 = ucont_c[k][j][ip2].x, uy_ip2 = ucont_c[k][j][ip2].y, uz_ip2 = ucont_c[k][j][ip2].z,
                  ux_im2 = ucont_c[k][j][im2].x, uy_im2 = ucont_c[k][j][im2].y, uz_im2 = ucont_c[k][j][im2].z, 
                  ux_ip3 = ucont_c[k][j][ip3].x, uy_ip3 = ucont_c[k][j][ip3].y, uz_ip3 = ucont_c[k][j][ip3].z, 
                  ux_im3 = ucont_c[k][j][im3].x, uy_im3 = ucont_c[k][j][im3].y, uz_im3 = ucont_c[k][j][im3].z,
                  ux_i   = ucont_c[k][j][i  ].x, uy_i   = ucont_c[k][j][i  ].y, uz_i   = ucont_c[k][j][i  ].z;

        if(nvert[k][j][i] + nvert[k][j][i+1] < 2.0 * solid)
        {
            // sixth order 
            if(i > 2 && i < mx-3) 
            {
                Fp_x = 10.0*(ux_ip1 - ux_i) -5.0*(ux_ip2 - ux_im1) + (ux_ip3 - ux_im2);
                Fm_x = 10.0*(ux_i - ux_im1) -5.0*(ux_ip1 - ux_im2) + (ux_ip2 - ux_im3);
                Fp_y = 10.0*(uy_ip1 - uy_i) -5.0*(uy_ip2 - uy_im1) + (uy_ip3 - uy_im2);
                Fm_y = 10.0*(uy_i - uy_im1) -5.0*(uy_ip1 - uy_im2) + (uy_ip2 - uy_im3);
                Fp_z = 10.0*(uz_ip1 - uz_i) -5.0*(uz_ip2 - uz_im1) + (uz_ip3 - uz_im2);
                Fm_z = 10.0*(uz_i - uz_im1) -5.0*(uz_ip1 - uz_im2) + (uz_ip2 - uz_im3);
            }
            // fourth order
            else if(i > 1 && i < mx-2) 
            {
                Fp_x = -(ux_ip2 - ux_im1) + 3.0*(ux_ip1 - ux_i);
                Fm_x = -(ux_ip1 - ux_im2) + 3.0*(ux_i - ux_im1);
                Fp_y = -(uy_ip2 - uy_im1) + 3.0*(uy_ip1 - uy_i);
                Fm_y = -(uy_ip1 - uy_im2) + 3.0*(uy_i - uy_im1);
                Fp_z = -(uz_ip2 - uz_im1) + 3.0*(uz_ip1 - uz_i);
                Fm_z = -(uz_ip1 - uz_im2) + 3.0*(uz_i - uz_im1);
            }
            // second order
            else
            {
                Fp_x = ux_ip1 - ux_i;
                Fm_x = ux_i - ux_im1;
                Fp_y = uy_ip1 - uy_i;
                Fm_y = uy_i - uy_im1;
                Fp_z = uz_ip1 - uz_i;
                Fm_z = uz_i - uz_im1;
            }

            if((ux_ip1 - ux_i) * Fp_x <= 0.0) Fp_x = 0.0;
            if((ux_i - ux_im1) * Fm_x <= 0.0) Fm_x = 0.0;
            if((uy_ip1 - uy_i) * Fp_y <= 0.0) Fp_y = 0.0;
            if((uy_i - uy_im1) * Fm_y <= 0.0) Fm_y = 0.0;
            if((uz_ip1 - uz_i) * Fp_z <= 0.0) Fp_z = 0.0;
            if((uz_i - uz_im1) * Fm_z <= 0.0) Fm_z = 0.0;

            rhs_c[k][j][i].x += eps_i * (Fp_x - Fm_x);
            rhs_c[k][j][i].y += eps_i * (Fp_y - Fm_y);
            rhs_c[k][j][i].z += eps_i * (Fp_z - Fm_z);
        }   

        // damping in j direction 
        PetscInt jm3, jm2, jm1, jp1, jp2, jp3;
        getCell2Cell7StencilEta(mesh, j, my, &jm3, &jm2, &jm1, &jp1, &jp2, &jp3);

        PetscReal ux_jp1 = ucont_c[k][jp1][i].x, uy_jp1 = ucont_c[k][jp1][i].y, uz_jp1 = ucont_c[k][jp1][i].z,
                  ux_jm1 = ucont_c[k][jm1][i].x, uy_jm1 = ucont_c[k][jm1][i].y, uz_jm1 = ucont_c[k][jm1][i].z,
                  ux_jp2 = ucont_c[k][jp2][i].x, uy_jp2 = ucont_c[k][jp2][i].y, uz_jp2 = ucont_c[k][jp2][i].z,
                  ux_jm2 = ucont_c[k][jm2][i].x, uy_jm2 = ucont_c[k][jm2][i].y, uz_jm2 = ucont_c[k][jm2][i].z, 
                  ux_jp3 = ucont_c[k][jp3][i].x, uy_jp3 = ucont_c[k][jp3][i].y, uz_jp3 = ucont_c[k][jp3][i].z, 
                  ux_jm3 = ucont_c[k][jm3][i].x, uy_jm3 = ucont_c[k][jm3][i].y, uz_jm3 = ucont_c[k][jm3][i].z,
                  ux_j   = ucont_c[k][j  ][i].x, uy_j   = ucont_c[k][j  ][i].y, uz_j   = ucont_c[k][j  ][i].z;

        if(nvert[k][j][i] + nvert[k][j+1][i] < 2.0 * solid)
        {
            // sixth order
            if(j > 2 && j < my-3)
            {
                Fp_x = 10.0*(ux_jp1 - ux_j) -5.0*(ux_jp2 - ux_jm1) + (ux_jp3 - ux_jm2);
                Fm_x = 10.0*(ux_j - ux_jm1) -5.0*(ux_jp1 - ux_jm2) + (ux_jp2 - ux_jm3);
                Fp_y = 10.0*(uy_jp1 - uy_j) -5.0*(uy_jp2 - uy_jm1) + (uy_jp3 - uy_jm2);
                Fm_y = 10.0*(uy_j - uy_jm1) -5.0*(uy_jp1 - uy_jm2) + (uy_jp2 - uy_jm3);
                Fp_z = 10.0*(uz_jp1 - uz_j) -5.0*(uz_jp2 - uz_jm1) + (uz_jp3 - uz_jm2);
                Fm_z = 10.0*(uz_j - uz_jm1) -5.0*(uz_jp1 - uz_jm2) + (uz_jp2 - uz_jm3);
            }
            // fourth order
            else if(j > 1 && j < my-2)
            {
                Fp_x = -(ux_jp2 - ux_jm1) + 3.0*(ux_jp1 - ux_j);
                Fm_x = -(ux_jp1 - ux_jm2) + 3.0*(ux_j - ux_jm1);
                Fp_y = -(uy_jp2 - uy_jm1) + 3.0*(uy_jp1 - uy_j);
                Fm_y = -(uy_jp1 - uy_jm2) + 3.0*(uy_j - uy_jm1);
                Fp_z = -(uz_jp2 - uz_jm1) + 3.0*(uz_jp1 - uz_j);
                Fm_z = -(uz_jp1 - uz_jm2) + 3.0*(uz_j - uz_jm1);
            }
            // second order
            else
            {
                Fp_x = ux_jp1 - ux_j;
                Fm_x = ux_j - ux_jm1;
                Fp_y = uy_jp1 - uy_j;
                Fm_y = uy_j - uy_jm1;
                Fp_z = uz_jp1 - uz_j;
                Fm_z = uz_j - uz_jm1;
            }

            if((ux_jp1 - ux_j) * Fp_x <= 0.0) Fp_x = 0.0;
            if((ux_j - ux_jm1) * Fm_x <= 0.0) Fm_x = 0.0;
            if((uy_jp1 - uy_j) * Fp_y <= 0.0) Fp_y = 0.0;
            if((uy_j - uy_jm1) * Fm_y <= 0.0) Fm_y = 0.0;
            if((uz_jp1 - uz_j) * Fp_z <= 0.0) Fp_z = 0.0;
            if((uz_j - uz_jm1) * Fm_z <= 0.0) Fm_z = 0.0;

            rhs_c[k][j][i].x += eps_j * (Fp_x - Fm_x);
            rhs_c[k][j][i].y += eps_j * (Fp_y - Fm_y);
            rhs_c[k][j][i].z += eps_j * (Fp_z - Fm_z);
        }

        // damping k direction
        PetscInt km3, km2, km1, kp1, kp2, kp3;
        getCell2Cell7StencilZet(mesh, k, mz, &km3, &km2, &km1, &kp1, &kp2, &kp3);

        PetscReal ux_kp1 = ucont_c[kp1][j][i].x, uy_kp1 = ucont_c[kp1][j][i].y, uz_kp1 = ucont_c[kp1][j][i].z,
                  ux_km1 = ucont_c[km1][j][i].x, uy_km1 = ucont_c[km1][j][i].y, uz_km1 = ucont_c[km1][j][i].z,
                  ux_kp2 = ucont_c[kp2][j][i].x, uy_kp2 = ucont_c[kp2][j][i].y, uz_kp2 = ucont_c[kp2][j][i].z,
                  ux_km2 = ucont_c[km2][j][i].x, uy_km2 = ucont_c[km2][j][i].y, uz_km2 = ucont_c[km2][j][i].z, 
                  ux_kp3 = ucont_c[kp3][j][i].x, uy_kp3 = ucont_c[kp3][j][i].y, uz_kp3 = ucont_c[kp3][j][i].z, 
                  ux_km3 = ucont_c[km3][j][i].x, uy_km3 = ucont_c[km3][j][i].y, uz_km3 = ucont_c[km3][j][i].z,
                  ux_k   = ucont_c[k  ][j][i].x, uy_k   = ucont_c[k  ][j][i].y, uz_k   = ucont_c[k  ][j][i].z;

        if(nvert[k][j][i] + nvert[k+1][j][i] < 2.0 * solid)
        {
            if(k > 2 && k < mz-3)
            {
                Fp_x = 10.0*(ux_kp1 - ux_k) -5.0*(ux_kp2 - ux_km1) + (ux_kp3 - ux_km2);
                Fm_x = 10.0*(ux_k - ux_km1) -5.0*(ux_kp1 - ux_km2) + (ux_kp2 - ux_km3);
                Fp_y = 10.0*(uy_kp1 - uy_k) -5.0*(uy_kp2 - uy_km1) + (uy_kp3 - uy_km2);
                Fm_y = 10.0*(uy_k - uy_km1) -5.0*(uy_kp1 - uy_km2) + (uy_kp2 - uy_km3);
                Fp_z = 10.0*(uz_kp1 - uz_k) -5.0*(uz_kp2 - uz_km1) + (uz_kp3 - uz_km2);
                Fm_z = 10.0*(uz_k - uz_km1) -5.0*(uz_kp1 - uz_km2) + (uz_kp2 - uz_km3);  
            }
            else if(k > 1 && k < mz-2)
            {
                Fp_x = -(ux_kp2 - ux_km1) + 3.0*(ux_kp1 - ux_k);
                Fm_x = -(ux_kp1 - ux_km2) + 3.0*(ux_k - ux_km1);
                Fp_y = -(uy_kp2 - uy_km1) + 3.0*(uy_kp1 - uy_k);
                Fm_y = -(uy_kp1 - uy_km2) + 3.0*(uy_k - uy_km1);
                Fp_z = -(uz_kp2 - uz_km1) + 3.0*(uz_kp1 - uz_k);
                Fm_z = -(uz_kp1 - uz_km2) + 3.0*(uz_k - uz_km1);
            }
            else
            {
                Fp_x = ux_kp1 - ux_k;
                Fm_x = ux_k - ux_km1;
                Fp_y = uy_kp1 - uy_k;
                Fm_y = uy_k - uy_km1;
                Fp_z = uz_kp1 - uz_k;
                Fm_z = uz_k - uz_km1;
            } 

            if((ux_kp1 - ux_k) * Fp_x <= 0.0) Fp_x = 0.0;
            if((ux_k - ux_km1) * Fm_x <= 0.0) Fm_x = 0.0;
            if((uy_kp1 - uy_k) * Fp_y <= 0.0) Fp_y = 0.0;
            if((uy_k - uy_km1) * Fm_y <= 0.0) Fm_y = 0.0;
            if((uz_kp1 - uz_k) * Fp_z <= 0.0) Fp_z = 0.0;
            if((uz_k - uz_km1) * Fm_z <= 0.0) Fm_z = 0.0;

            rhs_c[k][j][i].x += eps_k * (Fp_x - Fm_x);
            rhs_c[k][j][i].y += eps_k * (Fp_y - Fm_y);
            rhs_c[k][j][i].z += eps_k * (Fp_z - Fm_z);
        }

        // handle zerogradient and no penetration boundaries 
        if(!(mesh->i_periodic) && !(mesh->ii_periodic))
        {
            if (i == 1)
            {
                // generic zero gradient 
                rhs_c[k][j][i-1].x = rhs_c[k][j][i].x;
                rhs_c[k][j][i-1].y = rhs_c[k][j][i].y;
                rhs_c[k][j][i-1].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.iLeftPatchType == 1) // no penetration
                {
                    rhs_c[k][j][i-1].x = -rhs_c[k][j][i].x;
                }
            }
            else if(i == mx-2)
            {
                // generic zero gradient 
                rhs_c[k][j][i+1].x = rhs_c[k][j][i].x;
                rhs_c[k][j][i+1].y = rhs_c[k][j][i].y;
                rhs_c[k][j][i+1].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.iRightPatchType == 1) // no penetration
                {
                    rhs_c[k][j][i+1].x = -rhs_c[k][j][i].x;
                }
            }
        }

        if(!(mesh->j_periodic) && !(mesh->jj_periodic))
        {
            if (j == 1)
            {
                // generic zero gradient 
                rhs_c[k][j-1][i].x = rhs_c[k][j][i].x;
                rhs_c[k][j-1][i].y = rhs_c[k][j][i].y;
                rhs_c[k][j-1][i].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.jLeftPatchType == 1) // no penetration
                {
                    rhs_c[k][j-1][i].y = -rhs_c[k][j][i].y;
                }
            }
            else if(j == my-2)
            {
                // generic zero gradient 
                rhs_c[k][j+1][i].x = rhs_c[k][j][i].x;
                rhs_c[k][j+1][i].y = rhs_c[k][j][i].y;
                rhs_c[k][j+1][i].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.jRightPatchType == 1) // no penetration
                {
                    rhs_c[k][j+1][i].y = -rhs_c[k][j][i].y;
                }
            }
        }

        if(!(mesh->k_periodic) && !(mesh->kk_periodic))
        {
            if (k == 1)
            {
                // generic zero gradient 
                rhs_c[k-1][j][i].x = rhs_c[k][j][i].x;
                rhs_c[k-1][j][i].y = rhs_c[k][j][i].y;
                rhs_c[k-1][j][i].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.kLeftPatchType == 1) // no penetration
                {
                    rhs_c[k-1][j][i].z = -rhs_c[k][j][i].z;
                }
            }
            else if(k == mz-2)
            {
                // generic zero gradient 
                rhs_c[k+1][j][i].x = rhs_c[k][j][i].x;
                rhs_c[k+1][j][i].y = rhs_c[k][j][i].y;
                rhs_c[k+1][j][i].z = rhs_c[k][j][i].z;

                // no penetration 
                if(mesh->boundaryU.kRightPatchType == 1) // no penetration
                {
                    rhs_c[k+1][j][i].z = -rhs_c[k][j][i].z;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, Ucont_c, &ucont_c);
    DMDAVecRestoreArray(fda, Rhs_c, &rhs_c);
    DMLocalToLocalBegin(fda, Rhs_c, INSERT_VALUES, Rhs_c);
    DMLocalToLocalEnd  (fda, Rhs_c, INSERT_VALUES, Rhs_c);

    resetCellPeriodicFluxes(mesh, Rhs_c, Rhs_c, "vector", "localToLocal");

    // 3. interpolate back to faces and add to rhs
    DMDAVecGetArray(fda, Rhs,     &rhs);
    DMDAVecGetArray(fda, Rhs_c, &rhs_c);

    lxs = xs; if(lxs == 0) lxs=lxs+1;  lxe = xe; if(lxe == mx) lxe=lxe-1;
    lys = ys; if(lys == 0) lys=lys+1;  lye = ye; if(lye == my) lye=lye-1;
    lzs = zs; if(lzs == 0) lzs=lzs+1;  lze = ze; if(lze == mz) lze=lze-1;

    for(k = lzs; k < lze; k++)
    for(j = lys; j < lye; j++)
    for(i = lxs; i < lxe; i++)
    {
            rhs[k][j][i].x += 0.5 * (rhs_c[k][j][i+1].x + rhs_c[k][j][i].x);
            rhs[k][j][i].y += 0.5 * (rhs_c[k][j+1][i].y + rhs_c[k][j][i].y);
            rhs[k][j][i].z += 0.5 * (rhs_c[k+1][j][i].z + rhs_c[k][j][i].z);
    }

    DMDAVecRestoreArray(fda, Rhs_c, &rhs_c);
    DMDAVecRestoreArray(fda, Rhs,     &rhs);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    VecDestroy(&Ucont_c);
    VecDestroy(&Rhs_c);

    return(0);
}

//***************************************************************************************************************//

