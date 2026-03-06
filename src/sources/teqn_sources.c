//! \file  teqn_sources.c
//! \brief Temperature-equation source-term function implementations.

#include "teqn_sources.h"

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
    PetscReal *ltMean, *gtMean, *tD, *tM;
    PetscReal *src;

    if(abl->controllerActionT == "write")
    {
        PetscMalloc(sizeof(PetscReal) * (nLevels), &(ltMean));
        PetscMalloc(sizeof(PetscReal) * (nLevels), &(gtMean));
        PetscMalloc(sizeof(PetscReal) * (nLevels), &(tD));
        PetscMalloc(sizeof(PetscReal) * (nLevels), &(tM));

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
        else if((abl->controllerTypeT == "directProfileAssimilation") || (abl->controllerTypeT == "indirectProfileAssimilation") || (abl->controllerTypeT == "waveletProfileAssimilation"))
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

            // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
            if(abl->timeT[idx_2] < clock->time)
            {
                idx_1 = idx_2;
                idx_2 = idx_1 + 1;
            }

            if(abl->timeT[idx_1] > clock->time)
            {
                idx_2 = idx_1;
                idx_1 = idx_2 - 1;
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

            for(j=0; j<nLevels; j++)
            {
                // if cellLevels is below the lowest mesoscale data point use the source at the last available height
                if(abl->cellLevels[j] < abl->bottomSrcHtT)
                {
                    abl->tDes[j] = abl->tDes[abl->lowestIndT];
                    gtMean[j] = gtMean[abl->lowestIndT];
                    tD[j]     = tD[abl->lowestIndT];
                    tM[j]     = tM[abl->lowestIndT];
                }

                if(abl->cellLevels[j] > abl->hT[abl->numhT - 1])
                {
                    abl->tDes[j] = abl->tDes[abl->highestIndT];
                    gtMean[j] = gtMean[abl->highestIndT];
                    tD[j]     = tD[abl->highestIndT];
                    tM[j]     = tM[abl->highestIndT];
                }
            }

            if((abl->controllerTypeT == "indirectProfileAssimilation"))
            {
                //smooth the source based on the polynomial interpolation 
                for (PetscInt iCtr = 0; iCtr < nLevels; iCtr++) 
                {
                    PetscReal dotProduct1 = 0.0, dotProduct2 = 0.0;

                    for (PetscInt kCtr = 0; kCtr < nLevels; kCtr++) 
                    {
                        dotProduct1 += abl->polyCoeffT[iCtr][kCtr] * tD[kCtr];
                        dotProduct2 += abl->polyCoeffT[iCtr][kCtr] * tM[kCtr];
                    }
                    abl->tDes[iCtr] = dotProduct1;
                    gtMean[iCtr]    = dotProduct2;
                }  

                if(abl->flTypeT == "constantHeight")
                {
                    for(j=0; j<nLevels; j++)    
                    {
                        if(abl->cellLevels[j] < abl->bottomSrcHtT)
                        {
                            abl->tDes[j] = abl->tDes[abl->lowestIndT];
                            gtMean[j] = gtMean[abl->lowestIndT];
                        }
                    }
                }

                for(j=0; j<nLevels; j++)    
                {
                    if(abl->cellLevels[j] > abl->hT[abl->numhT - 1])
                    {
                        abl->tDes[j] = abl->tDes[abl->highestIndT];
                        gtMean[j] = gtMean[abl->highestIndT];
                    }
                }

            }

            if((abl->controllerTypeT == "waveletProfileAssimilation"))
            {
                if(abl->waveletBlend)
                {
                    PetscReal *tDesAbove, *gtMeanAbove;
                    PetscReal *tDesBelow, *gtMeanBelow;

                    PetscMalloc(sizeof(PetscReal) * (nLevels), &(tDesAbove));
                    PetscMalloc(sizeof(PetscReal) * (nLevels), &(gtMeanAbove));
                    PetscMalloc(sizeof(PetscReal) * (nLevels), &(tDesBelow));
                    PetscMalloc(sizeof(PetscReal) * (nLevels), &(gtMeanBelow));

                    #if USE_PYTHON
                        pywavedecScalar(abl, tD, tDesBelow, nLevels);
                        pywavedecScalar(abl, tM, gtMeanBelow, nLevels);
                    #endif

                    PetscInt waveLevel = abl->waveLevel;
                    
                    //hardcoded value for now
                    abl->waveLevel = abl->waveLevel - 4;
                    if(abl->waveLevel <= 0) abl->waveLevel = 1;

                    #if USE_PYTHON
                        pywavedecScalar(abl, tD, tDesAbove, nLevels);
                        pywavedecScalar(abl, tM, gtMeanAbove, nLevels);
                    #endif

                    abl->waveLevel = waveLevel;
                    
                    // blend between the below and above profiles 
                    PetscReal blendHeight = mesh->bounds.zmax-500.0, scaleFactor;
                    PetscReal ablgeoDelta = 300;
                    PetscReal blendTop = blendHeight + 0.5*ablgeoDelta;

                    for(j=0; j<nLevels; j++)
                    {
                        scaleFactor = scaleHyperTangTop(abl->cellLevels[j], blendTop, ablgeoDelta);
                        abl->tDes[j] = (1-scaleFactor)*tDesBelow[j] + scaleFactor*tDesAbove[j];                                                                     
                        gtMean[j] = (1-scaleFactor)*gtMeanBelow[j] + scaleFactor*gtMeanAbove[j];                                                                
                    }   

                    PetscFree(tDesAbove);
                    PetscFree(gtMeanAbove);
                    PetscFree(tDesBelow);
                    PetscFree(gtMeanBelow);
                }
                else 
                {
                    #if USE_PYTHON
                    pywavedecScalar(abl, tD, abl->tDes, nLevels);
                    pywavedecScalar(abl, tM, gtMean, nLevels);
                    #endif
                }


                if(abl->flTypeT == "constantHeight")
                {
                    for(j=0; j<nLevels; j++)    
                    {
                        if(abl->cellLevels[j] < abl->bottomSrcHtT)
                        {
                            abl->tDes[j] = abl->tDes[abl->lowestIndT];
                            gtMean[j] = gtMean[abl->lowestIndT];
                        }
                    }
                }

                for(j=0; j<nLevels; j++)    
                {
                    if(abl->cellLevels[j] > abl->hT[abl->numhT - 1])
                    {
                        abl->tDes[j] = abl->tDes[abl->highestIndT];
                        gtMean[j] = gtMean[abl->highestIndT];
                    }
                }
            }
        }
    }
    else if(abl->controllerActionT == "read")
    {
        if(abl->controllerTypeT == "timeHeightSeries")
        {
            PetscMalloc(sizeof(PetscReal) * (nLevels), &(src));
            PetscReal  tT1, tT2;
            PetscReal  wth1, wth2;
            PetscInt   idxh1, idxh2;

            PetscInt idx_1 = abl->currentCloseIdxT;
            PetscInt idx_2 = abl->currentCloseIdxT + 1;

            PetscInt lwrBound = 0;
            PetscInt uprBound = abl->nSourceTimes;

            // if past first iteration do the search on a subset to speed up the process
            if(clock->it > clock->itStart)
            {
                lwrBound = PetscMax(0, (abl->currentCloseIdxT - 50));
                uprBound = PetscMin(abl->nSourceTimes, (abl->currentCloseIdxT + 50));
            }

            // build error vector for the time search
            PetscReal  diff[abl->nSourceTimes];

            for(PetscInt i=lwrBound; i<uprBound; i++)
            {
                diff[i] = fabs(abl->timeHtSourcesT[i][0][0] - clock->time);
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
            if(abl->timeHtSourcesT[idx_2][0][0] < abl->timeHtSourcesT[idx_1][0][0])
            {
                PetscInt idx_tmp = idx_2;
                idx_2 = idx_1;
                idx_1 = idx_tmp;
            }

            // ensure time is bounded between idx_1 and idx_2 (necessary for non-uniform data)
            if(abl->timeHtSourcesT[idx_2][0][0] < clock->time)
            {
                idx_1 = idx_2;
                idx_2 = idx_1 + 1;
            }

            if(abl->timeHtSourcesT[idx_1][0][0] > clock->time)
            {
                idx_2 = idx_1;
                idx_1 = idx_2 - 1;
            }
            
            // find interpolation weights
            PetscReal idx = (idx_2 - idx_1) / (abl->timeHtSourcesT[idx_2][0][0] - abl->timeHtSourcesT[idx_1][0][0]) * (clock->time - abl->timeHtSourcesT[idx_1][0][0]) + idx_1;
            PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
            PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

                        
            if(print) PetscPrintf(mesh->MESH_COMM, "Correcting source terms T: selected time %lf for reading sources\n", w1 * abl->timeHtSourcesT[idx_1][0][0] + w2 * abl->timeHtSourcesT[idx_2][0][0]);
            if(print) PetscPrintf(mesh->MESH_COMM, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
            if(print) PetscPrintf(mesh->MESH_COMM, "                         closest avail. times : t1 = %lf, t2 = %lf\n", abl->timeHtSourcesT[idx_1][0][0], abl->timeHtSourcesT[idx_2][0][0]);

            // reset the closest index for nex iteration
            abl->currentCloseIdxT = idx_1;

            // get also the dt at the time the source was calculated
            double dtSource1 = std::max((abl->timeHtSourcesT[idx_1+1][0][0] - abl->timeHtSourcesT[idx_1][0][0]), 1e-5);
            double dtSource2 = std::max((abl->timeHtSourcesT[idx_2+1][0][0] - abl->timeHtSourcesT[idx_2][0][0]), 1e-5);

            // scale with time step for each height
            for(j=0; j<nLevels; j++)
            {
                if(abl->numhT != nLevels)
                {
                    //interpolating in height for each time
                    idxh1 = abl->tempInterpIdx[j][0];
                    idxh2 = abl->tempInterpIdx[j][1];

                    wth1  = abl->tempInterpWts[j][0];
                    wth2  = abl->tempInterpWts[j][1];
                    
                    tT1   = wth1 * abl->timeHtSourcesT[idx_1][idxh1][1] + wth2 * abl->timeHtSourcesT[idx_1][idxh2][1];
                    tT2   = wth1 * abl->timeHtSourcesT[idx_2][idxh1][1] + wth2 * abl->timeHtSourcesT[idx_2][idxh2][1];

                    // src[j] = (w1 * tT1 / dtSource1 + w2 * tT2 / dtSource2) * clock->dt;
                    src[j] = (w1 * tT1 + w2 * tT2);
                }
                else 
                {
                    // src[j] = (w1 * abl->timeHtSourcesT[idx_1][j][1] / dtSource1 + w2 * abl->timeHtSourcesT[idx_2][j][1] / dtSource2) * clock->dt;

                    src[j] = w1 * abl->timeHtSourcesT[idx_1][j][1] + w2 * abl->timeHtSourcesT[idx_2][j][1];

                } 
            }
        }
        else 
        {
            char error[512];
            sprintf(error, "unknown controllerType %s for T equation under read mode, available types are\n    1. timeHeightSeries\n   ", abl->controllerTypeT.c_str());
            fatalErrorInFunction("correctSourceTermT", error); 
        }
    }
    else
    {
        char error[512];
        sprintf(error, "unknown controllerAction %s for T equation, available actions are\n    1. write\n    2. read", abl->controllerAction.c_str());
        fatalErrorInFunction("correctSourceTermT", error);
    }

    if(abl->controllerActionT == "write")
    {
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
                    source[k][j][i]  = lsource[j-lys] * clock->dt;
                    lrelError += lsource[j-lys] / abl->tDes[j-1];
                    lcells++;
                }
            }
        }

        MPI_Allreduce(&lrelError, &grelError, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&lcells, &gcells, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

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
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e\t", width, abl->relax * (abl->tDes[j] - gtMean[j]));
                }

                PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                fclose(fp);
            }

            fileName = "postProcessing/temperatureSourceDPA_" + getStartTimeName(clock);
            fp = fopen(fileName.c_str(), "a");

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
                    PetscFPrintf(mesh->MESH_COMM, fp, "%*.5e\t", width, abl->relax * (tD[j] - tM[j]));
                }

                PetscFPrintf(mesh->MESH_COMM, fp, "\n");

                fclose(fp);
            }
        }
    }
    else if(abl->controllerActionT == "read")
    {
        if(abl->controllerTypeT == "timeHeightSeries")
        {
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        source[k][j][i]  = src[j-1] * clock->dt;
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da, teqn->sourceT, &source);
    DMDAVecRestoreArray(da, teqn->lTmprt,  &t);

    if(abl->controllerActionT == "write")
    {
        PetscFree(tD);
        PetscFree(tM);
        PetscFree(ltMean);
        PetscFree(gtMean);
    }
    else if(abl->controllerActionT == "read")
    {
        PetscFree(src);
    }

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


