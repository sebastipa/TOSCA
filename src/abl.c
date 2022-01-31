//! \file  abl.c
//! \brief Contains atmospheric boundary layer definition functions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"

//***************************************************************************************************************//

PetscErrorCode InitializeABL(abl_ *abl)
{
  if(abl != NULL)
  {
    mesh_         *mesh = abl->access->mesh;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscMPIInt   rank, nProcs;
    PetscInt      i, j, k, l;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Vec           Coor;

    PetscReal     lxmin, lymin, lzmin,       // local min coordinates
                  lxmax, lymax, lzmax,       // local max coordinates
                  xmin, ymin, zmin,          // global min coordinates
                  xmax, ymax, zmax;          // global max coordinates

    Cmpnts        ***coor, ***cent;          // point and cell center coordinates
    PetscReal     ***aj;

    char          dataLoc[256], fileName[500];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

    if(mesh->meshFileType != "cartesian")
    {
       char error[512];
        sprintf(error, "ABL capabilites only available for cartesian meshes\n");
        fatalErrorInFunction("InitializeABL",  error);
    }

    readDictDouble("ABLProperties.dat", "uTau",             &(abl->uTau));
    readDictDouble("ABLProperties.dat", "hRough",           &(abl->hRough));
    readDictDouble("ABLProperties.dat", "uRef",             &(abl->uRef));
    readDictDouble("ABLProperties.dat", "hRef",             &(abl->hRef));
    readDictDouble("ABLProperties.dat", "hInv",             &(abl->hInv));
    readDictDouble("ABLProperties.dat", "dInv",             &(abl->dInv));
    readDictDouble("ABLProperties.dat", "gInv",             &(abl->gInv));
    readDictDouble("ABLProperties.dat", "tRef",             &(abl->tRef));
    readDictDouble("ABLProperties.dat", "gTop",             &(abl->gTop));
    readDictDouble("ABLProperties.dat", "vkConst",          &(abl->vkConst));
    readDictDouble("ABLProperties.dat", "smearT",           &(abl->smear));
    readDictDouble("ABLProperties.dat", "controllerHeight", &(abl->controllerHeight));
    readDictDouble("ABLProperties.dat", "fCoriolis",        &(abl->fc));
    readDictWord  ("ABLProperties.dat", "controllerType",   &(abl->controllerType));

    // the vertical direction is the j direction in curvilinear coordinates
    PetscInt nLevels = my-2;

    // initialize some useful parameters used in fringe and velocity controller 'write'
    {
        PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->cellLevels));
        PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->totVolPerLevel));
        PetscMalloc(sizeof(PetscInt)  * nLevels, &(abl->totCelPerLevel));

        // initialize height levels for the velocity controller
        DMDAVecGetArray(fda, mesh->lCent, &cent);
        DMDAVecGetArray(da,  mesh->lAj,  &aj);

        std::vector<PetscReal> lLevels(nLevels);
        std::vector<PetscReal> gLevels(nLevels);
        std::vector<PetscReal> lVolumes(nLevels);
        std::vector<PetscReal> gVolumes(nLevels);
        std::vector<PetscInt>  lCells(nLevels);
        std::vector<PetscInt>  gCells(nLevels);

        for(l=0; l<nLevels; l++)
        {
            lLevels[l]  = 0.0;
            gLevels[l]  = 0.0;
            lVolumes[l] = 0.0;
            gVolumes[l] = 0.0;
            lCells[l]   = 0;
            gCells[l]   = 0;
        }

        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    lLevels[j-1]  += (cent[k][j][i].z - mesh->bounds.zmin);
                    lVolumes[j-1] += 1.0 / aj[k][j][i];
                    lCells[j-1]++;
                }
            }
        }

        MPI_Allreduce(&lLevels[0], &gLevels[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&lVolumes[0], &gVolumes[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
        MPI_Allreduce(&lCells[0], &gCells[0], nLevels, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

        for(l=0; l<nLevels; l++)
        {
            gLevels[l] = gLevels[l] / gCells[l];
        }

        DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
        DMDAVecRestoreArray(fda, mesh->lCent, &cent);

        for(l=0; l<nLevels; l++)
        {
            abl->cellLevels[l]     = gLevels[l];
            abl->totVolPerLevel[l] = gVolumes[l];
            abl->totCelPerLevel[l] = gCells[l];
        }

        std::vector<PetscReal> absLevelDelta(nLevels);

        for(l=0; l<nLevels; l++)
        {
            absLevelDelta[l] = std::fabs(abl->cellLevels[l] - abl->hRef);
        }

        for(PetscInt errI=0; errI<2; errI++)
        {
            PetscReal errMin   = 1e20;
            PetscReal errValue = 0.0;
            PetscInt  minLabel = 0;

            for(PetscInt errJ = errI; errJ < nLevels; errJ++)
            {
                if(absLevelDelta[errJ] < errMin)
                {
                    errValue = absLevelDelta[errJ];
                    minLabel = errJ;
                    errMin   = errValue;
                }
            }

            // exchange values so that elements are not ovwerwritten
            absLevelDelta[minLabel] = absLevelDelta[errI];

            // put the min value on the unchanged part at the last index of changed part
            absLevelDelta[errI] = errValue;

            // save the label adding one since DMDA labeling starts from physical ghost cells
            abl->closestLabels[errI] = minLabel + 1;
        }

        abl->levelWeights[0] = (abl->cellLevels[abl->closestLabels[1]-1]-abl->hRef) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);
        abl->levelWeights[1] = (abl->hRef-abl->cellLevels[abl->closestLabels[0]-1]) / (abl->cellLevels[abl->closestLabels[1]-1] - abl->cellLevels[abl->closestLabels[0]-1]);

        std::vector<PetscReal> ().swap(lLevels);
        std::vector<PetscReal> ().swap(gLevels);
        std::vector<PetscReal> ().swap(lVolumes);
        std::vector<PetscReal> ().swap(gVolumes);
        std::vector<PetscInt>  ().swap(lCells);
        std::vector<PetscInt>  ().swap(gCells);
        std::vector<PetscReal> ().swap(absLevelDelta);
    }

    // set cumulated sources to zero. They are needed for the integral part of the
    // controller if controllerType is set to 'write' or to store the average
    // sources if controllerType is set to 'average'
    abl->cumulatedSource.x = 0.0;
    abl->cumulatedSource.y = 0.0;
    abl->cumulatedSource.z = 0.0;

    // read the Rayleigh damping layer properties
    if(mesh->access->flags->isZDampingActive)
    {
        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingStart",   &(abl->zDampingStart));
        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingEnd",     &(abl->zDampingEnd));
        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingAlpha",   &(abl->zDampingAlpha));
        readSubDictInt   ("ABLProperties.dat", "zDampingProperties", "zDampingAlsoXY",  &(abl->zDampingAlsoXY));

        if(abl->zDampingAlsoXY)
        {
            // set average weight to zero
            abl->avgWeight = 0;

            // allocate memory for uBarMeanZ and set to zero
            PetscMalloc(my*sizeof(Cmpnts), &(abl->uBarMeanZ));

            for(j=0; j<my; j++)
            {
                abl->uBarMeanZ[j].x = 0.0;
                abl->uBarMeanZ[j].y = 0.0;
                abl->uBarMeanZ[j].z = 0.0;
            }
        }
    }

    // read the recycling fringe region properties
    if(mesh->access->flags->isXDampingActive)
    {
        readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingStart",            &(abl->xDampingStart));
        readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingEnd",              &(abl->xDampingEnd));
        readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingDelta",            &(abl->xDampingDelta));
        readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingAlpha",            &(abl->xDampingAlpha));
        readSubDictWord  ("ABLProperties.dat", "xDampingProperties", "xDampingAlphaControlType", &(abl->xDampingControlType));

        // if fringe controller is alphaOptimized read parameters
        if(abl->xDampingControlType == "alphaOptimized")
        {
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingLineSamplingYmin",  &(abl->xDampingLineSamplingYmin));
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingLineSamplingYmax",  &(abl->xDampingLineSamplingYmax));
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "xDampingTimeWindow",        &(abl->xDampingTimeWindow));
        }

        // read type of fringe region in uBarSelectionType
        // 1. periodized mapped
        // 2. interpolated mapped
        // 3. concurrent precursor (another simulation is started in parallel and we don't use inflow functions)

        readSubDictInt("ABLProperties.dat", "xDampingProperties", "uBarSelectionType", &(abl->xFringeUBarSelectionType));

        if(abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2)
        {
            // in x damping the uBar varies in time and it is actually the same as
            // the inflowFunction type 3 and 4, which correspond to uBarSelectionType
            // equal to 1 or 2 respectively. We use the same data structure as the inlet
            // function to store the informations, since it will never happen that both
            // inflow function and xDampingLayer are active (it must used with periodic
            // boundary conditions). Anyways, we make a check just to be sure.

            if(mesh->boundaryU.kLeft == "inletFunction")
            {
               char error[512];
                sprintf(error, "fringe region in x direction and inflowFunctions cannot be used together\n");
                fatalErrorInFunction("ABLInitialize", error);
            }

            // allocate memory for this patch inlet function data
            PetscMalloc(sizeof(inletFunctionTypes), &(mesh->inletF.kLeft));

            // set pointer to this inlet function type
            inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

            // set inflowFunction internal uBarSelectionType
            ifPtr->typeU = abl->xFringeUBarSelectionType;

            // unsteady mapped ubar
            if (ifPtr->typeU == 1)
            {
                // set tBarSelectionType the same as uBarSelectionType
                ifPtr->typeT  = 1;
                ifPtr->mapT   = 1;

                // nut value is dependent on the velocity (no need to map)
                ifPtr->mapNut = 0;

                readSubDictInt("ABLProperties.dat", "xDampingProperties", "n1Inflow",  &(ifPtr->n1));
                readSubDictInt("ABLProperties.dat", "xDampingProperties", "n2Inflow",  &(ifPtr->n2));
                readSubDictInt("ABLProperties.dat", "xDampingProperties", "n1Periods", &(ifPtr->prds1));
                readSubDictInt("ABLProperties.dat", "xDampingProperties", "n2Periods", &(ifPtr->prds2));

                // increase n1 and n2 accounting for side ghost cells
                ifPtr->n1wg = ifPtr->n1 + 2;
                ifPtr->n2wg = ifPtr->n2 + 2;

                // print mapping action information
                printInflowMappingAction(mesh, ifPtr);

                // initialize inflow data
                mappedInflowInitialize(ifPtr);
            }
            // unsteady interpolated ubar
            else if (ifPtr->typeU == 2)
            {
                // set tBarSelectionType the same as uBarSelectionType
                ifPtr->typeT  = 2;
                ifPtr->mapT   = 1;

                // nut value is dependent on the velocity (no need to map)
                ifPtr->mapNut = 0;

                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Inflow",   &(ifPtr->n1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Inflow",   &(ifPtr->n2));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Periods",  &(ifPtr->prds1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Periods",  &(ifPtr->prds2));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth1", &(ifPtr->width1));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth2", &(ifPtr->width2));

                // increase n1 and n2 accounting for side ghost cells
                ifPtr->n1wg = ifPtr->n1 + 2;
                ifPtr->n2wg = ifPtr->n2 + 2;

                // build interpolation weights and find periodized inflow cells indices
                SetInflowWeights(mesh, ifPtr);

                // initialize inflow data
                mappedInflowInitialize(ifPtr);
            }
            else
            {
               char error[512];
                sprintf(error, "unknown uBarSelectionType in xDampingLayer, available types are:\n        1 : unsteady mapped ubar from database\n        2 : unsteady interpolated uBar from database");
                fatalErrorInFunction("ABLInitialize",  error);
            }

            // allocate memory for uBarInstX and set to zero
            abl->uBarInstX = (Cmpnts **)malloc( sizeof(Cmpnts *) * my );

            for(j=0; j<my; j++)
            {
                abl->uBarInstX[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * mx );

                for(i=0; i<mx; i++)
                {
                    abl->uBarInstX[j][i].x = 0.0;
                    abl->uBarInstX[j][i].y = 0.0;
                    abl->uBarInstX[j][i].z = 0.0;
                }
            }

            // allocate memory for tBarInstX and set to zero
            abl->tBarInstX = (PetscReal **)malloc( sizeof(PetscReal *) * my );

            for(j=0; j<my; j++)
            {
                abl->tBarInstX[j] = (PetscReal *)malloc( sizeof(PetscReal) * mx );

                for(i=0; i<mx; i++)
                {
                    abl->tBarInstX[j][i] = 0.0;
                }
            }

            // allocate memory for nProcsKLine and compute processor influence
            abl->nProcsKLine = (PetscInt **)malloc(sizeof(PetscInt *) * my);

            // local vector
            std::vector<std::vector<PetscInt>> lnProcsKLine(my);

            for(j=0; j<my; j++)
            {
                abl->nProcsKLine[j] = (PetscInt *)malloc( sizeof(PetscInt) * mx );
                lnProcsKLine[j].resize(mx);

                for(i=0; i<mx; i++)
                {
                    abl->nProcsKLine[j][i] = 0;
                    lnProcsKLine[j][i] = 0;
                }
            }

            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    lnProcsKLine[j][i] = 1;
                }
            }

            // scatter to all processors
            for(j=0; j<my; j++)
            {
                // scatter on all processors
                MPI_Allreduce(&lnProcsKLine[j][0], &abl->nProcsKLine[j][0], mx, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);

                // clear memory
                std::vector<PetscInt> ().swap(lnProcsKLine[j]);
            }
        }
        if(abl->xFringeUBarSelectionType == 3)
        {
            concurrentPrecursorInitialize(abl);
        }
    }

    // read the side force region properties
    if(mesh->access->flags->isSideForceActive)
    {
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "xStartSideF",   &(abl->xStartSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "xEndSideF",     &(abl->xEndSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "zStartSideF",   &(abl->zStartSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "zEndSideF",     &(abl->zEndSideF));
    }

    // source terms are computed
    if(abl->controllerType == "write")
    {
        readDictDouble("ABLProperties.dat", "relaxPI",          &(abl->relax));
        readDictDouble("ABLProperties.dat", "alphaPI",          &(abl->alpha));
        readDictDouble("ABLProperties.dat", "timeWindowPI",     &(abl->timeWindow));
    }
    // source terms are read or averaged from available database
    else if(abl->controllerType=="read" || abl->controllerType=="average")
    {
        // use the vector class to append data in the vector (simpler)
        std::vector<std::vector<PetscReal>> preCompSourcesTmp;

        // word by word read
        char word[256];

        // buffer for read data
        PetscReal buffer;

        // time counter
        PetscInt ntimes;

        // file is located into the main folder
        std::ifstream indata;
        indata.open("inflowDatabase/momentumSource");

        if(!indata)
        {
           char error[512];
            sprintf(error, "cannot open file inflowDatabase/momentumSource\n");
            fatalErrorInFunction("ABLInitialize",  error);
        }
        else
        {
            std::string tmpStr;

            // read lines and get number of saved times
            for (ntimes = 0; std::getline(indata, tmpStr); ntimes++);

            // first line is header
            ntimes--;

            // save the number of times
            abl->nSourceTimes = ntimes;
            abl->currentCloseIdx = 0;

            // go back on top of file
            indata.close();
            indata.open("inflowDatabase/momentumSource");

            // skip header line
            std::getline(indata, tmpStr);

            // resize the source table
            preCompSourcesTmp.resize(ntimes);

            for(PetscInt t=0; t<ntimes; t++)
            {
                // read along the line: time | sourceX | sourceY | sourceZ
                for(PetscInt i=0; i<4; i++)
                {
                    indata >> word;
                    std::sscanf(word, "%lf", &buffer);

                    preCompSourcesTmp[t].push_back(buffer);
                }

            }

            indata.close();
        }

        // now store the source data into preCompSources and free the temporary variable
        PetscMalloc(sizeof(PetscReal) * ntimes, &(abl->preCompSources));
        for(PetscInt t=0; t<ntimes; t++)
        {
            PetscMalloc(sizeof(PetscReal) * 4, &(abl->preCompSources[t]));
        }

        for(PetscInt t=0; t<ntimes; t++)
        {
            for(PetscInt i=0; i<4; i++)
            {
               abl->preCompSources[t][i] =  preCompSourcesTmp[t][i];
            }
        }

        // clean the temporary variables
        for(PetscInt t=0; t<ntimes; t++)
        {
            std::vector<PetscReal> ().swap(preCompSourcesTmp[t]);
        }

        // if controllerType is average then average the preCompSources
        if(abl->controllerType=="average")
        {
            readDictDouble("ABLProperties.dat", "controllerAvgStartTime", &(abl->sourceAvgStartTime));

            // check that average start time is in the list
            if(abl->sourceAvgStartTime < abl->preCompSources[0][0])
            {
               char error[512];
                sprintf(error, "parameter 'controllerAvgStartTime' is lower than the first available time");
                fatalErrorInFunction("ABLInitialize",  error);
            }
            // check that more than 100 s of history are used to average
            else if(abl->sourceAvgStartTime > abl->preCompSources[ntimes-1][0] - 100.00)
            {
               char error[512];
                sprintf(error, "Lower 'controllerAvgStartTime' parameter. Average is too poor (less than 100 s)");
                fatalErrorInFunction("ABLInitialize",  error);
            }
            else
            {
                // warn if less then 1000 s of history are used to average
                if(abl->sourceAvgStartTime > abl->preCompSources[ntimes-1][0] - 1000.00)
                {
                   char error[512];
                    sprintf(error, "Lower 'controllerAvgStartTime' parameter. Average could be too poor (less than 1000 s)");
                    warningInFunction("ABLInitialize",  error);
                }

                // initialize average counter to zero
                PetscInt  nAvgSources = 0;
                PetscInt  timeOldSet  = 0;
                PetscReal timeOld;

                abl->avgTimeStep = 0.0;

                // average source terms
                for(PetscInt t=0; t<ntimes; t++)
                {
                    if(abl->preCompSources[t][0] > abl->sourceAvgStartTime)
                    {
                        if(!timeOldSet)
                        {
                            timeOld    = abl->preCompSources[t][0];
                            timeOldSet = 1;
                        }
                        
                        abl->cumulatedSource.x += abl->preCompSources[t][1];
                        abl->cumulatedSource.y += abl->preCompSources[t][2];
                        abl->cumulatedSource.z += abl->preCompSources[t][3];
                        abl->avgTimeStep       += (abl->preCompSources[t][0] - timeOld);
                        timeOld                =  abl->preCompSources[t][0];
                        nAvgSources++;
                    }
                }

                // divide by total number of sources used for the average
                abl->cumulatedSource.x = abl->cumulatedSource.x / nAvgSources;
                abl->cumulatedSource.y = abl->cumulatedSource.y / nAvgSources;
                abl->cumulatedSource.z = abl->cumulatedSource.z / nAvgSources;
                abl->avgTimeStep       = abl->avgTimeStep       / nAvgSources;

                PetscPrintf(mesh->MESH_COMM, "average driving sources = (%e %e %e), average time step = %lf\n\n", abl->cumulatedSource.x, abl->cumulatedSource.y, abl->cumulatedSource.z, abl->avgTimeStep);
            }
        }
    }
    else
    {
       char error[512];
        sprintf(error, "unknown controllerType, available types are:\n        1 : write\n        2 : read\n        3 : average");
        fatalErrorInFunction("ABLInitialize",  error);
    }
  }

  return(0);
}
