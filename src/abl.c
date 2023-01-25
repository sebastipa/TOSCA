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

    PetscInt nLevels = my-2;                 // the vertical direction is the j direction in curvilinear coordinates

    Cmpnts        ***coor, ***cent;          // point and cell center coordinates
    PetscReal     ***aj;

    char          dataLoc[256], fileName[500];

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

    DMGetCoordinatesLocal(mesh->da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    if(mesh->meshFileType != "cartesian")
    {
        char error[512];
        sprintf(error, "ABL capabilites only available for cartesian meshes\n");
        fatalErrorInFunction("InitializeABL",  error);
    }

    PetscPrintf(mesh->MESH_COMM, "Reading ABL properties...\n");

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
    readDictInt   ("ABLProperties.dat", "coriolisActive",   &(abl->coriolisActive));
    readDictInt   ("ABLProperties.dat", "controllerActive", &(abl->controllerActive));
    readDictInt   ("ABLProperties.dat", "controllerActiveT",&(abl->controllerActiveT));

    if(abl->coriolisActive)
    {
        readDictDouble("ABLProperties.dat", "fCoriolis", &(abl->fc));
    }

    // find friction velocity based on neutral log law
    abl->uTau = abl->uRef * abl->vkConst / std::log(abl->hRef / abl->hRough);

    // initialize mesh levels (used in fringe region and velocity controller type write)
    PetscPrintf(mesh->MESH_COMM, "   initializing levels\n");
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

        std::vector<PetscReal> ().swap(lLevels);
        std::vector<PetscReal> ().swap(gLevels);
        std::vector<PetscReal> ().swap(lVolumes);
        std::vector<PetscReal> ().swap(gVolumes);
        std::vector<PetscInt>  ().swap(lCells);
        std::vector<PetscInt>  ().swap(gCells);
    }

    if(abl->controllerActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading driving controller properties\n");

        if(abl->controllerActiveT)
        {
            PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->tDes));

            for(l=0; l<nLevels; l++)
            {
                abl->tDes[l] = 0.0;
            }

            // read proportional controller relaxation factor (same as the velocity one)
            readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
        }

        readSubDictDouble("ABLProperties.dat", "controllerProperties", "controllerMaxHeight", &(abl->controllerMaxHeight));
        readSubDictWord  ("ABLProperties.dat", "controllerProperties", "controllerType",   &(abl->controllerType));

        // set cumulated sources to zero. They are needed for the integral part of the
        // controller if controllerType is set to 'write' or to store the average
        // sources if controllerType is set to 'average'
        abl->cumulatedSource.x = 0.0;
        abl->cumulatedSource.y = 0.0;
        abl->cumulatedSource.z = 0.0;

        // source terms are computed
        if(abl->controllerType == "pressure" || abl->controllerType == "geostrophic")
        {
            // read PI controller properties
            readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
            readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaPI",          &(abl->alpha));
            readSubDictDouble("ABLProperties.dat", "controllerProperties", "timeWindowPI",     &(abl->timeWindow));

            PetscPrintf(mesh->MESH_COMM, "   -> controller type: %s\n", abl->controllerType.c_str());

            // calculating levels interpolation weights at reference height
            {
                PetscMalloc(sizeof(PetscInt)  * 2,       &(abl->closestLabels));
                PetscMalloc(sizeof(PetscReal) * 2,       &(abl->levelWeights));
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

                PetscPrintf(mesh->MESH_COMM, "   -> l1 = %lf, l2 = %lf, hRef = %lf\n", abl->cellLevels[abl->closestLabels[0]-1],abl->cellLevels[abl->closestLabels[1]-1], abl->hRef);
                PetscPrintf(mesh->MESH_COMM, "   -> sum of weights = %lf, w1 = %lf, w2 = %lf\n", abl->levelWeights[0]+abl->levelWeights[1], abl->levelWeights[0], abl->levelWeights[1]);

                std::vector<PetscReal> ().swap(absLevelDelta);

                if(abl->controllerType == "pressure")
                {
                    // read if geostrophic damping is active
                    readSubDictInt("ABLProperties.dat", "controllerProperties", "geostrophicDamping", &(abl->geostrophicDampingActive));

                    if(abl->geostrophicDampingActive)
                    {
                        if(!abl->coriolisActive)
                        {
                            readDictDouble("ABLProperties.dat", "fCoriolis", &(abl->fc));
                        }

                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingAlpha",      &(abl->geoDampAlpha));
                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingStartTime",  &(abl->geoDampStart));
                        readSubDictDouble("ABLProperties.dat", "controllerProperties", "geoDampingTimeWindow", &(abl->geoDampWindow));

                        abl->geoDampH     = abl->hInv + 0.5 * abl->dInv;
                        abl->geoDampDelta = abl->dInv;
                        abl->geoDampC     = 2.0*(2.0*abl->fc);
                        abl->geoDampUBar  = nSetZero();
						abl->geoDampAvgS  = nSetZero();
						abl->geoDampAvgDT = mesh->access->clock->dt;

						// allocate memory for filtered geostrophic velocity
                        PetscMalloc(sizeof(Cmpnts)*nLevels, &(abl->geoDampU));

						for(j=0; j<nLevels; j++)
						{
							abl->geoDampU[j] = nSetZero();
						}

                        // see if must read the average
                        std::stringstream stream;
                        stream << std::fixed << std::setprecision(mesh->access->clock->timePrecision) << mesh->access->clock->startTime;
                        word location = "./fields/" + mesh->meshName + "/" + stream.str();
                        word fileName = location + "/geostrophicDampingInfo";

                        FILE *fp=fopen(fileName.c_str(), "r");

                        if(fp==NULL)
                        {
                            // if start time > 0 should find the file
                            if(mesh->access->clock->startTime != 0.0)
                            {
                                char error[512];
                                sprintf(error, "cannot open file %s\n", fileName.c_str());
                                fatalErrorInFunction("ABLInitialize",  error);
                            }
                        }
                        else
                        {
                            fclose(fp);
                            readDictVector(fileName.c_str(), "filteredS",   &(abl->geoDampAvgS));
							readDictDouble(fileName.c_str(), "filteredDT",  &(abl->geoDampAvgDT));
							abl->geoDampUBar = nScale(1.0 / (2.0 * abl->fc * abl->geoDampAvgDT), nSetFromComponents(abl->geoDampAvgS.y, abl->geoDampAvgS.x, 0.0));
                            PetscPrintf(mesh->MESH_COMM, "   -> reading filtered geostrophic wind: Ug = (%.3lf, %.3lf, 0.000)\n",abl->geoDampUBar.x, abl->geoDampUBar.y);
                        }
                    }
                }
            }

            // calculating geostrophic wind and interpolation weight at geostrophic wind
            if(abl->controllerType == "geostrophic")
            {
                // read geosptrophic height
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "hGeo",     &(abl->hGeo));
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaGeo", &(abl->geoAngle));

                // initial parameters (should be correct at ABL convergence)
                abl->geoAngle = abl->geoAngle*M_PI/180;
                abl->hubAngle = 0.0;
                abl->omegaBar = 0.0;

                // compute geostrophic speed
                abl->uGeoBar  = nSetFromComponents(NieuwstadtGeostrophicWind(abl), 0.0, 0.0);

                // rotate according to initial angle
                Cmpnts uGeoBarTmp = nSetZero();
                uGeoBarTmp.x = std::cos(abl->geoAngle) * abl->uGeoBar.x - std::sin(abl->geoAngle) * abl->uGeoBar.y;
                uGeoBarTmp.y = std::sin(abl->geoAngle) * abl->uGeoBar.x + std::cos(abl->geoAngle) * abl->uGeoBar.y;
                mSet(abl->uGeoBar, uGeoBarTmp);

                // printf information
                PetscPrintf(mesh->MESH_COMM, "   -> Ug = (%f, %f, %f) m/s, UgMag = %f m/s\n", abl->uGeoBar.x, abl->uGeoBar.y, abl->uGeoBar.z, nMag(abl->uGeoBar));

                // compute levels at geostrophic height
                PetscMalloc(sizeof(PetscInt)  * 2,       &(abl->closestLabelsGeo));
                PetscMalloc(sizeof(PetscReal) * 2,       &(abl->levelWeightsGeo));
                std::vector<PetscReal> absLevelDelta(nLevels);

                for(l=0; l<nLevels; l++)
                {
                    absLevelDelta[l] = std::fabs(abl->cellLevels[l] - abl->hGeo);
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
                    abl->closestLabelsGeo[errI] = minLabel + 1;
                }

                abl->levelWeightsGeo[0] = (abl->cellLevels[abl->closestLabelsGeo[1]-1]-abl->hGeo) / (abl->cellLevels[abl->closestLabelsGeo[1]-1] - abl->cellLevels[abl->closestLabelsGeo[0]-1]);
                abl->levelWeightsGeo[1] = (abl->hGeo-abl->cellLevels[abl->closestLabelsGeo[0]-1]) / (abl->cellLevels[abl->closestLabelsGeo[1]-1] - abl->cellLevels[abl->closestLabelsGeo[0]-1]);

                PetscPrintf(mesh->MESH_COMM, "   -> l1 = %lf, l2 = %lf, hGeo = %lf\n", abl->cellLevels[abl->closestLabelsGeo[0]-1],abl->cellLevels[abl->closestLabelsGeo[1]-1], abl->hGeo);
                PetscPrintf(mesh->MESH_COMM, "   -> sum of weights = %lf, w1 = %lf, w2 = %lf\n", abl->levelWeightsGeo[0]+abl->levelWeightsGeo[1], abl->levelWeightsGeo[0], abl->levelWeightsGeo[1]);

                std::vector<PetscReal> ().swap(absLevelDelta);
            }
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
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "controllerAvgStartTime", &(abl->sourceAvgStartTime));

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

                    PetscPrintf(mesh->MESH_COMM, "   -> average driving sources = (%e %e %e), average time step = %lf\n", abl->cumulatedSource.x, abl->cumulatedSource.y, abl->cumulatedSource.z, abl->avgTimeStep);
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

    // read the recycling fringe region properties
    if(mesh->access->flags->isXDampingActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading x-damping properties\n");

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

            // calculating levels interpolation weights
            PetscMalloc(sizeof(PetscInt)  * 2,       &(abl->closestLabelsFringe));
            PetscMalloc(sizeof(PetscReal) * 2,       &(abl->levelWeightsFringe));
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
                abl->closestLabelsFringe[errI] = minLabel + 1;
            }

            abl->levelWeightsFringe[0] = (abl->cellLevels[abl->closestLabelsFringe[1]-1]-abl->hRef) / (abl->cellLevels[abl->closestLabelsFringe[1]-1] - abl->cellLevels[abl->closestLabelsFringe[0]-1]);
            abl->levelWeightsFringe[1] = (abl->hRef-abl->cellLevels[abl->closestLabelsFringe[0]-1]) / (abl->cellLevels[abl->closestLabelsFringe[1]-1] - abl->cellLevels[abl->closestLabelsFringe[0]-1]);

            PetscPrintf(mesh->MESH_COMM, "   -> l1 = %lf, l2 = %lf, hRef = %lf\n", abl->cellLevels[abl->closestLabelsFringe[0]-1],abl->cellLevels[abl->closestLabelsFringe[1]-1], abl->hRef);
            PetscPrintf(mesh->MESH_COMM, "   -> sum of weights = %lf, w1 = %lf, w2 = %lf\n", abl->levelWeightsFringe[0]+abl->levelWeightsFringe[1], abl->levelWeightsFringe[0], abl->levelWeightsFringe[1]);

            std::vector<PetscReal> ().swap(absLevelDelta);
        }

        // read advection term damping type (0: none, 1: LanzilaoMeyers2022 damping)
        readSubDictInt("ABLProperties.dat", "xDampingProperties", "advectionDampingType", &(abl->advectionDampingType));

        if(abl->advectionDampingType == 0)
        {
            // nothing to do: no advection damping
        }
        else if(abl->advectionDampingType == 1)
        {
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "advDampingStart",            &(abl->advDampingStart));
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "advDampingEnd",              &(abl->advDampingEnd));
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "advDampingDeltaStart",       &(abl->advDampingDeltaStart));
            readSubDictDouble("ABLProperties.dat", "xDampingProperties", "advDampingDeltaEnd",         &(abl->advDampingDeltaEnd));
        }
        else
        {
            char error[512];
            sprintf(error, "unknown advectionDampingType in xDampingLayer, available types are:\n        0 : no advection damping\n        1 : vertical advection damping");
            fatalErrorInFunction("ABLInitialize",  error);
        }

        // read type of fringe region in uBarSelectionType
        // 1. periodized mapped
        // 2. interpolated mapped
        // 3. concurrent precursor (another simulation is started in parallel and we don't use inflow functions)

        readSubDictInt("ABLProperties.dat", "xDampingProperties", "uBarSelectionType", &(abl->xFringeUBarSelectionType));

        if(abl->xFringeUBarSelectionType == 0 || abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2 || abl->xFringeUBarSelectionType == 4)
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

            // steady prescribed ubar with no wind veer
            if(ifPtr->typeU == 0)
            {
                PetscPrintf(mesh->MESH_COMM, "   -> damping type 0: prescribing constant uBar using log law\n");

                // set tBarSelectionType the same as uBarSelectionType
                ifPtr->typeT  = 0;

                // no need to map inflow
                ifPtr->mapT   = 0;
                ifPtr->mapNut = 0;

                // velocity properties
                readSubDictVector("ABLProperties.dat", "xDampingProperties", "directionU", &(ifPtr->Udir));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "hInversion", &(ifPtr->hInv));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "frictionU",  &(ifPtr->uTau));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "kRough",     &(ifPtr->roughness));
                mScale(1.0/nMag(ifPtr->Udir), ifPtr->Udir);

                // temperature properties: already stored in the abl struct
            }
            // steady prescribed uBar with wind veer
            else if(ifPtr->typeU == 4)
            {
                PetscPrintf(mesh->MESH_COMM, "   -> damping type 4: prescribing constant uBar with wind veer\n");

                // set tBarSelectionType the same as uBarSelectionType
                ifPtr->typeT  = 4;

                // no need to map inflow
                ifPtr->mapT   = 0;
                ifPtr->mapNut = 0;

                // velocity properties
                readSubDictVector("ABLProperties.dat", "xDampingProperties", "directionU", &(ifPtr->Udir));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "hInversion", &(ifPtr->hInv));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "hReference", &(ifPtr->Href));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "frictionU",  &(ifPtr->uTau));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "kRough",     &(ifPtr->roughness));
                readSubDictDouble("ABLProperties.dat", "xDampingProperties", "latitude",   &(ifPtr->latitude));
                ifPtr->fc =  2*7.292115e-5*std::sin(ifPtr->latitude * M_PI / 180.0);
                mScale(1.0/nMag(ifPtr->Udir), ifPtr->Udir);

                // temperature properties: already stored in the abl struct
            }
            // unsteady mapped ubar
            else if (ifPtr->typeU == 1)
            {
                PetscPrintf(mesh->MESH_COMM, "   -> damping type 1: unsteady uBar using one-to-one inflow mapping from database\n");

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

                PetscPrintf(mesh->MESH_COMM, "   -> averaging inflow at 10 top cells...");

                // top average to avoid top oscillations
                PetscMalloc(10*sizeof(Cmpnts),    &(ifPtr->uBarAvgTopX));
                PetscMalloc(10*sizeof(PetscReal), &(ifPtr->tBarAvgTopX));

                for(j=0; j<10; j++)
                {
                    mSetValue(ifPtr->uBarAvgTopX[j], 0.0);
                    ifPtr->tBarAvgTopX[j] = 0.0;
                }

                PetscReal heightTop = 0.5 * (abl->cellLevels[ifPtr->n1*ifPtr->prds1-1] + abl->cellLevels[ifPtr->n1*ifPtr->prds1]);
                PetscReal heightBot = 0.5 * (abl->cellLevels[ifPtr->n1*ifPtr->prds1-11] + abl->cellLevels[ifPtr->n1*ifPtr->prds1-10]);

                // height of the inflow database
                ifPtr->avgTopLength = heightTop;

                // width of the merging region
                ifPtr->avgTopDelta  = heightTop - heightBot;

                // 5 top points coordinates
                PetscMalloc(10*sizeof(PetscReal), &(ifPtr->avgTopPointCoords));
                i = 0;
                for (PetscInt ii=ifPtr->n1*ifPtr->prds1-11; ii<ifPtr->n1*ifPtr->prds1-1; ii++)
                {
                    ifPtr->avgTopPointCoords[i] = abl->cellLevels[ii];
                    i++;
                }

                // variable to store inflow function data
                std::vector<std::vector<Cmpnts>>    ucat_plane_tmp(ifPtr->n1wg);
                std::vector<std::vector<PetscReal>> t_plane_tmp(ifPtr->n1wg);

                // set to zero
                for(j=0; j<ifPtr->n1wg; j++)
                {
                    ucat_plane_tmp[j].resize(ifPtr->n2wg);
                    t_plane_tmp[j].resize(ifPtr->n2wg);

                    for(i=0; i<ifPtr->n2wg; i++)
                    {
                        mSetValue(ucat_plane_tmp[j][i], 0.0);
                        t_plane_tmp[j][i] = 0.0;
                    }
                }

                PetscInt ti, nAvg;
                PetscReal ntimes = std::min(ifPtr->inflowT.nInflowTimes, ifPtr->inflowU.nInflowTimes);
                word fname_U, fname_T;

                for(ti=0; ti<ntimes; ti++)
                {
                    fname_U = "inflowDatabase/U/" + getArbitraryTimeName(abl->access->clock, ifPtr->inflowU.inflowTimes[ti]);
                    fname_T = "inflowDatabase/T/" + getArbitraryTimeName(abl->access->clock, ifPtr->inflowT.inflowTimes[ti]);

                    // open the two inflow files and read
                    FILE *fp_U = fopen(fname_U.c_str(), "rb");
                    FILE *fp_T = fopen(fname_T.c_str(), "rb");

                    if(!fp_U || !fp_T)
                    {
                        char error[512];
                        sprintf(error, "cannot open files:\n    %s\n    %s\n", fname_U.c_str(), fname_T.c_str());
                        fatalErrorInFunction("ABLInitialize",  error);
                    }

                    for(j=0; j<ifPtr->n1wg; j++)
                    {
                        PetscInt err1, err2;
                        err1 = fread(&(ucat_plane_tmp[j][0]), sizeof(Cmpnts), ifPtr->n2wg, fp_U);
                        err2 = fread(&(t_plane_tmp[j][0]), sizeof(PetscReal), ifPtr->n2wg, fp_T);
                    }

                    fclose(fp_U);
                    fclose(fp_T);

                    // now average the top 10 cells (exclude ghosts)
                    PetscInt jAvg = 0;
                    for(j=ifPtr->n1wg-11; j<ifPtr->n1wg-1; j++)
                    {
                        for(i=1; i<ifPtr->n2; i++)
                        {
                            mSum(ifPtr->uBarAvgTopX[jAvg], ucat_plane_tmp[j][i]);
                            ifPtr->tBarAvgTopX[jAvg] += t_plane_tmp[j][i];
                        }

                        jAvg++;
                    }
                }

                // number of data summed per level (ntimes times n levels in direction 2)
                nAvg  = (ifPtr->n2-1) * ntimes;

                PetscPrintf(mesh->MESH_COMM, "done\n");

                // now average the top 10 cells (exclude ghosts)
                for(j=0; j<10; j++)
                {
                    mScale(1.0/nAvg, ifPtr->uBarAvgTopX[j]);
                    ifPtr->tBarAvgTopX[j] /= nAvg;
                    PetscPrintf(mesh->MESH_COMM, "   - Uavg[%ld] = (%.2f %.2f %.2f) m/s, thetaAvg = %.2f K\n", j, ifPtr->uBarAvgTopX[j].x, ifPtr->uBarAvgTopX[j].y, ifPtr->uBarAvgTopX[j].z, ifPtr->tBarAvgTopX[j]);
                }

                // temporary (basically forces zero gradient at the top for velocity)
                mSet(ifPtr->uBarAvgTopX[9], ifPtr->uBarAvgTopX[8]);

                for( j=0; j<ifPtr->n1wg; j++)
                {
                    std::vector<Cmpnts>   ().swap(ucat_plane_tmp[j]);
                    std::vector<PetscReal> ().swap(t_plane_tmp[j]);
                }
            }
            // unsteady interpolated ubar
            else if (ifPtr->typeU == 2)
            {
                PetscPrintf(mesh->MESH_COMM, "   -> damping type 2: unsteady uBar using bilinear inflow mapping from database\n");

                // set tBarSelectionType the same as uBarSelectionType
                ifPtr->typeT  = 2;
                ifPtr->mapT   = 1;

                // nut value is dependent on the velocity (no need to map)
                ifPtr->mapNut = 0;

                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Inflow",   &(ifPtr->n1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Inflow",   &(ifPtr->n2));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Periods",  &(ifPtr->prds1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Periods",  &(ifPtr->prds2));

                // read if source mesh is uniform or grading
                readSubDictWord("ABLProperties.dat", "xDampingProperties", "sourceType",  &(ifPtr->sourceType));

                if(ifPtr->sourceType == "uniform")
                {
                    PetscPrintf(mesh->MESH_COMM, "   -> using uniform source mesh type\n");

                    readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth1", &(ifPtr->width1));
                    readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth2", &(ifPtr->width2));

                    // height of the inflow database
                    ifPtr->avgTopLength = ifPtr->n1*ifPtr->prds1*ifPtr->width1;

                    // width of the merging region
                    ifPtr->avgTopDelta  = 10.0 * ifPtr->width1;

                    // 5 top points coordinates
                    PetscMalloc(10*sizeof(PetscReal), &(ifPtr->avgTopPointCoords));

                    for (i=0; i<10; i++)
                    {
                        ifPtr->avgTopPointCoords[i] = ifPtr->avgTopLength - (10-i)*ifPtr->width1 + 0.5*ifPtr->width1;
                    }

                }
                else if(ifPtr->sourceType == "grading")
                {
                    PetscPrintf(mesh->MESH_COMM, "   -> using grading source mesh type\n");

                    std::vector<PetscReal>  Zcart;

                    word pointsFileName     = "./inflowDatabase/inflowMesh.xyz";
                    FILE *meshFileID        = fopen(pointsFileName.c_str(), "r");

                    if(!meshFileID)
                    {
                        char error[512];
                        sprintf(error, "cannot open inflow points file %s\n", pointsFileName.c_str());
                        fatalErrorInFunction("SetInflowWeights", error);
                    }
                    else
                    {
                        // read the source mesh file in .xyz format
                        PetscReal bufferDouble;
                        PetscInt  npx, npy, npz;
                        PetscInt  error = fscanf(meshFileID, "%ld %ld %ld\n", &npx, &npy, &npz);

                        if(ifPtr->n1 != npz - 1 || ifPtr->n2 != npy - 1)
                        {
                            char error[512];
                            sprintf(error, "source mesh given in %s and expected number of cells do not match\n", pointsFileName.c_str());
                            fatalErrorInFunction("SetInflowWeights", error);
                        }

                        Zcart.resize(npz);

                        for (k = 0; k < npx; k++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
                        for (i = 0; i < npy; i++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
                        for (j = 0; j < npz; j++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &Zcart[j]);

                        fclose(meshFileID);

                        // height of the inflow database
                        ifPtr->avgTopLength = Zcart[npz-1] - Zcart[0];

                        // width of the merging region
                        ifPtr->avgTopDelta  = Zcart[npz-1] - Zcart[npz-11];

                        // 5 top points coordinates
                        PetscMalloc(10*sizeof(PetscReal), &(ifPtr->avgTopPointCoords));

                        for (i=0; i<10; i++)
                        {
                            ifPtr->avgTopPointCoords[i] = 0.5 * (Zcart[npz - 11 + i] + Zcart[npz - 11 + i + 1]);
                        }

                        // wipe vectors
                        std::vector<PetscReal> ().swap(Zcart);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "unknown sourceType in uBarSelectionType type 4, available types are\n    1: uniform\n    2: grading\n");
                    fatalErrorInFunction("ABLInitialize",  error);
                }

                // increase n1 and n2 accounting for side ghost cells
                ifPtr->n1wg = ifPtr->n1 + 2;
                ifPtr->n2wg = ifPtr->n2 + 2;

                // build interpolation weights and find periodized inflow cells indices
                SetInflowWeights(mesh, ifPtr);

                // initialize inflow data
                mappedInflowInitialize(ifPtr);

                PetscPrintf(mesh->MESH_COMM, "   -> averaging inflow at 10 top cells...");

                // top average to avoid top oscillations
                PetscMalloc(10*sizeof(Cmpnts),    &(ifPtr->uBarAvgTopX));
                PetscMalloc(10*sizeof(PetscReal), &(ifPtr->tBarAvgTopX));

                for(j=0; j<10; j++)
                {
                    mSetValue(ifPtr->uBarAvgTopX[j], 0.0);
                    ifPtr->tBarAvgTopX[j] = 0.0;
                }

                // variable to store inflow function data
                std::vector<std::vector<Cmpnts>>    ucat_plane_tmp(ifPtr->n1wg);
                std::vector<std::vector<PetscReal>> t_plane_tmp(ifPtr->n1wg);

                // set to zero
                for(j=0; j<ifPtr->n1wg; j++)
                {
                    ucat_plane_tmp[j].resize(ifPtr->n2wg);
                    t_plane_tmp[j].resize(ifPtr->n2wg);

                    for(i=0; i<ifPtr->n2wg; i++)
                    {
                        mSetValue(ucat_plane_tmp[j][i], 0.0);
                        t_plane_tmp[j][i] = 0.0;
                    }
                }

                PetscInt ti, nAvg;
                PetscReal ntimes = std::min(ifPtr->inflowT.nInflowTimes, ifPtr->inflowU.nInflowTimes);
                word fname_U, fname_T;

                for(ti=0; ti<ntimes; ti++)
                {
                    fname_U = "inflowDatabase/U/" + getArbitraryTimeName(abl->access->clock, ifPtr->inflowU.inflowTimes[ti]);
                    fname_T = "inflowDatabase/T/" + getArbitraryTimeName(abl->access->clock, ifPtr->inflowT.inflowTimes[ti]);

                    // open the two inflow files and read
                    FILE *fp_U = fopen(fname_U.c_str(), "rb");
                    FILE *fp_T = fopen(fname_T.c_str(), "rb");

                    if(!fp_U || !fp_T)
                    {
                        char error[512];
                        sprintf(error, "cannot open files:\n    %s\n    %s\n", fname_U.c_str(), fname_T.c_str());
                        fatalErrorInFunction("ABLInitialize",  error);
                    }

                    for(j=0; j<ifPtr->n1wg; j++)
                    {
                        PetscInt err1, err2;
                        err1 = fread(&(ucat_plane_tmp[j][0]), sizeof(Cmpnts), ifPtr->n2wg, fp_U);
                        err2 = fread(&(t_plane_tmp[j][0]), sizeof(PetscReal), ifPtr->n2wg, fp_T);
                    }

                    fclose(fp_U);
                    fclose(fp_T);

                    // now average the top 10 cells (exclude ghosts)
                    PetscInt jAvg = 0;
                    for(j=ifPtr->n1wg-11; j<ifPtr->n1wg-1; j++)
                    {
                        for(i=1; i<ifPtr->n2; i++)
                        {
                            mSum(ifPtr->uBarAvgTopX[jAvg], ucat_plane_tmp[j][i]);
                            ifPtr->tBarAvgTopX[jAvg] += t_plane_tmp[j][i];
                        }

                        jAvg++;
                    }
                }

                // number of data summed per level (ntimes times n levels in direction 2)
                nAvg  = (ifPtr->n2-1) * ntimes;

                PetscPrintf(mesh->MESH_COMM, "done\n");

                // now average the top 5 cells (exclude ghosts)
                for(j=0; j<10; j++)
                {
                    mScale(1.0/nAvg, ifPtr->uBarAvgTopX[j]);
                    ifPtr->tBarAvgTopX[j] /= nAvg;
                    PetscPrintf(mesh->MESH_COMM, "   - Uavg at %.2f m = (%.2f %.2f %.2f) m/s, thetaAvg = %.2f K\n", ifPtr->avgTopPointCoords[j], ifPtr->uBarAvgTopX[j].x, ifPtr->uBarAvgTopX[j].y, ifPtr->uBarAvgTopX[j].z, ifPtr->tBarAvgTopX[j]);
                }

                // temporary (basically forces zero gradient at the top for velocity)
                mSet(ifPtr->uBarAvgTopX[9], ifPtr->uBarAvgTopX[8]);

                for( j=0; j<ifPtr->n1wg; j++)
                {
                    std::vector<Cmpnts>    ().swap(ucat_plane_tmp[j]);
                    std::vector<PetscReal> ().swap(t_plane_tmp[j]);
                }
            }
            else
            {
                char error[512];
                sprintf(error, "unknown uBarSelectionType in xDampingLayer, available types are:\n        0 : steady mapped ubar from inflow function\n        1 : unsteady mapped ubar from database\n        2 : unsteady interpolated uBar from database");
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
        else if(abl->xFringeUBarSelectionType == 3)
        {
            concurrentPrecursorInitialize(abl);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown xFringeUBarSelectionType %ld, available types are:\n\t0: log law\n\t1: unsteadyMappedPeriodizedUniform\n\t2: unsteadyInterpPeriodizedUniform\n\t3: concurrentPrecursor\n\t4: Nieuwstadt model\n", abl->xFringeUBarSelectionType);
            fatalErrorInFunction("ABLInitialize",  error);
        }
    }

    // read the y damping layer properties
    if(mesh->access->flags->isYDampingActive)
    {
        // assume also x-damping is active (cross checked in SetSimulationFlags)

        // read parameters
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingStart",   &(abl->yDampingStart));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingEnd",     &(abl->yDampingEnd));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingDelta",   &(abl->yDampingDelta));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingAlpha",   &(abl->yDampingAlpha));

        if(abl->xFringeUBarSelectionType == 0 || abl->xFringeUBarSelectionType == 4)
        {
            // steady state y-damping

        }
        else if(abl->xFringeUBarSelectionType == 3)
        {
            // time-dependent y-damping

            // initialize arrays
            VecDuplicate(mesh->Cent,  &(abl->uBarInstY));      VecSet(abl->uBarInstY,    0.0);
            VecDuplicate(mesh->Nvert, &(abl->tBarInstY));      VecSet(abl->tBarInstY,    0.0);

            PetscMalloc(mz*sizeof(PetscInt*), &(abl->yFringeInterpIDs));
            PetscMalloc(mz*sizeof(PetscReal*), &(abl->yFringeInterpWeights));

            for(k=0; k<mz; k++)
            {
                PetscMalloc(2*sizeof(PetscInt), &(abl->yFringeInterpIDs[k]));
                PetscMalloc(2*sizeof(PetscReal), &(abl->yFringeInterpWeights[k]));

                for(i=0; i<2; i++)
                {
                    abl->yFringeInterpIDs[k][i]     = 0;
                    abl->yFringeInterpWeights[k][i] = 0.0;
                }
            }

            // create fictitious line
            std::vector<PetscReal> lxCoordinates, gxCoordinates;

            // selects those processors which are in the fringe
            if(abl->access->flags->isConcurrentPrecursorActive)
            {
                precursor_ *precursor = abl->precursor;
                domain_    *domain_p  = precursor->domain;
                mesh_      *mesh_p    = domain_p->mesh;

                DM            da_p   = mesh_p->da, fda_p = mesh_p->fda;
                DMDALocalInfo info_p = mesh->info;
                PetscInt      xs_p   = info_p.xs, xe_p = info_p.xs + info_p.xm;
                PetscInt      ys_p   = info_p.ys, ye_p = info_p.ys + info_p.ym;
                PetscInt      zs_p   = info_p.zs, ze_p = info_p.zs + info_p.zm;
                PetscInt      mx_p   = info_p.mx, my_p = info_p.my, mz_p = info_p.mz;

                Cmpnts        ***cent_p;

                PetscReal     gNperiods = 0, lNperiods = std::floor(mesh->bounds.Lx / mesh_p->bounds.Lx);

                PetscInt      k_p, r, lxs_p, lxe_p, lys_p, lye_p, lzs_p, lze_p;

                lxs_p = xs_p; lxe_p = xe_p; if (xs_p==0) lxs_p = xs_p+1; if (xe_p==mx_p) lxe_p = xe_p-1;
                lys_p = ys_p; lye_p = ye_p; if (ys_p==0) lys_p = ys_p+1; if (ye_p==my_p) lye_p = ye_p-1;
                lzs_p = zs_p; lze_p = ze_p; if (zs_p==0) lzs_p = zs_p+1; if (ze_p==mz_p) lze_p = ze_p-1;

                PetscMPIInt   rank_p, nProcs_p;

                MPI_Comm_rank(mesh_p->MESH_COMM, &rank_p);
                MPI_Comm_size(mesh_p->MESH_COMM, &nProcs_p);

                DMDAVecGetArray(fda_p, mesh_p->lCent,  &cent_p);
                DMDAVecGetArray(fda  , mesh->lCent,    &cent);



                DMDAVecRestoreArray(fda_p, mesh_p->lCent,  &cent_p);
                DMDAVecRestoreArray(fda  , mesh->lCent,    &cent);
            }
        }
        else
        {
            // y-damping not available for x-damping in which an unsteady term is spread in the fringe
            char error[512];
            sprintf(error, "yDampingLayer not available with xDampingLayer where uBarSelectionType is 1 or 2");
            fatalErrorInFunction("ABLInitialize",  error);
        }
    }

    // read the Rayleigh damping layer properties
    if(mesh->access->flags->isZDampingActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading z-damping properties\n");

        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingStart",   &(abl->zDampingStart));
        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingEnd",     &(abl->zDampingEnd));
        readSubDictDouble("ABLProperties.dat", "zDampingProperties", "zDampingAlpha",   &(abl->zDampingAlpha));
        readSubDictInt   ("ABLProperties.dat", "zDampingProperties", "zDampingAlsoXY",  &(abl->zDampingAlsoXY));

        if(abl->zDampingAlsoXY)
        {
            // read the type of XY computation for uBar
            readSubDictInt   ("ABLProperties.dat", "zDampingProperties", "zDampingXYType",  &(abl->zDampingXYType));

            if(abl->zDampingXYType == 2)
            {
                // check that xDampingLayer is activated and concurrent precursor is active
                if(!mesh->access->flags->isXDampingActive)
                {
                    char error[512];
                    sprintf(error, "inconsistent zDampingXYType = 2. Activate xDamping layer");
                    fatalErrorInFunction("ABLInitialize",  error);
                }
                else
                {
                    if(abl->xFringeUBarSelectionType != 3)
                    {
                        char error[512];
                        sprintf(error, "inconsistent zDampingXYType = 2. Set uBarSelectionType = 3 to activate concurrent precursor");
                        fatalErrorInFunction("ABLInitialize",  error);
                    }
                }
            }
            else if (abl->zDampingXYType != 1)
            {
                char error[512];
                sprintf(error, "unknown zDampingXYType. Available types are 1 or 2");
                fatalErrorInFunction("ABLInitialize",  error);
            }

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

    // read the side force region properties
    if(mesh->access->flags->isSideForceActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading side force properties\n");

        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "xStartSideF",   &(abl->xStartSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "xEndSideF",     &(abl->xEndSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "zStartSideF",   &(abl->zStartSideF));
        readSubDictDouble("ABLProperties.dat", "sideForceProperties", "zEndSideF",     &(abl->zEndSideF));
    }

    DMDAVecRestoreArray(fda, Coor, &coor);

    PetscPrintf(mesh->MESH_COMM, "done\n\n");

  }

  return(0);
}

//***************************************************************************************************************//

PetscReal NieuwstadtGeostrophicWind(abl_ *abl)
{
    PetscReal vkConstant = abl->vkConst, C, eta;
    complex   alpha, Const, Wg;

    C   = abl->fc * abl->hInv / (vkConstant * abl->uTau);

    // evaluate velocity above capping height
    eta    = 1.0;
    alpha  = 0.5 + 0.5 * std::sqrt(complex(1.0, 4.0*C));
    Const  = digamma(alpha + 1.0) + digamma(alpha - 1.0) - 2.0 * digamma(1.0);
    Wg     = 1.0 / vkConstant * (std::log(abl->hInv / abl->hRough) - Const);

    PetscReal U    = Wg.real() * abl->uTau,
              V    = Wg.imag() * abl->uTau;

    return(nMag(nSetFromComponents(U,V,0.0)));
}
