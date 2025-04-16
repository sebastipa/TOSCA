//! \file  abl.c
//! \brief Contains atmospheric boundary layer definition functions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"

#if USE_PYTHON
    #include <pybind11/embed.h>
    #include <pybind11/stl.h>
    namespace py = pybind11;
    static py::scoped_interpreter guard{};
#endif

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

    readDictDouble  ("ABLProperties.dat", "hRough",           &(abl->hRough));
    readDictVector2D("ABLProperties.dat", "uRef",             &(abl->uRef));
    readDictDouble  ("ABLProperties.dat", "hRef",             &(abl->hRef));
    readDictDouble  ("ABLProperties.dat", "hInv",             &(abl->hInv));
    readDictDouble  ("ABLProperties.dat", "dInv",             &(abl->dInv));
    readDictDouble  ("ABLProperties.dat", "gInv",             &(abl->gInv));
    readDictDouble  ("ABLProperties.dat", "tRef",             &(abl->tRef));
    readDictDouble  ("ABLProperties.dat", "gTop",             &(abl->gTop));
	readDictDouble  ("ABLProperties.dat", "gABL",             &(abl->gABL));
    readDictDouble  ("ABLProperties.dat", "vkConst",          &(abl->vkConst));
    readDictDouble  ("ABLProperties.dat", "smearT",           &(abl->smear));
    readDictInt     ("ABLProperties.dat", "coriolisActive",   &(abl->coriolisActive));
    readDictInt     ("ABLProperties.dat", "controllerActive", &(abl->controllerActive));
    readDictInt     ("ABLProperties.dat", "controllerActiveT",&(abl->controllerActiveT));
    readDictInt     ("ABLProperties.dat", "perturbations",    &(abl->perturbations));

    if(abl->coriolisActive)
    {
        readDictDouble("ABLProperties.dat", "fCoriolis", &(abl->fc));
    }

    PetscReal uRefMag = PetscSqrtReal(abl->uRef.x*abl->uRef.x + abl->uRef.y*abl->uRef.y);

    // find friction velocity based on neutral log law
    abl->uTau = uRefMag * abl->vkConst / std::log(abl->hRef / abl->hRough);

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
                    lLevels[j-1]  += (cent[k][j][i].z - mesh->grndLevel);
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

        readSubDictDouble("ABLProperties.dat", "controllerProperties", "controllerMaxHeight", &(abl->controllerMaxHeight));
        readSubDictWord  ("ABLProperties.dat", "controllerProperties", "controllerType",   &(abl->controllerType));
        readSubDictWord  ("ABLProperties.dat", "controllerProperties", "controllerAction",   &(abl->controllerAction));

        abl->cumulatedSource.x = 0.0;
        abl->cumulatedSource.y = 0.0;
        abl->cumulatedSource.z = 0.0;

        // source terms are computed
        if(abl->controllerAction == "write")
        {

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
                        readSubDictInt("ABLProperties.dat", "controllerProperties", "mesoScaleInput", &(abl->mesoScaleInputActive));

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
                        
                        if(abl->mesoScaleInputActive)
                        {
                            readMesoScaleVelocityData(abl);
                            findVelocityInterpolationWeightsOnePt(abl);
                        }
                    }
                }

                // calculating geostrophic wind and interpolation weight at geostrophic wind
                if(abl->controllerType == "geostrophic")
                {
                    // read geosptrophic height
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "hGeo",     &(abl->hGeo));
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaGeo", &(abl->geoAngle));

                    // read geostrophic speed
                    PetscReal geoWindMag;
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "uGeoMag", &geoWindMag);

                    // initial parameters (should be correct at ABL convergence)
                    abl->geoAngle = abl->geoAngle*M_PI/180;
                    abl->hubAngle = 0.0;
                    abl->omegaBar = 0.0;

                    // compute geostrophic speed
                    // abl->uGeoBar  = nSetFromComponents(NieuwstadtGeostrophicWind(abl), 0.0, 0.0);
                    abl->uGeoBar  = nSetFromComponents(geoWindMag, 0.0, 0.0);

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

                    // modify uTau to match exact uGeo above inversion (for proper initial condition)
                    abl->uTau = geoWindMag * abl->vkConst / std::log(abl->hInv / abl->hRough);
                }
            }
            else if( (abl->controllerType=="directProfileAssimilation") || (abl->controllerType=="indirectProfileAssimilation") || (abl->controllerType=="waveletProfileAssimilation"))
            {
                // read PI controller properties
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "alphaPI",          &(abl->alpha));
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "timeWindowPI",     &(abl->timeWindow));
                readSubDictInt   ("ABLProperties.dat", "controllerProperties", "avgSources",       &(abl->averageSource));
                readSubDictWord  ("ABLProperties.dat", "controllerProperties", "lowerLayerForcingType",   &(abl->flType));

                /* ABL Height Calculation*/
                //allocate memory for ABL height calculation 
                if(abl->flType == "ablHeight")
                {
                    PetscMalloc(sizeof(PetscReal) * nLevels,       &(abl->avgTotalStress));   
                    PetscMalloc(sizeof(PetscReal) * nLevels,       &(abl->avgHeatFlux));   
    
                    for(j=0; j<nLevels; j++)
                    {
                        abl->avgHeatFlux[j] = 0.0;
                        abl->avgTotalStress[j] = 0.0;
                    } 
    
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "hAverageTime",     &(abl->hAvgTime)); 
    
                    // create height directory, check for height file, if exists get ABL height
                    if
                    (
                        abl->access->clock->it == abl->access->clock->itStart && !rank
                    )
                    {
                        word postProcessFolder = "./postProcessing/";
                        errno = 0;
                        PetscInt dirRes = mkdir(postProcessFolder.c_str(), 0777);
                        if(dirRes != 0 && errno != EEXIST)
                        {
                            char error[512];
                            sprintf(error, "could not create %s directory",postProcessFolder.c_str());
                            fatalErrorInFunction("InitializeABL",  error);
                        }
                        else
                        {
                            word fileName = postProcessFolder + "/boundaryLayerHeight";
    
                            FILE *file = fopen(fileName.c_str(), "r");
    
                            if (file == NULL) 
                            {
                                
                            }
                            else
                            {
                                PetscReal closestAblHeight = 0.0;
                                PetscReal minTimeDiff = 1.0e20;
    
                                // Variables for parsing lines
                                char line[1024];
                                PetscReal time, ablHtAvg, ablHeight;
    
                                // Read the file line by line
                                while (fgets(line, sizeof(line), file)) 
                                {
                                    // Parse the line (assumes tab-separated values)
                                    PetscInt numFields = sscanf(line, "%lf\t%lf\t%lf",
                                                        &time, &ablHtAvg, &ablHeight);
    
                                    // Ensure the line has valid data
                                    if (numFields == 3) {
                                        // Calculate the time difference
                                        PetscReal timeDiff = fabs(time - abl->access->clock->startTime);
    
                                        // Update the closest match if this is closer
                                        if (timeDiff < minTimeDiff) {
                                            minTimeDiff = timeDiff;
                                            closestAblHeight = ablHtAvg;
                                        }
                                    }
                                }
    
                                fclose(file);
    
                                if(closestAblHeight == 0.0)
                                {
                                    abl->ablHt = abl->hRef;
                                }
                                else
                                {
                                    abl->ablHt = closestAblHeight;
                                }
                                
                            }
                        }
                    }    
                }
                else if(abl->flType == "constantHeight")
                {
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "lowestSrcHeight",  &(abl->lowestSrcHt));
                }
                else if(abl->flType == "mesoDataHeight")
                {

                }
                else
                {
                    char error[512];
                    sprintf(error, "unknown lower layer mesoscale forcing type, available types are:\n        1 : constantHeight\n        2 : ablHeight\n        3 : mesoDataHeight\n");
                    fatalErrorInFunction("ABLInitialize",  error);  
                }

                PetscPrintf(mesh->MESH_COMM, "   controller type velocity: %s\n", abl->controllerType.c_str());

                // allocate memory for the cumulated sources at all mesh heights
                PetscMalloc(sizeof(Cmpnts) * nLevels, &(abl->cumulatedSourceHt));

                for(PetscInt i = 0; i < nLevels; i++)
                {
                    abl->cumulatedSourceHt[i] = nSetZero();
                }

                readMesoScaleVelocityData(abl);

                if(abl->averageSource)
                {
                    PetscReal timeScaleMeso = 0.1 * (abl->timeV[abl->numtV - 1] - abl->timeV[0]);

                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "movingAvgWindow",       &(abl->tAvgWindow));

                    if(abl->tAvgWindow > timeScaleMeso)
                    {
                        char error[512];
                        sprintf(error, "moving average window for assimilation source averaging too high. reduce to less than %lf\n", timeScaleMeso);
                        fatalErrorInFunction("ABLInitialize",  error);
                    }

                    PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->avgsrc));

                    for(PetscInt i = 0; i < my-2; i++)
                    {
                        abl->avgsrc[i] = nSetZero();
                    }

                    abl->currAvgtime = 0.0;
                }

                //find the interpolation points and weights for the velocity and temperature fields from the available heights
                findVelocityInterpolationWeights(abl);

                // allocate memory for variables
                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->luMean));
                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->guMean));
                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->srcPA));

                if(abl->controllerType=="indirectProfileAssimilation")
                {
                    readSubDictInt   ("ABLProperties.dat", "controllerProperties", "polynomialOrder",   &(abl->polyOrder));

                    //precompute the polynomial coefficient matrix using least square regression method.
                    computeLSqPolynomialCoefficientMatrix(abl); 
                }

                if(abl->controllerType=="waveletProfileAssimilation")
                {
                    #if USE_PYTHON

                    #else
                        char error[512];
                        sprintf(error, "Wavelet profile assimilation requires python modules enabled. Set USE_PYTHON flag to 1 in makefile\n");
                        fatalErrorInFunction("ABLInitialize",  error);
                    #endif

                    readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletName", &(abl->waveName));
                    readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletTMethod", &(abl->waveTMethod));
                    readSubDictInt   ("ABLProperties.dat", "controllerProperties", "waveletDecompLevel", &(abl->waveLevel));
                    readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletExtn", &(abl->waveExtn));
                    readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletConvolution", &(abl->waveConv));
                    readSubDictInt   ("ABLProperties.dat", "controllerProperties", "waveletBlend", &(abl->waveletBlend));

                    PetscPrintf(PETSC_COMM_WORLD, "   wavelet profile assimilation using %s with wavelet filtering : %s and %ld levels of decomposition\n", abl->waveTMethod.c_str(), abl->waveName.c_str(), abl->waveLevel);

                }
            }
            else if(abl->controllerType == "geostrophicProfileAssimilation")
            {
                readSubDictWord  ("ABLProperties.dat", "controllerProperties", "lowerLayerForcingType",   &(abl->flType));
                if(abl->flType != "ablHeight")
                {
                    char error[512];
                    sprintf(error, "geostrophic profie assimilation requires the lower layer forcing type to be set to ablHeight based\n");
                    fatalErrorInFunction("ABLInitialize",  error);
                }

                PetscPrintf(mesh->MESH_COMM, "   controller type velocity: %s\n", abl->controllerType.c_str());

                //allocate memory for ABL height calculation 
                
                PetscMalloc(sizeof(PetscReal) * nLevels,       &(abl->avgTotalStress));   
                PetscMalloc(sizeof(PetscReal) * nLevels,       &(abl->avgHeatFlux));   

                for(j=0; j<nLevels; j++)
                {
                    abl->avgHeatFlux[j] = 0.0;
                    abl->avgTotalStress[j] = 0.0;
                } 

                readSubDictDouble("ABLProperties.dat", "controllerProperties", "hAverageTime",     &(abl->hAvgTime)); 

                // create height directory, check for height file, if exists get ABL height
                if
                (
                    abl->access->clock->it == abl->access->clock->itStart && !rank
                )
                {
                    word postProcessFolder = "./postProcessing/";
                    errno = 0;
                    PetscInt dirRes = mkdir(postProcessFolder.c_str(), 0777);
                    if(dirRes != 0 && errno != EEXIST)
                    {
                        char error[512];
                        sprintf(error, "could not create %s directory",postProcessFolder.c_str());
                        fatalErrorInFunction("InitializeABL",  error);
                    }
                    else
                    {
                        word fileName = postProcessFolder + "/boundaryLayerHeight";

                        FILE *file = fopen(fileName.c_str(), "r");

                        if (file == NULL) 
                        {

                        }
                        else
                        {
                            PetscReal closestAblHeight = 0.0;
                            PetscReal minTimeDiff = 1.0e20;

                            // Variables for parsing lines
                            char line[1024];
                            PetscReal time, ablHtAvg, ablHeight;

                            // Read the file line by line
                            while (fgets(line, sizeof(line), file)) 
                            {
                                // Parse the line (assumes tab-separated values)
                                PetscInt numFields = sscanf(line, "%lf\t%lf\t%lf",
                                                    &time, &ablHtAvg, &ablHeight);

                                // Ensure the line has valid data
                                if (numFields == 3) {
                                    // Calculate the time difference
                                    PetscReal timeDiff = fabs(time - abl->access->clock->startTime);

                                    // Update the closest match if this is closer
                                    if (timeDiff < minTimeDiff) {
                                        minTimeDiff = timeDiff;
                                        closestAblHeight = ablHtAvg;
                                    }
                                }
                            }

                            fclose(file);

                            abl->ablHt = closestAblHeight;
                        }
                    }
                }

                readMesoScaleVelocityData(abl);

                //find the interpolation points and weights for the velocity from the available heights 
                findVelocityInterpolationWeights(abl);

                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->srcPA));
                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->luMean));
                PetscMalloc(sizeof(Cmpnts) * (my-2), &(abl->guMean));
            }
            else
            {
                char error[512];
                sprintf(error, "unknown controllerType for controller action write, available types are:\n        1 : pressure\n        2 : geostrophic\n        3 : directProfileAssimilation\n");
                fatalErrorInFunction("ABLInitialize",  error);
            }
        }
        else if(abl->controllerAction == "read")
        {
            // source terms are read or averaged from available database
            if(abl->controllerType=="timeSeries" || abl->controllerType=="timeAverageSeries")
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

                    ntimes = 0;
                    for (int t = 0; std::getline(indata, tmpStr); t++)
                    {
                        if (!tmpStr.empty())
                        {
                            ntimes++;
                        }
                    }

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
                PetscMalloc(sizeof(PetscReal *) * ntimes, &(abl->preCompSources));
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

                // if controllerReadType is average then average the preCompSources
                if(abl->controllerType=="timeAverageSeries")
                {
                    readSubDictDouble("ABLProperties.dat", "controllerProperties", "controllerAvgStartTime", &(abl->sourceAvgStartTime));

                    // check that average start time is in the list
                    // if(abl->sourceAvgStartTime < abl->preCompSources[0][0])
                    // {
                    //     char error[512];
                    //     sprintf(error, "parameter 'controllerAvgStartTime' is lower than the first available time");
                    //     fatalErrorInFunction("ABLInitialize",  error);
                    // }
                    // check that more than 100 s of history are used to average
                    
                    if(abl->sourceAvgStartTime > abl->preCompSources[ntimes-1][0] - 100.00)
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

            else if(abl->controllerType=="timeHeightSeries")
            {
                std::vector<std::vector<PetscReal>> preCompSourcesTmp;

                // word by word read
                char word[256];

                // buffer for read data
                PetscReal buffer;

                // time counter
                PetscInt ntimes, numLevels;

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

                    ntimes = 0, numLevels = 0;
                    for (PetscInt t = 0; std::getline(indata, tmpStr); t++)
                    {
                        if (!tmpStr.empty())
                        {
                            ntimes++;
                        }
                    }

                    //first three lines are header
                    ntimes = ntimes-3;

                    // save the number of times
                    abl->nSourceTimes = ntimes;
                    abl->currentCloseIdx = 0;

                    // go back on top of file
                    indata.close();
                    indata.open("inflowDatabase/momentumSource");

                    // skip first line - levels
                    std::getline(indata, tmpStr);

                    //read the number of levels
                    indata >> word;
                    while(strcmp("time", word) !=0 )
                    {
                        numLevels++;
                        indata >> word;
                    }

                    //save number of levels in the original data
                    abl->numhV = numLevels;

                    // go back on top of file
                    indata.close();
                    indata.open("inflowDatabase/momentumSource");

                    // skip first line - levels
                    std::getline(indata, tmpStr);

                    PetscMalloc(sizeof(PetscReal) * numLevels, &(abl->hV));

                    //read each of the levels
                    for (PetscInt j = 0; j<numLevels; j++)
                    {
                        indata >> abl->hV[j];
                    }

                    // go back on top of file
                    indata.close();
                    indata.open("inflowDatabase/momentumSource");

                    //skip the first 3 lines
                    for (PetscInt t = 0; t<3; t++)
                    {
                        std::getline(indata, tmpStr);
                    }

                    // resize the source table
                    preCompSourcesTmp.resize(ntimes);

                    //x,y,z velocity per level + time
                    PetscInt wPerLine = numLevels * 3 + 1;

                    for (PetscInt t = 0; t<ntimes; t++)
                    {
                        for (PetscInt j = 0; j<wPerLine; j++)
                        {
                            indata >> word;
                            std::sscanf(word, "%lf", &buffer);

                            preCompSourcesTmp[t].push_back(buffer);
                        }
                    }

                    indata.close();
                }

                PetscMalloc(sizeof(PetscReal **) * ntimes, &(abl->timeHtSources));

                for(PetscInt t=0; t<ntimes; t++)
                {
                    PetscMalloc(sizeof(PetscReal*) * numLevels, &(abl->timeHtSources[t]));
                }

                for(PetscInt t=0; t<ntimes; t++)
                {
                    for(PetscInt j=0; j<numLevels; j++)
                    {
                        PetscMalloc(sizeof(PetscReal) * 4, &(abl->timeHtSources[t][j]));
                    }
                }

                //save the time for each height
                for(PetscInt t=0; t<ntimes; t++)
                {
                    for(PetscInt j=0; j<numLevels; j++)
                    {
                        abl->timeHtSources[t][j][0] = preCompSourcesTmp[t][0];
                    }
                }

                for(PetscInt t=0; t<ntimes; t++)
                {
                    for(PetscInt j=0; j<numLevels; j++)
                    {
                        abl->timeHtSources[t][j][1] = preCompSourcesTmp[t][3*j + 1];
                        abl->timeHtSources[t][j][2] = preCompSourcesTmp[t][3*j + 2];
                        abl->timeHtSources[t][j][3] = preCompSourcesTmp[t][3*j + 3];
                    }
                }

                // clean the temporary variables
                for(PetscInt t=0; t<ntimes; t++)
                {
                    std::vector<PetscReal> ().swap(preCompSourcesTmp[t]);
                }

                if(numLevels != my-2)
                {
                    //find interpolation ids and weights to interpolate from the original precursor mesh source data to the current mesh
                    findTimeHeightSeriesInterpolationWts(abl);
                }
            }
            // source terms are computed in concurrent precursor, they need to be transferred to successor. 
            else if(abl->controllerType=="timeSeriesFromPrecursor")
            {
                PetscMalloc(sizeof(PetscReal *) * 1, &(abl->preCompSources));
                PetscMalloc(sizeof(PetscReal) * 4, &(abl->preCompSources[0]));

                abl->preCompSources[0][0] = abl->access->clock->startTime;
                abl->preCompSources[0][1] = 0.0;
                abl->preCompSources[0][2] = 0.0;
                abl->preCompSources[0][3] = 0.0;  
            }
            else
            {
                char error[512];
                sprintf(error, "unknown controllerType for controller action read, available types are:\n        1 : timeSeries\n        2 : timeAverageSeries\n        3 : timeHeightSeries\n");
                fatalErrorInFunction("ABLInitialize",  error);
            }
        }
        else
        {
                char error[512];
                sprintf(error, "unknown controllerAction, available types are:\n        1 : write\n        2 : read\n");
                fatalErrorInFunction("ABLInitialize",  error);
        }
    }

    if(abl->controllerActiveT)
    {

        if(!(abl->access->flags->isTeqnActive))
        {
            char error[512];
            sprintf(error, "temperature controller cannot be used without potential temperature flag on\n");
            fatalErrorInFunction("ABLInitialize",  error);
        }
        
        PetscMalloc(sizeof(PetscReal) * nLevels, &(abl->tDes));

        for(l=0; l<nLevels; l++)
        {
            abl->tDes[l] = 0.0;
        }

        if(!abl->controllerActive)
        {
            char error[512];
            sprintf(error, "temperature controller is currently set only if velocity controller is active\n");
            fatalErrorInFunction("ABLInitialize",  error);
        }

        readDictWord     ("ABLProperties.dat", "controllerTypeT",    &(abl->controllerTypeT));

        if(abl->controllerTypeT=="indirectProfileAssimilation" || abl->controllerTypeT=="directProfileAssimilation" || abl->controllerTypeT=="waveletProfileAssimilation")
        {  
            // read proportional controller relaxation factor (same as the velocity one)
            readSubDictDouble("ABLProperties.dat", "controllerProperties", "relaxPI",          &(abl->relax));
            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "lowerLayerForcingTypeT",   &(abl->flTypeT));

            PetscPrintf(mesh->MESH_COMM, "   controller type temperature: %s\n", abl->controllerTypeT.c_str());

            if(abl->flTypeT == "constantHeight")
            {
                readSubDictDouble("ABLProperties.dat", "controllerProperties", "lowestSrcHeight",  &(abl->lowestSrcHt));
            }
            else if(abl->flTypeT == "mesoDataHeight")
            {

            }
            else
            {
                char error[512];
                sprintf(error, "unknown lower layer mesoscale forcing type, available types are:\n        1 : constantHeight\n      2 : mesoDataHeight\n");
                fatalErrorInFunction("ABLInitialize",  error);  
            } 

            readMesoScaleTemperatureData(abl);

            findTemperatureInterpolationWeights(abl);
        }

        if(abl->controllerTypeT=="indirectProfileAssimilation")
        {
            readSubDictInt   ("ABLProperties.dat", "controllerProperties", "polynomialOrderT",   &(abl->polyOrderT));

            //precompute the polynomial coefficient matrix using least square regression method.
            computeLSqPolynomialCoefficientMatrixT(abl); 
        }

        if(abl->controllerTypeT=="waveletProfileAssimilation")
        {
            #if USE_PYTHON
            #else
                char error[512];
                sprintf(error, "Wavelet profile assimilation requires python modules enabled. Set USE_PYTHON flag to 1 in makefile\n");
                fatalErrorInFunction("ABLInitialize",  error);
            #endif

            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletName", &(abl->waveName));
            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletTMethod", &(abl->waveTMethod));
            readSubDictInt   ("ABLProperties.dat", "controllerProperties", "waveletDecompLevel", &(abl->waveLevel));
            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletExtn", &(abl->waveExtn));
            readSubDictWord  ("ABLProperties.dat", "controllerProperties", "waveletConvolution", &(abl->waveConv));
            readSubDictInt   ("ABLProperties.dat", "controllerProperties", "waveletBlend", &(abl->waveletBlend));

            PetscPrintf(PETSC_COMM_WORLD, "   wavelet profile assimilation using %s with wavelet filtering : %s and %ld levels of decomposition\n", abl->waveTMethod.c_str(), abl->waveName.c_str(), abl->waveLevel);
        }
    }

    if(mesh->meshName == "overset")
    {
        // explicitly set damping to 0 for overset domains - ensure overset mesh dont intersect with damping region
        mesh->access->flags->isXDampingActive = 0;
        mesh->access->flags->isYDampingActive = 0;
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

                // read if source mesh is uniform or grading
                readSubDictWord("ABLProperties.dat", "xDampingProperties", "sourceType",  &(ifPtr->sourceType));

                // read interpolation method
                readSubDictWord("ABLProperties.dat", "xDampingProperties", "interpolation",  &(ifPtr->interpMethod));

                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Inflow",   &(ifPtr->n1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Inflow",   &(ifPtr->n2));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n1Periods",  &(ifPtr->prds1));
                readSubDictInt   ("ABLProperties.dat", "xDampingProperties", "n2Periods",  &(ifPtr->prds2));
                
                // overwrite remaining inputs
                ifPtr->merge1 = 1;
                ifPtr->shift2 = 0;

                if(ifPtr->sourceType == "uniform")
                {
                    readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth1", &(ifPtr->width1));
                    readSubDictDouble("ABLProperties.dat", "xDampingProperties", "cellWidth2", &(ifPtr->width2));
                    
                    // height of the inflow database (it is duplicate consider removing)
                    ifPtr->inflowHeigth = ifPtr->n1*ifPtr->prds1*ifPtr->width1;

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
                        ifPtr->inflowHeigth = Zcart[npz-1] - Zcart[0];

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
                    sprintf(error, "unknown sourceType in uBarSelectionType type 2, available types are\n    1: uniform\n    2: grading\n");
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
        PetscPrintf(mesh->MESH_COMM, "   reading y-damping properties\n");

        // read parameters
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingStart",   &(abl->yDampingStart));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingEnd",     &(abl->yDampingEnd));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingDelta",   &(abl->yDampingDelta));
        readSubDictDouble("ABLProperties.dat", "yDampingProperties", "yDampingAlpha",   &(abl->yDampingAlpha));
        readSubDictInt   ("ABLProperties.dat", "yDampingProperties", "yDampingNumPeriods",   &(abl->yDampingNumPeriods));

        if(abl->xFringeUBarSelectionType == 3)
        {
            //check if the streamwise length is a multiple of the number of periodization
            PetscReal xDampingLength = std::abs(abl->xDampingEnd - abl->xDampingStart);
            PetscReal error = mesh->bounds.Lx - xDampingLength*abl->yDampingNumPeriods;

            if(std::abs(error) > 1e-7)
            {
                char error[512];
                sprintf(error, "error in the yDampingNumPeriods, the sucessor domain length needs to be = xFringe length x yDampingNumPeriods to apply lateral damping \n");
                fatalErrorInFunction("ABLInitialize",  error);
            }

            if(abl->yDampingStart >= abl->yDampingEnd)
            {
                char error[512];
                sprintf(error, "yDampingStart defined greater than yDampingEnd. Correct it in ABLProperties.dat\n");
                fatalErrorInFunction("ABLInitialize",  error);
            }

            if(abl->yDampingNumPeriods == 1)
            {
                char error[512];
                sprintf(error, "sucessor domain length = xfringe length. Are you only simulating a fringe region?, if not increase the sucessor domain size.\n");
                warningInFunction("ABLInitialize",  error);
            }

            //x damping end must coincide with a mesh coordinate in the x axis, or will lead to issues when periodization of the cells in x axis
            PetscInt lEndcheckFlag = 0, gEndCheckFlag = 0;
            PetscInt lbegincheckFlag = 0, gbeginCheckFlag = 0;
            for (k=zs; k<lze; k++)
            {
                for (j=ys; j<lye; j++)
                {
                    for (i=xs; i<lxe; i++)
                    {
                        if(coor[k][j][i].x == abl->xDampingEnd)
                        {
                            lEndcheckFlag = 1;
                        }

                        if(coor[k][j][i].x == abl->xDampingStart)
                        {
                            lbegincheckFlag = 1;
                        }
                    }
                }
            }
            MPI_Allreduce(&lEndcheckFlag, &gEndCheckFlag, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lbegincheckFlag, &gbeginCheckFlag, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

            if(gEndCheckFlag == 0)
            {
                char error[512];
                sprintf(error, "x fringe region domain end does not coincide with a mesh co-ordinate in the x direction. This will lead to incomplete cells when periodizing along x\n");
                warningInFunction("ABLInitialize",  error);
            }

            if(gbeginCheckFlag == 0)
            {
                char error[512];
                sprintf(error, "x fringe region domain start does not coincide with a mesh co-ordinate in the x direction. This will lead to incomplete cells when periodizing along x\n");
                warningInFunction("ABLInitialize",  error);
            }

            //find the processors controlling the source damping region - this is intersection domain of x and y damping regions
            initializeYDampingMapping(abl);

            // initialize arrays
            VecDuplicate(mesh->lCent,  &(abl->uBarInstY));      VecSet(abl->uBarInstY,    0.0);
            VecDuplicate(mesh->lNvert, &(abl->tBarInstY));      VecSet(abl->tBarInstY,    0.0);
        }
        else
        {
            // y-damping not available for x-damping in which an unsteady term is spread in the fringe
            char error[512];
            sprintf(error, "yDampingLayer is only available with xDampingLayer with uBarSelectionType 3");
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

    DMDAVecRestoreArray(fda, Coor, &coor);

    // read the side force region properties
    if(mesh->access->flags->isCanopyActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading canopy properties\n");

        readSubDictDouble("ABLProperties.dat", "canopyProperties", "xStartCanopy",   &(abl->xStartCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "xEndCanopy",     &(abl->xEndCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "yStartCanopy",   &(abl->yStartCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "yEndCanopy",     &(abl->yEndCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "zStartCanopy",   &(abl->zStartCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "zEndCanopy",     &(abl->zEndCanopy));
        readSubDictDouble("ABLProperties.dat", "canopyProperties", "cftCanopy",      &(abl->cftCanopy));
        readSubDictVector("ABLProperties.dat", "canopyProperties", "diskDirCanopy",  &(abl->diskDirCanopy));
        mUnit(abl->diskDirCanopy);
    }

    // read kLeft Rayleigh damping properties
    if(mesh->access->flags->isKLeftRayleighDampingActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading kLeft-damping properties\n");

        readSubDictDouble("ABLProperties.dat", "kLeftDampingProperties", "kLeftPatchDist",      &(abl->kLeftPatchDist));
        readSubDictDouble("ABLProperties.dat", "kLeftDampingProperties", "kLeftDampingAlpha",   &(abl->kLeftDampingAlpha));
        readSubDictVector("ABLProperties.dat", "kLeftDampingProperties", "kLeftDampingUBar",    &(abl->kLeftDampingUBar));
        readSubDictDouble("ABLProperties.dat", "kLeftDampingProperties", "kLeftFilterHeight",   &(abl->kLeftDampingFilterHeight));
        readSubDictDouble("ABLProperties.dat", "kLeftDampingProperties", "kLeftFilterWidth",    &(abl->kLeftDampingFilterWidth));
    }

    // read kRight Rayleigh damping properties
    if(mesh->access->flags->isKRightRayleighDampingActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading kRigh-damping properties\n");

        readSubDictDouble("ABLProperties.dat", "kRighDampingProperties", "kRightPatchDist",      &(abl->kRightPatchDist));
        readSubDictDouble("ABLProperties.dat", "kRighDampingProperties", "kRightDampingAlpha",   &(abl->kRightDampingAlpha));
        readSubDictVector("ABLProperties.dat", "kRighDampingProperties", "kRightDampingUBar",    &(abl->kRightDampingUBar));
        readSubDictDouble("ABLProperties.dat", "kRighDampingProperties", "kRightFilterHeight",   &(abl->kRightDampingFilterHeight));
        readSubDictDouble("ABLProperties.dat", "kRighDampingProperties", "kRightFilterWidth",    &(abl->kRightDampingFilterWidth));
    }

    // read advection damping properties
    if(mesh->access->flags->isAdvectionDampingActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading advection damping properties\n");

        readSubDictDouble("ABLProperties.dat", "advectionDampingProperties", "advDampingStart",            &(abl->advDampingStart));
        readSubDictDouble("ABLProperties.dat", "advectionDampingProperties", "advDampingEnd",              &(abl->advDampingEnd));
        readSubDictDouble("ABLProperties.dat", "advectionDampingProperties", "advDampingDeltaStart",       &(abl->advDampingDeltaStart));
        readSubDictDouble("ABLProperties.dat", "advectionDampingProperties", "advDampingDeltaEnd",         &(abl->advDampingDeltaEnd));
    }

    if(mesh->access->flags->isAdvectionDampingYActive)
    {
        PetscPrintf(mesh->MESH_COMM, "   reading advection damping Y properties\n");

        readSubDictDouble("ABLProperties.dat", "advectionDampingYProperties", "advDampingStart",            &(abl->advDampingYStart));
        readSubDictDouble("ABLProperties.dat", "advectionDampingYProperties", "advDampingEnd",              &(abl->advDampingYEnd));
        readSubDictDouble("ABLProperties.dat", "advectionDampingYProperties", "advDampingDeltaStart",       &(abl->advDampingYDeltaStart));
        readSubDictDouble("ABLProperties.dat", "advectionDampingYProperties", "advDampingDeltaEnd",         &(abl->advDampingYDeltaEnd));
    }

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

//***************************************************************************************************************//

PetscErrorCode readMesoScaleVelocityData(abl_ *abl)
{
    word variableName;
    word fileName;
    word filePath = "inflowDatabase/mesoscaleData";

    PetscInt  dim1, dim2;

    PetscInt  numtV, numhV;

    std::vector<word> fileList;
    PetscInt nFiles = 0;

    getFileList(filePath.c_str(), fileList, nFiles);

    if(nFiles == 0)
    {
        char error[512];
        sprintf(error, "no file found in folder %s\n", filePath.c_str());
        fatalErrorInFunction("readMesoScaleVelocityData",  error);
    }

    PetscPrintf(PETSC_COMM_WORLD, "   reading mesoscale velocity data\n");

    for(PetscInt i=0; i<nFiles; i++)
    {
        fileName = filePath + "/" + fileList[i];

        std::ifstream indata;
        indata.open(fileName.c_str());

        if(!indata)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("readMesoScaleVelocityData",  error);
        }
        else
        {
            indata >> dim1 >> dim2;

            if(fileList[i] == "velTime")
            {
                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2;
                }
                else
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(abl->timeV));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> abl->timeV[j];
                }

                numtV = arrSize;

                abl->numtV = numtV;
            }

            if(fileList[i] == "velHeight")
            {

                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2;
                }
                else
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(abl->hV));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> abl->hV[j];
                }

                numhV = arrSize;

                abl->numhV = numhV;
            }

            indata.close();
        }
    }

    PetscInt memAllocatedV = 0;

    for(PetscInt i=0; i<nFiles; i++)
    {
        fileName = filePath + "/" + fileList[i];

        std::ifstream indata;
        indata.open(fileName.c_str());

        if(!indata)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("readMesoScaleVelocityData",  error);
        }
        else
        {
            indata >> dim1 >> dim2;

            if(fileList[i] == "uVel" || fileList[i] == "vVel")
            {
                if(memAllocatedV == 0)
                {
                    if(numhV == dim1 && numtV == dim2)
                    {
                        PetscMalloc(sizeof(Cmpnts *) * dim1, &(abl->uMeso));

                        for(PetscInt j=0; j<dim1; j++)
                        {
                            PetscMalloc(sizeof(Cmpnts) * dim2, &(abl->uMeso[j]));
                        }
                    }
                    else if(numhV == dim2 && numtV == dim1)
                    {
                        // in order to save in the array as variable[dist][time] irrespective of how it was written
                        PetscMalloc(sizeof(Cmpnts *) * dim2, &(abl->uMeso));

                        for(PetscInt j=0; j<dim2; j++)
                        {
                            PetscMalloc(sizeof(Cmpnts) * dim1, &(abl->uMeso[j]));
                        }
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "number of elements in the height array %ld or time array %ld and the mesoscale velocity field array [%ld x %ld] do not match. recheck file write\n", numhV, numtV, dim1, dim2);
                        fatalErrorInFunction("readMesoScaleVelocityData",  error);
                    }

                    memAllocatedV = 1;
                }

                //explicitly set the vertical velocity component to zero as it is not read currently
                if(numhV == dim1)
                {
                    // save the array
                    for(PetscInt j=0; j<dim1; j++)
                    {
                        for(PetscInt k=0; k<dim2; k++)
                        {
                            abl->uMeso[j][k].z = 0.;
                        }
                    }
                }
                else if(numhV == dim2)
                {
                    for(PetscInt j=0; j<dim2; j++)
                    {
                        for(PetscInt k=0; k<dim1; k++)
                        {
                            abl->uMeso[j][k].z = 0.;
                        }
                    }
                }

            }

            if(fileList[i] == "uVel")
            {
                PetscReal **tempArray;

                PetscMalloc(sizeof(PetscReal *) * dim1, &(tempArray));

                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscMalloc(sizeof(PetscReal) * dim2, &(tempArray[j]));
                }

                for(PetscInt j=0; j<dim1; j++)
                {
                    for(PetscInt k=0; k<dim2; k++)
                    {
                        indata >> tempArray[j][k];
                    }
                }

                if(numhV == dim1)
                {
                    // save the array
                    for(PetscInt j=0; j<dim1; j++)
                    {
                        for(PetscInt k=0; k<dim2; k++)
                        {
                            abl->uMeso[j][k].x = tempArray[j][k];
                        }
                    }
                }
                else if(numhV == dim2)
                {
                    for(PetscInt j=0; j<dim2; j++)
                    {
                        for(PetscInt k=0; k<dim1; k++)
                        {
                            abl->uMeso[j][k].x = tempArray[k][j];
                        }
                    }
                }

                //delete temp array
                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscFree(tempArray[j]);
                }

                PetscFree(tempArray); 
            }

            if(fileList[i] == "vVel")
            {
                PetscReal **tempArray;

                PetscMalloc(sizeof(PetscReal *) * dim1, &(tempArray));

                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscMalloc(sizeof(PetscReal) * dim2, &(tempArray[j]));
                }

                for(PetscInt j=0; j<dim1; j++)
                {
                    for(PetscInt k=0; k<dim2; k++)
                    {
                        indata >> tempArray[j][k];
                    }
                }

                if(numhV == dim1)
                {
                    // save the array
                    for(PetscInt j=0; j<dim1; j++)
                    {
                        for(PetscInt k=0; k<dim2; k++)
                        {
                            abl->uMeso[j][k].y = tempArray[j][k];
                        }
                    }
                }
                else if(numhV == dim2)
                {
                    for(PetscInt j=0; j<dim2; j++)
                    {
                        for(PetscInt k=0; k<dim1; k++)
                        {
                            abl->uMeso[j][k].y = tempArray[k][j];
                        }
                    }
                }

                //delete temp array
                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscFree(tempArray[j]);
                }

                PetscFree(tempArray); 
            }

            indata.close();
        }
    }
    return (0);
}

//***************************************************************************************************************//
PetscErrorCode readMesoScaleTemperatureData(abl_ *abl)
{
    word variableName;
    word fileName;
    word filePath = "inflowDatabase/mesoscaleData";

    PetscInt  dim1, dim2;

    PetscInt  numtT, numhT;

    std::vector<word> fileList;
    PetscInt nFiles = 0;

    getFileList(filePath.c_str(), fileList, nFiles);

    if(nFiles == 0)
    {
        char error[512];
        sprintf(error, "no file found in folder %s\n", filePath.c_str());
        fatalErrorInFunction("readMesoScaleTemperatureData",  error);
    }

    PetscPrintf(PETSC_COMM_WORLD, "   reading mesoscale temperature data\n");

    for(PetscInt i=0; i<nFiles; i++)
    {
        fileName = filePath + "/" + fileList[i];

        std::ifstream indata;
        indata.open(fileName.c_str());

        if(!indata)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("readMesoScaleTemperatureData",  error);
        }
        else
        {
            indata >> dim1 >> dim2;

            if(fileList[i] == "thetaTime")
            {

                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2;
                }
                else
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(abl->timeT));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> abl->timeT[j];
                }

                numtT = arrSize;
                abl->numtT = numtT;
            }

            if(fileList[i] == "thetaHeight")
            {

                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2;
                }
                else
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(abl->hT));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> abl->hT[j];
                }

                numhT = arrSize;

                abl->numhT = numhT;
            }

            indata.close();
        }
    }

    for(PetscInt i=0; i<nFiles; i++)
    {
        fileName = filePath + "/" + fileList[i];

        std::ifstream indata;
        indata.open(fileName.c_str());

        if(!indata)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("readMesoScaleTemperatureData",  error);
        }
        else
        {
            indata >> dim1 >> dim2;

            if(fileList[i] == "theta")
            {

                PetscReal **tempArray;

                PetscMalloc(sizeof(PetscReal *) * dim1, &(tempArray));

                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscMalloc(sizeof(PetscReal) * dim2, &(tempArray[j]));
                }

                for(PetscInt j=0; j<dim1; j++)
                {
                    for(PetscInt k=0; k<dim2; k++)
                    {
                        indata >> tempArray[j][k];
                    }
                }

                if(numhT == dim1)
                {
                    PetscMalloc(sizeof(PetscReal *) * dim1, &(abl->tMeso));

                    for(PetscInt j=0; j<dim1; j++)
                    {
                        PetscMalloc(sizeof(PetscReal) * dim2, &(abl->tMeso[j]));
                    }

                    // save the array
                    for(PetscInt j=0; j<dim1; j++)
                    {
                        for(PetscInt k=0; k<dim2; k++)
                        {
                            abl->tMeso[j][k] = tempArray[j][k];
                        }
                    }
                }
                else if(numhT == dim2)
                {
                    PetscMalloc(sizeof(PetscReal *) * dim2, &(abl->tMeso));

                    for(PetscInt j=0; j<dim2; j++)
                    {
                        PetscMalloc(sizeof(PetscReal) * dim1, &(abl->tMeso[j]));
                    }

                    for(PetscInt j=0; j<dim2; j++)
                    {
                        for(PetscInt k=0; k<dim1; k++)
                        {
                            abl->tMeso[j][k] = tempArray[k][j];
                        }
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "number of elements in the height array %ld and the mesoscale temperature field array [%ld x %ld] do not match. recheck file write\n", numhT, dim1, dim2);
                    fatalErrorInFunction("readMesoScaleTemperatureData",  error);
                }

                //delete temp array
                for(PetscInt j=0; j<dim1; j++)
                {
                    PetscFree(tempArray[j]);
                }

                PetscFree(tempArray); 
            }

            indata.close();
        }
    }
    return (0);
}

//***************************************************************************************************************//

PetscErrorCode findTimeHeightSeriesInterpolationWts(abl_ *abl)
{
    mesh_ *mesh = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j, nlevels;

    //local pointers to the abl cellLevels and available mesoscale data levels
    PetscReal *cellLevels, *hVel;

    cellLevels  = abl->cellLevels;
    hVel        = abl->hV;
    nlevels     = my-2;

    // allocate memory for abl->velInterpPts and abl->velInterpWts [nlevels x 2]
    PetscMalloc(sizeof(PetscInt *)  * (nlevels), &(abl->velInterpIdx));
    PetscMalloc(sizeof(PetscReal *) * (nlevels), &(abl->velInterpWts));

    for(j=0; j<nlevels; j++)
    {
        PetscMalloc(sizeof(PetscInt)  * 2, &(abl->velInterpIdx[j]));
        PetscMalloc(sizeof(PetscReal) * 2, &(abl->velInterpWts[j]));
    }

    //loop through the cell levels
    for(i=0; i<nlevels; i++)
    {
        //local variables
        PetscReal wt1, wt2, total;
        PetscInt  idx1, idx2;

        PetscReal currPt = cellLevels[i];

        idx1 = 0;
        idx2 = abl->numhV - 1;

        //loop through the mesoscale data points to find the closest points
        for(j=0; j<abl->numhV; j++)
        {
            if (hVel[j] <= currPt && hVel[j] > hVel[idx1])
            {
                idx1 = j;
            }
            if (hVel[j] >= currPt && hVel[j] < hVel[idx2])
            {
                idx2 = j;
            }
        }

        //ensure index are not same
        if(idx1 == idx2)
        {
            idx2 = idx1+1;
        }

        // Calculate interpolation weights
        wt1 = currPt - hVel[idx1];
        wt2 = hVel[idx2] - currPt;
        total = hVel[idx2] - hVel[idx1];
        abl->velInterpWts[i][0] = wt2 / total;
        abl->velInterpWts[i][1] = wt1 / total;

        if(idx2 > abl->numhV - 1) idx2 = abl->numhV - 1;

        abl->velInterpIdx[i][0] = idx1;
        abl->velInterpIdx[i][1] = idx2;
    }

    // find the lowest and highest cell levels for which data is available
    PetscInt lowestInd = 0, highestInd = nlevels-1;

    i = lowestInd;
    while(cellLevels[i] < hVel[0])
    {
        lowestInd = i;
        i++;
    }
    lowestInd = i;

    i = highestInd;
    while(cellLevels[i] > hVel[abl->numhV - 1])
    {
        highestInd = i;
        i--;
    }
    highestInd = i;

    //lowest and highest index based on the mesoscale input data
    abl->lMesoIndV = lowestInd;
    abl->hMesoIndV = highestInd;

    // for cellLevels outside the bounds of hVel, set the weights based on the lowest and highest available levels
    for(i=0; i<nlevels; i++)
    {
        if(cellLevels[i] < hVel[0])
        {
            abl->velInterpWts[i][0] = abl->velInterpWts[abl->lMesoIndV][0];
            abl->velInterpWts[i][1] = abl->velInterpWts[abl->lMesoIndV][1];

            abl->velInterpIdx[i][0] = abl->velInterpIdx[abl->lMesoIndV][0];
            abl->velInterpIdx[i][1] = abl->velInterpIdx[abl->lMesoIndV][1];
        }

        if(cellLevels[i] > hVel[abl->numhV - 1])
        {
            abl->velInterpWts[i][0] = abl->velInterpWts[abl->hMesoIndV][0];
            abl->velInterpWts[i][1] = abl->velInterpWts[abl->hMesoIndV][1];

            abl->velInterpIdx[i][0] = abl->velInterpIdx[abl->hMesoIndV][0];
            abl->velInterpIdx[i][1] = abl->velInterpIdx[abl->hMesoIndV][1];
        }

        // PetscPrintf(PETSC_COMM_WORLD, "cell %ld, height = %lf, interp ids = %ld, %ld, wts = %lf %lf\n", i, cellLevels[i], abl->velInterpIdx[i][0], abl->velInterpIdx[i][1], abl->velInterpWts[i][0], abl->velInterpWts[i][1]);

    }

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode findVelocityInterpolationWeightsOnePt(abl_ *abl)
{
    mesh_ *mesh = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j;

    //local pointers to the abl cellLevels and available mesoscale data levels
    PetscReal *cellLevels, *hVel;

    cellLevels  = abl->cellLevels;
    hVel        = abl->hV;

    // allocate memory for abl->velInterpPts and abl->velInterpWts 
    PetscMalloc(sizeof(PetscInt *)  * (1), &(abl->velInterpIdx));
    PetscMalloc(sizeof(PetscReal *) * (1), &(abl->velInterpWts));

    PetscMalloc(sizeof(PetscInt)  * 2, &(abl->velInterpIdx[0]));
    PetscMalloc(sizeof(PetscReal) * 2, &(abl->velInterpWts[0]));

    //local variables
    PetscReal wt1, wt2, total;
    PetscInt  idx1, idx2;

    PetscReal currPt = abl->hRef;

    idx1 = 0;
    idx2 = abl->numhV - 1;

    //loop through the mesoscale data points to find the closest points
    for(j=0; j<abl->numhV; j++)
    {
        if (hVel[j] <= currPt && hVel[j] > hVel[idx1])
        {
            idx1 = j;
        }
        if (hVel[j] >= currPt && hVel[j] < hVel[idx2])
        {
            idx2 = j;
        }
    }

    //ensure index are not same
    if(idx1 == idx2)
    {
        idx2 = idx1+1;
    }

    // Calculate interpolation weights
    wt1 = currPt - hVel[idx1];
    wt2 = hVel[idx2] - currPt;
    total = hVel[idx2] - hVel[idx1];
    abl->velInterpWts[0][0] = wt2 / total;
    abl->velInterpWts[0][1] = wt1 / total;

    if(idx2 > abl->numhV - 1) idx2 = abl->numhV - 1;

    abl->velInterpIdx[0][0] = idx1;
    abl->velInterpIdx[0][1] = idx2;

    // find the initial time assimilation index and weights
    PetscReal tim = clock->startTime;

    idx1 = 0;
    idx2 = abl->numtV - 1;

    if(tim < abl->timeV[0] || tim > abl->timeV[abl->numtV - 1])
    {
        char error[512];
        sprintf(error, "initial assimilation time out of bounds of the available mesoscale time data\n");
        fatalErrorInFunction("findVelocityInterpolationWeights",  error);
    }

    for(j=0; j<abl->numtV; j++)
    {
        if (abl->timeV[j] <= tim && abl->timeV[j] > abl->timeV[idx1])
        {
            idx1 = j;
        }
        if (abl->timeV[j] >= tim && abl->timeV[j] < abl->timeV[idx2])
        {
            idx2 = j;
        }
    }

    //ensure index are not same
    if(idx1 == idx2)
    {
        idx2 = idx1+1;
    }

    // Calculate interpolation weights
    wt1 = tim - abl->timeV[idx1];
    wt2 = abl->timeV[idx2] - tim;
    total = abl->timeV[idx2] - abl->timeV[idx1];
    abl->closestTimeWtV = wt2 / total;
    abl->closestTimeIndV = idx1;

    //other wt2 = 1-wt1, idx2 = idx1 + 1

    // PetscPrintf(PETSC_COMM_WORLD, "time interp id = %ld %ld, at time = %lf %lf, for time = %lf, wts = %lf %lf\n", abl->closestTimeIndV, (abl->closestTimeIndV) + 1, abl->timeV[idx1], abl->timeV[idx2], tim, abl->closestTimeWtV, 1 - (abl->closestTimeWtV) );

    return (0);
}
//***************************************************************************************************************//
PetscErrorCode findVelocityInterpolationWeights(abl_ *abl)
{
    mesh_ *mesh = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j, nlevels;

    //local pointers to the abl cellLevels and available mesoscale data levels
    PetscReal *cellLevels, *hVel;

    cellLevels  = abl->cellLevels;
    hVel        = abl->hV;
    nlevels     = my-2;

    // allocate memory for abl->velInterpPts and abl->velInterpWts [nlevels x 2]
    PetscMalloc(sizeof(PetscInt *)  * (nlevels), &(abl->velInterpIdx));
    PetscMalloc(sizeof(PetscReal *) * (nlevels), &(abl->velInterpWts));

    for(j=0; j<nlevels; j++)
    {
        PetscMalloc(sizeof(PetscInt)  * 2, &(abl->velInterpIdx[j]));
        PetscMalloc(sizeof(PetscReal) * 2, &(abl->velInterpWts[j]));
    }

    //loop through the cell levels
    for(i=0; i<nlevels; i++)
    {
        //local variables
        PetscReal wt1, wt2, total;
        PetscInt  idx1, idx2;

        PetscReal currPt = cellLevels[i];

        idx1 = 0;
        idx2 = abl->numhV - 1;

        //loop through the mesoscale data points to find the closest points
        for(j=0; j<abl->numhV; j++)
        {
            if (hVel[j] <= currPt && hVel[j] > hVel[idx1])
            {
                idx1 = j;
            }
            if (hVel[j] >= currPt && hVel[j] < hVel[idx2])
            {
                idx2 = j;
            }
        }

        //ensure index are not same
        if(idx1 == idx2)
        {
            idx2 = idx1+1;
        }

        // Calculate interpolation weights
        wt1 = currPt - hVel[idx1];
        wt2 = hVel[idx2] - currPt;
        total = hVel[idx2] - hVel[idx1];
        abl->velInterpWts[i][0] = wt2 / total;
        abl->velInterpWts[i][1] = wt1 / total;

        if(idx2 > abl->numhV - 1) idx2 = abl->numhV - 1;

        abl->velInterpIdx[i][0] = idx1;
        abl->velInterpIdx[i][1] = idx2;
    }

    if(abl->flType == "ablHeight" || abl->flType == "mesoDataHeight")
    {
        abl->bottomSrcHtV = hVel[0];
    }
    else if(abl->flType == "constantHeight")
    {
        if(abl->lowestSrcHt < hVel[0])
        {
            char error[512];
            sprintf(error, "the lowest interpolation height %lf is less than available mesoscale data height %lf\n", abl->lowestSrcHt, hVel[0]);
            fatalErrorInFunction("findVelocityInterpolationWeights",  error); 
        }
        abl->bottomSrcHtV = abl->lowestSrcHt;
    }

    i = 0;    
    while(cellLevels[i] < abl->bottomSrcHtV)
    {
        i++;
    }
    abl->lowestIndV = i;

    i = nlevels-1;    
    while(cellLevels[i] > hVel[abl->numhV - 1])
    {
        i--;
    }
    abl->highestIndV = i;

    // find the initial time assimilation index and weights
    PetscReal tim = clock->startTime;

    PetscReal wt1, wt2, total;
    PetscInt  idx1, idx2;

    idx1 = 0;
    idx2 = abl->numtV - 1;

    if(tim < abl->timeV[0] || tim > abl->timeV[abl->numtV - 1])
    {
        char error[512];
        sprintf(error, "initial assimilation time out of bounds of the available mesoscale time data\n");
        fatalErrorInFunction("findVelocityInterpolationWeights",  error);
    }

    for(j=0; j<abl->numtV; j++)
    {
        if (abl->timeV[j] <= tim && abl->timeV[j] > abl->timeV[idx1])
        {
            idx1 = j;
        }
        if (abl->timeV[j] >= tim && abl->timeV[j] < abl->timeV[idx2])
        {
            idx2 = j;
        }
    }

    //ensure index are not same
    if(idx1 == idx2)
    {
        idx2 = idx1+1;
    }

    // Calculate interpolation weights
    wt1 = tim - abl->timeV[idx1];
    wt2 = abl->timeV[idx2] - tim;
    total = abl->timeV[idx2] - abl->timeV[idx1];
    abl->closestTimeWtV = wt2 / total;
    abl->closestTimeIndV = idx1;

    //other wt2 = 1-wt1, idx2 = idx1 + 1

    // PetscPrintf(PETSC_COMM_WORLD, "time interp id = %ld %ld, at time = %lf %lf, for time = %lf, wts = %lf %lf\n", abl->closestTimeIndV, (abl->closestTimeIndV) + 1, abl->timeV[idx1], abl->timeV[idx2], tim, abl->closestTimeWtV, 1 - (abl->closestTimeWtV) );

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode findTemperatureInterpolationWeights(abl_ *abl)
{
    mesh_  *mesh  = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j, nlevels;

    //local pointers to the abl cellLevels and available mesoscale data levels
    PetscReal *cellLevels, *hT;

    cellLevels  = abl->cellLevels;
    hT          = abl->hT;
    nlevels     = my-2;

    // allocate memory for abl->velInterpPts and abl->velInterpWts [nlevels x 2]
    PetscMalloc(sizeof(PetscInt *)  * (nlevels), &(abl->tempInterpIdx));
    PetscMalloc(sizeof(PetscReal *) * (nlevels), &(abl->tempInterpWts));

    for(j=0; j<nlevels; j++)
    {
        PetscMalloc(sizeof(PetscInt)  * 2, &(abl->tempInterpIdx[j]));
        PetscMalloc(sizeof(PetscReal) * 2, &(abl->tempInterpWts[j]));
    }

    //loop through the cell levels
    for(i=0; i<nlevels; i++)
    {
        //local variables
        PetscReal wt1, wt2, total;
        PetscInt  idx1, idx2;

        PetscReal currPt = cellLevels[i];

        idx1 = 0;
        idx2 = abl->numhT - 1;

        //loop through the mesoscale data points to find the closest points
        for(j=0; j<abl->numhT; j++)
        {
            if (hT[j] <= currPt && hT[j] > hT[idx1])
            {
                idx1 = j;
            }
            if (hT[j] >= currPt && hT[j] < hT[idx2])
            {
                idx2 = j;
            }

        }

        //ensure index are not same
        if(idx1 == idx2)
        {
            idx2 = idx1+1;
        }

        // Calculate interpolation weights
        wt1 = currPt - hT[idx1];
        wt2 = hT[idx2] - currPt;
        total = hT[idx2] - hT[idx1];
        abl->tempInterpWts[i][0] = wt2 / total;
        abl->tempInterpWts[i][1] = wt1 / total;

        if(idx2 > abl->numhT - 1) idx2 = abl->numhT - 1;

        abl->tempInterpIdx[i][0] = idx1;
        abl->tempInterpIdx[i][1] = idx2;
    }

    if(abl->flTypeT == "mesoDataHeight")
    {
        abl->bottomSrcHtT = hT[0];
    }
    else if(abl->flTypeT == "constantHeight")
    {
        if(abl->lowestSrcHt < hT[0])
        {
            char error[512];
            sprintf(error, "the lowest interpolation height %lf is less than available mesoscale data height %lf\n", abl->lowestSrcHt, abl->hV[0]);
            fatalErrorInFunction("findTemperatureInterpolationWeights",  error); 
        }
        abl->bottomSrcHtT = abl->lowestSrcHt;
    }

    i = 0;    
    while(cellLevels[i] < abl->bottomSrcHtT)
    {
        i++;
    }
    abl->lowestIndT = i;

    i = nlevels-1;    
    while(cellLevels[i] > hT[abl->numhT - 1])
    {
        i--;
    }
    abl->highestIndT = i;

    //find initial time assimilation index and weights
    PetscReal tim = clock->startTime;

    PetscReal wt1, wt2, total;
    PetscInt  idx1, idx2;

    idx1 = 0;
    idx2 = abl->numtT - 1;

    if(tim < abl->timeT[0] || tim > abl->timeT[abl->numtT - 1])
    {
        char error[512];
        sprintf(error, "initial assimilation time out of bounds of the available mesoscale time data\n");
        fatalErrorInFunction("findTemperatureInterpolationWeights",  error);
    }

    for(j=0; j<abl->numtT; j++)
    {
        if (abl->timeT[j] <= tim && abl->timeT[j] > abl->timeT[idx1])
        {
            idx1 = j;
        }
        if (abl->timeT[j] >= tim && abl->timeT[j] < abl->timeT[idx2])
        {
            idx2 = j;
        }
    }

    //ensure index are not same
    if(idx1 == idx2)
    {
        idx2 = idx1+1;
    }

    // Calculate interpolation weights
    wt1 = tim - abl->timeT[idx1];
    wt2 = abl->timeT[idx2] - tim;
    total = abl->timeT[idx2] - abl->timeT[idx1];
    abl->closestTimeWtT = wt2 / total;
    abl->closestTimeIndT = idx1;

    //other wt2 = 1-wt1, idx2 = idx1 + 1

    // PetscPrintf(PETSC_COMM_WORLD, "time interp id = %ld %ld, at time = %lf %lf, for time = %lf, wts = %lf %lf\n", abl->closestTimeIndT, (abl->closestTimeIndT) + 1, abl->timeT[idx1], abl->timeT[idx2], tim, abl->closestTimeWtT, 1 - (abl->closestTimeWtT) );

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode initializeYDampingMapping(abl_ *abl)
{
    mesh_ *mesh = abl->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***cent, ***coor;
    Vec           Coor;
    boundingBox   sourceBound, destProcBound, yFringeBound;
    PetscInt      numSourceProc = 0, sourceBoundFlag = 0;
    PetscInt      numDestBounds = abl->yDampingNumPeriods;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, l, p;

    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    PetscMPIInt   globalSize;
    MPI_Comm_size(mesh->MESH_COMM, &globalSize);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    sourceBound.xmin = abl->xDampingStart;
    sourceBound.xmax = abl->xDampingEnd;
    sourceBound.ymin = abl->yDampingStart;
    sourceBound.ymax = abl->yDampingEnd;
    sourceBound.zmin = mesh->bounds.zmin;
    sourceBound.zmax = mesh->bounds.zmax;

    sourceBound.Lx = std::abs(abl->xDampingEnd - abl->xDampingStart);
    sourceBound.Ly = std::abs(abl->yDampingEnd - abl->yDampingStart);
    sourceBound.Lz = mesh->bounds.Lz;

    yFringeBound.xmin = abl->xDampingEnd;
    yFringeBound.xmax = mesh->bounds.xmax;
    yFringeBound.ymin = abl->yDampingStart;
    yFringeBound.ymax = abl->yDampingEnd;
    yFringeBound.zmin = mesh->bounds.zmin;
    yFringeBound.zmax = mesh->bounds.zmax;

    //find processors that control the source bounding box
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(isInsideBoundingBox(cent[k][j][i], sourceBound))
                {
                    sourceBoundFlag = 1;
                    break;
                }

            }
            if(sourceBoundFlag == 1) break;
        }
        if(sourceBoundFlag == 1) break;
    }

    //total source processors
    MPI_Allreduce(&sourceBoundFlag, &numSourceProc, 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

    abl->numSourceProc = numSourceProc;

    //allocate the memory for the abl variables
    PetscMalloc(numSourceProc * sizeof(PetscInt), &(abl->sourceProcList));
    PetscMalloc(numSourceProc * sizeof(PetscInt), &(abl->srcNumI));
    PetscMalloc(numSourceProc * sizeof(PetscInt), &(abl->srcNumJ));
    PetscMalloc(numSourceProc * sizeof(PetscInt), &(abl->srcNumK));

    PetscMalloc(numSourceProc * sizeof(PetscInt), &(abl->srcCommLocalRank));
    PetscMalloc(numSourceProc * sizeof(cellIds), &(abl->srcMinInd));
    PetscMalloc(numSourceProc * sizeof(cellIds), &(abl->srcMaxInd));
    PetscMalloc(numSourceProc * sizeof(MPI_Comm), &(abl->yDamp_comm));
    PetscMalloc(numSourceProc * sizeof(MPI_Request), &(abl->mapRequest));


    PetscMalloc(sizeof(PetscInt)  * numSourceProc, &(abl->isdestProc));
    PetscMalloc(sizeof(cellIds *)  * numSourceProc, &(abl->destMinInd));
    PetscMalloc(sizeof(cellIds *)  * numSourceProc, &(abl->destMaxInd));

    for(j=0; j<numSourceProc; j++)
    {
        PetscMalloc(sizeof(cellIds)  *  numDestBounds, &(abl->destMinInd[j]));
        PetscMalloc(sizeof(cellIds)  *  numDestBounds, &(abl->destMaxInd[j]));
    }

    //create min and max coordinate for each source processor
    Cmpnts lminCoor[numSourceProc], lmaxCoor[numSourceProc], gminCoor[numSourceProc], gmaxCoor[numSourceProc];
    cellIds lsrcMinInd[numSourceProc], lsrcMaxInd[numSourceProc];

    //create list of source processors
    PetscInt lsrcPrc[numSourceProc], gsrcPrc[numSourceProc];

    for (p = 0; p<numSourceProc; p++)
    {
        lminCoor[p] = nSetZero();
        lmaxCoor[p] = nSetZero();
        gminCoor[p] = nSetZero();
        gmaxCoor[p] = nSetZero();
        lsrcPrc[p]  = 0;
        gsrcPrc[p]  = 0;

        lsrcMinInd[p].i = 0;
        lsrcMinInd[p].j = 0;
        lsrcMinInd[p].k = 0;

        lsrcMaxInd[p].i = 0;
        lsrcMaxInd[p].j = 0;
        lsrcMaxInd[p].k = 0;
    }

    if(sourceBoundFlag == 0)
    {
        sourceBoundFlag = MPI_UNDEFINED;
    }

    //local communicator among the source processes
    MPI_Comm    SOURCE_BOUND_COMM;
    MPI_Comm_split(mesh->MESH_COMM, sourceBoundFlag, rank, &(SOURCE_BOUND_COMM));

    PetscMPIInt lsrcRoot = 0, gsrcRoot = 0;

    if(sourceBoundFlag == 1)
    {
        PetscMPIInt   srcRank;
        MPI_Comm_rank(SOURCE_BOUND_COMM, &srcRank);

        //find the co-ordinates and the range of indices that are within this subdomain
        //co-ordinates will be used to find processors to share with in other subdomains
        //range of indices will be used to share the data

        PetscReal xmin = 1e20, xmax = -1e20, ymin = 1e20, ymax = -1e20, zmin = 1e20, zmax = -1e20;

        //finding processor bounding co-ordinates - based on mesh coordinates
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    if(isInsideBoundingBox(cent[k][j][i], sourceBound))
                    {
                        if(cent[k][j][i].x <= xmin)
                        {
                            xmin = cent[k][j][i].x;
                        }

                        if(cent[k][j][i].y <= ymin)
                        {
                            ymin = cent[k][j][i].y;
                        }

                        if(cent[k][j][i].z <= zmin)
                        {
                            zmin = cent[k][j][i].z;
                        }

                        if(cent[k][j][i].x >= xmax)
                        {
                            xmax = cent[k][j][i].x;
                        }

                        if(cent[k][j][i].y >= ymax)
                        {
                            ymax = cent[k][j][i].y;
                        }

                        if(cent[k][j][i].z >= zmax)
                        {
                            zmax = cent[k][j][i].z;
                        }

                    }
                }
            }
        }

        lminCoor[srcRank] = nSetFromComponents(xmin, ymin, zmin);
        lmaxCoor[srcRank] = nSetFromComponents(xmax, ymax, zmax);
        lsrcPrc [srcRank] = rank;

        xmin = 1e20, xmax = -1e20, ymin = 1e20, ymax = -1e20, zmin = 1e20, zmax = -1e20;
        PetscInt imin, jmin, kmin, imax, jmax, kmax;

        //finding indexes - based on cell center coordinates
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    if(isInsideBoundingBox(cent[k][j][i], sourceBound))
                    {
                        if(cent[k][j][i].x <= xmin)
                        {
                            xmin = cent[k][j][i].x;
                            kmin = k;
                        }

                        if(cent[k][j][i].y <= ymin)
                        {
                            ymin = cent[k][j][i].y;
                            imin = i;
                        }

                        if(cent[k][j][i].z <= zmin)
                        {
                            zmin = cent[k][j][i].z;
                            jmin = j;
                        }

                        if(cent[k][j][i].x >= xmax)
                        {
                            xmax = cent[k][j][i].x;
                            kmax = k;
                        }

                        if(cent[k][j][i].y >= ymax)
                        {
                            ymax = cent[k][j][i].y;
                            imax = i;
                        }

                        if(cent[k][j][i].z >= zmax)
                        {
                            zmax = cent[k][j][i].z;
                            jmax = j;
                        }

                    }
                }
            }
        }

        lsrcMinInd[srcRank].i = imin;
        lsrcMinInd[srcRank].j = jmin;
        lsrcMinInd[srcRank].k = kmin;

        lsrcMaxInd[srcRank].i = imax;
        lsrcMaxInd[srcRank].j = jmax;
        lsrcMaxInd[srcRank].k = kmax;

        MPI_Reduce(&(lminCoor[0]), &(gminCoor[0]), numSourceProc * 3, MPIU_REAL, MPIU_SUM, 0, SOURCE_BOUND_COMM);
        MPI_Reduce(&(lmaxCoor[0]), &(gmaxCoor[0]), numSourceProc * 3, MPIU_REAL, MPIU_SUM,  0, SOURCE_BOUND_COMM);
        MPI_Reduce(&(lsrcPrc[0]), &(gsrcPrc[0]), numSourceProc, MPIU_INT, MPIU_SUM, 0, SOURCE_BOUND_COMM);

        if(srcRank == 0)
        {
            lsrcRoot = rank;
        }
    }

    // get the global rank of the sending processor
    MPI_Allreduce(&(lsrcRoot), &(gsrcRoot), 1, MPIU_INT, MPIU_SUM, mesh->MESH_COMM);

    // broadcast variables local to SOURCE_BOUND_COMM comm to world
    MPI_Bcast(&(gminCoor[0]), numSourceProc * 3, MPIU_REAL, gsrcRoot, mesh->MESH_COMM);
    MPI_Bcast(&(gmaxCoor[0]), numSourceProc * 3, MPIU_REAL, gsrcRoot, mesh->MESH_COMM);
    MPI_Bcast(&(gsrcPrc[0]), numSourceProc, MPIU_INT, gsrcRoot, mesh->MESH_COMM);

    //add buffer to min and max coordinates
    for(p=0;p<numSourceProc;p++)
    {
        gminCoor[p].x = gminCoor[p].x - 1e-7;
        gminCoor[p].y = gminCoor[p].y - 1e-7;
        gminCoor[p].z = gminCoor[p].z - 1e-7;
        gmaxCoor[p].x = gmaxCoor[p].x + 1e-7;
        gmaxCoor[p].y = gmaxCoor[p].y + 1e-7;
        gmaxCoor[p].z = gmaxCoor[p].z + 1e-7;
    }

    PetscPrintf(mesh->MESH_COMM, "   mapping from source processors:\n");

    for(p=0;p<numSourceProc;p++)
    {
        abl->sourceProcList[p] = gsrcPrc[p];
        PetscPrintf(mesh->MESH_COMM, "      processor %ld: %ld\n", p, gsrcPrc[p]);
    }

    PetscInt isdestProc[numSourceProc];

    //initialize local arrays
    for(p=0;p<numSourceProc;p++)
    {
        isdestProc[p] = 0;
    }

    for(p=0; p<numSourceProc; p++)
    {
        //add source processor
        if(rank == gsrcPrc[p])
        {
            isdestProc[p] = 1;
        }

        for(l=0; l<numDestBounds; l++)
        {
            destProcBound.xmin = gminCoor[p].x + (l)*sourceBound.Lx;
            destProcBound.xmax = gmaxCoor[p].x + (l)*sourceBound.Lx;
            destProcBound.ymin = gminCoor[p].y;
            destProcBound.ymax = gmaxCoor[p].y;
            destProcBound.zmin = gminCoor[p].z;
            destProcBound.zmax = gmaxCoor[p].z;

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        if(isInsideBoundingBox(cent[k][j][i], destProcBound))
                        {
                            isdestProc[p] = 1;
                            break;
                        }
                    }
                    if(isdestProc[p] == 1) break;
                }
                if(isdestProc[p] == 1) break;
            }

        }
    }

    // find the indexes in the destination processors where the mapped data needs to be set
    // this is local to a given processor. check for bounds when using.
    for(p=0; p<numSourceProc; p++)
    {
        for(l=0; l<numDestBounds; l++)
        {

            if(isdestProc[p] == 1)
            {
                destProcBound.xmin = gminCoor[p].x + (l)*sourceBound.Lx;
                destProcBound.xmax = gmaxCoor[p].x + (l)*sourceBound.Lx;
                destProcBound.ymin = gminCoor[p].y;
                destProcBound.ymax = gmaxCoor[p].y;
                destProcBound.zmin = gminCoor[p].z;
                destProcBound.zmax = gmaxCoor[p].z;

                PetscReal xmin = 1e20, xmax = -1e20, ymin = 1e20, ymax = -1e20, zmin = 1e20, zmax = -1e20;
                PetscInt imin, jmin, kmin, imax, jmax, kmax;

                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            if(isInsideBoundingBox(cent[k][j][i], destProcBound))
                            {
                                if(cent[k][j][i].x <= xmin)
                                {
                                    xmin = cent[k][j][i].x;
                                    kmin = k;
                                }

                                if(cent[k][j][i].y <= ymin)
                                {
                                    ymin = cent[k][j][i].y;
                                    imin = i;
                                }

                                if(cent[k][j][i].z <= zmin)
                                {
                                    zmin = cent[k][j][i].z;
                                    jmin = j;
                                }

                                if(cent[k][j][i].x >= xmax)
                                {
                                    xmax = cent[k][j][i].x;
                                    kmax = k;
                                }

                                if(cent[k][j][i].y >= ymax)
                                {
                                    ymax = cent[k][j][i].y;
                                    imax = i;
                                }

                                if(cent[k][j][i].z >= zmax)
                                {
                                    zmax = cent[k][j][i].z;
                                    jmax = j;
                                }

                            }
                        }
                    }
                }

                if(xmin == 1e20) kmin = -1;
                if(xmax == -1e20)kmax = -1;
                if(ymin == 1e20) imin = -1;
                if(ymax == -1e20)imax = -1;
                if(zmin == 1e20) jmin = -1;
                if(zmax == -1e20)jmax = -1;

                abl->destMinInd[p][l].i = imin;
                abl->destMinInd[p][l].j = jmin;
                abl->destMinInd[p][l].k = kmin;

                abl->destMaxInd[p][l].i = imax;
                abl->destMaxInd[p][l].j = jmax;
                abl->destMaxInd[p][l].k = kmax;

                // PetscPrintf(PETSC_COMM_SELF, "for processor %ld at periodization = %ld, connected to processor %d: coor = %lf %lf %lf %lf %lf %lf, index = %ld %ld %ld %ld %ld %ld \n", gsrcPrc[p], l, rank, xmin, xmax, ymin, ymax, zmin, zmax, kmin, kmax, imin, imax, jmin, jmax);
            }
        }
    }

    //create communicator between the source and dest processors for each source processor
    for(p=0;p<numSourceProc;p++)
    {
        abl->isdestProc[p] = isdestProc[p];

        if(isdestProc[p] == 0)
        {
            isdestProc[p] = MPI_UNDEFINED;
        }
    }

    for(p=0;p<numSourceProc;p++)
    {
        MPI_Comm_split(mesh->MESH_COMM, isdestProc[p], rank, &(abl->yDamp_comm[p]));
    }

    // make some data global to all the processors of the communicator

    for(p=0; p<numSourceProc; p++)
    {
        if(isdestProc[p] == 1)
        {
            PetscMPIInt lsrcLocal = 0, gsrcLocal = 0;
            PetscInt imin, imax, jmin, jmax, kmin, kmax;
            PetscInt lnumI=0, lnumJ=0, lnumK=0, gnumI = 0, gnumJ = 0, gnumK = 0;

            //making global the src rank for bradcasts later
            if(rank == abl->sourceProcList[p])
            {
                MPI_Comm_rank(abl->yDamp_comm[p], &lsrcLocal);

                imin = lsrcMinInd[p].i;
                imax = lsrcMaxInd[p].i;
                jmin = lsrcMinInd[p].j;
                jmax = lsrcMaxInd[p].j;
                kmin = lsrcMinInd[p].k;
                kmax = lsrcMaxInd[p].k;

                lnumI = imax - imin + 1;
                lnumJ = jmax - jmin + 1;
                lnumK = kmax - kmin + 1;
            }

            MPI_Allreduce(&(lsrcLocal), &(gsrcLocal), 1, MPI_INT, MPIU_SUM, abl->yDamp_comm[p]);
            MPI_Allreduce(&(lnumI), &(gnumI), 1, MPIU_INT, MPIU_SUM, abl->yDamp_comm[p]);
            MPI_Allreduce(&(lnumJ), &(gnumJ), 1, MPIU_INT, MPIU_SUM, abl->yDamp_comm[p]);
            MPI_Allreduce(&(lnumK), &(gnumK), 1, MPIU_INT, MPIU_SUM, abl->yDamp_comm[p]);
            MPI_Allreduce(&(lsrcMinInd[p]), &(abl->srcMinInd[p]), 3, MPIU_INT, MPIU_SUM, abl->yDamp_comm[p]);
            MPI_Allreduce(&(lsrcMaxInd[p]), &(abl->srcMaxInd[p]), 3, MPIU_INT, MPIU_SUM, abl->yDamp_comm[p]);

            abl->srcCommLocalRank[p] = gsrcLocal;
            abl->srcNumI[p] = gnumI;
            abl->srcNumJ[p] = gnumJ;
            abl->srcNumK[p] = gnumK;

            // PetscPrintf(PETSC_COMM_SELF, "for processor %d, local src rank = %d, num mapped elements = %ld %ld %ld \n", rank, abl->srcCommLocalRank[p], abl->srcNumI[p], abl->srcNumJ[p], abl->srcNumK[p]);
        }
    }


    /*allocate memory to the mapped one d array within both the src and destination processors
     if a destination processor gets data from two src processors, it allocates data for both the src processor velocity arrays depending on their
     respective size. */

    PetscMalloc(numSourceProc * sizeof(Cmpnts*), &(abl->velMapped));

    for(p=0; p<numSourceProc; p++)
    {
        PetscInt arraySize = abl->srcNumI[p]*abl->srcNumJ[p]*abl->srcNumK[p];

        if(isdestProc[p] == 1)
        {
            PetscMalloc(arraySize * sizeof(Cmpnts), &(abl->velMapped[p]));
        }
    }

    if(abl->access->flags->isTeqnActive)
    {
        PetscMalloc(numSourceProc * sizeof(PetscReal*), &(abl->tMapped));

        for(p=0; p<numSourceProc; p++)
        {
            PetscInt arraySize = abl->srcNumI[p]*abl->srcNumJ[p]*abl->srcNumK[p];

            if(isdestProc[p] == 1)
            {
                PetscMalloc(arraySize * sizeof(PetscReal), &(abl->tMapped[p]));
            }
        }
    }

    //create a fictitious mesh for interpolation in case the mesh is non-uniform
    setWeightsYDamping(abl);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, Coor, &coor);

    MPI_Comm_free(&SOURCE_BOUND_COMM);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode setWeightsYDamping(abl_ *abl)
{
    mesh_ *mesh = abl->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k, b;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscMPIInt   rank, size;

    // declare arrays
    Cmpnts        ***cent;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    PetscMalloc(sizeof(PetscInt *)  * (mz), &(abl->closestKCell));
    PetscMalloc(sizeof(PetscReal *)  * (mz), &(abl->wtsKCell));

    for(j=0; j<mz; j++)
    {
        PetscMalloc(sizeof(PetscInt)  *  2, &(abl->closestKCell[j]));
        PetscMalloc(sizeof(PetscReal)  *  2, &(abl->wtsKCell[j]));
    }

    //create fictitious k direction mesh based on the source
    std::vector <PetscInt>   lnum(size);                                                // local number of Acceptor cells in each processor
    std::vector <PetscInt>   gnum(size);

    PetscInt jInd = 1, iInd = 1, disp = 0, numkCells = 0;

    for (b=0; b < size; b++){
        lnum[b] = 0;
        gnum[b] = 0;
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(j==jInd && i==iInd)
                {
                    if(cent[k][jInd][iInd].x > abl->xDampingStart && cent[k][jInd][iInd].x < abl->xDampingEnd)
                    {
                        lnum[rank]++;
                    }
                }
            }
        }
    }


    MPI_Allreduce(&lnum[0], &gnum[0], size, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    for (b=0; b<rank; b++){
        disp +=gnum[b];
    }

    for (b=0; b< size; b++){
        numkCells +=gnum[b];
    }

    std::vector<PetscInt> lkSrcCellInd, gkSrcCellInd;
    std::vector<PetscReal> lkSrcCellCoor, gkSrcCellCoor;
    lkSrcCellInd.resize(numkCells);
    gkSrcCellInd.resize(numkCells);
    lkSrcCellCoor.resize(numkCells);
    gkSrcCellCoor.resize(numkCells);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(j==jInd && i==iInd)
                {
                    if(cent[k][jInd][iInd].x > abl->xDampingStart && cent[k][jInd][iInd].x < abl->xDampingEnd)
                    {
                        lkSrcCellInd[disp] = k;
                        lkSrcCellCoor[disp] = cent[k][jInd][iInd].x;
                        disp++;
                    }
                }
            }
        }
    }

    MPI_Allreduce(&lkSrcCellInd[0], &(gkSrcCellInd[0]), numkCells, MPIU_INT, MPI_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lkSrcCellCoor[0], &(gkSrcCellCoor[0]), numkCells, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);

    //find cell centers based on the periodization
    std::vector<PetscReal> kDirectionCenters;
    std::vector<PetscInt> kDirectionOrgIndex;

    PetscInt numK = numkCells*abl->yDampingNumPeriods;

    kDirectionCenters.resize(numK);
    kDirectionOrgIndex.resize(numK);

    for(j=0;j<abl->yDampingNumPeriods;j++)
    {
        for(i=0;i<numkCells;i++)
        {
            kDirectionCenters[j*numkCells + i] = gkSrcCellCoor[i] + j*(abl->xDampingEnd - abl->xDampingStart);
            kDirectionOrgIndex[j*numkCells + i] = gkSrcCellInd[i];

            // PetscPrintf(PETSC_COMM_WORLD, "num = %ld, kcenter = %lf, kOrgIndex = %ld\n", j*numkCells + i, kDirectionCenters[j*numkCells + i], kDirectionOrgIndex[j*numkCells + i]);
        }
    }

    // //now find the closest fictitious k cell centers to each cell
    for (k=lzs; k<lze; k++)
    {
        PetscReal  minDistMag = 1e20;
        PetscInt closestkInd = -1;
        PetscInt dir = 0;

        for(b=0;b<numK;b++)
        {
            PetscReal dist = cent[k][lys][lxs].x - kDirectionCenters[b];
            PetscReal distMag = std::abs(dist);

            if(distMag < minDistMag)
            {
                minDistMag = distMag;
                closestkInd = b;

                if(dist < 0)
                {
                    dir = -1;
                }
                else if(dist > 0)
                {
                    dir = 1;
                }
                else
                {
                    dir = 0;
                }
            }
        }

        if (closestkInd == -1)
        {
            PetscPrintf(PETSC_COMM_SELF, "Warning: No valid closest k-cell found for rank %d, k = %ld\n", rank, k);
            closestkInd = 0; // Assign default to avoid crashes
        }
        
        PetscInt    k1, k2;
        PetscReal   x1, x2;

        if(closestkInd == 0)
        {
            if(dir == -1)
            {
                k1 = closestkInd;
                k2 = closestkInd;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];
            }
            else
            {
                k1 = closestkInd;
                k2 = closestkInd+1;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];
            }
        }
        else if(closestkInd == numK-1)
        {
            if(dir == 1)
            {
                k1 = closestkInd;
                k2 = closestkInd;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];
            }
            else
            {
                k1 = numK-2;
                k2 = numK-1;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];
            }
        }
        else
        {
            if(dir == -1)
            {
                k1 = closestkInd-1;
                k2 = closestkInd;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];

                if(fabs(kDirectionOrgIndex[k1]-kDirectionOrgIndex[k2]) > 1)
                {
                    k1 = closestkInd;
                    k2 = closestkInd;
                    x1 = kDirectionCenters[k1];
                    x2 = kDirectionCenters[k2];
                }
            }
            else
            {
                k1 = closestkInd;
                k2 = closestkInd+1;
                x1 = kDirectionCenters[k1];
                x2 = kDirectionCenters[k2];

                if(fabs(kDirectionOrgIndex[k1]-kDirectionOrgIndex[k2]) > 1)
                {
                    k1 = closestkInd;
                    k2 = closestkInd;
                    x1 = kDirectionCenters[k1];
                    x2 = kDirectionCenters[k2];
                }
            }
        }

        PetscReal coeff = (x2 - x1);

        if(fabs(coeff) < 1e-10)
        {
            abl->wtsKCell[k][0] = 0.5;
            abl->wtsKCell[k][1] = 0.5;
        }
        else 
        {
            abl->wtsKCell[k][0] = (x2 - cent[k][lys][lxs].x)/coeff;
            abl->wtsKCell[k][1] = (cent[k][lys][lxs].x - x1)/coeff;
        }

        abl->closestKCell[k][0] = kDirectionOrgIndex[k1];
        abl->closestKCell[k][1] = kDirectionOrgIndex[k2];

        // if(lxs==1 && lys==1)
        // PetscPrintf(PETSC_COMM_SELF, "rank = %d, cell k= %ld, cell center = %lf, closest cells = %ld %ld, wts  = %lf %lf, closest cell = %ld, cell coord = %lf %lf, dir = %ld\n", rank, k, cent[k][lys][lxs].x, abl->closestKCell[k][0], abl->closestKCell[k][1], abl->wtsKCell[k][0], abl->wtsKCell[k][1], closestkInd, x1, x2, dir);
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    return (0);
}

//***************************************************************************************************************//

PetscErrorCode computeLSqPolynomialCoefficientMatrix(abl_ *abl)
{
    mesh_  *mesh  = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j,k, nlevels, col;

    PetscReal *cellLevels;
    PetscReal **Z, **Zt, **betaMat, **ZtWZ, **ZtWZInv, **W;
    cellLevels  = abl->cellLevels;
    nlevels     = my-2;
    col         = abl->polyOrder + 1;

    //allocate memory for the matrix of cell levels polynomials
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(Z));
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(abl->polyCoeffM));
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(W));

    PetscMalloc(sizeof(PetscReal *)  * (col), &(Zt));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(betaMat));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(ZtWZ));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(ZtWZInv));

    for(j=0; j<nlevels; j++)
    {
        PetscMalloc(sizeof(PetscReal)  * (col), &(Z[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(abl->polyCoeffM[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(W[j]));
    }

    for(j=0; j<col; j++)
    {
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(Zt[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(betaMat[j]));

        PetscMalloc(sizeof(PetscReal)  * (col), &(ZtWZ[j]));
        PetscMalloc(sizeof(PetscReal)  * (col), &(ZtWZInv[j]));
    }

    //set the polynomial matrix Z values
    for(i=0; i<nlevels; i++)
    {
        for(j=0; j<col; j++)
        {
            Z[i][j] = std::pow(cellLevels[i], j);
        }
    }

    //set the polynomial matrix Zt values
    for(i=0; i<col; i++)
    {
        for(j=0; j<nlevels; j++)
        {
            Zt[i][j] = Z[j][i];
        }
    }

    //set uniform weights for now
    abl->wtDist = "uniform";

    if(abl->wtDist == "uniform")
    {
        //define weight matrix
        for(i=0; i<nlevels; i++)
        {
            for(j=0; j<nlevels; j++)
            {
                // if((i == j) && (cellLevels[i] > abl->lowestSrcHt) && (cellLevels[i] < abl->highestSrcHt))
                if(i == j)
                {
                    W[i][j] = 1;
                }
                else
                {
                    W[i][j] = 0;
                }
            }
        }

        //find Z'W
        matMatProduct(Zt, W, betaMat, col, nlevels, nlevels, nlevels);

        //find Z'WZ
        matMatProduct(betaMat, Z, ZtWZ, col, nlevels, nlevels, col);

        //invert the matrix Z'WZ
        if (col == 4)
        {
            inv_4by4(ZtWZ, ZtWZInv, col);
        }
        else if (col == 3)
        {
            inv_3by3(ZtWZ, ZtWZInv, col);
        }
        else if (col == 2)
        {
            PetscReal det = ZtWZ[0][0]*ZtWZ[1][1] - ZtWZ[0][1]*ZtWZ[1][0];
            det = 1/det;

            ZtWZInv[0][0] =  ZtWZ[1][1] * det;
            ZtWZInv[0][1] = -ZtWZ[0][1] * det;
            ZtWZInv[1][0] = -ZtWZ[1][0] * det;
            ZtWZInv[1][1] =  ZtWZ[0][0] * det;

        }
        else
        {
            inverseMatrix(ZtWZ, ZtWZInv, col);
        }

        //find the beta coefficient matrix inv(Z'WZ)Z'
        matMatProduct(ZtWZInv, Zt, betaMat, col, col, col, nlevels);

        //inv(Z'WZ)Z'W
        matMatProduct(betaMat, W, Zt, col, nlevels, nlevels, nlevels);

        // Z inv(Z'WZ)Z'W
        matMatProduct(Z, Zt, abl->polyCoeffM, nlevels, col, col, nlevels);
    }

    //delete temp array
    for(j=0; j<nlevels; j++)
    {
        free(Z[j]);
        free(W[j]);
    }

    for(j=0; j<col; j++)
    {
        free(Zt[j]);
        free(betaMat[j]);
        free(ZtWZ[j]);
        free(ZtWZInv[j]);
    }

    free(W);
    free(Z);
    free(Zt);
    free(betaMat);
    free(ZtWZ);
    free(ZtWZInv);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode computeLSqPolynomialCoefficientMatrixT(abl_ *abl)
{
    mesh_  *mesh  = abl->access->mesh;
    clock_ *clock = abl->access->clock;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt     i,j,k, nlevels, col;

    PetscReal *cellLevels;
    PetscReal **Z, **Zt, **betaMat, **ZtWZ, **ZtWZInv, **W;
    cellLevels  = abl->cellLevels;
    nlevels     = my-2;
    col         = abl->polyOrder + 1;

    //allocate memory for the matrix of cell levels polynomials 
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(Z));
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(abl->polyCoeffT));
    PetscMalloc(sizeof(PetscReal *)  * (nlevels), &(W));

    PetscMalloc(sizeof(PetscReal *)  * (col), &(Zt));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(betaMat));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(ZtWZ));
    PetscMalloc(sizeof(PetscReal *)  * (col), &(ZtWZInv));

    for(j=0; j<nlevels; j++)
    {
        PetscMalloc(sizeof(PetscReal)  * (col), &(Z[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(abl->polyCoeffT[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(W[j]));
    }

    for(j=0; j<col; j++)
    {
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(Zt[j]));
        PetscMalloc(sizeof(PetscReal)  * (nlevels), &(betaMat[j]));

        PetscMalloc(sizeof(PetscReal)  * (col), &(ZtWZ[j]));
        PetscMalloc(sizeof(PetscReal)  * (col), &(ZtWZInv[j]));
    }

    //set the polynomial matrix Z values 
    for(i=0; i<nlevels; i++)
    {
        for(j=0; j<col; j++)
        {
            Z[i][j] = std::pow(cellLevels[i], j);
        }
    }

    //set the polynomial matrix Zt values 
    for(i=0; i<col; i++)
    {
        for(j=0; j<nlevels; j++)
        {
            Zt[i][j] = Z[j][i];
        }
    }

    //set uniform weights for now 
    abl->wtDist = "uniform";

    if(abl->wtDist == "uniform")
    {
        //define weight matrix 
        for(i=0; i<nlevels; i++)
        {
            for(j=0; j<nlevels; j++)
            {
                // if((i == j) && (cellLevels[i] > abl->lowestSrcHt) && (cellLevels[i] < abl->highestSrcHt))
                if(i == j)
                {
                    W[i][j] = 1;
                }
                else
                {
                    W[i][j] = 0;
                }
            }
        }

        //find Z'W
        matMatProduct(Zt, W, betaMat, col, nlevels, nlevels, nlevels);
        
        //find Z'WZ
        matMatProduct(betaMat, Z, ZtWZ, col, nlevels, nlevels, col);

        //invert the matrix Z'WZ
        if (col == 4)
        {
            inv_4by4(ZtWZ, ZtWZInv, col);
        }
        else if (col == 3)
        {
            inv_3by3(ZtWZ, ZtWZInv, col);
        }
        else if (col == 2)
        {
            PetscReal det = ZtWZ[0][0]*ZtWZ[1][1] - ZtWZ[0][1]*ZtWZ[1][0];
            det = 1/det;

            ZtWZInv[0][0] =  ZtWZ[1][1] * det;
            ZtWZInv[0][1] = -ZtWZ[0][1] * det;
            ZtWZInv[1][0] = -ZtWZ[1][0] * det;
            ZtWZInv[1][1] =  ZtWZ[0][0] * det;

        }
        else 
        {
            inverseMatrix(ZtWZ, ZtWZInv, col);
        }

        //find the beta coefficient matrix inv(Z'WZ)Z'
        matMatProduct(ZtWZInv, Zt, betaMat, col, col, col, nlevels);

        //inv(Z'WZ)Z'W
        matMatProduct(betaMat, W, Zt, col, nlevels, nlevels, nlevels);

        // Z inv(Z'WZ)Z'W
        matMatProduct(Z, Zt, abl->polyCoeffT, nlevels, col, col, nlevels);
    }

    //delete temp array 
    for(j=0; j<nlevels; j++)
    {
        free(Z[j]);
        free(W[j]);
    }

    for(j=0; j<col; j++)
    {
        free(Zt[j]);
        free(betaMat[j]);
        free(ZtWZ[j]);
        free(ZtWZInv[j]);
    }

    free(W); 
    free(Z); 
    free(Zt); 
    free(betaMat); 
    free(ZtWZ); 
    free(ZtWZInv);

    return (0);
}
//***************************************************************************************************************//

PetscErrorCode findABLHeight(abl_ *abl)
{
    mesh_  *mesh   = abl->access->mesh;
    ueqn_  *ueqn   = abl->access->ueqn;
    teqn_  *teqn   = abl->access->teqn;
    clock_ *clock  = abl->access->clock;
    les_   *les    = abl->access->les;

    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscMPIInt   rank, nProcs;
    PetscInt      i, j, k, l;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***ucat;                     // cartesian vel. and fluct. part
    PetscReal     ***tmprt, ***nut, ***nvert, ***meshTag;  // potential temp. and fluct. part and turb. visc.
    PetscReal     ***aj, ***ld_t;

    PetscReal     nu = abl->access->constants->nu;
    PetscReal     ablHeight;
    PetscReal     aN, m1, m2, n1, n2;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    MPI_Comm_size(mesh->MESH_COMM, &nProcs);

    DMDAVecGetArray(da,  mesh->lAj,       &aj);
    DMDAVecGetArray(da,  mesh->lNvert,    &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,  &meshTag);
    DMDAVecGetArray(fda, ueqn->lUcat,     &ucat);
    DMDAVecGetArray(fda, mesh->lCsi,      &csi);
    DMDAVecGetArray(fda, mesh->lEta,      &eta);
    DMDAVecGetArray(fda, mesh->lZet,      &zet);
    DMDAVecGetArray(da,  les->lNu_t,      &nut);
    DMDAVecGetArray(da,  les->lDiff_t,    &ld_t);

    DMDAVecGetArray(da,  teqn->lTmprt,    &tmprt);

	PetscInt nLevels = my-2;
    PetscInt stabilityFlag = 0, methodFlag = 0;
    PetscInt zMaxId, zMinId;
    PetscInt RiCrossId, tkeCrossId;
    PetscReal maxTempGrad;
    PetscReal ablHtConv, ablHtStable, ablHtShear, ablHtRiB;

    std::vector<Cmpnts> lVelocity(nLevels);
    std::vector<Cmpnts> gVelocity(nLevels);
    std::vector<PetscReal> lTemperature(nLevels);
    std::vector<PetscReal> gTemperature(nLevels);
    std::vector<PetscReal> luw(nLevels);
    std::vector<PetscReal> guw(nLevels);
    std::vector<PetscReal> lvw(nLevels);
    std::vector<PetscReal> gvw(nLevels);
    std::vector<PetscReal> lTw(nLevels);
    std::vector<PetscReal> gTw(nLevels);
    std::vector<PetscReal> lRxz(nLevels);
    std::vector<PetscReal> gRxz(nLevels);
    std::vector<PetscReal> lRyz(nLevels);
    std::vector<PetscReal> gRyz(nLevels);
    std::vector<PetscReal> lRTw(nLevels);
    std::vector<PetscReal> gRTw(nLevels); 

    PetscReal *RiB;
    PetscMalloc(sizeof(PetscReal *) * nLevels, &(RiB));

    for(l=0; l<nLevels; l++)
    {
        lVelocity[l].x = lVelocity[l].y = lVelocity[l].z = 0.0;
        gVelocity[l].x = gVelocity[l].y = gVelocity[l].z = 0.0;

        luw[l] = 0.0;
        guw[l] = 0.0;
        lvw[l] = 0.0;
        gvw[l] = 0.0;
        lRxz[l] = 0.0;
        gRxz[l] = 0.0;
        lRyz[l] = 0.0;
        gRyz[l] = 0.0;

        lTemperature[l] = 0.0;
        gTemperature[l] = 0.0;
        lTw[l]          = 0.0;
        gTw[l]          = 0.0;
        lRTw[l] = 0.0;
        gRTw[l] = 0.0;
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                lVelocity[j-1].x  += ucat[k][j][i].x / aj[k][j][i];
                lVelocity[j-1].y  += ucat[k][j][i].y / aj[k][j][i];
                lVelocity[j-1].z  += ucat[k][j][i].z / aj[k][j][i];

                lTemperature[j-1] += tmprt[k][j][i]  / aj[k][j][i];
            }
        }
    }

    MPI_Allreduce(&lVelocity[0], &gVelocity[0], 3*nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lTemperature[0], &gTemperature[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    for(l=0; l<nLevels; l++)
    {
        PetscReal totVolPerLevel = abl->totVolPerLevel[l];

        gVelocity[l].x   = gVelocity[l].x / totVolPerLevel;
        gVelocity[l].y   = gVelocity[l].y / totVolPerLevel;
        gVelocity[l].z   = gVelocity[l].z / totVolPerLevel;

        gTemperature[l]   = gTemperature[l] / totVolPerLevel;
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                luw[j-1]  += ((ucat[k][j][i].x - gVelocity[j-1].x) * (ucat[k][j][i].z - gVelocity[j-1].z)) / aj[k][j][i];
                lvw[j-1]  += ((ucat[k][j][i].y - gVelocity[j-1].y) * (ucat[k][j][i].z - gVelocity[j-1].z)) / aj[k][j][i];

                lTw[j-1]  += ((tmprt[k][j][i] - gTemperature[j-1]) * (ucat[k][j][i].z - gVelocity[j-1].z)) / aj[k][j][i];
            }
        }
    }

    MPI_Allreduce(&luw[0], &guw[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lvw[0], &gvw[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    MPI_Allreduce(&lTw[0], &gTw[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    PetscReal volCell;
    Cmpnts    uprimeCell;
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // pre-set base variables for speed
                PetscReal dudc, dvdc, dwdc,
                            dude, dvde, dwde,
                            dudz, dvdz, dwdz;
                PetscReal du_dx, du_dy, du_dz,
                            dv_dx, dv_dy, dv_dz,
                            dw_dx, dw_dy, dw_dz;

                PetscReal dtdc, dtde, dtdz;
                PetscReal dt_dx, dt_dy, dt_dz;

                PetscReal   ajc  = aj[k][j][i];

                PetscReal   csi0 = csi[k][j][i].x,
                            csi1 = csi[k][j][i].y,
                            csi2 = csi[k][j][i].z;
                PetscReal   eta0 = eta[k][j][i].x,
                            eta1 = eta[k][j][i].y,
                            eta2 = eta[k][j][i].z;
                PetscReal   zet0 = zet[k][j][i].x,
                            zet1 = zet[k][j][i].y,
                            zet2 = zet[k][j][i].z;

                volCell     = 1.0 / ajc;

                uprimeCell.x = ucat[k][j][i].x - gVelocity[j-1].x;
                uprimeCell.y = ucat[k][j][i].x - gVelocity[j-1].x;
                uprimeCell.z = ucat[k][j][i].x - gVelocity[j-1].x;

                Compute_du_center
                (
                    mesh,
                    i, j, k, mx, my, mz, ucat, nvert, meshTag, &dudc,
                    &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                );

                Compute_du_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                    zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                    dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                    &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                PetscReal nuEff = nut[k][j][i];
                lRxz[j-1] += ( - nuEff * (dw_dx + du_dz)) * volCell;
                lRyz[j-1] += ( - nuEff * (dv_dz + dw_dy)) * volCell;

                Compute_dscalar_center
                (
                    mesh,
                    i, j, k,
                    mx, my, mz,
                    tmprt, nvert, meshTag,
                    &dtdc, &dtde, &dtdz
                );

                // compute temperature derivative w.r.t. the cartesian coordinates
                Compute_dscalar_dxyz
                (
                    mesh,
                    csi0, csi1, csi2,
                    eta0, eta1, eta2,
                    zet0, zet1, zet2,
                    ajc,
                    dtdc, dtde, dtdz,
                    &dt_dx, &dt_dy, &dt_dz
                );

                PetscReal diffEff = ld_t[k][j][i];
                lRTw[j-1] += ( - diffEff * dt_dz) * volCell;

            }
        }
    }

    MPI_Allreduce(&lRxz[0], &gRxz[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lRyz[0], &gRyz[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    MPI_Allreduce(&lRTw[0], &gRTw[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    aN = 200.0;
    m1 = (aN)  / ((aN) + clock->dt);
    m2 = clock->dt / ((aN) + clock->dt);

    for(l=0; l<nLevels; l++)
    {
        PetscReal totVolPerLevel = abl->totVolPerLevel[l];

        guw[l]   = guw[l] / totVolPerLevel;
        gvw[l]   = gvw[l] / totVolPerLevel;
        gRxz[l]  = gRxz[l] / totVolPerLevel;
        gRyz[l]  = gRyz[l] / totVolPerLevel;

        gTw[l]  = gTw[l] / totVolPerLevel;
        gRTw[l] = gRTw[l] / totVolPerLevel;

        abl->avgTotalStress[l] = m1 * abl->avgTotalStress[l] + m2 * pow((pow(guw[l]+gRxz[l],2.0) + pow(gvw[l]+gRyz[l],2.0)),0.5);

        abl->avgHeatFlux[l]    = m1 * abl->avgHeatFlux[l]    + m2 * (gTw[l] + gRTw[l]);
    }

    //switching condition based on the value of the heat flux close to the surface - second cell above ground
    abl->heatFluxSwitch = abl->avgHeatFlux[0];

    //convective boundary layer - ABL height based on the heat flux minima after the maxima which occurs near surface
    //find the heat flux maximum
    zMaxId = 0, zMinId = -1;
    PetscReal minValue = 1.0e20;  

    for(l=1; l<nLevels; l++)
    {
        if(abl->avgHeatFlux[l] > abl->avgHeatFlux[zMaxId]) zMaxId = l;
    }

    for(PetscInt l = zMaxId + 1; l < nLevels-1; l++)
    {
        if(abl->avgHeatFlux[l] < abl->avgHeatFlux[l-1] && 
        abl->avgHeatFlux[l] < abl->avgHeatFlux[l+1])
        {
            // Track the deepest minima
            if(abl->avgHeatFlux[l] < minValue){
                minValue = abl->avgHeatFlux[l];
                zMinId = l;
            }
        }
    }

    if(zMinId != -1)
    {
        ablHtConv = abl->cellLevels[zMinId];
    }
    else 
    {
        // Use Potential temperature gradient Method
        maxTempGrad = -1e9;
        PetscInt gradId = 0;

        for(l=10; l<nLevels-1; l++) 
        {
            PetscReal grad = (gTemperature[l+1] - gTemperature[l-1])/(abl->cellLevels[l+1]-abl->cellLevels[l-1]);
            if(grad > maxTempGrad) 
            {
                maxTempGrad = grad;
                gradId = l;
            }
        }
        ablHtConv = abl->cellLevels[gradId];   
    }

    // stable and neutral boundary layer - shear stress threshold based 
    PetscReal tkeMax = -1e9;
    PetscInt  tkeMaxId = 0;
    tkeCrossId = -1;

    //find max shear stress in the surface layer
    for(l=0; l<nLevels; l++)
    {   
        if(abl->cellLevels[l] < 300.0)
        {
            if(abl->avgTotalStress[l] > tkeMax) 
            {
                tkeMax = abl->avgTotalStress[l];
                tkeMaxId = l;
            }
        }
    }

    for (l = 0; l < nLevels; l++) 
    {
        if(abl->cellLevels[l] > abl->cellLevels[tkeMaxId])
        {
            if (abl->avgTotalStress[l] <= 0.1*tkeMax) 
            {
                tkeCrossId = l; 
                break;
            }
        }
    }

    ablHtShear = abl->cellLevels[tkeCrossId];

    //based on the bulk richardson number
    PetscReal g = 9.81;
    PetscReal thetaRef = gTemperature[0];
    RiCrossId = -1;
    
    //find RiB
    for(l=0; l<nLevels; l++)
    {
        PetscReal dtheta = gTemperature[l] - thetaRef;
        PetscReal wind = sqrt(gVelocity[l].x*gVelocity[l].x + gVelocity[l].y*gVelocity[l].y + 1e-6);
        RiB[l] = (g/thetaRef) * (dtheta * abl->cellLevels[l]) / (wind * wind);
    }

    //stable when RiB crosses the critical value of 0.25
    for(l=0; l<nLevels; l++)
    {
        if(RiB[l] >= 0.25) 
        {
            RiCrossId = l;
            break;
        }
    }

    ablHtRiB = abl->cellLevels[RiCrossId]; 

    if(ablHtRiB < 100) ablHtRiB = 100;

    if(fabs(ablHtRiB - ablHtShear) < 500)
    {
        ablHtStable = ablHtShear;            
    }
    else 
    {
        ablHtStable = PetscMin(ablHtShear, ablHtRiB);

        if(ablHtStable == ablHtRiB)
            PetscPrintf(PETSC_COMM_WORLD, "RiB based ABL height used.\n");
    }

    if(abl->heatFluxSwitch > 0.02)
    {
        ablHeight = ablHtConv;
        stabilityFlag = 1;
    }
    else 
    {
        ablHeight = ablHtStable;
        stabilityFlag = 2;
    }

    aN = abl->hAvgTime;
    n1 = aN  / (aN + clock->dt);
    n2 = clock->dt / (aN + clock->dt);
    abl->ablHt = n1*abl->ablHt + n2*ablHeight;

    FILE *f;
    dataABL *ablStat = abl->access->acquisition->statisticsABL;
    word          fileName;

    // create the ABL averaging files
    if
    (
        !rank
    )
    {
        fileName = "./postProcessing/" + mesh->meshName + "/boundaryLayerHeight";
        f = fopen(fileName.c_str(), "a");

        if(!f)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("writeAveragingABL",  error);
        }
        else 
        {
            fprintf(f, "%.5lf\t", clock->time);
            fprintf(f, "%.5lf\t", abl->ablHt);
            fprintf(f, "%.5lf\t", ablHeight);
            fprintf(f, "\n");
            fclose(f);
        }
    }
    
    if(stabilityFlag == 1)
    {
            PetscPrintf(PETSC_COMM_WORLD, "Convective Boundary Layer, Using heat flux minima method, Avg ABL height = %lf, Inst. Ht = %lf\n", abl->ablHt, ablHeight);
            PetscPrintf(PETSC_COMM_WORLD, "SurfaceHeatFlux = %lf, Max heat flux = %lf at %lf, heat flux minima = %lf\n", abl->heatFluxSwitch, abl->avgHeatFlux[zMaxId], abl->cellLevels[zMaxId], abl->avgHeatFlux[zMinId]);

    }
    else if(stabilityFlag == 2)
    {
        if(abl->heatFluxSwitch < -0.005)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Stable Boundary Layer, Using turbulent shear stress threshold, Avg ABL height = %lf, Inst. Ht = %lf\n", abl->ablHt, ablHeight);
            PetscPrintf(PETSC_COMM_WORLD, "SurfaceHeatFlux = %lf, Max shear stress = %lf at %lf\n", abl->heatFluxSwitch, tkeMax, abl->cellLevels[tkeMaxId]);
        }
        else 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Neutral Boundary Layer, Using turbulent shear stress threshold, Avg ABL height = %lf, Inst. Ht = %lf\n", abl->ablHt, ablHeight);
            PetscPrintf(PETSC_COMM_WORLD, "SurfaceHeatFlux = %lf, Max shear stress = %lf at %lf\n", abl->heatFluxSwitch, tkeMax, abl->cellLevels[tkeMaxId]);
        }
    }
    

    std::vector<Cmpnts> ().swap(lVelocity);
    std::vector<Cmpnts> ().swap(gVelocity);
    std::vector<PetscReal> ().swap(lTemperature);
    std::vector<PetscReal> ().swap(gTemperature);
    std::vector<PetscReal> ().swap(luw);
    std::vector<PetscReal> ().swap(guw);
    std::vector<PetscReal> ().swap(lvw);
    std::vector<PetscReal> ().swap(gvw);
    std::vector<PetscReal> ().swap(lTw);
    std::vector<PetscReal> ().swap(gTw);
    std::vector<PetscReal> ().swap(lRxz);
    std::vector<PetscReal> ().swap(gRxz);
    std::vector<PetscReal> ().swap(lRyz);
    std::vector<PetscReal> ().swap(gRyz);
    std::vector<PetscReal> ().swap(lRTw);
    std::vector<PetscReal> ().swap(gRTw);

    PetscFree(RiB);

    DMDAVecRestoreArray(da,  mesh->lAj,       &aj);
    DMDAVecRestoreArray(da,  mesh->lNvert,    &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,  &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,     &ucat);
    DMDAVecRestoreArray(fda, mesh->lCsi,      &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,      &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,      &zet);
    DMDAVecRestoreArray(da,  les->lNu_t,      &nut);
    DMDAVecRestoreArray(da,  les->lDiff_t,    &ld_t);

    DMDAVecRestoreArray(da,  teqn->lTmprt,    &tmprt);
    return (0);
}

//***************************************************************************************************************//

PetscErrorCode waveletTransformContinuousVector(abl_ *abl, Cmpnts *srcPAIn, Cmpnts *srcPAOut, PetscInt sigLength)
{
    PetscReal    omega        = abl->omega;         // Central frequency of the wavelet
    PetscReal    sigma        = abl->sigma;         // Width of the wavelet
    PetscInt     kernelRadius = abl->kernelRadius;  // Radius of the wavelet kernel (affects smoothing extent)

    PetscInt    kernelLength  = 2 * kernelRadius + 1;
    PetscInt    paddedLength  = sigLength + 2 * kernelRadius;

    PetscReal   *kernel;
    Cmpnts      *paddedSignal;
    
    PetscMalloc(sizeof(PetscReal) * kernelLength, &kernel);
    PetscMalloc(sizeof(Cmpnts) * paddedLength, &paddedSignal);

    //add padding
    for (PetscInt i = 0; i < kernelRadius; i++) {
        paddedSignal[i] = nSetZero(); 
    }

    for (PetscInt i = 0; i < sigLength; i++) {
        paddedSignal[kernelRadius + i] = nSet(srcPAIn[i]);
    }

    for (PetscInt i = 0; i < kernelRadius; i++) {
        paddedSignal[kernelRadius + sigLength + i] = nSetZero(); 
    }

    // Create the wavelet kernel
    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        PetscReal t = i - kernelRadius;
        kernel[i] = exp(-t * t / (2 * sigma * sigma)) * cos(omega * t);
    }
 
    // Normalize the kernel to ensure energy conservation
    PetscReal kernelSum = 0.0;

    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        kernelSum += fabs(kernel[i]);
    }

    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        kernel[i] /= kernelSum;
    }

    // Perform convolution
    for (PetscInt i = 0; i < sigLength; i++) 
    {
        srcPAOut[i] = nSetZero();

        for (PetscInt j = -kernelRadius; j <= kernelRadius; j++) 
        {
            PetscInt paddedIndex = i + j + kernelRadius;

            srcPAOut[i].x += paddedSignal[paddedIndex].x * kernel[kernelRadius + j];
            srcPAOut[i].y += paddedSignal[paddedIndex].y * kernel[kernelRadius + j];
            srcPAOut[i].z += paddedSignal[paddedIndex].z * kernel[kernelRadius + j];
        }
    }

    PetscFree(kernel);
    PetscFree(paddedSignal);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode waveletTransformContinuousScalar(abl_ *abl, PetscReal *srcPAIn, PetscReal *srcPAOut, PetscInt sigLength)
{

    PetscReal    omega        = abl->omega;         // Central frequency of the wavelet
    PetscReal    sigma        = abl->sigma;         // Width of the wavelet
    PetscInt     kernelRadius = abl->kernelRadius;  // Radius of the wavelet kernel (affects smoothing extent)

    PetscInt    kernelLength  = 2 * kernelRadius + 1;
    PetscInt    paddedLength  = sigLength + 2 * kernelRadius;

    PetscReal   *kernel, *paddedSignal;

    PetscMalloc(sizeof(PetscReal) * kernelLength, &kernel);
    PetscMalloc(sizeof(PetscReal) * paddedLength, &paddedSignal);

    //add padding
    for (PetscInt i = 0; i < kernelRadius; i++) {
        paddedSignal[i] = 0.0; 
    }

    for (PetscInt i = 0; i < sigLength; i++) {
        paddedSignal[kernelRadius + i] = srcPAIn[i];
    }
    for (PetscInt i = 0; i < kernelRadius; i++) {
        paddedSignal[kernelRadius + sigLength + i] = 0.0; 
    }

    // Create the wavelet kernel
    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        PetscReal t = i - kernelRadius;
        kernel[i] = exp(-t * t / (2 * sigma * sigma)) * cos(omega * t);
    }
 
    // Normalize the kernel to ensure energy conservation
    PetscReal kernelSum = 0.0;

    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        kernelSum += fabs(kernel[i]);
    }

    for (PetscInt i = 0; i < kernelLength; i++) 
    {
        kernel[i] /= kernelSum;
    }

    // Perform convolution
    for (PetscInt i = 0; i < sigLength; i++) 
    {
        srcPAOut[i] = 0.0;

        for (PetscInt j = -kernelRadius; j <= kernelRadius; j++) 
        {
            PetscInt paddedIndex = i + j + kernelRadius;
            srcPAOut[i] += paddedSignal[paddedIndex] * kernel[kernelRadius + j];
        }
    }

    PetscFree(kernel);
    PetscFree(paddedSignal);

    return(0);
}

//*********************************************************************************************************************
#if USE_PYTHON
PetscErrorCode pywavedecVector(abl_ *abl, Cmpnts *srcPAIn, Cmpnts *srcPAOut, PetscInt sigLength)
{
    try 
    {
        // Import PyWavelets
        py::module pywt = py::module::import("pywt");

        // Convert C++ arrays to Python lists
        std::vector<double> signalX(sigLength), signalY(sigLength), signalZ(sigLength);
        for (PetscInt i = 0; i < sigLength; i++) 
        {
            signalX[i] = srcPAIn[i].x;
            signalY[i] = srcPAIn[i].y;
            signalZ[i] = srcPAIn[i].z;
        }

        PetscPrintf(PETSC_COMM_WORLD,"vector %s %ld\n", abl->waveName.c_str(), abl->waveLevel);
        PetscInt J = abl->waveLevel; // Decomposition level

        // Perform wavelet decomposition (wavedec)
        py::list coeffsX = pywt.attr("wavedec")(signalX, abl->waveName, py::arg("level") = J, py::arg("mode") = abl->waveExtn);
        py::list coeffsY = pywt.attr("wavedec")(signalY, abl->waveName, py::arg("level") = J, py::arg("mode") = abl->waveExtn);
        py::list coeffsZ = pywt.attr("wavedec")(signalZ, abl->waveName, py::arg("level") = J, py::arg("mode") = abl->waveExtn);

        // Get Approximation Coefficients (First element in the list)
        py::list pApproxCoeffX = coeffsX[0].cast<py::list>();
        py::list pApproxCoeffY = coeffsY[0].cast<py::list>();
        py::list pApproxCoeffZ = coeffsZ[0].cast<py::list>();

        // // Print Approximation Coefficients
        // std::cout << "Approximation Coefficients (X): ";
        // for (const auto &coeff : pApproxCoeffX) std::cout << coeff.cast<double>() << " ";
        // std::cout << std::endl;

        // Perform reconstruction using upcoef
        py::list pRecX = pywt.attr("upcoef")("a", pApproxCoeffX, abl->waveName, py::arg("level") = J, py::arg("take") = sigLength);
        py::list pRecY = pywt.attr("upcoef")("a", pApproxCoeffY, abl->waveName, py::arg("level") = J, py::arg("take") = sigLength);
        py::list pRecZ = pywt.attr("upcoef")("a", pApproxCoeffZ, abl->waveName, py::arg("level") = J, py::arg("take") = sigLength);

        // Convert reconstructed signals back to C++ array
        for (PetscInt i = 0; i < sigLength; i++) 
        {
            srcPAOut[i].x = pRecX[i].cast<double>();
            srcPAOut[i].y = pRecY[i].cast<double>();
            srcPAOut[i].z = pRecZ[i].cast<double>();
        }
    }
    catch (const py::error_already_set &e) 
    {
        std::cerr << "Python error: " << e.what() << std::endl;
        PetscFunctionReturn(PETSC_ERR_LIB);
    }

    return(0);
}

//*****************************************************************/
PetscErrorCode pywavedecScalar(abl_ *abl, PetscReal *srcPAIn, PetscReal *srcPAOut, PetscInt sigLength)
{

    try 
    {
        // Import PyWavelets
        py::module pywt = py::module::import("pywt");

        // Convert C++ arrays to Python lists
        std::vector<double> signal(sigLength);

        for (PetscInt i = 0; i < sigLength; i++) 
        {
            signal[i] = srcPAIn[i];
        }
        PetscPrintf(PETSC_COMM_WORLD,"scalar %s %ld\n", abl->waveName.c_str(), abl->waveLevel);

        PetscInt J = abl->waveLevel; // Decomposition level

        // Perform wavelet decomposition (wavedec)
        py::list coeffs = pywt.attr("wavedec")(signal, abl->waveName, py::arg("level") = J, py::arg("mode") = abl->waveExtn);

        py::list pApproxCoeff = coeffs[0].cast<py::list>();

        // Perform reconstruction using upcoef
        py::list pRec = pywt.attr("upcoef")("a", pApproxCoeff, abl->waveName, py::arg("level") = J, py::arg("take") = sigLength);

        // Convert reconstructed signals back to C++ array
        for (PetscInt i = 0; i < sigLength; i++) 
        {
            srcPAOut[i] = pRec[i].cast<double>();
        }
    }
    catch (const py::error_already_set &e) 
    {
        std::cerr << "Python error: " << e.what() << std::endl;
        PetscFunctionReturn(PETSC_ERR_LIB);
    }

    return(0);
}
#endif
