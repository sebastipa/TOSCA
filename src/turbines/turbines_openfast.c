#include "turbines_openfast.h"

//***************************************************************************************************************//

#if USE_OPENFAST

PetscErrorCode initOpenFAST(farm_ *farm)
{
    PetscInt nT = farm->size, nOpenFAST   = 0;
    
    // count how many turbines are coupled to OpenFAST
    for(PetscInt t=0; t<nT; t++)
    {
        windTurbine *wt = farm->wt[t];

        // check if this turbine uses OpenFAST 
        if(wt->useOpenFAST)
        {
            nOpenFAST++;    
        }
    }

    // set number of turbines that are coupled to OpenFAST
    farm->nOpenFAST = nOpenFAST;

    if(farm->nOpenFAST)
    {
        farm->FAST                = new fast::OpenFAST;

        farm->fi.comm             = farm->access->mesh->MESH_COMM;
        farm->fi.nTurbinesGlob    = (int)farm->nOpenFAST;
        farm->fi.dryRun           = (bool)0;
        farm->fi.debug            = (bool)farm->dbg;
        farm->fi.tStart           = (double)farm->access->clock->startTime;
        // for now restart is not available (need to check if the file is present)
        farm->fi.simStart         = fast::init; 
        farm->fi.restartFreq      = (int)100; // need to adjust
        farm->fi.outputFreq       = (int)farm->timeInterval; // if this is is seconds need to adjust
        farm->fi.tMax             = (double)farm->access->clock->endTime;
        farm->fi.dtFAST           = (double)(farm->access->clock->dt)/((PetscReal)farm->nFastSubSteps);
        farm->fi.dtDriver         = (double)farm->access->clock->dt;

        farm->fi.globTurbineData.resize(farm->fi.nTurbinesGlob);

        PetscMalloc(sizeof(PetscInt)*nT, &(farm->openfastIds));
        PetscInt count = 0;

        for(PetscInt t=0; t<nT; t++)
        {
            windTurbine *wt = farm->wt[t];

            // check if this turbine uses OpenFAST 
            if(wt->useOpenFAST)
            {
                farm->openfastIds[t] = count;

                // set this index here (dont'use wt because it's a local ptr)
                farm->wt[t]->openfastIndex = count;

                std::vector<float> baseLoc(3);
                baseLoc[0] = farm->base[t].x;
                baseLoc[1] = farm->base[t].y;
                baseLoc[2] = farm->base[t].z;

                std::vector<double> hubLoc(3);
                hubLoc[0] = wt->rotCenter.x;
                hubLoc[1] = wt->rotCenter.y;
                hubLoc[2] = wt->rotCenter.z;

                farm->fi.globTurbineData[count].TurbID              = (int)t;
                farm->fi.globTurbineData[count].FASTInputFileName   = (std::string)(wt->type + "." + std::to_string(t) + ".fst");
                farm->fi.globTurbineData[count].FASTRestartFileName = (std::string)(wt->type + "." + std::to_string(t) + ".chkpt");
                farm->fi.globTurbineData[count].TurbineBasePos      = baseLoc;
                farm->fi.globTurbineData[count].TurbineHubPos       = hubLoc;
                farm->fi.globTurbineData[count].numForcePtsBlade    = (int)wt->nBladeForcePtsOF;
                farm->fi.globTurbineData[count].numForcePtsTwr      = (int)wt->nTwrForcePtsOF;

                count++;
            }
            else 
            {
                farm->openfastIds[t] = -1;
            }
        }

        // send inputs to OpenFAST
        farm->FAST->setInputs(farm->fi);

        // assign OpenFAST instances to the writer rank of each turbine
        for(PetscInt t=0; t<nT; t++)
        {
            if(farm->wt[t]->useOpenFAST)
            {
                farm->FAST->setTurbineProcNo(farm->openfastIds[t], farm->wt[t]->writerRank);
            }
        }

        // initialize OpenFAST (reads OpenFAST input files and creates instance of OpenFAST)
        farm->FAST->init();

        PetscMPIInt   rank, thisTurbRank, nProcs;
        MPI_Comm_rank(farm->access->mesh->MESH_COMM, &rank);
        MPI_Comm_size(farm->access->mesh->MESH_COMM, &nProcs);

        for(PetscInt t=0; t<nT; t++)
        {
            if(farm->wt[t]->useOpenFAST)
            {
                windTurbine *wt = farm->wt[t];

                // set access pointers to the OpenFAST and fastInputs objects (we need them at the turbine level for parallel communication)
                wt->FAST = farm->FAST;
                wt->fi   = &(farm->fi);

                thisTurbRank = farm->FAST->get_procNo(wt->openfastIndex);
                
                // compute number of velocity points 
                PetscInt nVelPtsBlade = 0;
                if(rank == thisTurbRank)
                {
                    // when using OpenFAST, these are set later from the OpenFAST data
                    PetscInt nVelPtsAll   = wt->FAST->get_numVelPts(wt->openfastIndex);
                    PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node
                
                    // get the number of velocity points from OpenFAST 
                    for(PetscInt pi=0; pi<nVelPtsAll; pi++)
                    {
                        nodeType = wt->FAST->getVelNodeType(wt->openfastIndex, pi);
                        if(nodeType == 1) // blade node
                        {
                            nVelPtsBlade++;
                        }
                    }
                }

                // broadcast the number of blade velocity points to all processors
                MPI_Bcast(&nVelPtsBlade, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                if(nVelPtsBlade % wt->nBlades != 0)
                {
                    char error[512];
                    sprintf(error, "number of blade velocity points not multiple of number of blades\n");
                    fatalErrorInFunction("initOpenFAST",  error);
                }
                else 
                {
                    wt->nBladeVelPtsOF = nVelPtsBlade / wt->nBlades;
                }

                // allocate memory for velocity sample points 
                PetscMalloc(nVelPtsBlade*sizeof(Cmpnts), &(wt->velPts));
                PetscMalloc(nVelPtsBlade*sizeof(Cmpnts), &(wt->velVals));

                // compute number of force points 
                PetscInt nForcePtsBlade = 0;
                if(rank == thisTurbRank)
                {
                    // when using OpenFAST, these are set later from the OpenFAST data
                    PetscInt nForcePtsAll   = wt->FAST->get_numForcePts(wt->openfastIndex);
                    PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node
                
                    // get the number of velocity points from OpenFAST 
                    for(PetscInt pi=0; pi<nForcePtsAll; pi++)
                    {
                        nodeType = wt->FAST->getForceNodeType(wt->openfastIndex, pi);
                        if(nodeType == 1) // blade node
                        {
                            nForcePtsBlade++;
                        }
                    }
                }

                // broadcast the number of blade velocity points to all processors
                MPI_Bcast(&nForcePtsBlade, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                if(nForcePtsBlade % wt->nBlades != 0)
                {
                    char error[512];
                    sprintf(error, "number of blade force points not multiple of number of blades\n");
                    fatalErrorInFunction("getForcePtsOpenFAST",  error);
                }
                else 
                {
                    wt->nBladeForcePtsOF = nForcePtsBlade / wt->nBlades;
                }

                // allocate memory for velocity sample points 
                PetscMalloc(nForcePtsBlade*sizeof(Cmpnts), &(wt->forcePts));
                PetscMalloc(nForcePtsBlade*sizeof(Cmpnts), &(wt->forceVals));
            }
        }
    }
    
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode stepZeroOpenFAST(farm_ *farm)
{
    PetscInt nT = farm->size;

    PetscMPIInt   rank;
    MPI_Comm_rank(farm->fi.comm, &rank);

    
    // find which OpenFAST point this processor controls
    findControlledVelPointsOpenFAST(farm);

    // compute and send velocities to openfast 
    computeWindVectorsRotorOpenFAST(farm);

    if (farm->FAST->isTimeZero())
    {
        farm->FAST->solution0();
    }

    for(PetscInt t=0; t<nT; t++)
    {
        if(farm->wt[t]->useOpenFAST)
        {
            // update velocity points from OpenFAST
            getVelPtsBladeOpenFAST(farm->wt[t]);

            // update Force points from OpenFAST
            getForcePtsBladeOpenFAST(farm->wt[t]);

            // update global WT parameters from OpenFAST
            getGlobParamsOpenFAST(farm->wt[t]);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode stepOpenFAST(farm_ *farm)
{
    PetscInt nT = farm->size;

    PetscMPIInt   rank;
    MPI_Comm_rank(farm->fi.comm, &rank);

    for (PetscInt n = 0; n < farm->nFastSubSteps; n++)
    {
        farm->FAST->step((bool)0);
    }

    for(PetscInt t=0; t<nT; t++)
    {
        if(farm->wt[t]->useOpenFAST)
        {
            // update velocity points from OpenFAST
            getVelPtsBladeOpenFAST(farm->wt[t]);

            // update Force points from OpenFAST
            getForcePtsBladeOpenFAST(farm->wt[t]);

            // update global WT parameters from OpenFAST
            getGlobParamsOpenFAST(farm->wt[t]);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initControlledPointsOpenFAST(windTurbine *wt)
{
    // allocate arrays for vel points
    PetscMalloc(wt->nBlades*wt->nBladeVelPtsOF*sizeof(PetscInt), &(wt->thisVelPtControlled));
    PetscMalloc(wt->nBlades*wt->nBladeVelPtsOF*sizeof(cellIds), &(wt->closestVelCells));

    // allocate arrays for force points
    PetscMalloc(wt->nBlades*wt->nBladeForcePtsOF*sizeof(PetscInt), &(wt->thisForcePtControlled));
    PetscMalloc(wt->nBlades*wt->nBladeForcePtsOF*sizeof(cellIds), &(wt->closestForceCells));
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getGlobParamsOpenFAST(windTurbine *wt)
{
    PetscMPIInt   rank, thisTurbRank, nProcs;
    MPI_Comm_rank(wt->fi->comm, &rank);
    MPI_Comm_size(wt->fi->comm, &nProcs);
    thisTurbRank = wt->FAST->get_procNo(wt->openfastIndex);

    // update rotCenter
    std::vector<PetscReal> hubPosition(3,0.0);
    if(rank == thisTurbRank)
    {
        wt->FAST->getHubPos(hubPosition, wt->openfastIndex);
        wt->rotCenter.x = hubPosition[0];
        wt->rotCenter.y = hubPosition[1];
        wt->rotCenter.z = hubPosition[2];
    }
    MPI_Bcast(&wt->rotCenter, 3, MPIU_REAL, thisTurbRank, wt->fi->comm);

    // update rtrAxis
    std::vector<PetscReal> shaftOrientation(3,0.0);
    if(rank == thisTurbRank)
    {
        wt->FAST->getHubShftDir(shaftOrientation, wt->openfastIndex);
        wt->rtrAxis.x = shaftOrientation[0];
        wt->rtrAxis.y = shaftOrientation[1];
        wt->rtrAxis.z = shaftOrientation[2];
        mScale(-1.0, wt->rtrAxis); 
    }
    MPI_Bcast(&wt->rtrAxis, 3, MPIU_REAL, thisTurbRank, wt->fi->comm);
      
    // update rtrDir
    wt->rtrDir = nUnit(nSetFromComponents(wt->rtrAxis.x, wt->rtrAxis.y, 0.0));

    // update omega_hat (this is not really needed for OpenFAST, will compute for output)

    // update rtrOmega (this is not really needed for OpenFAST, will compute for output)

    // update twrDir (this is not really needed for OpenFAST, will compute for output)

    if(wt->dbg && rank == thisTurbRank) 
    {
        printf("\nTurbine %s, OpenFAST rank %d out of %d\n", wt->id.c_str(), rank, nProcs);
        printf(" - rotCenter: %f %f %f\n", wt->rotCenter.x, wt->rotCenter.y, wt->rotCenter.z);
        printf(" - rtrAxis:   %f %f %f\n", wt->rtrAxis.x,   wt->rtrAxis.y,   wt->rtrAxis.z);
        printf(" - rtrDir:    %f %f %f\n", wt->rtrDir.x,    wt->rtrDir.y,    wt->rtrDir.z);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getVelPtsBladeOpenFAST(windTurbine *wt)
{
    PetscMPIInt   rank, thisTurbRank, nProcs;
    MPI_Comm_rank(wt->fi->comm, &rank);
    MPI_Comm_size(wt->fi->comm, &nProcs);
    thisTurbRank = wt->writerRank;

    // temp storage for point location
    std::vector<double> pointLocation(3,0.0);
    PetscInt nVelPtsBlade = wt->nBlades * wt->nBladeVelPtsOF;

    // create local and global vectors for the blade sample points 
    std::vector<Cmpnts>  lVelPoints;
    std::vector<Cmpnts>  gVelPoints;
    lVelPoints.resize(nVelPtsBlade);
    gVelPoints.resize(nVelPtsBlade);

    if(rank == thisTurbRank)
    {
        PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node

        for(PetscInt i=0; i<wt->nBladeVelPtsOF; i++) 
        {
            for(PetscInt j=0; j<wt->nBlades; j++) 
            {
                // OpenFAST index: 1, rad_1_bld_1, rad_N_bld_1, rad_1_bld_2, ..., rad_N_bld_N
                PetscInt pi_of = 1 + i + j*wt->nBladeVelPtsOF; 
                // TOSCA index: rad1_bld1, ..., rad1_bldN, rad2_bld1, ..., radN_bldN
                PetscInt pi_ts = i*wt->nBlades + j;
                
                wt->FAST->getVelNodeCoordinates(pointLocation, pi_of, wt->openfastIndex);
                lVelPoints[pi_ts] = nSetFromComponents(pointLocation[0], pointLocation[1], pointLocation[2]);

                nodeType = wt->FAST->getVelNodeType(wt->openfastIndex, pi_of);
                if(nodeType != 1) // blade node
                {
                    char error[512];
                    sprintf(error, "wrong node type found when getting blade velocity points from OpenFAST\n");
                    fatalErrorInFunction("getVelPtsOpenFAST",  error);
                }
            }
        }
    }

    // scatter the local blade points to all processors
    MPI_Allreduce(&(lVelPoints[0]), &(gVelPoints[0]), nVelPtsBlade*3, MPIU_REAL, MPI_SUM, wt->fi->comm);

    if(rank == thisTurbRank && wt->dbg) 
    {
        printf("\nTurbine %s, OpenFAST rank %d out of %d\n", wt->id.c_str(), rank, nProcs);
        printf(" - total bld vel points: %ld\n", nVelPtsBlade);
        printf(" - vel points per blade: %ld\n", wt->nBladeVelPtsOF);
    }

    // assign them to the turbine structure
    for(PetscInt pi=0; pi<nVelPtsBlade; pi++)
    {
        mSet(wt->velPts[pi], gVelPoints[pi]);
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - bld vel point %ld: %f %f %f\n", pi, wt->velPts[pi].x, wt->velPts[pi].y, wt->velPts[pi].z);
    }

    // clean memory
    std::vector<Cmpnts> ().swap(lVelPoints);
    std::vector<Cmpnts> ().swap(gVelPoints);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getForcePtsBladeOpenFAST(windTurbine *wt)
{
    PetscMPIInt   rank, thisTurbRank, nProcs;
    MPI_Comm_rank(wt->fi->comm, &rank);
    MPI_Comm_size(wt->fi->comm, &nProcs);
    thisTurbRank = wt->writerRank;

    // temp storage for point location
    std::vector<double> pointLocation(3,0.0);
    PetscInt nForcePtsBlade = wt->nBlades * wt->nBladeForcePtsOF;

    // create local and global vectors for the blade sample points 
    std::vector<Cmpnts>  lForcePoints;
    std::vector<Cmpnts>  gForcePoints;
    lForcePoints.resize(nForcePtsBlade);
    gForcePoints.resize(nForcePtsBlade);

    if(rank == thisTurbRank)
    {
        PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node

        for(PetscInt i=0; i<wt->nBladeForcePtsOF; i++) 
        {
            for(PetscInt j=0; j<wt->nBlades; j++) 
            {
                // OpenFAST index: 1, rad_1_bld_1, rad_N_bld_1, rad_1_bld_2, ..., rad_N_bld_N
                PetscInt pi_of = 1 + i + j*wt->nBladeForcePtsOF; 
                // TOSCA index: rad1_bld1, ..., rad1_bldN, rad2_bld1, ..., radN_bldN
                PetscInt pi_ts = i*wt->nBlades + j;
                
                wt->FAST->getForceNodeCoordinates(pointLocation, pi_of, wt->openfastIndex);
                lForcePoints[pi_ts] = nSetFromComponents(pointLocation[0], pointLocation[1], pointLocation[2]);

                nodeType = wt->FAST->getForceNodeType(wt->openfastIndex, pi_of);
                if(nodeType != 1) // blade node
                {
                    char error[512];
                    sprintf(error, "wrong node type found when getting blade force points from OpenFAST\n");
                    fatalErrorInFunction("getForcePtsOpenFAST",  error);
                }
            }
        }
    }

    // scatter the local blade points to all processors
    MPI_Allreduce(&(lForcePoints[0]), &(gForcePoints[0]), nForcePtsBlade*3, MPIU_REAL, MPI_SUM, wt->fi->comm);

    if(rank == thisTurbRank && wt->dbg) 
    {
        printf("\nTurbine %s, OpenFAST rank %d out of %d\n", wt->id.c_str(), rank, nProcs);
        printf(" - total bld force points: %ld\n", nForcePtsBlade);
        printf(" - force points per blade: %ld\n", wt->nBladeForcePtsOF);
    }

    // assign them to the turbine structure
    for(PetscInt pi=0; pi<nForcePtsBlade; pi++)
    {
        mSet(wt->forcePts[pi], gForcePoints[pi]);
        
        mSet(wt->alm.points[pi], gForcePoints[pi]); // also set the ALM points
        
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - bld force point %ld: %f %f %f\n", pi, wt->forcePts[pi].x, wt->forcePts[pi].y, wt->forcePts[pi].z);
    }

    // clean memory
    std::vector<Cmpnts> ().swap(lForcePoints);
    std::vector<Cmpnts> ().swap(gForcePoints);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledVelPointsOpenFAST(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // create point perturbation range for breaking ties in determining the
    // owner processor. Note: the actual position of the points is not changed,
    // this is only to make sure that each velocity point is only controlled by a
    // single processor.

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // DO VELOCITY POINTS
        if(wt->turbineControlled)
        {
            // number of velocity points
            PetscInt npts_t = wt->nBlades*wt->nBladeVelPtsOF;

            // create temporary vectors
            std::vector<PetscReal> lminDist(npts_t);
            std::vector<PetscReal> gminDist(npts_t);
            std::vector<Cmpnts> perturb(npts_t);

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                // initialize min dists to a big value
                lminDist[p] = 1e20;
                gminDist[p] = 1e20;

                // set point perturbation
                perturb[p].x =  procContrib;
                perturb[p].y =  procContrib;
                perturb[p].z =  procContrib;

                // save this point locally for speed
                Cmpnts point_p = wt->velPts[p];

                // perturb the point position
                mSum(point_p, perturb[p]);

                // find the closest cell center
                PetscReal  r_c_minMag = 1e20;
                cellIds closestCell;

                // loop over the sphere cells
                for(c=0; c<wt->nControlled; c++)
                {
                    // cell indices
                    PetscInt i = wt->controlledCells[c].i,
                             j = wt->controlledCells[c].j,
                             k = wt->controlledCells[c].k;

                    // compute distance from mesh cell to AD point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    // compute magnitude
                    PetscReal r_c_mag = nMag(r_c);

                    if(r_c_mag < r_c_minMag)
                    {
                        r_c_minMag = r_c_mag;
                        closestCell.i = i;
                        closestCell.j = j;
                        closestCell.k = k;
                    }
                }

                // save closest cell indices
                wt->closestVelCells[p].i = closestCell.i;
                wt->closestVelCells[p].j = closestCell.j;
                wt->closestVelCells[p].k = closestCell.k;

                // save min dist
                lminDist[p] = r_c_minMag;
            }

            // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
            MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->nBlades*wt->nBladeVelPtsOF, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

            // now compare the lists, if they have the same min distance then the point
            // is controlled by this processor, otherwise not. Ties are break by making
            // the perturbation different among the processors:
            // distance was not equal: it can't become equal
            // distance was equal    : it is made different
            for(p=0; p<npts_t; p++)
            {
                // point is controlled
                if(lminDist[p] == gminDist[p])
                {
                    wt->thisVelPtControlled[p] = 1;
                }
                // point is not controlled
                else
                {
                    wt->thisVelPtControlled[p] = 0;
                }
            }

            // clean memory
            std::vector<PetscReal> ().swap(lminDist);
            std::vector<PetscReal> ().swap(gminDist);
            std::vector<Cmpnts> ().swap(perturb);
        }

        // DO FORCE POINTS
        if(wt->turbineControlled)
        {
            // number of velocity points
            PetscInt npts_t = wt->nBlades*wt->nBladeForcePtsOF;

            // create temporary vectors
            std::vector<PetscReal> lminDist(npts_t);
            std::vector<PetscReal> gminDist(npts_t);
            std::vector<Cmpnts> perturb(npts_t);

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                // initialize min dists to a big value
                lminDist[p] = 1e20;
                gminDist[p] = 1e20;

                // set point perturbation
                perturb[p].x =  procContrib;
                perturb[p].y =  procContrib;
                perturb[p].z =  procContrib;

                // save this point locally for speed
                Cmpnts point_p = wt->forcePts[p];

                // perturb the point position
                mSum(point_p, perturb[p]);

                // find the closest cell center
                PetscReal  r_c_minMag = 1e20;
                cellIds closestCell;

                // loop over the sphere cells
                for(c=0; c<wt->nControlled; c++)
                {
                    // cell indices
                    PetscInt i = wt->controlledCells[c].i,
                             j = wt->controlledCells[c].j,
                             k = wt->controlledCells[c].k;

                    // compute distance from mesh cell to AD point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    // compute magnitude
                    PetscReal r_c_mag = nMag(r_c);

                    if(r_c_mag < r_c_minMag)
                    {
                        r_c_minMag = r_c_mag;
                        closestCell.i = i;
                        closestCell.j = j;
                        closestCell.k = k;
                    }
                }

                // save closest cell indices
                wt->closestForceCells[p].i = closestCell.i;
                wt->closestForceCells[p].j = closestCell.j;
                wt->closestForceCells[p].k = closestCell.k;

                // save min dist
                lminDist[p] = r_c_minMag;
            }

            // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
            MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->nBlades*wt->nBladeForcePtsOF, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

            // now compare the lists, if they have the same min distance then the point
            // is controlled by this processor, otherwise not. Ties are break by making
            // the perturbation different among the processors:
            // distance was not equal: it can't become equal
            // distance was equal    : it is made different
            for(p=0; p<npts_t; p++)
            {
                // point is controlled
                if(lminDist[p] == gminDist[p])
                {
                    wt->thisForcePtControlled[p] = 1;
                }
                // point is not controlled
                else
                {
                    wt->thisForcePtControlled[p] = 0;
                }
            }

            // clean memory
            std::vector<PetscReal> ().swap(lminDist);
            std::vector<PetscReal> ().swap(gminDist);
            std::vector<Cmpnts> ().swap(perturb);
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeWindVectorsRotorOpenFAST(farm_ *farm)
{
    clock_           *clock = farm->access->clock;
    mesh_            *mesh  = farm->access->mesh;
    ueqn_            *ueqn  = farm->access->ueqn;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)
    Cmpnts           ***ucat;   // local vector (no ambiguity in this context)

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            // number of points in the sample mesh
            PetscInt npts_t = wt->nBlades*wt->nBladeVelPtsOF;

            // global inflow velocity for this processor
            std::vector<Cmpnts> lWind(npts_t);

            // global inflow velocity for this processor
            std::vector<Cmpnts> gWind(npts_t);

            // loop over the sample mesh points
            for(p=0; p<npts_t; p++)
            {
                // initialize to zero
                lWind[p].x = 0.0;
                lWind[p].y = 0.0;
                lWind[p].z = 0.0;

                gWind[p].x = 0.0;
                gWind[p].y = 0.0;
                gWind[p].z = 0.0;

                // save this point locally for speed
                Cmpnts point_p = wt->velPts[p];

                // get the closest cell center
                PetscInt i = wt->closestVelCells[p].i,
                         j = wt->closestVelCells[p].j,
                         k = wt->closestVelCells[p].k;

                if(wt->thisVelPtControlled[p])
                {
                    // find the point at which velocity must be sampled
                    Cmpnts sample = nSet(point_p);

                    // trilinear interpolate (do not go upstream for the sample)
                    vectorPointLocalVolumeInterpolation
                    (
                        mesh,
                        sample.x, sample.y, sample.z,
                        i, j, k,
                        cent, ucat, lWind[p]
                    );
                }
            }

            MPI_Allreduce(&(lWind[0]), &(gWind[0]), npts_t*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

            // send data to OpenFAST 
            if(rank == wt->writerRank)
            {
                // set the velocity at each point in OpenFAST
                std::vector<double> pointVelocity(3);
                for(PetscInt i=0; i<wt->nBladeVelPtsOF; i++)
                for(PetscInt j=0; j<wt->nBlades; j++)
                {
                    PetscInt pi_of = 1 + i + j * wt->nBladeVelPtsOF; // OpenFAST index
                    PetscInt pi_ts = i * wt->nBlades + j;            // TOSCA index
                    
                    // set velocity into TOSCA
                    wt->velVals[pi_ts] = nSet(gWind[pi_ts]);
                    
                    // set velocity into OpenFAST
                    pointVelocity[0] = wt->velVals[pi_ts].x;
                    pointVelocity[1] = wt->velVals[pi_ts].y;
                    pointVelocity[2] = wt->velVals[pi_ts].z;
                    wt->FAST->setVelocity(pointVelocity, pi_of, wt->openfastIndex);
                }
                std::vector<double> ().swap(pointVelocity);
            }

            PetscReal rtrAvgMagU = 0.0;
            PetscReal areaSum    = 0.0;

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                // save this point locally for speed
                Cmpnts point_p = wt->velVals[p];

                // this point position from COR
                Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                // build the blade reference frame based on rotation type:
                // 1. if counter-clockwise: z from tip to root,
                //    x as rotor axis from nacelle cone to back,
                //    y blade tangent, directed as wind due to rotation that blade sees.
                // 2. if clockwise: z from root to tip,
                //    x as rotor axis from nacelle cone to back,
                //    y blade tangent, directed as wind due to rotation that blade sees.
                // This is done in order to have the wind vector lying in the positive
                // quadrant in each case (it is a standard practice).

                Cmpnts xb_hat, yb_hat, zb_hat;

                if(wt->rotDir == "cw")
                {
                    zb_hat = nUnit(r_p);
                    xb_hat = nScale(-1.0, wt->rtrAxis);
                    yb_hat = nCross(zb_hat, xb_hat);
                }
                else if(wt->rotDir == "ccw")
                {
                    zb_hat = nUnit(r_p);
                                mScale(-1.0, zb_hat);
                    xb_hat = nScale(-1.0, wt->rtrAxis);
                    yb_hat = nCross(zb_hat, xb_hat);
                }

                // transform the velocity in the bladed reference frame and
                // remove radial (z) component
                PetscReal ub_x = nDot(gWind[p], xb_hat);
                PetscReal ub_y = nDot(gWind[p], yb_hat);

                // uniform weighting
                PetscReal dA   = 1.0;

                rtrAvgMagU += sqrt(ub_x*ub_x + ub_y*ub_y) * dA;

                areaSum += dA;
            }

            if((*farm->turbineModels[t]) == "ADM")
            {
                wt->adm.rtrAvgMagU = rtrAvgMagU / areaSum;
            }
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                wt->uadm.rtrAvgMagU = rtrAvgMagU / areaSum;
            }
            else if((*farm->turbineModels[t]) == "ALM")
            {
                wt->alm.rtrAvgMagU = rtrAvgMagU / areaSum;
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                wt->afm.U  = nScale(rtrAvgMagU / areaSum, wt->rtrDir);
            }
            
            // clean memory
            std::vector<Cmpnts> ().swap(lWind);
            std::vector<Cmpnts> ().swap(gWind);
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeBladeForceOpenFAST(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;

    // turbine and AD mesh point indices
    PetscInt t, p;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        Cmpnts rtrAxis  = farm->wt[t]->rtrAxis;

        // local vectors 
        std::vector<Cmpnts> lForceVals(wt->nBlades*wt->nBladeForcePtsOF);

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            // get data from OpenFAST 
            if(rank == wt->writerRank)
            {
                // set the velocity at each point in OpenFAST
                std::vector<double> pointForce(3);
                for(PetscInt i=0; i<wt->nBladeForcePtsOF; i++)
                for(PetscInt j=0; j<wt->nBlades; j++)
                {
                    PetscInt pi_of = 1 + i + j * wt->nBladeForcePtsOF; // OpenFAST index
                    PetscInt pi_ts = i * wt->nBlades + j;              // TOSCA index
                    wt->FAST->getForce(pointForce, pi_of, wt->openfastIndex);
                    lForceVals[pi_ts].x = pointForce[0];
                    lForceVals[pi_ts].y = pointForce[1];
                    lForceVals[pi_ts].z = pointForce[2];
                }
                std::vector<double> ().swap(pointForce);
            }
            // scatter to all ransk of this communicator 
            MPI_Allreduce(&(lForceVals[0]), &(wt->forceVals[0]), wt->nBlades*wt->nBladeForcePtsOF*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

            // number of points in the sample mesh
            PetscInt npts_t = wt->nBlades*wt->nBladeForcePtsOF;
            
            for(p=0; p<npts_t; p++)
            {
                // compute body force (flow on blade)
                wt->alm.B[p]  = nSet(wt->forceVals[p]);
            }
        }
    }

    return(0);  
}

#endif