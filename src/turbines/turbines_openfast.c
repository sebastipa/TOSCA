#include "turbines_openfast.h"

//***************************************************************************************************************//

#if USE_OPENFAST

PetscErrorCode initOpenFAST(farm_ *farm)
{
    PetscInt nT = farm->size, nOpenFAST = 0;
    
    // allocate memory for OpenFAST (different from TOSCA when not all turbines use OpenFAST)
    PetscMalloc(sizeof(PetscInt)*nT, &(farm->openfastIds));
    
    // count how many turbines are coupled to OpenFAST
    for(PetscInt t=0; t<nT; t++)
    {
        windTurbine *wt = farm->wt[t];

        // initialize global openfast index list to -1
        farm->openfastIds[t] = -1;

        // initialize openfast index to -1
        wt->openfastIndex    = -1;

        // check if this turbine uses OpenFAST 
        if(wt->useOpenFAST)
        {
            nOpenFAST++;    
        }
    }

    // set number of turbines that are coupled to OpenFAST
    farm->nOpenFAST = nOpenFAST;

    // print info 
    PetscPrintf(farm->access->mesh->MESH_COMM, "    > turbines coupled with OpenFAST: %ld/%ld\n", farm->nOpenFAST, farm->size);

    // there are coupled turbines 
    if(farm->nOpenFAST)
    {
        // initialize OpenFAST interface
        farm->FAST                = new fast::OpenFAST;

        // set OpenFAST inputs
        farm->fi.comm             = farm->access->mesh->MESH_COMM;
        farm->fi.nTurbinesGlob    = (int)farm->nOpenFAST;
        farm->fi.dryRun           = (bool)0;
        farm->fi.debug            = (bool)farm->dbg;
        farm->fi.tStart           = (double)farm->access->clock->startTime;
        farm->fi.simStart         = fast::init;                           // restart to be implemented
        farm->fi.restartFreq      = (int)100;                             // need to dynamically define it
        farm->fi.outputFreq       = (int)farm->timeInterval;              // need to dynamically define it
        farm->fi.tMax             = (double)farm->access->clock->endTime; 
        farm->fi.dtFAST           = (double)(farm->access->clock->dt)/((double)farm->nFastSubSteps);
        farm->fi.dtDriver         = (double)farm->access->clock->dt;

        // allocate initial OpenFAST iteration 
        farm->iterOpenFAST        = 0; // need to dynamically define it

        // allocate memory for turbine data
        farm->fi.globTurbineData.resize(farm->fi.nTurbinesGlob);
        
        PetscInt count = 0;
        for(PetscInt t=0; t<nT; t++)
        {
            windTurbine *wt = farm->wt[t];

            // check if this turbine uses OpenFAST 
            if(wt->useOpenFAST)
            {
                // set OpenFAST id for this turbine in global array (required for all interface functions)
                farm->openfastIds[t] = count;

                // set OpenFAST id for this turbine in turbine structure (required for all interface functions)
                wt->openfastIndex    = count;

                // get turbine base location 
                std::vector<float> baseLoc(3);
                baseLoc[0] = farm->base[t].x;
                baseLoc[1] = farm->base[t].y;
                baseLoc[2] = farm->base[t].z;

                // get turbine hub location 
                std::vector<double> hubLoc(3);
                hubLoc[0] = wt->rotCenter.x;
                hubLoc[1] = wt->rotCenter.y;
                hubLoc[2] = wt->rotCenter.z;

                // set input and checkpoint file names 
                std::string fastInputFile   = (std::string)("turbines/" + wt->type + "." + std::to_string(t) + ".fst");
                std::string fastRestartFile = (std::string)("turbines/" + wt->type + "." + std::to_string(t) + ".chkpt");

                farm->fi.globTurbineData[count].TurbID              = (int)t;
                farm->fi.globTurbineData[count].FASTInputFileName   = fastInputFile;
                farm->fi.globTurbineData[count].FASTRestartFileName = fastRestartFile;
                farm->fi.globTurbineData[count].TurbineBasePos      = baseLoc;
                farm->fi.globTurbineData[count].TurbineHubPos       = hubLoc;
                farm->fi.globTurbineData[count].numForcePtsBlade    = (int)wt->nBladeForcePtsOF;
                farm->fi.globTurbineData[count].numForcePtsTwr      = (int)wt->nTwrForcePtsOF;

                count++;
            }
        }

        // send inputs to OpenFAST (called by all)
        farm->FAST->setInputs(farm->fi);

        // do processor assignment: rank 0 of each TRB_COMM will have openfast control that turbine
        for(PetscInt t=0; t<nT; t++)
        {
            if(farm->wt[t]->useOpenFAST)
            {
                // called by all so that all know who is in control of this turbine
                farm->FAST->setTurbineProcNo(farm->openfastIds[t], farm->wt[t]->writerRank);
            }
        }

        // initialize OpenFAST (called by all)
        farm->FAST->init();

        // now we have to get from OpenFAST velocity and force points for rotor and tower
        PetscMPIInt   rank, thisTurbRank, nProcs;
        MPI_Comm_rank(farm->access->mesh->MESH_COMM, &rank);
        MPI_Comm_size(farm->access->mesh->MESH_COMM, &nProcs);

        // velocity points on rotor and tower 
        for(PetscInt t=0; t<nT; t++)
        {
            if(farm->wt[t]->useOpenFAST)
            {
                windTurbine *wt = farm->wt[t];

                // set access pointers to the OpenFAST at turbine level
                wt->FAST = farm->FAST;
                wt->fi   = &(farm->fi);

                // get rank that solves this turbine (called by all)
                thisTurbRank = farm->FAST->get_procNo(wt->openfastIndex);
                
                // compute number of velocity points 
                PetscInt nVelPtsBlade = 0;
                PetscInt nVelPtsTower = 0;
                if(rank == thisTurbRank)
                {
                    // get total number of velocity points (order is: 1 pt for nacelle, blade1, ... bladeN, tower)
                    PetscInt nVelPtsAll   = wt->FAST->get_numVelPts(wt->openfastIndex);
                    PetscInt nodeType     = 0; // 0: hub node, 1: blade node, 2: tower node
                
                    // discriminate between rotor and tower points 
                    for(PetscInt pi=0; pi<nVelPtsAll; pi++)
                    {
                        nodeType = wt->FAST->getVelNodeType(wt->openfastIndex, pi);
                        if(nodeType == 1)     // blade node
                        {
                            nVelPtsBlade++;
                        }
                        else if(nodeType == 2) // tower node
                        {
                            nVelPtsTower++;
                        }
                    }
                }

                // broadcast the number of blade and tower velocity points to all processors
                MPI_Bcast(&nVelPtsBlade, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                MPI_Bcast(&nVelPtsTower, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                if(nVelPtsBlade % wt->nBlades != 0)
                {
                    char error[512];
                    sprintf(error, "number of blade velocity points not multiple of number of blades\n");
                    fatalErrorInFunction("initOpenFAST",  error);
                }
                else 
                {
                    // set number of velocity points for one blade
                    wt->nBladeVelPtsOF = nVelPtsBlade / wt->nBlades;

                    // set number of tower velocity points
                    if(wt->includeTwr)
                    {
                        wt->nTwrVelPtsOF   = nVelPtsTower;
                    }
                }

                // allocate memory for rotor velocity sample points 
                PetscMalloc(nVelPtsBlade*sizeof(Cmpnts), &(wt->velPtsBlade));
                PetscMalloc(nVelPtsBlade*sizeof(Cmpnts), &(wt->velValsBlade));

                // allocate memory for tower velocity sample points 
                if(wt->includeTwr)
                {
                    PetscMalloc(nVelPtsTower*sizeof(Cmpnts), &(wt->velPtsTwr));
                    PetscMalloc(nVelPtsTower*sizeof(Cmpnts), &(wt->velValsTwr));
                }

                // compute number of force points 
                PetscInt nForcePtsBlade = 0;
                PetscInt nForcePtsTower = 0;
                if(rank == thisTurbRank)
                {
                    // get total number of force points (order is: 1 pt for nacelle, blade1, ... bladeN, tower)
                    PetscInt nForcePtsAll   = wt->FAST->get_numForcePts(wt->openfastIndex);
                    PetscInt nodeType       = 0; // 0: hub node, 1: blade node, 2: tower node
                
                    // discriminate between rotor and tower points
                    for(PetscInt pi=0; pi<nForcePtsAll; pi++)
                    {
                        nodeType = wt->FAST->getForceNodeType(wt->openfastIndex, pi);
                        if(nodeType == 1)      // blade node
                        {
                            nForcePtsBlade++;
                        }
                        else if(nodeType == 2) // tower node
                        {
                            nForcePtsTower++;
                        }
                    }
                }

                // broadcast the number of blade and tower force points to all processors
                MPI_Bcast(&nForcePtsBlade, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                MPI_Bcast(&nForcePtsTower, 1, MPIU_INT, thisTurbRank, wt->fi->comm);
                if(nForcePtsBlade % wt->nBlades != 0)
                {
                    char error[512];
                    sprintf(error, "number of blade force points not multiple of number of blades\n");
                    fatalErrorInFunction("getForcePtsOpenFAST",  error);
                }
                else 
                {
                    // set number of force points for one blade 
                    // here is a bit ugly: we reset it based on OpenFAST, but we initially passed the same value to OpenFAST
                    wt->nBladeForcePtsOF = nForcePtsBlade / wt->nBlades;

                    // set number of tower force points
                    // here is a bit ugly: we reset it based on OpenFAST, but we initially passed the same value to OpenFAST
                    if(wt->includeTwr)
                    {
                        wt->nTwrForcePtsOF   = nForcePtsTower;
                    }
                }

                // allocate memory for force sample points 
                PetscMalloc(nForcePtsBlade*sizeof(Cmpnts), &(wt->forcePtsBlade));
                PetscMalloc(nForcePtsBlade*sizeof(Cmpnts), &(wt->forceValsBlade));

                // allocate memory for tower force sample points 
                if(wt->includeTwr)
                {
                    PetscMalloc(nForcePtsTower*sizeof(Cmpnts), &(wt->forcePtsTwr));
                    PetscMalloc(nForcePtsTower*sizeof(Cmpnts), &(wt->forceValsTwr));
                } 
                
                // allocate arrays for blade vel points
                PetscMalloc(wt->nBlades*wt->nBladeVelPtsOF*sizeof(PetscInt), &(wt->thisBladeVelPtControlled));
                PetscMalloc(wt->nBlades*wt->nBladeVelPtsOF*sizeof(cellIds), &(wt->closestBladeVelCells));

                // allocate arrays for tower vel points
                if (wt->includeTwr)
                {
                    PetscMalloc(wt->nTwrVelPtsOF*sizeof(PetscInt), &(wt->thisTwrVelPtControlled));
                    PetscMalloc(wt->nTwrVelPtsOF*sizeof(cellIds), &(wt->closestTwrVelCells));
                }

                // allocate arrays for blade force points
                PetscMalloc(wt->nBlades*wt->nBladeForcePtsOF*sizeof(PetscInt), &(wt->thisBladeForcePtControlled));
                PetscMalloc(wt->nBlades*wt->nBladeForcePtsOF*sizeof(cellIds), &(wt->closestBladeForceCells));

                // allocate arrays for tower force points
                if (wt->includeTwr)
                {
                    PetscMalloc(wt->nTwrForcePtsOF*sizeof(PetscInt), &(wt->thisTwrForcePtControlled));
                    PetscMalloc(wt->nTwrForcePtsOF*sizeof(cellIds), &(wt->closestTwrForceCells));
                }
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

    // find which rotor point this processor controls 
    findControlledPointsRotorOpenFAST(farm);

    // find which tower point this processor controls 
    findControlledPointsTowerOpenFAST(farm);

    // compute and send velocities to openfast 
    computeWindVectorsRotorOpenFAST(farm);

    // compute and send tower velocities to openfast 
    computeWindVectorsTowerOpenFAST(farm);

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
            if(farm->wt[t]->includeTwr)
            {
                getVelPtsTwrOpenFAST(farm->wt[t]);
            }

            // update Force points from OpenFAST
            getForcePtsBladeOpenFAST(farm->wt[t]);
            if(farm->wt[t]->includeTwr)
            {
                getForcePtsTwrOpenFAST(farm->wt[t]);
            }

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
        bool writeFilesOF = false;
        farm->FAST->step(writeFilesOF);
    }

    for(PetscInt t=0; t<nT; t++)
    {
        if(farm->wt[t]->useOpenFAST)
        {
            // update velocity points from OpenFAST
            getVelPtsBladeOpenFAST(farm->wt[t]);
            if(farm->wt[t]->includeTwr)
            {
                getVelPtsTwrOpenFAST(farm->wt[t]);
            }

            // update Force points from OpenFAST
            getForcePtsBladeOpenFAST(farm->wt[t]);
            if(farm->wt[t]->includeTwr)
            {
                getForcePtsTwrOpenFAST(farm->wt[t]);
            }

            // update global WT parameters from OpenFAST
            getGlobParamsOpenFAST(farm->wt[t]);
        }
    }

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
                    fatalErrorInFunction("getVelPtsBladeOpenFAST",  error);
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
        mSet(wt->velPtsBlade[pi], gVelPoints[pi]);
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - bld vel point %ld: %f %f %f\n", pi, wt->velPtsBlade[pi].x, wt->velPtsBlade[pi].y, wt->velPtsBlade[pi].z);
    }

    // clean memory
    std::vector<Cmpnts> ().swap(lVelPoints);
    std::vector<Cmpnts> ().swap(gVelPoints);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getVelPtsTwrOpenFAST(windTurbine *wt)
{
    PetscMPIInt   rank, thisTurbRank, nProcs;
    MPI_Comm_rank(wt->fi->comm, &rank);
    MPI_Comm_size(wt->fi->comm, &nProcs);
    thisTurbRank = wt->writerRank;

    // temp storage for point location
    std::vector<double> pointLocation(3,0.0);
    PetscInt nVelPtsTower = wt->nTwrVelPtsOF;

    // create local and global vectors for the blade sample points 
    std::vector<Cmpnts>  lVelPoints;
    std::vector<Cmpnts>  gVelPoints;
    lVelPoints.resize(nVelPtsTower);
    gVelPoints.resize(nVelPtsTower);

    if(rank == thisTurbRank)
    {
        PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node

        for(PetscInt i=0; i<wt->nTwrVelPtsOF; i++) 
        {
            // OpenFAST index: 1, rad_1_bld_1, rad_N_bld_1, rad_1_bld_2, ..., rad_N_bld_N, tower1, ..., towerN
            PetscInt pi_of = 1 + i + wt->nBladeVelPtsOF*wt->nBlades; 
            // TOSCA index: tower1, ..., towerN
            PetscInt pi_ts = i;
            
            wt->FAST->getVelNodeCoordinates(pointLocation, pi_of, wt->openfastIndex);
            lVelPoints[pi_ts] = nSetFromComponents(pointLocation[0], pointLocation[1], pointLocation[2]);

            nodeType = wt->FAST->getVelNodeType(wt->openfastIndex, pi_of);
            if(nodeType != 2) // tower node
            {
                char error[512];
                sprintf(error, "wrong node type found when getting tower velocity points from OpenFAST\n");
                fatalErrorInFunction("getVelPtsTwrOpenFAST",  error);
            }
        }
    }

    // scatter the local blade points to all processors
    MPI_Allreduce(&(lVelPoints[0]), &(gVelPoints[0]), nVelPtsTower*3, MPIU_REAL, MPI_SUM, wt->fi->comm);

    // assign them to the turbine structure
    for(PetscInt pi=0; pi<wt->nTwrVelPtsOF; pi++)
    {
        mSet(wt->velPtsTwr[pi], gVelPoints[pi]);
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - twr vel point %ld: %f %f %f\n", pi, wt->velPtsTwr[pi].x, wt->velPtsTwr[pi].y, wt->velPtsTwr[pi].z);
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
                    fatalErrorInFunction("getForcePtsBladeOpenFAST",  error);
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
        mSet(wt->forcePtsBlade[pi], gForcePoints[pi]);
        
        mSet(wt->alm.points[pi], gForcePoints[pi]); // also set the ALM points
        
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - bld force point %ld: %f %f %f\n", pi, wt->forcePtsBlade[pi].x, wt->forcePtsBlade[pi].y, wt->forcePtsBlade[pi].z);
    }

    // clean memory
    std::vector<Cmpnts> ().swap(lForcePoints);
    std::vector<Cmpnts> ().swap(gForcePoints);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode getForcePtsTwrOpenFAST(windTurbine *wt)
{
    PetscMPIInt   rank, thisTurbRank, nProcs;
    MPI_Comm_rank(wt->fi->comm, &rank);
    MPI_Comm_size(wt->fi->comm, &nProcs);
    thisTurbRank = wt->writerRank;

    // temp storage for point location
    std::vector<double> pointLocation(3,0.0);
    PetscInt nForcePtsTower = wt->nTwrForcePtsOF;

    // create local and global vectors for the blade sample points 
    std::vector<Cmpnts>  lForcePoints;
    std::vector<Cmpnts>  gForcePoints;
    lForcePoints.resize(nForcePtsTower);
    gForcePoints.resize(nForcePtsTower);

    if(rank == thisTurbRank)
    {
        PetscInt nodeType = 0; // 0: hub node, 1: blade node, 2: tower node

        for(PetscInt i=0; i<wt->nTwrForcePtsOF; i++) 
        {
            // OpenFAST index: 1, rad_1_bld_1, rad_N_bld_1, rad_1_bld_2, ..., rad_N_bld_N, tower1, ..., towerN
            PetscInt pi_of = 1 + i + wt->nBladeForcePtsOF*wt->nBlades; 
            // TOSCA index: tower1, ..., towerN
            PetscInt pi_ts = i;
            
            wt->FAST->getForceNodeCoordinates(pointLocation, pi_of, wt->openfastIndex);
            lForcePoints[pi_ts] = nSetFromComponents(pointLocation[0], pointLocation[1], pointLocation[2]);

            nodeType = wt->FAST->getForceNodeType(wt->openfastIndex, pi_of);
            if(nodeType != 2) // tower node
            {
                char error[512];
                sprintf(error, "wrong node type found when getting tower force points from OpenFAST\n");
                fatalErrorInFunction("getForcePtsTwrOpenFAST",  error);
            }
        }
    }

    // scatter the local blade points to all processors
    MPI_Allreduce(&(lForcePoints[0]), &(gForcePoints[0]), nForcePtsTower*3, MPIU_REAL, MPI_SUM, wt->fi->comm);

    // assign them to the turbine structure
    for(PetscInt pi=0; pi<wt->nTwrForcePtsOF; pi++)
    {
        mSet(wt->forcePtsTwr[pi], gForcePoints[pi]);
        
        mSet(wt->twr.points[pi], gForcePoints[pi]); // also set the displaced tower points
        
        if(rank == thisTurbRank && wt->dbg) 
        printf(" - twr force point %ld: %f %f %f\n", pi, wt->forcePtsTwr[pi].x, wt->forcePtsTwr[pi].y, wt->forcePtsTwr[pi].z);
    }

    // clean memory
    std::vector<Cmpnts> ().swap(lForcePoints);
    std::vector<Cmpnts> ().swap(gForcePoints);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledPointsRotorOpenFAST(farm_ *farm)
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
                Cmpnts point_p = wt->velPtsBlade[p];

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
                wt->closestBladeVelCells[p].i = closestCell.i;
                wt->closestBladeVelCells[p].j = closestCell.j;
                wt->closestBladeVelCells[p].k = closestCell.k;

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
                    wt->thisBladeVelPtControlled[p] = 1;
                }
                // point is not controlled
                else
                {
                    wt->thisBladeVelPtControlled[p] = 0;
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
                Cmpnts point_p = wt->forcePtsBlade[p];

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
                wt->closestBladeForceCells[p].i = closestCell.i;
                wt->closestBladeForceCells[p].j = closestCell.j;
                wt->closestBladeForceCells[p].k = closestCell.k;

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
                    wt->thisBladeForcePtControlled[p] = 1;
                }
                // point is not controlled
                else
                {
                    wt->thisBladeForcePtControlled[p] = 0;
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

PetscErrorCode findControlledPointsTowerOpenFAST(farm_ *farm)
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

        if(wt->includeTwr)
        {
            // DO VELOCITY POINTS
            if(wt->twr.nControlled)
            {
                // number of velocity points
                PetscInt npts_t = wt->nTwrVelPtsOF;

                // create temporary vectors
                std::vector<PetscReal> lminDist(npts_t);
                std::vector<PetscReal> gminDist(npts_t);
                std::vector<Cmpnts> perturb(npts_t);

                // loop over the tower mesh points
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
                    Cmpnts point_p = wt->velPtsTwr[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->twr.nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->twr.controlledCells[c].i,
                                 j = wt->twr.controlledCells[c].j,
                                 k = wt->twr.controlledCells[c].k;

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
                    wt->closestTwrVelCells[p].i = closestCell.i;
                    wt->closestTwrVelCells[p].j = closestCell.j;
                    wt->closestTwrVelCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                // this call can be in the turbineControlled test as long as the communicator is TWR_COMM (will hang otherwise)
                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->nTwrVelPtsOF, MPIU_REAL, MPIU_MIN, wt->twr.TWR_COMM);

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
                        wt->thisTwrVelPtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->thisTwrVelPtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts> ().swap(perturb);
            }

            // DO FORCE POINTS
            if(wt->twr.nControlled)
            {
                // number of velocity points
                PetscInt npts_t = wt->nTwrForcePtsOF;

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
                    Cmpnts point_p = wt->forcePtsTwr[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->twr.nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->twr.controlledCells[c].i,
                                 j = wt->twr.controlledCells[c].j,
                                 k = wt->twr.controlledCells[c].k;

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
                    wt->closestTwrForceCells[p].i = closestCell.i;
                    wt->closestTwrForceCells[p].j = closestCell.j;
                    wt->closestTwrForceCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->nTwrForcePtsOF, MPIU_REAL, MPIU_MIN, wt->twr.TWR_COMM);

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
                        wt->thisTwrForcePtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->thisTwrForcePtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts> ().swap(perturb);
            }
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
                Cmpnts point_p = wt->velPtsBlade[p];

                // get the closest cell center
                PetscInt i = wt->closestBladeVelCells[p].i,
                         j = wt->closestBladeVelCells[p].j,
                         k = wt->closestBladeVelCells[p].k;

                if(wt->thisBladeVelPtControlled[p])
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
                    wt->velValsBlade[pi_ts] = nSet(gWind[pi_ts]);
                    
                    // set velocity into OpenFAST
                    pointVelocity[0] = wt->velValsBlade[pi_ts].x;
                    pointVelocity[1] = wt->velValsBlade[pi_ts].y;
                    pointVelocity[2] = wt->velValsBlade[pi_ts].z;
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
                Cmpnts point_p = wt->velValsBlade[p];

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

PetscErrorCode computeWindVectorsTowerOpenFAST(farm_ *farm)
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

        // test if tower is modeled in this turbine
        if(wt->includeTwr)
        {
            // test if this processor controls this tower
            if(wt->twr.nControlled)
            {
                // number of points in the sample mesh
                PetscInt npts_t = wt->nTwrVelPtsOF;

                // local velocity for this processor
                std::vector<Cmpnts> lU(npts_t);

                // loop over the sample mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize to zero
                    lU[p].x = 0.0;
                    lU[p].y = 0.0;
                    lU[p].z = 0.0;

                    // save this point locally for speed
                    Cmpnts point_p = wt->velPtsTwr[p];

                    if(wt->thisTwrVelPtControlled[p])
                    {
                        // get the closest cell center
                        PetscInt i = wt->closestTwrVelCells[p].i,
                                 j = wt->closestTwrVelCells[p].j,
                                 k = wt->closestTwrVelCells[p].k;
                        
                        // find the point at which velocity must be sampled
                        Cmpnts sample = nSet(point_p);

                        // trilinear interpolate (do not go upstream for the sample)
                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            sample.x, sample.y, sample.z,
                            i, j, k,
                            cent, ucat, lU[p]
                        );
                    }
                }

                MPI_Allreduce(&(lU[0]), &(wt->twr.U[0]), wt->twr.nPoints*3, MPIU_REAL, MPIU_SUM, wt->twr.TWR_COMM);

                // send data to OpenFAST 
                if(rank == wt->writerRank)
                {
                    // set the velocity at each point in OpenFAST
                    std::vector<double> pointVelocity(3);
                    for(PetscInt i=0; i<wt->nTwrVelPtsOF; i++)
                    {
                        PetscInt pi_of = 1 + i + wt->nBladeVelPtsOF*wt->nBlades; // OpenFAST index
                        PetscInt pi_ts = i;                                      // TOSCA index
                        
                        // set velocity into TOSCA
                        wt->velValsTwr[pi_ts] = nSet(wt->twr.U[pi_ts]);
                        
                        // set velocity into OpenFAST
                        pointVelocity[0] = wt->velValsTwr[pi_ts].x;
                        pointVelocity[1] = wt->velValsTwr[pi_ts].y;
                        pointVelocity[2] = wt->velValsTwr[pi_ts].z;
                        wt->FAST->setVelocity(pointVelocity, pi_of, wt->openfastIndex);
                    }
                    std::vector<double> ().swap(pointVelocity);
                }

                // clean memory
                std::vector<Cmpnts> ().swap(lU);
            }
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
            MPI_Allreduce(&(lForceVals[0]), &(wt->forceValsBlade[0]), wt->nBlades*wt->nBladeForcePtsOF*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

            // number of points in the sample mesh
            PetscInt npts_t = wt->nBlades*wt->nBladeForcePtsOF;
            
            for(p=0; p<npts_t; p++)
            {
                // compute body force (flow on blade)
                wt->alm.B[p]  = nSet(wt->forceValsBlade[p]);
            }
        }
    }

    return(0);  
}

#endif