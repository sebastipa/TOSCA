//! \file  initialization.c
//! \brief Contains simulation initialization function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/boundary.h"
#include "include/inflow.h"
#include "include/abl.h"
#include "include/turbines.h"
#include "include/initialization.h"     

//***************************************************************************************************************//

PetscErrorCode PrintOkWindLogo()
{
    PetscPrintf(PETSC_COMM_WORLD," \n\n");
    PetscPrintf(PETSC_COMM_WORLD," |========================================================|  TOSCA: The Toolbox fOr Stratified Convective Atmospheres\n");
    PetscPrintf(PETSC_COMM_WORLD," |                  \\\\                                    |  \n");
    PetscPrintf(PETSC_COMM_WORLD," |                   \\\\        T                          |  Copyright (c) 2023, The University of British Columbia. All rights reserved. This code is under\n");
    PetscPrintf(PETSC_COMM_WORLD," |                    \\\\    \\\\  O                         |  active development at UBC-Okanagan by Sebastiano Stipa, Arjun Ajay, Joshua Brinkerhoff.\n");
    PetscPrintf(PETSC_COMM_WORLD," |                     \\\\    \\\\  S                        |  \n");
    PetscPrintf(PETSC_COMM_WORLD," |                      \\\\    \\\\  C                       |  TOSCA is free software: you can redistribute it and/or modify it under the terms of the BSD 2-\n");
    PetscPrintf(PETSC_COMM_WORLD," |                       \\\\    \\\\  A                      |  Clause Simplified License.\n");
    PetscPrintf(PETSC_COMM_WORLD," |                        O ====\\\\=========               |  \n");
    PetscPrintf(PETSC_COMM_WORLD," |                       //|     \\\\                       |  Redistribution and use in source and binary forms, with or without modification, are permitted\n");
    PetscPrintf(PETSC_COMM_WORLD," |                      //||      O =============         |  provided that the following conditions are met:\n");
    PetscPrintf(PETSC_COMM_WORLD," |                     // ||     //| v2.0.0               |  1. Redistributions of source code must retain the above copyright notice, this list of conditi\n");
    PetscPrintf(PETSC_COMM_WORLD," |                    //  ||    //||                      |     ons and the following disclaimer.\n");
    PetscPrintf(PETSC_COMM_WORLD," |                   //   ||   // ||                      |  2. Redistributions in binary form must reproduce the above copyright notice, this list of cond\n");
    PetscPrintf(PETSC_COMM_WORLD," |                  //    ||  //  ||                      |     itions and the following disclaimer in the documentation and/or other materials provided wi\n");
    PetscPrintf(PETSC_COMM_WORLD," |                        || //   ||                      |     th the distribution.\n");
    PetscPrintf(PETSC_COMM_WORLD," |                        ||//    ||                      |  \n");
    PetscPrintf(PETSC_COMM_WORLD," |                        ||      ||                      |  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS -AS IS- AND ANY EXPRESS OR\n");
    PetscPrintf(PETSC_COMM_WORLD," |                        ||      ||                      |  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY A\n");
    PetscPrintf(PETSC_COMM_WORLD," |                        ||      ||                      |  ND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR \n");
    PetscPrintf(PETSC_COMM_WORLD," |________________________||______||______________________|  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENT\n");     
    PetscPrintf(PETSC_COMM_WORLD," | Copyright : University of British Columbia (Okanagan)  |  IAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS \n");    
    PetscPrintf(PETSC_COMM_WORLD," | Authors   : Stipa - Ajay - Haji - Brinkerhoff          |  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABIL\n");    
    PetscPrintf(PETSC_COMM_WORLD," |========================================================|  ITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISIN\n");        
    PetscPrintf(PETSC_COMM_WORLD," |     Toolbox for Stratified Convective Atmospheres      |  G IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAG\n");        
    PetscPrintf(PETSC_COMM_WORLD," |________________________________________________________|  E.\n\n\n\n");       

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode PrintNumberOfProcs()
{
    PetscMPIInt nProcs;

    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);
    PetscPrintf(PETSC_COMM_WORLD,"Job running with %ld processors\n", nProcs);
    PetscPrintf(PETSC_COMM_WORLD, "******************************************************************\n\n");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode simulationInitialize(domain_ **domainAddr, clock_ *clock, simInfo_ *info, flags_ *flags)
{
    // print logo
    PrintOkWindLogo();

    // print number of processors being used for this simulation
    PrintNumberOfProcs();

    // reads parent/child tree and allocates os memory
    SetDomainsAndAllocate(domainAddr, flags, info);

    domain_ *domain = *domainAddr;

    // timers 
    PetscReal timeStart, timeEnd;

    // set simulation start time
    SetStartTime(clock, domain,info);

    // domains 
    for(PetscInt d=0; d<info->nDomains; d++)
    {
        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);

        PetscTime(&timeStart);
        
        PetscPrintf(PETSC_COMM_WORLD, "\nInitializing domain %ld\n", d);
        PetscPrintf(PETSC_COMM_WORLD, "******************************************************************\n\n");

        // set pointer to time controls
        domain[d].clock = clock;

        domain[d].domainID = d;

        // read physical constants
        ReadPhysicalConstants(&domain[d]);

        // set access database pointers
        SetAccessPointers(&domain[d]);

        // set boundary conditions
        SetBoundaryConditions(domain[d].mesh);

        // initialize mesh
        InitializeMesh(domain[d].mesh);

        // initialize i/o controls and initialization type
        InitializeIO(domain[d].io);
        
        // set wall models
        SetWallModels(domain[d].ueqn);

        // initialize ibm
        InitializeIBM(domain[d].ibm);
        
        // set inflow functions
        SetInflowFunctions(domain[d].mesh);

        // initialize ABL parameters
        InitializeABL(domain[d].abl);

        // momentum equation initialize
        InitializeUEqn(domain[d].ueqn);

        // poisson equation initialize
        InitializePEqn(domain[d].peqn);

        // temperature equation initialize
        InitializeTEqn(domain[d].teqn);

        // LES model initialize
        InitializeLES(domain[d].les);

        // LES scalar model initialize
        InitializeLESScalar(domain[d].les);

        // initialize wind farm
        InitializeWindFarm(domain[d].farm);

        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);

        PetscTime(&timeEnd);

        PetscPrintf(PETSC_COMM_WORLD, "\nFinished initializing domain %ld: elapsed time = %lf s\n\n", d, timeEnd - timeStart); 
    }

    // acquisition system 
    {
        PetscPrintf(PETSC_COMM_WORLD, "\nInitializing acquisition system\n");
        PetscPrintf(PETSC_COMM_WORLD, "******************************************************************\n\n");


        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);

        PetscTime(&timeStart);

        // initialize acquisitions
        InitializeAcquisition(domain);

        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);

        PetscTime(&timeEnd);

        PetscPrintf(PETSC_COMM_WORLD, "Finished initializing acquisition system: elapsed time = %lf s\n\n", timeEnd - timeStart);
    }

    // overset & initial field
    if(info->nDomains == 1)
    {
        SetInitialField(&domain[0]);
    }
    else
    {
        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);
    
        // get timer 
        PetscTime(&timeStart);
    
        PetscPrintf(PETSC_COMM_WORLD, "\nOverset initialization:\n");
        PetscPrintf(PETSC_COMM_WORLD, "******************************************************************\n\n");

        // initialize overset
        InitializeOverset(domain);

        // sync processors
        MPI_Barrier(PETSC_COMM_WORLD);

        PetscTime(&timeEnd);

        PetscPrintf(PETSC_COMM_WORLD, "\nFinished initializing overset: elapsed time = %lf s\n", timeEnd - timeStart);
    }

    return(0);

}

//***************************************************************************************************************//

PetscErrorCode SetSimulationFlags(flags_ *flags)
{
    // read from control file
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    flags->isOversetActive               = 0;
    flags->isLesActive                   = 1;
    flags->isTeqnActive                  = 0;
    flags->isWindFarmActive              = 0;
    flags->isAquisitionActive            = 0;
    flags->isAblActive                   = 0;
    flags->isIBMActive                   = 0;
    flags->isZDampingActive              = 0;
    flags->isXDampingActive              = 0;
    flags->isYDampingActive              = 0;
    flags->isKLeftRayleighDampingActive  = 0;
    flags->isKRightRayleighDampingActive = 0;
    flags->isAdvectionDampingXActive     = 0;
    flags->isAdvectionDampingYActive     = 0;
    flags->isCanopyActive                = 0;
    flags->isConcurrentPrecursorActive   = 0;
    flags->isPvCatalystActive            = 0;
    flags->isGravityWaveModelingActive   = 0;
    flags->isNonInertialFrameActive      = 0;
    flags->isMeangradPForcingActive      = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-overset",          &(flags->isOversetActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-les",              &(flags->isLesActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-potentialT",       &(flags->isTeqnActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-abl",              &(flags->isAblActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-windplant",        &(flags->isWindFarmActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-ibm",              &(flags->isIBMActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-zDampingLayer",    &(flags->isZDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-xDampingLayer",    &(flags->isXDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-yDampingLayer",    &(flags->isYDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-canopy",           &(flags->isCanopyActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-adjustTimeStep",   &(flags->isAdjustableTime), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-pvCatalyst",       &(flags->isPvCatalystActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-kLeftRayleigh",    &(flags->isKLeftRayleighDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-kRightRayleigh",   &(flags->isKRightRayleighDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-advectionDampingX",&(flags->isAdvectionDampingXActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-advectionDampingY",&(flags->isAdvectionDampingYActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-nonInertial",      &(flags->isNonInertialFrameActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-agwModeling",      &(flags->isGravityWaveModelingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-meanGradPForce",   &(flags->isMeangradPForcingActive ), PETSC_NULL);

    // do some checks
    if(flags->isZDampingActive || flags->isXDampingActive || flags->isYDampingActive || flags->isKLeftRayleighDampingActive || flags->isKRightRayleighDampingActive || flags->isAdvectionDampingXActive)
    {
        // y damping only goes with x damping
        if(flags->isYDampingActive && !flags->isXDampingActive)
        {
            char error[512];
            sprintf(error, "xDampingLayer is required for yDampingLayer");
            fatalErrorInFunction("SetSimulationFlags", error);
        }

        if(!flags->isAblActive)
        {
            char error[512];
            sprintf(error, "abl flag is required for zDampingLayer, xDampingLayer, yDampingLayer, kLeftRayleigh, kRightRayleigh, advectionDamping");
            fatalErrorInFunction("SetSimulationFlags", error);
        }

        if(!flags->isTeqnActive)
        {
            char error[512];
            sprintf(error, "potentialT flag is required for zDampingLayer, xDampingLayer or yDampingLayer");
            fatalErrorInFunction("SetSimulationFlags", error);
        }

        // print a suggestion when running xdamping layer or kleft/kright rayleigh without advection damping
        if
        (
            (
                flags->isXDampingActive ||
                flags->isKLeftRayleighDampingActive ||
                flags->isKRightRayleighDampingActive
            )
            && !flags->isAdvectionDampingXActive
        )
        {
            char warning[512];
            sprintf(warning, "advection damping is suggested with inlet/outlet fringe layers and it is currently deactivated");
            warningInFunction("SetSimulationFlags", warning);
        }
    }

    // read acquisition flags
    PetscInt isProbesActive         = 0;
    PetscInt isSectionsActive       = 0;
    PetscInt isAverageABLActive     = 0;
    PetscInt isAverage3LMActive     = 0;
    PetscInt isPhaseAveragingActive = 0;
    PetscInt isAveragingActive      = 0;
    PetscInt isKEBudgetsActive      = 0;
    PetscInt isQCritActive          = 0;
    PetscInt isVgtQgActive          = 0;
    PetscInt isVgtRgActive          = 0;
    PetscInt isVgtQsActive          = 0;
    PetscInt isVgtRsActive          = 0;
    PetscInt isVgtQrActive          = 0;
    PetscInt isL2CritActive         = 0;
    PetscInt isSourcesActive        = 0;
    PetscInt isPerturbABLActive     = 0;
    PetscInt PvCatalystActive       = 0;
    PetscInt continuityActive       = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-probes",           &isProbesActive,         PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sections",         &isSectionsActive,       PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averageABL",       &isAverageABLActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-average3LM",       &isAverage3LMActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averaging",        &isAveragingActive,      PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-keBudgets",        &isKEBudgetsActive,      PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-phaseAveraging",   &isPhaseAveragingActive, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQ",         &isQCritActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeL2",        &isL2CritActive,         PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQg",        &isVgtQgActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeRg",        &isVgtRgActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQs",        &isVgtQsActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeRs",        &isVgtRsActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQr",        &isVgtQrActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeSources",   &isSourcesActive,        PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-perturbABL",       &isPerturbABLActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-pvCatalyst",       &PvCatalystActive,       PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeContinuity",&continuityActive,       PETSC_NULL);

    flags->isAquisitionActive
    =
    PetscMin(
    (
        isProbesActive + isSectionsActive + isAverageABLActive + isAverage3LMActive + isKEBudgetsActive +
        isAveragingActive + isPhaseAveragingActive + isQCritActive + isVgtQgActive + isVgtRgActive + isVgtQsActive + isVgtRsActive + isVgtQrActive + isL2CritActive + isSourcesActive + isPerturbABLActive + PvCatalystActive + continuityActive),
        1
    );

    return(0);

}

//***************************************************************************************************************//

PetscErrorCode SetSimulationInfo(simInfo_ *info)
{
    info->nDomains = 1;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetDomainsAndAllocate(domain_ **domainAddr, flags_ *flags, simInfo_ *info)
{
    // set the overset pointer
    if(flags->isOversetActive)
    {
        readDictInt("overset/oversetInput.dat", "MeshTotal", &(info->nDomains));
        readDictInt("overset/oversetInput.dat", "numHoleObjects", &(info->nHoleRegions));

        // allocate memory for the number of domains
        *domainAddr = new domain_[info->nDomains];

        domain_ *domain = *domainAddr;

        for(PetscInt d=0; d<info->nDomains; d++)
        {
            char userName[256];
            sprintf(userName, "Mesh%ld", d);

            domain[d].info  = *info;
            domain[d].flags = *flags;

            //set domain specific flags
            if(domain[d].flags.isIBMActive)
            {
                readSubDictInt("overset/oversetInput.dat", userName,"ibm", &(domain[d].flags.isIBMActive));
            }

            if(domain[d].flags.isWindFarmActive)
            {
                readSubDictInt("overset/oversetInput.dat", userName,"windplant", &(domain[d].flags.isWindFarmActive));
            }

            // allocate memory for domain objects
            SetDomainMemory(&(domain[d]));

            // allocate memory for overset
            domain[d].os = new overset_;

            overset_ *os = domain[d].os;


            readSubDictIntArray("overset/oversetInput.dat", userName, "parentMesh", os->parentMeshId);
            readSubDictIntArray("overset/oversetInput.dat", userName, "childMesh",  os->childMeshId);

            // set mesh name
            readSubDictWord("overset/oversetInput.dat", userName, "name",  &(domain[d].mesh->meshName));

        }
    }
    else
    {
        // allocate memory for the number of domains
        *domainAddr = new domain_;

        domain_ *domain = *domainAddr;

        domain->info  = *info;
        domain->flags = *flags;

        // allocate memory for domain objects
        SetDomainMemory(domain);

        // set mesh name
        domain->mesh->meshName = ".";
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ReadTimeControls(clock_ *clock)
{
    readDictWord  ("control.dat", "-startFrom", &(clock->startFrom));
    readDictDouble("control.dat", "-timeStep",  &(clock->dt));
    readDictDouble("control.dat", "-cfl",       &(clock->cfl));
    readDictDouble("control.dat", "-endTime",   &(clock->endTime));

    clock->it      = 0;
    clock->itStart = 0;
    clock->dxMin   = 1e10;
    clock->startDt = clock->dt;

    // get time precision (optional)
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    clock->timePrecision = 2;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-timePrecision",(PetscInt*)&(clock->timePrecision), PETSC_NULL);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetStartTime(clock_ *clock, domain_ *domain, simInfo_ *info)
{
    // get start time for each doman
    std::vector<PetscReal> startTime(info->nDomains);

    for(PetscInt d=0; d<info->nDomains; d++)
    {
        mesh_ *mesh = domain[d].mesh;

        if(clock->startFrom == "startTime")
        {
            readDictDouble("control.dat", "-startTime", &(startTime[d]));
        }
        else if(clock->startFrom == "latestTime")
        {
            word location = "./fields/" + mesh->meshName;
            std::vector<PetscReal>        timeSeries;
            PetscInt                      ntimes;
            getTimeList(location.c_str(), timeSeries, ntimes);

            startTime[d] = timeSeries[ntimes-1];

            std::vector<PetscReal> ().swap(timeSeries);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown parameter startFrom %s, known parameters are \n-startTime\n-latestTime", clock->startFrom.c_str());
            fatalErrorInFunction("SetStartTime",  error);
        }
    }

    // check that all start times are the same
    for(PetscInt d=1; d<info->nDomains; d++)
    {
        if(startTime[d] != startTime[0])
        {
            char error[512];
            sprintf(error, "available startTime value is not equal for each domain\n");
            fatalErrorInFunction("SetStartTime",  error);
        }
    }

    clock->startTime = startTime[0];

    clock->time    = clock->startTime;

    // clea nmemory
    std::vector<PetscReal> ().swap(startTime);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ReadPhysicalConstants(domain_ *domain)
{
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    // read kinematic viscosity
    readDictDouble("control.dat", "-nu", &(domain->constants.nu));

    // read flow density
    if(domain->flags.isWindFarmActive || domain->flags.isIBMActive)
    {
      readDictDouble("control.dat", "-rho", &(domain->constants.rho));
    }
    else
    {
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL,"-rho", &(domain->constants.rho), PETSC_NULL);
    }

    // read Prandtl number
    if(domain->flags.isTeqnActive)
    {
        readDictDouble("control.dat", "-Pr",  &(domain->constants.Pr));

        // read reference temperature if abl is not active
        if(!domain->flags.isAblActive)
        {
            readDictDouble("control.dat", "-tRef",  &(domain->constants.tRef));
        }
    }
    else
    {
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL,"-Pr", &(domain->constants.Pr), PETSC_NULL);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetDomainMemory(domain_ *domain)
{
    // initialize pointers to null
    domain->clock       = NULL;
    domain->mesh        = NULL;
    domain->io          = NULL;
    domain->ueqn        = NULL;
    domain->peqn        = NULL;
    domain->teqn        = NULL;
    domain->les         = NULL;
    domain->ibm         = NULL;
    domain->abl         = NULL;
    domain->os          = NULL;
    domain->farm        = NULL;
    domain->acquisition = NULL;

    // allocate pointers based on flags
    domain->mesh = new mesh_;
    domain->io   = new io_;
    domain->ueqn = new ueqn_;
    domain->peqn = new peqn_;

    if(domain->flags.isTeqnActive)       domain->teqn        = new teqn_;
    if(domain->flags.isLesActive)        domain->les         = new les_;
    if(domain->flags.isAblActive)        domain->abl         = new abl_;
    if(domain->flags.isWindFarmActive)   domain->farm        = new farm_;
    if(domain->flags.isAquisitionActive) domain->acquisition = new acquisition_;
    if(domain->flags.isIBMActive)        domain->ibm         = new ibm_;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetAccessPointers(domain_ *domain)
{
    // set pointer to clock
    domain->access.clock     = domain->clock;

    // set access to simulation flags
    domain->access.flags     = &(domain->flags);

    // set access to simulation info
    domain->access.info      = &(domain->info);

    // set access to domain id
    domain->access.domainID  = &(domain->domainID);

    // set access to physical constants
    domain->access.constants = &(domain->constants);

    // mesh two way access
    domain->access.mesh      = domain->mesh;
    domain->mesh->access     = &(domain->access);

    // io settings two way access
    domain->access.io        = domain->io;
    domain->io->access       = &(domain->access);

    // u equation two way access
    domain->access.ueqn      = domain->ueqn;
    domain->ueqn->access     = &(domain->access);

    // p equation two way access
    domain->access.peqn      = domain->peqn;
    domain->peqn->access     = &(domain->access);

    if(domain->flags.isTeqnActive)
    {
        // t equation two way access
        domain->access.teqn  = domain->teqn;
        domain->teqn->access = &(domain->access);
    }

    if(domain->flags.isLesActive)
    {
        // les model two way access
        domain->access.les       = domain->les;
        domain->les->access      = &(domain->access);
    }

    if(domain->flags.isWindFarmActive)
    {
        // farm model two way access
        domain->access.farm       = domain->farm;
        domain->farm->access      = &(domain->access);
    }

    if(domain->flags.isOversetActive)
    {
        // les model two way access
        domain->access.os         = domain->os;
        domain->os->access        = &(domain->access);
    }

    if(domain->flags.isAblActive)
    {
        // les model two way access
        domain->access.abl        = domain->abl;
        domain->abl->access       = &(domain->access);
    }

    if(domain->flags.isAquisitionActive)
    {
        // acquisition two way access
        domain->access.acquisition  = domain->acquisition;
        domain->acquisition->access = &(domain->access);
    }

    if(domain->flags.isIBMActive)
    {
        // ibm two way access
        domain->access.ibm  = domain->ibm;
        domain->ibm->access = &(domain->access);
    }

    return(0);
}

//***************************************************************************************************************//
