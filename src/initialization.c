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
    PetscPrintf(PETSC_COMM_WORLD,"\n\n   |==================================================|      This file is part of TOSCA.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |               \\\\        T                        |      TOSCA is free software: you can redistribute it and/or modify it\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                \\\\        O                       |      under the terms of the GNU General Public License as published by\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                 \\\\    \\\\  S                      |      the Free Software Foundation, either version 3 of the License, or\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                  \\\\    \\\\  C                     |      (at your option) any later version.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                   \\\\    \\\\  A                    |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                    O ====\\\\=======               |      TOSCA is distributed in the hope that it will be useful, but WITHOUT \n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                   //|     \\\\                     |      ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                  //||      O =============       |      FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                 // ||     //| v0122              |      for more details.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                //  ||    //||                    |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |               //   ||   // ||                    |      You should have received a copy of the GNU General Public License\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                    ||  //  ||                    |      along with TOSCA.  If not, see <http://www.gnu.org/licenses/>.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                    || //   ||                    |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                    ||      ||                    |      This software was developed by Sebastiano Stipa, Arjun Ajay and Mohammad\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |____________________||______||____________________|      Hadi during their time at the Okanagan CFD-LAB. It is a LES code in\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |            copyright : OK-CFD Lab                |      generalized curvilinear coordinates for the simulation of wind plants\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |            authors   : Stipa - Ajay - Haji       |      immersed in the atmospheric boundary layer. Momentum and pressure equations\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |==================================================|      are solved using operators splitting. If activated, potential temperature\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |  Toolbox for Stratified Convective Atmospheres   |      transport is also solved. Advanced actuator disk/line models and IBM with\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |__________________________________________________|      the moving least squares techinque are also implemented.\n\n\n\n");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode PrintNumberOfProcs()
{
    PetscMPIInt nProcs;

    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);
    PetscPrintf(PETSC_COMM_WORLD,"Simulation running with %ld processors\n", nProcs);
    PetscPrintf(PETSC_COMM_WORLD,"-----------------------------------\n\n");

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

    for(PetscInt d=0; d<info->nDomains; d++)
    {
        PetscPrintf(PETSC_COMM_WORLD, "\nDomain %ld\n", d);
        PetscPrintf(PETSC_COMM_WORLD, "--------\n");

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

        // set inflow functions
        SetInflowFunctions(domain[d].mesh);

        // set wall models
        SetWallModels(domain[d].ueqn);

        // initialize i/o controls and initialization type
        InitializeIO(domain[d].io);

        // initialize ABL parameters
        InitializeABL(domain[d].abl);

        // initialize ibm
        InitializeIBM(domain[d].ibm);

        // momentum equation initialize
        InitializeUEqn(domain[d].ueqn);

        // poisson equation initialize
        InitializePEqn(domain[d].peqn);

        // temperature equation initialize
        InitializeTEqn(domain[d].teqn);

        // LES model initialize
        InitializeLES(domain[d].les);

        // initialize wind farm
        InitializeWindFarm(domain[d].farm);

        PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");
    }

    // initialize acquisitions
    InitializeAcquisition(domain);

    // Set the initial field
    SetInitialField(domain);

    // initialize overset
    InitializeOverset(domain);

    // initialize ibm
    InitializeIBMInterpolation(domain);

    return(0);

}

//***************************************************************************************************************//

PetscErrorCode SetSimulationFlags(flags_ *flags)
{
    // read from control file
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    flags->isOversetActive    = 0;
    flags->isLesActive        = 1;
    flags->isTeqnActive       = 0;
    flags->isWindFarmActive   = 0;
    flags->isAquisitionActive = 0;
    flags->isAblActive        = 0;
    flags->isIBMActive        = 0;
    flags->isZDampingActive   = 0;
    flags->isXDampingActive   = 0;
    flags->isSideForceActive  = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-overset",       &(flags->isOversetActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-les",           &(flags->isLesActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-potentialT",    &(flags->isTeqnActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-abl",           &(flags->isAblActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-windplant",     &(flags->isWindFarmActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-ibm",           &(flags->isIBMActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-zDampingLayer", &(flags->isZDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-xDampingLayer", &(flags->isXDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sideForce",     &(flags->isSideForceActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-adjustTimeStep",&(flags->isAdjustableTime), PETSC_NULL);

    // read acquisition flags
    PetscInt isProbesActive         = 0;
    PetscInt isSectionsActive       = 0;
    PetscInt isAverageABLActive     = 0;
    PetscInt isAverage3LMActive     = 0;
    PetscInt isPhaseAveragingActive = 0;
    PetscInt isAveragingActive      = 0;
    PetscInt isKEBudgetsActive      = 0;
    PetscInt isQCritActive          = 0;
    PetscInt isL2CritActive         = 0;
    PetscInt isSourcesActive        = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-probes",        &isProbesActive,         PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sections",      &isSectionsActive,       PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averageABL",    &isAverageABLActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-average3LM",    &isAverage3LMActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-averaging",     &isAveragingActive,      PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-keBudgets",     &isKEBudgetsActive,      PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-phaseAveraging",&isPhaseAveragingActive, PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeQ",      &isQCritActive,          PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeL2",     &isL2CritActive,         PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeSources",&isSourcesActive,        PETSC_NULL);

    flags->isAquisitionActive
    =
    PetscMin(
    (
        isProbesActive + isSectionsActive + isAverageABLActive + isAverage3LMActive + isKEBudgetsActive +
        isAveragingActive + isPhaseAveragingActive + isQCritActive + isL2CritActive + isSourcesActive),
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
        readDictInt("Overset/OversetInput.dat", "MeshTotal", &(info->nDomains));

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
                readSubDictInt("Overset/OversetInput.dat", userName,"ibm", &(domain[d].flags.isIBMActive));
            }

            if(domain[d].flags.isWindFarmActive)
            {
                readSubDictInt("Overset/OversetInput.dat", userName,"windplant", &(domain[d].flags.isWindFarmActive));
            }

            // allocate memory for domain objects
            SetDomainMemory(&(domain[d]));

            // allocate memory for overset
            domain[d].os = new overset_;

            overset_ *os = domain[d].os;


            readSubDictIntArray("Overset/OversetInput.dat", userName, "parentMesh", os->parentMeshId);
            readSubDictIntArray("Overset/OversetInput.dat", userName, "childMesh",  os->childMeshId);

            // set mesh name
            readSubDictWord("Overset/OversetInput.dat", userName, "name",  &(domain[d].mesh->meshName));

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
    readDictDouble("control.dat", "-startTime", &(clock->startTime));
    readDictDouble("control.dat", "-timeStep",  &(clock->dt));
    readDictDouble("control.dat", "-cfl",       &(clock->cfl));
    readDictDouble("control.dat", "-endTime",   &(clock->endTime));

    clock->time    = clock->startTime;
    clock->it      = 0;
    clock->itStart = 0;
    clock->startDt = clock->dt;

    // get time precision (optional)
    PetscOptionsInsertFile(PETSC_COMM_WORLD, PETSC_NULL, "control.dat", PETSC_TRUE);

    clock->timePrecision = 2;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-timePrecision",(PetscInt*)&(clock->timePrecision), PETSC_NULL);

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
        domain->access.ibm  = domain->ibm;
        domain->ibm->access = &(domain->access);
    }

    return(0);
}

//***************************************************************************************************************//
