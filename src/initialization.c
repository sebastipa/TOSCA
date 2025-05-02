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
#include "include/vents.h"

//***************************************************************************************************************//

PetscErrorCode PrintOkWindLogo()
{
    PetscPrintf(PETSC_COMM_WORLD,"\n\n   |========================|      This file is part of TOSVA.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |            S           |      TOSVA is free software: you can redistribute it and/or modify it\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |         O  o  V        |      under the terms of the GNU General Public License as published by\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |       T o  |  o A      |      the Free Software Foundation, either version 3 of the License, or\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |          \\ | /         |      (at your option) any later version.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |       o -- O -- o      |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |          / | \\         |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |         o  |  o        |      TOSVA is distributed in the hope that it will be useful, but WITHOUT \n");
    PetscPrintf(PETSC_COMM_WORLD,"   |            o           |      ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |                        |      FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |  copyright: OK-CFD Lab |      for more details.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |  authors: Stipa - Ajay |\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |  Haji - Christianson   |      You should have received a copy of the GNU General Public License\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |========================|      along with TOSCA.  If not, see <http://www.gnu.org/licenses/>.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   | Toolbox fOr Simulating |      This software was developed by Sebastiano Stipa, Arjun Ajay, Mohammad\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |       Ventilation      |      Hadi and Cole Chrsitianson during their time at the Okanagan CFD-LAB.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |           and          |      It is a LES code in generalized curvilinear coordinates for the \n");
    PetscPrintf(PETSC_COMM_WORLD,"   |         Aerosols       |      simulation of ventilation and aerosols.\n");
    PetscPrintf(PETSC_COMM_WORLD,"   |========================|\n\n\n\n");

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

    // set simulation start time
    SetStartTime(clock, domain,info);

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

        // set wall models
        SetWallModels(domain[d].ueqn);

        // initialize ibm
        InitializeIBM(domain[d].ibm);

        // set inflow functions
        SetInflowFunctions(domain[d].mesh);

        // initialize i/o controls and initialization type
        InitializeIO(domain[d].io);

        // initialize ABL parameters
        InitializeABL(domain[d].abl);

        // momentum equation initialize
        InitializeUEqn(domain[d].ueqn);

        // poisson equation initialize
        InitializePEqn(domain[d].peqn);

        // temperature equation initialize
        InitializeTEqn(domain[d].teqn);

        InitializeSMObject(domain[d].smObject);

        // scalarMoment equations initialize
        for  (PetscInt  ii=0; ii < flags->isScalarMomentsActive; ii++)
        {
            InitializeSM(domain[d].smObject->sm[ii]);
        }

        // LES model initialize
        InitializeLES(domain[d].les);

        // initialize wind farm
        InitializeWindFarm(domain[d].farm);

        PetscPrintf(PETSC_COMM_WORLD, "------------------------------------------------------------------------\n");
    }

    PetscReal timeStart, timeEnd;

    // sync processors
    MPI_Barrier(PETSC_COMM_WORLD);

    PetscTime(&timeStart);

    // initialize acquisitions
    InitializeAcquisition(domain);

    // sync processors
    MPI_Barrier(PETSC_COMM_WORLD);

    PetscTime(&timeEnd);

    PetscPrintf(PETSC_COMM_WORLD, "Acquisition initialization time = %lf s\n", timeEnd - timeStart);

    // initialize vents
    initializeVents(domain);

    // Set the initial field
    SetInitialField(domain);

    // initialize overset
    InitializeOverset(domain);

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
    flags->isAdvectionDampingActive      = 0;
    flags->isCanopyActive                = 0;
    flags->isConcurrentPrecursorActive   = 0;
    flags->isPvCatalystActive            = 0;
    flags->isGravityWaveModelingActive   = 0;
    flags->isNonInertialFrameActive      = 0;
    flags->isVentsActive                 = 0;
    flags->isScalarMomentsActive         = 0;
    flags->isCoagSourceActive            = 0;
    flags->isGroSourceActive             = 0;
    flags->isDepoSourceActive            = 0;
    flags->isSediFluxActive              = 0;
    flags->isDeviFluxActive              = 0;

    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-overset",         &(flags->isOversetActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-les",             &(flags->isLesActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-potentialT",      &(flags->isTeqnActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-abl",             &(flags->isAblActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-windplant",       &(flags->isWindFarmActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-ibm",             &(flags->isIBMActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-zDampingLayer",   &(flags->isZDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-xDampingLayer",   &(flags->isXDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-yDampingLayer",   &(flags->isYDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-canopy",          &(flags->isCanopyActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-adjustTimeStep",  &(flags->isAdjustableTime), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-pvCatalyst",      &(flags->isPvCatalystActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-kLeftRayleigh",   &(flags->isKLeftRayleighDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-kRightRayleigh",  &(flags->isKRightRayleighDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-advectionDamping",&(flags->isAdvectionDampingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-nonInertial",     &(flags->isNonInertialFrameActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-agwModeling",     &(flags->isGravityWaveModelingActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-vents", &(flags->isVentsActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-scalarMoments", &(flags->isScalarMomentsActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-coagSource", &(flags->isCoagSourceActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-groSource", &(flags->isGroSourceActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-depoSource", &(flags->isDepoSourceActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-sediFlux", &(flags->isSediFluxActive), PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-deviFlux", &(flags->isDeviFluxActive), PETSC_NULL);

	// do some checks
	if(flags->isZDampingActive || flags->isXDampingActive || flags->isYDampingActive || flags->isKLeftRayleighDampingActive || flags->isKRightRayleighDampingActive || flags->isAdvectionDampingActive)
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
            && !flags->isAdvectionDampingActive
        )
        {
            char warning[512];
            sprintf(warning, "advection damping is suggested with inlet/outlet fringe layers and it is currently deactivated");
            warningInFunction("SetSimulationFlags", warning);
        }
	}

    if(flags->isScalarMomentsActive != 0 && flags->isScalarMomentsActive != 6)
    {
        printf("\n.......%li\n", flags->isScalarMomentsActive);
        char error[512];
        sprintf(error, "Scalar moment flag must be 0 or 6");
        fatalErrorInFunction("SetSimulationFlags", error);
    }

    // read acquisition flags
    PetscInt isProbesActive         = 0;
    PetscInt isSectionsActive       = 0;
    PetscInt isAverageABLActive     = 0;
    PetscInt isAverage3LMActive     = 0;
    PetscInt isPhaseAveragingActive = 0;
    PetscInt isAveragingActive      = 0;
    PetscInt isTKEActive            = 0;
    PetscInt isKEBudgetsActive      = 0;
    PetscInt isQCritActive          = 0;
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
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeSources",   &isSourcesActive,        PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-perturbABL",       &isPerturbABLActive,     PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-pvCatalyst",       &PvCatalystActive,       PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-computeContinuity",&continuityActive,       PETSC_NULL);
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-tke",           &isTKEActive,            PETSC_NULL);

    flags->isAquisitionActive
    =
    PetscMin(
    (
        isProbesActive + isSectionsActive + isAverageABLActive + isAverage3LMActive + isKEBudgetsActive +
        isAveragingActive + isPhaseAveragingActive + isQCritActive + isL2CritActive + isSourcesActive + isPerturbABLActive + PvCatalystActive + continuityActive),
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
        readDictDouble("control.dat", "-kCon",  &(domain->constants.kCon));
        readDictDouble("control.dat", "-alpha",  &(domain->constants.alpha));

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

    // read Scmidt number
    if(domain->flags.isScalarMomentsActive)
    {
        readDictDouble("control.dat", "-Sc",  &(domain->constants.Sc));
        readDictDouble("control.dat", "-ScT",  &(domain->constants.ScT));

        // read reference temperature if abl is not active
        if(!domain->flags.isAblActive)
        {
            readDictDouble("control.dat", "-tRef",  &(domain->constants.tRef));
        }
    }
    else
    {
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL,"-Sc", &(domain->constants.Sc), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL,"-ScT", &(domain->constants.ScT), PETSC_NULL);
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
    domain->vents       = NULL;
    domain->smObject    = NULL;

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
    if(domain->flags.isVentsActive)      domain->vents       = new vents_;
    if(domain->flags.isScalarMomentsActive)
    {
        domain->smObject       = new SMObj_;

        // allocate memory for each Scalar Moment object and the weights and abscissi
        PetscMalloc(domain->flags.isScalarMomentsActive * sizeof(sm_*), &(domain->smObject->sm));
        PetscMalloc(domain->flags.isScalarMomentsActive * sizeof(WandA*), &(domain->smObject->weightAbsc));

        for  (PetscInt  i=0; i < domain->flags.isScalarMomentsActive; i++)
        {
            PetscMalloc(sizeof(sm_), &(domain->smObject->sm[i]));
            PetscMalloc(sizeof(access_), &(domain->smObject->sm[i]->access));
        }

        for  (PetscInt  i=0; i < 3; i++)
        {
            PetscMalloc(sizeof(WandA), &(domain->smObject->weightAbsc[i]));
        }
    }


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

    if(domain->flags.isVentsActive)
    {
        domain->access.vents  = domain->vents;
        domain->vents->access = &(domain->access);
    }

    if(domain->flags.isScalarMomentsActive)
    {
        domain->access.smObject = domain->smObject; //oly 1 way access needed for smObject since each sm has 2 way accesss
        //domain->smObject->access = &(domain->access);

        for  (PetscInt  i=0; i < domain->flags.isScalarMomentsActive; i++)
        {
            domain->smObject->sm[i]->access->smObject = domain->smObject;
            domain->smObject->sm[i]->access = &(domain->access);
        }

    }

    return(0);
}

//***************************************************************************************************************//
