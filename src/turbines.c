//! \file  turbines.c
//! \brief Contains top to bottom level routines for the wind farm modeling.

// public definitions
#include "include/turbines.h"  

// private definitions 
#include "turbines/turbines_io.h"
#include "turbines/turbines_openfast.h"
#include "turbines/turbines_parallel.h"
#include "turbines/turbines_sampling.h"
#include "turbines/turbines_force.h"
#include "turbines/turbines_phys.h"

// Our ADM/ALM models are detailed actuator disk/line
// models with blade pitch and generator speed controls, blade and airfoil properties
// linearly interpolated where necessary. UniformADM is a simple AD version, where the
// thrust is weighted based on the rotor mesh element area. Yaw control
// can be activated for all models and is based upon 1D upstream velocity sampling. Turbines
// rotate with a prescribed nacelle rotation speed.
// We have fully optimized parallel communication: each turbine, tower and up-sampling
// rig (it s just a rotor defined 2.5D upstream of each rotor for velocty sampling)
// has its own communicator, and turbine solution is done in parallel by each communicator.
// Each turbine has an influence sphere of mesh cells, where the rotor can possibly
// be located in any yaw condition, and those cells define the turbine communicator based
// on the processor which they belong to. In this way communications are performed between
// max 8 processors (that is the maximum number of processors that can control a wind
// turbine, most of the time is 2 or 1). Turbine I/O is performed also in parallel:
// the master rank of each communicator is responsible for writing down turbine log
// and checkpoint file. In normal operation (no debug), a global communication between
// all turbines is never performed. There is only 1 communication of type MPI_Reduce
// for wind turbine yaw syncronization on the master node, which will then write
// the mesh to file. That is done only at global write time (not at turbine log write time).
// The only function call for turbine update is "PetscErrorCode UpdateWindTurbines(farm_ *farm)"".

PetscErrorCode UpdateWindTurbines(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    PetscReal ts, te;

    PetscTime(&ts);

    if(!farm->dbg) PetscPrintf(mesh->MESH_COMM,"Updating wind turbines, ");

    
#if USE_OPENFAST

    // find which force point this processor controls 
    findControlledVelPointsOpenFAST(farm); 

    // compute and send velocities to openfast 
    computeWindVectorsRotorOpenFAST(farm); 

    // advance openfast 
    stepOpenFAST(farm); 

    // find which AD points this processor controls
    findControlledPointsRotor(farm); 

    // find which sample points this processor controls
    findControlledPointsSample(farm); 

    // compute wind velocity at the sample mesh points
    computeWindVectorsSample(farm); 

    // compute wind velocity at the rotor mesh points
    computeWindVectorsRotor(farm); 

    // compute wind velocity at the tower mesh points
    computeWindVectorsTower(farm); 

    // compute wind velocity at the nacelle mesh point
    computeWindVectorsNacelle(farm); 

    // get force from openfast at current time 
    computeBladeForceOpenFAST(farm); 

#else 
    // solve rotor dynamics and compute rot speeds
    computeRotSpeed(farm);

    // apply the generator speed controller
    controlGenSpeed(farm);

    // apply the pitch PID controller
    controlBldPitch(farm);

    // apply the wind farm controller
    windFarmControl(farm);

    // apply the yaw controller
    controlNacYaw(farm);

    // find which AD points this processor controls
    findControlledPointsRotor(farm);

    // find which sample points this processor controls
    findControlledPointsSample(farm);

    // compute wind velocity at the sample mesh points
    computeWindVectorsSample(farm);

    // compute wind velocity at the rotor mesh points
    computeWindVectorsRotor(farm);

    // compute wind velocity at the tower mesh points
    computeWindVectorsTower(farm);

    // compute wind velocity at the nacelle mesh point
    computeWindVectorsNacelle(farm);
    
    // compute aerodynamic forces at the turbine mesh points
    computeBladeForce(farm);

#endif

    // project the wind turbine forces on the background mesh
    projectBladeForce(farm); 

    // project the tower forces on the background mesh
    projectTowerForce(farm); 

    // project the nacelle forces on the background mesh
    projectNacelleForce(farm); 

    // transform forces from cartesian to contravariant
    bodyForceCartesian2Contravariant(farm); 

    // write output if runTimeWrite == 1
    windTurbinesWrite(farm); 

    PetscTime(&te);

    if(!farm->dbg) PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %lf\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeWindFarm(farm_ *farm)
{
    if(farm != NULL)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscPrintf(mesh->MESH_COMM, "Initializing wind turbines...\n");

        // read farm propeties and initialize all parameters
        PetscPrintf(mesh->MESH_COMM, "   reading wind farm properties...\n");
        readFarmProperties(farm);
        PetscInt nT = farm->size;

        PetscPrintf(mesh->MESH_COMM, "   reading actuator model parameters...\n");
        for(PetscInt t=0; t<nT; t++)
        {
            // read actuator model parameters (this needs to be done before OpenFAST)
            readActuatorModelParameters((*farm->turbineModels[t]), farm->wt[t], mesh->meshName); 

            // initialize the tower model
            if(farm->wt[t]->includeTwr)
            {
                initTwrModel(farm->wt[t], farm->base[t]);
            }

            // initialize the tower model
            if(farm->wt[t]->includeNacelle)
            {
                initNacModel(farm->wt[t], farm->base[t]);
            }

            // initialize upstream velocity sampling rigs
            initSamplePoints(farm->wt[t], mesh->meshName);
        }

        PetscPrintf(mesh->MESH_COMM, "   pre-calculating influence sphere cells...\n");

        // initialize turbine and tower sphere cells and see what processor controls which wind turbine
        initControlledCells(farm);

        // ensure best practices with respect to mesh size and turbine dimension are followed
        checkTurbineMesh(farm);

        // initialize sample sphere cells and see what processor controls which sampling rig
        initSampleControlledCells(farm);

        // determine which points in each tower are controlled by which processor (also finds closest cell)
        findControlledPointsTower(farm);

        // determine which points in each nacelle are controlled by which processor (also finds closest cell)
        findControlledPointsNacelle(farm);

        // initialize OpenFAST 
        #if USE_OPENFAST
        PetscPrintf(mesh->MESH_COMM, "   sending inputs to OpenFAST...\n");
        initOpenFAST(farm);
        #endif

        // initialize actuator model meshes
        PetscPrintf(mesh->MESH_COMM, "   initializing actuator model meshes...\n");
        for(PetscInt t=0; t<nT; t++)
        {
            if((*farm->turbineModels[t]) == "ADM")
            {
                // initialize actuator disk model parameters
                meshADM(farm->wt[t]);
            }
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // initialize uniform actuator disk model parameters
                meshUADM(farm->wt[t]);
            }
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // initialize actuator line model parameters
                meshALM(farm->wt[t]);
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // initialize actuator farm model parameters
                meshAFM(farm->wt[t]);
            }

            #if USE_OPENFAST
            // initialize processor-controlled points
            initControlledPointsOpenFAST(farm->wt[t]);

            // compute initial condition 
            stepZeroOpenFAST(farm);
            #endif
        }

        // read the checkpoint file if present and prepare wind turbines for restart
        windTurbinesReadCheckpoint(farm);

        PetscPrintf(mesh->MESH_COMM, "done\n\n");
    }

    return(0);
}
