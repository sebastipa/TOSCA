//! \file  ibm.c
//! \brief Contains Immersed boundary method function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/wallfunctions.h"
#include "include/ibmInput.h"

//***************************************************************************************************************//

PetscErrorCode InitializeIBM(ibm_ *ibm)
{
    if(ibm != NULL)
    {
        mesh_ *mesh = ibm->access->mesh;
        VecDuplicate(mesh->lNvert, &(ibm->lNvertFixed));  VecSet(ibm->lNvertFixed,  0.0);

        //read ibm input file
        readIBMProperties(ibm);

        //compute the node to element reverse connectivity
        nodeElementConnectivity(ibm);

        //set the wall model properties
        setIBMWallModels(ibm);

        //find the ibm cartesian bounding box
        findBodyBoundingBox(ibm);

        //find the search cell dimensions from the average cell size
        findSearchCellDim(ibm);

        //create the ibm search cell list - ibm elements in each search cell
        createSearchCellList(ibm);

        //compute element normals and check that they point outwards
        computeIBMElementNormal(ibm);

        // create half edge data structure
        createHalfEdgeDataStructure(ibm);

        // ibm search algorithm
        PetscPrintf(mesh->MESH_COMM, "IBM search algorithm...");
        ibmSearch(ibm);
        PetscPrintf(mesh->MESH_COMM, "done\n");

        MPI_Barrier(mesh->MESH_COMM);

        //check that the ibm object has been detected
        if(ibm->dbg) checkIBMexists(ibm);

        //create local list of ibm fluid cells
        findIBMFluidCells(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        //find the processors that have ibm body in it - to parallelize
        findIBMControlledProcs(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        if(ibm->computeForce)
        {
            // divide the ibm elements into ibm processors based on the closest ibm fluid to the ibm element normal projection
            initElementProjectionProcs(ibm);

            MPI_Barrier(mesh->MESH_COMM);

        }

        //divide the ibm elements into ibm processors based on the closest ibm fluid to the element center
        initElementProcs(ibm);

        //find the closest normal projection element to every ibm fluid cell
        if(ibm->wallShearOn)
        {
            findClosestIBMElement2Solid(ibm);
        }
        else 
        {
            findClosestIBMElement(ibm);
        }

        if (ibm->IBInterpolationModel == "CURVIB")
        {
            //find the intereption point on the background grid
            if(ibm->curvibType == "CurvibTriangular")
            {
                findInterceptionPoint(ibm);
            }
        }

        MPI_Barrier(mesh->MESH_COMM);

        if (ibm->IBInterpolationModel == "MLS")
        {
            //find the fluid nodes within the support radius
            findFluidSupportNodes(ibm);

            //find the solid nodes within the support radius
            findIBMMeshSupportNodes(ibm);

            MPI_Barrier(mesh->MESH_COMM);

        }

    }

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode UpdateIBM(ibm_ *ibm)
{
    PetscReal solutionTimeStart, solutionTimeEnd;

    PetscTime(&solutionTimeStart);
    mesh_ *mesh = ibm->access->mesh;

    //copy current nvert before next time step
    VecCopy(mesh->lNvert, mesh->lNvert_o);
    DMLocalToLocalBegin(mesh->da, mesh->lNvert_o, INSERT_VALUES, mesh->lNvert_o);
	DMLocalToLocalEnd(mesh->da, mesh->lNvert_o, INSERT_VALUES, mesh->lNvert_o);

    PetscReal iterationTimeStart, iterationTimeEnd;

    //search to update the ibm cells and support nodes if dynamic ibm
    if(ibm->dynamic)
    {
        // destroy the lists created if dynamic before next iteration
        if(ibm->access->clock->it > ibm->access->clock->itStart)
        {
            destroyLists(ibm);
        }

        //reset nvert values to 0, they will be recomputed during ibm search
        VecSet(mesh->lNvert,0.);

        MPI_Barrier(mesh->MESH_COMM);

        UpdateIBMesh(ibm);

        findBodyBoundingBox(ibm);

        findSearchCellDim(ibm);

        createSearchCellList(ibm);

        ibmSearch(ibm);

        if(ibm->dbg) checkIBMexists(ibm);

        findIBMFluidCells(ibm);

        IBMElementProcessorTransfer(ibm);

        findClosestIBMElement(ibm);

        if (ibm->IBInterpolationModel == "CURVIB")
        {
            //find the intereption point on the background grid
            if(ibm->curvibType == "CurvibTriangular")
            {
                findInterceptionPoint(ibm);
            }
        }

        if (ibm->IBInterpolationModel == "MLS")
        {

            findFluidSupportNodes(ibm);

            MPI_Barrier(mesh->MESH_COMM);

            findIBMMeshSupportNodes(ibm);

            MPI_Barrier(mesh->MESH_COMM);

        }

    }

    //interpolate the ibm fluid cells
    if (ibm->IBInterpolationModel == "MLS")
    {
        MLSInterpolation(ibm);

        MPI_Barrier(mesh->MESH_COMM);
    }
    else if (ibm->IBInterpolationModel == "CURVIB")
    {
        if(ibm->wallShearOn)
        {
            CurvibInterpolationInternalCell(ibm);
        }
        else 
        {
            if(ibm->curvibType == "CurvibTrilinear")
            {
                if(ibm->curvibOrder == "linear")
                {
                    CurvibInterpolation(ibm);
                }
                else if(ibm->curvibOrder == "quadratic")
                {
                    CurvibInterpolationQuadratic(ibm);
                }
                else
                {
                    char error[512];
                    sprintf(error, "wrong interpolation order chosen. Available options are linear and quadratic\n");
                    fatalErrorInFunction("readIBMProperties",  error);
                }
            }
            else if(ibm->curvibType == "CurvibTriangular")
            {
                CurvibInterpolationTriangular(ibm);
            }
            else
            {
                char error[512];
                sprintf(error, "wrong curvib interpolation type\n");
                fatalErrorInFunction("UpdateIBM", error);
            }
        }

        MPI_Barrier(mesh->MESH_COMM);
    }

    //update the boundary condition around ibm cells
    UpdateImmersedBCs(ibm);

    contravariantToCartesian(ibm->access->ueqn);

    MPI_Barrier(mesh->MESH_COMM);

    PetscTime(&solutionTimeEnd);
    PetscPrintf(mesh->MESH_COMM, "IBM Update time = %lf s\n", solutionTimeEnd - solutionTimeStart);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode setIBMWallModels(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    flags_        *flags = ibm->access->flags;

    word          fileNameU = "./boundary/" + mesh->meshName + "/U";
    word          fileNameT = "./boundary/" + mesh->meshName + "/T";

    PetscPrintf(mesh->MESH_COMM, "\nReading IBM Boundary conditions...");

    for (PetscInt i=0; i < ibm->numBodies; i++)
    {
        char objectName[256];
        sprintf(objectName, "object%ld", i);

        ibmObject   *ibmBody = ibm->ibmBody[i];

        // read wall model properties
        readSubDictWord("./IBM/IBMProperties.dat", objectName, "velocityBCSetType", &(ibmBody->uBCSetType));

        //copy the wall model properties from the U boundary condition
        if(ibmBody->uBCSetType == "matchUiLeft")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.iLeft == "velocityWallFunction")
            {
                ibmBody->velocityBC = mesh->boundaryU.iLeft;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(ibmBody->wallFunctionTypeU));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelU));

                //cabot model
                if(ibmBody->wallFunctionTypeU == -1)
                {
                    PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallModelU->wmCabot));

                    Cabot *wm = ibmBody->ibmWallModelU->wmCabot;

                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -3)
                {
                    PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelU->wmShumann));

                    Shumann *wm = ibmBody->ibmWallModelU->wmShumann;

                    readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else if (mesh->boundaryU.iLeft == "slip")
            {
                ibmBody->velocityBC = mesh->boundaryU.iLeft;
            }
            else if (mesh->boundaryU.iLeft == "noSlip")
            {
                ibmBody->velocityBC = mesh->boundaryU.iLeft;
            }
            else
            {
                char error[512];
                sprintf(error, "Cannot match BC set in iLeft to IBM boundary, set BC in IBM file\n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else if(ibmBody->uBCSetType == "matchUiRight")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.iRight == "velocityWallFunction")
            {
                ibmBody->velocityBC = mesh->boundaryU.iRight;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(ibmBody->wallFunctionTypeU));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelU));

                //cabot model
                if(ibmBody->wallFunctionTypeU == -1)
                {
                    PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallModelU->wmCabot));

                    Cabot *wm = ibmBody->ibmWallModelU->wmCabot;

                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -3)
                {
                    PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelU->wmShumann));

                    Shumann *wm = ibmBody->ibmWallModelU->wmShumann;

                    readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else if (mesh->boundaryU.iRight == "slip")
            {
                ibmBody->velocityBC = mesh->boundaryU.iRight;
            }
            else if (mesh->boundaryU.iRight == "noSlip")
            {
                ibmBody->velocityBC = mesh->boundaryU.iRight;
            }
            else
            {
                char error[512];
                sprintf(error, "Cannot match BC set in iRight to IBM boundary, set BC in IBM file\n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else if(ibmBody->uBCSetType == "matchUjLeft")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.jLeft == "velocityWallFunction")
            {
                ibmBody->velocityBC = mesh->boundaryU.jLeft;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(ibmBody->wallFunctionTypeU));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelU));

                //cabot model
                if(ibmBody->wallFunctionTypeU == -1)
                {
                    PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallModelU->wmCabot));

                    Cabot *wm = ibmBody->ibmWallModelU->wmCabot;

                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -3)
                {
                    PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelU->wmShumann));

                    Shumann *wm = ibmBody->ibmWallModelU->wmShumann;

                    readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else if (mesh->boundaryU.jLeft == "slip")
            {
                ibmBody->velocityBC = mesh->boundaryU.jLeft;
            }
            else if (mesh->boundaryU.jLeft == "noSlip")
            {
                ibmBody->velocityBC = mesh->boundaryU.jLeft;
            }
            else
            {
                char error[512];
                sprintf(error, "Cannot match BC set in jLeft to IBM boundary, set BC in IBM file\n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else if(ibmBody->uBCSetType == "matchUjRight")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.jRight == "velocityWallFunction")
            {
                ibmBody->velocityBC = mesh->boundaryU.jRight;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(ibmBody->wallFunctionTypeU));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelU));

                //cabot model
                if(ibmBody->wallFunctionTypeU == -1)
                {
                    PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallModelU->wmCabot));

                    Cabot *wm = ibmBody->ibmWallModelU->wmCabot;

                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -3)
                {
                    PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelU->wmShumann));

                    Shumann *wm = ibmBody->ibmWallModelU->wmShumann;

                    readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                    readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else if (mesh->boundaryU.jRight == "slip")
            {
                ibmBody->velocityBC = mesh->boundaryU.jRight;
            }
            else if (mesh->boundaryU.jRight == "noSlip")
            {
                ibmBody->velocityBC = mesh->boundaryU.jRight;
            }
            else
            {
                char error[512];
                sprintf(error, "Cannot match BC set in jRight to IBM boundary, set BC in IBM file\n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else if (ibmBody->uBCSetType == "setHere")
        {
            readSubDictWord("./IBM/IBMProperties.dat", objectName, "velocityBC", &(ibmBody->velocityBC));

            if (ibmBody->velocityBC == "velocityWallFunction")
            {

                //set theta wall function type
                readSubDictInt("./IBM/IBMProperties.dat", objectName, "wallFunctionTypeU", &(ibmBody->wallFunctionTypeU));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelU));

                if (ibmBody->wallFunctionTypeU == -1)
                {
                    PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallModelU->wmCabot));

                    Cabot *wm = ibmBody->ibmWallModelU->wmCabot;

                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "roughness", &(wm->roughness));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",  &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -3)
                {
                    PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelU->wmShumann));

                    Shumann *wm = ibmBody->ibmWallModelU->wmShumann;

                    readSubDictWord  ("./IBM/IBMProperties.dat", objectName, "uStarEval", &(wm->wfEvalType));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",     &(wm->kappa));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "thetaRef",  &(wm->thetaRef));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "roughness", &(wm->roughness));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "gammaM",    &(wm->gammaM));
                }
                else if (ibmBody->wallFunctionTypeU == -4)
                {
                    PetscMalloc(sizeof(PowerLawAPG), &(ibmBody->ibmWallModelU->wmPowerLawAPG));

                    PowerLawAPG *wm = ibmBody->ibmWallModelU->wmPowerLawAPG;

                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "roughness", &(wm->roughness));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",  &(wm->kappa));
                }
                else if (ibmBody->wallFunctionTypeU == -5)
                {
                    PetscMalloc(sizeof(LogLawAPG), &(ibmBody->ibmWallModelU->wmLogLawAPG));

                    LogLawAPG *wm = ibmBody->ibmWallModelU->wmLogLawAPG;

                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "roughness", &(wm->roughness));
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",  &(wm->kappa));
                }
                else 
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error); 
                }
            }

            else if (ibmBody->velocityBC == "slip")
            {

            }
            else if (ibmBody->velocityBC == "noSlip")
            {

            }
            else
            {
                char error[512];
                sprintf(error, "invalid velocity boundary condition chosen. Please use option velocityWallFunction or slip/noSlip conditions\n");
                fatalErrorInFunction("setIBMWallModels", error);
            }
        }

        if(flags->isTeqnActive)
        {
            //read temperature BC set type
            readSubDictWord("./IBM/IBMProperties.dat", objectName, "temperatureBCSetType", &(ibmBody->tBCSetType));

            //copy the temperature BC type from the temperature boundary file
            if(ibmBody->tBCSetType == "matchTiLeft")
            {
                if (mesh->boundaryT.iLeft == "thetaWallFunction")
                {
                    if(ibmBody->velocityBC != "velocityWallFunction" && ibmBody->wallFunctionTypeU !=-3)
                    {
                        char error[512];
                        sprintf(error, "thetaWallFunction can be used currently only for velocity BC = velocityWallFunction and type = -3 (Schumann Model)\n");
                        fatalErrorInFunction("SetIBMWallModels", error);   
                    }

                    ibmBody->tempBC = mesh->boundaryT.iLeft;

                    //set theta wall function type
                    readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(ibmBody->wallFunctionTypeT));

                    // allocate memory and connect pointers
                    PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelT));

                    if (ibmBody->wallFunctionTypeT == -3)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                        char error[512];
                        sprintf(error, "thetaWallFunction type -3 is not yet implemented for IBM\n");
                        fatalErrorInFunction("SetIBMWallModels", error); 
                    }
                    else if(ibmBody->wallFunctionTypeT == -2)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;

                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "invalid wall model chosen. Please use option -3 or -2\n");
                        fatalErrorInFunction("SetWallModels", error);
                    } 
                }
                else if(mesh->boundaryT.iLeft == "zeroGradient")
                {
                    ibmBody->tempBC = mesh->boundaryT.iLeft;
                }
                else if(mesh->boundaryT.iLeft == "fixedValue")
                {
                    ibmBody->tempBC = mesh->boundaryT.iLeft;
                    ibmBody->fixedTemp = mesh->boundaryT.iLval;
                }
                else 
                {
                    char error[512];
                    sprintf(error, "Cannot match BC set in iLeft to IBM boundary, set BC in IBM file\n");
                    fatalErrorInFunction("SetWallModels", error);
                }
            }
            else if(ibmBody->tBCSetType == "matchTiRight")
            {
                if (mesh->boundaryT.iRight == "thetaWallFunction")
                {
                    ibmBody->tempBC = mesh->boundaryT.iRight;

                    if(ibmBody->velocityBC != "velocityWallFunction" && ibmBody->wallFunctionTypeU !=-3)
                    {
                        char error[512];
                        sprintf(error, "thetaWallFunction can be used currently only for velocity BC = velocityWallFunction and type = -3 (Schumann Model)\n");
                        fatalErrorInFunction("SetIBMWallModels", error);   
                    }

                    //set theta wall function type
                    readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(ibmBody->wallFunctionTypeT));

                    // allocate memory and connect pointers
                    PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelT));

                    if (ibmBody->wallFunctionTypeT == -3)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                        char error[512];
                        sprintf(error, "thetaWallFunction type -3 is not yet implemented for IBM\n");
                        fatalErrorInFunction("SetIBMWallModels", error); 
                    }
                    else if(ibmBody->wallFunctionTypeT == -2)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;

                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "invalid wall model chosen. Please use option -3 or -2\n");
                        fatalErrorInFunction("SetWallModels", error);
                    } 
                }
                else if(mesh->boundaryT.iRight == "zeroGradient")
                {
                    ibmBody->tempBC = mesh->boundaryT.iRight;
                }
                else if(mesh->boundaryT.iRight == "fixedValue")
                {
                    ibmBody->tempBC = mesh->boundaryT.iRight;
                    ibmBody->fixedTemp = mesh->boundaryT.iRval;
                }
                else 
                {
                    char error[512];
                    sprintf(error, "Cannot match BC set in iRight to IBM boundary, set BC in IBM file\n");
                    fatalErrorInFunction("SetWallModels", error);
                }
            }
            else if(ibmBody->tBCSetType == "matchTjLeft")
            {
                if (mesh->boundaryT.jLeft == "thetaWallFunction")
                {
                    ibmBody->tempBC = mesh->boundaryT.jLeft;

                    if(ibmBody->velocityBC != "velocityWallFunction" && ibmBody->wallFunctionTypeU !=-3)
                    {
                        char error[512];
                        sprintf(error, "thetaWallFunction can be used currently only for velocity BC = velocityWallFunction and type = -3 (Schumann Model)\n");
                        fatalErrorInFunction("SetIBMWallModels", error);   
                    }

                    //set theta wall function type
                    readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(ibmBody->wallFunctionTypeT));

                    // allocate memory and connect pointers
                    PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelT));

                    if (ibmBody->wallFunctionTypeT == -3)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                        char error[512];
                        sprintf(error, "thetaWallFunction type -3 is not yet implemented for IBM\n");
                        fatalErrorInFunction("SetIBMWallModels", error); 
                    }
                    else if(ibmBody->wallFunctionTypeT == -2)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;

                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
                    }
                    else if(ibmBody->wallFunctionTypeT == -4)
                    {

                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));

                        //read the surface temp and Obhukhov length 
                        readSurfaceTempData(wm);
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "invalid wall model chosen. Please use option -3, -4 or -2\n");
                        fatalErrorInFunction("SetWallModels", error);
                    } 
                }
                else if(mesh->boundaryT.jLeft == "zeroGradient")
                {
                    ibmBody->tempBC = mesh->boundaryT.jLeft;
                }
                else if(mesh->boundaryT.jLeft == "fixedValue")
                {
                    ibmBody->tempBC = mesh->boundaryT.jLeft;
                    ibmBody->fixedTemp = mesh->boundaryT.jLval;
                }
                else 
                {
                    char error[512];
                    sprintf(error, "Cannot match BC set in jLeft to IBM boundary, set BC in IBM file\n");
                    fatalErrorInFunction("SetWallModels", error);
                }
            }
            else if(ibmBody->tBCSetType == "matchTjRight")
            {
                if (mesh->boundaryT.jRight == "thetaWallFunction")
                {
                    ibmBody->tempBC = mesh->boundaryT.jRight;

                    if(ibmBody->velocityBC != "velocityWallFunction" && ibmBody->wallFunctionTypeU !=-3)
                    {
                        char error[512];
                        sprintf(error, "thetaWallFunction can be used currently only for velocity BC = velocityWallFunction and type = -3 (Schumann Model)\n");
                        fatalErrorInFunction("SetIBMWallModels", error);   
                    }

                    //set theta wall function type
                    readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(ibmBody->wallFunctionTypeT));

                    // allocate memory and connect pointers
                    PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelT));

                    if (ibmBody->wallFunctionTypeT == -3)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                        char error[512];
                        sprintf(error, "thetaWallFunction type -3 is not yet implemented for IBM\n");
                        fatalErrorInFunction("SetIBMWallModels", error); 
                    }
                    else if(ibmBody->wallFunctionTypeT == -2)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;

                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "invalid wall model chosen. Please use option -3 or -2\n");
                        fatalErrorInFunction("SetWallModels", error);
                    } 
                }
                else if(mesh->boundaryT.jRight == "zeroGradient")
                {
                    ibmBody->tempBC = mesh->boundaryT.jRight;
                }
                else if(mesh->boundaryT.jRight == "fixedValue")
                {
                    ibmBody->tempBC = mesh->boundaryT.jRight;
                    ibmBody->fixedTemp = mesh->boundaryT.jRval;
                }
                else 
                {
                    char error[512];
                    sprintf(error, "Cannot match BC set in jRight to IBM boundary, set BC in IBM file\n");
                    fatalErrorInFunction("SetWallModels", error);
                }
            }
            else if(ibmBody->tBCSetType == "setHere")
            {
                readSubDictWord("./IBM/IBMProperties.dat", objectName, "temperatureBC", &(ibmBody->tempBC));

                if (ibmBody->tempBC == "thetaWallFunction")
                {
                    if(ibmBody->velocityBC != "velocityWallFunction" && ibmBody->wallFunctionTypeU !=-3)
                    {
                        char error[512];
                        sprintf(error, "thetaWallFunction can be used currently only for velocity BC = velocityWallFunction and type = -3 (Schumann Model)\n");
                        fatalErrorInFunction("SetIBMWallModels", error);   
                    }

                    //set theta wall function type
                    readSubDictInt("./IBM/IBMProperties.dat", objectName, "wallFunctionTypeT", &(ibmBody->wallFunctionTypeT));

                    // allocate memory and connect pointers
                    PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallModelT));

                    if (ibmBody->wallFunctionTypeT == -3)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                        char error[512];
                        sprintf(error, "thetaWallFunction type -3 is not yet implemented for IBM\n");
                        fatalErrorInFunction("SetIBMWallModels", error); 
                    }
                    else if(ibmBody->wallFunctionTypeT == -2)
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;

                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
                    }
                    else if(ibmBody->wallFunctionTypeT == -4)
                    {

                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallModelT->wmShumann));

                        Shumann *wm = ibmBody->ibmWallModelT->wmShumann;
                        readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                        //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                        readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));

                        //read the surface temp and Obhukhov length 
                        readSurfaceTempData(wm);
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "invalid wall model chosen. Please use option -3, -4 or -2\n");
                        fatalErrorInFunction("SetWallModels", error);
                    } 
                }
                else if(ibmBody->tempBC == "zeroGradient")
                {
                }
                else if(ibmBody->tempBC == "fixedValue")
                {
                    readSubDictDouble("./IBM/IBMProperties.dat", objectName, "fixedValueT", &(ibmBody->fixedTemp));
                }
                else 
                {
                    char error[512];
                    sprintf(error, "Invalid temperature BC type, use temperatureWallFunction, zeroGradient or fixedValue\n");
                    fatalErrorInFunction("SetWallModels", error);
                }
            }
        }
    }

    PetscPrintf(mesh->MESH_COMM, "done\n\n");

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}
//***************************************************************************************************************//
PetscErrorCode ComputeForceMoment(ibm_ *ibm)
{

    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    peqn_         *peqn  = ibm->access->peqn;
    constants_    *cst   = ibm->access->constants;
    clock_        *clock = ibm->access->clock;
    io_           *io    = ibm->access->io;

    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys  = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs  = info.gzs, gze = info.gzs + info.gzm;

    cellIds       initCp;

    PetscInt      i1, j1, k1;
    PetscInt      b, c, s, e;
    PetscReal     sc, sb, uDot, ustar;
    Cmpnts        ***ucat, ***cent, uDiff, ibmPtVel, ibmPtVelPrev, bPtVel;
    Cmpnts        eN, eT1, eT2, dI, dF, eC;
    Cmpnts        ***csi,  ***eta,  ***zet;
    Cmpnts        checkPt, del, checkPtInit;                                // point where the pressure or shear stress needs to be interpolated
    PetscReal     refLength, eA;
    PetscReal     dx, dy, dz, area, diag;                                  // cell size in the x, y and z direction

    PetscInt      n1, n2, n3;
    PetscReal     ***lP, ***nvert, ***aj, bPtP;

    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, peqn->lP, &lP);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    //vectors to store the net force on each ibm
    std::vector<Cmpnts>    lPForce(ibm->numBodies);
    std::vector<Cmpnts>    gPForce(ibm->numBodies);
    std::vector<Cmpnts>    lVForce(ibm->numBodies);
    std::vector<Cmpnts>    gVForce(ibm->numBodies);

    std::vector<Cmpnts>    lMoment(ibm->numBodies);
    std::vector<Cmpnts>    gMoment(ibm->numBodies);

    for (b = 0; b < ibm->numBodies; b++)
    {
        lPForce[b] = nSetZero();
        gPForce[b] = nSetZero();

        lVForce[b] = nSetZero();
        gVForce[b] = nSetZero();

        lMoment[b] = nSetZero();
        gMoment[b] = nSetZero();
    }

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody = ibm->ibmBody[b];
        ibmMesh     *ibMsh   = ibmBody->ibMsh;
        ibmRotation *ibmRot  = ibmBody->ibmRot;

        //check if processor controls this ibm body
        if(ibmBody->ibmControlled)
        {
            PetscMPIInt   ibmnprocs; MPI_Comm_size(ibmBody->IBM_COMM, &ibmnprocs);
            PetscMPIInt   ibmrank;   MPI_Comm_rank(ibmBody->IBM_COMM, &ibmrank);

            PetscReal *gElemPressure = new PetscReal[ibMsh->elems];
            PetscReal *lElemPressure = new PetscReal[ibMsh->elems];

            //loop through the ibm mesh elements
            for(e = 0; e < ibMsh->elems; e++)
            {
                //set element pressure variable to 0
                gElemPressure[e] = 0.0;
                lElemPressure[e] = 0.0;

                ibmBody->ibmPForce[e] = nSetZero();

                //check that the element is inside the domain, if not then dont calculate force
                if(!isInsideBoundingBox(ibMsh->eCent[e], mesh->bounds)) continue;

                //this processor controls this ibm element
                if(ibmBody->thisPtControlled[e])
                {

                    eN = ibMsh->eN[e], eT1 = ibMsh->eT1[e], eT2 = ibMsh->eT2[e];
                    n1 = ibMsh->nID1[e], n2 = ibMsh->nID2[e], n3 = ibMsh->nID3[e];
                    eC = ibMsh->eCent[e];
                    eA = ibMsh->eA[e];

                    //current closest ibm fluid cell to point
                    PetscInt    ci, cj, ck;

                    ci = ibmBody->closestCells[e].i;
                    cj = ibmBody->closestCells[e].j;
                    ck = ibmBody->closestCells[e].k;

                    //reference length - square root of the element area vector magnitude
                    refLength = pow( 2.0 * ibMsh->eA[e], 1./2.);

                    //find the fluid cell size in each direction
                    area = nMag(csi[ck][cj][ci]); dy = 1.0/aj[ck][cj][ci]/area;
                    area = nMag(eta[ck][cj][ci]); dz = 1.0/aj[ck][cj][ci]/area;
                    area = nMag(zet[ck][cj][ci]); dx = 1.0/aj[ck][cj][ci]/area;

                    // diagonal main diagonal length
                    diag = pow( dx*dx + dy*dy + dz*dz, 1./2.);

                    //check that the IBM element is indeed close to the closest fluid mesh cell and not some IBM element inside the body (for multi-body overlapped IBM meshes)
                    PetscReal elemDist = nMag(nSub(cent[ck][cj][ci], eC));

                    if(elemDist > 2.0*diag)
                    {
                        if(ibm->dbg) PetscPrintf(PETSC_COMM_SELF,"inside IBM element, En = %lf %lf %lf\n", eN.x, eN.y, eN.z);
                        continue;
                    }

                    //initial point for interpolation
                    checkPt = nSum(nScale(refLength, eN), eC);

                    //find the closest neighbour of ck, cj, ci to checkPt
                    PetscInt    i1, j1, k1;
                    PetscReal   dmin = 10e6, d;
                    PetscInt    ic, jc, kc;

                    for (k1=ck-1; k1<ck+2; k1++)
                    for (j1=cj-1; j1<cj+2; j1++)
                    for (i1=ci-1; i1<ci+2; i1++)
                    {
                        if
                        (
                            (
                                k1>=1 && k1<mz-1 &&
                                j1>=1 && j1<my-1 &&
                                i1>=1 && i1<mx-1
                            ) && (!isIBMSolidCell(k1, j1, i1, nvert))
                        )
                        {
                            d = pow((checkPt.x - cent[k1][j1][i1].x), 2) +
                                pow((checkPt.y - cent[k1][j1][i1].y), 2) +
                                pow((checkPt.z - cent[k1][j1][i1].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                ic = i1;
                                jc = j1;
                                kc = k1;
                            }
                        }
                    }

                    PetscInt intId[6];
                    PetscInt intFlag = 1;

                    // get the trilinear interpolation cells
                    PointInterpolationCells
                    (
                            mesh,
                            checkPt.x, checkPt.y, checkPt.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );

                    // save the initial closest cell
                    initCp.i = ic; initCp.j = jc; initCp.k = kc;
                    checkPtInit = nSet(checkPt);

                    // max distance checkpoint has to move to be to the closest fluid trilinear interpolation box
                    // restrict this to 3*ref length
                    PetscReal sumDel = 0.0;

                    // if there is an ibm solid cell in the interpolation, increase reference length further
                    while ( (intFlag == 1) && isInsideBoundingBox(checkPt, mesh->bounds) && (sumDel <= 1.5*refLength))
                    {
                        PetscInt ibmCellCtr = 0;
                        PetscInt icc, jcc, kcc, setFlag = 0;

                        for (PetscInt kk = 0; kk<2; kk++)
                        for (PetscInt jj = 0; jj<2; jj++)
                        for (PetscInt ii = 0; ii<2; ii++)
                        {
                            if(isIBMSolidCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                            {
                                ibmCellCtr ++;
                            }

                        }

                        if (ibmCellCtr > 0)
                        {
                            del =  nScale(0.2 * refLength, eN);
                            mSum(checkPt, del);
                            dmin = 10e6;

                            sumDel += nMag(del);

                            for (k1=kc-1; k1<kc+2; k1++)
                            {
                                //check processor ghost bounds
                                if (k1 < gzs || k1 >= gze) {intFlag = 2; break;}
                                for (j1=jc-1; j1<jc+2; j1++)
                                {
                                    if (j1 < gys || j1 >= gye) {intFlag = 2; break;}
                                    for (i1=ic-1; i1<ic+2; i1++)
                                    {
                                        if (i1 < gxs || i1 >= gxe) {intFlag = 2; break;}

                                        if
                                        (
                                            (
                                                k1>=1 && k1<mz-1 &&
                                                j1>=1 && j1<my-1 &&
                                                i1>=1 && i1<mx-1
                                            ) && (!isIBMSolidCell(k1, j1, i1, nvert))
                                        )
                                        {
                                            d = pow((checkPt.x - cent[k1][j1][i1].x), 2) +
                                                pow((checkPt.y - cent[k1][j1][i1].y), 2) +
                                                pow((checkPt.z - cent[k1][j1][i1].z), 2);

                                            if
                                            (
                                                d < dmin
                                            )
                                            {
                                                dmin  = d;
                                                icc = i1;
                                                jcc = j1;
                                                kcc = k1;
                                                setFlag = 1;
                                            }
                                        }

                                    }
                                }
                            }

                            if(setFlag == 0)
                            {
                                //closest point not set, do not interpolate
                                intFlag = 2;
                            }

                            if( intFlag == 1)
                            {
                                kc = kcc; jc = jcc; ic = icc;

                                PointInterpolationCells
                                (
                                        mesh,
                                        checkPt.x, checkPt.y, checkPt.z,
                                        ic, jc, kc,
                                        cent,
                                        intId
                                );
                            }
                        }
                        else
                        {
                            intFlag = 0;

                            scalarPointLocalVolumeInterpolation
                            (
                                    mesh,
                                    checkPt.x, checkPt.y, checkPt.z,
                                    ic, jc, kc,
                                    cent,
                                    lP,
                                    lElemPressure[e]
                            );

                            vectorPointLocalVolumeInterpolation
                            (
                                    mesh,
                                    checkPt.x, checkPt.y, checkPt.z,
                                    ic, jc, kc,
                                    cent,
                                    ucat,
                                    bPtVel
                            );

                        }
                    }

                    // while loop fails
                    if(intFlag > 0)
                    {
                        lElemPressure[e] = lP[initCp.k][initCp.j][initCp.i];

                        bPtVel = nSet(ucat[initCp.k][initCp.j][initCp.i]);

                        checkPt = nSet(checkPtInit);
                    }

                    //Velocity of the ibm mesh element
                    ibmPtVel.x =  (ibMsh->nU[n1].x
                                 + ibMsh->nU[n2].x
                                 + ibMsh->nU[n3].x) / 3.0 ;

                    ibmPtVel.y =  (ibMsh->nU[n1].y
                                 + ibMsh->nU[n2].y
                                 + ibMsh->nU[n3].y) / 3.0 ;

                    ibmPtVel.z =  (ibMsh->nU[n1].z
                                 + ibMsh->nU[n2].z
                                 + ibMsh->nU[n3].z) / 3.0;

                    // interpolate the previous velocity
                    ibmPtVelPrev.x =   (ibMsh->nUPrev[n1].x
                                      + ibMsh->nUPrev[n2].x
                                      + ibMsh->nUPrev[n3].x) / 3.0 ;

                    ibmPtVelPrev.y =   (ibMsh->nUPrev[n1].y
                                      + ibMsh->nUPrev[n2].y
                                      + ibMsh->nUPrev[n3].y) / 3.0 ;

                    ibmPtVelPrev.z =   (ibMsh->nUPrev[n1].z
                                      + ibMsh->nUPrev[n2].z
                                      + ibMsh->nUPrev[n3].z) / 3.0;

                    // save the element acceleration term (dudt . elementNormal)
                    PetscReal bP = nDot(nScale(-1.0/clock->dtOld, nSub(ibmPtVel, ibmPtVelPrev)), eN);

                    // distance between interpolation point and ibm mesh element
                    sc    = nMag(nSub(checkPt, eC));

                    //set the ibm element pressure force
                    PetscReal pForce = (lElemPressure[e] - bP * sc) * eA * cst->rho;
                    ibmBody->ibmPForce[e] = nSet(nScale(-pForce, eN));

                    //find the friction velocity ustar
                    uDiff = nSub(bPtVel, ibmPtVel);

                    PetscReal  un, ut;              //normal and tangential component of udiff
                    Cmpnts     unV, utV;            // normal and tanential vector of udiff

                    // normal component to the wall
                    un = nDot(uDiff, eN);
                    unV = nScale(un,eN);

                    // tangential component to the wall
                    utV = nSub(uDiff, unV);
                    ut = nMag(utV);

                    Cmpnts tauWall = nSetZero();

                    if(ibmBody->velocityBC != "slip")
                    {
                        //friction velocity
                        ustar = uTauCabot(cst->nu, ut, sc, 0.01, 0);

                        if(ibmBody->velocityBC == "noSlip")
                        {
                            tauWall = nScale(eA * cst->rho * cst->nu * ut/sc, nUnit(utV));
                        }
                        else
                        {
                            //shear force on ibm mesh element
                            tauWall = nScale(eA * cst->rho * ustar * ustar, nUnit(utV));
                        }

                    }

                    //sum total pressure and viscous force
                    mSum(lPForce[b], ibmBody->ibmPForce[e]);

                    if(ibmBody->velocityBC != "Slip")
                    {
                        mSum(lVForce[b], tauWall);
                    }

                    if(ibmBody->bodyMotion == "rotation")
                    {
                        // find moment for rotating mesh

                        //radius vector from the rotation center to the ibm mesh coordinate
                        Cmpnts rvec   = nSub(eC, ibmRot->rotCenter);

                        Cmpnts pressureMoment = nCross(rvec, ibmBody->ibmPForce[e]), viscousMoment;

                        if(ibmBody->velocityBC != "slip")
                        {
                            viscousMoment  = nCross(rvec, tauWall);
                        }
                        else
                        {
                            viscousMoment  = nSetZero();
                        }

                        //net moment per processor
                        mSum(lMoment[b], nSum(pressureMoment, viscousMoment));

                    }


                }

            }

            MPI_Allreduce(&(lPForce[b]), &(gPForce[b]), 3, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);
            MPI_Allreduce(&(lVForce[b]), &(gVForce[b]), 3, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);
            MPI_Allreduce(lElemPressure, gElemPressure, ibMsh->elems, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);

            if(ibmBody->bodyMotion == "rotation")
            {
                MPI_Allreduce(&(lMoment[b]), &(gMoment[b]), 3, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);
            }

            //compute power
            PetscReal ibmPower = 0.0;
            PetscReal netMoment = 0.0;

            if(ibmBody->bodyMotion == "rotation")
            {
                ibmPower  = nDot(gMoment[b], nScale(ibmRot->angSpeed, ibmRot->rotAxis));
                netMoment = nDot(gMoment[b], ibmRot->rotAxis);
            }

            PetscPrintf(ibmBody->IBM_COMM, "%s: pforce(N) = %lf %lf %lf\n", ibmBody->bodyName.c_str(), gPForce[b].x, gPForce[b].y, gPForce[b].z);
            PetscPrintf(ibmBody->IBM_COMM, "%s: Moment(Nm) = %lf %lf %lf\n", ibmBody->bodyName.c_str(), gMoment[b].x, gMoment[b].y, gMoment[b].z);

            //write data
            if(ibmrank == 0)
            {
                writeIBMForceData(ibm, b, gElemPressure, netMoment, ibmPower, gPForce[b], gVForce[b], gMoment[b]);
            }

            delete[] gElemPressure;
            delete[] lElemPressure;
        }

    }

    std::vector<Cmpnts> ().swap(lPForce);
    std::vector<Cmpnts> ().swap(gPForce);
    std::vector<Cmpnts> ().swap(lVForce);
    std::vector<Cmpnts> ().swap(gVForce);

    std::vector<Cmpnts> ().swap(lMoment);
    std::vector<Cmpnts> ().swap(gMoment);

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(da, peqn->lP, &lP);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    MPI_Barrier(mesh->MESH_COMM);

    return 0;

}

//***************************************************************************************************************//
PetscErrorCode  writeIBMForceData(ibm_ *ibm, PetscInt b, PetscReal *gElemPressure, PetscReal netMoment, PetscReal ibmPower, Cmpnts PresForce, Cmpnts ViscForce, Cmpnts momentVector)
{
    clock_      *clock = ibm->access->clock;
    io_         *io    = ibm->access->io;
    mesh_       *mesh  = ibm->access->mesh;
    ibmObject   *ibmBody = ibm->ibmBody[b];
    ibmMesh     *ibMsh = ibmBody->ibMsh;

    // write flags for current time step
    PetscInt    writeNow         =  0;

    PetscReal   timeStart     = ibm->timeStart;
    PetscReal   timeInterval  = ibm->timeInterval;
    word        intervalType     = ibm->intervalType;

    PetscReal   epsilon          = 1e-6;

    word        timeName;
    word        elementForcePath, netForcePath;

    // set time folder name
    timeName = getTimeName(clock);
    elementForcePath     = "./postProcessing/" + mesh->meshName + "/IBM/" + getStartTimeName(clock) + "/" + ibmBody->bodyName.c_str() + "/elementForce/" + timeName;
    netForcePath         = "./postProcessing/" + mesh->meshName + "/IBM/" + getStartTimeName(clock) + "/" + ibmBody->bodyName.c_str() + "/netForce/";

    // create/initialize ibm write directory
    if
    (
        clock->it == clock->itStart
    )
    {
        errno = 0;
        word ibmFolder = "./postProcessing/" + mesh->meshName + "/IBM/";

        PetscInt dirRes = mkdir(ibmFolder.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create %s directory", ibmFolder.c_str());
            fatalErrorInFunction("writeElementPForce",  error);
        }

        word timeFolder = ibmFolder + getStartTimeName(clock);
        dirRes = mkdir(timeFolder.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create %s directory", timeFolder.c_str());
            fatalErrorInFunction("writeElementPForce",  error);
        }

        word ibmBodyFolder = timeFolder + "/" + ibmBody->bodyName.c_str();
        dirRes = mkdir(ibmBodyFolder.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create %s directory", ibmBodyFolder.c_str());
            fatalErrorInFunction("writeElementPForce",  error);
        }

        if(io->writePForce)
        {
            word elementForceFolder = ibmBodyFolder + "/elementForce/";
            dirRes = mkdir(elementForceFolder.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", elementForceFolder.c_str());
                fatalErrorInFunction("writeElementPForce",  error);
            }
        }

        word netForceFolder = ibmBodyFolder + "/netForce/";
        dirRes = mkdir(netForceFolder.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create %s directory", netForceFolder.c_str());
            fatalErrorInFunction("writeElementPForce",  error);
        }

        //create net force file
        word  fileName = netForcePath + "/" + ibmBody->bodyName.c_str() + "_netForce";

        FILE *fp;
        fp = fopen(fileName.c_str(), "w");

        if(!fp)
        {
           char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("ComputeForceMoment",  error);
        }
        else
        {
            word w1   = "time [s]";
            word w2   = "aeroTq [Nm]";
            word w3   = "aeroPwr [W]";
            word w4   = "omega [rpm]";
            word w5   = "FPx [N]";
            word w6   = "FPy [N]";
            word w7   = "FPz [N]";
            word w8   = "FVx [N]";
            word w9   = "FVy [N]";
            word w10   = "FVz [N]";
            word w11   = "Mx [Nm]";
            word w12   = "My [Nm]";
            word w13   = "Mz [Nm]";
            int width = -20;

            if(fp)
            {
                if(ibmBody->bodyMotion == "rotation")
                {
                    fprintf(fp, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n", width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(),  width, w5.c_str(),  width, w6.c_str(),  width, w7.c_str(), width, w8.c_str(), width, w9.c_str(), width, w10.c_str(), width, w11.c_str(), width, w12.c_str(), width, w13.c_str());
                }
                else if((ibmBody->bodyMotion == "sinusoidal") || (ibmBody->bodyMotion == "pitchingOscillation"))
                {
                    fprintf(fp, "%*s %*s %*s %*s %*s %*s %*s\n", width, w1.c_str(), width, w5.c_str(),  width, w6.c_str(),  width, w7.c_str(), width, w8.c_str(), width, w9.c_str(), width, w10.c_str());
                }
                else if (ibmBody->bodyMotion == "static")
                {
                    fprintf(fp, "%*s %*s %*s %*s %*s %*s %*s\n", width, w1.c_str(), width, w5.c_str(),  width, w6.c_str(),  width, w7.c_str(), width, w8.c_str(),  width, w9.c_str(),  width, w10.c_str());
                }
            }

            fclose(fp);
        }

    }

    // check if must accumulate
    if
    (
        (intervalType == "adjustableTime") &&
        (clock->time >= timeStart) &&
        (
            fabs((clock->time - timeStart) / timeInterval -
            std::floor
            (
                (clock->time - timeStart) / timeInterval + epsilon
            )) < 1e-6
        )
    )
    {
        writeNow = 1;
    }
    // write every "timeInterval" iterations
    else if
    (
        (clock->it > 0) &&
        (intervalType == "timeStep") &&
        (
            fabs(clock->it / timeInterval -
            std::floor
            (
                clock->it / timeInterval + epsilon
            )) < 1e-10
        )
    )
    {
        writeNow = 1;
    }

    if(writeNow)
    {
        word    fileName;
        FILE    *f;

        if(io->writePForce)
        {
            //write element force
            fileName = elementForcePath + "/" + ibmBody->bodyName.c_str() + "_pForce.dlo";
            PetscInt dirRes = mkdir(elementForcePath.c_str(), 0777);

            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory\n", elementForcePath.c_str());
                fatalErrorInFunction("writeElementPForce",  error);
            }

            f = fopen(fileName.c_str(), "w");

            if(!f)
            {
               char error[512];
                sprintf(error, "cannot open file %s\n", fileName.c_str());
                fatalErrorInFunction("writeElementPForce",  error);
            }
            else
            {
                fprintf(f, "** Pressure based on Esurface\n");

                if(ibmBody->bodyType == "surfaceBody")
                {

                    for(PetscInt e = 0; e < ibMsh->elems; e++)
                    {
                        for(PetscInt s = 0; s < ibmBody->numSurfaces; s++)
                        {
                            surface* ibmSurface = ibmBody->ibmSurface[s];

                            if(ibmSurface->surfaceFileType == "inp")
                            {
                                if(ibmSurface->surfaceId == ibMsh->eSurface[e])
                                {
                                    fprintf(f, "%ld, P, %lf \n", ibmBody->elementMapping[e], gElemPressure[e]);
                                }
                            }
                        }
                    }

                }
                else
                {
                    // use global element id for "inp" file format
                    if(ibmBody->fileType == "inp")
                    {
                        for(PetscInt e = 0; e < ibMsh->elems; e++)
                        {
                            fprintf(f, "%ld, P, %lf \n", ibmBody->elementMapping[e], gElemPressure[e]);
                        }
                    }
                    else
                    {
                        for(PetscInt e = 0; e < ibMsh->elems; e++)
                        {
                            fprintf(f, "%ld, P, %lf \n", e, gElemPressure[e]);
                        }
                    }
                }

            }

            fclose(f);

        }

        //write net force
        fileName = netForcePath + "/" + ibmBody->bodyName.c_str() + "_netForce";

        f = fopen(fileName.c_str(), "a");

        if(f)
        {
            int width = -20;

            if(ibmBody->bodyMotion == "rotation")
            {
                fprintf(f, "%*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f\n", width, clock->time, width, netMoment, width, ibmPower, width, ibmBody->ibmRot->angSpeed * 60/(2*M_PI), width, PresForce.x,  width, PresForce.y,  width, PresForce.z, width, ViscForce.x,  width, ViscForce.y,  width, ViscForce.z, width, momentVector.x, width, momentVector.y, width, momentVector.z);
            }
            else if((ibmBody->bodyMotion == "sinusoidal") || (ibmBody->bodyMotion == "pitchingOscillation"))
            {
                fprintf(f, "%*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f\n", width, clock->time, width, PresForce.x,  width, PresForce.y,  width, PresForce.z, width, ViscForce.x,  width, ViscForce.y,  width, ViscForce.z);
            }
            else if (ibmBody->bodyMotion == "static")
            {
                fprintf(f, "%*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f\n", width, clock->time, width, PresForce.x,  width, PresForce.y,  width, PresForce.z, width, ViscForce.x,  width, ViscForce.y,  width, ViscForce.z);
            }

        }

        fclose(f);
    }

    return (0);
}

//***************************************************************************************************************//
PetscErrorCode UpdateIBMesh(ibm_ *ibm)
{

    // loop through the ibm bodies
    for(PetscInt b = 0; b < ibm->numBodies; b++)
    {
        // update mesh only if not static body
        if(ibm->ibmBody[b]->bodyMotion != "static")
        {
            if(ibm->ibmBody[b]->bodyMotion == "rotation")
            {
                rotateIBMesh(ibm, b);
            }
            else if(ibm->ibmBody[b]->bodyMotion == "sinusoidal")
            {
                sineMotion(ibm, b);
            }
            else if(ibm->ibmBody[b]->bodyMotion == "pitchingOscillation")
            {
                pitchingMotion(ibm, b);
            }

            //write the new position into STL file
            if(ibm->writeSTL)
            {
                {
                    writeSTLFile(ibm, b);
                    MPI_Barrier(ibm->access->mesh->MESH_COMM);
                }
            }
        }

    }

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode recomputeIBMeshProperties(ibm_ *ibm, PetscInt b)
{
    ibmObject *ibmBody = ibm->ibmBody[b];
    ibmMesh   *ibMsh   = ibmBody->ibMsh;
    PetscInt  n1, n2, n3;
    Cmpnts    vec1, vec2, temp;

    for (PetscInt i=0; i<ibMsh->elems; i++)
    {
        PetscReal normMag;

        // get the element nodes
        n1 = ibMsh->nID1[i]; n2 = ibMsh->nID2[i]; n3 = ibMsh->nID3[i];

        vec1 = nSub(ibMsh->nCoor[n2], ibMsh->nCoor[n1]);

        vec2 = nSub(ibMsh->nCoor[n3], ibMsh->nCoor[n1]);

        // normal to the face is found as cross product of the edges vec1 and vec2
        ibMsh->eN[i] = nCross(vec1, vec2);
        normMag = nMag(ibMsh->eN[i]);
        mScale(1.0/normMag, ibMsh->eN[i]);

        //element area
        ibMsh->eA[i] = normMag/2.0;

        //element center
        temp = nSum(ibMsh->nCoor[n1], ibMsh->nCoor[n2]);
        ibMsh->eCent[i] = nSum( temp, ibMsh->nCoor[n3]);
        mScale(1/3.0, ibMsh->eCent[i]);
    }

    //use nearby elements to average the normal direction. this requires the node to element connectivity
    if(ibm->averageNormal)
    {
        ibmNode *ibmNodes  = ibMsh->ibmMeshNode;

        //temporary array to store the averaged element normals
        Cmpnts *tempEn = new Cmpnts[ibMsh->elems];

        for (PetscInt e=0; e<ibMsh->elems; e++)
        {
            PetscInt i, j;
            // get the element nodes
            n1 = ibMsh->nID1[e]; n2 = ibMsh->nID2[e]; n3 = ibMsh->nID3[e];

            //find the number of elements connected to each node
            PetscInt numElemn1 = ibmNodes[n1].numConnected;
            PetscInt numElemn2 = ibmNodes[n2].numConnected;
            PetscInt numElemn3 = ibmNodes[n3].numConnected;

            PetscInt mergedElemArray[numElemn1 + numElemn2 + numElemn3];

            //copy the first node elements into mergedElemArray
            for (i = 0; i < numElemn1; i++)
            {
                mergedElemArray[i] = ibmNodes[n1].elem[i];
            }

            //copy the second node elements into merged array
            for (j = 0; j < numElemn2; j++)
            {
                if (!isPresent(mergedElemArray, i+1, ibmNodes[n2].elem[j]))
                {
                    mergedElemArray[i++] = ibmNodes[n2].elem[j];
                }
            }

            //copy the third node elements into merged array
            for (j = 0; j < numElemn3; j++)
            {
                if (!isPresent(mergedElemArray, i+1, ibmNodes[n3].elem[j]))
                {
                    mergedElemArray[i++] = ibmNodes[n3].elem[j];
                }
            }

            //average the normal at element e from its surrounding elements

            PetscReal sumEnx = 0.0, sumEny = 0.0, sumEnz = 0.0;

            for (j = 0; j < i; j++)
            {
                sumEnx += ibMsh->eN[mergedElemArray[j]].x;
                sumEny += ibMsh->eN[mergedElemArray[j]].y;
                sumEnz += ibMsh->eN[mergedElemArray[j]].z;

            }

            // save the averaged element normals
            tempEn[e] = nSetFromComponents(sumEnx, sumEny, sumEnz);
            PetscReal normMag = nMag(tempEn[e]);
            mScale(1.0/normMag, tempEn[e]);
        }

        for (PetscInt e=0; e<ibMsh->elems; e++)
        {
            ibMsh->eN[e] = nSet(tempEn[e]);
        }

        delete[] tempEn;
    }

    for (PetscInt i=0; i<ibMsh->elems; i++)
    {
        // tangential to the face( eT1 and eT2)
        // eT1 = eN x k
        if (
            (((1.0 - ibMsh->eN[i].z ) <= 1e-6 ) && ((-1.0 + ibMsh->eN[i].z ) < 1e-6))
            ||
            (((ibMsh->eN[i].z + 1.0 ) <= 1e-6 ) && ((-1.0 - ibMsh->eN[i].z ) < 1e-6))
           )
        {
            ibMsh->eT1[i].x = 1.0;
            ibMsh->eT1[i].y = 0.0;
            ibMsh->eT1[i].z = 0.0;

            ibMsh->eT2[i].x = 0.0;
            ibMsh->eT2[i].y = 1.0;
            ibMsh->eT2[i].z = 0.0;
        }
        else
        {
            ibMsh->eT1[i].x =  ibMsh->eN[i].y/ sqrt(ibMsh->eN[i].x*ibMsh->eN[i].x + ibMsh->eN[i].y*ibMsh->eN[i].y);
            ibMsh->eT1[i].y = -ibMsh->eN[i].x/ sqrt(ibMsh->eN[i].x*ibMsh->eN[i].x + ibMsh->eN[i].y*ibMsh->eN[i].y);
            ibMsh->eT1[i].z = 0 ;

               // eT2 = eT2 x eN
            ibMsh->eT2[i].x = -ibMsh->eN[i].x*ibMsh->eN[i].z/ sqrt(ibMsh->eN[i].x*ibMsh->eN[i].x + ibMsh->eN[i].y*ibMsh->eN[i].y);
            ibMsh->eT2[i].y = -ibMsh->eN[i].y*ibMsh->eN[i].z/ sqrt(ibMsh->eN[i].x*ibMsh->eN[i].x + ibMsh->eN[i].y*ibMsh->eN[i].y);
            ibMsh->eT2[i].z = sqrt(ibMsh->eN[i].x*ibMsh->eN[i].x + ibMsh->eN[i].y*ibMsh->eN[i].y);
        }
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode pitchingMotion(ibm_ *ibm, PetscInt b)
{
    clock_         *clock   = ibm->access->clock;

    ibmObject      *ibmBody  = ibm->ibmBody[b];
    ibmMesh        *ibMsh    = ibmBody->ibMsh;                         // pointer to the ibm body mesh
    ibmPitchMotion *ibmPitch = ibmBody->ibmPitch;

    PetscReal      tCurrent = clock->time;

    PetscReal      delAlpha = ibmPitch->amplitude * ( cos(2*M_PI * ibmPitch->frequency * ibmPitch->tPrev) -  cos(2*M_PI * ibmPitch->frequency * tCurrent) );

    Cmpnts         rvec;
    Cmpnts         angVel = nScale(2*M_PI * ibmPitch->frequency * ibmPitch->amplitude * sin(2*M_PI * ibmPitch->frequency * ibmPitch->tPrev), ibmPitch->pitchAxis);

    // find the new ibm node co-ordinate after motion
    for(PetscInt n = 0; n < ibMsh->nodes; n++)
    {
        rvec    = nSub(ibMsh->nCoor[n], ibmPitch->pitchCenter);

        // rotate the vector to the required angular position
        mRot(ibmPitch->pitchAxis, rvec, delAlpha);

        // find the new co-ordinate
        ibMsh->nCoor[n]  = nSum(rvec, ibmPitch->pitchCenter);

        //save the old velocity
        ibMsh->nUPrev[n] = nSet(ibMsh->nU[n]);

        //find the new velocity
        ibMsh->nU[n]    = nCross(angVel, rvec);
    }

    ibmPitch->tPrev  = tCurrent;

    recomputeIBMeshProperties(ibm, b);

    return (0);
}

//***************************************************************************************************************//
PetscErrorCode sineMotion(ibm_ *ibm, PetscInt b)
{
    clock_        *clock   = ibm->access->clock;

    ibmObject     *ibmBody = ibm->ibmBody[b];
    ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
    ibmSineMotion *ibmSine  = ibmBody->ibmSine;

    PetscReal     tCurrent = clock->time;

    // find the new ibm node co-ordinate after motion
    for(PetscInt n = 0; n < ibMsh->nodes; n++)
    {
        PetscReal ct = ibmSine->amplitude * ( cos(2*M_PI * ibmSine->frequency * ibmSine->tPrev) -  cos(2*M_PI * ibmSine->frequency * tCurrent) );

        // find the new co-ordinate
        ibMsh->nCoor[n] = nSum(ibMsh->nCoor[n], nScale(ct, ibmSine->motionDir));

        //save the old velocity
        ibMsh->nUPrev[n] = nSet(ibMsh->nU[n]);

        //find the new velocity
        ibMsh->nU[n]    = nScale(2*M_PI*ibmSine->frequency*ibmSine->amplitude * sin(2*M_PI * ibmSine->frequency * tCurrent), ibmSine->motionDir);
    }

    ibmSine->tPrev  = tCurrent;

    recomputeIBMeshProperties(ibm, b);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode rotateIBMesh(ibm_ *ibm, PetscInt b)
{
    clock_        *clock   = ibm->access->clock;

    ibmObject     *ibmBody = ibm->ibmBody[b];
    ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
    ibmRotation   *ibmRot  = ibmBody->ibmRot;

    Cmpnts        rvec, rvecPar, angVel;

    // change in rotation angle
    PetscReal     dTheta   = ibmRot->angSpeed * clock->dt;

    // compute the new angular position of the body
    ibmRot->rotAngle += dTheta;

    // find the new ibm node co-ordinate after rotation
    for(PetscInt n = 0; n < ibMsh->nodes; n++)
    {
        // find the radius vector from the rotation center point
        rvec = nSub(ibMsh->nCoor[n], ibmRot->rotCenter);

        // rotate the vector
        mRot(ibmRot->rotAxis, rvec, dTheta);

        // find the new co-ordinate
        ibMsh->nCoor[n] = nSum(rvec, ibmRot->rotCenter);

        // set the node velocity
        angVel = nScale(ibmRot->angSpeed, ibmRot->rotAxis);

        //save the old velocity
        ibMsh->nUPrev[n] = nSet(ibMsh->nU[n]);

        //find the new velocity
        ibMsh->nU[n] = nCross(angVel, rvec);
    }

    //update the angular speed if body is accelerating
    ibmRot->angSpeed += (ibmRot->angAcc * clock->dt);

    recomputeIBMeshProperties(ibm, b);

    // write the current angular position of the ibm to a file
    writeIBMData(ibm, b);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode writeIBMData(ibm_ *ibm, PetscInt b)
{
    clock_      *clock = ibm->access->clock;
    io_         *io    = ibm->access->io;
    mesh_       *mesh  = ibm->access->mesh;

    word        timeName;
    word        path;

    PetscInt    t;

    PetscMPIInt rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set time folder name
    timeName = getTimeName(clock);
    path     = "./fields/" + mesh->meshName + "/ibm/" + timeName;

    // create/initialize ibm write directory
    if
    (
        clock->it == clock->itStart
    )
    {
        if(!rank)
        {

            errno = 0;
            word ibmfolder = "./fields/" + mesh->meshName + "/ibm/";
            PetscInt dirRes = mkdir(ibmfolder.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory", ibmfolder.c_str());
                fatalErrorInFunction("writeAngularPosition",  error);
            }

            // if directory already exist remove everything inside except the start time (safe)
            if(errno == EEXIST)
            {
                remove_subdirs_except(mesh->MESH_COMM, ibmfolder.c_str(), getStartTimeName(clock));
            }
        }
    }

    // test if must write at this time step
    if(io->runTimeWrite)
    {
        // creates time folder
        if(!rank)
        {
            PetscInt dirRes = mkdir(path.c_str(), 0777);

            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory\n", path.c_str());
                fatalErrorInFunction("writeAngularPosition",  error);
            }

            FILE *f;
            char fileName[80];
            sprintf(fileName, "%s/%s", path.c_str(), ibm->ibmBody[b]->bodyName.c_str());
            f = fopen(fileName, "w");

            ibmRotation   *ibmRot  = ibm->ibmBody[b]->ibmRot;

            if(!f)
            {
               char error[512];
                sprintf(error, "cannot open file %s\n", fileName);
                fatalErrorInFunction("writeAngularPosition",  error);
            }
            else
            {
                fprintf(f, "AngularPosition             %lf\n", ibmRot->rotAngle);
                fprintf(f, "AngularSpeed                %lf\n", ibmRot->angSpeed * 60/(2.0 * M_PI));
                fprintf(f, "AngularAcceleration         %lf\n", ibmRot->angAcc * 60/(2.0 * M_PI));
            }

            fclose(f);

        }

    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findInterceptionPoint(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    DM             da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo  info  = mesh->info;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscInt       i, j, k;
    PetscInt       c, nif;

    PetscInt	   ip[9], jp[9], kp[9];                        //indices for the background plane nine cells
    Cmpnts	       pc[9];

    Cmpnts      ***cent;

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmMesh   *ibMsh = ibm->ibmBody[ibF[c].bodyID]->ibMsh;
        PetscInt   cElem = ibF[c].closestElem;
        Cmpnts	   eNorm = ibMsh->eN[cElem];
        Cmpnts     pCoor = nSet(cent[k][j][i]);

        //check if the interception point is on the i-1 plane
        ip[0] = i-1; ip[1] = i-1; ip[2] = i-1;
        ip[3] = i-1; ip[4] = i-1; ip[5] = i-1;
        ip[6] = i-1; ip[7] = i-1; ip[8] = i-1;

        jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
        jp[1] = j;   jp[4] = j;   jp[7] = j;
        jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

        kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
        kp[3] = k;   kp[4] = k;   kp[5] = k;
        kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

        for (nif=0; nif<9; nif++)
        {
            pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
            pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
            pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
                ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
                ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
                ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
                break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }

        if (ibF[c].imode >= 0)
        {
            continue;
        }

        //check if the interception point is on the i+1 plane
        ip[0] = i+1; ip[1] = i+1; ip[2] = i+1;
        ip[3] = i+1; ip[4] = i+1; ip[5] = i+1;
        ip[6] = i+1; ip[7] = i+1; ip[8] = i+1;

        jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
        jp[1] = j;   jp[4] = j;   jp[7] = j;
        jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

        kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
        kp[3] = k;   kp[4] = k;   kp[5] = k;
        kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

        for (nif=0; nif<9; nif++)
        {
          pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
          pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
          pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
              ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
              ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }

        if (ibF[c].imode >=0)
        {
            continue;
        }

        //check if the interception point is on the j-1 plane
        ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
        ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
        ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

        jp[0] = j-1; jp[3] = j-1; jp[6] = j-1;
        jp[1] = j-1; jp[4] = j-1; jp[7] = j-1;
        jp[2] = j-1; jp[5] = j-1; jp[8] = j-1;

        kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
        kp[3] = k;   kp[4] = k;   kp[5] = k;
        kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

        for (nif=0; nif<9; nif++)
        {
          pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
          pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
          pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
              ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
              ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }

        if (ibF[c].imode >=0)
        {
            continue;
        }

        //check if the interception point is on the j+1 plane
        ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
        ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
        ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

        jp[0] = j+1; jp[3] = j+1; jp[6] = j+1;
        jp[1] = j+1; jp[4] = j+1; jp[7] = j+1;
        jp[2] = j+1; jp[5] = j+1; jp[8] = j+1;

        kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
        kp[3] = k;   kp[4] = k;   kp[5] = k;
        kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

        for (nif=0; nif<9; nif++)
        {
          pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
          pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
          pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
              ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
              ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }

        if (ibF[c].imode >=0)
        {
            continue;
        }

        //check if the interception point is on the k-1 plane
        ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
        ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
        ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

        jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
        jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
        jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;

        kp[0] = k-1; kp[1] = k-1; kp[2] = k-1;
        kp[3] = k-1; kp[4] = k-1; kp[5] = k-1;
        kp[6] = k-1; kp[7] = k-1; kp[8] = k-1;

        for (nif=0; nif<9; nif++)
        {
          pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
          pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
          pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
              ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
              ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7):
            {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }
        if (ibF[c].imode >=0)
        {
            continue;
        }

        //check if the interception point is on the k+1 plane
        ip[0] = i-1; ip[1] = i  ; ip[2] = i+1;
        ip[3] = i-1; ip[4] = i  ; ip[5] = i+1;
        ip[6] = i-1; ip[7] = i  ; ip[8] = i+1;

        jp[0] = j-1; jp[3] = j  ; jp[6] = j+1;
        jp[1] = j-1; jp[4] = j  ; jp[7] = j+1;
        jp[2] = j-1; jp[5] = j  ; jp[8] = j+1;

        kp[0] = k+1; kp[1] = k+1; kp[2] = k+1;
        kp[3] = k+1; kp[4] = k+1; kp[5] = k+1;
        kp[6] = k+1; kp[7] = k+1; kp[8] = k+1;

        for (nif=0; nif<9; nif++) {
          pc[nif].x = cent[kp[nif]][jp[nif]][ip[nif]].x;
          pc[nif].y = cent[kp[nif]][jp[nif]][ip[nif]].y;
          pc[nif].z = cent[kp[nif]][jp[nif]][ip[nif]].z;
        }

        interceptionPt(pCoor, pc, eNorm, &(ibF[c]));

        switch (ibF[c].imode)
        {
            case(0):
            {
              ibF[c].i1=ip[0]; ibF[c].j1 = jp[0]; ibF[c].k1 = kp[0];
              ibF[c].i2=ip[1]; ibF[c].j2 = jp[1]; ibF[c].k2 = kp[1];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (1):
            {
              ibF[c].i1=ip[1]; ibF[c].j1 = jp[1]; ibF[c].k1 = kp[1];
              ibF[c].i2=ip[2]; ibF[c].j2 = jp[2]; ibF[c].k2 = kp[2];
              ibF[c].i3=ip[4]; ibF[c].j3 = jp[4]; ibF[c].k3 = kp[4];
              break;
            }

            case (2):
            {
              ibF[c].i1=ip[2]; ibF[c].j1 = jp[2]; ibF[c].k1 = kp[2];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[5]; ibF[c].j3 = jp[5]; ibF[c].k3 = kp[5];
              break;
            }

            case (3):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[5]; ibF[c].j2 = jp[5]; ibF[c].k2 = kp[5];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (4):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[7]; ibF[c].j2 = jp[7]; ibF[c].k2 = kp[7];
              ibF[c].i3=ip[8]; ibF[c].j3 = jp[8]; ibF[c].k3 = kp[8];
              break;
            }

            case (5):
            {
              ibF[c].i1=ip[4]; ibF[c].j1 = jp[4]; ibF[c].k1 = kp[4];
              ibF[c].i2=ip[6]; ibF[c].j2 = jp[6]; ibF[c].k2 = kp[6];
              ibF[c].i3=ip[7]; ibF[c].j3 = jp[7]; ibF[c].k3 = kp[7];
              break;
            }

            case (6): {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[6]; ibF[c].j3 = jp[6]; ibF[c].k3 = kp[6];
              break;
            }

            case (7): {
              ibF[c].i1=ip[3]; ibF[c].j1 = jp[3]; ibF[c].k1 = kp[3];
              ibF[c].i2=ip[4]; ibF[c].j2 = jp[4]; ibF[c].k2 = kp[4];
              ibF[c].i3=ip[0]; ibF[c].j3 = jp[0]; ibF[c].k3 = kp[0];
              break;
            }
        }

        if (ibF[c].imode >=0)
        {
            continue;
        }

        if (ibF[c].imode < 0)
        {

            PetscPrintf(PETSC_COMM_SELF, "imode less than 0\n");

            Cmpnts ptmp;
            if (i==1 || i==mx-2 || j==1 || j==my-2)
            {
                if(eNorm.z > 0)
                {
                    ptmp = nSet(cent[k+1][j][i]);
                    ibF[c].dB = nMag(nSub(cent[k][j][i], ptmp));
                    ibF[c].cr1 = 1.;
                    ibF[c].cr2 = 0.;
                    ibF[c].cr3 = 0.;
                    ibF[c].i1 = i;
                    ibF[c].j1 = j;
                    ibF[c].k1 = k+1;

                    ibF[c].i2 = i;
                    ibF[c].j2 = j;
                    ibF[c].k2 = k+1;
                    ibF[c].i3 = i;
                    ibF[c].j3 = j;
                    ibF[c].k3 = k+1;
                }
                else
                {
                    ptmp = nSet(cent[k-1][j][i]);
                    ibF[c].dB = nMag(nSub(cent[k][j][i], ptmp));

                    ibF[c].cr1 = 1.;
                    ibF[c].cr2 = 0.;
                    ibF[c].cr3 = 0.;
                    ibF[c].i1 = i;
                    ibF[c].j1 = j;
                    ibF[c].k1 = k-1;

                    ibF[c].i2 = i;
                    ibF[c].j2 = j;
                    ibF[c].k2 = k-1;
                    ibF[c].i3 = i;
                    ibF[c].j3 = j;
                    ibF[c].k3 = k-1;
                }
            }
            else if (k==1 || k==mz-2)
            {
                ptmp = nSet(cent[k][j+1][i]);
                ibF[c].dB = nMag(nSub(cent[k][j][i], ptmp));

                ibF[c].cr1 = 1.;
                ibF[c].cr2 = 0.;
                ibF[c].cr3 = 0.;
                ibF[c].i1 = i;
                ibF[c].j1 = j+1;
                ibF[c].k1 = k;

                ibF[c].i2 = i;
                ibF[c].j2 = j+1;
                ibF[c].k2 = k;
                ibF[c].i3 = i;
                ibF[c].j3 = j+1;
                ibF[c].k3 = k;
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}
//***************************************************************************************************************//

PetscErrorCode interceptionPt(Cmpnts pCoor, Cmpnts pc[9], Cmpnts eNorm, ibmFluidCell *ibF)
{
    PetscInt 	triangles[3][8];
    Cmpnts   	p1, p2, p3;
    PetscReal	dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3, d;
    PetscReal	rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3;

    Cpt2D		pj1, pj2, pj3, pjp;
    PetscInt	cell, flag;

    PetscInt	i;
    Cmpnts	    pint;                                                // Interception point
    PetscReal	nfxt, nfyt, nfzt;

    ibF->imode = -100;

    triangles[0][0] = 0; triangles[1][0] = 1; triangles[2][0] = 4;
    triangles[0][1] = 1; triangles[1][1] = 2; triangles[2][1] = 4;
    triangles[0][2] = 2; triangles[1][2] = 4; triangles[2][2] = 5;
    triangles[0][3] = 4; triangles[1][3] = 5; triangles[2][3] = 8;
    triangles[0][4] = 4; triangles[1][4] = 7; triangles[2][4] = 8;
    triangles[0][5] = 4; triangles[1][5] = 6; triangles[2][5] = 7;
    triangles[0][6] = 3; triangles[1][6] = 4; triangles[2][6] = 6;
    triangles[0][7] = 3; triangles[1][7] = 4; triangles[2][7] = 0;

    for (i=0; i<8; i++)
    {
        p1 = pc[triangles[0][i]]; p2 = pc[triangles[1][i]], p3 = pc[triangles[2][i]];

        dx1 = pCoor.x - p1.x; dy1 = pCoor.y - p1.y; dz1 = pCoor.z - p1.z;
        dx2 = p2.x - p1.x; dy2 = p2.y - p1.y; dz2 = p2.z - p1.z;
        dx3 = p3.x - p1.x; dy3 = p3.y - p1.y; dz3 = p3.z - p1.z;

        d =  (eNorm.x * (dy2 * dz3 - dz2 * dy3)
            - eNorm.y * (dx2 * dz3 - dz2 * dx3)
            + eNorm.z * (dx2 * dy3 - dy2 * dx3));


        if (fabs(d) > 1.e-10)
        {
            d = -(dx1 * (dy2 * dz3 - dz2 * dy3)
                - dy1 * (dx2 * dz3 - dz2 * dx3)
                + dz1 * (dx2 * dy3 - dy2 * dx3)) / d;

            if (d>0)
            {
                pint.x = pCoor.x + d * eNorm.x;
                pint.y = pCoor.y + d * eNorm.y;
                pint.z = pCoor.z + d * eNorm.z;

                rx1 = p2.x - p1.x; ry1 = p2.y - p1.y; rz1 = p2.z - p1.z;
                rx2 = p3.x - p1.x; ry2 = p3.y - p1.y; rz2 = p3.z - p1.z;

                nfxt = ry1 * rz2 - rz1 * ry2;
                nfyt = -rx1 * rz2 + rz1 * rx2;
                nfzt = rx1 * ry2 - ry1 * rx2;

                flag = isPointInTriangle(pint, p1, p2, p3, eNorm);

                if (flag >= 0)
                {

                    cell = i;

                    if (fabs(nfxt) >= fabs(nfyt) && fabs(nfxt)>= fabs(nfzt))
                    {
                        pjp.x = pint.y; pjp.y = pint.z;
                        pj1.x = p1.y;   pj1.y = p1.z;
                        pj2.x = p2.y;   pj2.y = p2.z;
                        pj3.x = p3.y;   pj3.y = p3.z;
                        triangleIntpBg(pjp, pj1, pj2, pj3, ibF);
                    }
                    else if	(fabs(nfyt) >= fabs(nfxt) && fabs(nfyt)>= fabs(nfzt))
                    {
                      pjp.x = pint.x; pjp.y = pint.z;
                      pj1.x = p1.x;   pj1.y = p1.z;
                      pj2.x = p2.x;   pj2.y = p2.z;
                      pj3.x = p3.x;   pj3.y = p3.z;
                      triangleIntpBg(pjp, pj1, pj2, pj3, ibF);
                    }
                    else if (fabs(nfzt) >= fabs(nfyt) && fabs(nfzt)>= fabs(nfxt))
                    {
                      pjp.x = pint.y; pjp.y = pint.x;
                      pj1.x = p1.y;   pj1.y = p1.x;
                      pj2.x = p2.y;   pj2.y = p2.x;
                      pj3.x = p3.y;   pj3.y = p3.x;
                      triangleIntpBg(pjp, pj1, pj2, pj3, ibF);
                    }

                    ibF->dB = nMag(nSub(pint, pCoor));
                    ibF->imode = cell;
                    return(0);
                }
            }
        }
    }

    return(0);
}
//***************************************************************************************************************//

PetscErrorCode findClosestIBMElement2Solid(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      b, c;
    Cmpnts        ***cent, dist;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    PetscMPIInt        nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt        rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {
        // smallest bounding sphere algorithm
        elementBoundingSphere(ibm->ibmBody[b]);
    }

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        PetscInt      n1, n2, n3, vertexId;
        PetscInt      cellMin = -100;
        Cmpnts        dis, elemNorm, p1, p2, p3;
        PetscReal     d_center, dmin = 1.0e20, d, t, tmin;
        Cmpnts        pmin, po, pj;
        PetscReal     normProj;                             // normal projection of point to ibm mesh element
        PetscInt      bodyID;
        word          closestType;

        // loop through the ibm bodies
        for (b = 0; b < ibm->numBodies; b++)
        {
            ibmObject     *ibmBody = ibm->ibmBody[b];
            boundingBox   *ibBox   = ibmBody->bound;                         // bounding box of the ibm body
            ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
            elementBox    *eBox    = ibmBody->eBox;

            Cmpnts        *qvec   = ibMsh->eQVec;
            PetscReal     *rvec   = ibMsh->eRVec;

            PetscInt      *nv1   = ibMsh->nID1, *nv2  = ibMsh->nID2, *nv3  = ibMsh->nID3;
            Cmpnts        *eN    = ibMsh->eN;
            Cmpnts        *nCoor = ibMsh->nCoor;

            //check if processor controls this ibm body
            if(ibmBody->ibmControlled)
            {
                //loop through the IBM elements
                for(PetscInt e = 0; e < ibMsh->elems; e++)
                {
                    //this processor controls this ibm element
                    if(eBox->thisElemControlled[e])
                    {
                        n1 = nv1[e];
                        n2 = nv2[e];
                        n3 = nv3[e];

                        elemNorm = eN[e];

                        p1 = nCoor[n1];
                        p2 = nCoor[n2];
                        p3 = nCoor[n3];

                        dis = nSub(cent[k][j][i], p1);
                        normProj = nDot(dis, elemNorm);

                        if (fabs(normProj) < 1.e-10) normProj = 1.e-10;

                        if(normProj < 0)
                        {
                            //find the element whose bounding sphere is closest to cent[k][j][i]
                            d_center = nMag(nSub(cent[k][j][i], qvec[e]));

                            if(d_center - rvec[e] < dmin)
                            {
                                dmin = d_center - rvec[e];
                                cellMin = e;
                                bodyID  = b;
                            }
                        }

                    }
                }
            }
        }

        //no ibm elements within the processor buffer zone, do full search
        if (cellMin == -100)
        {
            char error[512];
            sprintf(error, "Rank %d, total ibm fluid cells = %ld, Nearest Cell Searching Error for cell %ld %ld %ld, coordinate %lf %lf %lf\n",rank, ibm->numIBMFluid, k, j, i, cent[k][j][i].x, cent[k][j][i].y, cent[k][j][i].z);
            fatalErrorInFunction("findClosestIBMElement",  error);
        }

        //for that element find if cent[k][j][i] is closest to the element face, edge or vertex
        ibmObject     *ibmBody = ibm->ibmBody[bodyID];
        boundingBox   *ibBox   = ibmBody->bound;                         // bounding box of the ibm body
        ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
        elementBox    *eBox    = ibmBody->eBox;

        PetscInt      *nv1   = ibMsh->nID1, *nv2  = ibMsh->nID2, *nv3  = ibMsh->nID3;
        Cmpnts        *eN    = ibMsh->eN;
        Cmpnts        *nCoor = ibMsh->nCoor;
        dmin = 1.0e20;

        //check if processor controls this ibm body
        if(ibmBody->ibmControlled)
        {
            //this processor controls this ibm element
            if(eBox->thisElemControlled[cellMin])
            {
                n1 = nv1[cellMin];
                n2 = nv2[cellMin];
                n3 = nv3[cellMin];

                elemNorm = eN[cellMin];

                p1 = nCoor[n1];
                p2 = nCoor[n2];
                p3 = nCoor[n3];

                dis = nSub(cent[k][j][i], p1);
                normProj = nDot(dis, elemNorm);

                if (fabs(normProj) < 1.e-10) normProj = 1.e-10;

                if(normProj < 0)
                {
                    dis = nScale(normProj, elemNorm);
                    pj  = nSub(cent[k][j][i], dis);

                    // The projected point is inside the triangle
                    if(isPointInTriangle(pj, p1, p2, p3, elemNorm) == 1)
                    {
                        pmin        = pj;
                        dmin        = fabs(normProj);
                        closestType = "face";
                    }
                    // The projected point is outside the triangle
                    else
                    {
                        // po is the approximated projected point and d is the distance to that point from cent[k][j][i]
                        disP2Line(cent[k][j][i], p1, p2, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n1;
                            }
                            else
                            {
                                vertexId = n2;
                            }
                        }

                        disP2Line(cent[k][j][i], p2, p3, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n2;
                            }
                            else
                            {
                                vertexId = n3;
                            }
                        }

                        disP2Line(cent[k][j][i], p3, p1, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n3;
                            }
                            else
                            {
                                vertexId = n1;
                            }
                        }

                        closestType = "edgeVertex";
                    }
                }
            }
        }

        //set the closest ibm mesh element to the ibm fluid cell
        ibF[c].closestElem = cellMin;
        ibF[c].pMin = pmin;
        ibF[c].minDist = dmin;
        ibF[c].bodyID = bodyID;

        // compute the element normal based on the closestType
        if(closestType == "face")
        {
            ibF[c].normal = nSet(ibMsh->eN[cellMin]);
        }
        else if(closestType == "edgeVertex")
        {
            //compute the angle averaged normal
            if(tmin <=0 || tmin >=1)
            {
                ibF[c].normal = nSet(computeVertexAverageNormal(ibMsh, vertexId));
            }
            else
            {
                ibF[c].normal = nSet(computeEdgeAverageNormal(ibMsh, vertexId, cellMin));
            }

        }
        else
        {
            char error[512];
            sprintf(error, "closestType not set. possible parallelization error");
            fatalErrorInFunction("findClosestIBMElement",  error);
        }

        Cpt2D     pjp, pj1, pj2, pj3;

        elemNorm  = ibMsh->eN[cellMin];
        n1        = ibMsh->nID1[cellMin];
        n2        = ibMsh->nID2[cellMin];
        n3        = ibMsh->nID3[cellMin];

        p1        = ibMsh->nCoor[n1];
        p2        = ibMsh->nCoor[n2];
        p3        = ibMsh->nCoor[n3];

        // to find the ibm interpolation coefficient of the projected point from the ibm element nodes
        if (fabs(elemNorm.x) >= fabs(elemNorm.y) && fabs(elemNorm.x) >= fabs(elemNorm.z))
        {
            pjp.x = pmin.y;
            pjp.y = pmin.z;
            pj1.x = p1.y;
            pj1.y = p1.z;
            pj2.x = p2.y;
            pj2.y = p2.z;
            pj3.x = p3.y;
            pj3.y = p3.z;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }
        else if (fabs(elemNorm.y) >= fabs(elemNorm.x) && fabs(elemNorm.y) >= fabs(elemNorm.z))
        {
            pjp.x = pmin.x;
            pjp.y = pmin.z;
            pj1.x = p1.x;
            pj1.y = p1.z;
            pj2.x = p2.x;
            pj2.y = p2.z;
            pj3.x = p3.x;
            pj3.y = p3.z;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }
        else if (fabs(elemNorm.z) >= fabs(elemNorm.y) && fabs(elemNorm.z) >= fabs(elemNorm.x))
        {
            pjp.x = pmin.y;
            pjp.y = pmin.x;
            pj1.x = p1.y;
            pj1.y = p1.x;
            pj2.x = p2.y;
            pj2.y = p2.x;
            pj3.x = p3.y;
            pj3.y = p3.x;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }

    }

    MPI_Barrier(mesh->MESH_COMM);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findClosestIBMElement(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      b, c;
    Cmpnts        ***cent, dist;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    PetscMPIInt        nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt        rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {
        // smallest bounding sphere algorithm
        elementBoundingSphere(ibm->ibmBody[b]);
    }

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        PetscInt      n1, n2, n3, vertexId;
        PetscInt      cellMin = -100;
        Cmpnts        dis, elemNorm, p1, p2, p3;
        PetscReal     d_center, dmin = 1.0e20, d, t, tmin;
        Cmpnts        pmin, po, pj;
        PetscReal     normProj;                             // normal projection of point to ibm mesh element
        PetscInt      bodyID;
        word          closestType;

        // loop through the ibm bodies
        for (b = 0; b < ibm->numBodies; b++)
        {
            ibmObject     *ibmBody = ibm->ibmBody[b];
            boundingBox   *ibBox   = ibmBody->bound;                         // bounding box of the ibm body
            ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
            elementBox    *eBox    = ibmBody->eBox;

            Cmpnts        *qvec   = ibMsh->eQVec;
            PetscReal     *rvec   = ibMsh->eRVec;

            PetscInt      *nv1   = ibMsh->nID1, *nv2  = ibMsh->nID2, *nv3  = ibMsh->nID3;
            Cmpnts        *eN    = ibMsh->eN;
            Cmpnts        *nCoor = ibMsh->nCoor;

            //check if processor controls this ibm body
            if(ibmBody->ibmControlled)
            {
                //loop through the IBM elements
                for(PetscInt e = 0; e < ibMsh->elems; e++)
                {
                    //this processor controls this ibm element
                    if(eBox->thisElemControlled[e])
                    {
                        n1 = nv1[e];
                        n2 = nv2[e];
                        n3 = nv3[e];

                        elemNorm = eN[e];

                        p1 = nCoor[n1];
                        p2 = nCoor[n2];
                        p3 = nCoor[n3];

                        dis = nSub(cent[k][j][i], p1);
                        normProj = nDot(dis, elemNorm);

                        if (fabs(normProj) < 1.e-10) normProj = 1.e-10;

                        if(normProj > 0)
                        {
                            //find the element whose bounding sphere is closest to cent[k][j][i]
                            d_center = nMag(nSub(cent[k][j][i], qvec[e]));

                            if(d_center - rvec[e] < dmin)
                            {
                                dmin = d_center - rvec[e];
                                cellMin = e;
                                bodyID  = b;
                            }
                        }

                    }
                }
            }
        }

        //no ibm elements within the processor buffer zone, do full search
        if (cellMin == -100)
        {
            char error[512];
            sprintf(error, "Rank %d, total ibm fluid cells = %ld, Nearest Cell Searching Error for cell %ld %ld %ld, coordinate %lf %lf %lf\n",rank, ibm->numIBMFluid, k, j, i, cent[k][j][i].x, cent[k][j][i].y, cent[k][j][i].z);
            fatalErrorInFunction("findClosestIBMElement",  error);
        }

        //for that element find if cent[k][j][i] is closest to the element face, edge or vertex
        ibmObject     *ibmBody = ibm->ibmBody[bodyID];
        boundingBox   *ibBox   = ibmBody->bound;                         // bounding box of the ibm body
        ibmMesh       *ibMsh   = ibmBody->ibMsh;                         // pointer to the ibm body mesh
        elementBox    *eBox    = ibmBody->eBox;

        PetscInt      *nv1   = ibMsh->nID1, *nv2  = ibMsh->nID2, *nv3  = ibMsh->nID3;
        Cmpnts        *eN    = ibMsh->eN;
        Cmpnts        *nCoor = ibMsh->nCoor;
        dmin = 1.0e20;

        //check if processor controls this ibm body
        if(ibmBody->ibmControlled)
        {
            //this processor controls this ibm element
            if(eBox->thisElemControlled[cellMin])
            {
                n1 = nv1[cellMin];
                n2 = nv2[cellMin];
                n3 = nv3[cellMin];

                elemNorm = eN[cellMin];

                p1 = nCoor[n1];
                p2 = nCoor[n2];
                p3 = nCoor[n3];

                dis = nSub(cent[k][j][i], p1);
                normProj = nDot(dis, elemNorm);

                if (fabs(normProj) < 1.e-10) normProj = 1.e-10;

                if(normProj > 0)
                {
                    dis = nScale(normProj, elemNorm);
                    pj  = nSub(cent[k][j][i], dis);

                    // The projected point is inside the triangle
                    if(isPointInTriangle(pj, p1, p2, p3, elemNorm) == 1)
                    {
                        pmin        = pj;
                        dmin        = normProj;
                        closestType = "face";
                    }
                    // The projected point is outside the triangle
                    else
                    {
                        // po is the approximated projected point and d is the distance to that point from cent[k][j][i]
                        disP2Line(cent[k][j][i], p1, p2, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n1;
                            }
                            else
                            {
                                vertexId = n2;
                            }
                        }

                        disP2Line(cent[k][j][i], p2, p3, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n2;
                            }
                            else
                            {
                                vertexId = n3;
                            }
                        }

                        disP2Line(cent[k][j][i], p3, p1, &po, &d, &t);

                        if (d < dmin)
                        {
                            dmin = d;
                            pmin = po;
                            tmin = t;

                            if(t < 1)
                            {
                                vertexId = n3;
                            }
                            else
                            {
                                vertexId = n1;
                            }
                        }

                        closestType = "edgeVertex";
                    }
                }
            }
        }

        //set the closest ibm mesh element to the ibm fluid cell
        ibF[c].closestElem = cellMin;
        ibF[c].pMin = pmin;
        ibF[c].minDist = dmin;
        ibF[c].bodyID = bodyID;

        // compute the element normal based on the closestType
        if(closestType == "face")
        {
            ibF[c].normal = nSet(ibMsh->eN[cellMin]);
        }
        else if(closestType == "edgeVertex")
        {
            //compute the angle averaged normal
            if(tmin <=0 || tmin >=1)
            {
                ibF[c].normal = nSet(computeVertexAverageNormal(ibMsh, vertexId));
            }
            else
            {
                ibF[c].normal = nSet(computeEdgeAverageNormal(ibMsh, vertexId, cellMin));
            }

        }
        else
        {
            char error[512];
            sprintf(error, "closestType not set. possible parallelization error");
            fatalErrorInFunction("findClosestIBMElement",  error);
        }

        Cpt2D     pjp, pj1, pj2, pj3;

        elemNorm  = ibMsh->eN[cellMin];
        n1        = ibMsh->nID1[cellMin];
        n2        = ibMsh->nID2[cellMin];
        n3        = ibMsh->nID3[cellMin];

        p1        = ibMsh->nCoor[n1];
        p2        = ibMsh->nCoor[n2];
        p3        = ibMsh->nCoor[n3];

        // to find the ibm interpolation coefficient of the projected point from the ibm element nodes
        if (fabs(elemNorm.x) >= fabs(elemNorm.y) && fabs(elemNorm.x) >= fabs(elemNorm.z))
        {
            pjp.x = pmin.y;
            pjp.y = pmin.z;
            pj1.x = p1.y;
            pj1.y = p1.z;
            pj2.x = p2.y;
            pj2.y = p2.z;
            pj3.x = p3.y;
            pj3.y = p3.z;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }
        else if (fabs(elemNorm.y) >= fabs(elemNorm.x) && fabs(elemNorm.y) >= fabs(elemNorm.z))
        {
            pjp.x = pmin.x;
            pjp.y = pmin.z;
            pj1.x = p1.x;
            pj1.y = p1.z;
            pj2.x = p2.x;
            pj2.y = p2.z;
            pj3.x = p3.x;
            pj3.y = p3.z;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }
        else if (fabs(elemNorm.z) >= fabs(elemNorm.y) && fabs(elemNorm.z) >= fabs(elemNorm.x))
        {
            pjp.x = pmin.y;
            pjp.y = pmin.x;
            pj1.x = p1.y;
            pj1.y = p1.x;
            pj2.x = p2.y;
            pj2.y = p2.x;
            pj3.x = p3.y;
            pj3.y = p3.x;
            triangleIntp(pjp, pj1, pj2, pj3, &(ibF[c]));
        }

    }

    MPI_Barrier(mesh->MESH_COMM);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode CurvibInterpolationTriangular(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    flags_        *flags = ibm->access->flags;
    constants_    *cst   = ibm->access->constants;
    clock_        *clock = ibm->access->clock;
    peqn_         *peqn  = ibm->access->peqn;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k, ii, kk, jj, k1, j1, i1;
    PetscInt      b, c, s;
    Cmpnts        ***ucat, ***lucat, ***cent;
    PetscReal     ***lt, ***temp, ***nvert, ***aj, dudt;
    Cmpnts        bPt;
    Cmpnts        bPtVel, ibmPtVel, ibmPtVelPrev;
    PetscReal     bPtTemp, ibmPtTemp, cellSize;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lAj, &aj);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
        DMDAVecGetArray(da, teqn->Tmprt, &temp);
    }

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    PetscReal    sc, sb, ustar;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmMesh   *ibMsh = ibm->ibmBody[ibF[c].bodyID]->ibMsh;
        PetscInt   cElem = ibF[c].closestElem;
        Cmpnts	   eNorm = ibMsh->eN[cElem];
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];

        PetscInt     ip1 = ibF[c].i1, jp1 = ibF[c].j1, kp1 = ibF[c].k1;
        PetscInt     ip2 = ibF[c].i2, jp2 = ibF[c].j2, kp2 = ibF[c].k2;
        PetscInt     ip3 = ibF[c].i3, jp3 = ibF[c].j3, kp3 = ibF[c].k3;

        // interpolate the velocity of the projected point on the IBM solid element from its nodes
        ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                     + ibMsh->nU[n2].x * ibF[c].cs2
                     + ibMsh->nU[n3].x * ibF[c].cs3;

        ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                     + ibMsh->nU[n2].y * ibF[c].cs2
                     + ibMsh->nU[n3].y * ibF[c].cs3;

        ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                     + ibMsh->nU[n2].z * ibF[c].cs2
                     + ibMsh->nU[n3].z * ibF[c].cs3;

        // interpolate the previous velocity
        ibmPtVelPrev.x =   ibMsh->nUPrev[n1].x * ibF[c].cs1
                         + ibMsh->nUPrev[n2].x * ibF[c].cs2
                         + ibMsh->nUPrev[n3].x * ibF[c].cs3;

        ibmPtVelPrev.y =   ibMsh->nUPrev[n1].y * ibF[c].cs1
                         + ibMsh->nUPrev[n2].y * ibF[c].cs2
                         + ibMsh->nUPrev[n3].y * ibF[c].cs3;

        ibmPtVelPrev.z =   ibMsh->nUPrev[n1].z * ibF[c].cs1
                         + ibMsh->nUPrev[n2].z * ibF[c].cs2
                         + ibMsh->nUPrev[n3].z * ibF[c].cs3;

         // save the element acceleration term (dudt . elementNormal)
         dudt = nDot(nScale(1.0/clock->dtOld, nSub(ibmPtVel, ibmPtVelPrev)), eNorm);

        bPtVel.x = lucat[kp1][jp1][ip1].x * ibF[c].cr1 + lucat[kp2][jp2][ip2].x * ibF[c].cr2 + lucat[kp3][jp3][ip3].x * ibF[c].cr3;
        bPtVel.y = lucat[kp1][jp1][ip1].y * ibF[c].cr1 + lucat[kp2][jp2][ip2].y * ibF[c].cr2 + lucat[kp3][jp3][ip3].y * ibF[c].cr3;
        bPtVel.z = lucat[kp1][jp1][ip1].z * ibF[c].cr1 + lucat[kp2][jp2][ip2].z * ibF[c].cr2 + lucat[kp3][jp3][ip3].z * ibF[c].cr3;

        if (flags->isTeqnActive)
        {
            if(flags->isAblActive)
            {
                ibmPtTemp = ibm->access->abl->tRef;
            }
            else
            {
                ibmPtTemp = ibm->access->constants->tRef;
            }
        }

        bPt = nScale(ibF[c].dB, eNorm);
        mSum(bPt, cent[k][j][i]);

        sb = nDot(nSub(cent[k][j][i], ibF[c].pMin), eNorm);
        sc = nDot(nSub(bPt, ibF[c].pMin), eNorm);

        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
        {
            if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
            {
                Cabot *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot;

                if(wm->roughness > 1.0e-12)
                {
                    wallFunctionCabotRoughness(cst->nu, wm->roughness, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
                }
                else
                {
                    wallFunctionCabot(cst->nu, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
                }
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
            {
                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;

                wallFunctionSchumann(cst->nu, sc, sb, wm->roughness, wm->kappa, ibmPtVel,
                                        bPtVel, &ucat[k][j][i], &ustar, eNorm);

            } 
        }
        else if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "slip")
        {
            slipBC(sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], eNorm);
        }
        else if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "noSlip")
        {
            ucat[k][j][i].x = (sb/sc) * bPtVel.x + (1.0 - (sb/sc)) * ibmPtVel.x;
            ucat[k][j][i].y= (sb/sc) * bPtVel.y + (1.0 - (sb/sc)) * ibmPtVel.y;
            ucat[k][j][i].z = (sb/sc) * bPtVel.z + (1.0 - (sb/sc)) * ibmPtVel.z;
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
    }

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    return (0);
}
//***************************************************************************************************************//

PetscErrorCode CurvibInterpolationInternalCell(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    peqn_         *peqn  = ibm->access->peqn;
    flags_        *flags = ibm->access->flags;
    constants_    *cst   = ibm->access->constants;
    clock_        *clock = ibm->access->clock;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys   = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs   = info.gzs, gze = info.gzs + info.gzm;

    PetscInt      i, j, k, ii, kk, jj;
    PetscInt      b, c, s;
    cellIds       initCp;
    Cmpnts        ***ucat, ***lucat, ***cent, ***csi, ***eta, ***zet;
    PetscReal     ***lt, ***temp, ***aj, ***nvert, ***lp, ***p;
    Cmpnts        bPt1, bPt2, bPtInit;
    Cmpnts        bPt1Vel, bPt2Vel, ibmPtVel, ibmPtVelPrev;
    PetscReal     nfMag, cellSize, bPt1Temp, bPt2Temp, ibmPtTemp, bPt1Pres, bPt2Pres, ibmPtPres, elemDist;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, peqn->lP, &lp);
    DMDAVecGetArray(da, peqn->P, &p);
    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
        DMDAVecGetArray(da, teqn->Tmprt, &temp);
    }

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    PetscReal    sc, sb, sd, ustar;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmMesh   *ibMsh = ibm->ibmBody[ibF[c].bodyID]->ibMsh;
        PetscInt   cElem = ibF[c].closestElem;
        Cmpnts	   eNorm = ibF[c].normal;
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal  roughness;

        /********************************************************************************************************
            Interpolation of the first background point field
        /*******************************************************************************************************/

        // background mesh projection point
        cellSize = ibm->interpDist;

        // distance from cell to the IBM element
        elemDist = fabs(nDot(nSub(ibF[c].pMin, cent[k][j][i]), eNorm));

        bPt1 = nScale(cellSize + elemDist, eNorm);
        mSum(bPt1, cent[k][j][i]);

        // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
        PetscInt    i1, j1, k1;
        PetscReal   dmin = 10e10, d;
        PetscInt    ic, jc, kc, setFlag = 0;

        // must be close so don't loop over all cells
        for (k1=k-1; k1<k+2; k1++)
        for (j1=j-1; j1<j+2; j1++)
        for (i1=i-1; i1<i+2; i1++)
        {

            if
            (
                (
                    k1>=1 && k1<mz-1 &&
                    j1>=1 && j1<my-1 &&
                    i1>=1 && i1<mx-1
                ) && (isFluidCell(k1, j1, i1, nvert))
            )
            {
                d = pow((bPt1.x - cent[k1][j1][i1].x), 2) +
                    pow((bPt1.y - cent[k1][j1][i1].y), 2) +
                    pow((bPt1.z - cent[k1][j1][i1].z), 2);

                if
                (
                    d < dmin
                )
                {
                    dmin  = d;
                    ic = i1;
                    jc = j1;
                    kc = k1;
                    setFlag = 1;
                }
            }
        }

        if(setFlag == 0)
        {
            char error[512];
            sprintf(error, "no closest fluid cell to ibm fluid %ld %ld %ld\n", k, j, i);
            fatalErrorInFunction("CurvibInterpolation",  error);
        }

        PetscInt intId[6];
        PetscInt intFlag = 1;

        // get the trilinear interpolation cells
        PointInterpolationCells
        (
                mesh,
                bPt1.x, bPt1.y, bPt1.z,
                ic, jc, kc,
                cent,
                intId
        );

        // save the initial closest cell and background point
        initCp.i = ic; initCp.j = jc; initCp.k = kc;
        bPtInit  = nSet(bPt1);

        // ensure none of the interpolation cells are solid cells.
        PetscReal sumDel = 0.0;

        while ( (intFlag == 1) && isInsideBoundingBox(bPt1, mesh->bounds) && (sumDel <= 1.5*cellSize))
        {
            PetscInt ibmCellCtr = 0;
            PetscInt icc, jcc, kcc, setFlag = 0;

            for (PetscInt kk = 0; kk<2; kk++)
            for (PetscInt jj = 0; jj<2; jj++)
            for (PetscInt ii = 0; ii<2; ii++)
            {
                if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                {
                    ibmCellCtr ++;
                }

            }

            if (ibmCellCtr > 0)
            {
                Cmpnts del =  nScale(0.2 * cellSize, eNorm);
                mSum(bPt1, del);
                dmin = 10e10;

                sumDel += nMag(del);

                for (k1=kc-1; k1<kc+2; k1++)
                {
                    //check processor ghost bounds
                    if (k1 < gzs || k1 >= gze) {intFlag = 2; break;}
                    for (j1=jc-1; j1<jc+2; j1++)
                    {
                        if (j1 < gys || j1 >= gye) {intFlag = 2; break;}
                        for (i1=ic-1; i1<ic+2; i1++)
                        {
                            if (i1 < gxs || i1 >= gxe) {intFlag = 2; break;}

                            if
                            (
                                (
                                    k1>=1 && k1<mz-1 &&
                                    j1>=1 && j1<my-1 &&
                                    i1>=1 && i1<mx-1
                                ) && (isFluidCell(k1, j1, i1, nvert))
                            )
                            {
                                d = pow((bPt1.x - cent[k1][j1][i1].x), 2) +
                                    pow((bPt1.y - cent[k1][j1][i1].y), 2) +
                                    pow((bPt1.z - cent[k1][j1][i1].z), 2);

                                if
                                (
                                    d < dmin
                                )
                                {
                                    dmin  = d;
                                    icc = i1;
                                    jcc = j1;
                                    kcc = k1;
                                    setFlag = 1;
                                }
                            }

                        }
                    }
                }

                if(setFlag == 0)
                {
                    //closest point not set, do not interpolate
                    intFlag = 2;
                }

                if(intFlag == 1)
                {
                    kc = kcc; jc = jcc; ic = icc;

                    PointInterpolationCells
                    (
                            mesh,
                            bPt1.x, bPt1.y, bPt1.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );
                }

            }
            else
            {
                intFlag = 0;

                // trilinear interpolate the velocity at this point
                vectorPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt1.x, bPt1.y, bPt1.z,
                        ic, jc, kc,
                        cent,
                        lucat,
                        bPt1Vel
                );

                if (flags->isTeqnActive)
                {
                    scalarPointLocalVolumeInterpolation
                    (
                            mesh,
                            bPt1.x, bPt1.y, bPt1.z,
                            ic, jc, kc,
                            cent,
                            lt,
                            bPt1Temp
                    );
                }

                scalarPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt1.x, bPt1.y, bPt1.z,
                        ic, jc, kc,
                        cent,
                        lp,
                        bPt1Pres
                );

            }
        }

        // while loop fails
        if(intFlag > 0)
        {
            bPt1Vel = nSet(lucat[initCp.k][initCp.j][initCp.i]);

            bPt1Pres = lp[initCp.k][initCp.j][initCp.i];

            if (flags->isTeqnActive)
            {
                bPt1Temp = lt[initCp.k][initCp.j][initCp.i];
            }

            bPt1 = nSet(bPtInit);

            ic = initCp.i;
            jc = initCp.j;
            kc = initCp.k;
        }

         /********************************************************************************************************
             Interpolation of the second background point field
         /*******************************************************************************************************/

         bPt2 = nScale(2.0 * cellSize, eNorm);
         mSum(bPt2, bPt1);

         // search for the closest background mesh cell center (it will be close to current closest fluid cell kc,jc,ic)
         dmin = 10e10; setFlag = 0;
         PetscInt    ic2, jc2, kc2;

         // must be close so don't loop over all cells
         for (k1=kc-1; k1<kc+2; k1++)
         for (j1=jc-1; j1<jc+2; j1++)
         for (i1=ic-1; i1<ic+2; i1++)
         {

             if
             (
                 (
                     k1>=1 && k1<mz-1 &&
                     j1>=1 && j1<my-1 &&
                     i1>=1 && i1<mx-1
                 ) && (isFluidCell(k1, j1, i1, nvert))
             )
             {
                 d = pow((bPt2.x - cent[k1][j1][i1].x), 2) +
                     pow((bPt2.y - cent[k1][j1][i1].y), 2) +
                     pow((bPt2.z - cent[k1][j1][i1].z), 2);

                 if
                 (
                     d < dmin
                 )
                 {
                     dmin  = d;
                     ic2 = i1;
                     jc2 = j1;
                     kc2 = k1;
                     setFlag = 1;
                 }
             }
         }

         if(setFlag == 0)
         {
             char error[512];
             sprintf(error, "no closest fluid cell to background point 2\n");
             fatalErrorInFunction("CurvibInterpolation",  error);
         }

         //check for the processor bound
         intFlag = 1;

         if (kc2 < gzs || kc2 >= gze) {intFlag = 2;}
         if (jc2 < gys || jc2 >= gye) {intFlag = 2;}
         if (ic2 < gxs || ic2 >= gxe) {intFlag = 2;}

         PetscInt ibmCellCtr = 0;

         if(intFlag == 1)
         {
             // get the trilinear interpolation cells
             PointInterpolationCells
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     intId
             );

             //check that the interpolation cells are not ibm solid
             for (PetscInt kk = 0; kk<2; kk++)
             for (PetscInt jj = 0; jj<2; jj++)
             for (PetscInt ii = 0; ii<2; ii++)
             {
                 if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                 {
                     ibmCellCtr ++;
                 }

             }
         }

         if(intFlag == 2)
         {
             bPt2Vel = nSet(bPt1Vel);

             bPt2Pres =  bPt1Pres;
             if (flags->isTeqnActive)
             {
                 bPt2Temp = bPt1Temp;
             }
         }
         else if(intFlag == 1 && ibmCellCtr > 0)
         {
             bPt2Vel = nSet(lucat[kc2][jc2][ic2]);
             
             bPt2Pres = lp[kc2][jc2][ic2];
             
             if (flags->isTeqnActive)
             {
                 bPt2Temp = lt[kc2][jc2][ic2];
             }
         }
         else
         {
             // trilinear interpolate the velocity at this point
             vectorPointLocalVolumeInterpolation
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     lucat,
                     bPt2Vel
             );

             scalarPointLocalVolumeInterpolation
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     lp,
                     bPt2Pres
             );

             if (flags->isTeqnActive)
             {
                 scalarPointLocalVolumeInterpolation
                 (
                         mesh,
                         bPt2.x, bPt2.y, bPt2.z,
                         ic2, jc2, kc2,
                         cent,
                         lt,
                         bPt2Temp
                 );
             }
         }

         /********************************************************************************************************
             Interpolation of the IBM surface field
         /*******************************************************************************************************/

        // interpolate the velocity of the projected point on the IBM solid element from its nodes
        ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                  + ibMsh->nU[n2].x * ibF[c].cs2
                  + ibMsh->nU[n3].x * ibF[c].cs3;

        ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                  + ibMsh->nU[n2].y * ibF[c].cs2
                  + ibMsh->nU[n3].y * ibF[c].cs3;

        ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                  + ibMsh->nU[n2].z * ibF[c].cs2
                  + ibMsh->nU[n3].z * ibF[c].cs3;

         // interpolate the previous velocity
         ibmPtVelPrev.x =   ibMsh->nUPrev[n1].x * ibF[c].cs1
                          + ibMsh->nUPrev[n2].x * ibF[c].cs2
                          + ibMsh->nUPrev[n3].x * ibF[c].cs3;

         ibmPtVelPrev.y =   ibMsh->nUPrev[n1].y * ibF[c].cs1
                          + ibMsh->nUPrev[n2].y * ibF[c].cs2
                          + ibMsh->nUPrev[n3].y * ibF[c].cs3;

         ibmPtVelPrev.z =   ibMsh->nUPrev[n1].z * ibF[c].cs1
                          + ibMsh->nUPrev[n2].z * ibF[c].cs2
                          + ibMsh->nUPrev[n3].z * ibF[c].cs3;

        /********************************************************************************************************
            Interpolation of the IBM fluid cell field
        /*******************************************************************************************************/
         // point placement b.......wall........c.......d
         //find the distance from the points b and c along the element normal
         sb = -elemDist;
         sc = nDot(nSub(bPt1, cent[k][j][i]), eNorm) - elemDist;
         sd = sc + nDot(nSub(bPt2, bPt1), eNorm); 

        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
        {
            if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
            {
                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;

                wallShearVelocityBCQuadratic(cst->nu, sd, sc, sb, wm->roughness, wm->kappa, ibmPtVel,
                            bPt2Vel, bPt1Vel, &ucat[k][j][i], &ustar, eNorm);
            }
            else if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
            {
                Cabot *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot;

                wallShearVelocityBCQuadratic(cst->nu, sd, sc, sb, wm->roughness, wm->kappa, ibmPtVel,
                            bPt2Vel, bPt1Vel, &ucat[k][j][i], &ustar, eNorm);
            }
        }
        else 
        {
            char error[512];
            sprintf(error, "wall shear model is on, however velocity BC is not velocityWallfunction \n");
            fatalErrorInFunction("CurvibInterpolationInternalCell", error);
        }

        
         //pressure and temperature boundary condition - neumann
         PetscReal aTemp, bTemp, cTemp, aP, bP, cP;

         bP = nDot(nScale(-1.0/clock->dtOld, nSub(ibmPtVel, ibmPtVelPrev)), eNorm);
         aP = (bPt2Pres - bPt1Pres - bP * cellSize)/(sd*sd - sc*sc);
         cP = bPt1Pres - aP*sc*sc - bP*sc;

        //  if(intFlag == 2 || (intFlag == 1 && ibmCellCtr > 0))
        //  {
        //      p[k][j][i] = bP * (sb - sc) + bPt1Pres;
        //  }
        //  else
        //  {
        //      p[k][j][i] = aP * sb * sb + bP * sb + cP;
        //  }

         p[k][j][i] = bP * (sb - sc) + bPt1Pres;
         
         if (flags->isTeqnActive)
         {

            if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "zeroGradient")
            {
                //apply zerogradient boundary condition. effect of boundary through heat flux bc
                temp[k][j][i] =  bPt1Temp;
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "fixedValue")
            {
                temp[k][j][i] =  ibm->ibmBody[ibF[c].bodyID]->fixedTemp;
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "thetaWallFunction")
            {
                if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT==-4)
                {
                    if(teqn->access->clock->it==teqn->access->clock->itStart)
                    {
                        temp[k][j][i] =  bPt1Temp;
                    }
                    else 
                    {
                        Shumann *wm  = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelT->wmShumann;

                        // find interpolation weights for surface temp
                        PetscReal w[2];
                        PetscInt  l[2];

                        findInterpolationWeights(w, l, wm->timeVec, wm->numT, clock->time);

                        PetscReal surfaceTemp = w[0] * wm->surfTemp[l[0]] + w[1] * wm->surfTemp[l[1]];

                        temp[k][j][i] = 2.0*surfaceTemp - bPt1Temp;
                    }
                }
                else 
                {
                    temp[k][j][i] =  bPt1Temp;
                }
            }
            //  bTemp = 0.0;
            //  aTemp = (bPt2Temp - bPt1Temp)/(sd*sd - sc*sc);
            //  cTemp = bPt1Temp - aTemp*sc*sc;

            //  if(intFlag == 2 || (intFlag == 1 && ibmCellCtr > 0))
            //  {
            //     temp[k][j][i] =  bPt1Temp;
            //  }
            //  else
            //  {
            //     temp[k][j][i] = aTemp * sb * sb + cTemp;
            //  }
         }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, peqn->lP, &lp);
    DMDAVecRestoreArray(da, peqn->P, &p);
    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
        DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    }

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd(da, peqn->P, INSERT_VALUES, peqn->lP);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CurvibInterpolation(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    peqn_         *peqn  = ibm->access->peqn;
    flags_        *flags = ibm->access->flags;
    constants_    *cst   = ibm->access->constants;
    clock_        *clock = ibm->access->clock;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys   = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs   = info.gzs, gze = info.gzs + info.gzm;

    PetscInt      i, j, k, ii, kk, jj;
    PetscInt      b, c, s;
    cellIds       initCp;
    Cmpnts        ***ucat, ***lucat, ***cent, ***csi, ***eta, ***zet;
    PetscReal     ***lt, ***temp, ***aj, ***nvert, ***meshTag, ***lp, ***p;
    Cmpnts        bPt, bPtInit;
    Cmpnts        bPtVel, ibmPtVel, ibmPtVelPrev;
    PetscReal     nfMag, cellSize, bPtTemp, bPtPres, elemDist;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da, peqn->lP, &lp);
    DMDAVecGetArray(da, peqn->P, &p);
    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
        DMDAVecGetArray(da, teqn->Tmprt, &temp);
    }

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    PetscReal    sc, sb, ustar, ustarMax = -1e10, ustarMin = 1e10;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmMesh   *ibMsh = ibm->ibmBody[ibF[c].bodyID]->ibMsh;
        PetscInt   cElem = ibF[c].closestElem;
        Cmpnts	   eNorm = ibF[c].normal;
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal  roughness;

        // background mesh projection point is taken 0.5 cell distance along the normal directional of the closest ibm mesh element
        elemDist = fabs(nDot(nSub(ibF[c].pMin, cent[k][j][i]), eNorm));
        cellSize = 1.0 * elemDist;
        bPt = nScale(cellSize, eNorm);
        mSum(bPt, cent[k][j][i]);

        // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
        PetscInt    i1, j1, k1;
        PetscReal   dmin = 10e10, d;
        PetscInt    ic, jc, kc, setFlag = 0;

        // must be close so don't loop over all cells
        for (k1=k-1; k1<k+2; k1++)
        for (j1=j-1; j1<j+2; j1++)
        for (i1=i-1; i1<i+2; i1++)
        {

            if
            (
                (
                    k1>=1 && k1<mz-1 &&
                    j1>=1 && j1<my-1 &&
                    i1>=1 && i1<mx-1
                ) && (!isIBMSolidCell(k1, j1, i1, nvert))
            )
            {
                d = pow((bPt.x - cent[k1][j1][i1].x), 2) +
                    pow((bPt.y - cent[k1][j1][i1].y), 2) +
                    pow((bPt.z - cent[k1][j1][i1].z), 2);

                if
                (
                    d < dmin
                )
                {
                    dmin  = d;
                    ic = i1;
                    jc = j1;
                    kc = k1;
                    setFlag = 1;
                }
            }
        }

        if(setFlag == 0)
        {
            char error[512];
            sprintf(error, "no closest fluid cell to ibm fluid %ld %ld %ld\n", k, j, i);
            fatalErrorInFunction("CurvibInterpolation",  error);
        }

        PetscInt intId[6];
        PetscInt intFlag = 1;

        // get the trilinear interpolation cells
        PointInterpolationCells
        (
                mesh,
                bPt.x, bPt.y, bPt.z,
                ic, jc, kc,
                cent,
                intId
        );

        // save the initial closest cell and background point
        initCp.i = ic; initCp.j = jc; initCp.k = kc;
        bPtInit  = nSet(bPt);

        // max distance checkpoint has to move to be to the closest fluid trilinear interpolation box
        // restrict this to 1.5*ref length
        PetscReal sumDel = 0.0;

        while ( (intFlag == 1) && isInsideBoundingBox(bPt, mesh->bounds) && (sumDel <= 1.5*cellSize))
        {
            PetscInt ibmCellCtr = 0;
            PetscInt icc, jcc, kcc, setFlag = 0;

            for (PetscInt kk = 0; kk<2; kk++)
            for (PetscInt jj = 0; jj<2; jj++)
            for (PetscInt ii = 0; ii<2; ii++)
            {
                if(isIBMSolidCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                {
                    ibmCellCtr ++;
                }

            }

            if (ibmCellCtr > 0)
            {
                Cmpnts del =  nScale(0.2 * cellSize, eNorm);
                mSum(bPt, del);
                dmin = 10e10;

                sumDel += nMag(del);

                for (k1=kc-1; k1<kc+2; k1++)
                {
                    //check processor ghost bounds
                    if (k1 < gzs || k1 >= gze) {intFlag = 2; break;}
                    for (j1=jc-1; j1<jc+2; j1++)
                    {
                        if (j1 < gys || j1 >= gye) {intFlag = 2; break;}
                        for (i1=ic-1; i1<ic+2; i1++)
                        {
                            if (i1 < gxs || i1 >= gxe) {intFlag = 2; break;}

                            if
                            (
                                (
                                    k1>=1 && k1<mz-1 &&
                                    j1>=1 && j1<my-1 &&
                                    i1>=1 && i1<mx-1
                                ) && (!isIBMSolidCell(k1, j1, i1, nvert))
                            )
                            {
                                d = pow((bPt.x - cent[k1][j1][i1].x), 2) +
                                    pow((bPt.y - cent[k1][j1][i1].y), 2) +
                                    pow((bPt.z - cent[k1][j1][i1].z), 2);

                                if
                                (
                                    d < dmin
                                )
                                {
                                    dmin  = d;
                                    icc = i1;
                                    jcc = j1;
                                    kcc = k1;
                                    setFlag = 1;
                                }
                            }

                        }
                    }
                }

                if(setFlag == 0)
                {
                    //closest point not set, do not interpolate
                    intFlag = 2;
                }

                if(intFlag == 1)
                {
                    kc = kcc; jc = jcc; ic = icc;

                    PointInterpolationCells
                    (
                            mesh,
                            bPt.x, bPt.y, bPt.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );
                }

            }
            else
            {
                intFlag = 0;

                // trilinear interpolate the velocity at this point
                vectorPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt.x, bPt.y, bPt.z,
                        ic, jc, kc,
                        cent,
                        lucat,
                        bPtVel
                );

                if (flags->isTeqnActive)
                {
                    scalarPointLocalVolumeInterpolation
                    (
                            mesh,
                            bPt.x, bPt.y, bPt.z,
                            ic, jc, kc,
                            cent,
                            lt,
                            bPtTemp
                    );
                }

                scalarPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt.x, bPt.y, bPt.z,
                        ic, jc, kc,
                        cent,
                        lp,
                        bPtPres
                );
            }
        }

        // while loop fails
        if(intFlag > 0)
        {
            bPtVel = nSet(lucat[initCp.k][initCp.j][initCp.i]);

            bPtPres = lp[initCp.k][initCp.j][initCp.i];

            if (flags->isTeqnActive)
            {
                bPtTemp = lt[initCp.k][initCp.j][initCp.i];
            }

            bPt = nSet(bPtInit);

            ic = initCp.i;
            jc = initCp.j;
            kc = initCp.k;
        }

        // interpolate the velocity of the projected point on the IBM solid element from its nodes
        ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                     + ibMsh->nU[n2].x * ibF[c].cs2
                     + ibMsh->nU[n3].x * ibF[c].cs3;

        ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                     + ibMsh->nU[n2].y * ibF[c].cs2
                     + ibMsh->nU[n3].y * ibF[c].cs3;

        ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                     + ibMsh->nU[n2].z * ibF[c].cs2
                     + ibMsh->nU[n3].z * ibF[c].cs3;

        // interpolate the previous velocity
        ibmPtVelPrev.x =   ibMsh->nUPrev[n1].x * ibF[c].cs1
                         + ibMsh->nUPrev[n2].x * ibF[c].cs2
                         + ibMsh->nUPrev[n3].x * ibF[c].cs3;

        ibmPtVelPrev.y =   ibMsh->nUPrev[n1].y * ibF[c].cs1
                         + ibMsh->nUPrev[n2].y * ibF[c].cs2
                         + ibMsh->nUPrev[n3].y * ibF[c].cs3;

        ibmPtVelPrev.z =   ibMsh->nUPrev[n1].z * ibF[c].cs1
                         + ibMsh->nUPrev[n2].z * ibF[c].cs2
                         + ibMsh->nUPrev[n3].z * ibF[c].cs3;

         //find the distance from the wall of points b and c along the element normal
         sb = elemDist;
         sc = sb + cellSize;

         if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
         {
            if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
            {
                Cabot *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot;

                if(roughness > 1.0e-12)
                {
                    wallFunctionCabotRoughness(cst->nu, wm->roughness, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
                }
                else
                {
                    wallFunctionCabot(cst->nu, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
                }
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -4 || ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -5)
            {
                PetscReal roughness, kappa;

                if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -4)
                {
                    roughness = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmPowerLawAPG->roughness;
                    kappa     = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmPowerLawAPG->kappa;
                }

                if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -5)
                {
                    roughness = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmLogLawAPG->roughness;
                    kappa     = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmLogLawAPG->kappa;
                }

                //find the pressure gradient in the tangential direction
                PetscReal dpdc, dpde, dpdz;
                PetscReal dp_dx, dp_dy, dp_dz;

                PetscReal csi0 = csi[k][j][i].x,
                        csi1 = csi[k][j][i].y,
                        csi2 = csi[k][j][i].z;
                PetscReal eta0 = eta[k][j][i].x,
                        eta1 = eta[k][j][i].y,
                        eta2 = eta[k][j][i].z;
                PetscReal zet0 = zet[k][j][i].x,
                        zet1 = zet[k][j][i].y,
                        zet2 = zet[k][j][i].z;
                PetscReal ajc  = aj[k][j][i];

                Compute_dscalar_center
                (
                    mesh,
                    i, j, k, mx, my, mz, lp, nvert, meshTag, &dpdc, &dpde, &dpdz
                );

                Compute_dscalar_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2,
                    zet0, zet1, zet2, ajc, dpdc, dpde, dpdz, &dp_dx,
                    &dp_dy, &dp_dz
                );

                if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -4)
                {
                    wallFunctionPowerlawAPG(cst->nu, sc, sb, roughness, kappa, ibmPtVel,
                                            bPtVel, &ucat[k][j][i], &ustar, eNorm,
                                            dp_dx, dp_dy, dp_dz);
                }

                if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -5)
                {
                    wallFunctionLogLawAPG(cst->nu, sc, sb, roughness, kappa, ibmPtVel,
                                            bPtVel, &ucat[k][j][i], &ustar, eNorm,
                                            dp_dx, dp_dy, dp_dz);
                }
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
            {
                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;

                wallFunctionSchumann(cst->nu, sc, sb, wm->roughness, wm->kappa, ibmPtVel,
                                bPtVel, &ucat[k][j][i], &ustar, eNorm);

            }
         }
         else if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "slip")
         {
             slipBC(sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], eNorm);
         }
         else if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "noSlip")
         {
             ucat[k][j][i].x = (sb/sc) * bPtVel.x + (1.0 - (sb/sc)) * ibmPtVel.x;
             ucat[k][j][i].y= (sb/sc) * bPtVel.y + (1.0 - (sb/sc)) * ibmPtVel.y;
             ucat[k][j][i].z = (sb/sc) * bPtVel.z + (1.0 - (sb/sc)) * ibmPtVel.z;
         }


         if (flags->isTeqnActive)
         {
            if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "zeroGradient" || ibm->ibmBody[ibF[c].bodyID]->tempBC == "thetaWallFunction")
            {
                //apply zerogradient boundary condition. effect of boundary through heat flux bc
                temp[k][j][i] =  bPtTemp;
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "fixedValue")
            {
                temp[k][j][i] =  ibm->ibmBody[ibF[c].bodyID]->fixedTemp;
            }
         }

         // save the element acceleration term (dudt . elementNormal)
         PetscReal bP = nDot(nScale(-1.0/clock->dtOld, nSub(ibmPtVel, ibmPtVelPrev)), eNorm);

         p[k][j][i] = bP * (sb - sc) + bPtPres;
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da, peqn->lP, &lp);
    DMDAVecRestoreArray(da, peqn->P, &p);
    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
        DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    }

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd(da, peqn->P, INSERT_VALUES, peqn->lP);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CurvibInterpolationQuadratic(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    peqn_         *peqn  = ibm->access->peqn;
    flags_        *flags = ibm->access->flags;
    constants_    *cst   = ibm->access->constants;
    clock_        *clock = ibm->access->clock;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys   = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs   = info.gzs, gze = info.gzs + info.gzm;

    PetscInt      i, j, k, ii, kk, jj;
    PetscInt      b, c, s;
    cellIds       initCp;
    Cmpnts        ***ucat, ***lucat, ***cent, ***csi, ***eta, ***zet;
    PetscReal     ***lt, ***temp, ***aj, ***nvert, ***lp, ***p;
    Cmpnts        bPt1, bPt2, bPtInit;
    Cmpnts        bPt1Vel, bPt2Vel, ibmPtVel, ibmPtVelPrev;
    PetscReal     nfMag, cellSize, bPt1Temp, bPt2Temp, ibmPtTemp, bPt1Pres, bPt2Pres, ibmPtPres, elemDist;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, peqn->lP, &lp);
    DMDAVecGetArray(da, peqn->P, &p);
    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
        DMDAVecGetArray(da, teqn->Tmprt, &temp);
    }

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    PetscReal    sc, sb, sd, ustar;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmMesh   *ibMsh = ibm->ibmBody[ibF[c].bodyID]->ibMsh;
        PetscInt   cElem = ibF[c].closestElem;
        Cmpnts	   eNorm = ibF[c].normal;
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal  roughness;

        /********************************************************************************************************
            Interpolation of the first background point field
        /*******************************************************************************************************/

        // distance from cell to the IBM element
        elemDist = fabs(nDot(nSub(ibF[c].pMin, cent[k][j][i]), eNorm));

        // background mesh projection point
        cellSize = 1.0 * elemDist;
        bPt1 = nScale(cellSize, eNorm);
        mSum(bPt1, cent[k][j][i]);

        // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
        PetscInt    i1, j1, k1;
        PetscReal   dmin = 10e10, d;
        PetscInt    ic, jc, kc, setFlag = 0;

        // must be close so don't loop over all cells
        for (k1=k-1; k1<k+2; k1++)
        for (j1=j-1; j1<j+2; j1++)
        for (i1=i-1; i1<i+2; i1++)
        {

            if
            (
                (
                    k1>=1 && k1<mz-1 &&
                    j1>=1 && j1<my-1 &&
                    i1>=1 && i1<mx-1
                ) && (isFluidCell(k1, j1, i1, nvert))
            )
            {
                d = pow((bPt1.x - cent[k1][j1][i1].x), 2) +
                    pow((bPt1.y - cent[k1][j1][i1].y), 2) +
                    pow((bPt1.z - cent[k1][j1][i1].z), 2);

                if
                (
                    d < dmin
                )
                {
                    dmin  = d;
                    ic = i1;
                    jc = j1;
                    kc = k1;
                    setFlag = 1;
                }
            }
        }

        if(setFlag == 0)
        {
            char error[512];
            sprintf(error, "no closest fluid cell to ibm fluid %ld %ld %ld\n", k, j, i);
            fatalErrorInFunction("CurvibInterpolation",  error);
        }

        PetscInt intId[6];
        PetscInt intFlag = 1;

        // get the trilinear interpolation cells
        PointInterpolationCells
        (
                mesh,
                bPt1.x, bPt1.y, bPt1.z,
                ic, jc, kc,
                cent,
                intId
        );

        // save the initial closest cell and background point
        initCp.i = ic; initCp.j = jc; initCp.k = kc;
        bPtInit  = nSet(bPt1);

        // ensure none of the interpolation cells are solid cells.
        PetscReal sumDel = 0.0;

        while ( (intFlag == 1) && isInsideBoundingBox(bPt1, mesh->bounds) && (sumDel <= 1.5*cellSize))
        {
            PetscInt ibmCellCtr = 0;
            PetscInt icc, jcc, kcc, setFlag = 0;

            for (PetscInt kk = 0; kk<2; kk++)
            for (PetscInt jj = 0; jj<2; jj++)
            for (PetscInt ii = 0; ii<2; ii++)
            {
                if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                {
                    ibmCellCtr ++;
                }

            }

            if (ibmCellCtr > 0)
            {
                Cmpnts del =  nScale(0.2 * cellSize, eNorm);
                mSum(bPt1, del);
                dmin = 10e10;

                sumDel += nMag(del);

                for (k1=kc-1; k1<kc+2; k1++)
                {
                    //check processor ghost bounds
                    if (k1 < gzs || k1 >= gze) {intFlag = 2; break;}
                    for (j1=jc-1; j1<jc+2; j1++)
                    {
                        if (j1 < gys || j1 >= gye) {intFlag = 2; break;}
                        for (i1=ic-1; i1<ic+2; i1++)
                        {
                            if (i1 < gxs || i1 >= gxe) {intFlag = 2; break;}

                            if
                            (
                                (
                                    k1>=1 && k1<mz-1 &&
                                    j1>=1 && j1<my-1 &&
                                    i1>=1 && i1<mx-1
                                ) && (isFluidCell(k1, j1, i1, nvert))
                            )
                            {
                                d = pow((bPt1.x - cent[k1][j1][i1].x), 2) +
                                    pow((bPt1.y - cent[k1][j1][i1].y), 2) +
                                    pow((bPt1.z - cent[k1][j1][i1].z), 2);

                                if
                                (
                                    d < dmin
                                )
                                {
                                    dmin  = d;
                                    icc = i1;
                                    jcc = j1;
                                    kcc = k1;
                                    setFlag = 1;
                                }
                            }

                        }
                    }
                }

                if(setFlag == 0)
                {
                    //closest point not set, do not interpolate
                    intFlag = 2;
                }

                if(intFlag == 1)
                {
                    kc = kcc; jc = jcc; ic = icc;

                    PointInterpolationCells
                    (
                            mesh,
                            bPt1.x, bPt1.y, bPt1.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );
                }

            }
            else
            {
                intFlag = 0;

                // trilinear interpolate the velocity at this point
                vectorPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt1.x, bPt1.y, bPt1.z,
                        ic, jc, kc,
                        cent,
                        lucat,
                        bPt1Vel
                );

                if (flags->isTeqnActive)
                {
                    scalarPointLocalVolumeInterpolation
                    (
                            mesh,
                            bPt1.x, bPt1.y, bPt1.z,
                            ic, jc, kc,
                            cent,
                            lt,
                            bPt1Temp
                    );
                }

                scalarPointLocalVolumeInterpolation
                (
                        mesh,
                        bPt1.x, bPt1.y, bPt1.z,
                        ic, jc, kc,
                        cent,
                        lp,
                        bPt1Pres
                );

            }
        }

        // while loop fails
        if(intFlag > 0)
        {
            bPt1Vel = nSet(lucat[initCp.k][initCp.j][initCp.i]);

            bPt1Pres = lp[initCp.k][initCp.j][initCp.i];

            if (flags->isTeqnActive)
            {
                bPt1Temp = lt[initCp.k][initCp.j][initCp.i];
            }

            bPt1 = nSet(bPtInit);

            ic = initCp.i;
            jc = initCp.j;
            kc = initCp.k;
        }

         /********************************************************************************************************
             Interpolation of the second background point field
         /*******************************************************************************************************/

         bPt2 = nScale(cellSize, eNorm);
         mSum(bPt2, bPt1);

         // search for the closest background mesh cell center (it will be close to current closest fluid cell kc,jc,ic)
         dmin = 10e10; setFlag = 0;
         PetscInt    ic2, jc2, kc2;

         // must be close so don't loop over all cells
         for (k1=kc-1; k1<kc+2; k1++)
         for (j1=jc-1; j1<jc+2; j1++)
         for (i1=ic-1; i1<ic+2; i1++)
         {

             if
             (
                 (
                     k1>=1 && k1<mz-1 &&
                     j1>=1 && j1<my-1 &&
                     i1>=1 && i1<mx-1
                 ) && (isFluidCell(k1, j1, i1, nvert))
             )
             {
                 d = pow((bPt2.x - cent[k1][j1][i1].x), 2) +
                     pow((bPt2.y - cent[k1][j1][i1].y), 2) +
                     pow((bPt2.z - cent[k1][j1][i1].z), 2);

                 if
                 (
                     d < dmin
                 )
                 {
                     dmin  = d;
                     ic2 = i1;
                     jc2 = j1;
                     kc2 = k1;
                     setFlag = 1;
                 }
             }
         }

         if(setFlag == 0)
         {
             char error[512];
             sprintf(error, "no closest fluid cell to background point 2\n");
             fatalErrorInFunction("CurvibInterpolation",  error);
         }

         //check for the processor bound
         intFlag = 1;

         if (kc2 < gzs || kc2 >= gze) {intFlag = 2;}
         if (jc2 < gys || jc2 >= gye) {intFlag = 2;}
         if (ic2 < gxs || ic2 >= gxe) {intFlag = 2;}

         PetscInt ibmCellCtr = 0;

         if(intFlag == 1)
         {
             // get the trilinear interpolation cells
             PointInterpolationCells
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     intId
             );

             //check that the interpolation cells are not ibm solid
             for (PetscInt kk = 0; kk<2; kk++)
             for (PetscInt jj = 0; jj<2; jj++)
             for (PetscInt ii = 0; ii<2; ii++)
             {
                 if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                 {
                     ibmCellCtr ++;
                 }

             }
         }

         if(intFlag == 2)
         {
             bPt2Vel = nSet(bPt1Vel);

             bPt2Pres =  bPt1Pres;
             if (flags->isTeqnActive)
             {
                 bPt2Temp = bPt1Temp;
             }
         }
         else if(intFlag == 1 && ibmCellCtr > 0)
         {
             bPt2Vel = nSet(lucat[kc2][jc2][ic2]);
             
             bPt2Pres = lp[kc2][jc2][ic2];
             
             if (flags->isTeqnActive)
             {
                 bPt2Temp = lt[kc2][jc2][ic2];
             }
         }
         else
         {
             // trilinear interpolate the velocity at this point
             vectorPointLocalVolumeInterpolation
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     lucat,
                     bPt2Vel
             );

             scalarPointLocalVolumeInterpolation
             (
                     mesh,
                     bPt2.x, bPt2.y, bPt2.z,
                     ic2, jc2, kc2,
                     cent,
                     lp,
                     bPt2Pres
             );

             if (flags->isTeqnActive)
             {
                 scalarPointLocalVolumeInterpolation
                 (
                         mesh,
                         bPt2.x, bPt2.y, bPt2.z,
                         ic2, jc2, kc2,
                         cent,
                         lt,
                         bPt2Temp
                 );
             }
         }

         /********************************************************************************************************
             Interpolation of the IBM surface field
         /*******************************************************************************************************/

        // interpolate the velocity of the projected point on the IBM solid element from its nodes
        ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                  + ibMsh->nU[n2].x * ibF[c].cs2
                  + ibMsh->nU[n3].x * ibF[c].cs3;

        ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                  + ibMsh->nU[n2].y * ibF[c].cs2
                  + ibMsh->nU[n3].y * ibF[c].cs3;

        ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                  + ibMsh->nU[n2].z * ibF[c].cs2
                  + ibMsh->nU[n3].z * ibF[c].cs3;

         // interpolate the previous velocity
         ibmPtVelPrev.x =   ibMsh->nUPrev[n1].x * ibF[c].cs1
                          + ibMsh->nUPrev[n2].x * ibF[c].cs2
                          + ibMsh->nUPrev[n3].x * ibF[c].cs3;

         ibmPtVelPrev.y =   ibMsh->nUPrev[n1].y * ibF[c].cs1
                          + ibMsh->nUPrev[n2].y * ibF[c].cs2
                          + ibMsh->nUPrev[n3].y * ibF[c].cs3;

         ibmPtVelPrev.z =   ibMsh->nUPrev[n1].z * ibF[c].cs1
                          + ibMsh->nUPrev[n2].z * ibF[c].cs2
                          + ibMsh->nUPrev[n3].z * ibF[c].cs3;

        /********************************************************************************************************
            Interpolation of the IBM fluid cell field
        /*******************************************************************************************************/

         //find the distance from the wall of points b and c along the element normal
         sb = elemDist;
         sc = sb + nMag(nSub(bPt1, cent[k][j][i]));
         sd = sc + nMag(nSub(bPt2, bPt1));

        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "noSlip")
        {
            //find quadratic interpolation coefficients a and b
            PetscReal a1, a2, a3, b1, b2, b3, denom;
            Cmpnts  aVel, bVel;
            denom = sd*sc*sc - sc*sd*sd;
            a1 = sd/denom;
            a2 = -sc/denom;
            a3 = -cellSize/denom;

            b1 = -sd*sd/denom;
            b2 = sc*sc/denom;
            b3 = (sd*sd - sd*sc + cellSize*sc)/denom;

            aVel  = nSum(nScale(a1, bPt1Vel), nScale(a2, bPt2Vel));
            mSum(aVel, nScale(a3, ibmPtVel));

            bVel  = nSum(nScale(b1, bPt1Vel), nScale(b2, bPt2Vel));
            mSum(bVel, nScale(b3, ibmPtVel));

            if(intFlag == 2 || (intFlag == 1 && ibmCellCtr > 0))
            {
                ucat[k][j][i].x = (sb/sc) * bPt1Vel.x + (1.0 - (sb/sc)) * ibmPtVel.x;
                ucat[k][j][i].y= (sb/sc) * bPt1Vel.y + (1.0 - (sb/sc)) * ibmPtVel.y;
                ucat[k][j][i].z = (sb/sc) * bPt1Vel.z + (1.0 - (sb/sc)) * ibmPtVel.z;
            }
            else
            {
                ucat[k][j][i].x =  aVel.x * sb * sb + bVel.x * sb + ibmPtVel.x;
                ucat[k][j][i].y =  aVel.y * sb * sb + bVel.y * sb + ibmPtVel.y;
                ucat[k][j][i].z =  aVel.z * sb * sb + bVel.z * sb + ibmPtVel.z;
            }

        }
        else if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
        {
            if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
            {
                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;

                wallFunctionSchumann(cst->nu, sc, sb, wm->roughness, wm->kappa, ibmPtVel,
                        bPt1Vel, &ucat[k][j][i], &ustar, eNorm);
            }
            else if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
            {
                Cabot *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot;

                if(roughness > 1.0e-12)
                {
                    wallFunctionCabotRoughness(cst->nu, wm->roughness, sc, sb, ibmPtVel, bPt1Vel, &ucat[k][j][i], &ustar, eNorm);
                }
                else
                {
                    wallFunctionCabot(cst->nu, sc, sb, ibmPtVel, bPt1Vel, &ucat[k][j][i], &ustar, eNorm);
                }
            }
        }


         //pressure and temperature boundary condition - neumann
         PetscReal aTemp, bTemp, cTemp, aP, bP, cP;

         bP = nDot(nScale(-1.0/clock->dtOld, nSub(ibmPtVel, ibmPtVelPrev)), eNorm);
         aP = (bPt2Pres - bPt1Pres - bP * cellSize)/(sd*sd - sc*sc);
         cP = bPt1Pres - aP*sc*sc - bP*sc;

         if(intFlag == 2 || (intFlag == 1 && ibmCellCtr > 0))
         {
             p[k][j][i] = bP * (sb - sc) + bPt1Pres;
         }
         else
         {
             p[k][j][i] = aP * sb * sb + bP * sb + cP;
         }

         if (flags->isTeqnActive)
         {
            if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "zeroGradient" || ibm->ibmBody[ibF[c].bodyID]->tempBC == "thetaWallFunction")
            {
                //apply zerogradient boundary condition. effect of boundary through heat flux bc
                temp[k][j][i] =  bPt1Temp;
            }
            else if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "fixedValue")
            {
                temp[k][j][i] =  ibm->ibmBody[ibF[c].bodyID]->fixedTemp;
            }
            // interpolated from the two background cells
            // else 
            // {
            //     bTemp = 0.0;
            //     aTemp = (bPt2Temp - bPt1Temp)/(sd*sd - sc*sc);
            //     cTemp = bPt1Temp - aTemp*sc*sc;

            //     if(intFlag == 2 || (intFlag == 1 && ibmCellCtr > 0))
            //     {
            //         temp[k][j][i] =  bPt1Temp;
            //     }
            //     else
            //     {
            //         temp[k][j][i] = aTemp * sb * sb + cTemp;
            //     }
            // }
         }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, peqn->lP, &lp);
    DMDAVecRestoreArray(da, peqn->P, &p);
    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
        DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    }

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd(da, peqn->P, INSERT_VALUES, peqn->lP);
    return(0);
}

//***************************************************************************************************************//
PetscErrorCode findIBMWallShearChester(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    les_          *les   = ibm->access->les;
    clock_        *clock = ibm->access->clock;
    constants_    *cst   = ibm->access->constants;
    flags_        *flags = ibm->access->flags;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys   = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs   = info.gzs, gze = info.gzs + info.gzm;

    PetscReal     du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

    PetscReal     dut1dn, dut2dn, dundn, dut1dt1, dut2dt1, dundt1, dut1dt2, dut2dt2, dundt2;

    PetscReal     dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz;
    PetscReal     dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords

    PetscReal     nu = cst->nu, nut, uTmag, ustar, cellSize, elemDist;
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc;      // surface area vectors components

    PetscReal     g11, g21, g31;                                             // metric tensor components
    PetscReal     r11, r21, r31,
                  r12, r22, r32,
                  r13, r23, r33;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i, j, k, c;
    PetscInt      i1, j1, k1;

    cellIds       initCp;
    Vec           Coor;
    Cmpnts        ***coor, ***ucat, ***cent, ibmPtVel;
    Cmpnts        ***visc1, ***visc2, ***visc3, ***viscT;

    Cmpnts	      ***icsi, ***ieta, ***izet;
	Cmpnts	      ***jcsi, ***jeta, ***jzet;
	Cmpnts	      ***kcsi, ***keta, ***kzet;

    Cmpnts        uC, uN, uT, uB;                          // background point velocity, its normal and tangential components
    Cmpnts        bPtVel, bPt, bPtInit;
    Cmpnts        eN, eT1, eT2;                          // local wall normal co-ordinate system
    PetscReal     ***nvert, ***iaj, ***jaj, ***kaj, ***aj, ***lt, ***lnu_t;
    PetscReal     bPtTemp;

    //local pointer for ibmFluidCells
    ibmFluidCell  *ibF = ibm->ibmFCells;

    //rotation tensor direction cosines
    PetscReal       a11, a21, a31,
                    a12, a22, a32,
                    a13, a23, a33;

    PetscReal       tau11, tau21, tau31,
                    tau12, tau22, tau32,
                    tau13, tau23, tau33;

    VecSet(ueqn->lViscIBM1,  0.0);
    VecSet(ueqn->lViscIBM2,  0.0);
    VecSet(ueqn->lViscIBM3,  0.0);

    if(flags->isTeqnActive)
    {
        VecSet(teqn->lViscIBMT,  0.0);
    }

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, ueqn->lViscIBM1, &visc1);
    DMDAVecGetArray(fda, ueqn->lViscIBM2, &visc2);
    DMDAVecGetArray(fda, ueqn->lViscIBM3, &visc3);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(fda, teqn->lViscIBMT, &viscT);
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
    }

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
	DMDAVecGetArray(fda, mesh->lIEta, &ieta);
	DMDAVecGetArray(fda, mesh->lIZet, &izet);

	DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
	DMDAVecGetArray(fda, mesh->lJEta, &jeta);
	DMDAVecGetArray(fda, mesh->lJZet, &jzet);

	DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
	DMDAVecGetArray(fda, mesh->lKEta, &keta);
	DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
	DMDAVecGetArray(da, mesh->lJAj, &jaj);
	DMDAVecGetArray(da, mesh->lKAj, &kaj);

    if(ibm->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
    }

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmObject     *ibmBody = ibm->ibmBody[ibF[c].bodyID];
        ibmMesh       *ibMsh = ibmBody->ibMsh;
        PetscInt      cElem = ibF[c].closestElem;
        Cmpnts	      eNorm = ibMsh->eN[cElem];
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal     tauWall, sb, qWall = 0;
        Cmpnts        qFlux;

        //find the distance from the wall where velocity is interpolated
        //a uniform distance from the wall is chosen for all ibm fluid points to prevent instabilities due to the stair step structure
        cellSize = ibm->interpDist;

        // distance from cell to the IBM element
        elemDist = fabs(nDot(nSub(ibF[c].pMin, cent[k][j][i]), eNorm));

        bPt = nScale(cellSize + elemDist, eNorm);
        mSum(bPt, cent[k][j][i]);

        //check the neighbourhood of current ibm fluid cell for fluid cells
        //find the fluid - ibm fluid cell interface faces
        for (PetscInt k11=k-1; k11<k+2; k11++)
        for (PetscInt j11=j-1; j11<j+2; j11++)
        for (PetscInt i11=i-1; i11<i+2; i11++)
        {
            PetscInt      isFace = 0;

            if
            (
                (
                    k11>=1 && k11<mz-1 &&
                    j11>=1 && j11<my-1 &&
                    i11>=1 && i11<mx-1
                ) && isFluidCell(k11, j11, i11, nvert)
            )
            {
                //define the faces
                // 1 - k:k+1
                // 2 - k-1:k
                // 3 - j:j+1
                // 4 - j-1:j
                // 5 - i:i+1
                // 6 - i-1:i

                if(k11 == k+1 && j11 == j && i11 == i)
                {
                    isFace = 1;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k-1 && j11 == j && i11 == i)
                {
                    isFace = 2;
                    k1 = k-1; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j+1 && i11 == i)
                {
                    isFace = 3;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j-1 && i11 == i)
                {
                    isFace = 4;
                    k1 = k; j1 = j-1; i1 = i;
                }
                else if (k11 == k && j11 == j && i11 == i+1)
                {
                    isFace = 5;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j && i11 == i-1)
                {
                    isFace = 6;
                    k1 = k; j1 = j; i1 = i-1;
                }

                // only for the cells that share a face, not diagonal cells
                if(isFace)
                {

                    // interpolate the velocity of the projected point on the IBM solid element from its nodes
                    ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                                 + ibMsh->nU[n2].x * ibF[c].cs2
                                 + ibMsh->nU[n3].x * ibF[c].cs3;

                    ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                                 + ibMsh->nU[n2].y * ibF[c].cs2
                                 + ibMsh->nU[n3].y * ibF[c].cs3;

                    ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                                 + ibMsh->nU[n2].z * ibF[c].cs2
                                 + ibMsh->nU[n3].z * ibF[c].cs3;

                    // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
                    PetscInt    ia, ja, ka;
                    PetscReal   dmin = 10e10, d;
                    PetscInt    ic, jc, kc, setFlag = 0;

                    // must be close so don't loop over all cells
                    for (ka=k-1; ka<k+2; ka++)
                    for (ja=j-1; ja<j+2; ja++)
                    for (ia=i-1; ia<i+2; ia++)
                    {

                        if
                        (
                            (
                                ka>=1 && ka<mz-1 &&
                                ja>=1 && ja<my-1 &&
                                ia>=1 && ia<mx-1
                            ) && isFluidCell(ka, ja, ia, nvert)
                        )
                        {
                            d = pow((bPt.x - cent[ka][ja][ia].x), 2) +
                                pow((bPt.y - cent[ka][ja][ia].y), 2) +
                                pow((bPt.z - cent[ka][ja][ia].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                ic = ia;
                                jc = ja;
                                kc = ka;
                                setFlag = 1;
                            }
                        }
                    }

                    if(setFlag == 0)
                    {
                        char error[512];
                        sprintf(error, "no closest fluid cell to ibm fluid %ld %ld %ld\n", k, j, i);
                        fatalErrorInFunction("findIBMWallShearChester",  error);
                    }  

                    PetscInt intId[6];
                    PetscInt intFlag = 1;

                    // get the trilinear interpolation cells
                    PointInterpolationCells
                    (
                            mesh,
                            bPt.x, bPt.y, bPt.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );

                    // save the initial closest cell and background point
                    initCp.i = ic; initCp.j = jc; initCp.k = kc;
                    bPtInit  = nSet(bPt);                    

                    PetscReal sumDel = 0.0;

                    while ( (intFlag == 1) && isInsideBoundingBox(bPt, mesh->bounds) && (sumDel <= 1.5*cellSize))
                    {
                        PetscInt ibmCellCtr = 0;
                        PetscInt icc, jcc, kcc, setFlag = 0;

                        for (PetscInt kk = 0; kk<2; kk++)
                        for (PetscInt jj = 0; jj<2; jj++)
                        for (PetscInt ii = 0; ii<2; ii++)
                        {
                            if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                            {
                                ibmCellCtr ++;
                            }

                        }

                        if (ibmCellCtr > 0)
                        {
                            Cmpnts del =  nScale(0.2 * cellSize, eNorm);
                            mSum(bPt, del);
                            dmin = 10e10;

                            sumDel += nMag(del);

                            for (ka=kc-1; ka<kc+2; ka++)
                            {
                                //check processor ghost bounds
                                if (ka < gzs || ka >= gze) {intFlag = 2; break;}
                                for (ja=jc-1; ja<jc+2; ja++)
                                {
                                    if (ja < gys || ja >= gye) {intFlag = 2; break;}
                                    for (ia=ic-1; ia<ic+2; ia++)
                                    {
                                        if (ia < gxs || ia >= gxe) {intFlag = 2; break;}

                                        if
                                        (
                                            (
                                                ka>=1 && ka<mz-1 &&
                                                ja>=1 && ja<my-1 &&
                                                ia>=1 && ia<mx-1
                                            ) && (isFluidCell(ka, ja, ia, nvert))
                                        )
                                        {
                                            d = pow((bPt.x - cent[ka][ja][ia].x), 2) +
                                                pow((bPt.y - cent[ka][ja][ia].y), 2) +
                                                pow((bPt.z - cent[ka][ja][ia].z), 2);

                                            if
                                            (
                                                d < dmin
                                            )
                                            {
                                                dmin  = d;
                                                icc = ia;
                                                jcc = ja;
                                                kcc = ka;
                                                setFlag = 1;
                                            }
                                        }

                                    }
                                }
                            }

                            if(setFlag == 0)
                            {
                                //closest point not set, do not interpolate
                                intFlag = 2;
                            }

                            if(intFlag == 1)
                            {
                                kc = kcc; jc = jcc; ic = icc;

                                PointInterpolationCells
                                (
                                        mesh,
                                        bPt.x, bPt.y, bPt.z,
                                        ic, jc, kc,
                                        cent,
                                        intId
                                );
                            }

                        }
                        else
                        {
                            intFlag = 0;

                            // trilinear interpolate the velocity at this point
                            vectorPointLocalVolumeInterpolation
                            (
                                    mesh,
                                    bPt.x, bPt.y, bPt.z,
                                    ic, jc, kc,
                                    cent,
                                    ucat,
                                    bPtVel
                            );

                            if(flags->isTeqnActive)
                            {
                                scalarPointLocalVolumeInterpolation
                                (
                                        mesh,
                                        bPt.x, bPt.y, bPt.z,
                                        ic, jc, kc,
                                        cent,
                                        lt,
                                        bPtTemp
                                );
                            }
                        }
                    }

                    // while loop fails
                    if(intFlag > 0)
                    {
                        bPtVel = nSet(ucat[initCp.k][initCp.j][initCp.i]);

                        if(flags->isTeqnActive)
                        {
                            bPtTemp = lt[initCp.k][initCp.j][initCp.i];
                        }

                        bPt = nSet(bPtInit);

                        ic = initCp.i;
                        jc = initCp.j;
                        kc = initCp.k;
                    }
                    
                    // distance of background point from the wall
                    sb = nDot(nSub(bPt, ibF[c].pMin), eNorm);

                    //velocity in local co-ordinate system
                    mSub(bPtVel, ibmPtVel);
                    uN = nScale(nDot(bPtVel, eNorm), eNorm);
                    uT = nSub(bPtVel, uN);
                    uTmag = nMag(uT);

                    if (flags->isTeqnActive)
                    {
                        if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "thetaWallFunction")
                        {
                            if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -3)
                            {

                            }
                            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -2)
                            {
                                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelT->wmShumann;

                                qWall  = wm->qWall;
                            }
                            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -4)
                            {

                                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelT->wmShumann;
                                
                                // find interpolation weights for surface temp
                                PetscReal w[2];
                                PetscInt  l[2];

                                findInterpolationWeights(w, l, wm->timeVec, wm->numT, clock->time);

                                //interpolate the surface temperature from closest available mesoscale times
                                PetscReal surfaceTemp = w[0] * wm->surfTemp[l[0]] + w[1] * wm->surfTemp[l[1]];
                                PetscReal surfaceL;

                                PetscReal deltaTheta = bPtTemp - surfaceTemp;
                                PetscReal phiM, phiH;
                                
                                qWallShumann
                                (
                                    uTmag, sb, wm->roughness,
                                    wm->gammaM, wm->gammaH, wm->alphaH,
                                    wm->thetaRef, deltaTheta, wm->kappa,
                                    qWall, ustar, phiM, phiH, surfaceL
                                );

                                // if(i  == 5 && k == 5 && j == 5)
                                //     PetscPrintf(PETSC_COMM_SELF, "surfTemp = %lf, L = %lf, bPtTemp = %lf, deltaTheta = %lf, qWall = %lf, ustar = %lf, utmag = %lf, walldist = %lf, wts = %lf %lf, ind = %ld\n", surfaceTemp, surfaceL, bPtTemp, deltaTheta, qWall, ustar, uTmag, sb, w[0], w[1], l[0]);
                            }
                        }
                    }

                    if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
                    {
                        if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
                        {
                            PetscReal roughness = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot->roughness;

                            if(roughness > 1.0e-12)
                            {
                                ustar = uTauCabotRoughness(nu, uTmag, sb, 0.01, 0, roughness);

                            }
                            else
                            {
                                ustar = uTauCabot(nu, uTmag, sb, 0.01, 0);
                            }
                        }
                        else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
                        {
                            Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;
                            PetscReal phiM, L;

                            uStarShumann
                            (
                                uTmag, sb, wm->roughness,
                                wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                                ustar, phiM, L
                            );

                            // if(i  == 5 && k == 5 && j == 5)
                            //     PetscPrintf(PETSC_COMM_SELF, "qWall = %lf, ustar = %lf, utmag = %lf, walldist = %lf\n", qWall, ustar, uTmag, sb);

                        }
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall shear model is on but velocity BC is not set to velocityWallFunction\n");
                        fatalErrorInFunction("findIBMWallShearChester", error);
                    }

                    //local co-ordinate system
                    eT1 = nUnit(uT);
                    eN = eNorm;
                    eT2 = nCross(eN, eT1);

                    tauWall = ustar * ustar;

                    //create the transformation vector for rotation to global axis.
                    a11 = eT1.x; a12 = eT2.x, a13 = eN.x;
                    a21 = eT1.y; a22 = eT2.y, a23 = eN.y;
                    a31 = eT1.z; a32 = eT2.z, a33 = eN.z;

                    //transform it to original co-ordinate system
                    tau11 = 2.0 * a11 * a13 * tauWall;
                    tau12 = (a13 * a21 + a11 * a23) * tauWall;
                    tau13 = (a13 * a31 + a11 * a33) * tauWall;

                    tau21 = (a23 * a11 + a21 * a13) * tauWall;
                    tau22 = 2.0 * a23 * a21 * tauWall;
                    tau23 = (a23 * a31 + a21 * a33) * tauWall;

                    tau31 = (a33 * a11 + a31 * a13) * tauWall;
                    tau32 = (a33 * a21 + a31 * a23) * tauWall;
                    tau33 = 2.0 * a33 * a31 * tauWall;

                    //save the wall shear force
                    if(isFace == 1 || isFace == 2)
                    {
                        visc3[k][j][i].x = tau11* kzet[k1][j1][i1].x + tau12 * kzet[k1][j1][i1].y + tau13 * kzet[k1][j1][i1].z;
                        visc3[k][j][i].y = tau21* kzet[k1][j1][i1].x + tau22 * kzet[k1][j1][i1].y + tau23 * kzet[k1][j1][i1].z;
                        visc3[k][j][i].z = tau31* kzet[k1][j1][i1].x + tau32 * kzet[k1][j1][i1].y + tau33 * kzet[k1][j1][i1].z;

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].z = qFlux.x * kzet[k1][j1][i1].x + qFlux.y * kzet[k1][j1][i1].y + qFlux.z * kzet[k1][j1][i1].z;
                        }
                    }
                    else if(isFace == 3 || isFace == 4)
                    {
                        visc2[k][j][i].x = tau11* jeta[k1][j1][i1].x + tau12 * jeta[k1][j1][i1].y + tau13 * jeta[k1][j1][i1].z;
                        visc2[k][j][i].y = tau21* jeta[k1][j1][i1].x + tau22 * jeta[k1][j1][i1].y + tau23 * jeta[k1][j1][i1].z;
                        visc2[k][j][i].z = tau31* jeta[k1][j1][i1].x + tau32 * jeta[k1][j1][i1].y + tau33 * jeta[k1][j1][i1].z;

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].y = qFlux.x * jeta[k1][j1][i1].x + qFlux.y * jeta[k1][j1][i1].y + qFlux.z * jeta[k1][j1][i1].z;
                        }
                    }
                    else if(isFace == 5 || isFace == 6)
                    {
                        visc1[k][j][i].x = tau11* icsi[k1][j1][i1].x + tau12 * icsi[k1][j1][i1].y + tau13 * icsi[k1][j1][i1].z;
                        visc1[k][j][i].y = tau21* icsi[k1][j1][i1].x + tau22 * icsi[k1][j1][i1].y + tau23 * icsi[k1][j1][i1].z;
                        visc1[k][j][i].z = tau31* icsi[k1][j1][i1].x + tau32 * icsi[k1][j1][i1].y + tau33 * icsi[k1][j1][i1].z;

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].x = qFlux.x * icsi[k1][j1][i1].x + qFlux.y * icsi[k1][j1][i1].y + qFlux.z * icsi[k1][j1][i1].z;
                        }
                    }

                }

            }
        }

    }

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    DMDAVecRestoreArray(fda, ueqn->lViscIBM1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lViscIBM2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lViscIBM3, &visc3);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
	DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
	DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

	DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
	DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
	DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

	DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
	DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
	DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
	DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
	DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    if(ibm->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(fda, teqn->lViscIBMT, &viscT);
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMLocalToLocalBegin(fda, teqn->lViscIBMT, INSERT_VALUES, teqn->lViscIBMT);
        DMLocalToLocalEnd  (fda, teqn->lViscIBMT, INSERT_VALUES, teqn->lViscIBMT);
    }

    DMLocalToLocalBegin(fda, ueqn->lViscIBM1, INSERT_VALUES, ueqn->lViscIBM1);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM1, INSERT_VALUES, ueqn->lViscIBM1);
    DMLocalToLocalBegin(fda, ueqn->lViscIBM2, INSERT_VALUES, ueqn->lViscIBM2);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM2, INSERT_VALUES, ueqn->lViscIBM2);
    DMLocalToLocalBegin(fda, ueqn->lViscIBM3, INSERT_VALUES, ueqn->lViscIBM3);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM3, INSERT_VALUES, ueqn->lViscIBM3);

    return (0);
}

//***************************************************************************************************************//
PetscErrorCode findIBMWallShear(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    les_          *les   = ibm->access->les;
    clock_        *clock = ibm->access->clock;
    constants_    *cst   = ibm->access->constants;
    flags_        *flags = ibm->access->flags;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys   = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs   = info.gzs, gze = info.gzs + info.gzm;

    PetscReal     du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

    PetscReal     dut1dn, dut2dn, dundn, dut1dt1, dut2dt1, dundt1, dut1dt2, dut2dt2, dundt2;

    PetscReal     dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz;
    PetscReal     dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords

    PetscReal     nu = cst->nu, nut, uTmag, ustar, cellSize, elemDist;
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc;      // surface area vectors components

    PetscReal     g11, g21, g31;                                             // metric tensor components
    PetscReal     r11, r21, r31,
                  r12, r22, r32,
                  r13, r23, r33;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i, j, k, c;
    PetscInt      i1, j1, k1;

    cellIds       initCp;
    Vec           Coor;
    Cmpnts        ***coor, ***ucat, ***cent, ibmPtVel;
    Cmpnts        ***visc1, ***visc2, ***visc3, ***viscT;

    Cmpnts	      ***icsi, ***ieta, ***izet;
	Cmpnts	      ***jcsi, ***jeta, ***jzet;
	Cmpnts	      ***kcsi, ***keta, ***kzet;

    Cmpnts        uC, uN, uT, uB;                          // background point velocity, its normal and tangential components
    Cmpnts        bPtVel, bPt, bPtInit;
    Cmpnts        eN, eT1, eT2;                          // local wall normal co-ordinate system
    PetscReal     ***nvert, ***meshTag, ***iaj, ***jaj, ***kaj, ***aj, ***lt, ***lnu_t;
    PetscReal     bPtTemp;

    //local pointer for ibmFluidCells
    ibmFluidCell  *ibF = ibm->ibmFCells;

    //rotation tensor direction cosines
    PetscReal       a11, a21, a31,
                    a12, a22, a32,
                    a13, a23, a33;

    PetscReal       tau11, tau21, tau31,
                    tau12, tau22, tau32,
                    tau13, tau23, tau33;

    VecSet(ueqn->lViscIBM1,  0.0);
    VecSet(ueqn->lViscIBM2,  0.0);
    VecSet(ueqn->lViscIBM3,  0.0);

    if(flags->isTeqnActive)
    {
        VecSet(teqn->lViscIBMT,  0.0);
    }

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(fda, ueqn->lViscIBM1, &visc1);
    DMDAVecGetArray(fda, ueqn->lViscIBM2, &visc2);
    DMDAVecGetArray(fda, ueqn->lViscIBM3, &visc3);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(fda, teqn->lViscIBMT, &viscT);
        DMDAVecGetArray(da, teqn->lTmprt, &lt);
    }

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
	DMDAVecGetArray(fda, mesh->lIEta, &ieta);
	DMDAVecGetArray(fda, mesh->lIZet, &izet);

	DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
	DMDAVecGetArray(fda, mesh->lJEta, &jeta);
	DMDAVecGetArray(fda, mesh->lJZet, &jzet);

	DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
	DMDAVecGetArray(fda, mesh->lKEta, &keta);
	DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
	DMDAVecGetArray(da, mesh->lJAj, &jaj);
	DMDAVecGetArray(da, mesh->lKAj, &kaj);

    if(ibm->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
    }

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        ibmObject     *ibmBody = ibm->ibmBody[ibF[c].bodyID];
        ibmMesh       *ibMsh = ibmBody->ibMsh;
        PetscInt      cElem = ibF[c].closestElem;
        Cmpnts	      eNorm = ibMsh->eN[cElem];
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal     tauWall, sb, qWall = 0;
        Cmpnts        qFlux;

        //find the distance from the wall where velocity is interpolated
        //a uniform distance from the wall is chosen for all ibm fluid points to prevent instabilities due to the stair step structure
        cellSize = ibm->interpDist;

        // distance from cell to the IBM element
        elemDist = fabs(nDot(nSub(ibF[c].pMin, cent[k][j][i]), eNorm));

        bPt = nScale(cellSize + elemDist, eNorm);        
        mSum(bPt, cent[k][j][i]);

        //check the neighbourhood of current ibm fluid cell for fluid cells
        //find the fluid - ibm fluid cell interface faces
        for (PetscInt k11=k-1; k11<k+2; k11++)
        for (PetscInt j11=j-1; j11<j+2; j11++)
        for (PetscInt i11=i-1; i11<i+2; i11++)
        {
            PetscInt      isFace = 0;

            if
            (
                (
                    k11>=1 && k11<mz-1 &&
                    j11>=1 && j11<my-1 &&
                    i11>=1 && i11<mx-1
                ) && isFluidCell(k11, j11, i11, nvert)
            )
            {
                //define the faces
                // 1 - k:k+1
                // 2 - k-1:k
                // 3 - j:j+1
                // 4 - j-1:j
                // 5 - i:i+1
                // 6 - i-1:i

                if(k11 == k+1 && j11 == j && i11 == i)
                {
                    isFace = 1;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k-1 && j11 == j && i11 == i)
                {
                    isFace = 2;
                    k1 = k-1; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j+1 && i11 == i)
                {
                    isFace = 3;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j-1 && i11 == i)
                {
                    isFace = 4;
                    k1 = k; j1 = j-1; i1 = i;
                }
                else if (k11 == k && j11 == j && i11 == i+1)
                {
                    isFace = 5;
                    k1 = k; j1 = j; i1 = i;
                }
                else if (k11 == k && j11 == j && i11 == i-1)
                {
                    isFace = 6;
                    k1 = k; j1 = j; i1 = i-1;
                }

                // only for the cells that share a face, not diagonal cells
                if(isFace)
                {

                    // interpolate the velocity of the projected point on the IBM solid element from its nodes
                    ibmPtVel.x =   ibMsh->nU[n1].x * ibF[c].cs1
                                 + ibMsh->nU[n2].x * ibF[c].cs2
                                 + ibMsh->nU[n3].x * ibF[c].cs3;

                    ibmPtVel.y =   ibMsh->nU[n1].y * ibF[c].cs1
                                 + ibMsh->nU[n2].y * ibF[c].cs2
                                 + ibMsh->nU[n3].y * ibF[c].cs3;

                    ibmPtVel.z =   ibMsh->nU[n1].z * ibF[c].cs1
                                 + ibMsh->nU[n2].z * ibF[c].cs2
                                 + ibMsh->nU[n3].z * ibF[c].cs3;

                    // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
                    PetscInt    ia, ja, ka;
                    PetscReal   dmin = 10e10, d;
                    PetscInt    ic, jc, kc, setFlag = 0;

                    // must be close so don't loop over all cells
                    for (ka=k-1; ka<k+2; ka++)
                    for (ja=j-1; ja<j+2; ja++)
                    for (ia=i-1; ia<i+2; ia++)
                    {

                        if
                        (
                            (
                                ka>=1 && ka<mz-1 &&
                                ja>=1 && ja<my-1 &&
                                ia>=1 && ia<mx-1
                            ) && isFluidCell(ka, ja, ia, nvert)
                        )
                        {
                            d = pow((bPt.x - cent[ka][ja][ia].x), 2) +
                                pow((bPt.y - cent[ka][ja][ia].y), 2) +
                                pow((bPt.z - cent[ka][ja][ia].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                ic = ia;
                                jc = ja;
                                kc = ka;
                                setFlag = 1;
                            }
                        }
                    }

                    if(setFlag == 0)
                    {
                        char error[512];
                        sprintf(error, "no closest fluid cell to ibm fluid %ld %ld %ld\n", k, j, i);
                        fatalErrorInFunction("findIBMWallShear",  error);
                    }  

                    PetscInt intId[6];
                    PetscInt intFlag = 1;

                    // get the trilinear interpolation cells
                    PointInterpolationCells
                    (
                            mesh,
                            bPt.x, bPt.y, bPt.z,
                            ic, jc, kc,
                            cent,
                            intId
                    );

                    // save the initial closest cell and background point
                    initCp.i = ic; initCp.j = jc; initCp.k = kc;
                    bPtInit  = nSet(bPt);                    

                    PetscReal sumDel = 0.0;

                    while ( (intFlag == 1) && isInsideBoundingBox(bPt, mesh->bounds) && (sumDel <= 1.5*cellSize))
                    {
                        PetscInt ibmCellCtr = 0;
                        PetscInt icc, jcc, kcc, setFlag = 0;

                        for (PetscInt kk = 0; kk<2; kk++)
                        for (PetscInt jj = 0; jj<2; jj++)
                        for (PetscInt ii = 0; ii<2; ii++)
                        {
                            if(isIBMCell(intId[kk], intId[jj+2], intId[ii+4], nvert))
                            {
                                ibmCellCtr ++;
                            }

                        }

                        if (ibmCellCtr > 0)
                        {
                            Cmpnts del =  nScale(0.2 * cellSize, eNorm);
                            mSum(bPt, del);
                            dmin = 10e10;

                            sumDel += nMag(del);

                            for (ka=kc-1; ka<kc+2; ka++)
                            {
                                //check processor ghost bounds
                                if (ka < gzs || ka >= gze) {intFlag = 2; break;}
                                for (ja=jc-1; ja<jc+2; ja++)
                                {
                                    if (ja < gys || ja >= gye) {intFlag = 2; break;}
                                    for (ia=ic-1; ia<ic+2; ia++)
                                    {
                                        if (ia < gxs || ia >= gxe) {intFlag = 2; break;}

                                        if
                                        (
                                            (
                                                ka>=1 && ka<mz-1 &&
                                                ja>=1 && ja<my-1 &&
                                                ia>=1 && ia<mx-1
                                            ) && (isFluidCell(ka, ja, ia, nvert))
                                        )
                                        {
                                            d = pow((bPt.x - cent[ka][ja][ia].x), 2) +
                                                pow((bPt.y - cent[ka][ja][ia].y), 2) +
                                                pow((bPt.z - cent[ka][ja][ia].z), 2);

                                            if
                                            (
                                                d < dmin
                                            )
                                            {
                                                dmin  = d;
                                                icc = ia;
                                                jcc = ja;
                                                kcc = ka;
                                                setFlag = 1;
                                            }
                                        }

                                    }
                                }
                            }

                            if(setFlag == 0)
                            {
                                //closest point not set, do not interpolate
                                intFlag = 2;
                            }

                            if(intFlag == 1)
                            {
                                kc = kcc; jc = jcc; ic = icc;

                                PointInterpolationCells
                                (
                                        mesh,
                                        bPt.x, bPt.y, bPt.z,
                                        ic, jc, kc,
                                        cent,
                                        intId
                                );
                            }

                        }
                        else
                        {
                            intFlag = 0;

                            // trilinear interpolate the velocity at this point
                            vectorPointLocalVolumeInterpolation
                            (
                                    mesh,
                                    bPt.x, bPt.y, bPt.z,
                                    ic, jc, kc,
                                    cent,
                                    ucat,
                                    bPtVel
                            );

                            if(flags->isTeqnActive)
                            {
                                scalarPointLocalVolumeInterpolation
                                (
                                        mesh,
                                        bPt.x, bPt.y, bPt.z,
                                        ic, jc, kc,
                                        cent,
                                        lt,
                                        bPtTemp
                                );
                            }
                        }
                    }

                    // while loop fails
                    if(intFlag > 0)
                    {
                        bPtVel = nSet(ucat[initCp.k][initCp.j][initCp.i]);

                        if(flags->isTeqnActive)
                        {
                            bPtTemp = lt[initCp.k][initCp.j][initCp.i];
                        }

                        bPt = nSet(bPtInit);

                        ic = initCp.i;
                        jc = initCp.j;
                        kc = initCp.k;
                    }

                    // distance of background point from the wall
                    sb = nDot(nSub(bPt, ibF[c].pMin), eNorm);

                    //velocity in local co-ordinate system
                    mSub(bPtVel, ibmPtVel);
                    uN = nScale(nDot(bPtVel, eNorm), eNorm);
                    uT = nSub(bPtVel, uN);
                    uTmag = nMag(uT);

                    if (flags->isTeqnActive)
                    {
                        if(ibm->ibmBody[ibF[c].bodyID]->tempBC == "thetaWallFunction")
                        {
                            if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -3)
                            {

                            }
                            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -2)
                            {
                                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelT->wmShumann;

                                qWall  = wm->qWall;
                            }
                            else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeT == -4)
                            {

                                Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelT->wmShumann;

                                // find interpolation weights for surface temp
                                PetscReal w[2];
                                PetscInt  l[2];

                                findInterpolationWeights(w, l, wm->timeVec, wm->numT, clock->time);

                                //interpolate the surface temperature from closest available mesoscale times
                                PetscReal surfaceTemp = w[0] * wm->surfTemp[l[0]] + w[1] * wm->surfTemp[l[1]];
                                PetscReal surfaceL;

                                PetscReal deltaTheta = bPtTemp - surfaceTemp;
                                PetscReal phiM, phiH;
                                
                                qWallShumann
                                (
                                    uTmag, sb, wm->roughness,
                                    wm->gammaM, wm->gammaH, wm->alphaH,
                                    wm->thetaRef, deltaTheta, wm->kappa,
                                    qWall, ustar, phiM, phiH, surfaceL
                                );
                            }
                        }
                    }

                    if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
                    {
                        if (ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
                        {
                            PetscReal roughness = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmCabot->roughness;

                            if(roughness > 1.0e-12)
                            {
                                ustar = uTauCabotRoughness(nu, uTmag, sb, 0.01, 0, roughness);

                            }
                            else
                            {
                                ustar = uTauCabot(nu, uTmag, sb, 0.01, 0);
                            }
                        }
                        else if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3)
                        {
                            Shumann *wm = ibm->ibmBody[ibF[c].bodyID]->ibmWallModelU->wmShumann;
                            PetscReal phiM, L;

                            uStarShumann
                            (
                                uTmag, sb, wm->roughness,
                                wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                                ustar, phiM, L
                            );
                        }
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall shear model is on but velocity BC is not set to velocityWallFunction\n");
                        fatalErrorInFunction("findIBMWallShear", error);
                    }

                    //local co-ordinate system
                    eT1 = nUnit(uT);
                    eN = eNorm;
                    eT2 = nCross(eN, eT1);

                    tauWall = ustar * ustar;

                    if(isFace == 5 || isFace == 6)
                    {
                        csi0 = icsi[k1][j1][i1].x, csi1 = icsi[k1][j1][i1].y, csi2 = icsi[k1][j1][i1].z;
                        eta0 = ieta[k1][j1][i1].x, eta1 = ieta[k1][j1][i1].y, eta2 = ieta[k1][j1][i1].z;
                        zet0 = izet[k1][j1][i1].x, zet1 = izet[k1][j1][i1].y, zet2 = izet[k1][j1][i1].z;
                        ajc  = iaj[k1][j1][i1];

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_i
                        (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert, meshTag,
                            &dudc, &dvdc, &dwdc,
                            &dude, &dvde, &dwde,
                            &dudz, &dvdz, &dwdz
                        );

                        // compute cartesian velocity derivatives w.r.t cartesian coords
                        Compute_du_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                            zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                            dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                            &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                        );

                        // transform cartesian velocity derivatives from cartesian coords to local wall normal coords
                        Compute_du_wmLocal
                        (
                            mesh,
                            eN, eT1, eT2,
                            du_dx, dv_dx, dw_dx,
                            du_dy, dv_dy, dw_dy,
                            du_dz, dv_dz, dw_dz,
                            &dut1dn, &dut2dn, &dundn,
                            &dut1dt1, &dut2dt1, &dundt1,
                            &dut1dt2, &dut2dt2, &dundt2
                        );

                        nut = lnu_t[k11][j11][i11];

                        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
                        {
                            if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3 || ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
                            {
                                //replace dut1dn
                                dut1dn = tauWall/(nu + nut) ;
                            }
                        }


                        //compute the transformation jacobian
                        Comput_JacobTensor_i
                        (
                            i1, j1, k1, mx, my, mz, coor,
                            &dxdc, &dxde, &dxdz,
                            &dydc, &dyde, &dydz,
                            &dzdc, &dzde, &dzdz
                        );

                        //transform from local wall model grid to computational grid
                        Compute_du_Compgrid
                        (
                            dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz,
                            eN.x, eN.y, eN.z, eT1.x, eT1.y, eT1.z, eT2.x, eT2.y, eT2.z,
                            dut1dn, dut2dn, dundn, dut1dt1, dut2dt1, dundt1, dut1dt2, dut2dt2, dundt2,
                            &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                        );

                        g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                        g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                        g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                        visc1[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nut + nu);
                        visc1[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nut + nu);
                        visc1[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nut + nu);

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].x = qFlux.x * icsi[k1][j1][i1].x + qFlux.y * icsi[k1][j1][i1].y + qFlux.z * icsi[k1][j1][i1].z;
                        }

                    }
                    else if(isFace == 3 || isFace == 4)
                    {
                        csi0 = jcsi[k1][j1][i1].x, csi1 = jcsi[k1][j1][i1].y, csi2 = jcsi[k1][j1][i1].z;
                        eta0 = jeta[k1][j1][i1].x, eta1 = jeta[k1][j1][i1].y, eta2 = jeta[k1][j1][i1].z;
                        zet0 = jzet[k1][j1][i1].x, zet1 = jzet[k1][j1][i1].y, zet2 = jzet[k1][j1][i1].z;
                        ajc =  jaj[k1][j1][i1];

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_j
                        (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert, meshTag, 
                            &dudc, &dvdc, &dwdc,
                            &dude, &dvde, &dwde,
                            &dudz, &dvdz, &dwdz
                        );

                        // compute cartesian velocity derivatives w.r.t cartesian coords
                        Compute_du_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                            zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                            dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                            &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                        );

                        // transform cartesian velocity derivatives from cartesian coords to local wall normal coords
                        Compute_du_wmLocal
                        (
                            mesh,
                            eN, eT1, eT2,
                            du_dx, dv_dx, dw_dx,
                            du_dy, dv_dy, dw_dy,
                            du_dz, dv_dz, dw_dz,
                            &dut1dn, &dut2dn, &dundn,
                            &dut1dt1, &dut2dt1, &dundt1,
                            &dut1dt2, &dut2dt2, &dundt2
                        );

                        nut = lnu_t[k11][j11][i11];

                        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
                        {
                            if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3 || ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
                            {
                                //replace dut1dn
                                dut1dn = tauWall/(nu + nut) ;
                            }
                        }

                        //compute the transformation jacobian
                        Comput_JacobTensor_j
                        (
                            i1, j1, k1, mx, my, mz, coor,
                            &dxdc, &dxde, &dxdz,
                            &dydc, &dyde, &dydz,
                            &dzdc, &dzde, &dzdz
                        );

                        //transform from local wall model grid to computational grid
                        Compute_du_Compgrid
                        (
                            dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz,
                            eN.x, eN.y, eN.z, eT1.x, eT1.y, eT1.z, eT2.x, eT2.y, eT2.z,
                            dut1dn, dut2dn, dundn, dut1dt1, dut2dt1, dundt1, dut1dt2, dut2dt2, dundt2,
                            &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                        );

                        g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
                        g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
                        g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                        visc2[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nut + nu);
                        visc2[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nut + nu);
                        visc2[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nut + nu);

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].y = qFlux.x * jeta[k1][j1][i1].x + qFlux.y * jeta[k1][j1][i1].y + qFlux.z * jeta[k1][j1][i1].z;
                        }

                    }
                    else if(isFace == 1 || isFace == 2)
                    {
                        csi0 = kcsi[k1][j1][i1].x, csi1 = kcsi[k1][j1][i1].y, csi2 = kcsi[k1][j1][i1].z;
                        eta0 = keta[k1][j1][i1].x, eta1 = keta[k1][j1][i1].y, eta2 = keta[k1][j1][i1].z;
                        zet0 = kzet[k1][j1][i1].x, zet1 = kzet[k1][j1][i1].y, zet2 = kzet[k1][j1][i1].z;
                        ajc = kaj[k1][j1][i1];

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_k
                        (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert, meshTag, 
                            &dudc, &dvdc, &dwdc,
                            &dude, &dvde, &dwde,
                            &dudz, &dvdz, &dwdz
                        );

                        // compute cartesian velocity derivatives w.r.t cartesian coords
                        Compute_du_dxyz
                        (
                            mesh,
                            csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                            zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                            dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                            &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                        );

                        // transform cartesian velocity derivatives from cartesian coords to local wall normal coords
                        Compute_du_wmLocal
                        (
                            mesh,
                            eN, eT1, eT2,
                            du_dx, dv_dx, dw_dx,
                            du_dy, dv_dy, dw_dy,
                            du_dz, dv_dz, dw_dz,
                            &dut1dn, &dut2dn, &dundn,
                            &dut1dt1, &dut2dt1, &dundt1,
                            &dut1dt2, &dut2dt2, &dundt2
                        );

                        nut = lnu_t[k11][j11][i11];

                        if(ibm->ibmBody[ibF[c].bodyID]->velocityBC == "velocityWallFunction")
                        {
                            if(ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -3 || ibm->ibmBody[ibF[c].bodyID]->wallFunctionTypeU == -1)
                            {
                                //replace dut1dn
                                dut1dn = tauWall/(nu + nut) ;
                            }
                        }

                        //compute the transformation jacobian
                        Comput_JacobTensor_k
                        (
                            i1, j1, k1, mx, my, mz, coor,
                            &dxdc, &dxde, &dxdz,
                            &dydc, &dyde, &dydz,
                            &dzdc, &dzde, &dzdz
                        );

                        //transform from local wall model grid to computational grid
                        Compute_du_Compgrid
                        (
                            dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz,
                            eN.x, eN.y, eN.z, eT1.x, eT1.y, eT1.z, eT2.x, eT2.y, eT2.z,
                            dut1dn, dut2dn, dundn, dut1dt1, dut2dt1, dundt1, dut1dt2, dut2dt2, dundt2,
                            &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                        );

                        g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                        g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                        g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;

                        visc3[k][j][i].x = (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nut + nu);
                        visc3[k][j][i].y = (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nut + nu);
                        visc3[k][j][i].z = (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nut + nu);

                        if (flags->isTeqnActive)
                        {
                            qFlux = nScale(-qWall, eNorm);
                            viscT[k][j][i].z = qFlux.x * kzet[k1][j1][i1].x + qFlux.y * kzet[k1][j1][i1].y + qFlux.z * kzet[k1][j1][i1].z;
                        }
                    }


                }

            }
        }

    }

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    DMDAVecRestoreArray(fda, ueqn->lViscIBM1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lViscIBM2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lViscIBM3, &visc3);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
	DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
	DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

	DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
	DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
	DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

	DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
	DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
	DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
	DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
	DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    if(ibm->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(fda, teqn->lViscIBMT, &viscT);
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMLocalToLocalBegin(fda, teqn->lViscIBMT, INSERT_VALUES, teqn->lViscIBMT);
        DMLocalToLocalEnd  (fda, teqn->lViscIBMT, INSERT_VALUES, teqn->lViscIBMT);
    }


    DMLocalToLocalBegin(fda, ueqn->lViscIBM1, INSERT_VALUES, ueqn->lViscIBM1);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM1, INSERT_VALUES, ueqn->lViscIBM1);
    DMLocalToLocalBegin(fda, ueqn->lViscIBM2, INSERT_VALUES, ueqn->lViscIBM2);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM2, INSERT_VALUES, ueqn->lViscIBM2);
    DMLocalToLocalBegin(fda, ueqn->lViscIBM3, INSERT_VALUES, ueqn->lViscIBM3);
    DMLocalToLocalEnd  (fda, ueqn->lViscIBM3, INSERT_VALUES, ueqn->lViscIBM3);

    return (0);
}

//***************************************************************************************************************//
PetscErrorCode ibmSearch(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    clock_        *clock = ibm->access->clock;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i, j, k, b, k1, j1, i1;
    PetscInt      ip, im, jp, jm, kp, km;
    PetscInt      ii, jj, kk;
    PetscMPIInt   rank, nProcs;

    cellIds       sCell;

    Vec           lCoor;

    PetscReal     ***nvert, ***nvert_o, ***gnvert, ***nvertFixed;
    Cmpnts        ***cent, ***coor;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    MPI_Comm_size(mesh->MESH_COMM, &nProcs);

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {

        boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
        ibmMesh       *ibMsh = ibm->ibmBody[b]->ibMsh;                         // pointer to the ibm body mesh

        searchBox *sBox           = &(ibm->sBox[b]);
        list      *searchCellList = ibm->ibmBody[b]->searchCellList;

        if(ibm->ibmBody[b]->thinBody)
        {
            DMDAVecGetArray(fda, mesh->lCent, &cent);
            DMDAVecGetArray(da, mesh->lNvert, &nvert);
            DMDAVecGetArray(da, mesh->lNvert_o, &nvert_o);

            // get mesh coordinates
            DMGetCoordinatesLocal(da, &lCoor);
            DMDAVecGetArray(fda, lCoor, &coor);

            for (k = lzs; k < lze; k++)
            for (j = lys; j < lye; j++)
            for (i = lxs; i < lxe; i++)
            {

                PetscInt checkCell = 0;

                //loop through the co-ordinates of the cell
                for (k1 = k-1; k1 < k+1; k1++)
                for (j1 = j-1; j1 < j+1; j1++)
                for (i1 = i-1; i1 < i+1; i1++)
                {
                    // Only if the fluid mesh cell center coordinates are in the IBM bounding box
                    // a fluid mesh cell can be inside the bounding box and still not be inside the IBM body due to its shape
                    if(   coor[k1][j1][i1].x > ibBox->xmin
                    && coor[k1][j1][i1].x < ibBox->xmax
                    && coor[k1][j1][i1].y > ibBox->ymin
                    && coor[k1][j1][i1].y < ibBox->ymax
                    && coor[k1][j1][i1].z > ibBox->zmin
                    && coor[k1][j1][i1].z < ibBox->zmax )
                    {
                        checkCell = 1;
                    }
                }

                if(checkCell == 1)
                {
                    // index of the cells neighbouring the cell i,j,k
                    ip = (i < mx - 2 ? (i + 1) : (i));
                    im = (i > 1 ? (i - 1) : (i));

                    jp = (j < my - 2 ? (j + 1) : (j));
                    jm = (j > 1 ? (j - 1) : (j));

                    kp = (k < mz - 2 ? (k + 1) : (k));
                    km = (k > 1 ? (k - 1) : (k));

                    PetscInt dosearch = 0, sign = 0;

                    //do a search if there is a change in the nvert value between the i,j,k cell and its neighbours due to IBM movement
                    for (kk = km; kk <= kp; kk++)
                    for (jj = jm; jj <= jp; jj++)
                    for (ii = im; ii <= ip; ii++)
                    {
                        PetscReal sign = nvert_o[k][j][i] - nvert_o[kk][jj][ii];
                        if (fabs(sign) > 1.e-6)
                        {
                        dosearch += 1;
                        }
                    }

                if ( (dosearch && ibm->ibmBody[b]->bodyMotion != "static") || (clock->it == clock->itStart))
                {

                        PetscReal val = 0, valPrev = 0;
                        PetscInt  ctr = 1;
                        PetscInt  bdaryCell = 0;

                        // do the ray casting test to check if a cell is inside or outside an IBM body
                        //ray casting for all the co-ordinates of the cell borders
                        //loop through the co-ordinates of the cell
                        for (k1 = k-1; k1 < k+1; k1++)
                        {
                            for (j1 = j-1; j1 < j+1; j1++)
                            {
                                for (i1 = i-1; i1 < i+1; i1++)
                                {
                                    // find the search cell were the fluid node is located
                                    sCell.i = floor((coor[k1][j1][i1].x - ibBox->xmin) / sBox->dcx);
                                    sCell.j = floor((coor[k1][j1][i1].y - ibBox->ymin) / sBox->dcy);
                                    sCell.k = floor((coor[k1][j1][i1].z - ibBox->zmin) / sBox->dcz);


                                    if(sCell.i < 0) sCell.i = 0;
                                    if(sCell.i > sBox->ncx-1) sCell.i = sBox->ncx-1;
                                    if(sCell.j < 0) sCell.j = 0;
                                    if(sCell.j > sBox->ncy-1) sCell.j = sBox->ncy-1;
                                    if(sCell.k < 0) sCell.k = 0;
                                    if(sCell.k > sBox->ncz-1) sCell.k = sBox->ncz-1;

                                    //save previous value
                                    valPrev = val;

                                    val = rayCastingTest(coor[k1][j1][i1], ibMsh, sCell, sBox, ibBox, searchCellList);

                                    //skip initial iteration of loop
                                    if(ctr != 1)
                                    {
                                        //if one of the co-ordinates is outside and other is inside the body, then it is a boundary cell
                                        if(valPrev!=val)
                                        {
                                            bdaryCell = 1;
                                            break;
                                        }
                                    }

                                    ctr++;
                                }

                                if(bdaryCell) break;
                            }

                            if(bdaryCell) break;
                        }


                        if(bdaryCell)
                        {
                            // find the search cell were the fluid node is located
                            sCell.i = floor((cent[k][j][i].x - ibBox->xmin) / sBox->dcx);
                            sCell.j = floor((cent[k][j][i].y - ibBox->ymin) / sBox->dcy);
                            sCell.k = floor((cent[k][j][i].z - ibBox->zmin) / sBox->dcz);


                            if(sCell.i < 0) sCell.i = 0;
                            if(sCell.i > sBox->ncx-1) sCell.i = sBox->ncx-1;
                            if(sCell.j < 0) sCell.j = 0;
                            if(sCell.j > sBox->ncy-1) sCell.j = sBox->ncy-1;
                            if(sCell.k < 0) sCell.k = 0;
                            if(sCell.k > sBox->ncz-1) sCell.k = sBox->ncz-1;

                            //check again now for the cell center
                            val = rayCastingTest(cent[k][j][i], ibMsh, sCell, sBox, ibBox, searchCellList);

                            if(val == 4.0)
                            {
                                //boundary cell is solid
                            }
                            else
                            {
                                //boundary cell is fluid. set to ibm fluid
                                val = 2.0;
                            }
                        }

                        nvert[k][j][i] = PetscMax(nvert[k][j][i], val);
                }
                else
                {   // set nvert to old value
                        nvert[k][j][i] = nvert_o[k][j][i];

                        if (PetscInt (nvert[k][j][i]+0.5) ==3)
                            nvert[k][j][i] = 4;             // only inside solid is set, here = 4
                        if (PetscInt (nvert[k][j][i]+0.5) ==1)
                            nvert[k][j][i] = 2;             // near boundary values are set to 0 and need to be recalculated

                    }
                }
            }

            DMDAVecRestoreArray(fda, mesh->lCent, &cent);
            DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
            DMDAVecRestoreArray(da, mesh->lNvert_o, &nvert_o);

            MPI_Barrier(mesh->MESH_COMM);

            DMLocalToLocalBegin(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
            DMLocalToLocalEnd(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);

            DMDAVecRestoreArray(fda, lCoor, &coor);

        }
        else
        {
            DMDAVecGetArray(fda, mesh->lCent, &cent);
            DMDAVecGetArray(da, mesh->lNvert, &nvert);
            DMDAVecGetArray(da, mesh->lNvert_o, &nvert_o);
            DMDAVecGetArray(da, ibm->lNvertFixed, &nvertFixed);

            for (k = lzs; k < lze; k++)
            for (j = lys; j < lye; j++)
            for (i = lxs; i < lxe; i++)
            {

            // Only if the fluid mesh cell center coordinates are in the IBM bounding box
            // a fluid mesh cell can be inside the bounding box and still not be inside the IBM body due to its shape
            if(cent[k][j][i].x > ibBox->xmin
                && cent[k][j][i].x < ibBox->xmax
                && cent[k][j][i].y > ibBox->ymin
                && cent[k][j][i].y < ibBox->ymax
                && cent[k][j][i].z > ibBox->zmin
                && cent[k][j][i].z < ibBox->zmax )
                {

                    // index of the cells neighbouring the cell i,j,k
                    ip = (i < mx - 2 ? (i + 1) : (i));
                    im = (i > 1 ? (i - 1) : (i));

                    jp = (j < my - 2 ? (j + 1) : (j));
                    jm = (j > 1 ? (j - 1) : (j));

                    kp = (k < mz - 2 ? (k + 1) : (k));
                    km = (k > 1 ? (k - 1) : (k));

                    PetscInt dosearch = 0, sign = 0;

                    //do a search if there is a change in the nvert value between the i,j,k cell and its neighbours due to IBM movement
                    for (kk = km; kk <= kp; kk++)
                    for (jj = jm; jj <= jp; jj++)
                    for (ii = im; ii <= ip; ii++)
                    {
                        PetscReal sign = nvert_o[k][j][i] - nvert_o[kk][jj][ii];
                        if (fabs(sign) > 1.e-6)
                        {
                        dosearch += 1;
                        }
                    }

                if ( (dosearch && ibm->ibmBody[b]->bodyMotion != "static") || (clock->it == clock->itStart))
                {
                        // find the search cell were the fluid node is located
                        sCell.i = floor((cent[k][j][i].x - ibBox->xmin) / sBox->dcx);
                        sCell.j = floor((cent[k][j][i].y - ibBox->ymin) / sBox->dcy);
                        sCell.k = floor((cent[k][j][i].z - ibBox->zmin) / sBox->dcz);

                        // do the ray casting test to check if a cell is inside or outside an IBM body
                        PetscReal val;
                        val = rayCastingTest(cent[k][j][i], ibMsh, sCell, sBox, ibBox, searchCellList);
                        nvert[k][j][i] = PetscMax(nvert[k][j][i], val);

                        //save the nvert for static bodies
                        if( (clock->it == clock->itStart) && (ibm->ibmBody[b]->bodyMotion == "static"))
                        {
                            nvertFixed[k][j][i] = PetscMax(nvertFixed[k][j][i], val);
                        }

                }
                else
                    {
                        if(ibm->ibmBody[b]->bodyMotion == "static")
                        {
                            nvert[k][j][i] = PetscMax(nvertFixed[k][j][i], nvert[k][j][i]);
                        }
                        else
                        {
                            nvert[k][j][i] = nvert_o[k][j][i];
                        }

                        if (PetscInt (nvert[k][j][i]+0.5) ==3)
                            nvert[k][j][i] = 4;             // only inside solid is set, here = 4
                        if (PetscInt (nvert[k][j][i]+0.5) ==1)
                            nvert[k][j][i] = 0;             // near boundary values are set to 0 and need to be recalculated

                }
                }
            }

            DMDAVecRestoreArray(fda, mesh->lCent, &cent);
            DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
            DMDAVecRestoreArray(da, mesh->lNvert_o, &nvert_o);
            DMDAVecRestoreArray(da, ibm->lNvertFixed, &nvertFixed);

            MPI_Barrier(mesh->MESH_COMM);

            DMLocalToLocalBegin(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
            DMLocalToLocalEnd(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
        }

        // set nvert at solid fluid intersection IB Nodes
        DMDAVecGetArray(da, mesh->lNvert, &nvert);

        if(ibm->wallShearOn)
        {
            for (k = lzs; k < lze; k++)
            for (j = lys; j < lye; j++)
            for (i = lxs; i < lxe; i++)
            {

                if (nvert[k][j][i] < 0)
                {
                    nvert[k][j][i] = 0;
                }

                ip = (i < mx - 1 ? (i + 1) : (i));
                im = (i > 0 ? (i - 1) : (i));

                jp = (j < my - 1 ? (j + 1) : (j));
                jm = (j > 0 ? (j - 1) : (j));

                kp = (k < mz - 1 ? (k + 1) : (k));
                km = (k > 0 ? (k - 1) : (k));

                //only cells that share a face with ibm solid are made ibm fluid
                if ((PetscInt) (nvert[k][j][i] + 0.5) != 0)
                {
                    for (kk = km; kk < kp + 1; kk++)
                    for (jj = jm; jj < jp + 1; jj++)
                    for (ii = im; ii < ip + 1; ii++)
                    {
                        if ((PetscInt) (nvert[kk][jj][ii] + 0.5) == 0)
                        {
                            if(kk == k+1 && jj == j && ii == i)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                            else if(kk == k-1 && jj == j && ii == i)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                            else if(kk == k && jj == j+1 && ii == i)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                            else if(kk == k && jj == j-1 && ii == i)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                            else if(kk == k && jj == j && ii == i+1)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                            else if(kk == k && jj == j && ii == i-1)
                            {
                                nvert[k][j][i] = 2.0;
                            }
                        }
                    }
                }
            }
        }
        else
        {
            for (k = lzs; k < lze; k++)
            for (j = lys; j < lye; j++)
            for (i = lxs; i < lxe; i++)
            {

                if (nvert[k][j][i] < 0)
                {
                    nvert[k][j][i] = 0;
                }

                ip = (i < mx - 1 ? (i + 1) : (i));
                im = (i > 0 ? (i - 1) : (i));

                jp = (j < my - 1 ? (j + 1) : (j));
                jm = (j > 0 ? (j - 1) : (j));

                kp = (k < mz - 1 ? (k + 1) : (k));
                km = (k > 0 ? (k - 1) : (k));

                //only cells that share a face with ibm solid are made ibm fluid
                if ((PetscInt) (nvert[k][j][i] + 0.5) != 4)
                {
                    for (kk = km; kk < kp + 1; kk++)
                    for (jj = jm; jj < jp + 1; jj++)
                    for (ii = im; ii < ip + 1; ii++)
                    {
                        if ((PetscInt) (nvert[kk][jj][ii] + 0.5) == 4)
                        {
                            if(kk == k+1 && jj == j && ii == i)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                            else if(kk == k-1 && jj == j && ii == i)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                            else if(kk == k && jj == j+1 && ii == i)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                            else if(kk == k && jj == j-1 && ii == i)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                            else if(kk == k && jj == j && ii == i+1)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                            else if(kk == k && jj == j && ii == i-1)
                            {
                                nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                            }
                        }
                    }
                }

                // if ((PetscInt) (nvert[k][j][i] + 0.5) != 4)
                // {
                //     for (kk = km; kk < kp + 1; kk++)
                //     for (jj = jm; jj < jp + 1; jj++)
                //     for (ii = im; ii < ip + 1; ii++)
                //     {
                //         if ((PetscInt) (nvert[kk][jj][ii] + 0.5) == 4)
                //         {
                //             nvert[k][j][i] = PetscMax(2.0, nvert[k][j][i]);
                //         }
                //     }
                // }
            }
        }

        DMDAVecGetArray(da, mesh->lNvert_o, &nvert_o);

        //check for cells that change phase (solid to fluid or viceversa) without being an ibm fluid cell
        if(b == ibm->numBodies-1)
        {
            for (k = lzs; k < lze; k++)
            for (j = lys; j < lye; j++)
            for (i = lxs; i < lxe; i++)
            {
                if (nvert_o[k][j][i] > 2.5 && nvert[k][j][i] < 0.5)
                {
                    //PetscPrintf(PETSC_COMM_SELF, "phase change for element %ld %ld %ld - nvert_o = %lf, nvert = %lf\n", k, j, i, nvert_o[k][j][i], nvert[k][j][i]);
                    nvert[k][j][i]=2;
                }
            }
        }

        DMDAVecRestoreArray(da, mesh->lNvert_o, &nvert_o);
        DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

        DMLocalToLocalBegin(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
        DMLocalToLocalEnd(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);

        // nvert values of 4 for solid and 2 for IB fluid nodes where used for the current body
        // to differentiate it from other bodies.
        // reset back to nvert values of 3 for solid and 1 for IB fluid
        DMDAVecGetArray(da, mesh->lNvert, &nvert);

        // Back to the old nvert 3 and 1
        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {
            if ((PetscInt) (nvert[k][j][i] + 0.5) == 2)
            {
                nvert[k][j][i] = 1;
            }
            if ((PetscInt) (nvert[k][j][i] + 0.5) == 4)
            {
                nvert[k][j][i] = 3;
            }

        }


        DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

        DMLocalToLocalBegin(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
        DMLocalToLocalEnd(da, mesh->lNvert, INSERT_VALUES, mesh->lNvert);
    }

    //ibm nvert cleanup - make solid, ibm fluid cells that dont have even one fluid cell around it
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->Nvert, &gnvert);

    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
        PetscInt ctr = 0;

        ip = (i < mx - 2 ? (i + 1) : (i));
        im = (i > 1 ? (i - 1) : (i));

        jp = (j < my - 2 ? (j + 1) : (j));
        jm = (j > 1 ? (j - 1) : (j));

        kp = (k < mz - 2 ? (k + 1) : (k));
        km = (k > 1 ? (k - 1) : (k));

        if (isIBMFluidCell(k, j, i, nvert))
        {

            for (kk = km; kk < kp + 1; kk++)
            for (jj = jm; jj < jp + 1; jj++)
            for (ii = im; ii < ip + 1; ii++)
            {
                if (isFluidCell(kk, jj, ii, nvert))
                {
                    ctr ++;
                }
            }

            if(ctr == 0)
            {
                nvert[k][j][i] = 3.0;
            }
        }

        //set the global nvert
        gnvert[k][j][i] = nvert[k][j][i];
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->Nvert, &gnvert);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    //apply periodicity or ghost node values
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->Nvert, &gnvert);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                PetscInt a=i, b=j, c=k, flag=0;

                if(i==0)
                {
                    if(mesh->i_periodic)       a=mx-2, flag=1;
                    else if(mesh->ii_periodic) a=-2, flag=1;
                    else                      a=1, flag=1;
                }
                if(i==mx-1)
                {
                    if(mesh->i_periodic)       a=1, flag=1;
                    else if(mesh->ii_periodic) a=mx+1, flag=1;
                    else                      a=mx-2, flag=1;
                }
                if(j==0)
                {
                    if(mesh->j_periodic)       b=my-2, flag=1;
                    else if(mesh->jj_periodic) b=-2, flag=1;
                    else                      b=1, flag=1;
                }
                if(j==my-1)
                {
                    if(mesh->j_periodic)       b=1, flag=1;
                    else if(mesh->jj_periodic) b=my+1, flag=1;
                    else                      b=my-2, flag=1;
                }
                if(k==0)
                {
                    if(mesh->k_periodic)       c=mz-2, flag=1;
                    else if(mesh->kk_periodic) c=-2, flag=1;
                    else                      c=1, flag=1;
                }
                if(k==mz-1)
                {
                    if(mesh->k_periodic)       c=1, flag=1;
                    else if(mesh->kk_periodic) c=mz+1, flag=1;
                    else                      c=mz-2, flag=1;
                }

                if(flag)
                {
                    gnvert[k][j][i] = nvert[c][b][a];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->Nvert, &gnvert);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    MPI_Barrier(mesh->MESH_COMM);
    
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode MLSInterpolation(ibm_ *ibm)
{
  mesh_         *mesh  = ibm->access->mesh;
  ueqn_         *ueqn  = ibm->access->ueqn;
  teqn_         *teqn  = ibm->access->teqn;
  flags_        *flags = ibm->access->flags;
  DM            da = mesh->da, fda = mesh->fda;
  DMDALocalInfo info = mesh->info;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      i, j, k, ii, kk, jj;
  PetscInt      b, c, s;
  PetscReal     normDist, ***nvert;

  Cmpnts        ***ucat, ***lucat, ***cent, dist;
  PetscReal     ***lt, ***temp;

  PetscInt      n = 0, nsupport=0;

  PetscReal     dis;
  PetscReal     *W, *PHI, **P, **B, *Ux, *Uy, *Uz, *tmp;
  PetscReal     **A, **inv_A;

  DMDAVecGetArray(fda, mesh->lCent, &cent);
  DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
  DMDAVecGetArray(fda, ueqn->Ucat, &ucat);

  if (flags->isTeqnActive)
  {
      DMDAVecGetArray(da, teqn->lTmprt, &lt);
      DMDAVecGetArray(da, teqn->Tmprt, &temp);
  }

  //local pointer for ibmFluidCells
  ibmFluidCell *ibF = ibm->ibmFCells;

  for(c = 0; c < ibm->numIBMFluid; c++)
  {
    //cell indices
    i = ibF[c].cellId.i;
    j = ibF[c].cellId.j;
    k = ibF[c].cellId.k;

    //total number of support cells - fluid and ibm mesh
    nsupport = ibF[c].numFl + ibF[c].numSl;

    // allocate local variables for MLS interpolation
    B = (PetscReal**) malloc(4*sizeof(PetscReal*));

    for (n=0;n<4;n++)
    {
        B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
    }

    P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));

    for (n=0;n<nsupport;n++)
    {
        P[n] = (PetscReal*) malloc(4*sizeof(PetscReal));
    }

    W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

    PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

    Ux = (PetscReal* ) malloc(nsupport*sizeof(PetscReal));
    Uy = (PetscReal* ) malloc(nsupport*sizeof(PetscReal));
    Uz = (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

    PetscMalloc(sizeof(PetscReal *)  * (4), &(A));
    PetscMalloc(sizeof(PetscReal *)  * (4), &(inv_A));

    for(n=0; n<4; n++)
    {
        PetscMalloc(sizeof(PetscReal)  * (4), &(A[n]));
        PetscMalloc(sizeof(PetscReal)  * (4), &(inv_A[n]));
    }

    if (flags->isTeqnActive)
    {
        tmp = (PetscReal* ) malloc(nsupport*sizeof(PetscReal));
    }

    n = -1;

    //build the fluid support polynomial basis
    //get access to the first fluid cell in the list
    cellNode *flList = ibF[c].flNodes.head;

    while(flList)
    {
      n++;

      ii = flList->Node.i;
      jj = flList->Node.j;
      kk = flList->Node.k;

      P[n][0] = 1.0;
      P[n][1] = (cent[kk][jj][ii].x - cent[k][j][i].x) / ibF[c].rad;
      P[n][2] = (cent[kk][jj][ii].y - cent[k][j][i].y) / ibF[c].rad;
      P[n][3] = (cent[kk][jj][ii].z - cent[k][j][i].z) / ibF[c].rad;

      Ux[n] = lucat[kk][jj][ii].x;
      Uy[n] = lucat[kk][jj][ii].y;
      Uz[n] = lucat[kk][jj][ii].z;

      if(flags->isTeqnActive)
      {
          tmp[n] = lt[kk][jj][ii];
      }

      // get normalized distance
      dist = nSub(cent[k][j][i], cent[kk][jj][ii]);
      normDist = nMag(dist)/ibF[c].rad;

      // get interpolation weights
      W[n]
      =
      (normDist < 0.5)
      ?
      2.0 / 3.0 - 4.0 * normDist * normDist + 4.0 * pow(normDist, 3.0)
      :
      4.0 / 3.0 - 4.0 * normDist + 4.0 * pow(normDist, 2.0) - 4.0 / 3.0 * pow(normDist, 3.0);

      flList = flList->next;
    }

    //build the ibm mesh node support polynomial basis
    //get access to the first ibm mesh node in the list
    node *slList = ibF[c].slNodes.head;

    // get access to the first ibm mesh body in the list
    node *bID = ibF[c].bID.head;

    while(slList)
    {
      n++;

      Cmpnts nCoor = ibm->ibmBody[bID->Node]->ibMsh->nCoor[slList->Node];
      Cmpnts uIBM = ibm->ibmBody[bID->Node]->ibMsh->nU[slList->Node];

      P[n][0]=1.0;

      P[n][1]=(nCoor.x-cent[k][j][i].x)/ibF[c].rad;
      P[n][2]=(nCoor.y-cent[k][j][i].y)/ibF[c].rad;
      P[n][3]=(nCoor.z-cent[k][j][i].z)/ibF[c].rad;

      Ux[n]=uIBM.x;
      Uy[n]=uIBM.y;
      Uz[n]=uIBM.z;

      if(flags->isTeqnActive)
      {
          tmp[n] = 20;//ibm->access->abl->tRef;
      }

      // get normalized distance
      dist = nSub(nCoor, cent[k][j][i]);
      normDist = nMag(dist)/ibF[c].rad;

      // get interpolation weights
      W[n]
      =
      (normDist < 0.5)
      ?
      2.0 / 3.0 - 4.0 * normDist * normDist + 4.0 * pow(normDist, 3.0)
      :
      4.0 / 3.0 - 4.0 * normDist + 4.0 * pow(normDist, 2.0) - 4.0 / 3.0 * pow(normDist, 3.0);

      slList = slList->next;
      bID = bID->next;
    }

    for (ii=0; ii<4; ii++)
    {
        for (jj=0; jj<nsupport; jj++)
        {
            B[ii][jj] = W[jj] * P[jj][ii];
        }
    }

    for (ii=0; ii<4; ii++)
    {
        for (jj=0; jj<4; jj++)
        {
            A[ii][jj] = 0.;

            for (kk=0; kk<nsupport; kk++)
            {
                A[ii][jj] += B[ii][kk] * P[kk][jj];
            }
        }
    }

    inv_4by4(A, inv_A, 4);

    // get the values of the shape function PHI on the support points
    mult_mats3_lin(inv_A, B, nsupport, PHI);

    ucat[k][j][i].x = ucat[k][j][i].y = ucat[k][j][i].z = 0.0;

    if(flags->isTeqnActive)
    {
        temp[k][j][i] = 0.0;
    }

    // interpolate with linear least squares
    for ( n = 0; n < nsupport; n++)
    {
        // velocity
        ucat[k][j][i].x += PHI[n]*Ux[n];
        ucat[k][j][i].y += PHI[n]*Uy[n];
        ucat[k][j][i].z += PHI[n]*Uz[n];

        // temperature
        if(flags->isTeqnActive)
        {
            temp[k][j][i] += PHI[n] * tmp[n];
        }
    }

    for( n = 0; n < 4; n++)
    {
        free(B[n]);
        free(A[n]);
        free(inv_A[n]);
    }

    for(n = 0; n < 4; n++)
    {
        free(P[n]);
    }

    free(B);
    free(P);
    free(W);
    free(PHI);
    free(Ux);
    free(Uy);
    free(Uz);
    free(A);
    free(inv_A);

    if(flags->isTeqnActive)
    {
        free(tmp);
    }


  }

  DMDAVecRestoreArray(fda, mesh->lCent, &cent);
  DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
  DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);

  if (flags->isTeqnActive)
  {
      DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
      DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
      DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
      DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

  }

  DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
  DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

  return 0;
}

//***************************************************************************************************************//

PetscErrorCode checkIBMexists(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, b;
    PetscInt      ip, im, jp, jm, kp, km, kk, jj, ii;
    PetscReal     ***nvert;
    PetscInt      lsum = 0., sum = 0.;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    //check if atleast one ibm found in the domain
    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
      if(isIBMCell(k, j, i, nvert))
      {
         lsum ++ ;
      }
    }

    MPI_Allreduce(&lsum, &sum, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    if(sum == 0)
    {
      char error[512];
      sprintf(error, "no ibm solid or fluid cell found in the mesh\n");
      fatalErrorInFunction("checkIBMexists", error);
    }

    //check that a solid ibm cell is definitely surrounded by ibm fluid cells
    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
        ip = (i < mx - 1 ? (i + 1) : (i));
        im = (i > 0 ? (i - 1) : (i));

        jp = (j < my - 1 ? (j + 1) : (j));
        jm = (j > 0 ? (j - 1) : (j));

        kp = (k < mz - 1 ? (k + 1) : (k));
        km = (k > 0 ? (k - 1) : (k));

        PetscInt sumCheck = 0;
        if(isIBMSolidCell(k, j, i, nvert))
        {
            for (kk = km; kk < kp + 1; kk++)
            for (jj = jm; jj < jp + 1; jj++)
            for (ii = im; ii < ip + 1; ii++)
            {
                if(!isIBMSolidCell(kk, jj, ii, nvert))
                {
                    if(kk == k+1 && jj == j && ii == i)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                    else if(kk == k-1 && jj == j && ii == i)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                    else if(kk == k && jj == j+1 && ii == i)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                    else if(kk == k && jj == j-1 && ii == i)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                    else if(kk == k && jj == j && ii == i+1)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                    else if(kk == k && jj == j && ii == i-1)
                    {
                        if(!isIBMFluidCell(kk, jj, ii, nvert)) sumCheck = 1;
                    }
                }

                if(sumCheck == 1)
                {
                    char error[512];
                    sprintf(error, "ibm fluid cell found next to ibm solid cell %ld %ld %ld\n", k, j, i);
                    fatalErrorInFunction("checkIBMexists", error);
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    //check ibm body bounding box is inside the ibm processor bounds

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody  = ibm->ibmBody[b];

        if(ibm->dynamic)
        {
            PetscReal xmin_b = ibmBody->procBoundCenter.x - 0.5 * ibmBody->procBoundSize.x,
            xmax_b = ibmBody->procBoundCenter.x + 0.5 * ibmBody->procBoundSize.x,
            ymin_b = ibmBody->procBoundCenter.y - 0.5 * ibmBody->procBoundSize.y,
            ymax_b = ibmBody->procBoundCenter.y + 0.5 * ibmBody->procBoundSize.y,
            zmin_b = ibmBody->procBoundCenter.z - 0.5 * ibmBody->procBoundSize.z,
            zmax_b = ibmBody->procBoundCenter.z + 0.5 * ibmBody->procBoundSize.z;

            PetscInt sumCheck = 0;

            if(ibmBody->bound->xmin < xmin_b)      sumCheck = 1;
            else if(ibmBody->bound->xmax > xmax_b) sumCheck = 1;
            else if(ibmBody->bound->ymin < ymin_b) sumCheck = 1;
            else if(ibmBody->bound->ymax > ymax_b) sumCheck = 1;
            else if(ibmBody->bound->zmin < zmin_b) sumCheck = 1;
            else if(ibmBody->bound->zmax > zmax_b) sumCheck = 1;

            if(sumCheck == 1)
            {
                char error[512];
                sprintf(error, "For body: %s ibm processor bounds set are smaller than its bounding box limits\n", ibmBody->bodyName.c_str());
                fatalErrorInFunction("checkIBMexists", error);
            }
        }
    }
    return 0;
}

//***************************************************************************************************************//
PetscErrorCode createHalfEdgeDataStructure(ibm_ *ibm)
{
    /*  this function creates the half edge data structure which provides efficient inter-access to an IBM element vertex, half edge and faces.
        for definition of vertex, half-edge or faces check ibmInput.h
        vertex stores its coordinates and pointer to one of the half-edges that originates from it
        half-edge stores pointers to its origin vertex, next half-edge, twin half-edge and the face it is part of
        face stores a pointer to one of the half-edge within it
        For introduction to understanding this data structure refer this link: https://jerryyin.info/geometry-processing-algorithms/half-edge/
    */
    PetscInt         i, b;
    PetscInt         nodes;                          // number of ib nodes
    PetscInt         elems;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {
        ibmMesh  *ibMesh = ibm->ibmBody[b]->ibMsh;

        nodes    =    ibMesh->nodes;
        elems    =    ibMesh->elems;

        // allocate memory for the half edge data structure
        PetscMalloc(nodes * sizeof(Vertex), &(ibMesh->vertices));
        PetscMalloc(3 * elems * sizeof(HalfEdge), &(ibMesh->halfEdges));
        PetscMalloc(elems * sizeof(Face), &(ibMesh->faces));

        for(PetscInt i = 0; i < nodes; i++)
        {
            // set the half edge structure vertices
            ibMesh->vertices[i].nCoor = nSet(ibMesh->nCoor[i]);
            ibMesh->vertices[i].edge  = NULL;
            ibMesh->vertices[i].vertexId = i;
        }

        for(PetscInt i = 0; i < elems; i++)
        {
            PetscInt id[3];
            id[0] = ibMesh->nID1[i];
            id[1] = ibMesh->nID2[i];
            id[2] = ibMesh->nID3[i];

            // Create half-edges for the current element
            for (PetscInt ed = 0; ed < 3; ed++) {
                PetscInt v1 = id[ed];
                PetscInt v2 = id[(ed + 1) % 3];

                ibMesh->halfEdges[3 * i + ed].origin = &ibMesh->vertices[v1];
                ibMesh->halfEdges[3 * i + ed].next = &ibMesh->halfEdges[3 * i + (ed + 1) % 3];
                ibMesh->halfEdges[3 * i + ed].twin = NULL;
                ibMesh->halfEdges[3 * i + ed].face = &ibMesh->faces[i];

                ibMesh->vertices[v1].edge = &ibMesh->halfEdges[3 * i + ed];
            }
        }

        // Form faces
        for (int i = 0; i < elems; i++) {
            ibMesh->faces[i].edge = &ibMesh->halfEdges[3 * i]; // Link the face to one of its half-edges
            ibMesh->faces[i].faceId = i;
        }

        // Link twins of the half-edges
        for (PetscInt i = 0; i < 3 * elems; i++) {
            Vertex *startVertex = ibMesh->halfEdges[i].origin;
            Vertex *endVertex = ibMesh->halfEdges[i].next->origin;

            // Search for the twin among the remaining half-edges
            PetscInt vertexId = startVertex->vertexId;

            // elements connected to startVertex
            PetscInt numElemConnected = ibMesh->ibmMeshNode[vertexId].numConnected;

            for (PetscInt j = 0; j < numElemConnected; j++)
            {
                PetscInt ele = ibMesh->ibmMeshNode[vertexId].elem[j];

                //loop through the half edges of the element
                HalfEdge  *start_he = ibMesh->faces[ele].edge;
                HalfEdge  *he = start_he;
                do
                {
                    if (he->origin == endVertex && he->next->origin == startVertex)
                    {
                        ibMesh->halfEdges[i].twin = he;
                        break;
                    }
                    he = he->next;
                }
                while(he !=start_he);

                if(ibMesh->halfEdges[i].twin!=NULL)
                {
                    break;
                }
            }

        }

        // for (PetscInt i = 0; i < 3 * elems; i++)
        // {
        //     printf("Half-edge %ld: Start vertex: (%lf, %lf, %lf)\n", i, ibMesh->halfEdges[i].origin->nCoor.x, ibMesh->halfEdges[i].origin->nCoor.y, ibMesh->halfEdges[i].origin->nCoor.z);
        //     printf("Half-edge %ld: second vertex: (%lf, %lf, %lf)\n", i, ibMesh->halfEdges[i].next->origin->nCoor.x, ibMesh->halfEdges[i].next->origin->nCoor.y, ibMesh->halfEdges[i].next->origin->nCoor.z);
        //     printf("Half-edge %ld: third vertex: (%lf, %lf, %lf)\n\n", i, ibMesh->halfEdges[i].next->next->origin->nCoor.x, ibMesh->halfEdges[i].next->next->origin->nCoor.y, ibMesh->halfEdges[i].next->next->origin->nCoor.z);
        //
        //     if(ibMesh->halfEdges[i].twin!=NULL)
        //     {
        //         printf("Half-edge twin of %ld: Start vertex: (%lf, %lf, %lf)\n", i, ibMesh->halfEdges[i].twin->origin->nCoor.x, ibMesh->halfEdges[i].twin->origin->nCoor.y, ibMesh->halfEdges[i].twin->origin->nCoor.z);
        //         printf("Half-edge twin of %ld: Start vertex: (%lf, %lf, %lf)\n\n", i, ibMesh->halfEdges[i].twin->next->origin->nCoor.x, ibMesh->halfEdges[i].twin->next->origin->nCoor.y, ibMesh->halfEdges[i].twin->next->origin->nCoor.z);
        //     }
        //     else
        //     {
        //         printf("Half-edge twin of %ld: NULL\n\n", i );
        //     }
        //
        // }


        for (PetscInt i = 0; i < 3 * elems; i++)
        {
            if(ibMesh->halfEdges[i].twin==NULL)
            {
                char error[512];
                sprintf(error, "IBM body: %s is not a closed body. If you are sure it is, then the issue could be that the element normals are not pointing outwards. Enable checkNormal in IBMProperties\n", ibm->ibmBody[b]->bodyName.c_str());
                fatalErrorInFunction("createHalfEdgeDataStructure",  error);
            }
        }

        PetscPrintf(PETSC_COMM_WORLD, "     Created half edge data structure for IBM body:%s \n", ibm->ibmBody[b]->bodyName.c_str());
        PetscPrintf(PETSC_COMM_WORLD, "     Num nodes: %ld, Num elements: %ld\n\n", nodes, elems);

    }
    return 0;
}
//***************************************************************************************************************//

PetscErrorCode computeIBMElementNormal(ibm_ *ibm)
{
    PetscInt      n1, n2, n3;                                     // nodes of a particular IBM element
    PetscInt      e, b;
    Cmpnts        vec1, vec2, temp, refPt, offsetVec;
    cellIds       sCell;
    PetscReal     normMag, minBound, offset;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {
        boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
        ibmMesh       *ibMesh = ibm->ibmBody[b]->ibMsh;                         // pointer to the ibm body mesh
        searchBox     *sBox = &(ibm->sBox[b]);
        list          *searchCellList = ibm->ibmBody[b]->searchCellList;

        //set initial element normal and element center
        for (e=0; e<ibMesh->elems; e++)
        {
            // get the element nodes
            n1 = ibMesh->nID1[e]; n2 = ibMesh->nID2[e]; n3 = ibMesh->nID3[e];

            vec1 = nSub(ibMesh->nCoor[n2], ibMesh->nCoor[n1]);

            vec2 = nSub(ibMesh->nCoor[n3], ibMesh->nCoor[n1]);

            // normal to the face is found as cross product of the edges vec1 and vec2
            ibMesh->eN[e] = nCross(vec1, vec2);
            normMag = nMag(ibMesh->eN[e]);
            mScale(1.0/normMag, ibMesh->eN[e]);

            //element center
            temp = nSum(ibMesh->nCoor[n1], ibMesh->nCoor[n2]);
            ibMesh->eCent[e] = nSum( temp, ibMesh->nCoor[n3]);
            mScale(1/3.0, ibMesh->eCent[e]);

            //element area
            ibMesh->eA[e] = normMag/2.0;

        }

        if(ibm->checkNormal)
        {
            PetscPrintf(PETSC_COMM_WORLD, "     Checking IBM element normal direction for body: %s...", ibm->ibmBody[b]->bodyName.c_str());
            // set offset distance
            minBound = PetscMin( PetscMin(ibBox->Lx, ibBox->Ly), ibBox->Lz);
            offset   = 1.0e-7;

            // check that the normal points outwards
            for (e=0; e<ibMesh->elems; e++)
            {
                //move reference distance from the element center
                offsetVec = nScale(offset, ibMesh->eN[e]);
                refPt = nSum(ibMesh->eCent[e], offsetVec);

                // find the search cell were the ref pt is located
                sCell.i = floor((refPt.x - ibBox->xmin) / sBox->dcx);
                sCell.j = floor((refPt.y - ibBox->ymin) / sBox->dcy);
                sCell.k = floor((refPt.z - ibBox->zmin) / sBox->dcz);

                //perform raycasting test to check if the point is inside or outside the body
                PetscReal val;

                val = rayCastingTest(refPt, ibMesh, sCell, sBox, ibBox, searchCellList);

                if (val > 0.1)
                {

                    // get the element nodes
                    n1 = ibMesh->nID1[e]; n2 = ibMesh->nID2[e]; n3 = ibMesh->nID3[e];

                    //reverse the element node order
                    ibMesh->nID2[e] = n3; ibMesh->nID3[e] = n2;

                    //save the normal with the correct orientation
                    mScale(-1.0, ibMesh->eN[e]);
                }

            }

            PetscPrintf(PETSC_COMM_WORLD, "done\n");
        }


        //use nearby elements to average the normal direction. this requires the node to element connectivity
        if(ibm->averageNormal)
        {
            PetscPrintf(PETSC_COMM_WORLD, "Averaging IBM element normals...");
            ibmNode *ibmNodes  = ibMesh->ibmMeshNode;

            //temporary array to store the averaged element normals
            Cmpnts *tempEn = new Cmpnts[ibMesh->elems];

            for (e=0; e<ibMesh->elems; e++)
            {
                PetscInt i, j;
                // get the element nodes
                n1 = ibMesh->nID1[e]; n2 = ibMesh->nID2[e]; n3 = ibMesh->nID3[e];

                //find the number of elements connected to each node
                PetscInt numElemn1 = ibmNodes[n1].numConnected;
                PetscInt numElemn2 = ibmNodes[n2].numConnected;
                PetscInt numElemn3 = ibmNodes[n3].numConnected;

                PetscInt mergedElemArray[numElemn1 + numElemn2 + numElemn3];

                //copy the first node elements into mergedElemArray
                for (i = 0; i < numElemn1; i++)
                {
                    mergedElemArray[i] = ibmNodes[n1].elem[i];

                    if(e == 100)
                    {
                        PetscPrintf(PETSC_COMM_WORLD, "NODE 1 elements = %ld\n", ibmNodes[n1].elem[i]);
                    }
                }

                //copy the second node elements into merged array
                for (j = 0; j < numElemn2; j++)
                {
                    if (!isPresent(mergedElemArray, i+1, ibmNodes[n2].elem[j]))
                    {
                        mergedElemArray[i++] = ibmNodes[n2].elem[j];
                    }

                    if(e == 100)
                    {
                        PetscPrintf(PETSC_COMM_WORLD, "NODE 2 elements = %ld\n", ibmNodes[n2].elem[j]);
                    }
                }

                //copy the third node elements into merged array
                for (j = 0; j < numElemn3; j++)
                {
                    if (!isPresent(mergedElemArray, i+1, ibmNodes[n3].elem[j]))
                    {
                        mergedElemArray[i++] = ibmNodes[n3].elem[j];
                    }

                    if(e == 100)
                    {
                        PetscPrintf(PETSC_COMM_WORLD, "NODE 3 elements = %ld\n", ibmNodes[n3].elem[j]);
                    }
                }

                //average the normal at element e from its surrounding elements

                PetscReal sumEnx = 0.0, sumEny = 0.0, sumEnz = 0.0;

                for (j = 0; j < i; j++)
                {
                    sumEnx += ibMesh->eN[mergedElemArray[j]].x;
                    sumEny += ibMesh->eN[mergedElemArray[j]].y;
                    sumEnz += ibMesh->eN[mergedElemArray[j]].z;

                    if(e == 100)
                    {
                        PetscPrintf(PETSC_COMM_WORLD, "merged elements = %ld\n", mergedElemArray[j]);
                    }
                }

                // save the averaged element normals
                tempEn[e] = nSetFromComponents(sumEnx, sumEny, sumEnz);

                PetscReal normMag = nMag(tempEn[e]);
                mScale(1.0/normMag, tempEn[e]);

                if(e == 100)
                {
                    PetscPrintf(PETSC_COMM_WORLD, "normal = %lf %lf %lf, avg normal = %lf %lf %lf\n", ibMesh->eN[e].x, ibMesh->eN[e].y, ibMesh->eN[e].z, tempEn[e].x, tempEn[e].y, tempEn[e].z);
                }
            }

            for (e=0; e<ibMesh->elems; e++)
            {
                ibMesh->eN[e] = nSet(tempEn[e]);
            }

            delete[] tempEn;
            PetscPrintf(PETSC_COMM_WORLD, "done\n");
        }

        //find tangential unit vectors for the ibm mesh elements
        for (e=0; e<ibMesh->elems; e++)
        {

            // tangential to the face( eT1 and eT2)
            // eT1 = eN x k
            if (
                (((1.0 - ibMesh->eN[e].z ) <= 1e-6 ) && ((-1.0 + ibMesh->eN[e].z ) < 1e-6))
                ||
                (((ibMesh->eN[e].z + 1.0 ) <= 1e-6 ) && ((-1.0 - ibMesh->eN[e].z ) < 1e-6))
            )
            {
                ibMesh->eT1[e].x = 1.0;
                ibMesh->eT1[e].y = 0.0;
                ibMesh->eT1[e].z = 0.0;

                ibMesh->eT2[e].x = 0.0;
                ibMesh->eT2[e].y = 1.0;
                ibMesh->eT2[e].z = 0.0;
            }
            else
            {
                ibMesh->eT1[e].x =  ibMesh->eN[e].y/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT1[e].y = -ibMesh->eN[e].x/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT1[e].z = 0 ;

                 // eT2 = eT1 x eN
                ibMesh->eT2[e].x = -ibMesh->eN[e].x*ibMesh->eN[e].z/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT2[e].y = -ibMesh->eN[e].y*ibMesh->eN[e].z/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT2[e].z = sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
            }

        }

        if(ibm->writeSTL)
        {
            {
                writeSTLFile(ibm, b);
                MPI_Barrier(ibm->access->mesh->MESH_COMM);
            }
        }

    }
    return (0);
}
//***************************************************************************************************************//

PetscErrorCode findIBMControlledProcs(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;

    DMDALocalInfo info = mesh->info;
    DM            da = mesh->da, fda = mesh->fda, sda = mesh->sda;

    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscInt       i, j, k, b;
    PetscInt       lxs, lxe, lys, lye, lzs, lze;

    Vec            Coor;
    Cmpnts         ***coor;

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    PetscMPIInt    rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        // define communicator color
        PetscInt    commColor = MPI_UNDEFINED;
        ibmObject   *ibmBody  = ibm->ibmBody[b];

        if(ibm->dynamic)
        {
            PetscReal xmin_b = ibmBody->procBoundCenter.x - 0.5 * ibmBody->procBoundSize.x,
                      xmax_b = ibmBody->procBoundCenter.x + 0.5 * ibmBody->procBoundSize.x,
                      ymin_b = ibmBody->procBoundCenter.y - 0.5 * ibmBody->procBoundSize.y,
                      ymax_b = ibmBody->procBoundCenter.y + 0.5 * ibmBody->procBoundSize.y,
                      zmin_b = ibmBody->procBoundCenter.z - 0.5 * ibmBody->procBoundSize.z,
                      zmax_b = ibmBody->procBoundCenter.z + 0.5 * ibmBody->procBoundSize.z;

            if
            (
                (
                  coor[zs ][ys ][xs ].x > xmax_b || coor[lze-1][lye-1][lxe-1].x < xmin_b
                ) ||
                (
                  coor[zs ][ys ][xs ].y > ymax_b || coor[lze-1][lye-1][lxe-1].y < ymin_b
                ) ||
                (
                  coor[zs ][ys ][xs ].z > zmax_b || coor[lze-1][lye-1][lxe-1].z < zmin_b
                )
            )
            {
                ibmBody->ibmControlled = 0;
            }
            else
            {
                // set this ibm body to be controlled by this processor
                ibmBody->ibmControlled = 1;

                // set the communicator color
                commColor = 1;
            }

        }
        else
        {
            if(ibm->numIBMFluid!=0)
            {
                // set this ibm body to be controlled by this processor
                ibmBody->ibmControlled = 1;

                // set the communicator color
                commColor = 1;
            }
            else
            {
                ibmBody->ibmControlled = 0;
            }
        }

        // create communicator
        MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(ibmBody->IBM_COMM));

        // the master rank of this IBM_COMM communicator will write this IBM I/O file
        PetscMPIInt thisIBMRank = 10;
        MPI_Comm_rank(ibmBody->IBM_COMM, &thisIBMRank);

        if(ibmBody->ibmControlled == 1 && ibm->dbg)
        {
            PetscPrintf(PETSC_COMM_SELF,"body %ld, controlling processor - global rank %ld (local rank %ld)\n", b, rank, thisIBMRank);
        }

    }

    DMDAVecRestoreArray(fda, Coor, &coor);

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode initElementProjectionProcs(ibm_ *ibm)
{
    mesh_              *mesh = ibm->access->mesh;
    DM                 da = mesh->da, fda = mesh->fda;

    PetscInt           i, j, k, c, b, e;
    ibmFluidCell       *ibF = ibm->ibmFCells;
    Cmpnts             eN;                                     // unit normal from the closest IBM mesh element of the IBM fluid cell
    Cmpnts             eC;
    Cmpnts             checkPt;                                // point where the pressure or shear stress needs to be interpolated
    PetscReal          refLength;                              // reference length set as the average cell size
    PetscReal          lminDist, gminDist;
    Cmpnts             perturb;

    Cmpnts             ***cent;

    PetscMPIInt        nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt        rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // max perturbation amplitude
    PetscReal          maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal          procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody  = ibm->ibmBody[b];
        ibmMesh     *ibMsh    = ibmBody->ibMsh;                         // pointer to the ibm body mesh

        //check if process controls this ibm body
        if(ibmBody->ibmControlled)
        {
            // number of points in the AD mesh
            PetscInt nptsIBM = ibMsh->elems;

            // create temporary vectors
            std::vector<PetscReal> lminDist(nptsIBM);
            std::vector<PetscReal> gminDist(nptsIBM);
            std::vector<Cmpnts>    perturb(nptsIBM);

            //loop through the ibm mesh elements
            for(e = 0; e < ibMsh->elems; e++)
            {
                //element normal
                eN = ibMsh->eN[e];

                //element center
                eC = ibMsh->eCent[e];

                //reference length - square root of the element area vector magnitude
                refLength = pow( 2.0 * ibMsh->eA[e], 1./2.);

                //initial point for interpolation
                checkPt = nSum(nScale(refLength, eN), eC);

                // initialize min dists to a big value
                lminDist[e] = 1e20;
                gminDist[e] = 1e20;

                // set point perturbation
                perturb[e].x  = procContrib;
                perturb[e].y  = procContrib;
                perturb[e].z  = procContrib;

                // perturb the point position
                mSum(checkPt, perturb[e]);

                // find the closest cell center
                PetscReal  distMinMag = 1e20;
                cellIds    closestCell;

                //loop through all the cells
                for(c = 0; c < ibm->numIBMFluid; c++)
                {
                    //cell indices
                    i = ibF[c].cellId.i;
                    j = ibF[c].cellId.j;
                    k = ibF[c].cellId.k;

                    // compute distance vector from check point to ibm fluid cell
                    Cmpnts  distVec = nSub(checkPt, cent[k][j][i]);

                    // compute magnitude
                    PetscReal dist  = nMag(distVec);

                    if(dist < distMinMag)
                    {
                        distMinMag = dist;
                        closestCell.i = i;
                        closestCell.j = j;
                        closestCell.k = k;
                    }
                }

                // save closest cell indices - local array
                ibmBody->closestCells[e].i = closestCell.i;
                ibmBody->closestCells[e].j = closestCell.j;
                ibmBody->closestCells[e].k = closestCell.k;

                // save min dist
                lminDist[e] = distMinMag;

            }

            // find the min distance among the ibm communicator processors
            MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), nptsIBM, MPIU_REAL, MPIU_MIN, ibmBody->IBM_COMM);

            PetscInt lsum = 0, gsum = 0;
            for(e = 0; e < ibMsh->elems; e++)
            {
                // point is controlled
                if(lminDist[e] == gminDist[e])
                {
                    ibmBody->thisPtControlled[e] = 1;

                    if(ibm->dbg)
                    {
                        lsum++;
                    }
                }
                // point is not controlled
                else
                {
                    ibmBody->thisPtControlled[e] = 0;
                }

                //initialize flag for ibm element processor flag to 0
                ibmBody->thisPtControlTransfer[e] = 0;
            }

            //for debug purpose
            if(ibm->dbg)
            {
                // the master rank of this IBM_COMM communicator will write this IBM I/O file
                PetscMPIInt thisIBMRank = 10;
                MPI_Comm_rank(ibmBody->IBM_COMM, &thisIBMRank);

                PetscPrintf(PETSC_COMM_SELF,"body %ld, controlling processor - global rank %ld (local rank %ld), number of elements = %ld\n", b, rank, thisIBMRank, lsum);
                MPI_Allreduce(&lsum, &gsum, 1, MPIU_INT, MPIU_SUM, ibmBody->IBM_COMM);
                PetscPrintf(ibmBody->IBM_COMM, "body %ld, total number of elements = %ld\n", b, gsum);

            }

            // clean memory
            std::vector<PetscReal> ().swap(lminDist);
            std::vector<PetscReal> ().swap(gminDist);
            std::vector<Cmpnts> ().swap(perturb);
        }

    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode IBMProjectionProcessorTransfer(ibm_ *ibm)
{
    mesh_            *mesh = ibm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, c, b, e;

    Cmpnts           eN;                                     // unit normal from the closest IBM mesh element of the IBM fluid cell
    Cmpnts           eC;
    Cmpnts           checkPt;                                // point where the pressure or shear stress needs to be interpolated
    PetscReal        refLength;

    Cmpnts           ***cent;
    PetscReal        ***nvert;

    PetscMPIInt        nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt        rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        //pointer to access the ibm body
        ibmObject   *ibmBody  = ibm->ibmBody[b];
        ibmMesh     *ibMsh    = ibmBody->ibMsh;                         // pointer to the ibm body mesh

        //check if process controls this ibm body
        if(ibmBody->ibmControlled)
        {
            // number of points in the AD mesh
            PetscInt nptsIBM = ibMsh->elems;

            // create temporary vectors
            std::vector<PetscInt> lptcontrolTransfer(nptsIBM);
            std::vector<cellIds>  lClosestCell(nptsIBM);

            //loop through the ibm mesh elements
            for(e = 0; e < ibMsh->elems; e++)
            {
                //this processor previosly controlled this ibm element
                if(ibmBody->thisPtControlled[e])
                {
                    /*find the new projection point of this element normal
                    The search for the ibm fluid cell that is closest to this
                    projection point is done in the neighbourhood of the old
                    closest point */

                    //element normal
                    eN = ibMsh->eN[e];

                    //element center
                    eC = ibMsh->eCent[e];

                    //reference length - square root of the element area vector magnitude
                    refLength = pow( 2.0 * ibMsh->eA[e], 1./2.);

                    //initial point for interpolation
                    checkPt = nSum(nScale(refLength, eN), eC);

                    //current closest ibm fluid cell to point
                    PetscInt    ci, cj, ck;
                    cellIds     currentClosest = ibmBody->closestCells[e];

                    ci = currentClosest.i;
                    cj = currentClosest.j;
                    ck = currentClosest.k;

                    //find the closest neighbour of ck,cj,ci to checkPt
                    PetscInt    i1, j1, k1;
                    PetscReal   dmin = 10e6, d;
                    PetscInt    ic, jc, kc;

                    for (k1=ck-2; k1<ck+3; k1++)
                    for (j1=cj-2; j1<cj+3; j1++)
                    for (i1=ci-2; i1<ci+3; i1++)
                    {
                        if
                        (
                            (
                                k1>=1 && k1<mz-1 &&
                                j1>=1 && j1<my-1 &&
                                i1>=1 && i1<mx-1
                            ) && (!isIBMSolidCell(k1, j1, i1, nvert))
                        )
                        {
                            d = pow((checkPt.x - cent[k1][j1][i1].x), 2) +
                                pow((checkPt.y - cent[k1][j1][i1].y), 2) +
                                pow((checkPt.z - cent[k1][j1][i1].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                ic = i1;
                                jc = j1;
                                kc = k1;
                            }
                        }
                    }

                    //update the closest cell id
                    lClosestCell[e].i = ic;
                    lClosestCell[e].j = jc;
                    lClosestCell[e].k = kc;

                    //set pt controlled if within the current processor limits
                    // else transfer

                    if (kc >= lzs && kc < lze && jc >= lys && jc < lye && ic >= lxs && ic < lxe)
                    {
                        //new closest cell is within the current processor
                        ibmBody->thisPtControlled[e] = 1;

                        //set flag for processor transfer to 0
                        lptcontrolTransfer[e] = 0;

                    }
                    else
                    {
                        //new closest cell belongs to the adjacent processor
                        ibmBody->thisPtControlled[e] = 0;

                        lptcontrolTransfer[e] = 1;

                    }

                }
                else
                {
                    lptcontrolTransfer[e] = 0;

                    lClosestCell[e].i = 0;
                    lClosestCell[e].j = 0;
                    lClosestCell[e].k = 0;
                }

            }

            //scatter the pt control transfer flag and new closest element
            MPI_Allreduce(&(lptcontrolTransfer[0]), &(ibmBody->thisPtControlTransfer[0]), nptsIBM, MPIU_INT, MPI_SUM, ibmBody->IBM_COMM);
            MPI_Allreduce(&(lClosestCell[0]), &(ibmBody->closestCells[0]), 3*nptsIBM, MPIU_INT, MPI_SUM, ibmBody->IBM_COMM);


            //clean the temporary vectors
            std::vector<PetscInt> ().swap(lptcontrolTransfer);
            std::vector<cellIds> ().swap(lClosestCell);

            //loop through the ibm mesh elements and transfer the ibm element control
            for(e = 0; e < ibMsh->elems; e++)
            {
                // now only check the elements which will be transferred
                if(ibmBody->thisPtControlTransfer[e])
                {
                    cellIds Id = ibmBody->closestCells[e];

                    if (Id.k >= lzs && Id.k < lze && Id.j >= lys && Id.j < lye && Id.i >= lxs && Id.i < lxe)
                    {
                        ibmBody->thisPtControlled[e] = 1;

                    }
                }

            }


            if(ibm->dbg)
            {
                MPI_Barrier(ibmBody->IBM_COMM);

                PetscInt lsum = 0, gsum = 0;
                for(e = 0; e < ibMsh->elems; e++)
                {
                    //for debug purpose
                    if(ibmBody->thisPtControlled[e] == 1 && ibm->dbg)
                    {
                        lsum++;
                    }
                }

                // the master rank of this IBM_COMM communicator will write this IBM I/O file
                PetscMPIInt thisIBMRank = 10;
                MPI_Comm_rank(ibmBody->IBM_COMM, &thisIBMRank);

                PetscPrintf(PETSC_COMM_SELF,"body %ld, controlling processor - global rank %ld (local rank %ld), number of elements = %ld\n", b, rank, thisIBMRank, lsum);
                MPI_Allreduce(&lsum, &gsum, 1, MPIU_INT, MPIU_SUM, ibmBody->IBM_COMM);
                PetscPrintf(ibmBody->IBM_COMM, "body %ld, total number of elements = %ld\n", b, gsum);

            }

        }

    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode createProcessorBufferZones(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    PetscInt      b;

    DMDALocalInfo info = mesh->info;
    DM            da = mesh->da, fda = mesh->fda, sda = mesh->sda;

    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscInt       lxs, lxe, lys, lye, lzs, lze;
    PetscReal      buffer = 0.0;
    Vec            Coor;
    Cmpnts         ***coor;

    PetscMPIInt   rank, size;
    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody  = ibm->ibmBody[b];
        elementBox  *eBox = ibmBody->eBox;

        if(ibmBody->ibmControlled)
        {
            //set the inner buffer zone - 2 cell widths from the cell boundary
            eBox->innerZone.xmin = coor[lzs+1 ][lys ][lxs ].x + buffer;
            eBox->innerZone.xmax = coor[lze-3 ][lye-1 ][lxe-1 ].x - buffer;
            eBox->innerZone.ymin = coor[lzs ][lys ][lxs+1 ].y + buffer;
            eBox->innerZone.ymax = coor[lze-1 ][lye-1 ][lxe-3 ].y - buffer;
            eBox->innerZone.zmin = coor[lzs ][lys+1 ][lxs ].z + buffer;
            eBox->innerZone.zmax = coor[lze-1 ][lye-3 ][lxe-1 ].z - buffer;

            eBox->innerZone.Lx = eBox->innerZone.xmax - eBox->innerZone.xmin;
            eBox->innerZone.Ly = eBox->innerZone.ymax - eBox->innerZone.ymin;
            eBox->innerZone.Lz = eBox->innerZone.zmax - eBox->innerZone.zmin;

            if(eBox->innerZone.Lx <= 0)
            {
                char error[512];
                sprintf(error, "processor domain size in coordinate x direction is 4 cells or less. Cannot create inner buffer zone for ibm elements\n");
                fatalErrorInFunction("createProcessorBufferZones",  error);
            }

            if(eBox->innerZone.Ly <= 0)
            {
                char error[512];
                sprintf(error, "processor domain size in coordinate y direction is 4 cells or less. Cannot create inner buffer zone for ibm elements\n");
                fatalErrorInFunction("createProcessorBufferZones",  error);
            }

            if(eBox->innerZone.Lz <= 0)
            {
                char error[512];
                sprintf(error, "processor domain size in coordinate z direction is 4 cells or less. Cannot create inner buffer zone for ibm elements\n");
                fatalErrorInFunction("createProcessorBufferZones",  error);
            }

            //set the outer buffer zone
            //min values
            if(lzs == 1)
            {
                eBox->outerZone.xmin = coor[lzs-1 ][lys ][lxs].x - buffer;
            }
            else
            {
                eBox->outerZone.xmin = coor[lzs-3 ][ys ][xs].x - buffer;
            }

            if(lxs == 1)
            {
                eBox->outerZone.ymin = coor[lzs ][lys ][lxs-1].y - buffer;
            }
            else
            {
                eBox->outerZone.ymin = coor[lzs ][lys ][lxs-3].y - buffer;
            }

            if(lys == 1)
            {
                eBox->outerZone.zmin = coor[lzs ][lys-1 ][lxs].z - buffer;
            }
            else
            {
                eBox->outerZone.zmin = coor[lzs ][lys-3 ][lxs].z - buffer;
            }

            //max values
            if(lze == mz-1)
            {
                eBox->outerZone.xmax = coor[lze-1 ][lye-1 ][lxe-1].x + buffer;
            }
            else
            {
                eBox->outerZone.xmax = coor[lze+1 ][lye-1 ][lxe-1].x + buffer;
            }

            if(lxe == mx-1)
            {
                eBox->outerZone.ymax = coor[lze-1 ][lye-1 ][lxe-1].y + buffer;
            }
            else
            {
                eBox->outerZone.ymax = coor[lze-1 ][lye-1 ][lxe+1].y + buffer;
            }

            if(lye == my-1)
            {
                eBox->outerZone.zmax = coor[lze-1 ][lye-1 ][lxe-1].z + buffer;
            }
            else
            {
                eBox->outerZone.zmax = coor[lze-1 ][lye+1 ][lxe-1].z + buffer;
            }

            eBox->outerZone.Lx = eBox->outerZone.xmax - eBox->outerZone.xmin;
            eBox->outerZone.Ly = eBox->outerZone.ymax - eBox->outerZone.ymin;
            eBox->outerZone.Lz = eBox->outerZone.zmax - eBox->outerZone.zmin;

            if(ibm->dbg)
            {
                PetscPrintf(PETSC_COMM_SELF, "Rank = %d, inner box dimensions:\n   xmin = %lf, xmax = %lf\n   ymin = %lf, ymax = %lf\n   zmin = %lf, zmax = %lf\n   Lx = %lf, Ly = %lf, Lz = %lf\n", rank, eBox->innerZone.xmin, eBox->innerZone.xmax, eBox->innerZone.ymin,
                                eBox->innerZone.ymax, eBox->innerZone.zmin, eBox->innerZone.zmax, eBox->innerZone.Lx, eBox->innerZone.Ly, eBox->innerZone.Lz);

                PetscPrintf(PETSC_COMM_SELF, "Rank = %d, outer box dimensions:\n   xmin = %lf, xmax = %lf\n   ymin = %lf, ymax = %lf\n   zmin = %lf, zmax = %lf\n   Lx = %lf, Ly = %lf, Lz = %lf\n", rank, eBox->outerZone.xmin, eBox->outerZone.xmax, eBox->outerZone.ymin,
                                                eBox->outerZone.ymax, eBox->outerZone.zmin, eBox->outerZone.zmax, eBox->outerZone.Lx, eBox->outerZone.Ly, eBox->outerZone.Lz);

            }
        }
        else
        {
            //initialize the bounding box to extreme values
            eBox->outerZone.xmin =  1.e23;
            eBox->outerZone.xmax = -1.e23;
            eBox->outerZone.ymin =  1.e23;
            eBox->outerZone.ymax = -1.e23;
            eBox->outerZone.zmin =  1.e23;
            eBox->outerZone.zmax = -1.e23;

            eBox->outerZone.Lx   =  0.;
            eBox->outerZone.Ly   =  0.;
            eBox->outerZone.Lz   =  0.;

        }
    }

    DMDAVecRestoreArray(fda, Coor, &coor);

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode initElementProcs(ibm_ *ibm)
{
    Cmpnts             p1, p2, p3;
    PetscInt           b, e;
    mesh_              *mesh  = ibm->access->mesh;

    //create the processor inner and outer buffer zones
    createProcessorBufferZones(ibm);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {

        ibmObject   *ibmBody  = ibm->ibmBody[b];
        ibmMesh     *ibMsh    = ibmBody->ibMsh;                         // pointer to the ibm body mesh
        PetscInt    *nv1      = ibMsh->nID1,
                    *nv2      = ibMsh->nID2,
                    *nv3      = ibMsh->nID3;
        Cmpnts      *nCoor    = ibMsh->nCoor;

        PetscMPIInt   ibmRank, ibmProcSize;
        MPI_Comm_size(ibmBody->IBM_COMM, &ibmProcSize);
        MPI_Comm_rank(ibmBody->IBM_COMM, &ibmRank);

        //check if processor controls this ibm body
        if(ibmBody->ibmControlled)
        {
            elementBox  *eBox = ibmBody->eBox;
            PetscInt    sum = 0;

            //loop through the ibm mesh elements
            for(e = 0; e < ibMsh->elems; e++)
            {
                //element node coordinates
                p1 = nCoor[nv1[e]];
                p2 = nCoor[nv2[e]];
                p3 = nCoor[nv3[e]];

                //check if the element node coordinates is in the outer buffer zone of the processor
                if(    isInsideBoundingBox(p1, eBox->outerZone)
                    || isInsideBoundingBox(p2, eBox->outerZone)
                    || isInsideBoundingBox(p3, eBox->outerZone)
                )
                {
                    eBox->thisElemControlled[e] = 1;
                }
                else
                {
                    eBox->thisElemControlled[e] = 0;
                }

                //initialize flag for ibm element transfer to 0
                eBox->thisElemTransfered[e] = 0;

                if(ibm->dbg && eBox->thisElemControlled[e] == 1)
                {
                    sum++;
                }
            }

            if(ibm->dbg)
            {
                PetscPrintf(PETSC_COMM_SELF, "rank = %d, number of ibm elements = %ld\n", ibmRank, sum);
            }

            if(ibmRank == 0)
            {
                PetscPrintf(PETSC_COMM_SELF, "IBM body: %s, num of controlling processors = %d\n", ibmBody->bodyName.c_str(), ibmProcSize);
            }
        }

    }

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

/*the elements of the ibm mesh are divided among the ibm processors based on the location of the element center.
 2 buffer zones are created around the processor boundary - the inner zone and the outer zone
 these buffer zones determine if the processor controls this element, if needs to be transferred to another processor
 or removed from the current processor.
 Based on this, an element can belong to multiple processors at the same time.

 Based on the element location w.r.t processor boundary and the buffer zone, the processor can decide to:
 a) Control this element - use it in its calculations
 b) Transfer it to another processor - transfer does not mean elimination from the current processor
 c) Delete it from the current processor
*/

PetscErrorCode IBMElementProcessorTransfer(ibm_ *ibm)
{
    mesh_            *mesh = ibm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         c, b, e;

    Cmpnts           p1, p2, p3;

    Cmpnts           ***cent;
    PetscReal        ***nvert;

    PetscMPIInt      globalRank, ProcSize;
    MPI_Comm_size(mesh->MESH_COMM, &ProcSize);
    MPI_Comm_rank(mesh->MESH_COMM, &globalRank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody  = ibm->ibmBody[b];
        ibmMesh     *ibMsh    = ibmBody->ibMsh;
        PetscInt    *nv1      = ibMsh->nID1,
                    *nv2      = ibMsh->nID2,
                    *nv3      = ibMsh->nID3;
        Cmpnts      *nCoor    = ibMsh->nCoor;


        //check if process controls this ibm body
        if(ibmBody->ibmControlled)
        {
            elementBox  *eBox     = ibmBody->eBox;

            //find the current ibm processor rank and total ibm processors
            PetscMPIInt   ibmRank, ibmProcSize;
            MPI_Comm_size(ibmBody->IBM_COMM, &ibmProcSize);
            MPI_Comm_rank(ibmBody->IBM_COMM, &ibmRank);

            // create temporary vectors
            std::vector <PetscInt>   lNumTransfers(ibmProcSize);
            std::vector <PetscInt>   gNumTransfers(ibmProcSize);

            std::vector <PetscInt>   lTransferElem;
            std::vector <PetscInt>   gTransferElem;

            // total number of transfer elements, start index within the global element transfer array
            PetscInt   totalNumTransfers = 0, startInd = 0;

            //initialize the temporary arrays
            for (c=0; c < ibmProcSize; c++){
                lNumTransfers[c] = 0;
                gNumTransfers[c] = 0;
            }

            //loop through the ibm mesh elements
            for(e = 0; e < ibMsh->elems; e++)
            {
                //reset to 0
                eBox->thisElemTransfered[e] = 0;

                //this processor previosly controlled this ibm element
                if(eBox->thisElemControlled[e])
                {
                    //element node coordinates
                    p1 = nCoor[nv1[e]];
                    p2 = nCoor[nv2[e]];
                    p3 = nCoor[nv3[e]];

                    //check the element location w.r.t processor buffer zones
                    if(    isInsideBoundingBox(p1, eBox->outerZone)
                        || isInsideBoundingBox(p2, eBox->outerZone)
                        || isInsideBoundingBox(p3, eBox->outerZone)
                    )
                    {
                        eBox->thisElemControlled[e] = 1;
                        eBox->thisElemTransfered[e] = 1;

                        // if fully inside the inner buffer zone dont transfer, else transfer
                        if(    isInsideBoundingBox(p1, eBox->innerZone)
                            && isInsideBoundingBox(p2, eBox->innerZone)
                            && isInsideBoundingBox(p3, eBox->innerZone)
                        )
                        {
                            eBox->thisElemTransfered[e] = 0;
                        }

                    }
                    else
                    {
                        eBox->thisElemControlled[e] = 0;
                        eBox->thisElemTransfered[e] = 1;
                    }

                    if(eBox->thisElemTransfered[e])
                    {
                        lNumTransfers[ibmRank]++;
                    }

                }

            }

            //scatter the number of transfers for each processor
            MPI_Allreduce(&(lNumTransfers[0]), &(gNumTransfers[0]), ibmProcSize, MPIU_INT, MPI_SUM, ibmBody->IBM_COMM);

            //find total number of transfer elements
            for (c=0; c < ibmProcSize; c++){
                totalNumTransfers += gNumTransfers[c];
            }

            lTransferElem.resize(totalNumTransfers);
            gTransferElem.resize(totalNumTransfers);

            //find the start index for each processor in lTransferElem array
            for (c=0; c< ibmRank; c++){
                startInd += gNumTransfers[c];
            }

            //loop through the elements and fill the lTransferElem array
            for(e = 0; e < ibMsh->elems; e++)
            {
                if(eBox->thisElemTransfered[e])
                {
                    lTransferElem[startInd] = e;
                    startInd++;
                }
            }

            //scatter the lTransferElem array among all the ibm processors
            MPI_Allreduce(&lTransferElem[0], &(gTransferElem[0]), totalNumTransfers, MPIU_INT, MPI_SUM, ibmBody->IBM_COMM);

            //loop through the transferred elements and see if they are in the processor buffer zone for transfer
            for(c = 0; c < totalNumTransfers; c++)
            {
                e = gTransferElem[c];

                //element center
                p1 = nCoor[nv1[e]];
                p2 = nCoor[nv2[e]];
                p3 = nCoor[nv3[e]];

                if(    isInsideBoundingBox(p1, eBox->outerZone)
                    || isInsideBoundingBox(p2, eBox->outerZone)
                    || isInsideBoundingBox(p3, eBox->outerZone)
                )
                {
                    eBox->thisElemControlled[e] = 1;
                }

            }

            if(ibm->dbg)
            {
                PetscInt sum = 0;
                for(e = 0; e < ibMsh->elems; e++)
                {
                    if(eBox->thisElemControlled[e])
                    {
                        sum++;
                    }

                }
                PetscPrintf(PETSC_COMM_SELF, "rank = %d, ibm elements controlled= %ld, transferred = %ld\n", ibmRank, sum, lNumTransfers[ibmRank]);
            }

            //clean the temporary vectors
            std::vector<PetscInt> ().swap(lNumTransfers);
            std::vector<PetscInt> ().swap(gNumTransfers);

            std::vector<PetscInt> ().swap(lTransferElem);
            std::vector<PetscInt> ().swap(gTransferElem);

        }

    }

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode elementBoundingSphere(ibmObject *ibmBody)
{
    ibmMesh     *ibMesh = ibmBody->ibMsh;
    elementBox  *eBox = ibmBody->eBox;
    PetscInt    n1, n2, n3;
    Cmpnts      p1, p2, p3;
    Cmpnts      pa, pb, pc, pf, pu, pv, pd, pt;
    PetscReal   l12, l23, l31, lu, lv;
    PetscReal   gama, lamda;

    Cmpnts      *qvec = ibMesh->eQVec;
    PetscReal   *rvec = ibMesh->eRVec;

    //check if processor controls this ibm body
    if(ibmBody->ibmControlled)
    {
        for (PetscInt i=0; i<ibMesh->elems; i++)
        {
            //this processor controls this ibm element
            if(eBox->thisElemControlled[i])
            {
                // get the element nodes
                n1 = ibMesh->nID1[i]; n2 = ibMesh->nID2[i]; n3 = ibMesh->nID3[i];

                p1 = ibMesh->nCoor[n1]; p2 = ibMesh->nCoor[n2]; p3 = ibMesh->nCoor[n3];

                l12 = nMag(nSub(p1, p2));
                l23 = nMag(nSub(p2, p3));
                l31 = nMag(nSub(p3, p1));

                //Find the longest edge and assign the corresponding two vertices to pa and pb
                if (l12 > l23)
                {
                    if (l12 > l31)
                    {
                        pa = p1;
                        pb = p2;
                        pc = p3;
                    }
                    else
                    {
                        pa = p3;
                        pb = p1;
                        pc = p2;
                    }
                }
                else
                {
                    if (l31 < l23)
                    {
                        pa = p2;
                        pb = p3;
                        pc = p1;
                    }
                    else
                    {
                        pa = p3;
                        pb = p1;
                        pc = p2;
                    }
                }

                pf = nSum(pa, pb);
                mScale(0.5, pf);

                // u = a - f; v = c - f;
                pu = nSub(pa, pf);
                pv = nSub(pc, pf);

                // d = (u X v) X u;
                pt = nCross(pu, pv);
                pd = nCross(pt, pu);

                // gama = (v^2 - u^2) / (2 d \dot (v - u));
                lu   = nMag(pu);
                lv   = nMag(pv);

                gama = lv * lv - lu * lu ;

                pt    = nSub(pv, pu);
                lamda = 2.0 * nDot(pd, pt);

                gama  /= lamda;

                if (gama < 0)
                {
                    lamda = 0;
                }
                else
                {
                    lamda = gama;
                }


                qvec[i] = nSet(pd);
                mScale(lamda, qvec[i]);
                mSum(qvec[i], pf);

                rvec[i] = nMag(nSub(qvec[i], pa));
            }

        }

    }

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode findBodyBoundingBox(ibm_ *ibm)
{

  boundingBox    ibBox;

  //loop through the IBM bodies
  for (PetscInt i = 0; i < ibm->numBodies; i++)
  {

    //initialize the bounding box
    ibBox.xmin =  1.e23;
    ibBox.xmax = -1.e23;
    ibBox.ymin =  1.e23;
    ibBox.ymax = -1.e23;
    ibBox.zmin =  1.e23;
    ibBox.zmax = -1.e23;

    ibBox.Lx   =  0.;
    ibBox.Ly   =  0.;
    ibBox.Lz   =  0.;

    // get pointer to the current body mesh
    ibmMesh *ibMsh = ibm->ibmBody[i]->ibMsh;

    // loop through the nodes of the IBM Mesh
    for (PetscInt n=0; n < ibMsh->nodes; n++)
    {
      ibBox.xmin = PetscMin(ibBox.xmin, ibMsh->nCoor[n].x);
      ibBox.xmax = PetscMax(ibBox.xmax, ibMsh->nCoor[n].x);

      ibBox.ymin = PetscMin(ibBox.ymin, ibMsh->nCoor[n].y);
      ibBox.ymax = PetscMax(ibBox.ymax, ibMsh->nCoor[n].y);

      ibBox.zmin = PetscMin(ibBox.zmin, ibMsh->nCoor[n].z);
      ibBox.zmax = PetscMax(ibBox.zmax, ibMsh->nCoor[n].z);
    }

    // set a buffer for the bounding box
    ibBox.xmin -= 0.05;
    ibBox.xmax += 0.05;
    ibBox.ymin -= 0.05;
    ibBox.ymax += 0.05;
    ibBox.zmin -= 0.05;
    ibBox.zmax += 0.05;

    ibBox.Lx = ibBox.xmax - ibBox.xmin;
    ibBox.Ly = ibBox.ymax - ibBox.ymin;
    ibBox.Lz = ibBox.zmax - ibBox.zmin;

    //save the bounding box for this body
    ibm->ibmBody[i]->bound->xmin = ibBox.xmin;
    ibm->ibmBody[i]->bound->xmax = ibBox.xmax;
    ibm->ibmBody[i]->bound->ymin = ibBox.ymin;
    ibm->ibmBody[i]->bound->ymax = ibBox.ymax;
    ibm->ibmBody[i]->bound->zmin = ibBox.zmin;
    ibm->ibmBody[i]->bound->zmax = ibBox.zmax;

    ibm->ibmBody[i]->bound->Lx = ibBox.Lx;
    ibm->ibmBody[i]->bound->Ly = ibBox.Ly;
    ibm->ibmBody[i]->bound->Lz = ibBox.Lz;

    if(ibm->dbg)
    {
      PetscPrintf(PETSC_COMM_WORLD, "Bounding box for ibm object: %s\n\n", ibm->ibmBody[i]->bodyName.c_str());
      PetscPrintf(PETSC_COMM_WORLD, "xmin = %lf\n", ibBox.xmin);
      PetscPrintf(PETSC_COMM_WORLD, "xmax = %lf\n", ibBox.xmax);
      PetscPrintf(PETSC_COMM_WORLD, "ymin = %lf\n", ibBox.ymin);
      PetscPrintf(PETSC_COMM_WORLD, "ymax = %lf\n", ibBox.ymax);
      PetscPrintf(PETSC_COMM_WORLD, "zmin = %lf\n", ibBox.zmin);
      PetscPrintf(PETSC_COMM_WORLD, "zmax = %lf\n", ibBox.zmax);

      PetscPrintf(PETSC_COMM_WORLD, "Lx = %lf\n", ibBox.Lx);
      PetscPrintf(PETSC_COMM_WORLD, "Ly = %lf\n", ibBox.Ly);
      PetscPrintf(PETSC_COMM_WORLD, "Lz = %lf\n", ibBox.Lz);
    }
  }
  return 0;
}

//***************************************************************************************************************//

PetscErrorCode findFluidSupportNodes(ibm_ *ibm)
{
  mesh_         *mesh = ibm->access->mesh;
  DM            da = mesh->da, fda = mesh->fda;
  DMDALocalInfo info = mesh->info;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;
  PetscInt      gxs = info.gxs, gxe = info.gxs + info.gxm;
  PetscInt      gys = info.gys, gye = info.gys + info.gym;
  PetscInt      gzs = info.gzs, gze = info.gzs + info.gzm;

  PetscInt      i, j, k, ii, kk, jj, c;
  PetscReal     beta = 2.1, normDist, ***nvert;
  Cmpnts        ***cent, dist;

  DMDAVecGetArray(fda, mesh->lCent, &cent);
  DMDAVecGetArray(da, mesh->lNvert, &nvert);

  //local pointer for ibmFluidCells
  ibmFluidCell *ibF = ibm->ibmFCells;

  for(c = 0; c < ibm->numIBMFluid; c++)
  {
    //cell indices
    i = ibF[c].cellId.i;
    j = ibF[c].cellId.j;
    k = ibF[c].cellId.k;

    //initialize the fluid support List of the ibm fluid cell i, j, k
    initCellList(&(ibF[c].flNodes));

    // set the support radius for the least square method
    if (k==mz-2 || j==my-2 || i==mx-2)
    {
      dist = nSub(cent[k][j][i], cent[k-1][j-1][i-1]);
      ibF[c].rad = nMag(dist);
    }
    else if (k==1 || j==1 || i==1)
    {
      dist = nSub(cent[k+1][j+1][i+1], cent[k][j][i]);
      ibF[c].rad = nMag(dist);
    }
    else
    {
      dist = nSub(cent[k+1][j+1][i+1], cent[k-1][j-1][i-1]);
      ibF[c].rad = 0.5 * nMag(dist);
    }

    //scaled support radius
    ibF[c].rad = beta * ibF[c].rad;

    //initialize the number of fluid support cells around ibm fluid cell i,j,k
    ibF[c].numFl = 0;

    for (kk = k-floor(beta); kk <= k+floor(beta); kk++)
    {
      //do not include cells outside the ghost cells of the processor - will cause segmentation fault
      if (kk < gzs || kk >= gze || kk >= mz-2 || kk <= 1 ) continue;

      for (jj = j-floor(beta); jj <= j+floor(beta); jj++)
      {
        if (jj < gys || jj >= gye || jj >= my-2 || jj <= 1) continue;

        for (ii = i-floor(beta); ii <= i+floor(beta); ii++)
        {
          if (ii < gxs || ii >= gxe || ii >= mx-2 || ii <= 1) continue;

          if(isFluidCell(kk, jj, ii, nvert))
          {
            //normalized distance from the ibm fluid node
            dist = nSub(cent[k][j][i],cent[kk][jj][ii]);
            normDist = nMag(dist);
            normDist = normDist / ibF[c].rad;

            //check that cell is within the support and not the ibm fluid cell
            if (normDist <= 1 && !(ii==i && jj==j && kk==k))
            {
              cellIds tmp;
              tmp.i = ii; tmp.j = jj; tmp.k = kk;

              //insert into the fluid support list
              insertCellNode1(&(ibF[c].flNodes), tmp);
              ibF[c].numFl++;
            }

          }

        }
      }
    }

  }

  DMDAVecRestoreArray(fda, mesh->lCent, &cent);
  DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

  return 0;
}

//***************************************************************************************************************//

PetscErrorCode findIBMFluidCells(ibm_ *ibm)
{
  mesh_         *mesh = ibm->access->mesh;
  DM            da = mesh->da, fda = mesh->fda;
  DMDALocalInfo info = mesh->info;
  PetscInt      xs = info.xs, xe = info.xs + info.xm;
  PetscInt      ys = info.ys, ye = info.ys + info.ym;
  PetscInt      zs = info.zs, ze = info.zs + info.zm;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  PetscInt      i, j, k, ctr;
  PetscReal     ***nvert;

  lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
  lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
  lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

  DMDAVecGetArray(da, mesh->lNvert, &nvert);

  ctr = 0;

  for (k = lzs; k < lze; k++)
  for (j = lys; j < lye; j++)
  for (i = lxs; i < lxe; i++)
  {
    if(isIBMFluidCell(k, j, i, nvert))
    {
      ctr++;
    }
  }
  ibm->numIBMFluid = ctr;

  PetscMalloc(ctr * sizeof(ibmFluidCell), &(ibm->ibmFCells));
  ibmFluidCell *ibmFluid = ibm->ibmFCells;

  ctr = 0;

  for (k = lzs; k < lze; k++)
  for (j = lys; j < lye; j++)
  for (i = lxs; i < lxe; i++)
  {
        if(isIBMFluidCell(k, j, i, nvert))
        {
          ibmFluid[ctr].cellId.i = i;
          ibmFluid[ctr].cellId.j = j;
          ibmFluid[ctr].cellId.k = k;

          ctr++;
        }
  }

  DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

  return 0;
}

//***************************************************************************************************************//

PetscErrorCode findIBMMeshSupportNodes(ibm_ *ibm)
{
  mesh_         *mesh = ibm->access->mesh;
  DM            da = mesh->da, fda = mesh->fda;
  DMDALocalInfo info = mesh->info;
  PetscInt      mx = info.mx, my = info.my, mz = info.mz;

  PetscInt      i, j, k, ii, kk, jj;
  PetscInt      b, c, s;
  PetscReal     normDist, ***nvert;
  Cmpnts        ***cent, dist;
  cellList      localSearchList;

  DMDAVecGetArray(fda, mesh->lCent, &cent);
  DMDAVecGetArray(da, mesh->lNvert, &nvert);

  //local pointer for ibmFluidCells
  ibmFluidCell *ibF = ibm->ibmFCells;

  for(c = 0; c < ibm->numIBMFluid; c++)
  {
      //cell indices
      i = ibF[c].cellId.i;
      j = ibF[c].cellId.j;
      k = ibF[c].cellId.k;

    // initialize the ibm node support list
    initlist(&(ibF[c].slNodes));

    // initialize the list of ibm body each node in the support list belongs to
    initlist(&(ibF[c].bID));

    //initialize number of support cells to 0
    ibF[c].numSl = 0;

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {

        boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
        searchBox     *sBox  = &(ibm->sBox[b]);

        //initialize the local search list for search cells within the support radius
        initCellList(&localSearchList);

        // loop through the searchCellList
        for(ii = 0; ii < sBox->ncx; ii++)
        for(jj = 0; jj < sBox->ncy; jj++)
        for(kk = 0; kk < sBox->ncz; kk++)
        {
          cellIds tmp;
          tmp.i = ii;
          tmp.j = jj;
          tmp.k = kk;

          if(isSearchCellSupport(ibF[c].rad, cent[k][j][i], tmp, ibBox, sBox))
          {
            insertCellNode1(&localSearchList, tmp);
          }
        }

        //get access to the first search cell in the list
        cellNode *currentCell = localSearchList.head;

        // loop through the search cells in the support radius
        while(currentCell)
        {
          // get the searchcell index
          cellIds cellInd;
          cellInd.i = currentCell->Node.i;
          cellInd.j = currentCell->Node.j;
          cellInd.k = currentCell->Node.k;

          insertIBMSupportNodes(cent[k][j][i], cellInd, &ibF[c], ibm->ibmBody[b], sBox);

          currentCell = currentCell->next;
        }

        destroyCellList(&localSearchList);


    }

  }

  DMDAVecRestoreArray(fda, mesh->lCent, &cent);
  DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

  return 0;
}

//***************************************************************************************************************//

PetscErrorCode findSearchCellDim(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***aj, lcellSize = 0., cellSize = 0.;
    PetscInt      lsum = 0., sum = 0.;
    PetscReal     searchCellsize;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->lAj, &aj);

    //find average cellsize
    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
      lcellSize += pow( 1.0/aj[k][j][i], 1./3.);
      lsum ++ ;
    }

    MPI_Allreduce(&lcellSize, &cellSize, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lsum, &sum, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    cellSize /= sum;

    for (PetscInt b = 0; b < ibm->numBodies; b++)
    {
        //compare the searchCellsize to the body bounding box
        boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
        searchBox     *sBox  = &(ibm->sBox[b]);
        //search cell size taken to be n times the average cell size
        searchCellsize = ibm->ibmBody[b]->searchCellRatio * cellSize;

        sBox->ncx = ceil(ibBox->Lx/searchCellsize);
        sBox->ncy = ceil(ibBox->Ly/searchCellsize);
        sBox->ncz = ceil(ibBox->Lz/searchCellsize);

        //set a default value (this needs optimization)
        // if(b == 1)
        // {
            sBox->ncx = 10; sBox->ncy = 10; sBox->ncz = 10;
        // }

        // cell size for the search algorithm
        sBox->dcx = (ibBox->Lx) / (sBox->ncx);
        sBox->dcy = (ibBox->Ly) / (sBox->ncy);
        sBox->dcz = (ibBox->Lz) / (sBox->ncz);

        if(ibm->dbg)
        {
          PetscPrintf(mesh->MESH_COMM, "\n ncx = %ld, ncy = %ld, ncz = %ld\n", sBox->ncx, sBox->ncy ,sBox->ncz);
          PetscPrintf(mesh->MESH_COMM, "\n dcx = %lf, dcy = %lf, dcz = %lf\n", sBox->dcx, sBox->dcy ,sBox->dcz);
        }
    }

    DMDAVecRestoreArray(da, mesh->lAj, &aj);

  return(0);
}

//***************************************************************************************************************//

// creates the search cell list
// each element of the search list has an array of the IBM elements within than search cell
// when the ray cast algorithm is used, need to go only through the search cells the ray is casted along
PetscErrorCode createSearchCellList(ibm_ *ibm)
{
    PetscReal     xv_min, yv_min, zv_min, xv_max, yv_max, zv_max; // min and max co-ordinate of the 3 vertices of an element
    PetscInt      iv_min, iv_max, jv_min, jv_max, kv_min, kv_max; // min and max search cell index of an IBM element
    PetscInt      n1, n2, n3;                                     // nodes of a particular IBM element
    PetscInt      i, j, k, e, b;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {

        boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
        ibmMesh       *ibMsh = ibm->ibmBody[b]->ibMsh;                         // pointer to the ibm body mesh

        searchBox *sBox = &(ibm->sBox[b]);

        // allocate memory for the search cell list
        PetscMalloc(sBox->ncz * sBox->ncy * sBox->ncx * sizeof(list), &(ibm->ibmBody[b]->searchCellList));
        list  *searchCellList = ibm->ibmBody[b]->searchCellList;

        //initialize each array element of the search cell list to null pointer
        for (k = 0; k < sBox->ncz; k++) {
          for (j = 0; j < sBox->ncy; j++) {
              for (i = 0; i < sBox->ncx; i++) {
                  initlist(&searchCellList[k * sBox->ncx * sBox->ncy + j * sBox->ncx + i]);
              }
          }
        }

        //insert the ibm body triangular mesh elements into the search cell list based on its position
        //loop through the ibm mesh elements

        for(e = 0; e < ibMsh->elems; e++)
        {
            // 3 vertices of the element
            n1 = ibMsh->nID1[e];
            n2 = ibMsh->nID2[e];
            n3 = ibMsh->nID3[e];

            PetscReal solutionTimeStart, solutionTimeEnd;

            if(e%1000 == 1 && ibm->dbg)
            {
                PetscTime(&solutionTimeStart);
                PetscPrintf(PETSC_COMM_WORLD, "element = %ld\n", e);
            }

            // min and max coordinate value of the element vertices
            xv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].x, ibMsh->nCoor[n2].x ), ibMsh->nCoor[n3].x );
            xv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].x, ibMsh->nCoor[n2].x ), ibMsh->nCoor[n3].x );

            yv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].y, ibMsh->nCoor[n2].y ), ibMsh->nCoor[n3].y );
            yv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].y, ibMsh->nCoor[n2].y ), ibMsh->nCoor[n3].y );

            zv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].z, ibMsh->nCoor[n2].z ), ibMsh->nCoor[n3].z );
            zv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].z, ibMsh->nCoor[n2].z ), ibMsh->nCoor[n3].z );

            // min and max index of the search cell where this element is located
            iv_min = floor((xv_min - ibBox->xmin) / sBox->dcx);          // find the search cell index of xmin coordinate of the element
            iv_max = floor((xv_max - ibBox->xmin) / sBox->dcx);

            jv_min = floor((yv_min - ibBox->ymin) / sBox->dcy); //
            jv_max = floor((yv_max - ibBox->ymin) / sBox->dcy);

            kv_min = floor((zv_min - ibBox->zmin) / sBox->dcz); //
            kv_max = floor((zv_max - ibBox->zmin) / sBox->dcz);

            //ensure that the search cell indices of the current element are bounded between 0 and max index
            iv_min = (iv_min < 0) ? 0 : iv_min;
            iv_max = (iv_max > sBox->ncx-1) ? sBox->ncx-1 : iv_max;

            jv_min = (jv_min < 0) ? 0 : jv_min;
            jv_max = (jv_max > sBox->ncy-1) ? sBox->ncy-1 : jv_max;

            kv_min = (kv_min < 0) ? 0 : kv_min;
            kv_max = (kv_max > sBox->ncz-1) ? sBox->ncz-1 : kv_max;

            //insert element into search cell
            for (k = kv_min; k <= kv_max; k++) {
                for (j = jv_min; j <= jv_max; j++) {
                    for (i = iv_min; i <= iv_max; i++) {
                        insertnode(&(searchCellList[k * sBox->ncx * sBox->ncy + j * sBox->ncx + i]), e);
                    }
                }
            }

            if(e%1000 == 0 && ibm->dbg)
            {
                PetscTime(&solutionTimeEnd);
                PetscPrintf(PETSC_COMM_WORLD, "Total ibm search cell time = %lf s\n", solutionTimeEnd - solutionTimeStart);
            }
        }
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode destroyLists(ibm_ *ibm)
{
    PetscInt      i, j, k, b, c;

    // loop through the ibm bodies and destroy the searchCellList
    for(b = 0; b < ibm->numBodies; b++)
    {

        PetscInt ncx = ibm->sBox[b].ncx;
        PetscInt ncy = ibm->sBox[b].ncy;
        PetscInt ncz = ibm->sBox[b].ncz;

        for (k = 0; k < ncz; k++) {
            for (j = 0; j < ncy; j++) {
                for (i = 0; i < ncx; i++) {
                    destroy(&(ibm->ibmBody[b]->searchCellList[k * ncx * ncy + j * ncx + i]));
                }
            }
        }

        PetscFree(ibm->ibmBody[b]->searchCellList);

    }


    if (ibm->IBInterpolationModel == "MLS")
    {
        //local pointer for ibmFluidCells
        ibmFluidCell *ibF = ibm->ibmFCells;

        for(c = 0; c < ibm->numIBMFluid; c++)
        {

            destroyCellList(&(ibF[c].flNodes));

            destroy(&(ibF[c].slNodes));

            destroy(&(ibF[c].bID));
        }

    }

    PetscFree(ibm->ibmFCells);

    return (0);
}

//***************************************************************************************************************//

inline void insertIBMSupportNodes(Cmpnts pt, cellIds sCell, ibmFluidCell *ibF, ibmObject *ibBody, searchBox *sBox)
{
  PetscInt       elemID, i;
  PetscInt       n[3];
  Cmpnts    dist;
  PetscReal    normDist;
  ibmMesh   *ibMsh = ibBody->ibMsh;

  // get access to the first IBM element in a search cell
  node *currentElem = ibBody->searchCellList[sCell.k * sBox->ncx * sBox->ncy + sCell.j * sBox->ncx + sCell.i].head;

  //loop through the IBM elements within the current search cell
  while(currentElem)
  {
    elemID = currentElem->Node;
    n[0]   = ibMsh->nID1[elemID];
    n[1]   = ibMsh->nID2[elemID];
    n[2]   = ibMsh->nID3[elemID];

    for(i = 0; i < 3; i++)
    {
      dist = nSub(pt, ibMsh->nCoor[n[i]]);
      normDist = nMag(dist)/ibF->rad;

      if (normDist<=1)
      {
        //insert node if not previously added
        if(insertnode(&(ibF->slNodes), n[i]))
        {
          //list of the body id of the node support ibm node inserted
          insertnode1(&(ibF->bID), ibBody->bodyID);
          ibF->numSl++;
        }
      }
    }

    currentElem = currentElem->next;
  }
  return;
}

//***************************************************************************************************************//
PetscErrorCode rayCastLocal(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, ibmMesh *ibMsh, cellIds sCell, searchBox *sBox, boundingBox *ibBox, list *searchCellList, PetscInt &intersect)
{
    PetscInt	i, j, k, e;
    PetscBool   cutElement;
    PetscInt    ks, js, is;
    PetscInt    kp, jp, ip;

    node        *current;
    PetscBool	*Element_Searched;

    kp = sCell.k;
    jp = sCell.j;
    ip = sCell.i;

    ks=PetscMax(kp-1,0);
    js=PetscMax(jp-1,0);
    is=PetscMax(ip-1,0);

    Element_Searched = new PetscBool[ibMsh->elems];

    for (e = 0; e < ibMsh->elems; e++)
    {
        Element_Searched[e] = PETSC_FALSE;
    }

    for (k=ks; k<kp+2 && k<sBox->ncz; k++)
    {
        for (j=js; j<jp+2 && j<sBox->ncy; j++)
        {
            for (i=is; i<ip+2 && i<sBox->ncx; i++)
            {
                current = searchCellList[k * sBox->ncx * sBox->ncy + j * sBox->ncx + i].head;

                while(current)
                {
                    // element currently pointed to
                    e = current->Node;

                    if (!Element_Searched[e])
                    {
                        Element_Searched[e] = PETSC_TRUE;

                        cutElement = isLineTriangleInt(p, p1, ibMsh, e);

                        if (cutElement)
                        {
                            break;
                        }

                        cutElement = isLineTriangleInt(p, p2, ibMsh, e);

                        if (cutElement)
                        {
                            break;
                        }

                        cutElement = isLineTriangleInt(p, p3, ibMsh, e);

                        if (cutElement)
                        {
                            break;
                        }

                        cutElement = isLineTriangleInt(p, p4, ibMsh, e);

                        if (cutElement)
                        {
                            break;
                        }
                    }

                    current = current->next;
                }
                if (cutElement)
                {
                    break;
                }
            }
            if (cutElement)
            {
                break;
            }
        }
        if (cutElement)
        {
            break;
        }
    }

    if (cutElement)
    {
        intersect = 2;
    }
    else
    {
        intersect = 0;
    }

    PetscFree(Element_Searched);
    return (0);
}

//***************************************************************************************************************//

PetscReal rayCastingTest(Cmpnts p, ibmMesh *ibMsh, cellIds sCell, searchBox *sBox, boundingBox *ibBox, list *searchCellList)
{
  PetscBool     *Element_Searched;                    //bool to indicate if an element has been searched or not
  PetscBool      NotDecided = PETSC_TRUE;
  PetscBool      Singularity = PETSC_FALSE;

  PetscInt       searchtimes = 0, nvertLoc, numInt;                     //number of times search is performed - each time a new ray is cast
  Cmpnts         dir, nor, dnn[1000];
  node           *current;                           //pointer to the first element in the current searchCellList node
  PetscInt            e, i, j, k;
  PetscInt            n1, n2, n3;                         // vertices of an ibm mesh element
  PetscReal      epsilon = 1.e-25;
  PetscReal      dt[1000], ndotn, dirdotn;
  PetscReal      t, u, v;

  PetscMalloc(ibMsh->elems * sizeof(PetscBool), &Element_Searched);

  j = sCell.j;
  i = sCell.i;

  while (NotDecided)
  {
    searchtimes++;
    numInt = 0;
    Singularity = PETSC_FALSE;

    for (e = 0; e < ibMsh->elems; e++)
    {
        Element_Searched[e] = PETSC_FALSE;
    }

    //get a random direction from the fluid mesh point in the search box z direction
    dir = randomdirection(p, sCell, ibBox, sBox, searchtimes);

    // ray cast is done along the z direction of the search cells
    // as the ray is cast in the positive z cell direction, need to search only between point p.z and ncz
    for (k = sCell.k; k < sBox->ncz; k++)
    {

      current = searchCellList[k * sBox->ncx * sBox->ncy + j * sBox->ncx + i].head;

      // loop through the elements in the current search cell
      while (current)
      {
        // element currently pointed to
        e = current->Node;

        if (!Element_Searched[e])
        {
          Element_Searched[e] = PETSC_TRUE;
          n1 = ibMsh->nID1[e];
          n2 = ibMsh->nID2[e];
          n3 = ibMsh->nID3[e];
          nor = ibMsh->eN[e];

          dirdotn = nDot(dir, nor);

          // t: is the magnitude of the ray (positive ore negative based on the intersection side)
          // u and v are a and b of the Borazjani, Ge, Sotiropulos 2008 paper
          nvertLoc = intsectElement(p, dir, ibMsh->nCoor[n1], ibMsh->nCoor[n2], ibMsh->nCoor[n3], &t, &u, &v);

          // if ray intersects element
          if (nvertLoc > 0 && t > 0)
          {
            //store the element normal and dist from point to element at current intersection
            dt[numInt] = t;
            dnn[numInt] = nSet(nor);

            numInt++;

            // if 2 intersections points are at same distance from point p
            // ensure no singularity. If singularity point, recast and find number of intersections again
            for( PetscInt tmp = 0; tmp < numInt - 1; tmp++)
            {
              ndotn = nDot(dnn[tmp], nor);

              if(fabs(t - dt[tmp]) < epsilon && ndotn > -0.99)
              {
                Singularity = PETSC_TRUE;
                break;
              }
            }

          }

        }
        if(Singularity)
        {
          break;
        }
        else
        {
          current = current->next;
        }
      }
      if (Singularity)
      {
          break;
      }
    }

    if (!Singularity)
    {
        // The interception point number is odd, inside body else outise the IBM
        if (numInt % 2)
        {

            PetscFree(Element_Searched);
            return 4.;
        }
        else
        {
            PetscFree(Element_Searched);
            return 0.;
        }
    }
    // if code control reaches here, then singularity has occured so recast ray to repeat until no singularity
  }
  PetscFree(Element_Searched);

  return 0.;
}

//***************************************************************************************************************//

//returns true if search cell intersects with the support radius of the point pt
inline bool isSearchCellSupport(PetscReal rad, Cmpnts pt, cellIds cID, boundingBox *ibBox, searchBox *sBox)
{

    Cmpnts Bmin, Bmax;
    PetscReal dist = rad*rad;

    // co-ordinates of the min and max of the cell box
    Bmin.x = ibBox->xmin + (cID.i)*sBox->dcx;
    Bmin.y = ibBox->ymin + (cID.j)*sBox->dcy;
    Bmin.z = ibBox->zmin + (cID.k)*sBox->dcz;

    Bmax.x = ibBox->xmin + (cID.i+1)*sBox->dcx;
    Bmax.y = ibBox->ymin + (cID.j+1)*sBox->dcy;
    Bmax.z = ibBox->zmin + (cID.k+1)*sBox->dcz;

    if (pt.x < Bmin.x)      dist = dist - (pt.x - Bmin.x)*(pt.x - Bmin.x);
    else if (pt.x > Bmax.x) dist = dist - (pt.x - Bmax.x)*(pt.x - Bmax.x);

    if (pt.y < Bmin.y)      dist = dist - (pt.y - Bmin.y)*(pt.y - Bmin.y);
    else if (pt.y > Bmax.y) dist = dist - (pt.y - Bmax.y)*(pt.y - Bmax.y);

    if (pt.z < Bmin.z)      dist = dist - (pt.z - Bmin.z)*(pt.z - Bmin.z);
    else if (pt.z > Bmax.z) dist = dist - (pt.z - Bmax.z)*(pt.z - Bmax.z);

    return dist > 0;
}

//***************************************************************************************************************//

inline Cmpnts randomdirection(Cmpnts p, cellIds sCell, boundingBox *ibBox, searchBox *sBox, PetscInt seed)
{
  Cmpnts endpoint, dir;
  PetscReal s;
  PetscReal xpc, ypc;

//Find the co-ordinate of the center of the search cell (xc, yc) containing the current node p
  xpc = sBox->dcx * (sCell.i + 0.5) + ibBox->xmin;
  ypc = sBox->dcy * (sCell.j + 0.5) + ibBox->ymin;

  // seeded differently each search time to get a new direction
  srand(seed);

// A ray is cast in the positive z cell direction from p such that its endpoint lies
// between (xc + 0.5 and xc - 0.5, yc + 0.5, yc - 0.5)

  s = rand() / ((PetscReal) RAND_MAX + 1) - 0.5;

  endpoint.x = xpc + s * sBox->dcx;
  endpoint.y = ypc + s * sBox->dcy;
  endpoint.z = ibBox->zmax + 0.2;

  dir = nSub(endpoint, p);

  // unit vector along the cast ray from p
  s = nMag(dir);
  mScale(1.0/s, dir);

  return dir;
}

//***************************************************************************************************************//

//to check if a ray from point p intersects an IBM elemtent. return 1 - success
inline PetscInt intsectElement(Cmpnts p, Cmpnts dir, Cmpnts node1, Cmpnts node2, Cmpnts node3, PetscReal *t, PetscReal *u, PetscReal *v)
{
  Cmpnts    edge1, edge2, tvec, pvec, qvec;
  PetscReal    det, invDet;
  PetscReal epsilon = 1.e-15;
  edge1 = nSub(node2, node1);
  edge2 = nSub(node3, node1);

  pvec = nCross(dir, edge2);
  det = nDot(edge1, pvec);

  if(det > -epsilon && det < epsilon)
  {
    return 0;
  }

  invDet = 1.0 / det;

  tvec = nSub(p, node1);

  *u = nDot(tvec, pvec) * invDet;

  if (*u < 0.0 || *u > 1.0)
  {
      return 0;
  }

  qvec = nCross(tvec, edge1);

  *v = nDot(dir, qvec) *invDet;

  if (*v < 0.0 || *u + *v > 1.0)
  {
      return 0;
  }

  *t = nDot(edge2, qvec) * invDet;
  return 1;
}

//***************************************************************************************************************//

inline PetscBool isLineTriangleInt(Cmpnts p1, Cmpnts p2, ibmMesh *ibMsh, PetscInt e)
{
    PetscInt  cutthrough = 0;
    PetscInt  n1e, n2e, n3e;

    Cmpnts    pe1, pe2, pe3, nor;
    Cmpnts    dp1, p2p1, pe2pe1, pe3pe1, pint;

    PetscReal tf1, tf2, d;

    n1e = ibMsh->nID1[e];
    n2e = ibMsh->nID2[e];
    n3e = ibMsh->nID3[e];

    pe1 = ibMsh->nCoor[n1e];
    pe2 = ibMsh->nCoor[n2e];
    pe3 = ibMsh->nCoor[n3e];

    nor = ibMsh->eN[e];

    tf1 = nDot(nSub(p1, pe1), nor);
    tf2 = nDot(nSub(p2, pe1), nor);

    if(tf1 * tf2 < 0)
    {
        //p1 & p2 on different sides of triangle
        dp1 = nSub(p1, pe1);
        p2p1 = nSub(p2, p1);

        pe2pe1 = nSub(pe2, pe1);
        pe3pe1 = nSub(pe3, pe1);

        PetscReal dx1, dy1, dz1, dx2, dy2, dz2, dx3, dy3, dz3;

        PetscReal nftx, nfty, nftz;

        dx1 = dp1.x, dy1 = dp1.y, dz1 = dp1.z;

        nftx = p2p1.x; nfty = p2p1.y; nftz = p2p1.z;

        dx2 = pe2pe1.x; dy2 = pe2pe1.y; dz2 = pe2pe1.z;
        dx3 = pe3pe1.x; dy3 = pe3pe1.y; dz3 = pe3pe1.z;

        d = - (dx1 * (dy2 * dz3 - dz2 * dy3) -
           dy1 * (dx2 * dz3 - dz2 * dx3) +
           dz1 * (dx2 * dy3 - dy2 * dx3)) /
             (nftx * (dy2 * dz3 - dz2 * dy3) -
          nfty * (dx2 * dz3 - dz2 * dx3) +
          nftz * (dx2 * dy3 - dy2 * dx3));

        pint = nSum(p1, nScale(d, p2p1));

        cutthrough = isPointInTriangle(pint, pe1, pe2, pe3, nor);

        if (cutthrough == 1)
        {
            return (PETSC_TRUE);
        }
    }

    return (PETSC_FALSE);
}

//***************************************************************************************************************//

inline PetscInt isPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts norm)
{
  PetscInt   flag;
  Cpt2D	     pj, pj1, pj2, pj3;
  PetscReal	 epsilon = 1.e-10, area;

  if ( fabs(norm.z) >= fabs(norm.x) && fabs(norm.z) >= fabs(norm.y) )
  {
    pj.x  = p.x;  pj.y  = p.y;
    pj1.x = p1.x; pj1.y = p1.y;
    pj2.x = p2.x; pj2.y = p2.y;
    pj3.x = p3.x; pj3.y = p3.y;
  }
  else if (fabs(norm.x) >= fabs(norm.y) && fabs(norm.x) >= fabs(norm.z))
  {
    pj.x  = p.z;  pj.y  = p.y;
    pj1.x = p1.z; pj1.y = p1.y;
    pj2.x = p2.z; pj2.y = p2.y;
    pj3.x = p3.z; pj3.y = p3.y;
  }
  else {
    pj.x  = p.x;  pj.y  = p.z;
    pj1.x = p1.x; pj1.y = p1.z;
    pj2.x = p2.x; pj2.y = p2.z;
    pj3.x = p3.x; pj3.y = p3.z;
  }

  area = pj1.x * (pj2.y - pj3.y) + pj2.x * (pj3.y - pj1.y) + pj3.x * (pj1.y - pj2.y);

  if(fabs(area) < epsilon)
  {
      flag = -1;
  }
  else
  {
      flag = ISInsideTriangle2D(pj, pj1, pj2, pj3);
  }

  return(flag);
}

//***************************************************************************************************************//

inline PetscInt ISInsideTriangle2D(Cpt2D p, Cpt2D pa, Cpt2D pb, Cpt2D pc)
{
  // Check if point p and p3 is located on the same side of line p1p2
  PetscInt 	ls;

  ls = ISSameSide2D(p, pa, pb, pc);

  if (ls < 0) {
    return (ls);
  }

  ls = ISSameSide2D(p, pb, pc, pa);

  if (ls < 0) {
    return (ls);
  }

  ls = ISSameSide2D(p, pc, pa, pb);

  if (ls <0) {
    return(ls);
  }

  return (ls);
}

//***************************************************************************************************************//

/*  Check whether 2D point p is located on the same side of line p1p2
    with point p3. Returns:
    -1	different side
     1	same side (including the case when p is located right on the line)

    If p and p3 is located on the same side to line p1p2, then
    the (p-p1) X (p2-p1) and (p3-p1) X (p2-p1) should have the same sign
*/

inline PetscInt ISSameSide2D(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3)
{
  PetscReal t1, t2, t3;
  PetscReal	epsilon = 1.e-10;
  PetscReal lt;

  t1 = (p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x);
  t2 = (p3.x - p1.x) * (p2.y - p1.y) - (p3.y - p1.y) * (p2.x - p1.x);

  lt = sqrt((p1.x - p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));

  if (fabs(t1/lt) < epsilon)
  { // Point is located along the line of p1p2
    return(1);
  }

  if (t1 > 0)
  {
    if (t2 > 0) return (1); // same side
    else return(-1);  // not
  }
  else
  {
    if (t2 < 0) return(1); // same side
    else return(-1);
  }
}

//***************************************************************************************************************//

inline void disP2Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d, PetscReal *t)
{
  PetscReal	dmin;
  PetscReal	dx21, dy21, dz21, dx31, dy31, dz31;

  dx21 = p2.x - p1.x; dy21 = p2.y - p1.y; dz21 = p2.z - p1.z;
  dx31 = p.x  - p1.x; dy31 = p.y  - p1.y; dz31 = p.z  - p1.z;

  *t = (dx31 * dx21 + dy31 * dy21 + dz31 * dz21) / (dx21*dx21 + dy21*dy21 +
       dz21 * dz21);

  if ((*t)<0)
  { // The closest point is p1
    po->x = p1.x; po->y = p1.y; po->z = p1.z;
    *d = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
  }
  else if ((*t)>1)
  { // The closest point is p2
    po->x = p2.x; po->y = p2.y; po->z = p2.z;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
         (p.z - po->z) * (p.z - po->z));
  }
  else
  { // The closest point lies between p1 & p2
    po->x = p1.x + (*t) * dx21; po->y = p1.y + (*t) * dy21; po->z = p1.z + (*t) *dz21;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
         (p.z - po->z) * (p.z - po->z));
  }

  return;
}

//***************************************************************************************************************//
inline Cmpnts computeVertexAverageNormal(ibmMesh *ibMsh, PetscInt vertexId)
{
    Cmpnts    avgNormal = nSetZero(), eNorm, n;
    HalfEdge  *start_he = ibMsh->vertices[vertexId].edge;
    HalfEdge  *he_1 = start_he;    //first half edge
    HalfEdge  *he_2 = NULL;        //second half edge

    Cmpnts    p1, p2, p3, pc;
    Cmpnts    v1, v2;
    PetscReal alpha, normMag;
    PetscInt  faceId;

    do
    {
        he_2 = he_1->twin->next;

        faceId = he_1->face->faceId;
        eNorm  = ibMsh->eN[faceId];

        p1 = nSet(he_1->origin->nCoor);
        p2 = nSet(he_1->next->origin->nCoor);
        p3 = nSet(he_2->next->origin->nCoor);
        v1 = nSub(p2, p1);
        v2 = nSub(p3, p1);

        pc = nCross(v1,v2);
        n  = nUnit(pc);

        alpha = std::atan2(nDot(n, pc), nDot(v1,v2));

        mScale(alpha, eNorm);
        mSum(avgNormal, eNorm);

        he_1 = he_2;
    }
    while (he_1 != start_he);

    mUnit(avgNormal);

    return avgNormal;
}
//***************************************************************************************************************//
inline Cmpnts computeEdgeAverageNormal(ibmMesh *ibMsh, PetscInt vertexId, PetscInt faceId)
{
    Cmpnts    avgNormal = nSetZero(), eN1, eN2;
    HalfEdge  *start_he = ibMsh->vertices[vertexId].edge;
    HalfEdge  *he_1 = start_he;    //first half edge
    HalfEdge  *he_2 = NULL;        //second half edge

    do
    {
        if(he_1->face->faceId == faceId)
        {
            break;
        }
        else
        {
            he_1 = he_1->twin->next;
        }
    }
    while (he_1 != start_he);

    eN1 = ibMsh->eN[he_1->face->faceId];

    he_2 = he_1->twin;
    eN2 = ibMsh->eN[he_2->face->faceId];

    avgNormal = nSum(eN1, eN2);
    mUnit(avgNormal);

    return avgNormal;
}
//***************************************************************************************************************//

inline void triangleIntp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, ibmFluidCell *ibF)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;

  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a   = x13 * y23 - x23 * y13;

  ibF->cs1 = (y23 * xp3 - x23 * yp3) / a;
  ibF->cs2 = (-y13 * xp3 + x13 * yp3) / a;
  ibF->cs3 = 1. - ibF->cs1 - ibF->cs2;

  return;
}

//***************************************************************************************************************//

inline void triangleIntpBg(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, ibmFluidCell *ibF)
{
  PetscReal  x13, y13, x23, y23, xp3, yp3, a;

  x13 = p1.x - p3.x; y13 = p1.y - p3.y;
  x23 = p2.x - p3.x; y23 = p2.y - p3.y;
  xp3 = p.x - p3.x; yp3 = p.y - p3.y;
  a   = x13 * y23 - x23 * y13;

  ibF->cr1 = (y23 * xp3 - x23 * yp3) / a;
  ibF->cr2 = (-y13 * xp3 + x13 * yp3) / a;
  ibF->cr3 = 1. - ibF->cr1 - ibF->cr2;

  return;
}
