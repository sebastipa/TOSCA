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

        readIBMProperties(ibm);

        setIBMWallModels(ibm);

        findBodyBoundingBox(ibm);

        findSearchCellDim(ibm);

        createSearchCellList(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        computeIBMElementNormal(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        // ibm search algorithm
        PetscPrintf(PETSC_COMM_WORLD, "IBM search algorithm...");
        ibmSearch(ibm);
        PetscPrintf(PETSC_COMM_WORLD, "done\n");

        MPI_Barrier(mesh->MESH_COMM);

        checkIBMexists(ibm);

        findIBMFluidCells(ibm);

        findClosestIBMElement(ibm);

        if(ibm->computeForce)
        {
            findIBMControlledProcs(ibm);

            initFindIBMElementProcs(ibm);
        }

        if (ibm->IBInterpolationModel == "MLS")
        {
            findFluidSupportNodes(ibm);

            findIBMMeshSupportNodes(ibm);

            MPI_Barrier(mesh->MESH_COMM);

        }

    }

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode UpdateIBM(ibm_ *ibm)
{
    mesh_ *mesh = ibm->access->mesh;

    PetscReal ibmTimeStart, ibmTimeEnd;
    PetscTime(&ibmTimeStart);

    //copy current nvert before next time step
    VecCopy(mesh->Nvert, mesh->Nvert_o);
    DMGlobalToLocalBegin(mesh->da, mesh->Nvert_o, INSERT_VALUES, mesh->lNvert_o);
    DMGlobalToLocalEnd(mesh->da, mesh->Nvert_o, INSERT_VALUES, mesh->lNvert_o);

    //search to update the ibm cells and support nodes if dynamic ibm
    if(ibm->dynamic)
    {
        // destroy the lists created if dynamic before next iteration
        if(ibm->access->clock->it > ibm->access->clock->itStart)
        {
            destroyLists(ibm);
        }

        //reset nvert values to 0, they will be recomputed during ibm search
        VecSet(mesh->Nvert,0.);
        VecSet(mesh->lNvert,0.);

        UpdateIBMesh(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        findBodyBoundingBox(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        findSearchCellDim(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        createSearchCellList(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        ibmSearch(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        findIBMFluidCells(ibm);

        MPI_Barrier(mesh->MESH_COMM);

        findClosestIBMElement(ibm);

        MPI_Barrier(mesh->MESH_COMM);

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
        CurvibInterpolation(ibm);

        MPI_Barrier(mesh->MESH_COMM);
    }

    //update the boundary condition around ibm cells
    UpdateImmersedBCs(ibm);

    contravariantToCartesian(ibm->access->ueqn);

    MPI_Barrier(mesh->MESH_COMM);

    PetscTime(&ibmTimeEnd);
    PetscPrintf(mesh->MESH_COMM, "ibm update time = %lf s\n", ibmTimeEnd - ibmTimeStart);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode setIBMWallModels(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    flags_        *flags = ibm->access->flags;

    word          fileNameU = "./boundary/" + mesh->meshName + "/U";

    PetscPrintf(mesh->MESH_COMM, "Reading IBM wall model...");

    for (PetscInt i=0; i < ibm->numBodies; i++)
    {
        char objectName[256];
        sprintf(objectName, "object%ld", i);

        ibmObject   *ibmBody = ibm->ibmBody[i];

        // read wall model
        readSubDictWord("./IBM/IBMProperties.dat", objectName, "wallModelU", &(ibmBody->wallModelU));

        // read wall model properties
        readSubDictWord("./IBM/IBMProperties.dat", objectName, "wallModelUProp", &(ibmBody->wallModelUProp));

        //copy the wall model properties from the U boundary condition
        if(ibmBody->wallModelUProp == "matchUiLeft")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.iLeft == "velocityWallFunction")
            {
                PetscInt wallModelType;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(wallModelType));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallmodel));

                if(wallModelType == -1)
                {
                    if(ibmBody->wallModelU == "Cabot")
                    {
                        PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallmodel->wmCabot));

                        Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));

                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U iLeft bc is not Cabot model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else if (wallModelType == -3)
                {
                    if(ibmBody->wallModelU == "Schumann")
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallmodel->wmShumann));

                        Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                        readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U iLeft bc is not Schumann model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else
            {
                char error[512];
                sprintf(error, "wall model prop = %s. But U iLeft bc is not a wall model\n", ibmBody->wallModelUProp.c_str());
                fatalErrorInFunction("setIBMWallModels", error);
            }
        }
        else if(ibmBody->wallModelUProp == "matchUiRight")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.iRight == "velocityWallFunction")
            {
                PetscInt wallModelType;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(wallModelType));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallmodel));

                if(wallModelType == -1)
                {
                    if(ibmBody->wallModelU == "Cabot")
                    {
                        PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallmodel->wmCabot));

                        Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));

                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U iRight bc is not Cabot model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else if (wallModelType == -3)
                {
                    if(ibmBody->wallModelU == "Schumann")
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallmodel->wmShumann));

                        Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                        readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U iRight bc is not Schumann model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else
            {
                char error[512];
                sprintf(error, "wall model prop = %s. But U iRight bc is not a wall model\n", ibmBody->wallModelUProp.c_str());
                fatalErrorInFunction("setIBMWallModels", error);
            }
        }
        else if(ibmBody->wallModelUProp == "matchUjLeft")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.jLeft == "velocityWallFunction")
            {
                PetscInt wallModelType;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(wallModelType));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallmodel));

                if(wallModelType == -1)
                {
                    if(ibmBody->wallModelU == "Cabot")
                    {
                        PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallmodel->wmCabot));

                        Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));

                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U jLeft bc is not Cabot model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else if (wallModelType == -3)
                {
                    if(ibmBody->wallModelU == "Schumann")
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallmodel->wmShumann));

                        Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                        readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U jLeft bc is not Schumann model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else
            {
                char error[512];
                sprintf(error, "wall model prop = %s. But U jLeft bc is not a wall model\n", ibmBody->wallModelUProp.c_str());
                fatalErrorInFunction("setIBMWallModels", error);
            }
        }
        else if(ibmBody->wallModelUProp == "matchUjRight")
        {
            //check that the iLeft has a wall model bc
            if (mesh->boundaryU.jRight == "velocityWallFunction")
            {
                PetscInt wallModelType;

                readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(wallModelType));

                // allocate memory and connect pointers
                PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallmodel));

                if(wallModelType == -1)
                {
                    if(ibmBody->wallModelU == "Cabot")
                    {
                        PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallmodel->wmCabot));

                        Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",   &(wm->kappa));

                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U jRight bc is not Cabot model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else if (wallModelType == -3)
                {
                    if(ibmBody->wallModelU == "Schumann")
                    {
                        PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallmodel->wmShumann));

                        Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                        readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                        readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                    }
                    else
                    {
                        char error[512];
                        sprintf(error, "wall model type = %s. But U jRight bc is not Schumann model\n", ibmBody->wallModelUProp.c_str());
                        fatalErrorInFunction("setIBMWallModels", error);
                    }
                }
                else
                {
                    char error[512];
                    sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                    fatalErrorInFunction("setIBMWallModels", error);
                }

            }
            else
            {
                char error[512];
                sprintf(error, "wall model prop = %s. But U jRight bc is not a wall model\n", ibmBody->wallModelUProp.c_str());
                fatalErrorInFunction("setIBMWallModels", error);
            }
        }
        else if (ibmBody->wallModelUProp == "setHere")
        {
            // allocate memory and connect pointers
            PetscMalloc(sizeof(wallModel), &(ibmBody->ibmWallmodel));

            if(ibmBody->wallModelU == "Cabot")
            {
                PetscMalloc(sizeof(Cabot), &(ibmBody->ibmWallmodel->wmCabot));

                Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "roughness", &(wm->roughness));
                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",  &(wm->kappa));

            }
            else if (ibmBody->wallModelU == "Schumann")
            {
                PetscMalloc(sizeof(Shumann), &(ibmBody->ibmWallmodel->wmShumann));

                Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                readSubDictWord  ("./IBM/IBMProperties.dat", objectName, "uStarEval", &(wm->wfEvalType));
                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kappa",     &(wm->kappa));
                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "thetaRef",  &(wm->thetaRef));
                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "kRough",    &(wm->roughness));
                readSubDictDouble("./IBM/IBMProperties.dat", objectName, "gammaM",    &(wm->gammaM));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option Cabot or Schumann \n");
                fatalErrorInFunction("setIBMWallModels", error);
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
    PetscInt           gxs   = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt           gys   = info.gys, gye = info.gys + info.gym;
    PetscInt           gzs   = info.gzs, gze = info.gzs + info.gzm;

    cellIds       initCp;

    PetscInt      i1, j1, k1;
    PetscInt      b, c, s, e;
    PetscReal     sc, sb, uDot, ustar;
    Cmpnts        ***ucat, ***cent, uDiff, ibmPtVel, bPtVel;
    Cmpnts        eN, eT1, eT2, dI, dF, eC;
    Cmpnts        checkPt, del;                                // point where the pressure or shear stress needs to be interpolated
    PetscReal     refLength, eA;

    PetscInt      n1, n2, n3;
    PetscReal     ***lP, ***nvert, ***aj, bPtP;

    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, peqn->lP, &lP);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lAj, &aj);

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    //vectors to store the net force on each ibm
    std::vector<Cmpnts>    lForce(ibm->numBodies);
    std::vector<Cmpnts>    gForce(ibm->numBodies);

    std::vector<Cmpnts>    lMoment(ibm->numBodies);
    std::vector<Cmpnts>    gMoment(ibm->numBodies);

    for (b = 0; b < ibm->numBodies; b++)
    {
        lForce[b] = nSetZero();
        gForce[b] = nSetZero();

        lMoment[b] = nSetZero();
        gMoment[b] = nSetZero();
    }

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        ibmObject   *ibmBody = ibm->ibmBody[b];
        ibmMesh     *ibMsh   = ibmBody->ibMsh;
        ibmRotation *ibmRot  = ibmBody->ibmRot;

        //check if process controls this ibm body
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
                    refLength = pow( 1./aj[ck][cj][ci], 1./3.);

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

                    // max distance checkpoint has to move to be to the closest fluid trilinear interpolation box
                    // restrict this to 3*ref length
                    PetscReal sumDel = 0.0;

                    // if there is an ibm solid cell in the interpolation, increase reference length further
                    while ( (intFlag == 1) && isInsideSimDomain(checkPt, mesh->bounds) && (sumDel <= 3*refLength))
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
                            del =  nScale(0.5 * refLength, eN);
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
                    }

                    //set the ibm element pressure force
                    PetscReal pForce = lElemPressure[e] * eA * cst->rho;
                    ibmBody->ibmPForce[e] = nSet(nScale(-pForce, eN));

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

                    //find the friction velocity ustar
                    uDiff = nSub(bPtVel, ibmPtVel);

                    // distance between interpolation point and ibm mesh element
                    sc    = nMag(nSub(checkPt, eC));

                    PetscReal  un, ut;              //normal and tangential component of udiff
                    Cmpnts     unV, utV;            // normal and tanential vector of udiff

                    // normal component to the wall
                    un = nDot(uDiff, eN);
                    unV = nScale(un,eN);

                    // tangential component to the wall
                    utV = nSub(uDiff, unV);
                    ut = nMag(utV);

                    //friction velocity
                    ustar = uTauCabot(cst->nu, ut, sc, 0.01, 0);

                    //shear force on ibm mesh element
                    Cmpnts tauWall = nScale(eA * cst->rho * ustar * ustar, nUnit(utV));

                    //sum total pressure and viscous force
                    mSum(lForce[b], ibmBody->ibmPForce[e]);
                    mSum(lForce[b], tauWall);

                    if(ibmBody->bodyMotion == "rotation")
                    {
                        // find moment for rotating mesh

                        //radius vector from the rotation center to the ibm mesh coordinate
                        Cmpnts rvec   = nSub(eC, ibmRot->rotCenter);

                        Cmpnts pressureMoment = nCross(rvec, ibmBody->ibmPForce[e]);
                        Cmpnts viscousMoment  = nCross(rvec, tauWall);

                        //net moment per processor
                        mSum(lMoment[b], nSum(pressureMoment, viscousMoment));

                    }


                }

            }

            MPI_Allreduce(&(lForce[b]), &(gForce[b]), 3, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);
            MPI_Allreduce(&(lMoment[b]), &(gMoment[b]), 3, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);
            MPI_Allreduce(lElemPressure, gElemPressure, ibMsh->elems, MPIU_REAL, MPIU_SUM, ibmBody->IBM_COMM);

            //compute power
            PetscReal ibmPower = 0.0;
            PetscReal netMoment = 0.0;

            if(ibmBody->bodyMotion == "rotation")
            {
                ibmPower  = nDot(gMoment[b], nScale(ibmRot->angSpeed, ibmRot->rotAxis));
                netMoment = nDot(gMoment[b], ibmRot->rotAxis);
            }

            //write data
            if(ibmrank == 0)
            {
                writeIBMForceData(ibm, b, gElemPressure, netMoment, ibmPower, gForce[b]);
            }

            delete[] gElemPressure;
            delete[] lElemPressure;
        }

    }

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(da, peqn->lP, &lP);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    MPI_Barrier(mesh->MESH_COMM);

    return 0;

}

//***************************************************************************************************************//
PetscErrorCode  writeIBMForceData(ibm_ *ibm, PetscInt b, PetscReal *gElemPressure, PetscReal netMoment, PetscReal ibmPower, Cmpnts netForce)
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

    PetscReal   epsilon          = 1e-8;

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
            word w5   = "Fx [N]";
            word w6   = "Fy [N]";
            word w7   = "Fz [N]";
            int width = -20;

            if(fp)
            {
                if(ibmBody->bodyMotion == "rotation")
                {
                    fprintf(fp, "%*s %*s %*s %*s %*s %*s %*s\n", width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(),  width, w5.c_str(),  width, w6.c_str(),  width, w7.c_str());
                }
                else if (ibmBody->bodyMotion == "static")
                {
                    fprintf(fp, "%*s %*s %*s %*s\n", width, w1.c_str(), width, w5.c_str(),  width, w6.c_str(),  width, w7.c_str());
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
            (clock->time - timeStart) / timeInterval -
            std::floor
            (
                (clock->time - timeStart) / timeInterval + epsilon
            ) < 1e-10
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
            clock->it / timeInterval -
            std::floor
            (
                clock->it / timeInterval
            ) < 1e-10
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
                fprintf(f, "%*.5f %*.5f %*.5f %*.5f %*.5f %*.5f %*.5f\n", width, clock->time, width, netMoment, width, ibmPower, width, ibmBody->ibmRot->angSpeed * 60/(2*M_PI), width, netForce.x,  width, netForce.y,  width, netForce.z);
            }
            else if (ibmBody->bodyMotion == "static")
            {
                fprintf(f, "%*.5f %*.5f %*.5f %*.5f\n", width, clock->time, width, netForce.x,  width, netForce.y,  width, netForce.z);
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
        }
    }

    return(0);
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

        ibMsh->nU[n] = nCross(angVel, rvec);
    }

    //update the angular speed if body is accelerating
    ibmRot->angSpeed += (ibmRot->angAcc * clock->dt);

    // recompute the ibm mesh properties - as rigid body motion the elements do not change
    PetscInt  n1, n2, n3;
    PetscReal normMag;
    Cmpnts    vec1, vec2, temp;

    for (PetscInt i=0; i<ibMsh->elems; i++)
    {
        // get the element nodes
        n1 = ibMsh->nID1[i]; n2 = ibMsh->nID2[i]; n3 = ibMsh->nID3[i];

        vec1 = nSub(ibMsh->nCoor[n2], ibMsh->nCoor[n1]);

        vec2 = nSub(ibMsh->nCoor[n3], ibMsh->nCoor[n1]);

        // normal to the face is found as cross product of the edges vec1 and vec2
        ibMsh->eN[i] = nCross(vec1, vec2);
        normMag = nMag(ibMsh->eN[i]);
        mScale(1.0/normMag, ibMsh->eN[i]);

        // flip normal if pointing inwards
        if(ibMsh->flipNormal[i] == 1)
        {
            mScale(-1, ibMsh->eN[i]);
        }

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

        //element area
        ibMsh->eA[i] = normMag/2.0;

        //element center
        temp = nSum(ibMsh->nCoor[n1], ibMsh->nCoor[n2]);
        ibMsh->eCent[i] = nSum( temp, ibMsh->nCoor[n3]);
        mScale(1/3.0, ibMsh->eCent[i]);
    }

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

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;

    // loop through the ibm bodies
    for(b = 0; b < ibm->numBodies; b++)
    {
            // smallest bounding sphere algorithm
        elementBoundingSphere(ibm->ibmBody[b]->ibMsh);
    }

    for(c = 0; c < ibm->numIBMFluid; c++)
    {
        //cell indices
        i = ibF[c].cellId.i;
        j = ibF[c].cellId.j;
        k = ibF[c].cellId.k;

        PetscInt      n1, n2, n3;
        PetscInt      cellMin = -100;
        Cmpnts        dis, elemNorm, p1, p2, p3;
        PetscReal     d_center, dmin = 1.0e20, d;
        Cmpnts        pmin, po, pj;
        PetscReal     normProj;                             // normal projection of point to ibm mesh element
        PetscInt      bodyID;

        // loop through the ibm bodies
        for (b = 0; b < ibm->numBodies; b++)
        {

            boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
            ibmMesh       *ibMsh = ibm->ibmBody[b]->ibMsh;                         // pointer to the ibm body mesh

            Cmpnts        *qvec   = ibMsh->eQVec;
            PetscReal     *rvec   = ibMsh->eRVec;

            PetscInt      *nv1   = ibMsh->nID1, *nv2  = ibMsh->nID2, *nv3  = ibMsh->nID3;
            Cmpnts        *eN    = ibMsh->eN;
            Cmpnts        *nCoor = ibMsh->nCoor;

            //loop through the IBM elements
            for(PetscInt e = 0; e < ibMsh->elems; e++)
            {

                dis = nSub(cent[k][j][i], qvec[e]);
                d_center = nMag(dis);

                // find the ibm mesh element whose bounding sphere center is closest to cent[k][j][i]
                if(d_center - rvec[e] < dmin)
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

                    // cent[k][j][i] located on the positive side of surface triangle
                    if(normProj >= 0)
                    {
                        dis = nScale(normProj, elemNorm);
                        pj  = nSub(cent[k][j][i], dis);

                        // The projected point is inside the triangle
                        if(isPointInTriangle(pj, p1, p2, p3, elemNorm) == 1)
                        {
                            if(normProj < dmin)
                            {
                                dmin    = normProj;
                                pmin    = pj;
                                cellMin = e;
                                bodyID  = b;
                            }
                        }
                        // The projected point is on the triangle line
                        else
                        {
                            // po is the projected point on the line and d is the distance to that point from cent[k][j][i]
                            disP2Line(cent[k][j][i], p1, p2, &po, &d);

                            if (d < dmin)
                            {
                                dmin = d;
                                pmin = po;
                                cellMin = e;
                                bodyID  = b;
                            }

                            disP2Line(cent[k][j][i], p2, p3, &po, &d);

                            if (d < dmin)
                            {
                                dmin = d;
                                pmin = po;
                                cellMin = e;
                                bodyID  = b;
                            }

                            disP2Line(cent[k][j][i], p3, p1, &po, &d);

                            if (d < dmin)
                            {
                                dmin = d;
                                pmin = po;
                                cellMin = e;
                                bodyID  = b;
                            }
                        }
                    }
                }

            }
        }

        if (cellMin == -100)
        {
           char error[512];
            sprintf(error, "Nearest Cell Searching Error for cell %ld %ld %ld\n", k, j, i);
            fatalErrorInFunction("findClosestIBMElement",  error);
        }

        //set the closest ibm mesh element to the ibm fluid cell
        ibF[c].closestElem = cellMin;
        ibF[c].pMin = pmin;
        ibF[c].minDist = dmin;
        ibF[c].bodyID = bodyID;

        Cpt2D     pjp, pj1, pj2, pj3;
        ibmMesh   *ibMsh = ibm->ibmBody[bodyID]->ibMsh;                         // pointer to the ibm body mesh

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

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CurvibInterpolation(ibm_ *ibm)
{
    mesh_         *mesh  = ibm->access->mesh;
    ueqn_         *ueqn  = ibm->access->ueqn;
    teqn_         *teqn  = ibm->access->teqn;
    flags_        *flags = ibm->access->flags;
    constants_    *cst   = ibm->access->constants;

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
    Cmpnts        ***ucat, ***lucat, ***cent;
    PetscReal     ***lt, ***temp, ***aj, ***nvert;
    Cmpnts        bPt, bPtInit;
    Cmpnts        bPtVel, ibmPtVel;
    PetscReal     nfMag, cellSize, bPtTemp, ibmPtTemp;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

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
        Cmpnts	   eNorm = ibMsh->eN[cElem];
        PetscInt      n1 = ibMsh->nID1[cElem], n2 = ibMsh->nID2[cElem], n3 = ibMsh->nID3[cElem];
        PetscReal  roughness = ibm->ibmBody[ibF[c].bodyID]->ibmWallmodel->wmCabot->roughness;
        // background mesh projection point is taken 1 cell distance along the normal directional of the closest ibm mesh element
        cellSize = pow( 1./aj[k][j][i], 1./3.);
        bPt = nScale(cellSize, eNorm);
        mSum(bPt, cent[k][j][i]);

        // search for the closest background mesh cell center (it will be close to current ibm fluid cell k,j,i)
        PetscInt    i1, j1, k1;
        PetscReal   dmin = 10e10, d;
        PetscInt    ic, jc, kc;

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
                }
            }
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
        // restrict this to 3*ref length
        PetscReal sumDel = 0.0;

        while ( (intFlag == 1) && isInsideSimDomain(bPt, mesh->bounds) && (sumDel <= 3*cellSize))
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

            }
        }

        // while loop fails
        if(intFlag > 0)
        {
            bPtVel = nSet(lucat[initCp.k][initCp.j][initCp.i]);

            if (flags->isTeqnActive)
            {
                bPtTemp = lt[initCp.k][initCp.j][initCp.i];
            }

            bPt = nSet(bPtInit);
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

         if (flags->isTeqnActive)
         {
             ibmPtTemp = 320;//ibm->access->abl->tRef;
         }

         // cabot wall model to interpolate the velocity at the ibm fluid node
         sb = nDot(nSub(cent[k][j][i], ibF[c].pMin), eNorm);
         sc = nDot(nSub(bPt, ibF[c].pMin), eNorm);

         // ucat[k][j][i].x = (sb/sc) * bPtVel.x + (1.0 - (sb/sc)) * ibmPtVel.x;
         // ucat[k][j][i].y= (sb/sc) * bPtVel.y + (1.0 - (sb/sc)) * ibmPtVel.y;
         // ucat[k][j][i].z = (sb/sc) * bPtVel.z + (1.0 - (sb/sc)) * ibmPtVel.z;

         if (roughness > 1.e-12)
         {
             wallFunctionCabotRoughness(cst->nu, roughness, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
         }
         else
         {
             wallFunctionCabot(cst->nu, sc, sb, ibmPtVel, bPtVel, &ucat[k][j][i], &ustar, eNorm);
         }

         if (flags->isTeqnActive)
         {
             temp[k][j][i] = (1.0 - (sb/sc)) * ibmPtTemp + (sb/sc) * bPtTemp;
         }

         findIBMWallShear(ibm, c, ustar);
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
        DMDAVecRestoreArray(da, teqn->Tmprt, &temp);
        DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    }

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    return(0);
}

//***************************************************************************************************************//
PetscErrorCode findIBMWallShear(ibm_ *ibm, PetscInt c, PetscReal ustar)
{
    mesh_         *mesh  = ibm->access->mesh;
    clock_        *clock = ibm->access->clock;
    constants_    *cst   = ibm->access->constants;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i, j, k;
    PetscInt      i1, j1, k1;

    Vec           Coor;
    Cmpnts        ***coor, faceCent;
    PetscReal     ***nvert;

    //local pointer for ibmFluidCells
    ibmFluidCell  *ibF = ibm->ibmFCells;

    ibmObject     *ibmBody = ibm->ibmBody[ibF[c].bodyID];
    ibmMesh       *ibMsh = ibmBody->ibMsh;
    PetscInt      cElem = ibF[c].closestElem;
    Cmpnts	      eNorm = ibMsh->eN[cElem];
    PetscReal     yFace, tau13;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // current ibm fluid cell index
    i = ibF[c].cellId.i;
    j = ibF[c].cellId.j;
    k = ibF[c].cellId.k;

    PetscInt      isFace = 0;
    //check the neighbourhood of current ibm fluid cell for fluid cells
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
            ) && isFluidCell(k1, j1, i1, nvert)
        )
        {
            //save the face center co-ordinate to interpolate
            if(k1 > k && j1 == j && i1 == i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k][j][i], nSum(coor[k][j][i-1], nSum(coor[k][j-1][i], coor[k][j-1][i-1])))));
                isFace = 1;
            }
            else if (k1 < k && j1 == j && i1 == i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k-1][j][i], nSum(coor[k-1][j][i-1], nSum(coor[k-1][j-1][i], coor[k-1][j-1][i-1])))));
                isFace = 1;
            }
            else if (k1 == k && j1 > j && i1 == i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k][j][i], nSum(coor[k][j][i-1], nSum(coor[k-1][j][i], coor[k-1][j][i-1])))));
                isFace = 1;
            }
            else if (k1 == k && j1 < j && i1 == i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k][j-1][i], nSum(coor[k][j-1][i-1], nSum(coor[k-1][j-1][i], coor[k-1][j-1][i-1])))));
                isFace = 1;
            }
            else if (k1 == k && j1 == j && i1 > i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k][j][i], nSum(coor[k-1][j][i], nSum(coor[k][j-1][i], coor[k-1][j-1][i])))));
                isFace = 1;
            }
            else if (k1 == k && j1 == j && i1 < i)
            {
                faceCent = nSet(nScale(1.0/4.0, nSum(coor[k][j][i-1], nSum(coor[k-1][j][i-1], nSum(coor[k][j-1][i-1], coor[k-1][j-1][i-1])))));
                isFace = 1;
            }

            // only for the cells that share a face, not diagonal cells
            if(isFace)
            {
                //find distance from the closest ibm mesh cell to this face
                yFace = nDot(nSub(faceCent, ibF[c].pMin), eNorm);

                if(ibmBody->wallModelU == "Cabot")
                {
                    Cabot *wm = ibmBody->ibmWallmodel->wmCabot;

                    tau13 = cst->nu * ustar / (wm->kappa * yFace);
                }
                else if(ibmBody->wallModelU == "Schumann")
                {
                    Shumann *wm = ibmBody->ibmWallmodel->wmShumann;

                    tau13 = cst->nu * ustar / (wm->kappa * yFace);
                }

                //find the shear stress tensor at this face

                //transform it to wall normal co-ordinate system

                //replace its tau13 component

                //transform back to the original co-ordinate system

                //save this tensor
            }

        }
    }

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

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

  PetscInt      i, j, k, b;
  PetscInt      ip, im, jp, jm, kp, km;
  PetscInt      ii, jj, kk;

  cellIds       sCell;

  PetscReal     ***nvert, ***nvert_o, ***gnvert;
  Cmpnts        ***cent;

  lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
  lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
  lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

  // loop through the ibm bodies
  for(b = 0; b < ibm->numBodies; b++)
  {

    boundingBox   *ibBox = ibm->ibmBody[b]->bound;                         // bounding box of the ibm body
    ibmMesh       *ibMsh = ibm->ibmBody[b]->ibMsh;                         // pointer to the ibm body mesh

    searchBox *sBox           = &(ibm->sBox[b]);
    list      *searchCellList = ibm->ibmBody[b]->searchCellList;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, mesh->Nvert, &nvert);
    DMDAVecGetArray(da, mesh->lNvert_o, &nvert_o);

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

           }
           else
            {   // set nvert to old value
                nvert[k][j][i] = nvert_o[k][j][i];

                if (PetscInt (nvert[k][j][i]+0.5) ==3)
                    nvert[k][j][i] = 4;             // only inside solid is set, here = 4
                if (PetscInt (nvert[k][j][i]+0.5) ==1)
                    nvert[k][j][i] = 0;             // near boundary values are set to 0 and need to be recalculated

           }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->Nvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lNvert_o, &nvert_o);

    MPI_Barrier(mesh->MESH_COMM);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    // set nvert at solid fluid intersection IB Nodes
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->Nvert, &gnvert);

    for (k = zs; k < ze; k++)
    for (j = ys; j < ye; j++)
    for (i = xs; i < xe; i++)
    {
        if (nvert[k][j][i] < 0)
        {
          gnvert[k][j][i] = 0;
        }

        ip = (i < mx - 1 ? (i + 1) : (i));
        im = (i > 0 ? (i - 1) : (i));

        jp = (j < my - 1 ? (j + 1) : (j));
        jm = (j > 0 ? (j - 1) : (j));

        kp = (k < mz - 1 ? (k + 1) : (k));
        km = (k > 0 ? (k - 1) : (k));

        //if not a solid node
        if ((PetscInt) (nvert[k][j][i] + 0.5) != 4)
        {
            for (kk = km; kk < kp + 1; kk++)
            for (jj = jm; jj < jp + 1; jj++)
            for (ii = im; ii < ip + 1; ii++)
            {
                if ((PetscInt) (nvert[kk][jj][ii] + 0.5) == 4)
                { //if next to a solid node
                    gnvert[k][j][i] = PetscMax(2, nvert[k][j][i]);
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->Nvert, &gnvert);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    // nvert values of 4 for solid and 2 for IB fluid nodes where used for the current body
    // to differentiate it from other bodies.
    // reset back to nvert values of 3 for solid and 1 for IB fluid
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->Nvert, &gnvert);

    // Back to the old nvert 3 and 1
    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
        if ((PetscInt) (nvert[k][j][i] + 0.5) == 2)
        {
            gnvert[k][j][i] = 1;
        }
        if ((PetscInt) (nvert[k][j][i] + 0.5) == 4)
        {
            gnvert[k][j][i] = 3;
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->Nvert, &gnvert);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    //ibm nvert cleanup - make solid ibm fluid cells that dont have even one fluid cell around it
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->Nvert, &gnvert);

    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
        PetscInt ctr = 0;

        ip = (i < mx - 1 ? (i + 1) : (i));
        im = (i > 0 ? (i - 1) : (i));

        jp = (j < my - 1 ? (j + 1) : (j));
        jm = (j > 0 ? (j - 1) : (j));

        kp = (k < mz - 1 ? (k + 1) : (k));
        km = (k > 0 ? (k - 1) : (k));

        if (isIBMFluidCell(k, j, i, gnvert))
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
                gnvert[k][j][i] = 3.0;
            }
        }

    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->Nvert, &gnvert);

    DMGlobalToLocalBegin(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);
    DMGlobalToLocalEnd(da, mesh->Nvert, INSERT_VALUES, mesh->lNvert);

    MPI_Barrier(mesh->MESH_COMM);

  }

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
  PetscReal     A[4][4]={0.}, inv_A[4][4]={0.};
  PetscInt      n = 0, nsupport=0;

  PetscReal     dis;
  PetscReal     *W, *PHI, **P, **B, *Ux, *Uy, *Uz, *tmp;

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

    inv_4by4(A, inv_A);

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
    PetscReal     ***nvert;
    PetscInt      lsum = 0., sum = 0.;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {

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
              sprintf(error, "IBMProperties file indicates that mesh [%s] is linked to ibm body [%s]. But it is not found in the mesh. Check ibm body base location or ensure that the correct nvert field has been read or mesh is too coarse.\n", mesh->meshName.c_str(), ibm->ibmBody[b]->bodyName.c_str());
              fatalErrorInFunction("checkIBMexists", error);
          }

    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

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
            PetscPrintf(PETSC_COMM_WORLD, "Checking IBM element normal direction for body: %s...", ibm->ibmBody[b]->bodyName.c_str());
            // set offset distance
            minBound = PetscMin( PetscMin(ibBox->Lx, ibBox->Ly), ibBox->Lz);
            offset   = 0.001 * minBound;

            // check that the normal points outwards
            for (e=0; e<ibMesh->elems; e++)
            {
                if(e%1000 == 0 && ibm->dbg)
                PetscPrintf(PETSC_COMM_WORLD, "element = %ld\n", e);

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
                    //normal points inverts, reverse the direction
                    mScale(-1.0, ibMesh->eN[e]);
                    ibMesh->flipNormal[e] = 1;
                }

            }

            PetscPrintf(PETSC_COMM_WORLD, "done\n");
        }

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

                 // eT2 = eT2 x eN
                ibMesh->eT2[e].x = -ibMesh->eN[e].x*ibMesh->eN[e].z/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT2[e].y = -ibMesh->eN[e].y*ibMesh->eN[e].z/ sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
                ibMesh->eT2[e].z = sqrt(ibMesh->eN[e].x*ibMesh->eN[e].x + ibMesh->eN[e].y*ibMesh->eN[e].y);
            }

        }

        if(ibm->writeSTL)
        {
            {
                writeSTLFile(ibm->ibmBody[b]);
                MPI_Barrier(ibm->access->mesh->MESH_COMM);
            }
        }

    }
    return (0);
}
//***************************************************************************************************************//

PetscErrorCode findIBMControlledProcs(ibm_ *ibm)
{
    PetscInt      b, c;
    mesh_         *mesh = ibm->access->mesh;

    PetscMPIInt      rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    //local pointer for ibmFluidCells
    ibmFluidCell *ibF = ibm->ibmFCells;


    // loop through the ibm bodies
    for (b = 0; b < ibm->numBodies; b++)
    {
        // define communicator color
        PetscInt    commColor = MPI_UNDEFINED;
        ibmObject   *ibmBody  = ibm->ibmBody[b];

        if(ibm->numIBMFluid!=0)
        {
            // loop through the ibm fluid cells
            for(c = 0; c < ibm->numIBMFluid; c++)
            {
                if( b == ibF[c].bodyID)
                {
                    // set this ibm body to be controlled by this processor
                    ibmBody->ibmControlled = 1;

                    // set the communicator color
                    commColor = 1;

                    //atleast one ibm fluid cell close to this body found for this processor, break out of loop
                    break;
                }
                else
                {
                    ibmBody->ibmControlled = 0;
                }
            }
        }
        else
        {
            ibmBody->ibmControlled = 0;
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

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode initFindIBMElementProcs(ibm_ *ibm)
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

                //loop through the ibm fluid cells
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

            for(e = 0; e < ibMsh->elems; e++)
            {
                // point is controlled
                if(lminDist[e] == gminDist[e])
                {
                    ibmBody->thisPtControlled[e] = 1;
                }
                // point is not controlled
                else
                {
                    ibmBody->thisPtControlled[e] = 0;
                }

                //initialize flag for ibm element processor flag to 0
                ibmBody->thisPtControlTransfer[e] = 0;
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

PetscErrorCode checkIBMElementControlledProcessor(ibm_ *ibm)
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

        }

    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    MPI_Barrier(mesh->MESH_COMM);

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode elementBoundingSphere(ibmMesh *ibMesh)
{
    PetscInt    n1, n2, n3;
    Cmpnts      p1, p2, p3, tmp1, tmp2, tmp3;
    Cmpnts      pa, pb, pc, pf, pu, pv, pd, pt;
    PetscReal   l12, l23, l31, lu, lv;
    PetscReal   gama, lamda;

    Cmpnts      *qvec = ibMesh->eQVec;
    PetscReal   *rvec = ibMesh->eRVec;

    for (PetscInt i=0; i<ibMesh->elems; i++)
    {
        // get the element nodes
        n1 = ibMesh->nID1[i]; n2 = ibMesh->nID2[i]; n3 = ibMesh->nID3[i];

        p1 = ibMesh->nCoor[n1]; p2 = ibMesh->nCoor[n2]; p3 = ibMesh->nCoor[n3];

        tmp1 = nSum(p1, p2);
        tmp2 = nSum(p2, p3);
        tmp3 = nSum(p3, p1);

        l12 = nMag(tmp1);
        l23 = nMag(tmp2);
        l31 = nMag(tmp3);

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

        tmp1 = nSub(qvec[i], pa);
        rvec[i] = nMag(tmp1);
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
      lcellSize += pow( 1./aj[k][j][i], 1./3.);
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

        // //set a default value
        // sBox->ncx = 1; sBox->ncy = 1; sBox->ncz = 1;

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

            if(e%1000 == 0 && ibm->dbg)
                PetscPrintf(PETSC_COMM_WORLD, "element = %ld\n", e);

            // min and max coordinate value of the element vertices
            xv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].x, ibMsh->nCoor[n2].x ), ibMsh->nCoor[n3].x );
            xv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].x, ibMsh->nCoor[n2].x ), ibMsh->nCoor[n3].x );

            yv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].y, ibMsh->nCoor[n2].y ), ibMsh->nCoor[n3].y );
            yv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].y, ibMsh->nCoor[n2].y ), ibMsh->nCoor[n3].y );

            zv_min = PetscMin( PetscMin( ibMsh->nCoor[n1].z, ibMsh->nCoor[n2].z ), ibMsh->nCoor[n3].z );
            zv_max = PetscMax( PetscMax( ibMsh->nCoor[n1].z, ibMsh->nCoor[n2].z ), ibMsh->nCoor[n3].z );

            // min and max index of the search cell where this element is located
            iv_min = floor((xv_min - ibBox->xmin) / sBox->dcx);          // find the search cell index of xmin coordinate of the element
            iv_max = floor((xv_max - ibBox->xmin) / sBox->dcx) + 1;      // +1 ensures that iv_min and iv_max are not the same

            jv_min = floor((yv_min - ibBox->ymin) / sBox->dcy); //
            jv_max = floor((yv_max - ibBox->ymin) / sBox->dcy) + 1;

            kv_min = floor((zv_min - ibBox->zmin) / sBox->dcz); //
            kv_max = floor((zv_max - ibBox->zmin) / sBox->dcz) + 1;

            //ensure that the search cell indices of the current element are bounded between 0 and max index
            iv_min = (iv_min < 0) ? 0 : iv_min;
            iv_max = (iv_max > sBox->ncx) ? sBox->ncx : iv_max;

            jv_min = (jv_min < 0) ? 0 : jv_min;
            jv_max = (jv_max > sBox->ncx) ? sBox->ncy : jv_max;

            kv_min = (kv_min < 0) ? 0 : kv_min;
            kv_max = (kv_max > sBox->ncz) ? sBox->ncz : kv_max;

            //insert element into search cell
            for (k = kv_min; k < kv_max; k++) {
                for (j = jv_min; j < jv_max; j++) {
                    for (i = iv_min; i < iv_max; i++) {
                        insertnode(&(searchCellList[k * sBox->ncx * sBox->ncy + j * sBox->ncx + i]), e);
                    }
                }
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

inline PetscInt isPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts norm)
{
  PetscInt   flag;
  Cpt2D	     pj, pj1, pj2, pj3;

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

  flag = ISInsideTriangle2D(pj, pj1, pj2, pj3);

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

inline void disP2Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d)
{
  PetscReal	dmin;
  PetscReal	dx21, dy21, dz21, dx31, dy31, dz31, t;

  dx21 = p2.x - p1.x; dy21 = p2.y - p1.y; dz21 = p2.z - p1.z;
  dx31 = p.x  - p1.x; dy31 = p.y  - p1.y; dz31 = p.z  - p1.z;

  t = (dx31 * dx21 + dy31 * dy21 + dz31 * dz21) / (dx21*dx21 + dy21*dy21 +
       dz21 * dz21);

  if (t<0)
  { // The closet point is p1
    po->x = p1.x; po->y = p1.y; po->z = p1.z;
    *d = sqrt(dx31*dx31 + dy31*dy31 + dz31*dz31);
  }
  else if (t>1)
  { // The closet point is p2
    po->x = p2.x; po->y = p2.y; po->z = p2.z;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
         (p.z - po->z) * (p.z - po->z));
  }
  else
  { // The closet point lies between p1 & p2
    po->x = p1.x + t * dx21; po->y = p1.y + t * dy21; po->z = p1.z + t*dz21;
    *d = sqrt((p.x - po->x)*(p.x - po->x)+(p.y - po->y) * (p.y - po->y) +
         (p.z - po->z) * (p.z - po->z));
  }

  return;
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
