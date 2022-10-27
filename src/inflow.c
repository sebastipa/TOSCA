//! \file  initialization.c
//! \brief Contains inflow boundary condition function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/boundary.h"
#include "include/inflow.h"

//***************************************************************************************************************//

PetscErrorCode SetInflowFunctions(mesh_ *mesh)
{
    PetscInt inletFunctionsAllocated = 0;

    // read inlet functions subdictionaries (if present)
    if(mesh->boundaryU.kLeft == "inletFunction")
    {
        word           fileName, location, field;

        PetscPrintf(mesh->MESH_COMM, "Creating U inflow boundary data...");

        // allocate memory for this patch inlet function data and init types
        PetscMalloc(sizeof(inletFunctionTypes), &mesh->inletF.kLeft);
        mesh->inletF.kLeft->typeU   = -1;
        mesh->inletF.kLeft->typeT   = -1;
        mesh->inletF.kLeft->typeNut = -1;
        inletFunctionsAllocated ++;

        // set local pointer to this inlet function type
        inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

        location = "./boundary/" + mesh->meshName + "/";
        field    = "U";
        fileName = location + field;

        readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeU));

        // power law profile
        if (ifPtr->typeU == 1)
        {
            readSubDictVector(fileName.c_str(), "inletFunction", "Uref", &(ifPtr->Uref));
            readSubDictDouble(fileName.c_str(), "inletFunction", "Href", &(ifPtr->Href));
            readSubDictDouble(fileName.c_str(), "inletFunction", "uPrimeRMS", &(ifPtr->uPrimeRMS));
        }
        // log law profile
        else if (ifPtr->typeU == 2)
        {
            readSubDictVector(fileName.c_str(), "inletFunction", "directionU", &(ifPtr->Udir));
            readSubDictDouble(fileName.c_str(), "inletFunction", "hInversion", &(ifPtr->hInv));
            readSubDictDouble(fileName.c_str(), "inletFunction", "frictionU",  &(ifPtr->uTau));
            readSubDictDouble(fileName.c_str(), "inletFunction", "kRough",     &(ifPtr->roughness));
            mScale(1.0/nMag(ifPtr->Udir), ifPtr->Udir);
        }
        // unsteady mapped inflow
        else if (ifPtr->typeU == 3)
        {
            readSubDictInt(fileName.c_str(), "inletFunction", "n1Inflow",  &(ifPtr->n1));
            readSubDictInt(fileName.c_str(), "inletFunction", "n2Inflow",  &(ifPtr->n2));
            readSubDictInt(fileName.c_str(), "inletFunction", "n1Periods", &(ifPtr->prds1));
            readSubDictInt(fileName.c_str(), "inletFunction", "n2Periods", &(ifPtr->prds2));

            // increase n1 and n2 accounting for side ghost cells
            ifPtr->n1wg = ifPtr->n1 + 2;
            ifPtr->n2wg = ifPtr->n2 + 2;

            // flags for T and nut
            ifPtr->mapT = 0;
            ifPtr->mapNut = 0;

            // print mapping action information
            printInflowMappingAction(mesh, ifPtr);

            // check also temperature specification for this inlet function
            if(mesh->access->flags->isTeqnActive)
            {
                if(mesh->boundaryT.kLeft == "inletFunction")
                {
                    field = "T";
                    fileName = location + field;
                    readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeT));

                    if(ifPtr->typeT == 3)
                    {
                        if(ifPtr->typeT != ifPtr->typeU)
                        {
                            // throw error
                            char error[512];
                            sprintf(error, "unsteadyMappedInflow must be applied to U, T if temperatureTransport is active\n");
                            fatalErrorInFunction("SetInflowFunctions", error);
                        }
                        else
                        {
                            ifPtr->mapT = 1;

                            // read parameters and check that are the same as U
                            PetscInt n1, n2, prds1, prds2;
                            readSubDictInt(fileName.c_str(), "inletFunction", "n1Inflow",  &n1);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n2Inflow",  &n2);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n1Periods", &prds1);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n2Periods", &prds2);

                            if(n1!=ifPtr->n1 || n2 != ifPtr->n2 || prds1!=ifPtr->prds1 || prds2!=ifPtr->prds2)
                            {
                                char error[512];
                                sprintf(error, "inletFunction type 3 parameters in boundary/T must match boundary/U");
                                fatalErrorInFunction("SetInflowFunctions",  error);
                            }
                        }
                    }
                }
            }

            // check also nut specifiation for this inlet function
            if(mesh->access->flags->isLesActive)
            {
                if(mesh->boundaryNut.kLeft == "inletFunction")
                {
                    field = "nut";
                    fileName = location + field;
                    readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeNut));

                    if(ifPtr->typeNut == 3)
                    {
                        if(ifPtr->typeNut == ifPtr->typeU)
                        {
                            ifPtr->mapNut = 1;

                            // read parameters and check that are the same as U
                            PetscInt n1, n2, prds1, prds2;
                            readSubDictInt(fileName.c_str(), "inletFunction", "n1Inflow",  &n1);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n2Inflow",  &n2);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n1Periods", &prds1);
                            readSubDictInt(fileName.c_str(), "inletFunction", "n2Periods", &prds2);
                            if(n1!=ifPtr->n1 || n2 != ifPtr->n2 || prds1!=ifPtr->prds1 || prds2!=ifPtr->prds2)
                            {
                                char error[512];
                                sprintf(error, "inletFunction type 3 parameters in boundary/nut must match boundary/U");
                                fatalErrorInFunction("SetInflowFunctions",  error);
                            }
                        }
                    }
                }
            }

            // initialize inflow data
            mappedInflowInitialize(ifPtr);
        }

        // unsteady interpolated inflow
        else if (ifPtr->typeU == 4)
        {
            field = "U";
            fileName = location + field;

            // read if source mesh is uniform or grading
            readSubDictWord(fileName.c_str(), "inletFunction", "sourceType",  &(ifPtr->sourceType));

            // read number of cells and periods
            readSubDictInt   (fileName.c_str(), "inletFunction", "n1Inflow",   &(ifPtr->n1));
            readSubDictInt   (fileName.c_str(), "inletFunction", "n2Inflow",   &(ifPtr->n2));
            readSubDictInt   (fileName.c_str(), "inletFunction", "n1Periods",  &(ifPtr->prds1));
            readSubDictInt   (fileName.c_str(), "inletFunction", "n2Periods",  &(ifPtr->prds2));

            if(ifPtr->sourceType == "uniform")
            {
                readSubDictDouble(fileName.c_str(), "inletFunction", "cellWidth1", &(ifPtr->width1));
                readSubDictDouble(fileName.c_str(), "inletFunction", "cellWidth2", &(ifPtr->width2));

                ifPtr->inflowHeigth = ifPtr->n1*ifPtr->prds1*ifPtr->width1;
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

                    for (PetscInt k = 0; k < npx; k++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
                    for (PetscInt i = 0; i < npy; i++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
                    for (PetscInt j = 0; j < npz; j++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &Zcart[j]);

                    fclose(meshFileID);

                    // height of the inflow database
                    ifPtr->inflowHeigth = Zcart[npz-1] - Zcart[0];

                    // wipe vectors
                    std::vector<PetscReal> ().swap(Zcart);
                }
            }
            else
            {
                char error[512];
                sprintf(error, "unknown sourceType in inletFunction type 4, available types are\n    1: uniform\n    2: grading\n");
                fatalErrorInFunction("SetInflowFunctions",  error);
            }

            // increase n1 and n2 accounting for side ghost cells
            ifPtr->n1wg = ifPtr->n1 + 2;
            ifPtr->n2wg = ifPtr->n2 + 2;

            // flags for T and nut
            ifPtr->mapT   = 0;
            ifPtr->mapNut = 0;

            // check also temperature specification for this inlet function
            if(mesh->access->flags->isTeqnActive)
            {
                if(mesh->boundaryT.kLeft == "inletFunction")
                {
                    field = "T";
                    fileName = location + field;

                    readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeT));

                    if(ifPtr->typeT == 4)
                    {
                        ifPtr->mapT = 1;

                        // read parameters and check that are the same as U
                        PetscInt  n1, n2, prds1, prds2;
                        PetscReal width1, width2;
                        word      sourceTypeT;
                        readSubDictWord(fileName.c_str(),   "inletFunction", "sourceType", &sourceTypeT);
                        readSubDictInt(fileName.c_str(),    "inletFunction", "n1Inflow",   &n1);
                        readSubDictInt(fileName.c_str(),    "inletFunction", "n2Inflow",   &n2);
                        readSubDictInt(fileName.c_str(),    "inletFunction", "n1Periods",  &prds1);
                        readSubDictInt(fileName.c_str(),    "inletFunction", "n2Periods",  &prds2);
                        if
                        (
                            n1!=ifPtr->n1 || n2 != ifPtr->n2 ||
                            prds1!=ifPtr->prds1 || prds2!=ifPtr->prds2 ||
                            sourceTypeT!=ifPtr->sourceType
                        )
                        {
                            char error[512];
                            sprintf(error, "inletFunction type 4 parameters in boundary/T must match boundary/U");
                            fatalErrorInFunction("SetInflowFunctions",  error);
                        }
                    }
                    else
                    {
                        // throw error
                        char error[512];
                        sprintf(error, "unsteadyMappedInflow must be applied to U, T if temperatureTransport is active\n");
                        fatalErrorInFunction("SetInflowFunctions", error);
                    }
                }
            }

            // check also nut specifiation for this inlet function
            if(mesh->access->flags->isLesActive)
            {
                if(mesh->boundaryNut.kLeft == "inletFunction")
                {
                    field = "nut";
                    fileName = location + field;

                    readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeNut));

                    // if only nut is mapped throw error
                    if(ifPtr->typeNut == 4)
                    {
                        if(ifPtr->typeNut == ifPtr->typeU)
                        {
                            ifPtr->mapNut = 1;

                            // read parameters and check that are the same as U
                            PetscInt n1, n2, prds1, prds2;
                            PetscReal width1, width2;
                            word      sourceTypeNut;
                            readSubDictWord(fileName.c_str(),   "inletFunction", "sourceType", &sourceTypeNut);
                            readSubDictInt(fileName.c_str(),    "inletFunction", "n1Inflow",   &n1);
                            readSubDictInt(fileName.c_str(),    "inletFunction", "n2Inflow",   &n2);
                            readSubDictInt(fileName.c_str(),    "inletFunction", "n1Periods",  &prds1);
                            readSubDictInt(fileName.c_str(),    "inletFunction", "n2Periods",  &prds2);
                            if
                            (
                                n1!=ifPtr->n1 || n2 != ifPtr->n2 ||
                                prds1!=ifPtr->prds1 || prds2!=ifPtr->prds2 ||
                                sourceTypeNut!=ifPtr->sourceType
                            )
                            {
                                char error[512];
                                sprintf(error, "inletFunction type 4 parameters in boundary/nut must match boundary/U");
                                fatalErrorInFunction("SetInflowFunctions",  error);
                            }
                        }
                    }
                    else
                    {
                        // throw error
                        char error[512];
                        sprintf(error, "unsteadyMappedInflow must be applied also on nut\n");
                        fatalErrorInFunction("SetInflowFunctions", error);
                    }
                }
            }

            // build interpolation weights and find periodized inflow cells indices
            SetInflowWeights(mesh, ifPtr);

            // initialize inflow data
            mappedInflowInitialize(ifPtr);
        }
        // Nieuwstadt
        else if (ifPtr->typeU == 5)
        {
            readSubDictVector(fileName.c_str(), "inletFunction", "directionU", &(ifPtr->Udir));
            readSubDictDouble(fileName.c_str(), "inletFunction", "hInversion", &(ifPtr->hInv));
            readSubDictDouble(fileName.c_str(), "inletFunction", "hReference", &(ifPtr->Href));
            readSubDictDouble(fileName.c_str(), "inletFunction", "frictionU",  &(ifPtr->uTau));
            readSubDictDouble(fileName.c_str(), "inletFunction", "kRough",     &(ifPtr->roughness));
            readSubDictDouble(fileName.c_str(), "inletFunction", "latitude",   &(ifPtr->latitude));
            mScale(1.0/nMag(ifPtr->Udir), ifPtr->Udir);
            ifPtr->fc =  2*7.292115e-5*std::sin(ifPtr->latitude * M_PI / 180.0);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown inflow profile on k-left boundary, available profiles are:\n        1 : power law (alpha = 0.107027)\n        2 : log law according to ABLProperties.dat\n        3 : unsteady mapped inflow from database\n        4 : unsteady interpolated inflow from database\n        5 : Nieuwstadt inflow (with veer)");
            fatalErrorInFunction("SetInflowFunctions",  error);
        }

        PetscPrintf(mesh->MESH_COMM, "done\n\n");
    }

    // check also temperature specification
    if(mesh->access->flags->isTeqnActive)
    {
        if(mesh->boundaryT.kLeft == "inletFunction")
        {
            word           fileName, location, field;

            PetscPrintf(mesh->MESH_COMM, "Creating T inflow boundary data...");

            if(!inletFunctionsAllocated)
            {
                PetscMalloc(sizeof(inletFunctionTypes), &(mesh->inletF.kLeft));
                mesh->inletF.kLeft->typeU   = -1;
                mesh->inletF.kLeft->typeT   = -1;
                mesh->inletF.kLeft->typeNut = -1;
                inletFunctionsAllocated++;
            }

            // set local pointer to this inlet function type
            inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

            location = "./boundary/" + mesh->meshName + "/";
            field    = "T";
            fileName = location + field;

            readSubDictInt(fileName.c_str(), "inletFunction", "type", &(ifPtr->typeT));

            // RampanelliZardi stratification model
            if (ifPtr->typeT == 2)
            {
                PetscReal hInv;
                readSubDictDouble(fileName.c_str(), "inletFunction", "hInv",  &hInv);
                readSubDictDouble(fileName.c_str(), "inletFunction", "dInv",  &(ifPtr->dInv));
                readSubDictDouble(fileName.c_str(), "inletFunction", "gInv",  &(ifPtr->gInv));
                readSubDictDouble(fileName.c_str(), "inletFunction", "tRef",  &(ifPtr->tRef));
                readSubDictDouble(fileName.c_str(), "inletFunction", "gTop",  &(ifPtr->gTop));
                readSubDictDouble(fileName.c_str(), "inletFunction", "smearT",&(ifPtr->smear));

                // check input
                if((ifPtr->typeU == 2 || ifPtr->typeU == 5) && hInv != ifPtr->hInv)
                {
                    char error[512];
                    sprintf(error, "inletFunction type 2 requires the same hInv if active in both U and T");
                    fatalErrorInFunction("SetInflowFunctions",  error);
                }
            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");
        }
    }

    MPI_Barrier(mesh->MESH_COMM);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInflowWeights(mesh_ *mesh, inletFunctionTypes *ifPtr)
{
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs  = info.xs, xe = info.xs + info.xm;
    PetscInt      ys  = info.ys, ye = info.ys + info.ym;
    PetscInt      zs  = info.zs, ze = info.zs + info.zm;
    PetscInt      mx  = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    // declare arrays
    Cmpnts        ***cent;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    // allocate space for weights and closest cell ids
    ifPtr->closestCells  = (cellIds***)malloc(my*sizeof(cellIds**));
    ifPtr->inflowWeights = (PetscReal***)malloc(my*sizeof(PetscReal**));

    for(j=0; j<my; j++)
    {
        ifPtr->closestCells[j]  = (cellIds**)malloc(mx*sizeof(cellIds*));
        ifPtr->inflowWeights[j] = (PetscReal**)malloc(mx*sizeof(PetscReal*));

        for (i=0; i<mx; i++)
        {
            ifPtr->closestCells[j][i]  = (cellIds*)malloc(4*sizeof(cellIds));
            ifPtr->inflowWeights[j][i] = (PetscReal*)malloc(4*sizeof(PetscReal));

            /*
            ifPtr->closestCells[j][i][0].i = ifPtr->closestCells[j][i][0].j = ifPtr->closestCells[j][i][0].k = 0;
            ifPtr->closestCells[j][i][1].i = ifPtr->closestCells[j][i][1].j = ifPtr->closestCells[j][i][1].k = 0;
            ifPtr->closestCells[j][i][2].i = ifPtr->closestCells[j][i][2].j = ifPtr->closestCells[j][i][2].k = 0;
            ifPtr->closestCells[j][i][3].i = ifPtr->closestCells[j][i][3].j = ifPtr->closestCells[j][i][3].k = 0;
            */
        }
    }

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // due to marblles convention on block structured meshes:
    // i = y direction
    // j = z direction
    // k = x direction (useles for kLeftPatch)

    PetscInt jN;
    PetscInt iN;

    // create fictitious inflow plane (it has the same size of the actual mesh inflow)
    std::vector<std::vector<Cmpnts>> inflowCellCenters;

    // uniform mesh: just fill the plane using provided widths
    if(ifPtr->sourceType == "uniform")
    {
        // add one point because cell indexing starts from 1, ghost cells are
        // unnecessary, we only need the first for indexing the loops
        jN = std::floor(mesh->bounds.Lz / ifPtr->width1)+1;
        iN = std::floor(mesh->bounds.Ly / ifPtr->width2)+1;

        inflowCellCenters.resize(jN);

        for(PetscInt jif=1; jif<jN; jif++)
        {
            inflowCellCenters[jif].resize(iN);

            for(PetscInt iif=1; iif<iN; iif++)
            {
                inflowCellCenters[jif][iif].x = 0.0;
                inflowCellCenters[jif][iif].y = mesh->bounds.ymin + 0.5 * ifPtr->width2 + (iif-1) * ifPtr->width2;
                inflowCellCenters[jif][iif].z = mesh->bounds.zmin + 0.5 * ifPtr->width1 + (jif-1) * ifPtr->width1;
            }
        }
    }
    else if(ifPtr->sourceType == "grading")
    {
        // we have to actually read the points, numbering goes on as TOSCA indexing
        // so 1 = j, 2 = i for k-normal slices.

        std::vector<PetscReal>  Ycart, Zcart;
        std::vector<PetscReal>  points1(ifPtr->n1), points2(ifPtr->n2);

        word pointsFileName     = "./inflowDatabase/inflowMesh.xyz";
        FILE *meshFileID        = fopen(pointsFileName.c_str(), "r");

        PetscPrintf(mesh->MESH_COMM, "   -> reading source mesh %s\n",pointsFileName.c_str());

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

            Ycart.resize(npy);
            Zcart.resize(npz);

            for (k = 0; k < npx; k++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble, &bufferDouble);
            for (i = 0; i < npy; i++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble,  &Ycart[i], &bufferDouble);
            for (j = 0; j < npz; j++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble,  &Zcart[j]);

            fclose(meshFileID);

            // compute cell centers recentered so that the slice has one corner at zero
            for (j = 0; j < ifPtr->n1; j++) points1[j] = 0.5*(Zcart[j+1] + Zcart[j]) - Zcart[0];
            for (i = 0; i < ifPtr->n2; i++) points2[i] = 0.5*(Ycart[i+1] + Ycart[i]) - Ycart[0];

            // now compute fictitious mesh
            PetscReal Lyf,
                      Lzf;
            PetscReal Ly_inflow = Ycart[npy-1] - Ycart[0],
                      Lz_inflow = Zcart[npz-1] - Zcart[0];
            PetscInt  nTimes1 = 0, nTimes2 = 0;

            do{nTimes2++; Lyf = nTimes2*Ly_inflow;} while(Lyf < mesh->bounds.Ly);
            do{nTimes1++; Lzf = nTimes1*Lz_inflow;} while(Lzf < mesh->bounds.Lz);

            PetscPrintf(mesh->MESH_COMM, "   -> source mesh multipliers dir1 = %ld, dir2 = %ld\n",nTimes1, nTimes2);

            jN = nTimes1 * ifPtr->n1 + 1;
            iN = nTimes2 * ifPtr->n2 + 1;

            inflowCellCenters.resize(jN);

            for(PetscInt jif=1; jif<jN; jif++)
            {
                inflowCellCenters[jif].resize(iN);

                for(PetscInt iif=1; iif<iN; iif++)
                {
                    PetscInt sourceI = (iif-1) % ifPtr->n2,
                             sourceJ = (jif-1) % ifPtr->n1;
                    inflowCellCenters[jif][iif].x = 0.0;
                    inflowCellCenters[jif][iif].y = mesh->bounds.ymin + points2[sourceI] + std::floor((iif-1) / ifPtr->n2)*Ly_inflow;
                    inflowCellCenters[jif][iif].z = mesh->bounds.zmin + points1[sourceJ] + std::floor((jif-1) / ifPtr->n1)*Lz_inflow;
                }
            }

            // wipe vectors
            std::vector<PetscReal> ().swap(Ycart);
            std::vector<PetscReal> ().swap(Zcart);
            std::vector<PetscReal> ().swap(points2);
            std::vector<PetscReal> ().swap(points1);
        }
    }

    PetscPrintf(mesh->MESH_COMM, "   -> calculating inflow bilinear interpolation weigths...");
    // now do the search and, for each cell center on the kLeftPatch, find its
    // 4 closest cell centers belonging to the fictitious inflow plane

    // the k is free so all processors in the k line build the data

    // loop over internal cells on kLeftPatch
    for(j=lys; j<lye; j++)
    {
        for(i=lxs; i<lxe; i++)
        {
            PetscReal  minDistMag = 1e20;
            cellIds closestCell;

            // loop over cells on fictitious inflow plane
            for(PetscInt jif=1; jif<jN; jif++)
            {
                for(PetscInt iif=1; iif<iN; iif++)
                {
                    // compute distance
                    Cmpnts dist    = nSub(cent[lzs][j][i], inflowCellCenters[jif][iif]);

                    // set x component to zero
                    dist.x         = 0.0;

                    PetscReal distMag = nMag(dist);

                    if(distMag < minDistMag)
                    {
                        minDistMag = distMag;
                        closestCell.i = iif;
                        closestCell.j = jif;
                        closestCell.k = 0;
                    }
                }
            }

            PetscInt    j_n, j_s, i_w, i_e;
            PetscReal z1, z2, y1, y2;

            if(closestCell.i == 1)
            {
                i_w = 1;
                i_e = 2;
            }
            else if(closestCell.i == iN-1)
            {
                i_w = iN-2;
                i_e = iN-1;
            }
            else
            {
                // find surrounding points for bilinear interpolation
                PetscReal dyl = fabs(cent[lzs][j][i].y - inflowCellCenters[closestCell.j][closestCell.i-1].y);
                PetscReal dyr = fabs(inflowCellCenters[closestCell.j][closestCell.i+1].y - cent[lzs][j][i].y);

                if(dyr<=dyl)
                {
                    i_e = closestCell.i+1;
                    i_w = closestCell.i;
                }
                else if (dyr>dyl)
                {
                    i_e = closestCell.i;
                    i_w = closestCell.i-1;
                }
            }

            y1 = inflowCellCenters[closestCell.j][i_w].y;
            y2 = inflowCellCenters[closestCell.j][i_e].y;

            if(closestCell.j == 1)
            {
                j_s = 1;
                j_n = 2;
            }
            else if(closestCell.j == jN-1)
            {
                j_s = jN-2;
                j_n = jN-1;
            }
            else
            {
                // find surrounding points for bilinear interpolation
                PetscReal dzl = fabs(cent[lzs][j][i].z - inflowCellCenters[closestCell.j-1][closestCell.i].z);
                PetscReal dzr = fabs(inflowCellCenters[closestCell.j+1][closestCell.i].z - cent[lzs][j][i].z);

                if(dzr<=dzl)
                {
                    j_n = closestCell.j+1;
                    j_s = closestCell.j;
                }
                else if (dzr>dzl)
                {
                    j_n = closestCell.j;
                    j_s = closestCell.j-1;
                }
            }

            z1 = inflowCellCenters[j_s][closestCell.i].z;
            z2 = inflowCellCenters[j_n][closestCell.i].z;

            // remove singularity if points are coincident
            PetscReal coeff = (y2 - y1) * (z2 - z1);

            // compute weights
            ifPtr->inflowWeights[j][i][0]  = (y2 - cent[lzs][j][i].y) * (z2 - cent[lzs][j][i].z) / coeff;
            ifPtr->inflowWeights[j][i][1]  = (y2 - cent[lzs][j][i].y) * (cent[lzs][j][i].z - z1) / coeff;
            ifPtr->inflowWeights[j][i][2]  = (cent[lzs][j][i].y - y1) * (z2 - cent[lzs][j][i].z) / coeff;
            ifPtr->inflowWeights[j][i][3]  = (cent[lzs][j][i].y - y1) * (cent[lzs][j][i].z - z1) / coeff;

            /*
            if(ifPtr->inflowWeights[j][i][0]+ifPtr->inflowWeights[j][i][1]+ifPtr->inflowWeights[j][i][2]+ifPtr->inflowWeights[j][i][3] != 1.0)
            {
                PetscPrintf(mesh->MESH_COMM, " --> Warning in function: non-normal weights at P = (0.0, %lf, %lf)\n", cent[lzs][j][i].y, cent[lzs][j][i].z);
                PetscPrintf(mesh->MESH_COMM, "W1 = %lf, W2 = %lf, W3 = %lf, W4 = %lf, SUM = %lf\n", ifPtr->inflowWeights[j][i][0], ifPtr->inflowWeights[j][i][1], ifPtr->inflowWeights[j][i][2], ifPtr->inflowWeights[j][i][3],  ifPtr->inflowWeights[j][i][0]+ifPtr->inflowWeights[j][i][1]+ifPtr->inflowWeights[j][i][2]+ifPtr->inflowWeights[j][i][3]);
            }
            */

            // compute indices with periodicization

            // j south index periodicization
            if(j_s<=ifPtr->n1*ifPtr->prds1)
            {
                ifPtr->closestCells[j][i][0].j = j_s % ifPtr->n1 == 0 ? ifPtr->n1 : j_s % ifPtr->n1;
                ifPtr->closestCells[j][i][2].j = j_s % ifPtr->n1 == 0 ? ifPtr->n1 : j_s % ifPtr->n1;
            }
            else
            {
                ifPtr->closestCells[j][i][0].j = ifPtr->n1;
                ifPtr->closestCells[j][i][2].j = ifPtr->n1;
            }

            // j north index periodicization
            if(j_n<=ifPtr->n1*ifPtr->prds1)
            {
                ifPtr->closestCells[j][i][1].j = j_n % ifPtr->n1 == 0 ? ifPtr->n1 : j_n % ifPtr->n1;
                ifPtr->closestCells[j][i][3].j = j_n % ifPtr->n1 == 0 ? ifPtr->n1 : j_n % ifPtr->n1;
            }
            else
            {
                ifPtr->closestCells[j][i][1].j = ifPtr->n1;
                ifPtr->closestCells[j][i][3].j = ifPtr->n1;
            }

            // i west index periodicization
            if(i_w<=ifPtr->n2*ifPtr->prds2)
            {
                ifPtr->closestCells[j][i][0].i = i_w % ifPtr->n2 == 0 ? ifPtr->n2 : i_w % ifPtr->n2;
                ifPtr->closestCells[j][i][1].i = i_w % ifPtr->n2 == 0 ? ifPtr->n2 : i_w % ifPtr->n2;
            }
            else
            {
                ifPtr->closestCells[j][i][0].i = ifPtr->n2;
                ifPtr->closestCells[j][i][1].i = ifPtr->n2;
            }

            // i east index periodicization
            if(i_e<=ifPtr->n2*ifPtr->prds2)
            {
                ifPtr->closestCells[j][i][2].i = i_e % ifPtr->n2 == 0 ? ifPtr->n2 : i_e % ifPtr->n2;
                ifPtr->closestCells[j][i][3].i = i_e % ifPtr->n2 == 0 ? ifPtr->n2 : i_e % ifPtr->n2;
            }
            else
            {
                ifPtr->closestCells[j][i][2].i = ifPtr->n2;
                ifPtr->closestCells[j][i][3].i = ifPtr->n2;
            }
        }
    }

    for(PetscInt jif=0; jif<jN; jif++)
    {
        std::vector<Cmpnts> ().swap(inflowCellCenters[jif]);
    }

    PetscPrintf(mesh->MESH_COMM, "done\n");

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode mappedInflowInitialize(inletFunctionTypes *ifPtr)
{
    PetscInt           m1 = ifPtr->n1wg, m2 = ifPtr->n2wg;

    PetscInt           i,j,k;

    // allocate memory for temporary inflow data
    ifPtr->ucat_plane = (Cmpnts **)malloc( sizeof(Cmpnts *) * m1 );

    for(j=0; j<m1; j++)
    {
        ifPtr->ucat_plane[j] = (Cmpnts *)malloc( sizeof(Cmpnts) * m2 );
    }

    // Velocity: read the time series in inflowDatabase/U folder and store available times, only
    // master process do the search and reads the files.

    // get file names (physical times at which the slice has been saved)
    DIR *dir; struct dirent *diread;

    std::vector<PetscReal> timeSeries;

    char dataLoc[256];
    sprintf(dataLoc, "inflowDatabase/U/");

    PetscInt     ntimes = 0;

    if ((dir = opendir(dataLoc)) != nullptr)
    {
        while ((diread = readdir(dir)) != nullptr)
        {
            char* timeName = diread->d_name;
            if
            (
                strcmp(timeName, ".") !=0 &&
                strcmp(timeName, "..") !=0
            )
            {
                PetscReal timeValue;
                std::sscanf(timeName, "%lf", &timeValue);
                timeSeries.push_back(timeValue);
                ntimes++;
            }
        }
        closedir (dir);
    }
    else
    {
        char error[512];
        sprintf(error, "could not access ./inflowDatabase/U directory\n");
        fatalErrorInFunction("mappedInflowInitialize", error);
    }

    // sort the timeSeries
    for(PetscInt i=0; i<ntimes; i++)
    {
        PetscReal min   = 1e20;
        PetscReal value = 0.;
        PetscInt label    = 0;

        for(PetscInt s=i; s<ntimes; s++)
        {
            if(timeSeries[s] < min)
            {
                value = timeSeries[s];
                label = s;
                min   = value;
            }
        }
        // exchange values so that elements are not lost
        timeSeries[label] = timeSeries[i];
        // put the min value on the unchanged part at the last index of changed part
        timeSeries[i] = value;
    }

    // store into inflowData struct
    PetscMalloc(sizeof(PetscReal)*ntimes, &(ifPtr->inflowU.inflowTimes));

    for(PetscInt i=0; i<ntimes; i++)
    {
        ifPtr->inflowU.inflowTimes[i]  = timeSeries[i];
    }

    ifPtr->inflowU.nInflowTimes = ntimes;
    ifPtr->inflowU.currentCloseIdx = 0;

    // free memory
    std::vector<PetscReal> ().swap(timeSeries);

    // Temperature: read the time series in inflowDatabase/T folder and store available times, only
    // master process do the search and reads the files.
    if(ifPtr->mapT)
    {
        // allocate memory for temporary inflow data
        ifPtr->t_plane = (PetscReal **)malloc( sizeof(PetscReal *) * m1 );

        for(j=0; j<m1; j++)
        {
            ifPtr->t_plane[j] = (PetscReal *)malloc( sizeof(PetscReal) * m2 );
        }

        // get file names (physical times at which the slice has been saved)
        DIR *dir; struct dirent *diread;

        std::vector<PetscReal> timeSeries;

        char dataLoc[256];
        sprintf(dataLoc, "inflowDatabase/T/");

        PetscInt     ntimes = 0;

        if ((dir = opendir(dataLoc)) != nullptr)
        {
            while ((diread = readdir(dir)) != nullptr)
            {
                char* timeName = diread->d_name;
                if
                (
                    strcmp(timeName, ".") !=0 &&
                    strcmp(timeName, "..") !=0
                )
                {
                    PetscReal timeValue;
                    std::sscanf(timeName, "%lf", &timeValue);
                    timeSeries.push_back(timeValue);
                    ntimes++;
                }
            }
            closedir (dir);
        }
        else
        {
           char error[512];
            sprintf(error, "could not access ./inflowDatabase/T directory\n");
            fatalErrorInFunction("mappedInflowInitialize", error);
        }

        // sort the timeSeries
        for(PetscInt i=0; i<ntimes; i++)
        {
            PetscReal min   = 1e20;
            PetscReal value = 0.;
            PetscInt label    = 0;

            for(PetscInt s=i; s<ntimes; s++)
            {
                if(timeSeries[s] < min)
                {
                    value = timeSeries[s];
                    label = s;
                    min   = value;
                }
            }
            // exchange values so that elements are not lost
            timeSeries[label] = timeSeries[i];
            // put the min value on the unchanged part at the last index of changed part
            timeSeries[i] = value;
        }

        // store into inflowData struct
        PetscMalloc(sizeof(PetscReal)*ntimes, &(ifPtr->inflowT.inflowTimes));

        for(PetscInt i=0; i<ntimes; i++)
        {
            ifPtr->inflowT.inflowTimes[i]  = timeSeries[i];
        }

        ifPtr->inflowT.nInflowTimes = ntimes;
        ifPtr->inflowT.currentCloseIdx = 0;

        // free memory
        std::vector<PetscReal> ().swap(timeSeries);
    }

    // initialize nut inflow sections if the inflow data is mapped
    if(ifPtr->mapNut)
    {
        // allocate memory for temporary inflow data
        ifPtr->nut_plane = (PetscReal **)malloc( sizeof(PetscReal *) * m1 );

        for(j=0; j<m1; j++)
        {
            ifPtr->nut_plane[j] = (PetscReal *)malloc( sizeof(PetscReal) * m2 );
        }

        // Nut: read the time series in inflowDatabase/nut folder and store available times, only
        // master process do the search and reads the files.

        // get file names (physical times at which the slice has been saved)
        DIR *dir; struct dirent *diread;

        std::vector<PetscReal> timeSeries;

        char dataLoc[256];
        sprintf(dataLoc, "inflowDatabase/nut/");

        PetscInt     ntimes = 0;

        if ((dir = opendir(dataLoc)) != nullptr)
        {
            while ((diread = readdir(dir)) != nullptr)
            {
                char* timeName = diread->d_name;
                if
                (
                    strcmp(timeName, ".") !=0 &&
                    strcmp(timeName, "..") !=0
                )
                {
                    PetscReal timeValue;
                    std::sscanf(timeName, "%lf", &timeValue);
                    timeSeries.push_back(timeValue);
                    ntimes++;
                }
            }
            closedir (dir);
        }
        else
        {
           char error[512];
            sprintf(error, "could not access ./inflowDatabase/nut directory\n");
            fatalErrorInFunction("mappedInflowInitialize", error);
        }

        // sort the timeSeries
        for(PetscInt i=0; i<ntimes; i++)
        {
            PetscReal min   = 1e20;
            PetscReal value = 0.;
            PetscInt label    = 0;

            for(PetscInt s=i; s<ntimes; s++)
            {
                if(timeSeries[s] < min)
                {
                    value = timeSeries[s];
                    label = s;
                    min   = value;
                }
            }
            // exchange values so that elements are not lost
            timeSeries[label] = timeSeries[i];
            // put the min value on the unchanged part at the last index of changed part
            timeSeries[i] = value;
        }

        // store into inflowData struct
        PetscMalloc(sizeof(PetscReal)*ntimes, &(ifPtr->inflowNut.inflowTimes));

        for(PetscInt i=0; i<ntimes; i++)
        {
            ifPtr->inflowNut.inflowTimes[i]  = timeSeries[i];
        }

        ifPtr->inflowNut.nInflowTimes = ntimes;
        ifPtr->inflowNut.currentCloseIdx = 0;

        // free memory
        std::vector<PetscReal> ().swap(timeSeries);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode printInflowMappingAction(mesh_ *mesh, inletFunctionTypes *ifPtr)
{
    // check input consistency
    word actionOnJ = "mismatch";
    word actionOnI = "mismatch";

    // action along J
    if((mesh->info.my - 2) - ifPtr->prds1 * ifPtr->n1 == 0)
    {
        if(ifPtr->prds1 == 1)      actionOnJ  = "matchmap";
        else                       actionOnJ  = "periodize";
    }
    if((mesh->info.my - 2) - ifPtr->prds1 * ifPtr->n1 > 0)
    {
        if(ifPtr->prds1 == 0)      actionOnJ  = "mismatch";
        else if(ifPtr->prds1 == 1) actionOnJ  = "extrapolate";
        else                       actionOnJ  = "periodize and extrapolate";
    }

    // action along I
    if((mesh->info.mx - 2) - ifPtr->prds2 * ifPtr->n2 == 0)
    {
        if(ifPtr->prds2 == 1)      actionOnI  = "matchmap";
        else                       actionOnI  = "periodize";
    }
    if((mesh->info.mx - 2) - ifPtr->prds2 * ifPtr->n2 > 0)
    {
        if(ifPtr->prds2 == 0)      actionOnI  = "mismatch";
        else if(ifPtr->prds2 == 1) actionOnI  = "extrapolate";
        else                       actionOnI  = "periodize and extrapolate";
    }

    if(actionOnJ == "mismatch" || actionOnI == "mismatch")
    {
       char error[512];
        sprintf(error, "inflow cannot correctly periodized (j check = %s, i check = %s)", actionOnJ.c_str(), actionOnI.c_str());
        fatalErrorInFunction("printInflowMappingAction",  error);
    }

    PetscPrintf(mesh->MESH_COMM, "   Inflow data mapping info:\n");
    PetscPrintf(mesh->MESH_COMM, "      type 3 on kLeft:\n      - action along j = %s\n      - action along i = %s\n", actionOnJ.c_str(), actionOnI.c_str());

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readInflowU(inletFunctionTypes *ifPtr, clock_ *clock)
{
    PetscInt           m1 = ifPtr->n1wg, m2 = ifPtr->n2wg;

    PetscInt           i,j,k;

    PetscReal     ts, te;

    // initialize velocity temporary local planes in each process
    std::vector<std::vector<Cmpnts>> ucat_plane_tmp_1(m1);
    std::vector<std::vector<Cmpnts>> ucat_plane_tmp_2(m1);

    // set to zero
    for( j=0; j<m1; j++)
    {
        ucat_plane_tmp_1[j].resize(m2);
        ucat_plane_tmp_2[j].resize(m2);

        for(i=0; i<m2; i++)
        {
            ucat_plane_tmp_1[j][i].x = 0;
            ucat_plane_tmp_1[j][i].y = 0;
            ucat_plane_tmp_1[j][i].z = 0;

            ucat_plane_tmp_2[j][i].x = 0;
            ucat_plane_tmp_2[j][i].y = 0;
            ucat_plane_tmp_2[j][i].z = 0;
        }
    }

    // do the search on inflowDatabase/U
    PetscReal nxtItrTime = clock->time;

    PetscInt ntimes   = ifPtr->inflowU.nInflowTimes;

    PetscInt lwrBound = 0;
    PetscInt uprBound = ntimes;

    // if past first iteration do the search on a subset to speed up the process
    if(clock->it > clock->itStart)
    {
        lwrBound = PetscMax(0, (ifPtr->inflowU.currentCloseIdx - 50));
        uprBound = PetscMin(ntimes, (ifPtr->inflowU.currentCloseIdx + 50));
    }

    // get the 2 time values closest to the current time
    PetscInt idx_1 = ifPtr->inflowU.currentCloseIdx;
    PetscInt idx_2 = ifPtr->inflowU.currentCloseIdx + 1;

    PetscReal  diff[ntimes];

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        diff[i] = fabs(ifPtr->inflowU.inflowTimes[i] - nxtItrTime);
    }

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        if(diff[i] < diff[idx_1])
        {
            idx_2 = idx_1;
            idx_1 = i;
        }
        if(diff[i] < diff[idx_2] && i != idx_1)
        {
            idx_2 = i;
        }
    }

    // always put the lower time at idx_1 and higher at idx_2
    if
    (
        ifPtr->inflowU.inflowTimes[idx_2]
        <
        ifPtr->inflowU.inflowTimes[idx_1]
    )
    {
        PetscInt idx_tmp = idx_2;
        idx_2 = idx_1;
        idx_1 = idx_tmp;
    }

    // find interpolation weights
    PetscReal idx = (idx_2 - idx_1) / (ifPtr->inflowU.inflowTimes[idx_2] - ifPtr->inflowU.inflowTimes[idx_1]) * (nxtItrTime - ifPtr->inflowU.inflowTimes[idx_1]) + idx_1;
    PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
    PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

    // get the file name to read
    word fname_1 = "inflowDatabase/U/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_1]);
    word fname_2 = "inflowDatabase/U/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_2]);

    // open the two inflow files and read
    FILE *fp_1 = fopen(fname_1.c_str(), "rb");
    FILE *fp_2 = fopen(fname_2.c_str(), "rb");

    if(!fp_1 || !fp_2)
    {
       char error[512];
        sprintf(error, "cannot open files:\n    %s\n    %s\n", fname_1.c_str(), fname_2.c_str());
        fatalErrorInFunction("read_inflow_data",  error);
    }

    for(j=0; j<m1; j++)
    {
        PetscInt err1, err2;
        err1 = fread(&(ucat_plane_tmp_1[j][0]), sizeof(Cmpnts), m2, fp_1);
        err2 = fread(&(ucat_plane_tmp_2[j][0]), sizeof(Cmpnts), m2, fp_2);
    }

    fclose(fp_1);
    fclose(fp_2);

    // linearly interpolate
    for(j=0; j<m1; j++)
    {
        for(i=0; i<m2; i++)
        {
            ifPtr->ucat_plane[j][i].x
            =
            w1 * ucat_plane_tmp_1[j][i].x + w2 * ucat_plane_tmp_2[j][i].x;

            ifPtr->ucat_plane[j][i].y
            =
            w1 * ucat_plane_tmp_1[j][i].y + w2 * ucat_plane_tmp_2[j][i].y;

            ifPtr->ucat_plane[j][i].z
            =
            w1 * ucat_plane_tmp_1[j][i].z + w2 * ucat_plane_tmp_2[j][i].z;
        }
    }

    // update closest indices
    ifPtr->inflowU.currentCloseIdx   = idx_1;

    // print information
    // PetscPrintf(PETSC_COMM_WORLD, "Selected for reading inflow: Time = %lf, Elapsed Time = %lf s\n", w1 * ifPtr->inflowU.inflowTimes[idx_1] + w2 * ifPtr->inflowU.inflowTimes[idx_2], te-ts);
    // PetscPrintf(PETSC_COMM_WORLD, "                         interpolation weights: w1 = %lf, w2 = %lf\n", w1, w2);
    // PetscPrintf(PETSC_COMM_WORLD, "                         closest avail. times : t1 = %lf, t2 = %lf\n", ifPtr->inflowU.inflowTimes[idx_1], ifPtr->inflowU.inflowTimes[idx_2]);

    // free memory
    for( j=0; j<m1; j++)
    {
        // velocity
        std::vector<Cmpnts> ().swap(ucat_plane_tmp_1[j]);
        std::vector<Cmpnts> ().swap(ucat_plane_tmp_2[j]);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readInflowT(inletFunctionTypes *ifPtr, clock_ *clock)
{
    PetscInt           m1 = ifPtr->n1wg, m2 = ifPtr->n2wg;

    PetscInt           i,j,k;

    PetscReal     ts, te;

    // initialize temperature temporary local planes in each process
    std::vector<std::vector<PetscReal>> t_plane_tmp_1(m1);
    std::vector<std::vector<PetscReal>> t_plane_tmp_2(m1);

    // set to zero
    for( j=0; j<m1; j++)
    {
        t_plane_tmp_1[j].resize(m2);
        t_plane_tmp_2[j].resize(m2);

        for(i=0; i<m2; i++)
        {
            t_plane_tmp_1[j][i] = 0;
            t_plane_tmp_2[j][i] = 0;
        }
    }

    // do the search on inflowDatabase/T
    PetscReal nxtItrTime = clock->time;

    PetscInt ntimes   = ifPtr->inflowT.nInflowTimes;

    PetscInt lwrBound = 0;
    PetscInt uprBound = ntimes;

    // if past first iteration do the search on a subset to speed up the process
    if(clock->it > clock->itStart)
    {
        lwrBound = PetscMax(0, (ifPtr->inflowT.currentCloseIdx - 50));
        uprBound = PetscMin(ntimes, (ifPtr->inflowT.currentCloseIdx + 50));
    }

    // get the 2 time values closest to the current time
    PetscInt idx_1 = ifPtr->inflowT.currentCloseIdx;
    PetscInt idx_2 = ifPtr->inflowT.currentCloseIdx + 1;

    PetscReal  diff[ntimes];

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        diff[i] = fabs(ifPtr->inflowT.inflowTimes[i] - nxtItrTime);
    }

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        if(diff[i] < diff[idx_1])
        {
            idx_2 = idx_1;
            idx_1 = i;
        }
        if(diff[i] < diff[idx_2] && i != idx_1)
        {
            idx_2 = i;
        }
    }

    // always put the lower time at idx_1 and higher at idx_2
    if
    (
        ifPtr->inflowT.inflowTimes[idx_2]
        <
        ifPtr->inflowT.inflowTimes[idx_1]
    )
    {
        PetscInt idx_tmp = idx_2;
        idx_2 = idx_1;
        idx_1 = idx_tmp;
    }

    // find interpolation weights
    PetscReal idx = (idx_2 - idx_1) / (ifPtr->inflowT.inflowTimes[idx_2] - ifPtr->inflowT.inflowTimes[idx_1]) * (nxtItrTime - ifPtr->inflowT.inflowTimes[idx_1]) + idx_1;
    PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
    PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

    // get the file name to read
    word fname_1 = "inflowDatabase/T/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_1]);
    word fname_2 = "inflowDatabase/T/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_2]);

    // open the two inflow files and read
    FILE *fp_1 = fopen(fname_1.c_str(), "rb");
    FILE *fp_2 = fopen(fname_2.c_str(), "rb");

    if(!fp_1 || !fp_2)
    {
        char error[512];
        sprintf(error, "cannot open files:\n    %s\n    %s\n", fname_1.c_str(), fname_2.c_str());
        fatalErrorInFunction("read_inflow_data",  error);
    }

    for(j=0; j<m1; j++)
    {
        PetscInt err1, err2;
        err1 = fread(&(t_plane_tmp_1[j][0]), sizeof(PetscReal), m2, fp_1);
        err2 = fread(&(t_plane_tmp_2[j][0]), sizeof(PetscReal), m2, fp_2);
    }

    fclose(fp_1);
    fclose(fp_2);

    // linearly interpolate
    for(j=0; j<m1; j++)
    {
        for(i=0; i<m2; i++)
        {
            ifPtr->t_plane[j][i]
            =
            w1 * t_plane_tmp_1[j][i] + w2 * t_plane_tmp_2[j][i];
        }
    }

    // update closest indices
    ifPtr->inflowT.currentCloseIdx   = idx_1;

    // free memory
    for( j=0; j<m1; j++)
    {
        // temeperature
        std::vector<PetscReal> ().swap(t_plane_tmp_1[j]);
        std::vector<PetscReal> ().swap(t_plane_tmp_2[j]);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readInflowNut(inletFunctionTypes *ifPtr, clock_ *clock)
{
    PetscInt           m1 = ifPtr->n1wg, m2 = ifPtr->n2wg;

    PetscInt           i,j,k;

    PetscReal     ts, te;

    // initialize nut temporary local planes in each process
    std::vector<std::vector<PetscReal>> nut_plane_tmp_1(m1);
    std::vector<std::vector<PetscReal>> nut_plane_tmp_2(m1);

    // set to zero
    for( j=0; j<m1; j++)
    {
        nut_plane_tmp_1[j].resize(m2);
        nut_plane_tmp_2[j].resize(m2);

        for(i=0; i<m2; i++)
        {
            nut_plane_tmp_1[j][i] = 0;
            nut_plane_tmp_2[j][i] = 0;
        }
    }

    // do the search on inflowDatabase/Nut
    PetscReal nxtItrTime = clock->time;

    PetscInt ntimes   = ifPtr->inflowNut.nInflowTimes;

    PetscInt lwrBound = 0;
    PetscInt uprBound = ntimes;

    // if past first iteration do the search on a subset to speed up the process
    if(clock->it > clock->itStart)
    {
        lwrBound = PetscMax(0, (ifPtr->inflowNut.currentCloseIdx - 50));
        uprBound = PetscMin(ntimes, (ifPtr->inflowNut.currentCloseIdx + 50));
    }

    // get the 2 time values closest to the current time
    PetscInt idx_1 = ifPtr->inflowNut.currentCloseIdx;
    PetscInt idx_2 = ifPtr->inflowNut.currentCloseIdx + 1;

    PetscReal  diff[ntimes];

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        diff[i] = fabs(ifPtr->inflowNut.inflowTimes[i] - nxtItrTime);
    }

    for(PetscInt i=lwrBound; i<uprBound; i++)
    {
        if(diff[i] < diff[idx_1])
        {
            idx_2 = idx_1;
            idx_1 = i;
        }
        if(diff[i] < diff[idx_2] && i != idx_1)
        {
            idx_2 = i;
        }
    }

    // always put the lower time at idx_1 and higher at idx_2
    if
    (
        ifPtr->inflowNut.inflowTimes[idx_2]
        <
        ifPtr->inflowNut.inflowTimes[idx_1]
    )
    {
        PetscInt idx_tmp = idx_2;
        idx_2 = idx_1;
        idx_1 = idx_tmp;
    }

    // find interpolation weights
    PetscReal idx = (idx_2 - idx_1) / (ifPtr->inflowNut.inflowTimes[idx_2] - ifPtr->inflowNut.inflowTimes[idx_1]) * (nxtItrTime - ifPtr->inflowNut.inflowTimes[idx_1]) + idx_1;
    PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
    PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

    // get the file name to read
    word fname_1 = "inflowDatabase/nut/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_1]);
    word fname_2 = "inflowDatabase/nut/" + getArbitraryTimeName(clock, ifPtr->inflowU.inflowTimes[idx_2]);

    // open the two inflow files and read
    FILE *fp_1 = fopen(fname_1.c_str(), "rb");
    FILE *fp_2 = fopen(fname_2.c_str(), "rb");

    if(!fp_1 || !fp_2)
    {
       char error[512];
        sprintf(error, "cannot open files:\n    %s\n    %s\n", fname_1.c_str(), fname_2.c_str());
        fatalErrorInFunction("read_inflow_data",  error);
    }

    for(j=0; j<m1; j++)
    {
        PetscInt err1, err2;
        err1 = fread(&(nut_plane_tmp_1[j][0]), sizeof(PetscReal), m2, fp_1);
        err2 = fread(&(nut_plane_tmp_2[j][0]), sizeof(PetscReal), m2, fp_2);
    }

    fclose(fp_1);
    fclose(fp_2);

    // linearly interpolate
    for(j=0; j<m1; j++)
    {
        for(i=0; i<m2; i++)
        {
            ifPtr->nut_plane[j][i]
            =
            w1 * nut_plane_tmp_1[j][i] + w2 * nut_plane_tmp_2[j][i];
        }
    }

    // update closest indices
    ifPtr->inflowNut.currentCloseIdx = idx_1;

    // free memory
    for( j=0; j<m1; j++)
    {
        // nut
        std::vector<PetscReal> ().swap(nut_plane_tmp_1[j]);
        std::vector<PetscReal> ().swap(nut_plane_tmp_2[j]);
    }

    return(0);
}

//***************************************************************************************************************//

Cmpnts NieuwstadtInflowEvaluate(inletFunctionTypes *ifPtr, PetscReal h)
{
    PetscReal vkConstant = 0.4, C, eta, U = 0.0, V = 0.0;
    complex   alpha, Const, Wg, deltaW, W;

    C   = ifPtr->fc * ifPtr->hInv / (vkConstant * ifPtr->uTau);

    if(h < ifPtr->hInv)
    {
        // evaluate velocity at wanted height
        eta    = h / ifPtr->hInv;
        alpha  = 0.5 + 0.5 * std::sqrt(complex(1.0, 4.0*C));
        Const  = digamma(alpha + 1.0) + digamma(alpha - 1.0) - 2.0 * digamma(1.0);
        Wg     = 1.0 / vkConstant * (std::log(ifPtr->hInv / ifPtr->roughness) - Const);
        deltaW = complex(0.0, -1.0/vkConstant) * alpha * alpha * gamma(alpha) * gamma(alpha) / (C * gamma(2.0*alpha)) * pow(1.0-eta,alpha-1.0)*hypergeom(alpha+1.0, alpha-1.0,2.0*alpha, 1.0-eta);
        W      = Wg - deltaW;

        PetscReal U_p    = W.real() * ifPtr->uTau,
                  V_p    = W.imag() * ifPtr->uTau;

        // evaluate velocity at reference height
        eta    = ifPtr->Href / ifPtr->hInv;
        alpha  = 0.5 + 0.5 * std::sqrt(complex(1.0, 4.0*C));
        Const  = digamma(alpha + 1.0) + digamma(alpha - 1.0) - 2.0 * digamma(1.0);
        Wg     = 1.0 / vkConstant * (std::log(ifPtr->hInv / ifPtr->roughness) - Const);
        deltaW = complex(0.0, -1.0/vkConstant) * alpha * alpha * gamma(alpha) * gamma(alpha) / (C * gamma(2.0*alpha)) * pow(1.0-eta,alpha-1.0)*hypergeom(alpha+1.0, alpha-1.0,2.0*alpha, 1.0-eta);
        W      = Wg - deltaW;

        PetscReal U_r    = W.real() * ifPtr->uTau,
                  V_r    = W.imag() * ifPtr->uTau;

        // normalize
        U_r = U_r / sqrt(U_r*U_r + V_r*V_r);
        V_r = V_r / sqrt(U_r*U_r + V_r*V_r);

        // find rotation angle of the wind profile
        PetscReal theta = std::acos(ifPtr->Udir.x * U_r + ifPtr->Udir.y * V_r);

        // rotate components
        U = std::cos(theta) * U_p - std::sin(theta) * V_p;
        V = std::sin(theta) * U_p + std::cos(theta) * V_p;
    }
    else
    {
        // evaluate velocity above capping height
        eta    = 1.0;
        alpha  = 0.5 + 0.5 * std::sqrt(complex(1.0, 4.0*C));
        Const  = digamma(alpha + 1.0) + digamma(alpha - 1.0) - 2.0 * digamma(1.0);
        Wg     = 1.0 / vkConstant * (std::log(ifPtr->hInv / ifPtr->roughness) - Const);
        deltaW = complex(0.0, -1.0/vkConstant) * alpha * alpha * gamma(alpha) * gamma(alpha) / (C * gamma(2.0*alpha)) * pow(1.0-eta,alpha-1.0)*hypergeom(alpha+1.0, alpha-1.0,2.0*alpha, 1.0-eta);
        W      = Wg - deltaW;

        PetscReal U_p    = W.real() * ifPtr->uTau,
                  V_p    = W.imag() * ifPtr->uTau;

        // evaluate velocity at reference height
        eta    = ifPtr->Href / ifPtr->hInv;
        alpha  = 0.5 + 0.5 * std::sqrt(complex(1.0, 4.0*C));
        Const  = digamma(alpha + 1.0) + digamma(alpha - 1.0) - 2.0 * digamma(1.0);
        Wg     = 1.0 / vkConstant * (std::log(ifPtr->hInv / ifPtr->roughness) - Const);
        deltaW = complex(0.0, -1.0/vkConstant) * alpha * alpha * gamma(alpha) * gamma(alpha) / (C * gamma(2.0*alpha)) * pow(1.0-eta,alpha-1.0)*hypergeom(alpha+1.0, alpha-1.0,2.0*alpha, 1.0-eta);
        W      = Wg - deltaW;

        PetscReal U_r    = W.real() * ifPtr->uTau,
                  V_r    = W.imag() * ifPtr->uTau;

        // normalize
        U_r = U_r / sqrt(U_r*U_r + V_r*V_r);
        V_r = V_r / sqrt(U_r*U_r + V_r*V_r);

        // find rotation angle of the wind profile
        PetscReal theta = std::acos(ifPtr->Udir.x * U_r + ifPtr->Udir.y * V_r);

        // rotate components
        U = std::cos(theta) * U_p - std::sin(theta) * V_p;
        V = std::sin(theta) * U_p + std::cos(theta) * V_p;
    }

    return(nSetFromComponents(U,V,0.0));
}
