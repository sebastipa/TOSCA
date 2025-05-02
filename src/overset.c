//! \file  overset.c
//! \brief overset function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/ibm.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/initialField.h"
#include "include/overset.h"
#include "include/ibmInput.h"

// Comments on TOSCA's overset method: 
// Overset is a mesh refinement technique where two disconnected meshes are overalpped. A finer mesh is fully contained 
// whithin a coarser mesh. At the ghost nodes of the finer mesh, we interpolate the values from the coarser mesh. For 
// the coarser mesh, we introduce a blanking region where the fields are not solved. This blanking region is treated as an
// IBM body and must be offset by 2/3 cells inward with respect to the finer mesh. The boundary conditions for the coarser mesh at the  
// boundaries of the blanking regions are taken from the finer mesh and applied at the interface cells, as done for the IBM method. 
// In order to preserve overall continuity in the finer domain while being able to prescribe fixedValue-type BCs at all 
// boundaries, we need a way to correct for mass imbalance. In TOSCA, while the values of the cartesian velocity at the ghost nodes 
// are always interpolated from the coarser mesh, the contravariant fluxes (the variable that is actually being solved) are
// treated like a fixedValue BC if the flow is entering the domain, while they are solved, similarly to the zeroGradient BC, if the flow is 
// leaving the domain. Mass imbalance is corrected on a per-cell basis, only at those cells where the flow is outgoing. 
// This ensures that the solution remains "attached" to that of the outer domain through the formation of the RHS (which requires the 
// cartesian velocity), but also that the solver is somewhat free to develop its own outflow contravariant fluxes when the flow
// is leaving the domain. This is particularly important in those cases where the mean flow switches direction over the course of the 
// simulation. In fact, when the switch in flow direction happens, the BC has to switch from zeroGradient to fixedValue. If a pure zero 
// gradient would be applied for the outgoing flow, the flow close to the bopundary would have no knowledge of the flow in the outer domain 
// and in some cases they could be largely different when the switch is performed, causing instability. With our method instead, the solution 
// is maintained close to the coarse solution by applying the interpolated cartesian velocity at the ghost nodes all while the flow has been exiting. 
// This makes the transition from zeroGradient to fixedValue less hard on the solver and more stable. 

//***************************************************************************************************************//
// Acceptor cells refer to the cells where the fields are interpolated from a donor or set of donor cells. There are two kinds of acceptor cells possible for an overset domain.
// Domain boundary cells which are at the boundary of the overset mesh. Additionally, if the overset mesh has other overset meshes within it, they a hole cut region is 
// required, which allows to interpolate from the second level mesh to the current mesh. This hole cut region interface cells are the other kind of 
// acceptor cells. The overset parent child system works like a family tree. Starting with the highest level, we recursively move through each lower level of the 
//line until we reach the last generation of a line, before moving to the next branch at the previous level. 

// Note: this overset implementation does not handle intersecting domains at the same level. Overset meshes need to be either cascaded (within another domain) or non intersecting if at same level
//       Use a buffer of alteast 2 cells of the coarser mesh between the donor and acceptor meshes. 

//! \brief Initialize overset coupling 
PetscErrorCode InitializeOverset(domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;

    // Read all hole objects from the input file
    std::vector<HoleObject> holeObjects;
    readHoleObjects(holeObjects, domain[0].info.nHoleRegions);

    // Set initial fields for all domains
    for (PetscInt d = 0; d < nDomains; d++)
    {
        // Only process domains that are not already processed as children
        PetscBool isChild = PETSC_FALSE;
        for (PetscInt other = 0; other < nDomains; other++)
        {
            if (other != d && domain[other].os != NULL)
            {
                for (PetscInt ci = 0; ci < domain[other].os->childMeshId.size(); ci++)
                {
                    if (domain[other].os->childMeshId[ci] == d)
                    {
                        isChild = PETSC_TRUE;
                        break;
                    }
                }
            }
            if (isChild) break;
        }

        if (!isChild)
        {
            SetInitialFieldsOverset(d, domain);
        }
    }

    // find acceptor cells starting from the top level
    for (PetscInt d = 0; d < nDomains; d++)
    {
        // Only process domains that are not already processed as children
        PetscBool isChild = PETSC_FALSE;
        for (PetscInt other = 0; other < nDomains; other++)
        {
            if (other != d && domain[other].os != NULL)
            {
                for (PetscInt ci = 0; ci < domain[other].os->childMeshId.size(); ci++)
                {
                    if (domain[other].os->childMeshId[ci] == d)
                    {
                        isChild = PETSC_TRUE;
                        break;
                    }
                }
            }
            if (isChild) break;
        }

        if (!isChild)
        {
            PetscPrintf(domain[d].mesh->MESH_COMM, "\nStarted recursive acceptor search from domain %ld:\n", d);
            findAcceptorCells(d, domain, 0, holeObjects);
        }
    }

    //it is important to compute all the acceptor cells at a level - from Domain boundary cells, multiple hole cut regions
    //before moving to find the closest donors for the acceptor cells. 

    for (PetscInt d = 0; d < nDomains; d++)
    {
        // Only process domains that are not already processed as children
        PetscBool isChild = PETSC_FALSE;
        for (PetscInt other = 0; other < nDomains; other++)
        {
            if (other != d && domain[other].os != NULL)
            {
                for (PetscInt ci = 0; ci < domain[other].os->childMeshId.size(); ci++)
                {
                    if (domain[other].os->childMeshId[ci] == d)
                    {
                        isChild = PETSC_TRUE;
                        break;
                    }
                }
            }
            if (isChild) break;
        }

        if (!isChild)
        {
            PetscPrintf(domain[d].mesh->MESH_COMM, "\nStarted recursive donor search from domain %ld:\n", d);
            findClosestDomainDonors(d, domain, 0, holeObjects);
        }
    }
    
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode UpdateOversetInterpolation(domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;
    PetscReal ts, te;
    PetscTime(&ts);

    // Main loop: Update interpolation for top-level domains
    for (PetscInt d = 0; d < nDomains; d++)
    {
        overset_ *os   = domain[d].os;

        // Only process domains that are not already processed as children
        PetscBool isChild = PETSC_FALSE;
        for (PetscInt other = 0; other < nDomains; other++)
        {
            if (other != d && domain[other].os != NULL)
            {
                for (PetscInt ci = 0; ci < domain[other].os->childMeshId.size(); ci++)
                {
                    if (domain[other].os->childMeshId[ci] == d)
                    {
                        isChild = PETSC_TRUE;
                        break;
                    }
                }
            }
            if (isChild) break;
        }

        if (!isChild)
        {
            UpdateDomainInterpolation(d, domain, 0);
        }

        // Update acceptor coordinates for dynamic overset
        if (os->dynamicOverset)
        {
            updateAcceptorCoordinates(os);
        }
    }

    PetscTime(&te);
    PetscPrintf(PETSC_COMM_WORLD, "Overset Interpolation Elapsed Time = %lf\n", te - ts);
    return 0;
}

//***************************************************************************************************************//
//! \brief Update overset interpolation for all domains.

// function to update interpolation for a single domain and its dependencies recursively
PetscErrorCode UpdateDomainInterpolation(PetscInt d, domain_ *domain, PetscInt level)
{
    PetscInt nDomains = domain[0].info.nDomains;
    if (d < 0 || d >= nDomains) return 0; 
    if (domain[d].os == NULL) return 0;   

    overset_ *os    = domain[d].os;
    mesh_    *mesh  = domain[d].mesh;

    // Branch 1: Update interpolation from parent meshes to this domain
    for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
    {
        if (os->parentMeshId[pi] != -1)
        {
            mesh_ *parentMesh = domain[os->parentMeshId[pi]].mesh;

            interpolateACellTrilinearP2C(parentMesh, mesh);

            // ABL source term handling
            if (os->access->flags->isAblActive)
            {
                abl_ *ablP = parentMesh->access->abl;
                if (ablP->controllerActive && ablP->controllerAction == "read" && 
                    ablP->controllerType == "timeSeriesFromPrecursor")
                {
                    abl_ *ablC = mesh->access->abl;
                    ablC->preCompSources[0][0] = ablP->preCompSources[0][0];
                    ablC->preCompSources[0][1] = ablP->preCompSources[0][1];
                    ablC->preCompSources[0][2] = ablP->preCompSources[0][2];
                    ablC->preCompSources[0][3] = ablP->preCompSources[0][3];
                }
            }
            
            MPI_Barrier(mesh->MESH_COMM);
        }
    }

    // Update boundary conditions
    UpdateCartesianBCs(domain[d].ueqn);
    UpdateContravariantBCs(domain[d].ueqn);

    // Branch 2: Update interpolation from child meshes to this domain
    for (PetscInt ci = 0; ci < os->childMeshId.size(); ci++)
    {
        if (os->childMeshId[ci] != -1)
        {
            mesh_ *childMesh = domain[os->childMeshId[ci]].mesh;
            PetscInt childId = os->childMeshId[ci];

            interpolateACellTrilinearC2P(childMesh, mesh, childId);

            MPI_Barrier(mesh->MESH_COMM);

            // Update boundary conditions
            UpdateCartesianBCs(domain[d].ueqn);
            UpdateContravariantBCs(domain[d].ueqn);

            // Recursively update child domain
            UpdateDomainInterpolation(childId, domain, level + 1);
        }
    }
    return 0;
}

//***************************************************************************************************************//

// set initial field for a domain and its children recursively
PetscErrorCode SetInitialFieldsOverset(PetscInt d, domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;
    if (d < 0 || d >= nDomains) return 0;
    if (domain[d].os == NULL)   return 0;

    overset_ *os    = domain[d].os;
    mesh_    *mesh  = domain[d].mesh;
    
    SetInitialField(&domain[d]);

    for (PetscInt ci = 0; ci < os->childMeshId.size(); ci++)
    {
        if (os->childMeshId[ci] != -1)
        {
            mesh_ *childMesh = domain[os->childMeshId[ci]].mesh;
            PetscInt childId = os->childMeshId[ci];

            // Recursively initialize child domain
            SetInitialFieldsOverset(childId, domain);
        }
    }
    return 0;
}

//***************************************************************************************************************//

// function to find acceptors cells for a domain and its children recursively
PetscErrorCode findAcceptorCells(PetscInt d, domain_ *domain, PetscInt level, 
                            const std::vector<HoleObject> &holeObjects)
{
    PetscInt nDomains = domain[0].info.nDomains;
    if (d < 0 || d >= nDomains) return 0;
    if (domain[d].os == NULL) return 0;

    flags_ flags = domain[d].flags;
    overset_ *os = domain[d].os;
    mesh_ *mesh = domain[d].mesh;

    // timers 
    PetscReal timeStart, timeEnd;

    //Read overset properties for the current domain
    readOversetProperties(os);

    // Branch 1: create acceptor cells for domain boundary of the overset mesh
    for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
    {
        if (os->parentMeshId[pi] != -1)
        {
            mesh_ *parentMesh = domain[os->parentMeshId[pi]].mesh;

            PetscPrintf(mesh->MESH_COMM, "Creating ghost acceptor cells from %s to %s (level %ld)\n", parentMesh->meshName.c_str(), mesh->meshName.c_str(), level);
        
            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeStart);

            createAcceptorCellOverset(os);

            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeEnd);

            PetscPrintf(mesh->MESH_COMM, "     Elapsed time = %lf\n", timeEnd - timeStart); 
        }
    }

    // Branch 2: create acceptor cells for hole cut boundary
    for (PetscInt ci = 0; ci < os->childMeshId.size(); ci++)
    {
        if (os->childMeshId[ci] != -1)
        {
            mesh_ *childMesh = domain[os->childMeshId[ci]].mesh;
            PetscInt childId = os->childMeshId[ci];

            PetscPrintf(mesh->MESH_COMM, "Creating hole cutting acceptor cells from %s to %s (level %ld)\n", childMesh->meshName.c_str(), mesh->meshName.c_str(), level);
            
            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeStart);

            // Find the hole object for this parent-child pair
            char *holeObjectName = NULL;
            FindHoleObject(holeObjects, d, childId, &holeObjectName);

            if (holeObjectName != NULL)
            {
                PetscPrintf(mesh->MESH_COMM, "     Reading hole object: %s\n", holeObjectName);
                readBlankingIBMObject(os, &domain[d], holeObjectName, holeObjects);   
            }

            createAcceptorCellBackground(os, childId);

            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeEnd);

            PetscPrintf(mesh->MESH_COMM, "     Elapsed time = %lf\n", timeEnd - timeStart); 

            // Recursively initialize child domain
            findAcceptorCells(childId, domain, level + 1, holeObjects);
        }
    }
    return 0;
}

//***************************************************************************************************************//

// function to compute the closest domain for a domain and its children recursively
PetscErrorCode findClosestDomainDonors(PetscInt d, domain_ *domain, PetscInt level, 
                              const std::vector<HoleObject> &holeObjects)
{
    PetscInt nDomains = domain[0].info.nDomains;
    if (d < 0 || d >= nDomains) return 0; 
    if (domain[d].os == NULL) return 0;   

    flags_   flags = domain[d].flags;
    overset_ *os   = domain[d].os;
    mesh_    *mesh = domain[d].mesh;

    // timers 
    PetscReal timeStart, timeEnd;

    // Branch 1: find closest donor cells for domain boundary acceptor cells
    for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
    {
        if (os->parentMeshId[pi] != -1)
        {
            mesh_ *parentMesh = domain[os->parentMeshId[pi]].mesh;

            PetscPrintf(mesh->MESH_COMM, "Creating ghost donor cells from %s to %s (level %ld):\n", parentMesh->meshName.c_str(), mesh->meshName.c_str(), level);
            
            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeStart);

            PetscPrintf(mesh->MESH_COMM, "     Finding closest donor from parent to child...\n");
            findClosestDonorP2C(parentMesh, mesh);
            
            PetscPrintf(mesh->MESH_COMM, "     Interpolating fields...\n");
            interpolateACellTrilinearP2C(parentMesh, mesh);
            
            // sync processors 
            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeEnd);

            PetscPrintf(mesh->MESH_COMM, "     Elapsed time = %lf\n", timeEnd - timeStart); 
        }
    }

    // Branch 2: find closest donor cells for hole cut boundary acceptor cells
    for (PetscInt ci = 0; ci < os->childMeshId.size(); ci++)
    {
        if (os->childMeshId[ci] != -1)
        {
            mesh_ *childMesh = domain[os->childMeshId[ci]].mesh;
            PetscInt childId = os->childMeshId[ci];

            PetscPrintf(mesh->MESH_COMM, "Creating hole cutting donor cells from %s to %s (level %ld):\n", childMesh->meshName.c_str(), mesh->meshName.c_str(), level);

            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeStart);

            PetscPrintf(mesh->MESH_COMM, "     Finding closest donor from child to parent...\n");
            findClosestDonorC2P(childMesh, mesh, childId);
            
            PetscPrintf(mesh->MESH_COMM, "     Interpolating fields...\n");
            interpolateACellTrilinearC2P(childMesh, mesh, childId);

            // sync processors
            MPI_Barrier(mesh->MESH_COMM);
            PetscTime(&timeEnd);

            PetscPrintf(mesh->MESH_COMM, "     Elapsed time = %lf\n", timeEnd - timeStart); 

            // recursively initialize child domain
            findClosestDomainDonors(childId, domain, level + 1, holeObjects);
        }
    }

    return 0;
}
//***************************************************************************************************************//

// function to find the hole object for a parent-child pair
PetscErrorCode FindHoleObject(const std::vector<HoleObject> &holeObjects, 
                                PetscInt parentId, PetscInt childId, char **holeObjectName)
{
    for (const auto &hole : holeObjects)
    {
        if (hole.ownerMesh == parentId && hole.donorMesh == childId)
        {
            *holeObjectName = (char *)hole.bodyName.c_str();
            return 0; // Found the matching hole object
        }
    }

    // If no match found, set to NULL and warn
    *holeObjectName = NULL;
    
    PetscPrintf(PETSC_COMM_WORLD, "Warning: No hole object found for parent %d and child %d\n", parentId, childId);
    
    return 1;
}

//***************************************************************************************************************//

// Function to read hole objects from the input file
PetscErrorCode readHoleObjects(std::vector<HoleObject> &holeObjects, PetscInt numHoleObjects)
{
    PetscErrorCode ierr;
    char objectName[256];
    
    // Clear any existing hole objects
    holeObjects.clear();

    if (numHoleObjects < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "Error: Invalid numHoleObjects (%d) in oversetInput.dat\n", numHoleObjects);
        return 1; // Non-zero error code
    }
    if (numHoleObjects == 0) {
        PetscPrintf(PETSC_COMM_WORLD, "Warning: No hole objects specified in oversetInput.dat\n");
    }

    // Loop through the specified number of hole objects
    for (PetscInt holeIndex = 0; holeIndex < numHoleObjects; holeIndex++)
    {
        // Construct the sub-dictionary name (e.g., "holeObject0")
        sprintf(objectName, "holeObject%ld", holeIndex);

        // Create a new HoleObject
        HoleObject hole;

        // Read bodyName
        word bodyName;
        ierr = readSubDictWord("overset/oversetInput.dat", objectName, "bodyName", &bodyName);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read bodyName for %s\n", objectName);
            return ierr;
        }
        hole.bodyName = bodyName;

        // Read ownerMesh
        PetscReal ownerMeshReal;
        ierr = readSubDictDouble("overset/oversetInput.dat", objectName, "ownerMesh", &ownerMeshReal);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read ownerMesh for %s\n", objectName);
            return ierr;
        }
        hole.ownerMesh = static_cast<PetscInt>(ownerMeshReal);

        // Read donorMesh
        PetscReal donorMeshReal;
        ierr = readSubDictDouble("overset/oversetInput.dat", objectName, "donorMesh", &donorMeshReal);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read donorMesh for %s\n", objectName);
            return ierr;
        }
        hole.donorMesh = static_cast<PetscInt>(donorMeshReal);

        // Read fileType
        word fileType;
        ierr = readSubDictWord("overset/oversetInput.dat", objectName, "fileType", &fileType);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read fileType for %s\n", objectName);
            return ierr;
        }
        hole.fileType = fileType;

        // Read baseLocation
        Cmpnts baseLocation;
        ierr = readSubDictVector("overset/oversetInput.dat", objectName, "baseLocation", &baseLocation);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read baseLocation for %s\n", objectName);
            return ierr;
        }
        hole.baseLocation = baseLocation;

        // Read searchCellRatio
        PetscReal searchCellRatioReal;
        ierr = readSubDictDouble("overset/oversetInput.dat", objectName, "searchCellRatio", &searchCellRatioReal);
        if (ierr != 0) 
        {
            PetscPrintf(PETSC_COMM_WORLD, "Error: Failed to read searchCellRatio for %s\n", objectName);
            return ierr;
        }
        hole.searchCellRatio = static_cast<PetscInt>(searchCellRatioReal);

        // Add the hole object to the vector
        holeObjects.push_back(hole);
    }

    PetscPrintf(PETSC_COMM_WORLD, "\nRead %d hole objects from oversetInput.dat\n\n", (PetscInt)holeObjects.size());
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode readBlankingIBMObject(overset_ *os, domain_ *domain, char *holeObjectName, const std::vector<HoleObject> &holeObjects)
{
    PetscPrintf(PETSC_COMM_WORLD, "     Hole cutting mesh: %s\n", domain->mesh->meshName.c_str());

    //read blank region for background mesh 
    os->oibm = new ibm_;

    //set flags for the ibm obsect 
    os->oibm->dbg           = 0;
    os->oibm->dynamic       = 0;
    os->oibm->computeForce  = 0;
    os->oibm->checkNormal   = 1;
    os->oibm->averageNormal = 0;
    os->oibm->wallShearOn   = 0;
    os->oibm->ibmABL        = 0;
    os->oibm->writeSTL      = 0;
    os->oibm->numBodies     = 1;    

    //set access pointer 
    os->oibm->access = &(domain->access);

    os->oibm->ibmBody = new ibmObject*[1];
    os->oibm->sBox    = new searchBox[1];

    os->oibm->ibmBody[0] = new ibmObject;
    ibmObject   *ibmBody  = os->oibm->ibmBody[0];

    // set pointers to null
    ibmBody->bound          = NULL;
    ibmBody->searchCellList = NULL;
    ibmBody->ibMsh          = NULL;
    ibmBody->ibmRot         = NULL;

    // allocate memory for the IBM mesh of the object
    ibmBody->ibMsh = new ibmMesh;

    // allocate memory for the bounding box  of the object
    ibmBody->bound = new boundingBox;

    //allocate memory for the local ibm elements box
    ibmBody->eBox = new elementBox;

    for (const auto &hole : holeObjects)
    {
        if(hole.bodyName.c_str() == holeObjectName)
        {
            ibmBody->bodyName = hole.bodyName;
            ibmBody->baseLocation = nSet(hole.baseLocation);
            ibmBody->searchCellRatio = hole.searchCellRatio;
            ibmBody->fileType = hole.fileType;
        }
    }

    readIBMBodyFileUCD(ibmBody);

    ibmMesh       *ibMesh  = ibmBody->ibMsh;

    // allocate memory for the element normal, area and center coordinate
    PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eN));
    PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eT1));
    PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eT2));

    PetscMalloc(ibMesh->elems * sizeof(PetscReal), &(ibMesh->eA));
    PetscMalloc(ibMesh->elems * sizeof(Cmpnts), &(ibMesh->eCent));

    //find the ibm cartesian bounding box
    findBodyBoundingBox(os->oibm);

    //find the search cell dimensions from the average cell size
    findSearchCellDim(os->oibm);

    //create the ibm search cell list - ibm elements in each search cell
    createSearchCellList(os->oibm);

    //compute element normals and check that they point outwards
    computeOversetIBMElementNormal(os->oibm);

    oversetIbmSearch(os->oibm);

    return (0);
}
//***************************************************************************************************************//
// overset simulation properties
PetscErrorCode readOversetProperties(overset_ *os)
{
    // to set dynamic overset on
    readDictInt("overset/oversetInput.dat", "dynamicOverset", &(os->dynamicOverset));

    // read the interpolation type
    new(&(os->interpolationType)) word{};
    readDictWord("overset/oversetInput.dat", "interpolationType", &(os->interpolationType));

    // search radius size = cell size X cell factors
    if( (os->interpolationType == "LS1") || (os->interpolationType == "LS2") || (os->interpolationType == "LS3"))
    {
        readDictDouble("overset/oversetInput.dat", "cellFactor", &(os->cellFactor));
    }

    //allocate memory for the overset motion struc if dynamicOverset is on
    if(os->dynamicOverset)
    {
        PetscMalloc(sizeof(oversetMotion), &(os->oMotion));

        oversetMotion *osetMotion = os->oMotion;

        // read the prescribed motion switch
        readSubDictInt("overset/oversetInput.dat", "oversetMotion", "setMotion", &(osetMotion->setMotion));

        if(osetMotion->setMotion)
        {
        // read the interpolation type
        readSubDictWord("overset/oversetInput.dat", "oversetMotion", "motionType", &(osetMotion->motionType));

        if(osetMotion->motionType == "Translation")
        {
            readSubDictVector("overset/oversetInput.dat", "oversetMotion", "prescribedVel", &(osetMotion->prescribedVel));
        }
        }
        else
        {
            // motion not prescribed. motion will be prescribed from IBM motion
            osetMotion->ibmAttached = 1;
        }

    }

    return 0;
}

//*************************************************************** */
PetscErrorCode interpolateACellTrilinearP2C(mesh_ *meshD, mesh_ *meshA)
{
    overset_         *os     = meshA->access->os;
    ueqn_            *ueqnA  = meshA->access->ueqn;
    ueqn_            *ueqnD  = meshD->access->ueqn;
    teqn_            *teqnA  = meshA->access->teqn;
    teqn_            *teqnD  = meshD->access->teqn;
    flags_           *flags  = meshA->access->flags;
    DM               daA     = meshA->da, fdaA = meshA->fda;
    DMDALocalInfo    infoA   = meshA->info;
    DM               daD     = meshD->da, fdaD = meshD->fda;
    DMDALocalInfo    infoD   = meshD->info;

    PetscInt         xs = infoA.xs, xe = infoA.xs + infoA.xm;
    PetscInt         ys = infoA.ys, ye = infoA.ys + infoA.ym;
    PetscInt         zs = infoA.zs, ze = infoA.zs + infoA.zm;
    PetscInt         mx = infoA.mx, my = infoA.my, mz = infoA.mz;

    PetscInt         i, j, k, ic, kc, jc, b, n, m;

    Cmpnts           ***lucatD, ***ucatA, ***cent, ucart;
    PetscReal        ***ltempD, ***tempA, Temp;

    Cmpnts           pCoor;

    PetscMPIInt      rankA, sizeA, rankD, sizeD;
    PetscInt         sum_ind1 = 0;

    MPI_Comm_size(meshA->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshA->MESH_COMM, &rankA);

    MPI_Comm_size(meshD->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshD->MESH_COMM, &rankD);

    std::vector<Acell> aCell = os->aCellDb;
    std::vector<Dcell> dCell = os->closestDonorDb;
    std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMatDb;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcDb;

    DMDAVecGetArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecGetArray(fdaA, ueqnA->Ucat, &ucatA);
    DMDAVecGetArray(fdaD, meshD->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecGetArray(daA, teqnA->Tmprt, &tempA);
    }

    // loop through the ranks
    for(n = 0; n < sizeA; n++)
    {

        if(NumAcellPerProc[n]!=0)
        {

            if(AcellProcMat[n][rankD] !=MPI_UNDEFINED)
            {

                // loop through the aCell cells of a given processor n
                for(b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++)
                {
                    // aCell cell index
                    i = aCell[b].indi;
                    j = aCell[b].indj;
                    k = aCell[b].indk;

                    pCoor.x = aCell[b].coorx;
                    pCoor.y = aCell[b].coory;
                    pCoor.z = aCell[b].coorz;

                    ucart.x = 0.0;
                    ucart.y = 0.0;
                    ucart.z = 0.0;
                    Temp = 0.0;

                    if (rankD == dCell[b].rank)
                    {
                        ic = dCell[b].indi;
                        jc = dCell[b].indj;
                        kc = dCell[b].indk;

                        vectorPointLocalVolumeInterpolation
                        (
                                meshD,
                                pCoor.x, pCoor.y, pCoor.z,
                                ic, jc, kc,
                                cent,
                                lucatD,
                                ucart
                        );

                        if (flags->isTeqnActive)
                        {
                            scalarPointLocalVolumeInterpolation
                            (
                                    meshD,
                                    pCoor.x, pCoor.y, pCoor.z,
                                    ic, jc, kc,
                                    cent,
                                    ltempD,
                                    Temp
                            );
                        }

                        MPI_Send(&ucart, 3, MPIU_REAL, aCell[b].rank, 0, meshD->MESH_COMM);
                        MPI_Send(&Temp, 1, MPIU_REAL, aCell[b].rank, 1, meshD->MESH_COMM);

                        // if(k == 25 && j == 25 && i == 20)
                        //     PetscPrintf(PETSC_COMM_SELF, "donor = %ld %ld %ld, ucatD = %lf %lf %lf, ucatA = %lf %lf %lf\n", kc, jc, ic, lucatD[kc][jc][ic].x, lucatD[kc][jc][ic].y, lucatD[kc][jc][ic].z, ucart.x, ucart.y, ucart.z );
                    }

                    if (rankA == aCell[b].rank)
                    {

                        MPI_Recv(&ucart, 3, MPIU_REAL, dCell[b].rank, 0, meshD->MESH_COMM, MPI_STATUS_IGNORE);
                        MPI_Recv(&Temp, 1, MPIU_REAL, dCell[b].rank, 1, meshD->MESH_COMM, MPI_STATUS_IGNORE);

                        ucatA[k][j][i].x = ucart.x;
                        ucatA[k][j][i].y = ucart.y;
                        ucatA[k][j][i].z = ucart.z;

                        if (flags->isTeqnActive)
                        {
                            tempA[k][j][i] = Temp;
                        }

                    }
                }
            }

            sum_ind1 +=NumAcellPerProc[n];
        }

    }

    std::vector<Acell> ().swap(aCell);
    std::vector<Dcell> ().swap(dCell);
    std::vector<std::vector<PetscInt>> ().swap(AcellProcMat);
    std::vector<PetscInt> ().swap(NumAcellPerProc);

    DMDAVecRestoreArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecRestoreArray(fdaA, ueqnA->Ucat, &ucatA);
    DMDAVecRestoreArray(fdaD, meshD->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecRestoreArray(daA, teqnA->Tmprt, &tempA);

        DMGlobalToLocalBegin(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
        DMGlobalToLocalEnd(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
    }

    DMGlobalToLocalBegin(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);
    DMGlobalToLocalEnd(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);

    setBackgroundBC(meshA);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode interpolateACellTrilinearC2P(mesh_ *meshD, mesh_ *meshA, PetscInt donorId)
{
    overset_         *os     = meshA->access->os;
    ueqn_            *ueqnA  = meshA->access->ueqn;
    ueqn_            *ueqnD  = meshD->access->ueqn;
    teqn_            *teqnA  = meshA->access->teqn;
    teqn_            *teqnD  = meshD->access->teqn;
    flags_           *flags  = meshA->access->flags;
    DM               daA     = meshA->da, fdaA = meshA->fda;
    DMDALocalInfo    infoA   = meshA->info;
    DM               daD     = meshD->da, fdaD = meshD->fda;
    DMDALocalInfo    infoD   = meshD->info;

    PetscInt         xs = infoA.xs, xe = infoA.xs + infoA.xm;
    PetscInt         ys = infoA.ys, ye = infoA.ys + infoA.ym;
    PetscInt         zs = infoA.zs, ze = infoA.zs + infoA.zm;
    PetscInt         mx = infoA.mx, my = infoA.my, mz = infoA.mz;

    PetscInt         i, j, k, b, n;

    Cmpnts           ***lucatD, ***ucatA, ***cent;
    PetscReal        ***ltempD, ***tempA;

    Cmpnts           pCoor, ucart;
    PetscReal        Temp;

    PetscMPIInt      rankA, sizeA, rankD, sizeD;
    PetscInt         sum_ind1 = 0;

    MPI_Comm_size(meshA->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshA->MESH_COMM, &rankA);
    MPI_Comm_size(meshD->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshD->MESH_COMM, &rankD);

    std::vector<Acell> aCell = os->aCellHc;
    std::vector<Dcell> dCell = os->closestDonorHc;
    std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMatHc;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcHc;

    DMDAVecGetArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecGetArray(fdaA, ueqnA->Ucat, &ucatA);
    DMDAVecGetArray(fdaD, meshD->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecGetArray(daA, teqnA->Tmprt, &tempA);
    }

    // Map to store interpolated values for each parent cell
    std::map<PetscInt, std::vector<Cmpnts>> vertexVelocities; // parentCellId -> list of velocities
    std::map<PetscInt, std::vector<PetscReal>> vertexTemps;   // parentCellId -> list of temperatures
    std::map<PetscInt, std::tuple<PetscInt, PetscInt, PetscInt>> cellIndices; // parentCellId -> (i, j, k)

    // Loop through the ranks
    for (n = 0; n < sizeA; n++)
    {
        if (NumAcellPerProc[n] != 0)
        {
            if (AcellProcMat[n][rankD] != MPI_UNDEFINED)
            {
                // Loop through the aCell entries (vertices) of a given processor n
                for (b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++)
                {
                    if (aCell[b].donorId == donorId)
                    {
                        pCoor.x = aCell[b].coorx;
                        pCoor.y = aCell[b].coory;
                        pCoor.z = aCell[b].coorz;

                        ucart.x = 0.0;
                        ucart.y = 0.0;
                        ucart.z = 0.0;
                        Temp = 0.0;

                        if (rankD == dCell[b].rank)
                        {
                            // Interpolate velocity at the vertex
                            PetscInt ic = dCell[b].indi;
                            PetscInt jc = dCell[b].indj;
                            PetscInt kc = dCell[b].indk;

                            vectorPointLocalVolumeInterpolation
                            (
                                meshD,
                                pCoor.x, pCoor.y, pCoor.z,
                                ic, jc, kc,
                                cent,
                                lucatD,
                                ucart
                            );

                            if (flags->isTeqnActive)
                            {
                                scalarPointLocalVolumeInterpolation
                                (
                                    meshD,
                                    pCoor.x, pCoor.y, pCoor.z,
                                    ic, jc, kc,
                                    cent,
                                    ltempD,
                                    Temp
                                );
                            }

                            // Send interpolated values to acceptor processor
                            MPI_Send(&ucart, 3, MPIU_REAL, aCell[b].rank, b, meshD->MESH_COMM);
                            if (flags->isTeqnActive)
                            {
                                MPI_Send(&Temp, 1, MPIU_REAL, aCell[b].rank, b + sizeA, meshD->MESH_COMM);
                            }
                        }

                        if (rankA == aCell[b].rank)
                        {
                            // Receive interpolated values
                            MPI_Recv(&ucart, 3, MPIU_REAL, dCell[b].rank, b, meshD->MESH_COMM, MPI_STATUS_IGNORE);
                            vertexVelocities[aCell[b].parentCellId].push_back(ucart);

                            if (flags->isTeqnActive)
                            {
                                MPI_Recv(&Temp, 1, MPIU_REAL, dCell[b].rank, b + sizeA, meshD->MESH_COMM, MPI_STATUS_IGNORE);
                                vertexTemps[aCell[b].parentCellId].push_back(Temp);
                            }

                            // Store cell indices for this parentCellId (only once per cell)
                            if (cellIndices.find(aCell[b].parentCellId) == cellIndices.end())
                            {
                                cellIndices[aCell[b].parentCellId] = {aCell[b].indi, aCell[b].indj, aCell[b].indk};
                            }
                        }
                    }
                }
            }

            sum_ind1 += NumAcellPerProc[n];
        }
    }

    //perform averaging of the fields in any processor that has the acceptor cells 
    bool isAcceptorProcessor = aCell.empty();
    for (const auto& cell : aCell)
    {
        if (cell.donorId == donorId && cell.rank == rankA)
        {
            isAcceptorProcessor = true;
            break;
        }
    }

    // Average the velocities and temperatures for each parent cell
    if (isAcceptorProcessor) // Ensure only acceptor processor processes
    {
        for (const auto& [parentCellId, velocities] : vertexVelocities)
        {
            if (velocities.size() == 8) // Ensure all 8 vertices are present
            {
                Cmpnts avgVelocity = {0.0, 0.0, 0.0};
                for (const auto& v : velocities)
                {
                    avgVelocity.x += v.x / 8.0;
                    avgVelocity.y += v.y / 8.0;
                    avgVelocity.z += v.z / 8.0;
                }

                // Get cell indices
                auto [i, j, k] = cellIndices[parentCellId];

                // Store averaged velocity
                ucatA[k][j][i].x = avgVelocity.x;
                ucatA[k][j][i].y = avgVelocity.y;
                ucatA[k][j][i].z = avgVelocity.z;

                // Handle temperature if active
                if (flags->isTeqnActive)
                {
                    PetscReal avgTemp = 0.0;
                    for (const auto& t : vertexTemps[parentCellId])
                    {
                        avgTemp += t / 8.0;
                    }
                    tempA[k][j][i] = avgTemp;
                }
            }
        }
    }

    // Clean up
    std::vector<Acell>().swap(aCell);
    std::vector<Dcell>().swap(dCell);
    std::vector<std::vector<PetscInt>>().swap(AcellProcMat);
    std::vector<PetscInt>().swap(NumAcellPerProc);

    DMDAVecRestoreArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecRestoreArray(fdaA, ueqnA->Ucat, &ucatA);
    DMDAVecRestoreArray(fdaD, meshD->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecRestoreArray(daA, teqnA->Tmprt, &tempA);

        DMGlobalToLocalBegin(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
        DMGlobalToLocalEnd(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
    }

    DMGlobalToLocalBegin(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);
    DMGlobalToLocalEnd(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);

    setBackgroundBC(meshA);

    return 0;
}

//***************************************************************************************************************//
PetscErrorCode setBackgroundBC(mesh_ *meshA)
{
    ueqn_         *ueqn = meshA->access->ueqn;
    teqn_         *teqn = meshA->access->teqn;
    flags_        *flags= meshA->access->flags;
    DM            da    = meshA->da, fda = meshA->fda;
    DMDALocalInfo info  = meshA->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    PetscReal        ***meshTag, ***Temp;
    PetscReal        ucx, ucy, ucz;

    Cmpnts           ***lucat, ***ucat, ***ucont, ***icsi, ***jeta, ***kzet;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(fda, meshA->lICsi, &icsi);
    DMDAVecGetArray(fda, meshA->lJEta, &jeta);
    DMDAVecGetArray(fda, meshA->lKZet, &kzet);
    DMDAVecGetArray(da, meshA->lmeshTag, &meshTag);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->Tmprt, &Temp);
    }

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        if (isZeroedIFace(k, j, i, i+1, meshTag))
        {
            ucont[k][j][i].x = 0;
        }

        if (isZeroedJFace(k, j, i, j+1, meshTag))
        {
            ucont[k][j][i].y = 0;
        }

        if (isZeroedKFace(k, j, i, k+1, meshTag))
        {
            ucont[k][j][i].z = 0;
        }

        if (isInterpolatedIFace(k,j,i, i+1, meshTag)) 
        {
            ucx = (lucat[k][j][i].x + lucat[k][j][i+1].x) * 0.5;
            ucy = (lucat[k][j][i].y + lucat[k][j][i+1].y) * 0.5;
            ucz = (lucat[k][j][i].z + lucat[k][j][i+1].z) * 0.5;

            ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
        }

        if (isInterpolatedJFace(k,j,i, j+1, meshTag)) 
        {
            ucx = (lucat[k][j+1][i].x + lucat[k][j][i].x) * 0.5;
            ucy = (lucat[k][j+1][i].y + lucat[k][j][i].y) * 0.5;
            ucz = (lucat[k][j+1][i].z + lucat[k][j][i].z) * 0.5;

            ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
        }

        if (isInterpolatedKFace(k,j,i, k+1, meshTag))
        {
            ucx = (lucat[k+1][j][i].x + lucat[k][j][i].x) * 0.5;
            ucy = (lucat[k+1][j][i].y + lucat[k][j][i].y) * 0.5;
            ucz = (lucat[k+1][j][i].z + lucat[k][j][i].z) * 0.5;

            ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
        }

        if(isZeroedCell(k, j, i, meshTag))
        {
            mSetValue(ucat[k][j][i], 0);

            if (flags->isTeqnActive)
            {
                PetscReal tRef;

                if(flags->isAblActive) tRef = teqn->access->abl->tRef;
                else                   tRef = teqn->access->constants->tRef;
                
                Temp[k][j][i] = tRef;
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(fda, meshA->lICsi, &icsi);
    DMDAVecRestoreArray(fda, meshA->lJEta, &jeta);
    DMDAVecRestoreArray(fda, meshA->lKZet, &kzet);
    DMDAVecRestoreArray(da,  meshA->lmeshTag, &meshTag);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->Tmprt, &Temp);

        DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    }

    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    return (0);
}

//***************************************************************************************************************//
PetscErrorCode createAcceptorCellOverset(overset_ *os)
{
    mesh_         *mesh = os->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, b;

    PetscReal     ***aj;
    Cmpnts        ***cent, ***csi, ***eta, ***zet;

    PetscMPIInt   rank, size;
    PetscInt      localCount = 0, globalCount = 0;

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::vector<PetscInt> localNum(size, 0);
    std::vector<PetscInt> globalNum(size, 0);
    std::vector<Acell> localCells;

    lxs = xs; lxe = xe; if (xs == 0) lxs = xs + 1; if (xe == mx) lxe = xe - 1;
    lys = ys; lye = ye; if (ys == 0) lys = ys + 1; if (ye == my) lye = ye - 1;
    lzs = zs; lze = ze; if (zs == 0) lzs = zs + 1; if (ze == mz) lze = ze - 1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    for (k = zs; k < ze; k++) 
    {
        for (j = ys; j < ye; j++) 
        {
            for (i = xs; i < xe; i++) 
            {
                if (isOnCornerCellCenters(i, j, k, info)) continue;
                if ((k == 0) || (k == mz - 1) || (j == 0) || (j == my - 1) || (i == 0) || (i == mx - 1)) 
                {
                    localCount++;
                }
            }
        }
    }

    localNum[rank] = localCount;
    MPI_Allreduce(&localNum[0], &globalNum[0], size, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    os->NumAcellPerProcDb = globalNum;

    PetscInt offset = 0;
    for (b = 0; b < rank; b++) 
    {
        offset += globalNum[b];
    }

    globalCount = 0;
    for (b = 0; b < size; b++) 
    {
        globalCount += globalNum[b];
    }

    os->aCellDb.resize(globalCount);
    for (PetscInt idx = 0; idx < globalCount; idx++) 
    {
        os->aCellDb[idx].indi = 0;
        os->aCellDb[idx].indj = 0;
        os->aCellDb[idx].indk = 0;
        os->aCellDb[idx].coorx = 0.0;
        os->aCellDb[idx].coory = 0.0;
        os->aCellDb[idx].coorz = 0.0;
        os->aCellDb[idx].rank = -1;
        os->aCellDb[idx].cell_size = 0.0;
        os->aCellDb[idx].face = 0;
        os->aCellDb[idx].donorId = 0;
        os->aCellDb[idx].parentCellId = -1;
    }

    localCells.resize(localCount);
    PetscInt localIndex = 0;
    for (k = zs; k < ze; k++) 
    {
        for (j = ys; j < ye; j++) 
        {
            for (i = xs; i < xe; i++) 
            {
                if (isOnCornerCellCenters(i, j, k, info)) continue;
                if ((k == 0) || (k == mz - 1) || (j == 0) || (j == my - 1) || (i == 0) || (i == mx - 1)) 
                {
                    Acell &cell = localCells[localIndex];
                    cell.indi = i;
                    cell.indj = j;
                    cell.indk = k;
                    cell.coorx = cent[k][j][i].x;
                    cell.coory = cent[k][j][i].y;
                    cell.coorz = cent[k][j][i].z;
                    cell.rank = rank;
                    cell.cell_size = 0;
                    cell.face = 0;
                    cell.donorId = 0;
                    cell.parentCellId = -1;
                    localIndex++;
                }
            }
        }
    }

    for (PetscInt idx = 0; idx < localCount; idx++) 
    {
        os->aCellDb[offset + idx] = localCells[idx];
    }

    MPI_Datatype mpi_Acell;
    defineStruct_Acell(&mpi_Acell);
    MPI_Op sumstruct;

    MPI_Op_create(sum_struct_Acell, 1, &sumstruct);

    if (globalCount > 0) 
    {
        MPI_Allreduce(MPI_IN_PLACE, &os->aCellDb[0], globalCount, mpi_Acell, sumstruct, mesh->MESH_COMM);
    }

    MPI_Op_free(&sumstruct);
    MPI_Type_free(&mpi_Acell);
    localCells.clear();
    localNum.clear();
    globalNum.clear();

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return 0;
}

//***************************************************************************************************************//
//! \brief Create the list of background acceptor cells 

PetscErrorCode createAcceptorCellBackground(overset_ *os, PetscInt donorMeshId)
{
    mesh_         *mesh = os->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, b;
    
    Vec           Coor;
    PetscReal     ***aj, ***meshTag;
    Cmpnts        ***coor, ***cent, ***nvert; // Added nvert for vertex coordinates

    PetscMPIInt   rank, size;
    PetscInt      localCount = 0, globalCount = 0, prevGlobalCount = os->aCellHc.size();

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::vector<PetscInt> localNum(size, 0);
    std::vector<PetscInt> globalNum(size, 0);
    std::vector<Acell> localNewCells;

    lxs = xs; lxe = xe; if (xs == 0) lxs = xs + 1; if (xe == mx) lxe = xe - 1;
    lys = ys; lye = ye; if (ys == 0) lys = ys + 1; if (ye == my) lye = ye - 1;
    lzs = zs; lze = ze; if (zs == 0) lzs = zs + 1; if (ze == mz) lze = ze - 1;

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // Counter for unique parent cell IDs (local to this rank)
    PetscInt localCellId = 0;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if (isInterpolatedCell(k, j, i, meshTag)) 
                {
                    bool isNewCell = true;
                    for (const auto& cell : os->aCellHc) 
                    {
                        if (cell.indi == i && cell.indj == j && cell.indk == k && cell.rank == rank) 
                        {
                            isNewCell = false;
                            break;
                        }
                    }

                    if (isNewCell) 
                    {
                        // Compute the 8 vertex coordinates for the cell (i, j, k)

                        PetscReal vx[8], vy[8], vz[8];
                        
                        // Vertex indices: (i-1,i), (j,j-1), (k,k-1)

                        vx[0] = coor[k-1][j-1][i-1].x; // (i-1, j-1, k-1)
                        vy[0] = coor[k-1][j-1][i-1].y;
                        vz[0] = coor[k-1][j-1][i-1].z;
                        
                        vx[1] = coor[k-1][j-1][i].x;   // (i, j-1, k-1)
                        vy[1] = coor[k-1][j-1][i].y;
                        vz[1] = coor[k-1][j-1][i].z;
                        
                        vx[2] = coor[k-1][j][i-1].x;   // (i-1, j, k-1)
                        vy[2] = coor[k-1][j][i-1].y;
                        vz[2] = coor[k-1][j][i-1].z;
                        
                        vx[3] = coor[k-1][j][i].x;     // (i, j, k-1)
                        vy[3] = coor[k-1][j][i].y;
                        vz[3] = coor[k-1][j][i].z;
                        
                        vx[4] = coor[k][j-1][i-1].x;   // (i-1, j-1, k)
                        vy[4] = coor[k][j-1][i-1].y;
                        vz[4] = coor[k][j-1][i-1].z;
                        
                        vx[5] = coor[k][j-1][i].x;     // (i, j-1, k)
                        vy[5] = coor[k][j-1][i].y;
                        vz[5] = coor[k][j-1][i].z;
                        
                        vx[6] = coor[k][j][i-1].x;     // (i-1, j, k)
                        vy[6] = coor[k][j][i-1].y;
                        vz[6] = coor[k][j][i-1].z;
                        
                        vx[7] = coor[k][j][i].x;       // (i, j, k)
                        vy[7] = coor[k][j][i].y;
                        vz[7] = coor[k][j][i].z;
                        

                        // Assign a unique parentCellId for this acceptor cell
                        PetscInt parentCellId = rank * 1000000 + localCellId; 
                        localCellId++;

                        // Add 8 vertices as Acell entries
                        for (PetscInt v = 0; v < 8; v++) 
                        {
                            Acell newCell;
                            newCell.indi = i;
                            newCell.indj = j;
                            newCell.indk = k;
                            newCell.coorx = vx[v]; 
                            newCell.coory = vy[v];
                            newCell.coorz = vz[v];
                            newCell.rank = rank;
                            newCell.cell_size = 0;
                            newCell.face = 0;
                            newCell.donorId = donorMeshId;
                            newCell.parentCellId = parentCellId; 
                            localNewCells.push_back(newCell);
                            localCount++;
                        }
                    }
                }
            }
        }
    }

    localNum[rank] = localCount;
    MPI_Allreduce(&localNum[0], &globalNum[0], size, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    if (os->NumAcellPerProcHc.empty()) 
    {
        os->NumAcellPerProcHc = globalNum;
    } 
    else 
    {
        for (PetscInt b = 0; b < size; b++) 
        {
            os->NumAcellPerProcHc[b] += globalNum[b];
        }
    }

    PetscInt offset = prevGlobalCount;
    for (PetscInt b = 0; b < rank; b++) 
    {
        offset += globalNum[b];
    }
    globalCount = 0;

    for (PetscInt b = 0; b < size; b++) 
    {
        globalCount += globalNum[b];
    }

    os->aCellHc.resize(prevGlobalCount + globalCount);
    for (PetscInt idx = prevGlobalCount; idx < prevGlobalCount + globalCount; idx++) 
    {
        os->aCellHc[idx].indi = 0;
        os->aCellHc[idx].indj = 0;
        os->aCellHc[idx].indk = 0;
        os->aCellHc[idx].coorx = 0.0;
        os->aCellHc[idx].coory = 0.0;
        os->aCellHc[idx].coorz = 0.0;
        os->aCellHc[idx].rank = -1;
        os->aCellHc[idx].cell_size = 0.0;
        os->aCellHc[idx].face = 0;
        os->aCellHc[idx].donorId = 0;
        os->aCellHc[idx].parentCellId = -1; // Initialize parentCellId
    }

    for (PetscInt idx = 0; idx < localCount; idx++) 
    {
        os->aCellHc[offset + idx] = localNewCells[idx];
    }

    MPI_Datatype mpi_Acell;
    defineStruct_Acell(&mpi_Acell); // Must include parentCellId in MPI datatype
    MPI_Op sumstruct;
    MPI_Op_create(sum_struct_Acell, 1, &sumstruct);

    if (globalCount > 0) 
    {
        MPI_Allreduce(MPI_IN_PLACE, &os->aCellHc[prevGlobalCount], globalCount, mpi_Acell, sumstruct, mesh->MESH_COMM);
    }

    // Sort os->aCellHc by rank and parentCellId for easier grouping
    if (!os->aCellHc.empty()) 
    {
        std::sort(os->aCellHc.begin(), os->aCellHc.end(), 
                  [](const Acell& a, const Acell& b) { 
                      if (a.rank == b.rank) return a.parentCellId < b.parentCellId;
                      return a.rank < b.rank; 
                  });

        // Update NumAcellPerProcHc
        std::fill(os->NumAcellPerProcHc.begin(), os->NumAcellPerProcHc.end(), 0);

        for (const auto& cell : os->aCellHc) 
        {
            if (cell.rank >= 0 && cell.rank < size) 
            {
                os->NumAcellPerProcHc[cell.rank]++;
            }
        }
    }

    MPI_Op_free(&sumstruct);
    MPI_Type_free(&mpi_Acell);

    localNewCells.clear();
    localNum.clear();
    globalNum.clear();

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, Coor, &coor);

    return 0;
}

// experimental octree structure for donor cells
struct OctreeNode 
{
    std::vector<cellIds> donorCells; 
    Cmpnts minBounds;  // Minimum bounds of the node
    Cmpnts maxBounds;  // Maximum bounds of the node
    OctreeNode* children[8]; // Pointers to child nodes (8 for octree)

    OctreeNode(Cmpnts minB, Cmpnts maxB) : minBounds(minB), maxBounds(maxB) 
    {
        for (int i = 0; i < 8; i++) children[i] = nullptr;
    }

    ~OctreeNode() 
    {
        for (int i = 0; i < 8; i++) 
        {
            if (children[i]) delete children[i];
        }
    }
};

void buildOctree(OctreeNode* node, Cmpnts*** donorCells, 
    PetscInt lxs, PetscInt lxe, PetscInt lys, PetscInt lye, PetscInt lzs, PetscInt lze,
    PetscInt maxDepth, PetscInt maxCellsPerNode) 
{
    // Count the number of donor cells in the current node's bounds
    PetscInt cellCount = 0;
    for (PetscInt k = lzs; k < lze; k++) 
    for (PetscInt j = lys; j < lye; j++) 
    for (PetscInt i = lxs; i < lxe; i++) 
    {
        Cmpnts centroid = donorCells[k][j][i];
        if
        (
            centroid.x >= node->minBounds.x && centroid.x < node->maxBounds.x &&
            centroid.y >= node->minBounds.y && centroid.y < node->maxBounds.y &&
            centroid.z >= node->minBounds.z && centroid.z < node->maxBounds.z
        ) 
        {
            cellCount++;
        }
    }

    // ff the number of cells is below the threshold or max depth is reached, store the cells in this node
    if (cellCount <= maxCellsPerNode || maxDepth == 0) 
    {
        for (PetscInt k = lzs; k < lze; k++) 
        for (PetscInt j = lys; j < lye; j++) 
        for (PetscInt i = lxs; i < lxe; i++) 
        {
            Cmpnts centroid = donorCells[k][j][i];

            if
            (
                centroid.x >= node->minBounds.x && centroid.x < node->maxBounds.x &&
                centroid.y >= node->minBounds.y && centroid.y < node->maxBounds.y &&
                centroid.z >= node->minBounds.z && centroid.z < node->maxBounds.z
            ) 
            {
                // Store the cell in this node
                cellIds cell;
                cell.i = i;
                cell.j = j;
                cell.k = k;
                node->donorCells.push_back(cell);   
            }
        }
        return;
    }

    // Compute midpoints of the current node's bounds
    Cmpnts mid;
    mid.x = (node->minBounds.x + node->maxBounds.x) / 2.0;
    mid.y = (node->minBounds.y + node->maxBounds.y) / 2.0;
    mid.z = (node->minBounds.z + node->maxBounds.z) / 2.0;

    // Create child nodes
    for (int i = 0; i < 8; i++) 
    {
        Cmpnts childMin = node->minBounds;
        Cmpnts childMax = mid;

        if (i & 1) childMin.x = mid.x, childMax.x = node->maxBounds.x;
        if (i & 2) childMin.y = mid.y, childMax.y = node->maxBounds.y;
        if (i & 4) childMin.z = mid.z, childMax.z = node->maxBounds.z;

        node->children[i] = new OctreeNode(childMin, childMax);
    }

    // Recursively build child nodes
    for (int i = 0; i < 8; i++) 
    {
        buildOctree(node->children[i], donorCells, lxs, lxe, lys, lye, lzs, lze, maxDepth - 1, maxCellsPerNode);
    }
}

Dcell searchOctree(OctreeNode* node, PetscReal procContrib, const Cmpnts& acceptorCoord, Cmpnts*** centroids, PetscReal& minDist) 
{
    Dcell closestDonor;
    closestDonor.rank   = -1;
    closestDonor.indi   = -1;
    closestDonor.indj   = -1;
    closestDonor.indk   = -1;
    closestDonor.dist2p = minDist; 


    // check if the node is null
    if (node == nullptr) 
    {
        PetscPrintf(PETSC_COMM_SELF, "Error: Octree node is null.\n");
        return closestDonor;
    }

    // check if the acceptor is outside the node bounds (this is essential to narrow down the search)
    if
    (
        acceptorCoord.x < node->minBounds.x || acceptorCoord.x >= node->maxBounds.x ||
        acceptorCoord.y < node->minBounds.y || acceptorCoord.y >= node->maxBounds.y ||
        acceptorCoord.z < node->minBounds.z || acceptorCoord.z >= node->maxBounds.z
    ) 
    {
        return closestDonor;
    }

    // check donor cells in this node
    for(PetscInt c=0; c<node->donorCells.size(); c++)
    {
        cellIds cell = node->donorCells[c];

        // Calculate distance to acceptor
        Cmpnts centroid = centroids[cell.k][cell.j][cell.i];
        PetscReal dist = sqrt(pow(centroid.x - acceptorCoord.x - procContrib, 2) +
                              pow(centroid.y - acceptorCoord.y - procContrib, 2) +
                              pow(centroid.z - acceptorCoord.z - procContrib, 2));

        if (dist < minDist) 
        {
            minDist             = dist + procContrib;
            closestDonor.indi   = cell.i;
            closestDonor.indj   = cell.j;
            closestDonor.indk   = cell.k;
            closestDonor.dist2p = dist;
            closestDonor.rank   = 1;
        }
    }

    // recursively search child nodes
    for (int i = 0; i < 8; i++) 
    {
        if (node->children[i] != nullptr) 
        {
            PetscReal childMinDist = minDist;
            Dcell childClosest = searchOctree(node->children[i], procContrib, acceptorCoord, centroids, childMinDist);

            if (childClosest.rank != -1 && childClosest.dist2p < minDist) 
            {
                closestDonor = childClosest;
                minDist = childClosest.dist2p;
            }
        }
    }

    return closestDonor;
}

//***************************************************************************************************************//

PetscErrorCode findClosestDonorC2P_Bins(mesh_ *meshDonor, mesh_ *meshAcceptor, PetscInt donorId)
{
    overset_         *os  = meshAcceptor->access->os;
    DM               da   = meshDonor->da, fda = meshDonor->fda;
    DMDALocalInfo    info = meshDonor->info;
    PetscMPIInt      rankD, sizeD, rankA, sizeA;

    MPI_Comm_size(meshDonor->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshDonor->MESH_COMM, &rankD);

    MPI_Comm_size(meshAcceptor->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshAcceptor->MESH_COMM, &rankA);

    std::vector<Acell> aCell = os->aCellHc;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcHc;
    std::vector<std::vector<PetscInt>> lAcellProcMat(sizeA);

    Cmpnts  ***cent;
    PetscReal ***aj;
    DMDAVecGetArray(fda, meshDonor->lCent, &cent);
    DMDAVecGetArray(da, meshDonor->lAj, &aj);

    PetscInt lxs = info.xs, lxe = info.xs + info.xm;
    PetscInt lys = info.ys, lye = info.ys + info.ym;
    PetscInt lzs = info.zs, lze = info.zs + info.zm;
    PetscInt mx  = info.mx, my = info.my, mz = info.mz;

    if (lxs == 0) lxs++; if (lxe == mx) lxe--;
    if (lys == 0) lys++; if (lye == my) lye--;
    if (lzs == 0) lzs++; if (lze == mz) lze--;

    // Estimate local average donor cell spacing
    PetscReal dxSum = 0.0, dySum = 0.0, dzSum = 0.0;
    PetscInt count = 0;
    for (PetscInt k = lzs; k < lze; k++)
    for (PetscInt j = lys; j < lye; j++)
    for (PetscInt i = lxs; i < lxe; i++) 
    {
        PetscReal cellSize = pow(1.0/aj[k][j][i], 1.0/3.0);
        dxSum += cellSize;
        dySum += cellSize;
        dzSum += cellSize;
        count++;
    }

    PetscReal localSpacing[3] = {count > 0 ? dxSum/count : 0.0, 
                                count > 0 ? dySum/count : 0.0, 
                                count > 0 ? dzSum/count : 0.0};
    PetscReal globalSpacing[3];

    // Gather global maximum spacing
    MPI_Allreduce(&localSpacing[0], &globalSpacing[0], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);
    MPI_Allreduce(&localSpacing[1], &globalSpacing[1], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);
    MPI_Allreduce(&localSpacing[2], &globalSpacing[2], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);

    PetscReal binSize = 2.0 * PetscMax(PetscMax(globalSpacing[0], globalSpacing[1]), globalSpacing[2]);

    // Build bins
    std::unordered_map<BinIndex, std::vector<std::tuple<PetscInt, PetscInt, PetscInt>>> bins;

    for (PetscInt k = lzs; k < lze; k++)
    for (PetscInt j = lys; j < lye; j++)
    for (PetscInt i = lxs; i < lxe; i++) {
        BinIndex idx;
        idx.ix = floor(cent[k][j][i].x / binSize);
        idx.iy = floor(cent[k][j][i].y / binSize);
        idx.iz = floor(cent[k][j][i].z / binSize);
        bins[idx].emplace_back(k, j, i);
    }

    // Resize vectors
    os->closestDonorHc.resize(aCell.size());
    os->AcellProcMatHc.resize(sizeA);


    for(PetscInt b = 0; b < sizeA; b++) 
    {
        lAcellProcMat[b].resize(sizeD);
        os->AcellProcMatHc[b].resize(sizeD);
    }

    // initial local processor matrix to previous value
    for (PetscInt b = 0; b < sizeA; b++) 
    {
        if (!os->AcellProcMatHc[b].empty()) 
        {
            for (PetscInt m = 0; m < sizeD; m++) 
            {
                if (os->AcellProcMatHc[b][m] != MPI_UNDEFINED) 
                {
                    lAcellProcMat[b][m] = os->AcellProcMatHc[b][m];
                }
            }
        }
    }

    // Find closest donor for each acceptor cell
    for(PetscInt b = 0; b < aCell.size(); b++) 
    {
        if(os->aCellHc[b].donorId == donorId)
        {
            Dcell dCell;
            dCell.rank = -1;
            PetscReal maxPerturb = 1e-10;
            PetscReal procContrib = maxPerturb * ((PetscReal)rankD + 1) / (PetscReal)sizeD;
            PetscReal lminDist = 1e20;
            std::vector<PetscInt> indices {0, 0, 0};
    
            BinIndex query;
            query.ix = floor(aCell[b].coorx / binSize); // Remove procContrib from bin index
            query.iy = floor(aCell[b].coory / binSize);
            query.iz = floor(aCell[b].coorz / binSize);
    
            // Search neighboring bins (3x3x3 for robustness)
            for (int dx = -1; dx <= 1; dx++)
            for (int dy = -1; dy <= 1; dy++)
            for (int dz = -1; dz <= 1; dz++) 
            {
                BinIndex neighbor = {query.ix + dx, query.iy + dy, query.iz + dz};
                if (bins.find(neighbor) != bins.end()) 
                {
                    for (auto& idx : bins[neighbor]) {
                        PetscInt k = std::get<0>(idx);
                        PetscInt j = std::get<1>(idx);
                        PetscInt i = std::get<2>(idx);
                        PetscReal ds = sqrt(
                            (cent[k][j][i].x - aCell[b].coorx - procContrib) * (cent[k][j][i].x - aCell[b].coorx - procContrib) +
                            (cent[k][j][i].y - aCell[b].coory - procContrib) * (cent[k][j][i].y - aCell[b].coory - procContrib) +
                            (cent[k][j][i].z - aCell[b].coorz - procContrib) * (cent[k][j][i].z - aCell[b].coorz - procContrib)
                        );
                        if (ds < lminDist) {
                            lminDist = ds + procContrib;
                            indices[0] = k;
                            indices[1] = j;
                            indices[2] = i;
                        }
                    }
                }
            }
    
            PetscReal gminDist;
            MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPI_MIN, meshDonor->MESH_COMM);
    
            if (lminDist == gminDist && lminDist < 1e19) 
            {
                dCell.indi = indices[2];
                dCell.indj = indices[1];
                dCell.indk = indices[0];
                dCell.dist2p = gminDist;
                dCell.rank = rankD;
                os->closestDonorHc[b] = dCell;
                lAcellProcMat[aCell[b].rank][rankD] = 1;
            }
    
            PetscReal lClosestSize = (lminDist < 1e19) ? pow(1./aj[indices[0]][indices[1]][indices[2]], 1./3.) : 0.0;
            MPI_Allreduce(&lClosestSize, &os->aCellHc[b].cell_size, 1, MPIU_REAL, MPI_SUM, meshDonor->MESH_COMM);
            MPI_Allreduce(&dCell.rank, &os->closestDonorHc[b].rank, 1, MPI_INT, MPI_MAX, meshDonor->MESH_COMM);
    
            // Check if no valid donor was found
            if (os->aCellHc[b].cell_size == 0.0 && rankD == 0) 
            {
                PetscPrintf(PETSC_COMM_SELF, "Warning: No donor cell found for acceptor cell %d\n", b);
            }
        }
    }

    // Finalize processor matrix
    for(PetscInt b = 0; b < sizeA; b++) 
    {
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMatHc[b][0], sizeD, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    }

    if (sizeA != sizeD) {
        char error[512];
        sprintf(error, "Meshes must have same number of processors.\n");
        fatalErrorInFunction("findClosestDonor", error);
    }

    for(PetscInt b = 0; b < sizeA; b++) 
    {
        if (NumAcellPerProc[b] != 0) os->AcellProcMatHc[b][b] = 1;
        for(PetscInt m = 0; m < sizeD; m++) 
        {
            if (os->AcellProcMatHc[b][m] == 0) os->AcellProcMatHc[b][m] = MPI_UNDEFINED;
        }
    }

    // Cleanup
    std::vector<Acell>().swap(aCell);
    std::vector<PetscInt>().swap(NumAcellPerProc);
    std::vector<std::vector<PetscInt>>().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, meshDonor->lCent, &cent);
    DMDAVecRestoreArray(da, meshDonor->lAj, &aj);

    return 0;
}

PetscErrorCode findClosestDonorC2P(mesh_ *meshDonor, mesh_ *meshAcceptor, PetscInt donorId)
{
    overset_         *os  = meshAcceptor->access->os;
    DM               da   = meshDonor->da, fda = meshDonor->fda;
    DMDALocalInfo info = meshDonor->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Vec           Coor;

    Cmpnts        ***cent;
    Cmpnts        ***coor;
    PetscReal     ***aj;

    PetscMPIInt   rankD, sizeD, rankA, sizeA;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_size(meshDonor->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshDonor->MESH_COMM, &rankD);

    MPI_Comm_size(meshAcceptor->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshAcceptor->MESH_COMM, &rankA);

    std::vector<Acell> aCell = os->aCellHc;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcHc;
    std::vector<std::vector<PetscInt>> lAcellProcMat(sizeA);

    DMDAVecGetArray(fda, meshDonor->lCent, &cent);
    DMDAVecGetArray(da, meshDonor->lAj, &aj);
    DMGetCoordinatesLocal(meshDonor->da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // there has to be not gap between adjacent octrees otherwise we are loosing cells
    PetscInt bsz = zs; if(zs!=0) bsz = bsz - 1;
    PetscInt bsy = ys; if(ys!=0) bsy = bsy - 1;
    PetscInt bsx = xs; if(xs!=0) bsx = bsx - 1;

    // find min and max bounds for this processor (in terms of points coordinates)
    Cmpnts minBounds = {coor[bsz  ][bsy  ][bsx  ].x, coor[bsz  ][bsy  ][bsx  ].y, coor[bsz  ][bsy  ][bsx  ].z};
    Cmpnts maxBounds = {coor[lze-1][lye-1][lxe-1].x, coor[lze-1][lye-1][lxe-1].y, coor[lze-2][lye-1][lxe-1].z};

    // build the octree
    PetscInt maxDepth        = 10;   // Maximum depth of the octree
    PetscInt maxCellsPerNode = 1000; // Maximum cells per leaf node
    OctreeNode *root         = new OctreeNode(minBounds, maxBounds);
    buildOctree(root, cent, lxs, lxe, lys, lye, lzs, lze, maxDepth, maxCellsPerNode);

    // Resize vectors
    os->closestDonorHc.resize(aCell.size());
    os->AcellProcMatHc.resize(sizeA);

    for(PetscInt b = 0; b < sizeA; b++) 
    {
        lAcellProcMat[b].resize(sizeD);
        os->AcellProcMatHc[b].resize(sizeD);
    }

    // initial local processor matrix to previous value
    for (PetscInt b = 0; b < sizeA; b++) 
    {
        if (!os->AcellProcMatHc[b].empty()) 
        {
            for (PetscInt m = 0; m < sizeD; m++) 
            {
                if (os->AcellProcMatHc[b][m] != MPI_UNDEFINED) 
                {
                    lAcellProcMat[b][m] = os->AcellProcMatHc[b][m];
                }
            }
        }
    }

    // Find closest donor for each acceptor cell
    for(PetscInt b = 0; b < aCell.size(); b++) 
    {
        if(os->aCellHc[b].donorId == donorId)
        {
            Dcell dCell, dCellLocal;
            PetscReal maxPerturb   = 1e-10;
            PetscReal procContrib  = maxPerturb * ((PetscReal)rankD + 1) / (PetscReal)sizeD;
            
            // initialize to huge
            PetscReal lminDist     = 1e20;
    
            // Search the octree for the closest donor cell
            Cmpnts acceptorCoord   = nSetFromComponents(aCell[b].coorx, aCell[b].coory, aCell[b].coorz);

            // exclude acceptor cell outside of this processor bounds (ocree search is useless)
            if
            (
                acceptorCoord.x >= root->minBounds.x && acceptorCoord.x < root->maxBounds.x &&
                acceptorCoord.y >= root->minBounds.y && acceptorCoord.y < root->maxBounds.y &&
                acceptorCoord.z >= root->minBounds.z && acceptorCoord.z < root->maxBounds.z
            )
            {
                dCellLocal           = searchOctree(root, procContrib, acceptorCoord, cent, lminDist);
                lminDist             = dCellLocal.dist2p;
            }

            PetscReal gminDist;
            MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPI_MIN, meshDonor->MESH_COMM);
            
            if (lminDist == gminDist && lminDist < 1e19) 
            {
                dCell.indi            = dCellLocal.indi;
                dCell.indj            = dCellLocal.indj;
                dCell.indk            = dCellLocal.indk;  
                dCell.dist2p          = gminDist;
                dCell.rank            = rankD;

                os->closestDonorHc[b] = dCell;
                lAcellProcMat[aCell[b].rank][rankD] = 1;
            }
            else 
            {
                dCell.rank            = -1;
            }
    
            PetscReal lClosestSize = (lminDist < 1e19) ? pow(1./aj[dCellLocal.indk][dCellLocal.indj][dCellLocal.indi], 1./3.) : 0.0;
            MPI_Allreduce(&lClosestSize, &os->aCellHc[b].cell_size, 1, MPIU_REAL, MPI_SUM, meshDonor->MESH_COMM);
            MPI_Allreduce(&dCell.rank, &os->closestDonorHc[b].rank, 1, MPI_INT, MPI_MAX, meshDonor->MESH_COMM);
    
            // Check if no valid donor was found
            if (os->aCellHc[b].cell_size == 0.0 && rankD == 0) 
            {
                PetscPrintf(PETSC_COMM_SELF, "Warning: No donor cell found for acceptor cell %d\n", b);
            }
        }
    }

    // Finalize processor matrix
    for(PetscInt b = 0; b < sizeA; b++) 
    {
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMatHc[b][0], sizeD, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    }

    if (sizeA != sizeD) 
    {
        char error[512];
        sprintf(error, "Meshes must have same number of processors.\n");
        fatalErrorInFunction("findClosestDonor", error);
    }

    for(PetscInt b = 0; b < sizeA; b++) 
    {
        if (NumAcellPerProc[b] != 0) os->AcellProcMatHc[b][b] = 1;
        for(PetscInt m = 0; m < sizeD; m++) 
        {
            if (os->AcellProcMatHc[b][m] == 0) os->AcellProcMatHc[b][m] = MPI_UNDEFINED;
        }
    }

    // Cleanup
    delete root;
    std::vector<Acell>().swap(aCell);
    std::vector<PetscInt>().swap(NumAcellPerProc);
    std::vector<std::vector<PetscInt>>().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(fda, meshDonor->lCent, &cent);
    DMDAVecRestoreArray(da, meshDonor->lAj, &aj);

    return 0;
}

// using octree to find closest donor cell for acceptor cells
PetscErrorCode findClosestDonorP2C(mesh_ *meshDonor, mesh_ *meshAcceptor)
{
    overset_      *os  = meshAcceptor->access->os;
    DM            da   = meshDonor->da, fda = meshDonor->fda;
    DMDALocalInfo info = meshDonor->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Vec           Coor;

    Cmpnts        ***cent;
    Cmpnts        ***coor;
    PetscReal     ***aj;

    PetscMPIInt   rankD, sizeD, rankA, sizeA;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_size(meshDonor->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshDonor->MESH_COMM, &rankD);

    MPI_Comm_size(meshAcceptor->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshAcceptor->MESH_COMM, &rankA);

    std::vector<Acell> aCell = os->aCellDb;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcDb;
    std::vector<std::vector<PetscInt>> lAcellProcMat(sizeA);

    DMDAVecGetArray(fda, meshDonor->lCent, &cent);
    DMDAVecGetArray(da, meshDonor->lAj, &aj);
    DMGetCoordinatesLocal(meshDonor->da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // there has to be not gap between adjacent octrees otherwise we are loosing cells
    PetscInt bsz = zs; if(zs!=0) bsz = bsz - 1;
    PetscInt bsy = ys; if(ys!=0) bsy = bsy - 1;
    PetscInt bsx = xs; if(xs!=0) bsx = bsx - 1;

    // find min and max bounds for this processor (in terms of points coordinates)
    Cmpnts minBounds = {coor[bsz  ][bsy  ][bsx  ].x, coor[bsz  ][bsy  ][bsx  ].y, coor[bsz  ][bsy  ][bsx  ].z};
    Cmpnts maxBounds = {coor[lze-1][lye-1][lxe-1].x, coor[lze-1][lye-1][lxe-1].y, coor[lze-2][lye-1][lxe-1].z};

    // build the octree
    PetscInt maxDepth        = 10;   // Maximum depth of the octree
    PetscInt maxCellsPerNode = 1000; // Maximum cells per leaf node
    OctreeNode *root         = new OctreeNode(minBounds, maxBounds);
    buildOctree(root, cent, lxs, lxe, lys, lye, lzs, lze, maxDepth, maxCellsPerNode);

    // Resize vectors
    os->closestDonorDb.resize(aCell.size());
    os->AcellProcMatDb.resize(sizeA);

    for (PetscInt b = 0; b < sizeA; b++) {
        lAcellProcMat[b].resize(sizeD);
        os->AcellProcMatDb[b].resize(sizeD);
    }

    // Find closest donor for each acceptor cell using the octree
    for (PetscInt b = 0; b < aCell.size(); b++) 
    {
        Dcell dCell, dCellLocal;
        PetscReal maxPerturb   = 1e-10;
        PetscReal procContrib  = maxPerturb * ((PetscReal)rankD + 1) / (PetscReal)sizeD;
        
        // initialize to huge
        PetscReal lminDist     = 1e20;

        // Search the octree for the closest donor cell
        Cmpnts acceptorCoord   = nSetFromComponents(aCell[b].coorx, aCell[b].coory, aCell[b].coorz);

        // exclude acceptor cell outside of this processor bounds (ocree search is useless)
        if
        (
            acceptorCoord.x >= root->minBounds.x && acceptorCoord.x < root->maxBounds.x &&
            acceptorCoord.y >= root->minBounds.y && acceptorCoord.y < root->maxBounds.y &&
            acceptorCoord.z >= root->minBounds.z && acceptorCoord.z < root->maxBounds.z
        )
        {
            dCellLocal           = searchOctree(root, procContrib, acceptorCoord, cent, lminDist);
            lminDist             = dCellLocal.dist2p;
        }
        
        PetscReal gminDist;
        MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPI_MIN, meshDonor->MESH_COMM);

        if (lminDist == gminDist && lminDist < 1e19) 
        {
            dCell.indi            = dCellLocal.indi;
            dCell.indj            = dCellLocal.indj;
            dCell.indk            = dCellLocal.indk;  
            dCell.dist2p          = gminDist;
            dCell.rank            = rankD;

            os->closestDonorDb[b] = dCell;
            lAcellProcMat[aCell[b].rank][rankD] = 1;
        }
        else 
        {
            dCell.rank            = -1;
        }

        PetscReal lClosestSize = (lminDist < 1e19) ? pow(1./aj[dCellLocal.indk][dCellLocal.indj][dCellLocal.indi], 1./3.) : 0.0;
        
        MPI_Allreduce(&lClosestSize, &os->aCellDb[b].cell_size, 1, MPIU_REAL, MPI_SUM, meshDonor->MESH_COMM);
        MPI_Allreduce(&dCell.rank, &os->closestDonorDb[b].rank, 1, MPI_INT, MPI_MAX, meshDonor->MESH_COMM);

        // check if no valid donor was found
        if (os->aCellDb[b].cell_size == 0.0 && rankD == 0) 
        {
            PetscPrintf(PETSC_COMM_SELF, "Warning: No donor cell found for acceptor cell %d\n", b);
        }
    }

    // Finalize processor matrix
    for (PetscInt b = 0; b < sizeA; b++) 
    {
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMatDb[b][0], sizeD, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    }

    if (sizeA != sizeD) 
    {
        char error[512];
        sprintf(error, "Meshes must have the same number of processors.\n");
        fatalErrorInFunction("findClosestDonorP2C", error);
    }

    for (PetscInt b = 0; b < sizeA; b++) 
    {
        if (NumAcellPerProc[b] != 0) os->AcellProcMatDb[b][b] = 1;
        for (PetscInt m = 0; m < sizeD; m++) {
            if (os->AcellProcMatDb[b][m] == 0) os->AcellProcMatDb[b][m] = MPI_UNDEFINED;
        }
    }

    // Cleanup
    delete root;
    std::vector<Acell>().swap(aCell);
    std::vector<PetscInt>().swap(NumAcellPerProc);
    std::vector<std::vector<PetscInt>>().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(fda, meshDonor->lCent, &cent);
    DMDAVecRestoreArray(da, meshDonor->lAj, &aj);

    return 0;
}


PetscErrorCode findClosestDonorP2C_Bins(mesh_ *meshDonor, mesh_ *meshAcceptor)
{
    overset_         *os  = meshAcceptor->access->os;
    DM               da   = meshDonor->da, fda = meshDonor->fda;
    DMDALocalInfo    info = meshDonor->info;
    PetscMPIInt      rankD, sizeD, rankA, sizeA;

    MPI_Comm_size(meshDonor->MESH_COMM, &sizeD);
    MPI_Comm_rank(meshDonor->MESH_COMM, &rankD);

    MPI_Comm_size(meshAcceptor->MESH_COMM, &sizeA);
    MPI_Comm_rank(meshAcceptor->MESH_COMM, &rankA);

    std::vector<Acell> aCell = os->aCellDb;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProcDb;
    std::vector<std::vector<PetscInt>> lAcellProcMat(sizeA);

    Cmpnts  ***cent;
    PetscReal ***aj;
    DMDAVecGetArray(fda, meshDonor->lCent, &cent);
    DMDAVecGetArray(da, meshDonor->lAj, &aj);

    PetscInt lxs = info.xs, lxe = info.xs + info.xm;
    PetscInt lys = info.ys, lye = info.ys + info.ym;
    PetscInt lzs = info.zs, lze = info.zs + info.zm;
    PetscInt mx  = info.mx, my = info.my, mz = info.mz;

    if (lxs == 0) lxs++; if (lxe == mx) lxe--;
    if (lys == 0) lys++; if (lye == my) lye--;
    if (lzs == 0) lzs++; if (lze == mz) lze--;

    // Estimate local average donor cell spacing
    PetscReal dxSum = 0.0, dySum = 0.0, dzSum = 0.0;
    PetscInt count = 0;
    for (PetscInt k = lzs; k < lze; k++)
    for (PetscInt j = lys; j < lye; j++)
    for (PetscInt i = lxs; i < lxe; i++) 
    {
        PetscReal cellSize = pow(1.0/aj[k][j][i], 1.0/3.0);
        dxSum += cellSize;
        dySum += cellSize;
        dzSum += cellSize;
        count++;
    }

    PetscReal localSpacing[3] = {count > 0 ? dxSum/count : 0.0, 
                                count > 0 ? dySum/count : 0.0, 
                                count > 0 ? dzSum/count : 0.0};
    PetscReal globalSpacing[3];

    // Gather global maximum spacing
    MPI_Allreduce(&localSpacing[0], &globalSpacing[0], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);
    MPI_Allreduce(&localSpacing[1], &globalSpacing[1], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);
    MPI_Allreduce(&localSpacing[2], &globalSpacing[2], 1, MPIU_REAL, MPI_MAX, meshDonor->MESH_COMM);

    PetscReal binSize = 2.0 * PetscMax(PetscMax(globalSpacing[0], globalSpacing[1]), globalSpacing[2]);

    // Build bins
    std::unordered_map<BinIndex, std::vector<std::tuple<PetscInt, PetscInt, PetscInt>>> bins;

    for (PetscInt k = lzs; k < lze; k++)
    for (PetscInt j = lys; j < lye; j++)
    for (PetscInt i = lxs; i < lxe; i++) 
    {
        BinIndex idx;
        idx.ix = floor(cent[k][j][i].x / binSize);
        idx.iy = floor(cent[k][j][i].y / binSize);
        idx.iz = floor(cent[k][j][i].z / binSize);
        bins[idx].emplace_back(k, j, i);
    }

    // Resize vectors
    os->closestDonorDb.resize(aCell.size());
    os->AcellProcMatDb.resize(sizeA);

    for(PetscInt b = 0; b < sizeA; b++) 
    {
        lAcellProcMat[b].resize(sizeD);
        os->AcellProcMatDb[b].resize(sizeD);
    }

    // Find closest donor for each acceptor cell
    for(PetscInt b = 0; b < aCell.size(); b++) 
    {
        Dcell dCell;
        dCell.rank = -1;
        PetscReal maxPerturb = 1e-10;
        PetscReal procContrib = maxPerturb * ((PetscReal)rankD + 1) / (PetscReal)sizeD;
        PetscReal lminDist = 1e20;
        std::vector<PetscInt> indices {0, 0, 0};

        BinIndex query;
        query.ix = floor(aCell[b].coorx / binSize); // Remove procContrib from bin index
        query.iy = floor(aCell[b].coory / binSize);
        query.iz = floor(aCell[b].coorz / binSize);

        // Search neighboring bins (5x5x5 for robustness)
        for (int dx = -1; dx <= 1; dx++)
        for (int dy = -1; dy <= 1; dy++)
        for (int dz = -1; dz <= 1; dz++) 
        {
            BinIndex neighbor = {query.ix + dx, query.iy + dy, query.iz + dz};
            if (bins.find(neighbor) != bins.end()) 
            {
                for (auto& idx : bins[neighbor]) 
                {
                    PetscInt k = std::get<0>(idx);
                    PetscInt j = std::get<1>(idx);
                    PetscInt i = std::get<2>(idx);
                    PetscReal ds = sqrt(
                        (cent[k][j][i].x - aCell[b].coorx - procContrib) * (cent[k][j][i].x - aCell[b].coorx - procContrib) +
                        (cent[k][j][i].y - aCell[b].coory - procContrib) * (cent[k][j][i].y - aCell[b].coory - procContrib) +
                        (cent[k][j][i].z - aCell[b].coorz - procContrib) * (cent[k][j][i].z - aCell[b].coorz - procContrib)
                    );
                    if (ds < lminDist) 
                    {
                        lminDist = ds + procContrib;
                        indices[0] = k;
                        indices[1] = j;
                        indices[2] = i;
                    }
                }
            }
        }

        PetscReal gminDist;
        MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPI_MIN, meshDonor->MESH_COMM);

        if (lminDist == gminDist && lminDist < 1e19) 
        {
            dCell.indi = indices[2];
            dCell.indj = indices[1];
            dCell.indk = indices[0];
            dCell.dist2p = gminDist;
            dCell.rank = rankD;
            os->closestDonorDb[b] = dCell;
            lAcellProcMat[aCell[b].rank][rankD] = 1;
        }

        printf("Rank %d: Acceptor cell %ld, closest donor cell (%ld, %ld, %ld) with distance %f\n", dCell.rank, b, dCell.indi, dCell.indj, dCell.indk, gminDist);

        PetscReal lClosestSize = (lminDist < 1e19) ? pow(1./aj[indices[0]][indices[1]][indices[2]], 1./3.) : 0.0;
        MPI_Allreduce(&lClosestSize, &os->aCellDb[b].cell_size, 1, MPIU_REAL, MPI_SUM, meshDonor->MESH_COMM);
        MPI_Allreduce(&dCell.rank, &os->closestDonorDb[b].rank, 1, MPI_INT, MPI_MAX, meshDonor->MESH_COMM);

        // Check if no valid donor was found
        if (os->aCellDb[b].cell_size == 0.0 && rankD == 0) {
            PetscPrintf(PETSC_COMM_SELF, "Warning: No donor cell found for acceptor cell %d\n", b);
        }
    }

    // Finalize processor matrix
    for(PetscInt b = 0; b < sizeA; b++) 
    {
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMatDb[b][0], sizeD, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    }

    if (sizeA != sizeD) {
        char error[512];
        sprintf(error, "Meshes must have same number of processors.\n");
        fatalErrorInFunction("findClosestDonor", error);
    }

    for(PetscInt b = 0; b < sizeA; b++) {
        if (NumAcellPerProc[b] != 0) os->AcellProcMatDb[b][b] = 1;
        for(PetscInt m = 0; m < sizeD; m++) {
            if (os->AcellProcMatDb[b][m] == 0) os->AcellProcMatDb[b][m] = MPI_UNDEFINED;
        }
    }

    // Cleanup
    std::vector<Acell>().swap(aCell);
    std::vector<PetscInt>().swap(NumAcellPerProc);
    std::vector<std::vector<PetscInt>>().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, meshDonor->lCent, &cent);
    DMDAVecRestoreArray(da, meshDonor->lAj, &aj);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode updateAcceptorCoordinates(overset_ *os)
{
  oversetMotion *osetMotion = os->oMotion;

  if(osetMotion->setMotion)
  {
    if(osetMotion->motionType == "Translation")
    {
      oversetMeshTranslation(os);
    }
    else if (osetMotion->motionType == "Rotation")
    {
      PetscPrintf(os->access->mesh->MESH_COMM, "Rotating the overset mesh\n");
      exit(0);
    }
    else
    {
     char error[512];
      sprintf(error, "Only translation and rotation motion is available presently \n");
      fatalErrorInFunction("readOversetParameters",  error);
    }
  }
  else
  {
    // check if ibm is moving from ibm motion
    // temporary variable for now to be added in IBM
    PetscInt ibmMotion = 0;

    if(!ibmMotion)
    {
     char error[512];
      sprintf(error, "IBM is not moving. Use static overset as dynamic overset motion not required.\n");
      fatalErrorInFunction("readOversetParameters",  error);
    }
  }
  return(0);
}

//***************************************************************************************************************//

PetscErrorCode oversetMeshTranslation(overset_ *os)
{
  mesh_           *mesh = os->access->mesh;
  DM              da = mesh->da, fda = mesh->fda;
  DMDALocalInfo	  info = mesh->info;
  PetscInt	      xs = info.xs, xe = info.xs + info.xm;
  PetscInt        ys = info.ys, ye = info.ys + info.ym;
  PetscInt	      zs = info.zs, ze = info.zs + info.zm;
  PetscInt	      mx = info.mx, my = info.my, mz = info.mz;
  PetscInt        lxs, lxe, lys, lye, lzs, lze;

  PetscInt             k, j, i;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  Vec             Coor;
  Cmpnts          ***coor, presVel;

  oversetMotion *osetMotion = os->oMotion;

  DMGetCoordinates(da, &Coor);
  DMDAVecGetArray(fda, Coor, &coor);

  //obtain the prescribed velocity
  presVel = nSet(osetMotion->prescribedVel);

  for (k=zs; k<lze; k++)
      for (j=ys; j<lye; j++)
          for (i=xs; i<lxe; i++){
            // move mesh
          }
  DMDAVecRestoreArray(fda, Coor, &coor);

  PetscPrintf(mesh->MESH_COMM, "Dynamic overset translation: This function is not complete. Exiting ...\n");
  exit(0);
  return 0;
}

//***************************************************************************************************************//

PetscErrorCode oversetIbmSearch(ibm_ *ibm)
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

    PetscReal     ***meshTag, ***gmeshTag;
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

        DMDAVecGetArray(fda, mesh->lCent, &cent);
        DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

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

                // find the search cell were the fluid node is located
                sCell.i = floor((cent[k][j][i].x - ibBox->xmin) / sBox->dcx);
                sCell.j = floor((cent[k][j][i].y - ibBox->ymin) / sBox->dcy);
                sCell.k = floor((cent[k][j][i].z - ibBox->zmin) / sBox->dcz);

                // do the ray casting test to check if a cell is inside or outside an IBM body
                PetscReal val;
                val = rayCastingTest(cent[k][j][i], ibMsh, sCell, sBox, ibBox, searchCellList);
                meshTag[k][j][i] = PetscMax(meshTag[k][j][i], val);
            }
        }

        DMDAVecRestoreArray(fda, mesh->lCent, &cent);
        DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

        MPI_Barrier(mesh->MESH_COMM);

        DMLocalToLocalBegin(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);
        DMLocalToLocalEnd(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);

        // set nvert at solid fluid intersection IB Nodes
        DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {

            if (meshTag[k][j][i] < 0)
            {
                meshTag[k][j][i] = 0;
            }

            ip = (i < mx - 1 ? (i + 1) : (i));
            im = (i > 0 ? (i - 1) : (i));

            jp = (j < my - 1 ? (j + 1) : (j));
            jm = (j > 0 ? (j - 1) : (j));

            kp = (k < mz - 1 ? (k + 1) : (k));
            km = (k > 0 ? (k - 1) : (k));

            if ((PetscInt) (meshTag[k][j][i] + 0.5) != 4)
            {
                for (kk = km; kk < kp + 1; kk++)
                for (jj = jm; jj < jp + 1; jj++)
                for (ii = im; ii < ip + 1; ii++)
                {
                    if ((PetscInt) (meshTag[kk][jj][ii] + 0.5) == 4)
                    {
                        meshTag[k][j][i] = PetscMax(2.0, meshTag[k][j][i]);
                    }
                }
            }
        }

        DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

        DMLocalToLocalBegin(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);
        DMLocalToLocalEnd(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);

        // nvert values of 4 for solid and 2 for IB fluid nodes where used for the current body
        // to differentiate it from other bodies.
        // reset back to nvert values of 3 for solid and 1 for IB fluid
        DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

        // Back to the old nvert 3 and 1
        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {
            if ((PetscInt) (meshTag[k][j][i] + 0.5) == 2)
            {
                meshTag[k][j][i] = 1;
            }
            if ((PetscInt) (meshTag[k][j][i] + 0.5) == 4)
            {
                meshTag[k][j][i] = 3;
            }

        }

        DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

        DMLocalToLocalBegin(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);
        DMLocalToLocalEnd(da, mesh->lmeshTag, INSERT_VALUES, mesh->lmeshTag);
    }

    //ibm nvert cleanup - make solid, ibm fluid cells that dont have even one fluid cell around it
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da, mesh->meshTag, &gmeshTag);

    for (k = lzs; k < lze; k++)
    for (j = lys; j < lye; j++)
    for (i = lxs; i < lxe; i++)
    {
        //set the global nvert
        gmeshTag[k][j][i] = meshTag[k][j][i];
    }

    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da, mesh->meshTag, &gmeshTag);

    DMGlobalToLocalBegin(da, mesh->meshTag, INSERT_VALUES, mesh->lmeshTag);
    DMGlobalToLocalEnd(da, mesh->meshTag, INSERT_VALUES, mesh->lmeshTag);

    MPI_Barrier(mesh->MESH_COMM);
    
    return 0;
}

PetscErrorCode computeOversetIBMElementNormal(ibm_ *ibm)
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

void defineStruct_Acell(MPI_Datatype *tstype) {
    const PetscInt    count = 11; 
    int               blocklens[count];
    MPI_Aint          disps[count];

    for (PetscInt i = 0; i < count; i++) {
        blocklens[i] = 1; // Each field is a single element
    }

    MPI_Datatype types[count] = {
        MPIU_INT,   // indi
        MPIU_INT,   // indj
        MPIU_INT,   // indk
        MPIU_REAL,  // coorx
        MPIU_REAL,  // coory
        MPIU_REAL,  // coorz
        MPIU_INT,   // rank (PetscMPIInt)
        MPIU_REAL,  // cell_size
        MPIU_INT,   // face
        MPIU_INT,   // donorId
        MPIU_INT    // parentCellId (new)
    };

    disps[0] = offsetof(Acell, indi);
    disps[1] = offsetof(Acell, indj);
    disps[2] = offsetof(Acell, indk);
    disps[3] = offsetof(Acell, coorx);
    disps[4] = offsetof(Acell, coory);
    disps[5] = offsetof(Acell, coorz);
    disps[6] = offsetof(Acell, rank);
    disps[7] = offsetof(Acell, cell_size);
    disps[8] = offsetof(Acell, face);
    disps[9] = offsetof(Acell, donorId);
    disps[10] = offsetof(Acell, parentCellId); // New offset

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
    return;
}

//***************************************************************************************************************//

//MPI operation function to find the sum the elements of the vector of structs
void sum_struct_Acell(void *in, void *inout, int *len, MPI_Datatype *type) {
    Acell *invals    = (Acell*)in;
    Acell *inoutvals = (Acell*)inout;

    for (PetscInt i = 0; i < *len; i++) {
        // If invals[i] is a valid cell, copy it to inoutvals[i]
        if (invals[i].rank >= 0) {
            inoutvals[i].indi = invals[i].indi;
            inoutvals[i].indj = invals[i].indj;
            inoutvals[i].indk = invals[i].indk;
            inoutvals[i].coorx = invals[i].coorx;
            inoutvals[i].coory = invals[i].coory;
            inoutvals[i].coorz = invals[i].coorz;
            inoutvals[i].rank = invals[i].rank;
            inoutvals[i].cell_size = invals[i].cell_size;
            inoutvals[i].face = invals[i].face;
            inoutvals[i].donorId = invals[i].donorId;
            inoutvals[i].parentCellId = invals[i].parentCellId; // New field
        }
        // Otherwise, keep inoutvals unchanged (it should already contain the valid cell or be default)
    }

    return;
}

// ************************************************************************************************* //
//Deprecated functions for other interpolations. Variables names need to be updated for integration.

// PetscErrorCode interpolateACellInvD(mesh_ *meshP, mesh_ *mesh)
// {
//     overset_         *os    = mesh->access->os;
//     ueqn_            *ueqn  = mesh->access->ueqn;
//     ueqn_            *ueqnP = meshP->access->ueqn;
//     teqn_            *teqn  = mesh->access->teqn;
//     teqn_            *teqnP = meshP->access->teqn;

//     DM               da1    = mesh->da, fda1 = mesh->fda;
//     DMDALocalInfo    info1  = mesh->info;
//     DM               da0    = meshP->da, fda0 = meshP->fda;
//     DMDALocalInfo    info0  = meshP->info;

//     PetscInt         xs = info1.xs, xe = info1.xs + info1.xm;
//     PetscInt         ys = info1.ys, ye = info1.ys + info1.ym;
//     PetscInt         zs = info1.zs, ze = info1.zs + info1.zm;
//     PetscInt         mx = info1.mx, my = info1.my, mz = info1.mz;

//     PetscInt         lxs, lxe, lys, lye, lzs, lze;
//     PetscInt         i, j, k, b, m, n;
//     PetscInt         ii, jj, kk;
//     PetscInt         ip, im, jp, jm, kp, km;

//     Cmpnts           ***lucat0, ***ucat1, lucart, gucart;

//     PetscReal        dist;
//     PetscReal        ***ltemp0, ***temp1, lT, gT;

//     PetscMPIInt      rank, size, rankP, sizeP;
//     PetscInt         sum_ind1 = 0;

//     MPI_Comm_size(mesh->MESH_COMM, &size);
//     MPI_Comm_rank(mesh->MESH_COMM, &rank);

//     MPI_Comm_size(meshP->MESH_COMM, &sizeP);
//     MPI_Comm_rank(meshP->MESH_COMM, &rankP);

//     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
//     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
//     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

//     std::vector<Acell> aCell = os->aCell;
//     std::vector<std::vector<Dcell>> dCell = os->dCell;
//     std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMat;
//     std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;
//     std::vector<std::vector<PetscReal>> DWeights = os->DWeights;

//     DMDAVecGetArray(fda0, ueqnP->lUcat, &lucat0);
//     DMDAVecGetArray(fda1, ueqn->Ucat, &ucat1);

//     if (mesh->access->flags->isTeqnActive)
//     {
//         DMDAVecGetArray(da0, teqnP->lTmprt, &ltemp0);
//         DMDAVecGetArray(da1, teqn->Tmprt, &temp1);
//     }

//     // to set the local rank within the new communicator
//     std::vector<PetscInt> local_rank_Vec;
//     PetscInt aCellLocalrank = 0;

//     // loop through the ranks
//     for(n = 0; n < size; n++){

//         for(m = 0; m < sizeP; m++){
//             if(AcellProcMat[n][m] == 1)
//                 local_rank_Vec.push_back(m);
//         }

//         // find the rank of aCell cells within the new communicator
//         for(m = 0; m < local_rank_Vec.size(); m++){
//             if(n == local_rank_Vec[m])
//                 aCellLocalrank = m;
//         }

//         std::vector<PetscInt> ().swap(local_rank_Vec);

//         if(NumAcellPerProc[n]!=0){

//             if(AcellProcMat[n][rankP] !=MPI_UNDEFINED){

//                 // loop through the aCell cells of a given processor n
//                 for(b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++){

//                     PetscReal   lsumwt = 0.;
//                     PetscReal   gsumwt = 0.;

//                     lucart.x = 0.; gucart.x = 0.;
//                     lucart.y = 0.; gucart.y = 0.;
//                     lucart.z = 0.; gucart.z = 0.;
//                     lT = 0;        gT = 0;

//                     // aCell cell index
//                     i = aCell[b].indi;
//                     j = aCell[b].indj;
//                     k = aCell[b].indk;

//                     // loop through the donor cells within the current processor - gnumdCell[b][rank]
//                     for(m = 0; m < dCell[b].size(); m++){

//                         kk = dCell[b][m].indk;
//                         jj = dCell[b][m].indj;
//                         ii = dCell[b][m].indi;

//                         dist = dCell[b][m].dist2p;

//                         lucart.x += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].x;
//                         lucart.y += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].y;
//                         lucart.z += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].z;

//                         if (mesh->access->flags->isTeqnActive)
//                         {
//                             lT += (1.0/(PetscMax(dist, 1e-10))) * ltemp0[kk][jj][ii];
//                         }

//                         lsumwt += 1.0/(PetscMax(dist, 1e-10));

//                     }

//                     // reduce the contribution of all valid processors to the local rank aCellLocalrank within the new communicator
//                     MPI_Reduce(&lucart, &gucart, 3, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
//                     MPI_Reduce(&lT, &gT, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
//                     MPI_Reduce(&lsumwt, &gsumwt, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);

//                     gucart.x /= gsumwt;
//                     gucart.y /= gsumwt;
//                     gucart.z /= gsumwt;
//                     gT /= gsumwt;

//                     if (rank == aCell[b].rank){

//                        if(aCell[b].face == 0)
//                        {
//                          ucat1[k][j][i].x = gucart.x;
//                          ucat1[k][j][i].y = gucart.y;
//                          ucat1[k][j][i].z = gucart.z;

//                          if (mesh->access->flags->isTeqnActive)
//                          {
//                              temp1[k][j][i] = gT;
//                          }
//                        }
//                        else
//                        {
//                          oversetContravariantBC(mesh, i, j, k, gucart, aCell[b].face);
//                        }

//                     }

//                 }
//             }

//             sum_ind1 +=NumAcellPerProc[n];

//         }

//     }

//     std::vector<Acell> ().swap(aCell);
//     std::vector<std::vector<Dcell>> ().swap(dCell);
//     std::vector<std::vector<PetscReal>> ().swap(DWeights);
//     std::vector<std::vector<PetscInt>> ().swap(AcellProcMat);
//     std::vector<PetscInt> ().swap(NumAcellPerProc);

//     DMDAVecRestoreArray(fda0, ueqnP->lUcat, &lucat0);
//     DMDAVecRestoreArray(fda1, ueqn->Ucat, &ucat1);

//     if (mesh->access->flags->isTeqnActive)
//     {
//         DMDAVecRestoreArray(da0, teqnP->lTmprt, &ltemp0);
//         DMDAVecRestoreArray(da1, teqn->Tmprt, &temp1);

//         DMGlobalToLocalBegin(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
//         DMGlobalToLocalEnd(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
//     }

//     DMGlobalToLocalBegin(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
//     DMGlobalToLocalEnd(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

//     DMGlobalToLocalBegin(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
//     DMGlobalToLocalEnd(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

//     return 0;
// }

// //***************************************************************************************************************//

// // split the global aCell vector based on processors. For the cells within each processor
// // interpolate from the donor cell connectivity and using the processor matrix (to know which processors have the donor cells)

// PetscErrorCode interpolateACellLS(mesh_ *meshP, mesh_ *mesh)
// {
//     overset_         *os    = mesh->access->os;
//     ueqn_            *ueqn  = mesh->access->ueqn;
//     ueqn_            *ueqnP = meshP->access->ueqn;
//     teqn_            *teqn  = mesh->access->teqn;
//     teqn_            *teqnP = meshP->access->teqn;
//     flags_           *flags = mesh->access->flags;
//     DM               da1    = mesh->da, fda1 = mesh->fda;
//     DMDALocalInfo    info1  = mesh->info;
//     DM               da0    = meshP->da, fda0 = meshP->fda;
//     DMDALocalInfo    info0  = meshP->info;

//     PetscInt         xs = info1.xs, xe = info1.xs + info1.xm;
//     PetscInt         ys = info1.ys, ye = info1.ys + info1.ym;
//     PetscInt         zs = info1.zs, ze = info1.zs + info1.zm;
//     PetscInt         mx = info1.mx, my = info1.my, mz = info1.mz;

//     PetscInt         lxs, lxe, lys, lye, lzs, lze;
//     PetscInt         i, j, k, b, m, n;
//     PetscInt         ii, jj, kk;

//     Cmpnts           ***lucat0, ***ucat1, lucart, gucart;

//     PetscReal        dist, ds;
//     PetscReal        ***nvert, ***ltemp0, ***temp1, lT, gT;
//     PetscMPIInt      rank, size, rankP, sizeP;
//     PetscInt         sum_ind1 = 0;

//     MPI_Comm_size(mesh->MESH_COMM, &size);
//     MPI_Comm_rank(mesh->MESH_COMM, &rank);

//     MPI_Comm_size(meshP->MESH_COMM, &sizeP);
//     MPI_Comm_rank(meshP->MESH_COMM, &rankP);

//     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
//     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
//     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

//     std::vector<Acell> aCell = os->aCell;
//     std::vector<std::vector<Dcell>> dCell = os->dCell;
//     std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMat;
//     std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;
//     std::vector<std::vector<PetscReal>> DWeights = os->DWeights;

//     DMDAVecGetArray(fda0, ueqnP->lUcat, &lucat0);
//     DMDAVecGetArray(fda1, ueqn->Ucat, &ucat1);

//     if (flags->isTeqnActive)
//     {
//         DMDAVecGetArray(da0, teqnP->lTmprt, &ltemp0);
//         DMDAVecGetArray(da1, teqn->Tmprt, &temp1);
//     }

//     // to set the local rank within the new communicator
//     std::vector<PetscInt> local_rank_Vec;
//     PetscInt aCellLocalrank = 0;

//     // loop through the ranks
//     for(n = 0; n < size; n++){

//         for(m = 0; m < sizeP; m++){
//             if(AcellProcMat[n][m] == 1)
//                 local_rank_Vec.push_back(m);
//         }

//         // find the rank of aCell cells within the new communicator
//         for(m = 0; m < local_rank_Vec.size(); m++){
//             if(n == local_rank_Vec[m])
//                 aCellLocalrank = m;
//         }

//         std::vector<PetscInt> ().swap(local_rank_Vec);

//         // proceed only if there are non zero acceptor cells in processor n
//         if(NumAcellPerProc[n]!=0){

//             // proceed only if the current processor has donor cells to processor n
//             if(AcellProcMat[n][rankP] !=MPI_UNDEFINED){

//                 // now data is exchanged only between the acceptor donor processors through the new communicator
//                 // loop through the aCell cells of a given processor n
//                 for(b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++){

//                     lucart.x = 0.; gucart.x = 0.;
//                     lucart.y = 0.; gucart.y = 0.;
//                     lucart.z = 0.; gucart.z = 0.;
//                     lT = 0;        gT = 0;

//                     // aCell cell index
//                     i = aCell[b].indi;
//                     j = aCell[b].indj;
//                     k = aCell[b].indk;

//                     // loop through the donor cells of the current acceptor
//                     for(m = 0; m < dCell[b].size(); m++){

//                         kk = dCell[b][m].indk;
//                         jj = dCell[b][m].indj;
//                         ii = dCell[b][m].indi;

//                         lucart.x += DWeights[b][m] * lucat0[kk][jj][ii].x;
//                         lucart.y += DWeights[b][m] * lucat0[kk][jj][ii].y;
//                         lucart.z += DWeights[b][m] * lucat0[kk][jj][ii].z;

//                         if (flags->isTeqnActive)
//                         {
//                             lT += DWeights[b][m] * ltemp0[kk][jj][ii];
//                         }

//                    }

//                     // reduce the contribution of all valid processors to the local rank aCellLocalrank within the new communicator
//                     MPI_Reduce(&lucart, &gucart, 3, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
//                     MPI_Reduce(&lT, &gT, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);

//                     if (rank == aCell[b].rank){

//                         if(aCell[b].face == 0)
//                         {
//                             ucat1[k][j][i].x = gucart.x;
//                             ucat1[k][j][i].y = gucart.y;
//                             ucat1[k][j][i].z = gucart.z;

//                             if (flags->isTeqnActive)
//                             {
//                                 temp1[k][j][i] = gT;
//                             }
//                         }
//                         else
//                         {
//                             oversetContravariantBC(mesh, i, j, k, gucart, aCell[b].face);
//                         }

//                     }

//                 }
//             }

//             sum_ind1 +=NumAcellPerProc[n];

//         }

//     }

//     std::vector<Acell> ().swap(aCell);
//     std::vector<std::vector<Dcell>> ().swap(dCell);
//     std::vector<std::vector<PetscReal>> ().swap(DWeights);
//     std::vector<std::vector<PetscInt>> ().swap(AcellProcMat);
//     std::vector<PetscInt> ().swap(NumAcellPerProc);

//     DMDAVecRestoreArray(fda0, ueqnP->lUcat, &lucat0);
//     DMDAVecRestoreArray(fda1, ueqn->Ucat, &ucat1);

//     if (flags->isTeqnActive)
//     {
//         DMDAVecRestoreArray(da0, teqnP->lTmprt, &ltemp0);
//         DMDAVecRestoreArray(da1, teqn->Tmprt, &temp1);

//         DMGlobalToLocalBegin(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
//         DMGlobalToLocalEnd(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
//     }

//     DMGlobalToLocalBegin(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
//     DMGlobalToLocalEnd(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

//     DMGlobalToLocalBegin(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
//     DMGlobalToLocalEnd(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

//     return 0;
// }

// //***************************************************************************************************************//
// PetscErrorCode acell1Dcell0Connectivity(mesh_ *meshP, mesh_ *mesh)
// {

//     overset_         *os  = mesh->access->os;
//     DM               da   = meshP->da, fda = meshP->fda;
//     DMDALocalInfo    info = meshP->info;
//     PetscInt         xs   = info.xs, xe = info.xs + info.xm;
//     PetscInt         ys   = info.ys, ye = info.ys + info.ym;
//     PetscInt         zs   = info.zs, ze = info.zs + info.zm;
//     PetscInt         mx   = info.mx, my = info.my, mz = info.mz;

//     PetscInt         lxs, lxe, lys, lye, lzs, lze;
//     PetscInt         i, j, k, b, m;

//     PetscMPIInt      rankP, sizeP, rank, size;
//     PetscReal        dist = 0., ds;

//     Cmpnts           ***cent;

//     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
//     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
//     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

//     MPI_Comm_size(meshP->MESH_COMM, &sizeP);
//     MPI_Comm_rank(meshP->MESH_COMM, &rankP);

//     MPI_Comm_size(mesh->MESH_COMM, &size);
//     MPI_Comm_rank(mesh->MESH_COMM, &rank);

//     DMDAVecGetArray(fda, meshP->lCent, &cent);

//     std::vector<Acell> aCell = os->aCell;

//     // vector that stores the number of aCell cells in each processor
//     std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;

//     std::vector<Dcell> dCellVec;
//     Dcell              dCell;

//     std::vector<std::vector<PetscInt>> lAcellProcMat(size);
//     os->AcellProcMat.resize(size);

//     os->oset_comm.resize(size);

//     for(b = 0; b < size; b++){
//         lAcellProcMat[b].resize(sizeP);
//         os->AcellProcMat[b].resize(sizeP);
//     }

//     os->dCell.resize( aCell.size() );   // stores the dCell elements within a given radius of each element of aCell vector

//     for(b = 0; b < aCell.size(); b++){

//         dist = aCell[b].cell_size * os->cellFactor;

//         for (k=lzs; k<lze; k++)
//             for (j=lys; j<lye; j++)
//                 for (i=lxs; i<lxe; i++){

//                     ds = sqrt(pow(cent[k][j][i].x - aCell[b].coorx,2.)
//                             +pow(cent[k][j][i].y - aCell[b].coory,2.)
//                             +pow(cent[k][j][i].z - aCell[b].coorz,2.));

//                     if (ds < dist){

//                         lAcellProcMat[aCell[b].rank][rankP] = 1;

//                         dCell.indi = i;
//                         dCell.indj = j;
//                         dCell.indk = k;

//                         dCell.rank = rank;
//                         dCell.dist2p = ds;

//                         dCellVec.push_back(dCell);

//                     }

//                 }

//         // store the connectivity for all the acceptor cells
//         os->dCell[b] = dCellVec;
//         std::vector<Dcell> ().swap(dCellVec);

//     }

//     // reduce the local processor matrix to get the global processor matrix
//     // this is a matrix of dim - size x size - rows - processors of aCell cells, columns- processors of dCell cells
//     for(b = 0; b < size; b++){
//         MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMat[b][0], sizeP, MPIU_INT, MPI_SUM, meshP->MESH_COMM);
//     }

//     // ensure the communicator colours are set right - include the diagonal element and set all that are 0 - MPI_UNDEFINED
//     // communicators will be created for each row of the processor matrix

//     if(size != sizeP)
//     {
//       char error[512];
//       sprintf(error, "current implementation requires the 2 meshes to have same number of processors. Recheck acell1Dcell0Connectivity function\n");
//       fatalErrorInFunction("acell1Dcell0Connectivity", error);
//     }

//     for(b = 0; b < size; b++){
//         if(NumAcellPerProc[b] != 0)
//             os->AcellProcMat[b][b] = 1; // set to 1  if the aCell cells are within the given processor
//         for(m = 0; m < sizeP; m++){
//             if (os->AcellProcMat[b][m] == 0)
//                 os->AcellProcMat[b][m] = MPI_UNDEFINED;
//         }
//     }

//     //create communicator for each row of the AcellProcMat that are non 0
//     for(b = 0; b < size; b++)
//     {
//       if(NumAcellPerProc[b]!=0)
//       {
//         MPI_Comm_split(PETSC_COMM_WORLD, os->AcellProcMat[b][rankP], rankP, &(os->oset_comm[b]));
//       }
//     }

//     std::vector<Acell> ().swap(aCell);

//     // vector that stores the number of aCell cells in each processor
//     std::vector<PetscInt> ().swap(NumAcellPerProc);

//     std::vector<std::vector<PetscInt>> ().swap(lAcellProcMat);

//     DMDAVecRestoreArray(fda, meshP->lCent, &cent);

//     return 0;
// }

// //***************************************************************************************************************//

// PetscErrorCode getLSWeights(mesh_ *meshP, mesh_ *mesh){

//     overset_         *os  = mesh->access->os;

//     DM               da1   = mesh->da, fda1 = mesh->fda;
//     DMDALocalInfo    info1 = mesh->info;
//     DM               da0   = meshP->da, fda0 = meshP->fda;
//     DMDALocalInfo    info0 = meshP->info;

//     PetscInt         i, j, k, ii, jj, kk, b, n, m;
//     PetscInt         nsupport = 0;
//     PetscReal        *W, *PHI, **P, **B, **A, **inv_A;
//     PetscReal        lA1[16] = {0.}, A1[16] = {0.};

//     Cmpnts           ***cent1, ***cent0;
//     PetscReal        ***nvert, ds, dist;
//     cellIds          supportNode;

//     std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
//     std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

//     os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

//     DMDAVecGetArray(fda1, mesh->lCent, &cent1);
//     DMDAVecGetArray(fda0, meshP->lCent, &cent0);

//     for(b = 0; b < aCell.size(); b++){

//         os->DWeights[b].resize(dCell[b].size());

//         // aCell cell index
//         i = aCell[b].indi;
//         j = aCell[b].indj;
//         k = aCell[b].indk;

//         // support radius size
//         ds = aCell[b].cell_size * os->cellFactor;

//         nsupport = dCell[b].size();

//         // allocate local variables for MLS interpolation
//         B = (PetscReal**) malloc(4*sizeof(PetscReal*));

//         for (n=0;n<4;n++)
//         {
//             B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
//         }

//         P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
//         for (n=0;n<nsupport;n++)
//         {
//             P[n] = (PetscReal*) malloc(4*sizeof(PetscReal));
//         }

//         W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         PetscMalloc(sizeof(PetscReal *)  * (4), &(A));
//         PetscMalloc(sizeof(PetscReal *)  * (4), &(inv_A));

//         for(n=0; n<4; n++)
//         {
//             PetscMalloc(sizeof(PetscReal)  * (4), &(A[n]));
//             PetscMalloc(sizeof(PetscReal)  * (4), &(inv_A[n]));
//         }

//         //initialise matrix and vectors
//         for (ii=0; ii<4; ii++)
//         {
//             for (jj=0; jj<4; jj++)
//             {
//                 A[ii][jj] = 0.;
//             }
//         }
//         for (ii=0; ii<16; ii++){
//             A1[ii] = 0.;
//             lA1[ii] = 0.;
//         }

//         for (n=0; n<nsupport; n++){
//             PHI[n] = 0.;
//             os->DWeights[b][n] = 0.;
//         }

//         // loop through the donor cells within the current processor - gnumdCell[b][rank]
//         for(m = 0; m < nsupport; m++){

//             kk = dCell[b][m].indk;
//             jj = dCell[b][m].indj;
//             ii = dCell[b][m].indi;

//             P[m][0] = 1.0;
//             P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
//             P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
//             P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;

//             // get normalized distance
//             dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
//                     +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
//                     +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

//             // get interpolation weights
//             W[m]
//               =
//                       (dist < 0.5)
//                       ?
//                               2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
//             :
//                               4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

//             // set B
//             for (ii=0; ii<4; ii++){
//                 B[ii][m] = W[m] * P[m][ii];
//             }


//             // set A
//             for (ii=0; ii<4; ii++)
//             {
//                 for (jj=0; jj<4; jj++)
//                 {

//                     A[ii][jj] += B[ii][m] * P[m][jj];
//                 }
//             }

//         }

//         // set A1 elements to scatter A matrix elements
//         for (ii=0; ii<4; ii++)
//         {
//             for (jj=0; jj<4; jj++)
//             {
//                 lA1[4*ii + jj] = A[ii][jj];
//             }
//         }

//         // reduce and scatter A matrix to all processors
//         MPI_Allreduce(&lA1, &A1, 16, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

//         for (ii=0; ii<4; ii++)
//         {
//             for (jj=0; jj<4; jj++)
//             {
//                 A[ii][jj] = A1[4*ii + jj];
//             }
//         }

//         //invert matrix A
//         inv_4by4(A, inv_A, 4);

//         for(m = 0; m < nsupport; m++){
//             for (ii=0; ii<4; ii++)
//             {
//                 PHI[m] += inv_A[0][ii]*B[ii][m];
//             }

//             os->DWeights[b][m] = PHI[m];
//         }


//         // free all the vectors
//         for(PetscInt ind=0; ind<4; ind++)
//         {
//             free(B[ind]);
//             free(A[ind]);
//             free(inv_A[ind]);
//         }

//         for(PetscInt ind=0; ind<nsupport; ind++)
//         {
//             free(P[ind]);
//         }

//         free(B);
//         free(A);
//         free(inv_A);
//         free(P);
//         free(W);
//         free(PHI);
//     }


//     DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
//     DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
//     return 0;
// }

// //***************************************************************************************************************//

// PetscErrorCode getLSWeights_2nd(mesh_ *meshP, mesh_ *mesh)
// {
//     overset_         *os  = mesh->access->os;

//     DM               da1   = mesh->da, fda1 = mesh->fda;
//     DMDALocalInfo    info1 = mesh->info;
//     DM               da0   = meshP->da, fda0 = meshP->fda;
//     DMDALocalInfo    info0 = meshP->info;

//     PetscInt         i, j, k, ii, jj, kk, b, n, m;
//     PetscInt         nsupport = 0;
//     PetscReal        *W, *PHI, **P, **B;
//     PetscReal        A[10][10]={0.}, inv_A[10][10]={0.}, lA1[100] = {0.}, A1[100] = {0.};

//     Cmpnts           ***cent1, ***cent0;
//     PetscReal        ***nvert, ds, dist;
//     cellIds          supportNode;

//     std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
//     std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

//     os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

//     DMDAVecGetArray(fda1, mesh->lCent, &cent1);
//     DMDAVecGetArray(fda0, meshP->lCent, &cent0);

//     for(b = 0; b < aCell.size(); b++){

//         os->DWeights[b].resize(dCell[b].size());

//         // aCell cell index
//         i = aCell[b].indi;
//         j = aCell[b].indj;
//         k = aCell[b].indk;

//         // support radius size
//         ds = aCell[b].cell_size * os->cellFactor;

//         nsupport = dCell[b].size();

//         // allocate local variables for MLS interpolation
//         B = (PetscReal**) malloc(10*sizeof(PetscReal*));
//         for (n=0;n<10;n++)
//         {
//             B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
//         }

//         P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
//         for (n=0;n<nsupport;n++)
//         {
//             P[n] = (PetscReal*) malloc(10*sizeof(PetscReal));
//         }

//         W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         //initialise matrix and vectors
//         for (ii=0; ii<10; ii++)
//         {
//             for (jj=0; jj<10; jj++)
//             {
//                 A[ii][jj] = 0.;
//             }
//         }
//         for (ii=0; ii<100; ii++){
//             A1[ii] = 0.;
//             lA1[ii] = 0.;
//         }

//         for (n=0; n<nsupport; n++){
//             PHI[n] = 0.;
//             os->DWeights[b][n] = 0.;
//         }

//         // loop through the donor cells within the current processor - gnumdCell[b][rank]
//         for(m = 0; m < nsupport; m++){

//             kk = dCell[b][m].indk;
//             jj = dCell[b][m].indj;
//             ii = dCell[b][m].indi;

//             P[m][0] = 1.0;
//             P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
//             P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
//             P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;
//             P[m][4]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
//             P[m][5]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
//             P[m][6]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
//             P[m][7]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
//             P[m][8]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
//             P[m][9]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;

//             // get normalized distance
//             dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
//                     +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
//                     +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

//             // get interpolation weights
//             W[m]
//               =
//                       (dist < 0.5)
//                       ?
//                               2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
//             :
//                               4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

//             // set B
//             for (ii=0; ii<10; ii++){
//                 B[ii][m] = W[m] * P[m][ii];
//             }


//             // set A
//             for (ii=0; ii<10; ii++)
//             {
//                 for (jj=0; jj<10; jj++)
//                 {

//                     A[ii][jj] += B[ii][m] * P[m][jj];
//                 }
//             }

//         }

//         // set A1 elements to scatter A matrix elements
//         for (ii=0; ii<10; ii++)
//         {
//             for (jj=0; jj<10; jj++)
//             {
//                 lA1[10*ii + jj] = A[ii][jj];
//             }
//         }

//         // reduce and scatter A matrix to all processors
//         MPI_Allreduce(&lA1, &A1, 100, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

//         for (ii=0; ii<10; ii++)
//         {
//             for (jj=0; jj<10; jj++)
//             {
//                 A[ii][jj] = A1[10*ii + jj];
//             }
//         }

//         //invert matrix A
//         inv_10(A,inv_A,10);

//         for(m = 0; m < nsupport; m++){
//             for (ii=0; ii<10; ii++)
//             {
//                 PHI[m] += inv_A[0][ii]*B[ii][m];
//             }

//             os->DWeights[b][m] = PHI[m];
//         }


//         // free all the vectors
//         for(PetscInt ind=0; ind<10; ind++)
//         {
//             free(B[ind]);
//         }

//         for(PetscInt ind=0; ind<nsupport; ind++)
//         {
//             free(P[ind]);
//         }

//         free(B);
//         free(P);
//         free(W);
//         free(PHI);
//     }


//     DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
//     DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
//     return 0;
// }

// //***************************************************************************************************************//
// // third order least square method for finding the weights for the donor cells of the background mesh to interpolate the overset mesh acceptor cells
// PetscErrorCode getLSWeights_3rd(mesh_ *meshP, mesh_ *mesh)
// {
//     overset_         *os  = mesh->access->os;

//     DM               da1   = mesh->da, fda1 = mesh->fda;
//     DMDALocalInfo    info1 = mesh->info;
//     DM               da0   = meshP->da, fda0 = meshP->fda;
//     DMDALocalInfo    info0 = meshP->info;

//     PetscInt         ii, jj, kk, b, n, m;
//     PetscInt         nsupport = 0;
//     PetscReal        *W, *PHI, **P, **B;
//     PetscReal        A[20][20]={0.}, inv_A[20][20]={0.}, lA1[400] = {0.}, A1[400] = {0.};

//     Cmpnts           ***cent1, ***cent0;
//     PetscReal        ***nvert, ds, dist;
//     cellIds          supportNode;

//     std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
//     std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

//     os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

//     DMDAVecGetArray(fda1, mesh->lCent, &cent1);
//     DMDAVecGetArray(fda0, meshP->lCent, &cent0);

//     for(b = 0; b < aCell.size(); b++){

//         os->DWeights[b].resize(dCell[b].size());

//         // support radius size
//         ds = aCell[b].cell_size * os->cellFactor;

//         nsupport = dCell[b].size();

//         // allocate local variables for MLS interpolation
//         B = (PetscReal**) malloc(20*sizeof(PetscReal*));
//         for (n=0;n<20;n++)
//         {
//             B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
//         }

//         P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
//         for (n=0;n<nsupport;n++)
//         {
//             P[n] = (PetscReal*) malloc(20*sizeof(PetscReal));
//         }

//         W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

//         //initialise matrix and vectors
//         for (ii=0; ii<20; ii++)
//         {
//             for (jj=0; jj<20; jj++)
//             {
//                 A[ii][jj] = 0.;
//             }
//         }
//         for (ii=0; ii<400; ii++){
//             A1[ii] = 0.;
//             lA1[ii] = 0.;
//         }

//         for (n=0; n<nsupport; n++){
//             PHI[n] = 0.;
//             os->DWeights[b][n] = 0.;
//         }

//         // loop through the donor cells within the current processor - gnumdCell[b][rank]
//         for(m = 0; m < nsupport; m++){

//             kk = dCell[b][m].indk;
//             jj = dCell[b][m].indj;
//             ii = dCell[b][m].indi;

//             P[m][0] = 1.0;
//             P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
//             P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
//             P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;
//             P[m][4]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
//             P[m][5]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
//             P[m][6]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
//             P[m][7]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
//             P[m][8]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
//             P[m][9]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
//             P[m][10]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
//             P[m][11]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
//             P[m][12]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
//             P[m][13]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
//             P[m][14]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
//             P[m][15]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
//             P[m][16]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
//             P[m][17]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
//             P[m][18]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
//             P[m][19]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;


//             // get normalized distance
//             dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
//                     +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
//                     +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

//             // get interpolation weights
//             W[m]
//               =
//                       (dist < 0.5)
//                       ?
//                               2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
//             :
//                               4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

//             // set B
//             for (ii=0; ii<20; ii++){
//                 B[ii][m] = W[m] * P[m][ii];
//             }


//             // set A
//             for (ii=0; ii<20; ii++)
//             {
//                 for (jj=0; jj<20; jj++)
//                 {

//                     A[ii][jj] += B[ii][m] * P[m][jj];
//                 }
//             }

//         }

//         // set A1 elements to scatter A matrix elements
//         for (ii=0; ii<20; ii++)
//         {
//             for (jj=0; jj<20; jj++)
//             {
//                 lA1[20*ii + jj] = A[ii][jj];
//             }
//         }

//         // reduce and scatter A matrix to all processors
//         MPI_Allreduce(&lA1, &A1, 400, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);

//         for (ii=0; ii<20; ii++)
//         {
//             for (jj=0; jj<20; jj++)
//             {
//                 A[ii][jj] = A1[20*ii + jj];
//             }
//         }

//         //invert matrix A
//         inv_20(A,inv_A,20);

//         for(m = 0; m < nsupport; m++){
//             for (ii=0; ii<20; ii++)
//             {
//                 PHI[m] += inv_A[0][ii]*B[ii][m];
//             }

//             os->DWeights[b][m] = PHI[m];
//         }


//         // free all the vectors
//         for(PetscInt ind=0; ind<20; ind++)
//         {
//             free(B[ind]);
//         }

//         for(PetscInt ind=0; ind<nsupport; ind++)
//         {
//             free(P[ind]);
//         }

//         free(B);
//         free(P);
//         free(W);
//         free(PHI);
//     }


//     DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
//     DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
//     return 0;
// }

// //***************************************************************************************************************//

// PetscErrorCode oversetContravariantBC(mesh_ *mesh, PetscInt i, PetscInt j, PetscInt k, Cmpnts ucart, PetscInt face)
// {

//     ueqn_         *ueqn = mesh->access->ueqn;
//     DM            da    = mesh->da, fda = mesh->fda;
//     DMDALocalInfo info  = mesh->info;
//     PetscInt      xs    = info.xs, xe = info.xs + info.xm;
//     PetscInt      ys    = info.ys, ye = info.ys + info.ym;
//     PetscInt      zs    = info.zs, ze = info.zs + info.zm;
//     PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

//     Cmpnts           ***ucont, ***icsi, ***jeta, ***kzet;

//     DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
//     DMDAVecGetArray(fda, mesh->lICsi, &icsi);
//     DMDAVecGetArray(fda, mesh->lJEta, &jeta);
//     DMDAVecGetArray(fda, mesh->lKZet, &kzet);

//     // flux bc based on the face where ucart is interpolated
//     if(face == 1 && !(mesh->i_periodic) && !(mesh->ii_periodic))
//     {
//         ucont[k][j][i-1].x = (ucart.x * icsi[k][j][i-1].x + ucart.y * icsi[k][j][i-1].y + ucart.z * icsi[k][j][i-1].z);
//     }

//     if(face == 2 && !(mesh->i_periodic) && !(mesh->ii_periodic))
//     {
//         ucont[k][j][i].x = (ucart.x * icsi[k][j][i].x + ucart.y * icsi[k][j][i].y + ucart.z * icsi[k][j][i].z);
//     }

//     if(face == 3 && !(mesh->j_periodic) && !(mesh->jj_periodic))
//     {
//         ucont[k][j-1][i].y = (ucart.x * jeta[k][j-1][i].x + ucart.y * jeta[k][j-1][i].y + ucart.z * jeta[k][j-1][i].z);
//     }

//     if(face == 4 && !(mesh->j_periodic) && !(mesh->jj_periodic))
//     {
//         ucont[k][j][i].y = (ucart.x * jeta[k][j][i].x + ucart.y * jeta[k][j][i].y + ucart.z * jeta[k][j][i].z);
//     }

//     if(face == 5 && !(mesh->k_periodic) && !(mesh->kk_periodic))
//     {
//         ucont[k-1][j][i].z = (ucart.x * kzet[k-1][j][i].x + ucart.y * kzet[k-1][j][i].y + ucart.z * kzet[k-1][j][i].z );
//     }

//     if(face == 6 && !(mesh->k_periodic) && !(mesh->kk_periodic))
//     {
//         ucont[k][j][i].z = (ucart.x * kzet[k][j][i].x + ucart.y * kzet[k][j][i].y + ucart.z * kzet[k][j][i].z );
//     }

//     DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
//     DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
//     DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
//     DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

//     return 0;
// }