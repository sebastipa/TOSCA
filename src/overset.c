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

//***************************************************************************************************************//

PetscErrorCode InitializeOverset(domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;

    // Read all hole objects from the input file
    std::vector<HoleObject> holeObjects;
    readHoleObjects(holeObjects, domain[0].info.nHoleRegions);

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
            PetscPrintf(domain[d].mesh->MESH_COMM, "Started recursive acceptor search from domain %ld:\n", d);
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

        // poisson equation initialize
        InitializePEqn(domain[d].peqn);
    }

    // Re-run a clean interpolation cycle with all domains fully initialized
    PetscPrintf(PETSC_COMM_WORLD, "\nOverset: initializing interface cells after field initialization...\n");
    UpdateOversetInterpolation(domain);

    // Save old-field vectors so that Ucont_o and Tmprt_o reflect the corrected interface state after interpolation  
    PetscInt nDomainsFinal = domain[0].info.nDomains;
    for (PetscInt d = 0; d < nDomainsFinal; d++)
    {
        VecCopy(domain[d].ueqn->Ucont, domain[d].ueqn->Ucont_o);

        if (domain[d].flags.isTeqnActive)
        {
            VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
        }
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode UpdateOversetInterpolation(domain_ *domain)
{
    PetscPrintf(PETSC_COMM_WORLD, "\n");

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
    PetscPrintf(PETSC_COMM_WORLD, "OS interpolation, Elapsed Time = %lf\n", te - ts);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode UpdateDomainInterpolation(PetscInt d, domain_ *domain, PetscInt level)
{
    // function to update interpolation for a single domain and its dependencies recursively

    PetscInt nDomains = domain[0].info.nDomains;
    if (d < 0 || d >= nDomains) return 0; 
    if (domain[d].os == NULL) return 0;   

    overset_ *os    = domain[d].os;
    mesh_    *mesh  = domain[d].mesh;

    // 1: update interpolation from parent meshes to this domain: updates interpolated domain boundaries
    for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
    {
        if (os->parentMeshId[pi] != -1)
        {
            mesh_ *parentMesh = domain[os->parentMeshId[pi]].mesh;

            // print interpolation info for debugging
            PetscPrintf(mesh->MESH_COMM, "OS Interpolation, %s > %s (acceptor is lvl %ld)\n", parentMesh->meshName.c_str(), mesh->meshName.c_str(), level);

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

    // 2: update boundary conditions: updates remaining domain boundaries
    UpdateCartesianBCs(domain[d].ueqn);
    UpdateContravariantBCs(domain[d].ueqn);
    UpdatePressureBCs(domain[d].peqn);

    if (domain[d].flags.isTeqnActive)
    {
        UpdateTemperatureBCs(domain[d].teqn);
    }

    // 3: update interpolation from child meshes to this domain: updates hole boundaries
    for (PetscInt ci = 0; ci < os->childMeshId.size(); ci++)
    {
        if (os->childMeshId[ci] != -1)
        {
            mesh_ *childMesh = domain[os->childMeshId[ci]].mesh;
            PetscInt childId = os->childMeshId[ci];

            // print interpolation info for debugging
            PetscPrintf(mesh->MESH_COMM, "OS Interpolation, %s > %s (acceptor is lvl %ld)\n", childMesh->meshName.c_str(),  mesh->meshName.c_str(), level);

            interpolateACellTrilinearC2P(childMesh, mesh, childId);

            MPI_Barrier(mesh->MESH_COMM);

            // Recursively update child domain
            UpdateDomainInterpolation(childId, domain, level + 1);
        }
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode SyncPressureAcrossDomains(domain_ *domain)
{
    PetscInt nDomains = domain[0].info.nDomains;

    // Single domain: shift entire pressure field so that p[1][1][1] = 0
    if (nDomains == 1)
    {
        mesh_         *mesh = domain[0].mesh;
        DM             da   = mesh->da;
        DMDALocalInfo  info = mesh->info;
        PetscInt       xs   = info.xs, xe = info.xs + info.xm;
        PetscInt       ys   = info.ys, ye = info.ys + info.ym;
        PetscInt       zs   = info.zs, ze = info.zs + info.zm;

        PetscReal ***p;
        DMDAVecGetArray(da, domain[0].peqn->P, &p);

        // read p[1][1][1] on whichever rank owns it and send it to all other ranks
        PetscReal localRef = 0.0;

        if (xs <= 1 && 1 < xe && ys <= 1 && 1 < ye && zs <= 1 && 1 < ze)
        {
            localRef  = p[1][1][1];
        }

        PetscReal globalRef = 0.0;
        MPI_Allreduce(&localRef, &globalRef, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        // shift all pressure values 
        for (PetscInt k = zs; k < ze; k++)
            for (PetscInt j = ys; j < ye; j++)
                for (PetscInt i = xs; i < xe; i++)
                    p[k][j][i] -= globalRef;

        DMDAVecRestoreArray(da, domain[0].peqn->P, &p);

        DMGlobalToLocalBegin(da, domain[0].peqn->P, INSERT_VALUES, domain[0].peqn->lP);
        DMGlobalToLocalEnd  (da, domain[0].peqn->P, INSERT_VALUES, domain[0].peqn->lP);

        PetscPrintf(mesh->MESH_COMM, "Gauge pressure removal, shift is %.6e\n", -globalRef);

        return 0;
    }

    // Multiple domains: first find the root parent, then for each non root shift its
    // pressure to match the parent's value at the cell closest to the child's cent[1][1][1]. 
    // This ensures that the gauge is consistent across all domains and that the solution 
    // remains attached across the overset interfaces (it doesn't matter currently, but we do not
    // want pressure shifts). 

    // Identify the root parent domain
    PetscInt rootId = -1;
    for (PetscInt d = 0; d < nDomains; d++)
    {
        PetscBool isChild = PETSC_FALSE;
        for (PetscInt other = 0; other < nDomains && !isChild; other++)
        {
            if (other == d || domain[other].os == NULL) continue;
            for (PetscInt ci = 0; ci < (PetscInt)domain[other].os->childMeshId.size(); ci++)
            {
                if (domain[other].os->childMeshId[ci] == d)
                {
                    isChild = PETSC_TRUE;
                    break;
                }
            }
        }
        if (!isChild) { rootId = d; break; }
    }

    if (rootId < 0)
    {
        char error[512];
        sprintf(error, "no root parent domain found for pressure synchronization\n");
        fatalErrorInFunction("SyncPressureAcrossDomains", error);
    }

    // Shift pressure so that p[1][1][1] = 0 in the root parent domain
    {
        mesh_        *rootMesh = domain[rootId].mesh;
        DM            rootDa   = rootMesh->da;
        DMDALocalInfo rInfo    = rootMesh->info;
        PetscInt      rxs = rInfo.xs, rxe = rInfo.xs + rInfo.xm;
        PetscInt      rys = rInfo.ys, rye = rInfo.ys + rInfo.ym;
        PetscInt      rzs = rInfo.zs, rze = rInfo.zs + rInfo.zm;

        PetscReal ***rootP;
        DMDAVecGetArray(rootDa, domain[rootId].peqn->P, &rootP);

        PetscReal localRef = 0.0;
        if (rxs <= 1 && 1 < rxe && rys <= 1 && 1 < rye && rzs <= 1 && 1 < rze)
            localRef = rootP[1][1][1];

        PetscReal globalRef = 0.0;
        MPI_Allreduce(&localRef, &globalRef, 1, MPIU_REAL, MPIU_SUM, rootMesh->MESH_COMM);

        for (PetscInt k = rzs; k < rze; k++)
            for (PetscInt j = rys; j < rye; j++)
                for (PetscInt i = rxs; i < rxe; i++)
                    rootP[k][j][i] -= globalRef;

        DMDAVecRestoreArray(rootDa, domain[rootId].peqn->P, &rootP);

        // sync lP so children read the shifted values below
        DMGlobalToLocalBegin(rootDa, domain[rootId].peqn->P, INSERT_VALUES, domain[rootId].peqn->lP);
        DMGlobalToLocalEnd  (rootDa, domain[rootId].peqn->P, INSERT_VALUES, domain[rootId].peqn->lP);

        PetscPrintf(rootMesh->MESH_COMM, "OS gauge pressure removal: root domain %ld shift is %.6e\n", rootId, -globalRef);
    }

    // Shift each child domain so that its closest cell to the parent's cent[1][1][1] has the same pressure as the parent at that cell
    // Note: use multiple passes to handle any domain numbering order: after at most
    //       nDomains-1 passes every level in a telescopic chain is processed.
    //       A domain is processed only once its direct parent's lP is already finalized.
    std::vector<PetscBool> done(nDomains, PETSC_FALSE);
    done[rootId] = PETSC_TRUE;

    for (PetscInt pass = 0; pass < nDomains; pass++)
    {
        for (PetscInt d = 0; d < nDomains; d++)
        {
            if (done[d]) continue;

            // skip domains with no overset structure or no parent
            if (domain[d].os == NULL) { done[d] = PETSC_TRUE; continue; }
            if (domain[d].os->parentMeshId.empty()) { done[d] = PETSC_TRUE; continue; }

            PetscInt parentId = domain[d].os->parentMeshId[0];
            if (parentId < 0 || parentId >= nDomains) { done[d] = PETSC_TRUE; continue; }

            // only process this domain once its direct parent is finalized
            if (!done[parentId]) continue;

            mesh_        *childMesh  = domain[d].mesh;
            mesh_        *parentMesh = domain[parentId].mesh;
            DM            childDa    = childMesh->da;
            DM            parentDa   = parentMesh->da;
            DM            parentFda  = parentMesh->fda;

            DMDALocalInfo  cInfo  = childMesh->info;
            DMDALocalInfo  pInfo  = parentMesh->info;

            PetscReal ***childP, ***parentP;
            Cmpnts    ***parentCent;

            // Step 1: broadcast cent[1][1][1] of the child 
            PetscReal localChildX = 0.0, localChildY = 0.0, localChildZ = 0.0;
            {
                DM     cFda = childMesh->fda;
                Cmpnts ***childCent;

                DMDAVecGetArray(cFda, childMesh->lCent, &childCent);

                PetscInt cxs = cInfo.xs, cxe = cInfo.xs + cInfo.xm;
                PetscInt cys = cInfo.ys, cye = cInfo.ys + cInfo.ym;
                PetscInt czs = cInfo.zs, cze = cInfo.zs + cInfo.zm;

                if (cxs <= 1 && 1 < cxe && cys <= 1 && 1 < cye && czs <= 1 && 1 < cze)
                {
                    localChildX = childCent[1][1][1].x;
                    localChildY = childCent[1][1][1].y;
                    localChildZ = childCent[1][1][1].z;
                }

                DMDAVecRestoreArray(cFda, childMesh->lCent, &childCent);
            }

            PetscReal childRefX = 0.0, childRefY = 0.0, childRefZ = 0.0;
            MPI_Allreduce(&localChildX, &childRefX, 1, MPIU_REAL, MPIU_SUM, childMesh->MESH_COMM);
            MPI_Allreduce(&localChildY, &childRefY, 1, MPIU_REAL, MPIU_SUM, childMesh->MESH_COMM);
            MPI_Allreduce(&localChildZ, &childRefZ, 1, MPIU_REAL, MPIU_SUM, childMesh->MESH_COMM);

            // Step 2: find closest cell in parent and read its p value
            DMDAVecGetArray(parentFda, parentMesh->lCent, &parentCent);
            DMDAVecGetArray(parentDa,  domain[parentId].peqn->lP, &parentP);

            PetscInt  pxs = pInfo.xs, pxe = pInfo.xs + pInfo.xm;
            PetscInt  pys = pInfo.ys, pye = pInfo.ys + pInfo.ym;
            PetscInt  pzs = pInfo.zs, pze = pInfo.zs + pInfo.zm;
            PetscInt  pmx = pInfo.mx, pmy = pInfo.my, pmz = pInfo.mz;

            PetscReal localMinDist = 1.0e30;
            PetscReal localParentP = 0.0;

            for (PetscInt k = PetscMax(pzs,1); k < PetscMin(pze, pmz-1); k++)
            {
                for (PetscInt j = PetscMax(pys,1); j < PetscMin(pye, pmy-1); j++)
                {
                    for (PetscInt i = PetscMax(pxs,1); i < PetscMin(pxe, pmx-1); i++)
                    {
                        PetscReal dx = parentCent[k][j][i].x - childRefX;
                        PetscReal dy = parentCent[k][j][i].y - childRefY;
                        PetscReal dz = parentCent[k][j][i].z - childRefZ;
                        PetscReal dist = dx*dx + dy*dy + dz*dz;
                        if (dist < localMinDist)
                        {
                            localMinDist = dist;
                            localParentP = parentP[k][j][i];
                        }
                    }
                }
            }

            DMDAVecRestoreArray(parentFda, parentMesh->lCent, &parentCent);
            DMDAVecRestoreArray(parentDa,  domain[parentId].peqn->lP, &parentP);

            // find global minimum and corresponding parent p value
            PetscReal globalMinDist = 0.0;
            MPI_Allreduce(&localMinDist, &globalMinDist, 1, MPIU_REAL, MPIU_MIN, parentMesh->MESH_COMM);

            PetscReal localContrib = 0.0;
            if (localMinDist == globalMinDist)
                localContrib = localParentP;

            PetscReal parentRefP = 0.0;
            MPI_Allreduce(&localContrib, &parentRefP, 1, MPIU_REAL, MPIU_SUM, parentMesh->MESH_COMM);

            // Step 3: read child p[1][1][1] 
            DMDAVecGetArray(childDa, domain[d].peqn->P, &childP);

            PetscReal localChildP = 0.0;
            {
                PetscInt cxs = cInfo.xs, cxe = cInfo.xs + cInfo.xm;
                PetscInt cys = cInfo.ys, cye = cInfo.ys + cInfo.ym;
                PetscInt czs = cInfo.zs, cze = cInfo.zs + cInfo.zm;
                if (cxs <= 1 && 1 < cxe && cys <= 1 && 1 < cye && czs <= 1 && 1 < cze)
                    localChildP = childP[1][1][1];
            }
            PetscReal globalChildP = 0.0;
            MPI_Allreduce(&localChildP, &globalChildP, 1, MPIU_REAL, MPIU_SUM, childMesh->MESH_COMM);

            // Step 4: shift child pressure by the difference 
            PetscReal shift = parentRefP - globalChildP;

            PetscInt cxs = cInfo.xs, cxe = cInfo.xs + cInfo.xm;
            PetscInt cys = cInfo.ys, cye = cInfo.ys + cInfo.ym;
            PetscInt czs = cInfo.zs, cze = cInfo.zs + cInfo.zm;

            for (PetscInt k = czs; k < cze; k++)
                for (PetscInt j = cys; j < cye; j++)
                    for (PetscInt i = cxs; i < cxe; i++)
                        childP[k][j][i] += shift;

            DMDAVecRestoreArray(childDa, domain[d].peqn->P, &childP);

            DMGlobalToLocalBegin(childDa, domain[d].peqn->P, INSERT_VALUES, domain[d].peqn->lP);
            DMGlobalToLocalEnd  (childDa, domain[d].peqn->P, INSERT_VALUES, domain[d].peqn->lP);

            PetscPrintf(childMesh->MESH_COMM,
                "OS gauge pressure removal: child (%ld) to parent (%ld) shift is %+.6e\n",
                d, parentId, shift);

            done[d] = PETSC_TRUE;
        }
    }

    PetscPrintf(PETSC_COMM_WORLD, "\n");

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

    SetInitialField(&domain[d]);

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

PetscErrorCode FindHoleObject(const std::vector<HoleObject> &holeObjects, PetscInt parentId, PetscInt childId, char **holeObjectName)
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

//***************************************************************************************************************//

PetscErrorCode interpolateACellTrilinearP2C(mesh_ *meshD, mesh_ *meshA)
{
    overset_         *os     = meshA->access->os;
    ueqn_            *ueqnA  = meshA->access->ueqn;
    ueqn_            *ueqnD  = meshD->access->ueqn;
    teqn_            *teqnA  = meshA->access->teqn;
    teqn_            *teqnD  = meshD->access->teqn;
    peqn_            *peqnA  = meshA->access->peqn;
    peqn_            *peqnD  = meshD->access->peqn;
    flags_           *flags  = meshA->access->flags;
    DM               daA     = meshA->da, fdaA = meshA->fda;
    DM               daD     = meshD->da, fdaD = meshD->fda;

    Cmpnts           ***lucatD, ***ucatA, ***cent;
    PetscReal        ***ltempD = NULL, ***tempA = NULL;
    PetscReal        ***lpressD, ***pressA;

    DMDAVecGetArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecGetArray(fdaA, ueqnA->Ucat,  &ucatA);
    DMDAVecGetArray(fdaD, meshD->lCent, &cent);
    DMDAVecGetArray(daD,  peqnD->lP,    &lpressD);
    DMDAVecGetArray(daA,  peqnA->P,     &pressA);
    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecGetArray(daA, teqnA->Tmprt,  &tempA);
    }

    // Donor (root) side: pack interpolated Ux, Uy, Uz, T, p per slot
    const PetscInt stride = 5;
    PetscInt nRoots = os->nRootsDb;
    std::vector<PetscReal> rootBuf((size_t)nRoots * stride, 0.0);
    const std::vector<PetscInt>  &flatI = os->rootSlotIdxDb;
    const std::vector<PetscReal> &flatC = os->rootSlotCoordsDb;

    for (PetscInt s = 0; s < nRoots; s++)
    {
        PetscInt ic = flatI[3*s + 0];
        PetscInt jc = flatI[3*s + 1];
        PetscInt kc = flatI[3*s + 2];
        PetscReal px = flatC[3*s + 0];
        PetscReal py = flatC[3*s + 1];
        PetscReal pz = flatC[3*s + 2];

        Cmpnts ucart;
        vectorPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, lucatD, ucart);

        rootBuf[s*stride + 0] = ucart.x;
        rootBuf[s*stride + 1] = ucart.y;
        rootBuf[s*stride + 2] = ucart.z;

        if (flags->isTeqnActive)
        {
            PetscReal Temp = 0.0;
            scalarPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, ltempD, Temp);
            rootBuf[s*stride + 3] = Temp;
        }

        PetscReal Pres = 0.0;
        scalarPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, lpressD, Pres);
        rootBuf[s*stride + 4] = Pres;
    }

    // Acceptor (leaf) side: one slot per local acceptor
    PetscInt nLeaves = (PetscInt)os->localAcceptorsDb.size();
    std::vector<PetscReal> leafBuf((size_t)nLeaves * stride, 0.0);

    // PetscSFBcast as 5 contiguous doubles per slot via a derived type
    MPI_Datatype t5;
    MPI_Type_contiguous(5, MPIU_REAL, &t5);
    MPI_Type_commit(&t5);

    PetscSFBcastBegin(os->sfP2C, t5, rootBuf.data(), leafBuf.data(), MPI_REPLACE);
    PetscSFBcastEnd  (os->sfP2C, t5, rootBuf.data(), leafBuf.data(), MPI_REPLACE);

    MPI_Type_free(&t5);

    // Scatter the received slot data into each acceptor's (indi, indj, indk) cell
    for (PetscInt b = 0; b < nLeaves; b++)
    {
        if (os->localDonorMapDb[b].rank < 0) continue;
        const Acell &ac = os->localAcceptorsDb[b];
        PetscInt i = ac.indi, j = ac.indj, k = ac.indk;
        ucatA[k][j][i].x = leafBuf[b*stride + 0];
        ucatA[k][j][i].y = leafBuf[b*stride + 1];
        ucatA[k][j][i].z = leafBuf[b*stride + 2];
        pressA[k][j][i]  = leafBuf[b*stride + 4];
        if (flags->isTeqnActive)
        {
            tempA[k][j][i] = leafBuf[b*stride + 3];
        }
    }

    DMDAVecRestoreArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecRestoreArray(fdaA, ueqnA->Ucat,  &ucatA);
    DMDAVecRestoreArray(fdaD, meshD->lCent, &cent);
    DMDAVecRestoreArray(daD,  peqnD->lP,    &lpressD);
    DMDAVecRestoreArray(daA,  peqnA->P,     &pressA);
    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecRestoreArray(daA, teqnA->Tmprt,  &tempA);
        DMGlobalToLocalBegin(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
        DMGlobalToLocalEnd  (daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
    }
    DMGlobalToLocalBegin(daA,  peqnA->P,    INSERT_VALUES, peqnA->lP);
    DMGlobalToLocalEnd  (daA,  peqnA->P,    INSERT_VALUES, peqnA->lP);
    DMGlobalToLocalBegin(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);
    DMGlobalToLocalEnd  (fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);

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
    peqn_            *peqnA  = meshA->access->peqn;
    peqn_            *peqnD  = meshD->access->peqn;
    flags_           *flags  = meshA->access->flags;
    DM               daA     = meshA->da, fdaA = meshA->fda;
    DM               daD     = meshD->da, fdaD = meshD->fda;

    Cmpnts           ***lucatD, ***ucatA, ***cent;
    PetscReal        ***ltempD = NULL, ***tempA = NULL;
    PetscReal        ***lpressD, ***pressA;

    // Nothing to do if no acceptors are registered for this donor mesh
    if (os->localAcceptorsHc.find(donorId) == os->localAcceptorsHc.end() ||
        os->sfC2P.find(donorId) == os->sfC2P.end())
    {
        return 0;
    }

    DMDAVecGetArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecGetArray(fdaA, ueqnA->Ucat,  &ucatA);
    DMDAVecGetArray(fdaD, meshD->lCent, &cent);
    DMDAVecGetArray(daD,  peqnD->lP,    &lpressD);
    DMDAVecGetArray(daA,  peqnA->P,     &pressA);
    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecGetArray(daA, teqnA->Tmprt,  &tempA);
    }

    const PetscInt stride = 5;
    PetscInt nRoots = os->nRootsHc[donorId];
    std::vector<PetscReal> rootBuf((size_t)nRoots * stride, 0.0);
    const std::vector<PetscInt>  &flatI = os->rootSlotIdxHc[donorId];
    const std::vector<PetscReal> &flatC = os->rootSlotCoordsHc[donorId];

    for (PetscInt s = 0; s < nRoots; s++)
    {
        PetscInt ic = flatI[3*s + 0];
        PetscInt jc = flatI[3*s + 1];
        PetscInt kc = flatI[3*s + 2];
        PetscReal px = flatC[3*s + 0];
        PetscReal py = flatC[3*s + 1];
        PetscReal pz = flatC[3*s + 2];

        Cmpnts ucart;
        vectorPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, lucatD, ucart);
        rootBuf[s*stride + 0] = ucart.x;
        rootBuf[s*stride + 1] = ucart.y;
        rootBuf[s*stride + 2] = ucart.z;

        if (flags->isTeqnActive)
        {
            PetscReal Temp = 0.0;
            scalarPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, ltempD, Temp);
            rootBuf[s*stride + 3] = Temp;
        }

        PetscReal Pres = 0.0;
        scalarPointLocalVolumeInterpolation(meshD, px, py, pz, ic, jc, kc, cent, lpressD, Pres);
        rootBuf[s*stride + 4] = Pres;
    }

    const std::vector<Acell> &acceptors = os->localAcceptorsHc[donorId];
    PetscInt nLeaves = (PetscInt)acceptors.size();
    std::vector<PetscReal> leafBuf((size_t)nLeaves * stride, 0.0);

    MPI_Datatype t5;
    MPI_Type_contiguous(5, MPIU_REAL, &t5);
    MPI_Type_commit(&t5);
    PetscSFBcastBegin(os->sfC2P[donorId], t5, rootBuf.data(), leafBuf.data(), MPI_REPLACE);
    PetscSFBcastEnd  (os->sfC2P[donorId], t5, rootBuf.data(), leafBuf.data(), MPI_REPLACE);
    MPI_Type_free(&t5);

    // Acceptor side: group the 8 vertices per parentCellId and average to the cell centre
    const std::vector<Dcell> &donorMap = os->localDonorMapHc[donorId];
    std::map<PetscInt, std::vector<Cmpnts>>    vertexVel;
    std::map<PetscInt, std::vector<PetscReal>> vertexTemp;
    std::map<PetscInt, std::vector<PetscReal>> vertexPres;
    std::map<PetscInt, std::tuple<PetscInt, PetscInt, PetscInt>> cellIdx;

    for (PetscInt b = 0; b < nLeaves; b++)
    {
        if (donorMap[b].rank < 0)
        {
            PetscPrintf(PETSC_COMM_SELF,
                        "Warning: C2P interp no donor for vertex of parent cell (%ld,%ld,%ld) parentCellId=%ld\n",
                        acceptors[b].indi, acceptors[b].indj, acceptors[b].indk, acceptors[b].parentCellId);
            continue;
        }
        Cmpnts ucart;
        ucart.x = leafBuf[b*stride + 0];
        ucart.y = leafBuf[b*stride + 1];
        ucart.z = leafBuf[b*stride + 2];
        vertexVel[acceptors[b].parentCellId].push_back(ucart);
        if (flags->isTeqnActive) vertexTemp[acceptors[b].parentCellId].push_back(leafBuf[b*stride + 3]);
        vertexPres[acceptors[b].parentCellId].push_back(leafBuf[b*stride + 4]);
        if (cellIdx.find(acceptors[b].parentCellId) == cellIdx.end())
        {
            cellIdx[acceptors[b].parentCellId] = {acceptors[b].indi, acceptors[b].indj, acceptors[b].indk};
        }
    }

    for (const auto& kv : vertexVel)
    {
        PetscInt parentCellId = kv.first;
        const std::vector<Cmpnts> &vels = kv.second;
        if (vels.size() != 8)
        {
            PetscPrintf(PETSC_COMM_SELF,
                        "Warning: parentCellId %ld has %lu vertices instead of 8 skipping\n",
                        (long)parentCellId, (unsigned long)vels.size());
            continue;
        }
        Cmpnts avgV = {0.0, 0.0, 0.0};
        for (const auto& v : vels) { avgV.x += v.x/8.0; avgV.y += v.y/8.0; avgV.z += v.z/8.0; }
        auto [i, j, k] = cellIdx[parentCellId];
        ucatA[k][j][i].x = avgV.x;
        ucatA[k][j][i].y = avgV.y;
        ucatA[k][j][i].z = avgV.z;
        PetscReal avgP = 0.0;
        for (PetscReal p : vertexPres[parentCellId]) avgP += p/8.0;
        pressA[k][j][i] = avgP;
        if (flags->isTeqnActive)
        {
            PetscReal avgT = 0.0;
            for (PetscReal t : vertexTemp[parentCellId]) avgT += t/8.0;
            tempA[k][j][i] = avgT;
        }
    }

    DMDAVecRestoreArray(fdaD, ueqnD->lUcat, &lucatD);
    DMDAVecRestoreArray(fdaA, ueqnA->Ucat,  &ucatA);
    DMDAVecRestoreArray(fdaD, meshD->lCent, &cent);
    DMDAVecRestoreArray(daD,  peqnD->lP,    &lpressD);
    DMDAVecRestoreArray(daA,  peqnA->P,     &pressA);
    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(daD, teqnD->lTmprt, &ltempD);
        DMDAVecRestoreArray(daA, teqnA->Tmprt,  &tempA);
        DMGlobalToLocalBegin(daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
        DMGlobalToLocalEnd  (daA, teqnA->Tmprt, INSERT_VALUES, teqnA->lTmprt);
    }
    DMGlobalToLocalBegin(daA,  peqnA->P,    INSERT_VALUES, peqnA->lP);
    DMGlobalToLocalEnd  (daA,  peqnA->P,    INSERT_VALUES, peqnA->lP);
    DMGlobalToLocalBegin(fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);
    DMGlobalToLocalEnd  (fdaA, ueqnA->Ucat, INSERT_VALUES, ueqnA->lUcat);

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
    DM            fda   = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***cent;
    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // Each rank rebuilds its own local acceptor list from scratch
    std::vector<Acell> &local = os->localAcceptorsDb;
    local.clear();

    for (PetscInt k = zs; k < ze; k++)
    {
        for (PetscInt j = ys; j < ye; j++)
        {
            for (PetscInt i = xs; i < xe; i++)
            {
                if (isOnCornerCellCenters(i, j, k, info)) continue;
                if ((k == 0) || (k == mz - 1) || (j == 0) || (j == my - 1) || (i == 0) || (i == mx - 1))
                {
                    Acell c;
                    c.indi  = i; c.indj = j; c.indk = k;
                    c.coorx = cent[k][j][i].x;
                    c.coory = cent[k][j][i].y;
                    c.coorz = cent[k][j][i].z;
                    c.rank  = rank;
                    c.cell_size    = 0.0;
                    c.face         = 0;
                    c.donorId      = 0;
                    c.parentCellId = -1;
                    local.push_back(c);
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode createAcceptorCellBackground(overset_ *os, PetscInt donorMeshId)
{
    mesh_         *mesh = os->access->mesh;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs = xs, lxe = xe; if (xs == 0) lxs = xs + 1; if (xe == mx) lxe = xe - 1;
    PetscInt lys = ys, lye = ye; if (ys == 0) lys = ys + 1; if (ye == my) lye = ye - 1;
    PetscInt lzs = zs, lze = ze; if (zs == 0) lzs = zs + 1; if (ze == mz) lze = ze - 1;

    Vec           Coor;
    PetscReal     ***meshTag;
    Cmpnts        ***coor;
    PetscMPIInt   rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // A cell may only have one donor: skip any cell already claimed by a previously registered donor mesh
    auto alreadyTaken = [&](PetscInt i, PetscInt j, PetscInt k) -> bool
    {
        for (const auto &kv : os->localAcceptorsHc)
        {
            for (const Acell &c : kv.second)
            {
                if (c.indi == i && c.indj == j && c.indk == k && c.rank == rank) return true;
            }
        }
        return false;
    };

    std::vector<Acell> &dst = os->localAcceptorsHc[donorMeshId];
    PetscInt localCellId = (PetscInt)(dst.size() / 8);

    for (PetscInt k = lzs; k < lze; k++)
    {
        for (PetscInt j = lys; j < lye; j++)
        {
            for (PetscInt i = lxs; i < lxe; i++)
            {
                if (!isInterpolatedCell(k, j, i, meshTag)) continue;
                if (alreadyTaken(i, j, k)) continue;

                PetscReal vx[8], vy[8], vz[8];
                vx[0] = coor[k-1][j-1][i-1].x; vy[0] = coor[k-1][j-1][i-1].y; vz[0] = coor[k-1][j-1][i-1].z;
                vx[1] = coor[k-1][j-1][i  ].x; vy[1] = coor[k-1][j-1][i  ].y; vz[1] = coor[k-1][j-1][i  ].z;
                vx[2] = coor[k-1][j  ][i-1].x; vy[2] = coor[k-1][j  ][i-1].y; vz[2] = coor[k-1][j  ][i-1].z;
                vx[3] = coor[k-1][j  ][i  ].x; vy[3] = coor[k-1][j  ][i  ].y; vz[3] = coor[k-1][j  ][i  ].z;
                vx[4] = coor[k  ][j-1][i-1].x; vy[4] = coor[k  ][j-1][i-1].y; vz[4] = coor[k  ][j-1][i-1].z;
                vx[5] = coor[k  ][j-1][i  ].x; vy[5] = coor[k  ][j-1][i  ].y; vz[5] = coor[k  ][j-1][i  ].z;
                vx[6] = coor[k  ][j  ][i-1].x; vy[6] = coor[k  ][j  ][i-1].y; vz[6] = coor[k  ][j  ][i-1].z;
                vx[7] = coor[k  ][j  ][i  ].x; vy[7] = coor[k  ][j  ][i  ].y; vz[7] = coor[k  ][j  ][i  ].z;

                // parentCellId only needs to be unique within this rank (lists are local)
                PetscInt parentCellId = localCellId++;

                for (PetscInt v = 0; v < 8; v++)
                {
                    Acell c;
                    c.indi = i; c.indj = j; c.indk = k;
                    c.coorx = vx[v]; c.coory = vy[v]; c.coorz = vz[v];
                    c.rank  = rank;
                    c.cell_size    = 0.0;
                    c.face         = 0;
                    c.donorId      = donorMeshId;
                    c.parentCellId = parentCellId;
                    dst.push_back(c);
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, Coor, &coor);
    return 0;
}

// experimental octree structure for donor cells
struct OctreeNode
{
    PetscInt hasCells;
    PetscInt imin, imax, jmin, jmax, kmin, kmax; // node bounding indices
    Cmpnts minBounds;                            // node minimum physical bounds
    Cmpnts maxBounds;                            // node maximum physical bounds
    OctreeNode* children[8];                     // pointed to the 8 child nodes

    // initialize pointers to null 
    OctreeNode(Cmpnts minB, Cmpnts maxB) : minBounds(minB), maxBounds(maxB) 
    {
        for (int i = 0; i < 8; i++) children[i] = nullptr;
    }

    // destructor to delete child nodes
    ~OctreeNode() 
    {
        for (int i = 0; i < 8; i++) 
        {
            if (children[i]) delete children[i];
        }
    }
};

void buildOctree
(
    OctreeNode* node, Cmpnts*** donorCells, 
    PetscInt lxs, PetscInt lxe, PetscInt lys, PetscInt lye, PetscInt lzs, PetscInt lze,
    PetscInt maxDepth, PetscInt maxCellsPerNode
) 
{
    // initialize minmax node indices 
    PetscInt imin = 1000000, imax = 0;
    PetscInt jmin = 1000000, jmax = 0;
    PetscInt kmin = 1000000, kmax = 0;

    // count how many cells are in this node 
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

            // set minmax node indices 
            if(i < imin) imin = i;
            if(i > imax) imax = i;
            if(j < jmin) jmin = j;
            if(j > jmax) jmax = j;
            if(k < kmin) kmin = k;
            if(k > kmax) kmax = k;
        }
    }

    // if the number of cells is below the threshold or max depth is reached, this is a leaf node 
    if (cellCount <= maxCellsPerNode || maxDepth == 0) 
    {
        node->hasCells = 1; 
        node->imin     = imin;
        node->imax     = imax;
        node->jmin     = jmin;
        node->jmax     = jmax;
        node->kmin     = kmin;
        node->kmax     = kmax;
        
        return;
    }

    node->hasCells = 0; 

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

Dcell searchOctree
(
    OctreeNode* node, PetscReal procContrib, const Cmpnts& acceptorCoord, Cmpnts*** centroids, PetscReal& minDist,
    PetscInt lxs, PetscInt lxe, PetscInt lys, PetscInt lye, PetscInt lzs, PetscInt lze
) 
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
    
    // we are in a leaf node 
    if(node->hasCells)
    {
        for (PetscInt k = node->kmin; k <= node->kmax; k++)
        for (PetscInt j = node->jmin; j <= node->jmax; j++)
        for (PetscInt i = node->imin; i <= node->imax; i++)
        {
            // Calculate distance to acceptor
            Cmpnts centroid = centroids[k][j][i];

            PetscReal dist = sqrt(pow(centroid.x - acceptorCoord.x - procContrib, 2) +
                                  pow(centroid.y - acceptorCoord.y - procContrib, 2) +
                                  pow(centroid.z - acceptorCoord.z - procContrib, 2));

            if (dist < minDist) 
            {
                minDist             = dist + procContrib;
                closestDonor.indi   = i;
                closestDonor.indj   = j;
                closestDonor.indk   = k;
                closestDonor.dist2p = dist;
                closestDonor.rank   = 1;
            }
        }
    }

    // recursively search child nodes
    for (int i = 0; i < 8; i++) 
    {
        if (node->children[i] != nullptr) 
        {
            PetscReal childMinDist = minDist;
            Dcell childClosest = searchOctree(node->children[i], procContrib, acceptorCoord, centroids, childMinDist, lxs, lxe, lys, lye, lzs, lze);

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

//! \brief Per-rank axis-aligned bounding box of donor mesh cell centres.
// local struct and using static only visible in overset.c (may be changed based on the scope)
struct OversetBbox
{
    PetscReal xmin, ymin, zmin;
    PetscReal xmax, ymax, zmax;
};

//! \brief Returns true if (x,y,z) lies inside the bounding box with a small tolerance.
static inline bool bboxContains(const OversetBbox &bb, PetscReal x, PetscReal y, PetscReal z, PetscReal tol)
{
    return x >= bb.xmin - tol && x <= bb.xmax + tol &&
           y >= bb.ymin - tol && y <= bb.ymax + tol &&
           z >= bb.zmin - tol && z <= bb.zmax + tol;
}

//! \brief Query packet: acceptor rank asks candidate donor ranks for its closest donor.
struct OversetQuery
{
    PetscInt  originIdx;
    PetscReal x, y, z;
};

//! \brief Reply packet: donor rank reports its best candidate slot and distance back.
struct OversetReply
{
    PetscInt  originIdx;
    PetscInt  slotIdx;
    PetscReal dist;
};

//! \brief Builds the PetscSF (and matching local donor map) for one donor/acceptor mesh pair.
static PetscErrorCode BuildOversetSF(
                                        mesh_                    *meshDonor,
                                        mesh_                    *meshAcceptor,
                                        const std::vector<Acell> &localAcceptors,
                                        std::vector<Dcell>       &localDonorMap,
                                        PetscInt                 &nRoots,
                                        std::vector<PetscInt>    &rootI,
                                        std::vector<PetscInt>    &rootJ,
                                        std::vector<PetscInt>    &rootK,
                                        std::vector<PetscReal>   &rootX,
                                        std::vector<PetscReal>   &rootY,
                                        std::vector<PetscReal>   &rootZ,
                                        PetscSF                  &sf
                                    )
{
    DM            da   = meshDonor->da, fda = meshDonor->fda;
    DMDALocalInfo info = meshDonor->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs = xs, lxe = xe; if (xs == 0) lxs = xs + 1; if (xe == mx) lxe = xe - 1;
    PetscInt lys = ys, lye = ye; if (ys == 0) lys = ys + 1; if (ye == my) lye = ye - 1;
    PetscInt lzs = zs, lze = ze; if (zs == 0) lzs = zs + 1; if (ze == mz) lze = ze - 1;

    MPI_Comm    comm = meshDonor->MESH_COMM;
    PetscMPIInt rankD, sizeD, sizeA;
    MPI_Comm_size(comm, &sizeD);
    MPI_Comm_rank(comm, &rankD);
    MPI_Comm_size(meshAcceptor->MESH_COMM, &sizeA);

    if (sizeA != sizeD)
    {
        char error[512];
        sprintf(error, "Donor and acceptor meshes must share the same MPI communicator size.\n");
        fatalErrorInFunction("BuildOversetSF", error);
    }

    Vec           Coor;
    Cmpnts        ***cent, ***coor;
    PetscReal     ***aj;
    DMDAVecGetArray(fda, meshDonor->lCent, &cent);
    DMDAVecGetArray(da, meshDonor->lAj, &aj);
    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    // Expand the donor bounding box by one cell so neighbouring ranks overlap
    PetscInt bsz = zs; if (zs != 0) bsz = bsz - 1;
    PetscInt bsy = ys; if (ys != 0) bsy = bsy - 1;
    PetscInt bsx = xs; if (xs != 0) bsx = bsx - 1;
    Cmpnts procDeltas = {fabs(coor[lze-1][lye-1][lxe-1].x - coor[bsz][bsy][bsx].x),
                         fabs(coor[lze-1][lye-1][lxe-1].y - coor[bsz][bsy][bsx].y),
                         fabs(coor[lze-1][lye-1][lxe-1].z - coor[bsz][bsy][bsx].z)};
    Cmpnts oMin = {coor[bsz][bsy][bsx].x - procDeltas.x,
                   coor[bsz][bsy][bsx].y - procDeltas.y,
                   coor[bsz][bsy][bsx].z - procDeltas.z};
    Cmpnts oMax = {coor[lze-1][lye-1][lxe-1].x + procDeltas.x,
                   coor[lze-1][lye-1][lxe-1].y + procDeltas.y,
                   coor[lze-1][lye-1][lxe-1].z + procDeltas.z};

    // Build the local donor-side octree
    OctreeNode *root = new OctreeNode(oMin, oMax);
    buildOctree(root, cent, lxs, lxe, lys, lye, lzs, lze, 15, 1000);

    // Allgather donor bboxes across all ranks
    OversetBbox myBb;
    myBb.xmin = oMin.x; myBb.xmax = oMax.x;
    myBb.ymin = oMin.y; myBb.ymax = oMax.y;
    myBb.zmin = oMin.z; myBb.zmax = oMax.z;
    std::vector<OversetBbox> allBb(sizeD);
    MPI_Allgather(&myBb, sizeof(OversetBbox), MPI_BYTE,
                  allBb.data(), sizeof(OversetBbox), MPI_BYTE, comm);

    // Count candidate donor ranks per local acceptor (those whose bbox contains the acceptor)
    PetscInt nLocal = (PetscInt)localAcceptors.size();
    std::vector<int> sCount(sizeD, 0), sDispl(sizeD, 0);
    std::vector<int> rCount(sizeD, 0), rDispl(sizeD, 0);

    for (PetscInt b = 0; b < nLocal; b++)
    {
        PetscReal x = localAcceptors[b].coorx;
        PetscReal y = localAcceptors[b].coory;
        PetscReal z = localAcceptors[b].coorz;
        for (PetscMPIInt d = 0; d < sizeD; d++)
        {
            if (bboxContains(allBb[d], x, y, z, 0.0)) sCount[d]++;
        }
    }
    for (PetscMPIInt d = 1; d < sizeD; d++) sDispl[d] = sDispl[d-1] + sCount[d-1];
    PetscInt totalSend = (sizeD > 0) ? sDispl[sizeD-1] + sCount[sizeD-1] : 0;

    // Pack the queries grouped by destination rank
    std::vector<OversetQuery> sendQ(totalSend);
    std::vector<PetscInt> cursor(sizeD, 0);
    for (PetscInt b = 0; b < nLocal; b++)
    {
        PetscReal x = localAcceptors[b].coorx;
        PetscReal y = localAcceptors[b].coory;
        PetscReal z = localAcceptors[b].coorz;
        for (PetscMPIInt d = 0; d < sizeD; d++)
        {
            if (bboxContains(allBb[d], x, y, z, 0.0))
            {
                PetscInt pos = sDispl[d] + cursor[d];
                sendQ[pos].originIdx = b;
                sendQ[pos].x = x; sendQ[pos].y = y; sendQ[pos].z = z;
                cursor[d]++;
            }
        }
    }

    // Exchange counts, then queries, via Alltoall/Alltoallv
    MPI_Alltoall(sCount.data(), 1, MPI_INT, rCount.data(), 1, MPI_INT, comm);
    for (PetscMPIInt d = 1; d < sizeD; d++) rDispl[d] = rDispl[d-1] + rCount[d-1];
    PetscInt totalRecv = (sizeD > 0) ? rDispl[sizeD-1] + rCount[sizeD-1] : 0;

    MPI_Datatype qtype;
    {
        int blocklens[2] = {1, 3};
        MPI_Aint disp[2];
        OversetQuery tmp;
        MPI_Aint base;
        MPI_Get_address(&tmp, &base);
        MPI_Get_address(&tmp.originIdx, &disp[0]);
        MPI_Get_address(&tmp.x, &disp[1]);
        disp[0] -= base; disp[1] -= base;
        MPI_Datatype types[2] = {MPIU_INT, MPIU_REAL};
        MPI_Type_create_struct(2, blocklens, disp, types, &qtype);
        MPI_Type_commit(&qtype);
    }

    std::vector<OversetQuery> recvQ(totalRecv);
    MPI_Alltoallv(sendQ.data(), sCount.data(), sDispl.data(), qtype,
                  recvQ.data(), rCount.data(), rDispl.data(), qtype, comm);

    // Donor side: search the octree for each incoming query and allocate a root slot
    rootI.clear(); rootJ.clear(); rootK.clear();
    rootX.clear(); rootY.clear(); rootZ.clear();
    rootI.reserve(totalRecv); rootJ.reserve(totalRecv); rootK.reserve(totalRecv);
    rootX.reserve(totalRecv); rootY.reserve(totalRecv); rootZ.reserve(totalRecv);
    std::vector<OversetReply> reply(totalRecv);

    PetscInt slotCounter = 0;
    for (PetscInt q = 0; q < totalRecv; q++)
    {
        Cmpnts pt = nSetFromComponents(recvQ[q].x, recvQ[q].y, recvQ[q].z);
        PetscReal lminDist = 1e20;
        Dcell res = searchOctree(root, 0.0, pt, cent, lminDist, lxs, lxe, lys, lye, lzs, lze);

        reply[q].originIdx = recvQ[q].originIdx;
        if (res.rank == -1)
        {
            reply[q].slotIdx = -1;
            reply[q].dist    = 1e20;
        }
        else
        {
            reply[q].slotIdx = slotCounter;
            reply[q].dist    = res.dist2p;
            rootI.push_back(res.indi);
            rootJ.push_back(res.indj);
            rootK.push_back(res.indk);
            // Cache acceptor query coords for trilinear interpolation at every step
            rootX.push_back(recvQ[q].x);
            rootY.push_back(recvQ[q].y);
            rootZ.push_back(recvQ[q].z);
            slotCounter++;
        }
    }
    nRoots = slotCounter;

    // Send replies back; send and recv counts are swapped for the return trip
    MPI_Datatype rtype;
    {
        int blocklens[3] = {1, 1, 1};
        MPI_Aint disp[3];
        OversetReply tmp;
        MPI_Aint base;
        MPI_Get_address(&tmp, &base);
        MPI_Get_address(&tmp.originIdx, &disp[0]);
        MPI_Get_address(&tmp.slotIdx,   &disp[1]);
        MPI_Get_address(&tmp.dist,      &disp[2]);
        disp[0] -= base; disp[1] -= base; disp[2] -= base;
        MPI_Datatype types[3] = {MPIU_INT, MPIU_INT, MPIU_REAL};
        MPI_Type_create_struct(3, blocklens, disp, types, &rtype);
        MPI_Type_commit(&rtype);
    }
    std::vector<OversetReply> recvReply(totalSend);
    MPI_Alltoallv(reply.data(),    rCount.data(), rDispl.data(), rtype,
                  recvReply.data(), sCount.data(), sDispl.data(), rtype, comm);

    // Pick the minimum-distance donor for each local acceptor
    std::vector<PetscReal> bestDist(nLocal, 1e20);
    std::vector<PetscInt>  bestRank(nLocal, -1);
    std::vector<PetscInt>  bestSlot(nLocal, -1);
    for (PetscMPIInt d = 0; d < sizeD; d++)
    {
        PetscInt off = sDispl[d];
        PetscInt cnt = sCount[d];
        for (PetscInt k = 0; k < cnt; k++)
        {
            const OversetReply &r = recvReply[off + k];
            if (r.slotIdx < 0) continue;
            if (r.dist < bestDist[r.originIdx])
            {
                bestDist[r.originIdx] = r.dist;
                bestRank[r.originIdx] = d;
                bestSlot[r.originIdx] = r.slotIdx;
            }
        }
    }

    // Build the local donor map and the PetscSF graph
    localDonorMap.assign(nLocal, Dcell());
    PetscInt nLeaves = 0;
    for (PetscInt b = 0; b < nLocal; b++) if (bestRank[b] >= 0) nLeaves++;

    std::vector<PetscInt>    ilocal(nLeaves);
    std::vector<PetscSFNode> iremote(nLeaves);
    PetscInt leafCount = 0;
    for (PetscInt b = 0; b < nLocal; b++)
    {
        Dcell &dc = localDonorMap[b];
        if (bestRank[b] >= 0)
        {
            ilocal[leafCount] = b;
            iremote[leafCount].rank  = (PetscMPIInt)bestRank[b];
            iremote[leafCount].index = bestSlot[b];
            leafCount++;
            dc.rank   = (PetscMPIInt)bestRank[b];
            dc.dist2p = bestDist[b];
            dc.indi   = -1; dc.indj = -1; dc.indk = -1;
        }
        else
        {
            dc.rank   = -1;
            dc.dist2p = 1e20;
            dc.indi   = -1; dc.indj = -1; dc.indk = -1;
        }
    }

    PetscSFCreate(comm, &sf);
    PetscSFSetGraph(sf, nRoots, nLeaves,
                    ilocal.data(),  PETSC_COPY_VALUES,
                    iremote.data(), PETSC_COPY_VALUES);
    PetscSFSetUp(sf);

    MPI_Type_free(&qtype);
    MPI_Type_free(&rtype);
    delete root;

    DMDAVecRestoreArray(fda, Coor, &coor);
    DMDAVecRestoreArray(fda, meshDonor->lCent, &cent);
    DMDAVecRestoreArray(da, meshDonor->lAj, &aj);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode findClosestDonorC2P(mesh_ *meshDonor, mesh_ *meshAcceptor, PetscInt donorId)
{
    overset_ *os = meshAcceptor->access->os;

    // No acceptors registered for this donor mesh: install an empty PetscSF and return
    auto it = os->localAcceptorsHc.find(donorId);
    if (it == os->localAcceptorsHc.end())
    {
        if (os->sfC2P.count(donorId) && os->sfC2P[donorId]) PetscSFDestroy(&os->sfC2P[donorId]);
        PetscSFCreate(meshDonor->MESH_COMM, &os->sfC2P[donorId]);
        PetscSFSetGraph(os->sfC2P[donorId], 0, 0, NULL, PETSC_COPY_VALUES, NULL, PETSC_COPY_VALUES);
        os->nRootsHc[donorId] = 0;
        os->rootSlotIdxHc[donorId].clear();
        os->rootSlotCoordsHc[donorId].clear();
        os->localDonorMapHc[donorId].clear();
        return 0;
    }

    std::vector<PetscInt>  rootI, rootJ, rootK;
    std::vector<PetscReal> rootX, rootY, rootZ;
    PetscInt               nRoots = 0;
    PetscSF                sf     = NULL;

    // Build the PetscSF for this (donor, acceptor) pair
    PetscErrorCode ierr = BuildOversetSF(meshDonor, meshAcceptor,
                                         it->second,
                                         os->localDonorMapHc[donorId],
                                         nRoots,
                                         rootI, rootJ, rootK,
                                         rootX, rootY, rootZ,
                                         sf);
    CHKERRQ(ierr);

    os->nRootsHc[donorId] = nRoots;
    std::vector<PetscInt>  &flatI = os->rootSlotIdxHc[donorId];
    std::vector<PetscReal> &flatC = os->rootSlotCoordsHc[donorId];
    flatI.resize(3 * nRoots);
    flatC.resize(3 * nRoots);
    for (PetscInt s = 0; s < nRoots; s++)
    {
        flatI[3*s + 0] = rootI[s];
        flatI[3*s + 1] = rootJ[s];
        flatI[3*s + 2] = rootK[s];
        flatC[3*s + 0] = rootX[s];
        flatC[3*s + 1] = rootY[s];
        flatC[3*s + 2] = rootZ[s];
    }
    if (os->sfC2P.count(donorId) && os->sfC2P[donorId]) PetscSFDestroy(&os->sfC2P[donorId]);
    os->sfC2P[donorId] = sf;

    PetscInt nMissing = 0;
    for (const auto &dc : os->localDonorMapHc[donorId]) if (dc.rank < 0) nMissing++;
    PetscInt gMissing = 0;
    MPI_Allreduce(&nMissing, &gMissing, 1, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    if (gMissing > 0)
    {
        PetscPrintf(meshDonor->MESH_COMM,
                    "     C2P donor search for donorId %ld left %ld vertices without a donor\n",
                    (long)donorId, (long)gMissing);
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode findClosestDonorP2C(mesh_ *meshDonor, mesh_ *meshAcceptor)
{
    overset_ *os = meshAcceptor->access->os;

    std::vector<PetscInt>  rootI, rootJ, rootK;
    std::vector<PetscReal> rootX, rootY, rootZ;
    PetscInt               nRoots = 0;
    PetscSF                sf     = NULL;

    PetscErrorCode ierr = BuildOversetSF(meshDonor, meshAcceptor,
                                         os->localAcceptorsDb,
                                         os->localDonorMapDb,
                                         nRoots,
                                         rootI, rootJ, rootK,
                                         rootX, rootY, rootZ,
                                         sf);
    CHKERRQ(ierr);

    os->nRootsDb = nRoots;
    // Pack (i,j,k) and (x,y,z) into flat stride-3 arrays
    os->rootSlotIdxDb.resize(3 * nRoots);
    os->rootSlotCoordsDb.resize(3 * nRoots);
    for (PetscInt s = 0; s < nRoots; s++)
    {
        os->rootSlotIdxDb[3*s + 0] = rootI[s];
        os->rootSlotIdxDb[3*s + 1] = rootJ[s];
        os->rootSlotIdxDb[3*s + 2] = rootK[s];
        os->rootSlotCoordsDb[3*s + 0] = rootX[s];
        os->rootSlotCoordsDb[3*s + 1] = rootY[s];
        os->rootSlotCoordsDb[3*s + 2] = rootZ[s];
    }
    // Destroy any previous PetscSF before overwriting
    if (os->sfP2C) PetscSFDestroy(&os->sfP2C);
    os->sfP2C = sf;

    // Report acceptors that did not find a donor (likely a mesh-overlap issue)
    PetscMPIInt rankD;
    MPI_Comm_rank(meshDonor->MESH_COMM, &rankD);
    PetscInt nMissing = 0;
    for (const auto &dc : os->localDonorMapDb) if (dc.rank < 0) nMissing++;
    PetscInt gMissing = 0;
    MPI_Allreduce(&nMissing, &gMissing, 1, MPIU_INT, MPI_SUM, meshDonor->MESH_COMM);
    if (gMissing > 0)
    {
        PetscPrintf(meshDonor->MESH_COMM,
                    "     P2C donor search left %ld acceptors without a donor check overlap of meshes\n",
                    (long)gMissing);
    }

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
  DMDALocalInfo   info = mesh->info;
  PetscInt        xs = info.xs, xe = info.xs + info.xm;
  PetscInt        ys = info.ys, ye = info.ys + info.ym;
  PetscInt        zs = info.zs, ze = info.zs + info.zm;
  PetscInt        mx = info.mx, my = info.my, mz = info.mz;
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
            if
            (
                cent[k][j][i].x > ibBox->xmin
                && cent[k][j][i].x < ibBox->xmax
                && cent[k][j][i].y > ibBox->ymin
                && cent[k][j][i].y < ibBox->ymax
                && cent[k][j][i].z > ibBox->zmin
                && cent[k][j][i].z < ibBox->zmax 
            )
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

    // set boundary conditions for meshTag: 
    // if the boundary is set to overset interpolate then set it to 1 
    DMDAVecGetArray(da, mesh->meshTag, &gmeshTag);
    for (k = zs; k < ze; k++)
    for (j = ys; j < ye; j++)
    for (i = xs; i < xe; i++)
    {
        if(i==0 && mesh->boundaryU.iLeft == "oversetInterpolate")
        {
            gmeshTag[k][j][i] = 1;
        }
        if(i==mx-1 && mesh->boundaryU.iRight == "oversetInterpolate")
        {
            gmeshTag[k][j][i] = 1;
        }
        if(j==0 && mesh->boundaryU.jLeft == "oversetInterpolate")
        {
            gmeshTag[k][j][i] = 1;
        }
        if(j==my-1 && mesh->boundaryU.jRight == "oversetInterpolate")   
        {
            gmeshTag[k][j][i] = 1;
        }
        if(k==0 && mesh->boundaryU.kLeft == "oversetInterpolate")
        {
            gmeshTag[k][j][i] = 1;
        }
        if(k==mz-1 && mesh->boundaryU.kRight == "oversetInterpolate")
        {
            gmeshTag[k][j][i] = 1;
        }
    }

    DMDAVecRestoreArray(da, mesh->meshTag, &gmeshTag);

    DMGlobalToLocalBegin(da, mesh->meshTag, INSERT_VALUES, mesh->lmeshTag);
    DMGlobalToLocalEnd(da, mesh->meshTag, INSERT_VALUES, mesh->lmeshTag);
           
    return 0;
}

//***************************************************************************************************************//

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

// DEPRECATED FUNCTION: for other interpolation methods. Variables names need to be updated for integration.
//***************************************************************************************************************//

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
