//! \file  mesh.c
//! \brief Contains mesh definition functions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

//***************************************************************************************************************//

PetscErrorCode InitializeMesh(mesh_ *mesh)
{
    // set distributed arrays
    SetDistributedArrays(mesh);

    // bounding box initialize
    SetBoundingBox(mesh);

    // deform mesh for BL displacement if gravity waves are modelled
    if(mesh->access->flags->isGravityWaveModelingActive)
    {
        DeformMeshBasedOnBLDisp(mesh);
    }

    // set curvilinear coordinates metrics
    SetMeshMetrics(mesh);

    MPI_Barrier(mesh->MESH_COMM);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetDistributedArrays(mesh_ *mesh)
{
    PetscInt             i, j, k, l;
    word                 meshFileName;
    PetscReal            bufferDouble;
    PetscInt             bufferInt;
    char                 bufferChar;
    std::vector<PetscReal>  Ycart, Zcart, Xcart;
    PetscInt             m, n, p, s = 3;
    PetscSubcomm         MeshComm;
    PetscInt             error;
    PetscMPIInt          rank, nProcs;
    Vec                  lCoor, gCoor;
    Cmpnts               ***lcoor, ***gcoor;

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

    // read type of mesh file provided
    readDictWord("control.dat", "-meshFileType", &(mesh->meshFileType));

    FILE *meshFileID;

    // only specify axes
    if(mesh->meshFileType == "cartesian")
    {
        if(mesh->meshName == ".") meshFileName ="mesh.xyz";
        else                      meshFileName = mesh->meshName + ".xyz";

        SetPeriodicConnectivity(mesh, meshFileName);

        meshFileID = fopen(meshFileName.c_str(), "r");

        if(!meshFileID)
        {
            char error[512];
            sprintf(error, "cannot open mesh file %s\n", meshFileName.c_str());
            fatalErrorInFunction("SetDistributedArrays", error);
        }
        else
        {
            // skip lines indicating periodic connectivity
            for(l=0; l<mesh->access->info->periodic; l++)
            {
                error = fscanf(meshFileID, "%s %ld", &bufferChar, &bufferInt);
            }

            // first line contains number of nodes in x,y,z directions
            PetscInt npx, npy, npz;
            error = fscanf(meshFileID, "%ld %ld %ld\n", &npx, &npy, &npz);

            Xcart.resize(npx);
            Ycart.resize(npy);
            Zcart.resize(npz);

            for (k = 0; k < npx; k++) error = fscanf(meshFileID, "%le %le %le\n", &Xcart[k], &bufferDouble, &bufferDouble);
            for (i = 0; i < npy; i++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble,  &Ycart[i], &bufferDouble);
            for (j = 0; j < npz; j++) error = fscanf(meshFileID, "%le %le %le\n", &bufferDouble, &bufferDouble,  &Zcart[j]);

            mesh->IM = npy + 1;
            mesh->JM = npz + 1;
            mesh->KM = npx + 1;
        }

        fclose(meshFileID);
    }
    // specify all points
    else if(mesh->meshFileType == "curvilinear")
    {
        if(mesh->meshName == ".") meshFileName ="mesh.grid";
        else                      meshFileName = mesh->meshName + ".grid";

        SetPeriodicConnectivity(mesh, meshFileName);

        meshFileID = fopen(meshFileName.c_str(), "r");

        if(!meshFileID)
        {
           char error[512];
            sprintf(error, "cannot open mesh file %s\n", meshFileName.c_str());
            fatalErrorInFunction("SetDistributedArrays", error);
        }
        else
        {
            // skip lines indicating periodic connectivity
            for(l=0; l<mesh->access->info->periodic; l++)
            {
                error = fscanf(meshFileID, "%s %ld", &bufferChar, &bufferInt);
            }

            // first line contains number of nodes in x,y,z directions
            PetscInt npx, npy, npz;
            error = fscanf(meshFileID, "%ld %ld %ld\n", &npz, &npx, &npy);

            mesh->IM = npy + 1;
            mesh->JM = npz + 1;
            mesh->KM = npx + 1;
        }

    }

    // initialize domain boundary to ghost cells
    DMBoundaryType bx = DM_BOUNDARY_GHOSTED, by = DM_BOUNDARY_GHOSTED, bz = DM_BOUNDARY_GHOSTED;

    // correct domain boundary if periodic
    if (mesh->ii_periodic) bx = DM_BOUNDARY_PERIODIC;
    if (mesh->jj_periodic) by = DM_BOUNDARY_PERIODIC;
    if (mesh->kk_periodic) bz = DM_BOUNDARY_PERIODIC;

    m = n = p = PETSC_DECIDE;

    if (mesh->i_periodic) m = 1;
    if (mesh->j_periodic) n = 1;
    if (mesh->k_periodic) p = 1;

    // create new communicator
    PetscSubcommCreate(PETSC_COMM_WORLD, &(MeshComm));
    PetscSubcommSetTypeGeneral(MeshComm, 1, rank);
    mesh->MESH_COMM = PetscSubcommChild(MeshComm);

    // create da: distributed array data structure for scalars
    DMDACreate3d(mesh->MESH_COMM, bx, by, bz,                   // boundary type
            DMDA_STENCIL_BOX,                                   // stencil type
            mesh->IM, mesh->JM, mesh->KM,                       // global points in x,y,z
            m, n, p,                                            // processes in x,y,z
            1,                                                  // dofs
            s,                                                  // ghost stencil width
            PETSC_NULL, PETSC_NULL, PETSC_NULL,                 // lx,ly,lz
            &(mesh->da));

    // set up da
    DMSetUp(mesh->da);

    // set coordinates on da (xmin,xmax,ymin,ymax,zmin,zmax)
    DMDASetUniformCoordinates(mesh->da, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);

    // create fda: distributed array data structure for vectors
    DMGetCoordinateDM(mesh->da, &(mesh->fda));

    // get information about da and this processor location in it
    DMDAGetLocalInfo(mesh->da, &(mesh->info));

    // create da: distributed array data structure for symmetric tensors
    DMDACreate3d(mesh->MESH_COMM, bx, by, bz,                   // boundary type
            DMDA_STENCIL_BOX,                                   // stencil type
            mesh->IM, mesh->JM, mesh->KM,                       // global points in x,y,z
            m, n, p,                                            // processes in x,y,z
            6,                                                  // dofs
            s,                                                  // ghost stencil width
            PETSC_NULL, PETSC_NULL, PETSC_NULL,                 // lx,ly,lz
            &(mesh->sda));

    // set up sda
    DMSetUp(mesh->sda);

    MPI_Barrier(mesh->MESH_COMM);

    PetscPrintf(mesh->MESH_COMM, "Reading mesh...");

    DMDALocalInfo info = mesh->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs, lxe, lys, lye, lzs, lze;

    lxe = xe; if (xe == mx) lxe = xe - 1;
    lye = ye; if (ye == my) lye = ye - 1;
    lze = ze; if (ze == mz) lze = ze - 1;

    // get initial global and local co-ordinates from the da
    DMGetCoordinates(mesh->da, &gCoor);
    DMGetCoordinatesLocal(mesh->da, &lCoor);

    //reset the global co-ordinates to 0 and set them from the mesh file
    VecSet(gCoor, 0.);

    DMDAVecGetArray(mesh->fda, gCoor, &gcoor);

    // read x coords - loop over k,i,j with global indexing

    if (mesh->meshFileType == "cartesian")
    {
        for (k = 0; k < mesh->KM; k++)
        {
            for (j = 0; j < mesh->JM; j++)
            {
                for (i = 0; i < mesh->IM; i++)
                {

                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].x = Xcart[k];
                    }
                }
            }
        }

        // read y coords - loop over k,i,j with global indexing
        for (k = 0; k < mesh->KM; k++)
        {
            for (j = 0; j < mesh->JM; j++)
            {
                for (i = 0; i < mesh->IM; i++)
                {

                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].y = Ycart[i];
                    }
                }
            }
        }

        // read z coords - loop over k,i,j with global indexing
        for (k = 0; k < mesh->KM; k++)
        {
            for (j = 0; j < mesh->JM; j++)
            {
                for (i = 0; i < mesh->IM; i++)
                {
                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].z = Zcart[j];
                    }
                }
            }
        }

    }
    else if (mesh->meshFileType == "curvilinear")
    {
        for (i = 0; i < mesh->IM-1; i++)
        {
            for (k = 0; k < mesh->KM-1; k++)
            {
                for (j = 0; j < mesh->JM-1; j++)
                {
                    error = fscanf(meshFileID, "%le", &bufferDouble);

                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].x = bufferDouble;
                    }
                }
            }
        }

            // read y coords - loop over k,i,j with global indexing
        for (i = 0; i < mesh->IM-1; i++)
        {
            for (k = 0; k < mesh->KM-1; k++)
            {
                for (j = 0; j < mesh->JM-1; j++)
                {
                    error = fscanf(meshFileID, "%le", &bufferDouble);

                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].y = bufferDouble;
                    }
                }
            }
        }

        // read z coords - loop over k,i,j with global indexing
        for (i = 0; i < mesh->IM-1; i++)
        {
            for (k = 0; k < mesh->KM-1; k++)
            {
                for (j = 0; j < mesh->JM-1; j++)
                {
                    error = fscanf(meshFileID, "%le", &bufferDouble);

                    if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
                    {
                        gcoor[k][j][i].z = bufferDouble;
                    }
                }
            }
        }

        // for (i = 0; i < mesh->IM-1; i++)
        // {
        //     for (k = mesh->KM-2; k >= 0; k--)
        //     {
        //         for (j = 0; j < mesh->JM-1; j++)
        //         {
        //             error = fscanf(meshFileID, "%le", &bufferDouble);
        //
        //             if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
        //             {
        //                 gcoor[k][j][i].x = bufferDouble;
        //             }
        //         }
        //     }
        // }
        //
        //     // read y coords - loop over k,i,j with global indexing
        // for (i = 0; i < mesh->IM-1; i++)
        // {
        //     for (k = mesh->KM-2; k >= 0; k--)
        //     {
        //         for (j = 0; j < mesh->JM-1; j++)
        //         {
        //             error = fscanf(meshFileID, "%le", &bufferDouble);
        //
        //             if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
        //             {
        //                 gcoor[k][j][i].y = bufferDouble;
        //             }
        //         }
        //     }
        // }
        //
        // // read z coords - loop over k,i,j with global indexing
        // for (i = 0; i < mesh->IM-1; i++)
        // {
        //     for (k = mesh->KM-2; k >= 0; k--)
        //     {
        //         for (j = mesh->JM-2; j >= 0; j--)
        //         {
        //             error = fscanf(meshFileID, "%le", &bufferDouble);
        //
        //             if (k >= zs && k < ze && j >= ys && j < ye && i >= xs && i < xe)
        //             {
        //                 gcoor[k][j][i].z = bufferDouble;
        //             }
        //         }
        //     }
        // }

        fclose(meshFileID);
    }

    DMDAVecRestoreArray(mesh->fda, gCoor, &gcoor);

    //scatter the global co-ordinates to local co-ordinates
    // DMGlobalToLocalBegin(mesh->fda, gCoor, INSERT_VALUES, lCoor);
    // DMGlobalToLocalEnd(mesh->fda, gCoor, INSERT_VALUES, lCoor);

    // set coordinates into da
    DMSetCoordinates(mesh->da, gCoor);
    //DMSetCoordinatesLocal(mesh->da, lCoor);

    MPI_Barrier(mesh->MESH_COMM);

    PetscPrintf(mesh->MESH_COMM, "done\n\n");

    // delete the local mesh co-ordinate vectors
    std::vector<PetscReal> ().swap(Xcart);
    std::vector<PetscReal> ().swap(Ycart);
    std::vector<PetscReal> ().swap(Zcart);

    // initialize mesh and metric local and global fields
    DMCreateGlobalVector(mesh->fda, &(mesh->Cent));
    DMCreateLocalVector (mesh->fda, &(mesh->lCent));

    VecDuplicate(mesh->lCent, &(mesh->lCsi));
    VecDuplicate(mesh->lCent, &(mesh->lEta));
    VecDuplicate(mesh->lCent, &(mesh->lZet));

    VecDuplicate(mesh->lCent, &(mesh->lICsi));
    VecDuplicate(mesh->lCent, &(mesh->lIEta));
    VecDuplicate(mesh->lCent, &(mesh->lIZet));
    VecDuplicate(mesh->lCent, &(mesh->lJCsi));
    VecDuplicate(mesh->lCent, &(mesh->lJEta));
    VecDuplicate(mesh->lCent, &(mesh->lJZet));
    VecDuplicate(mesh->lCent, &(mesh->lKCsi));
    VecDuplicate(mesh->lCent, &(mesh->lKEta));
    VecDuplicate(mesh->lCent, &(mesh->lKZet));

    VecDuplicate(mesh->Cent, &(mesh->fluxLimiter));

    DMCreateGlobalVector(mesh->da,  &(mesh->Nvert));
    DMCreateLocalVector (mesh->da,  &(mesh->lAj));

    VecDuplicate(mesh->lAj, &(mesh->lIAj));
    VecDuplicate(mesh->lAj, &(mesh->lJAj));
    VecDuplicate(mesh->lAj, &(mesh->lKAj));
    VecDuplicate(mesh->lAj, &(mesh->lNvert));
    VecDuplicate(mesh->lAj, &(mesh->lNvert_o));

    VecDuplicate(mesh->Nvert, &(mesh->Nvert_o));

    VecDuplicate(mesh->Nvert, &(mesh->ventMarkers));

    VecSet(mesh->Nvert, 0.0);
    VecSet(mesh->lNvert, 0.0);
    VecSet(mesh->Nvert_o, 0.0);
    VecSet(mesh->lNvert_o, 0.0);
    VecSet(mesh->ventMarkers, 0.0);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DeformMeshBasedOnBLDisp(mesh_ *mesh)
{
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs  = info.xs, xe = info.xs + info.xm;
    PetscInt      ys  = info.ys, ye = info.ys + info.ym;
    PetscInt      zs  = info.zs, ze = info.zs + info.zm;
    PetscInt      gxs = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs = info.gzs, gze = info.gzs + info.gzm;
    PetscInt      mx  = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      error;

    Vec           gCoor, lCoor;

    Cmpnts        ***coor;                     // point coordinates

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // read BL displacement
    FILE       *fileID = fopen("atmosphericSources/eta", "r");

    if(!fileID)
    {
        char error[512];
        sprintf(error, "cannot open atmosphericSources/eta file\n");
        fatalErrorInFunction("DeformMeshBasedOnBLDisp", error);
    }

    PetscReal  dxSources, dySources,
               xsSources, ysSources;
    PetscInt   nxSources, nySources;
    char       bufferChar;

    // read number of points
    error = fscanf(fileID, "%s %ld", &bufferChar, &nxSources);
    error = fscanf(fileID, "%s %ld", &bufferChar, &nySources);

    // read cell size
    error = fscanf(fileID, "%s %lf", &bufferChar, &dxSources);
    error = fscanf(fileID, "%s %lf", &bufferChar, &dySources);

    // read starting coords
    error = fscanf(fileID, "%s %lf", &bufferChar, &xsSources);
    error = fscanf(fileID, "%s %lf", &bufferChar, &ysSources);

    // allocate local vectors
    std::vector<std::vector<PetscReal>>  etaSources(nxSources);

    // read values
    for (i=0; i<nxSources; i++)
    {
        etaSources[i].resize(nySources);

        for (j=0; j<nySources; j++)
        {
            error = fscanf(fileID, "%lf", &(etaSources[i][j]));
        }
    }

    // close file
    fclose(fileID);

    // get global co-ordinates from the da
    DMGetCoordinates(mesh->da, &gCoor);

    DMDAVecGetArray(mesh->fda, gCoor, &coor);

    // only displace z coordinate by eta
    for (k=zs; k<lze; k++)
    {
        for (i=xs; i<lxe; i++)
        {
            // the j index is dummy as the mesh is assumed cartesian
            PetscReal xQuery = coor[k][lys][i].x,
                      yQuery = coor[k][lys][i].y;

            // find bilinear interpolation parameters
            PetscInt kClose = std::floor((xQuery-xsSources)/dxSources);
            PetscInt iClose = std::floor((yQuery-ysSources)/dySources);

            PetscInt kl = kClose,
                     kr = kClose + 1,
                     ib = iClose,
                     it = iClose + 1;

            PetscReal wl = 0.0,
                      wr = 0.0,
                      wb = 0.0,
                      wt = 0.0;

            // out of right boundary
            if(kr>=nxSources)
            {
                kr = nxSources-1;
                kl = nxSources-1;
            }
            // out of left boundary
            if(kl<=0)
            {
                kr = 0;
                kl = 0;
            }
            // out of top boundary
            if(it>=nySources)
            {
                it = nySources-1;
                ib = nySources-1;
            }
            // out of bottom boundary
            if(ib<=0)
            {
                it = 0;
                ib = 0;
            }

            // set left and right weights
            if(kl!=kr)
            {
                PetscReal xr = xsSources + kr*dxSources,
                          xl = xsSources + kl*dxSources;
                wl = (xr - xQuery) / dxSources;
                wr = (xQuery - xl) / dxSources;
            }
            else
            {
                wl = 0.5;
                wr = 0.5;
            }

            // set top and bottom weights
            if(it!=ib)
            {
                PetscReal yt = ysSources + it*dySources,
                          yb = ysSources + ib*dySources;
                wb = (yt - yQuery) / dySources;
                wt = (yQuery - yb) / dySources;
            }
            else
            {
                wt = 0.5;
                wb = 0.5;
            }

            // interpolate vertical mesh deformation
            PetscReal etaTop = wl*etaSources[kl][it] + wr*etaSources[kr][it],
                      etaBot = wl*etaSources[kl][ib] + wr*etaSources[kr][ib],
                      etaPt  = wt*etaTop + wb*etaBot;

            // distribute the vertical mesh deformation
            for (j=ys; j<ye; j++)
            {

                PetscReal fraction  = (coor[k][j][i].z - mesh->bounds.zmin)/mesh->bounds.Lz;
                coor[k][j][i].z     =  coor[k][j][i].z + fraction*etaPt;
            }
        }
    }

    DMDAVecRestoreArray(mesh->fda, gCoor, &coor);

    // set global coordinates into da
    DMSetCoordinates(mesh->da, gCoor);

    // clean local vectors
    for (i=0; i<nxSources; i++)
    {
        std::vector<PetscReal> ().swap(etaSources[i]);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetMeshMetrics(mesh_ *mesh)
{
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs  = info.xs, xe = info.xs + info.xm;
    PetscInt      ys  = info.ys, ye = info.ys + info.ym;
    PetscInt      zs  = info.zs, ze = info.zs + info.zm;
    PetscInt      gxs = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs = info.gzs, gze = info.gzs + info.gzm;
    PetscInt      mx  = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    Vec           Csi, Eta, Zet;
    Vec           ICsi, IEta, IZet;
    Vec           JCsi, JEta, JZet;
    Vec           KCsi, KEta, KZet;
    Vec           Aj;
    Vec           IAj, JAj, KAj;

    Vec           lCoor;

    Vec           Centx, Centy, Centz;

    Cmpnts        ***coor;
    Cmpnts        ***cent, ***centx, ***centy, ***centz;
    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;
    PetscReal     ***aj;
    PetscReal     ***iaj, ***jaj, ***kaj;

    PetscReal     xcp, ycp, zcp, xcm, ycm, zcm;

    // indices for internal cells
    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // create global vectors
    VecDuplicate(mesh->lCent,  &Centx);
    VecDuplicate(mesh->lCent,  &Centy);
    VecDuplicate(mesh->lCent,  &Centz);

    VecDuplicate(mesh->Cent,  &(Csi));
    VecDuplicate(mesh->Cent,  &(Eta));
    VecDuplicate(mesh->Cent,  &(Zet));

    VecDuplicate(mesh->Cent,  &(ICsi));
    VecDuplicate(mesh->Cent,  &(IEta));
    VecDuplicate(mesh->Cent,  &(IZet));
    VecDuplicate(mesh->Cent,  &(JCsi));
    VecDuplicate(mesh->Cent,  &(JEta));
    VecDuplicate(mesh->Cent,  &(JZet));
    VecDuplicate(mesh->Cent,  &(KCsi));
    VecDuplicate(mesh->Cent,  &(KEta));
    VecDuplicate(mesh->Cent,  &(KZet));

    VecDuplicate(mesh->Nvert, &(Aj));
    VecDuplicate(mesh->Nvert, &(IAj));
    VecDuplicate(mesh->Nvert, &(JAj));
    VecDuplicate(mesh->Nvert, &(KAj));

    // get global empty metrics
    DMDAVecGetArray(fda, Csi, &csi);
    DMDAVecGetArray(fda, Eta, &eta);
    DMDAVecGetArray(fda, Zet, &zet);
    DMDAVecGetArray(da,  Aj,  &aj);

    // get mesh coordinates
    DMGetCoordinatesLocal(da, &lCoor);
    DMDAVecGetArray(fda, lCoor, &coor);

    // ------------------------------------------------------------------------ //
    // create cell centers
    // ------------------------------------------------------------------------ //

    DMDAVecGetArray(fda, mesh->Cent, &cent);

    // loop over internal cells
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                cent[k][j][i].x = 0.125 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                    coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x +
                    coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x +
                    coor[k-1][j  ][i-1].x + coor[k-1][j-1][i-1].x
                );
                cent[k][j][i].y = 0.125 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                    coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y +
                    coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y +
                    coor[k-1][j  ][i-1].y + coor[k-1][j-1][i-1].y
                );
                cent[k][j][i].z = 0.125 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                    coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z +
                    coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z +
                    coor[k-1][j  ][i-1].z + coor[k-1][j-1][i-1].z
                );

            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->Cent, &cent);

    // scatter the cell centers coords
    DMGlobalToLocalBegin(fda, mesh->Cent, INSERT_VALUES, mesh->lCent);
    DMGlobalToLocalEnd  (fda, mesh->Cent, INSERT_VALUES, mesh->lCent);

    // ------------------------------------------------------------------------ //
    // create cell centered contrav. basis and transformation jacobian
    // ------------------------------------------------------------------------ //

    // create dummy contrav. basis to be interpolated at cell centers
    Vec     lCsitmp,   lEtatmp,   lZettmp;
    Cmpnts  ***csitmp, ***etatmp, ***zettmp;

    DMCreateLocalVector(fda, &lCsitmp);
    DMCreateLocalVector(fda, &lEtatmp);
    DMCreateLocalVector(fda, &lZettmp);

    DMDAVecGetArray(fda, lCsitmp, &csitmp);
    DMDAVecGetArray(fda, lEtatmp, &etatmp);
    DMDAVecGetArray(fda, lZettmp, &zettmp);

    // when computing derivatives, cell index increment is assumed to
    // be the change in the curvilinear coord in that direction, hence
    // the dcsi, deta, dzet is always 1. Cell indices actually are the
    // curvilinear coordinates.

    // calculating i direction metrics: csi
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                // csi = de X dz
                PetscReal dxde = 0.5 *
                (
                    coor[k][j  ][i].x + coor[k-1][j  ][i].x -
                    coor[k][j-1][i].x - coor[k-1][j-1][i].x
                );
                PetscReal dyde = 0.5 *
                (
                    coor[k][j  ][i].y + coor[k-1][j  ][i].y -
                    coor[k][j-1][i].y - coor[k-1][j-1][i].y
                );
                PetscReal dzde = 0.5 *
                (
                    coor[k][j  ][i].z + coor[k-1][j  ][i].z -
                    coor[k][j-1][i].z - coor[k-1][j-1][i].z
                );

                PetscReal dxdz = 0.5 *
                (
                    coor[k  ][j-1][i].x + coor[k  ][j][i].x -
                    coor[k-1][j-1][i].x - coor[k-1][j][i].x
                );
                PetscReal dydz = 0.5 *
                (
                    coor[k  ][j-1][i].y + coor[k  ][j][i].y -
                    coor[k-1][j-1][i].y - coor[k-1][j][i].y
                );
                PetscReal dzdz = 0.5 *
                (
                    coor[k  ][j-1][i].z + coor[k  ][j][i].z -
                    coor[k-1][j-1][i].z - coor[k-1][j][i].z
                );

                csitmp[k][j][i].x = dyde * dzdz - dzde * dydz;
                csitmp[k][j][i].y =-dxde * dzdz + dzde * dxdz;
                csitmp[k][j][i].z = dxde * dydz - dyde * dxdz;
            }
        }
    }

    // calculating j direction metrics: eta
    for (k=lzs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // eta = dz X dc
                PetscReal dxdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
                    coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x
                );
                PetscReal dydc = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
                    coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y
                );
                PetscReal dzdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
                    coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z
                );
                PetscReal dxdz = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                    coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x
                );
                PetscReal dydz = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                    coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y
                );
                PetscReal dzdz = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                    coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z
                );

                etatmp[k][j][i].x = dydz * dzdc - dzdz * dydc;
                etatmp[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
                etatmp[k][j][i].z = dxdz * dydc - dydz * dxdc;
            }
        }
    }

    // calculating k direction metrics: zet
    for (k=zs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                // zet = dc X de
                PetscReal dxdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
                    coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x
                );
                PetscReal dydc = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
                    coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y
                );
                PetscReal dzdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
                    coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z
                );
                PetscReal dxde = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                    coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x
                );
                PetscReal dyde = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                    coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y
                );
                PetscReal dzde = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                    coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z
                );

                zettmp[k][j][i].x = dydc * dzde - dzdc * dyde;
                zettmp[k][j][i].y =-dxdc * dzde + dzdc * dxde;
                zettmp[k][j][i].z = dxdc * dyde - dydc * dxde;
            }
        }
    }

    DMDAVecRestoreArray(fda, lCsitmp, &csitmp);
    DMDAVecRestoreArray(fda, lEtatmp, &etatmp);
    DMDAVecRestoreArray(fda, lZettmp, &zettmp);

    // scatter dummy basis
    DMLocalToLocalBegin(fda, lCsitmp, INSERT_VALUES, lCsitmp);
    DMLocalToLocalEnd  (fda, lCsitmp, INSERT_VALUES, lCsitmp);
    DMLocalToLocalBegin(fda, lEtatmp, INSERT_VALUES, lEtatmp);
    DMLocalToLocalEnd  (fda, lEtatmp, INSERT_VALUES, lEtatmp);
    DMLocalToLocalBegin(fda, lZettmp, INSERT_VALUES, lZettmp);
    DMLocalToLocalEnd  (fda, lZettmp, INSERT_VALUES, lZettmp);

    DMDAVecGetArray(fda, lCsitmp, &csitmp);
    DMDAVecGetArray(fda, lEtatmp, &etatmp);
    DMDAVecGetArray(fda, lZettmp, &zettmp);

    // interpolate at cell centers from dummy basis
    // ------------------------------------
    // csi = 0.5 * ( csi_i + csi_i-1 )
    // eta = 0.5 * ( eta_j + eta_j-1 )
    // zet = 0.5 * ( zet_k + zet_k-1 )
    // ------------------------------------
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                AxByC
                (
                    0.5, csitmp[k][j][i],
                    0.5, csitmp[k][j][i-1],
                    &csi[k][j][i]
                );
                AxByC
                (
                    0.5, etatmp[k][j][i],
                    0.5, etatmp[k][j-1][i],
                    &eta[k][j][i]
                );
                AxByC
                (
                    0.5, zettmp[k][j][i],
                    0.5, zettmp[k-1][j][i],
                    &zet[k][j][i]
                );
            }
        }
    }

    DMDAVecRestoreArray(fda, lCsitmp, &csitmp);
    DMDAVecRestoreArray(fda, lEtatmp, &etatmp);
    DMDAVecRestoreArray(fda, lZettmp, &zettmp);

    // destroy dummy basis
    VecDestroy(&lCsitmp);
    VecDestroy(&lEtatmp);
    VecDestroy(&lZettmp);

    // create the Jacobian of the transformation
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if (i>0 && j>0 && k>0)
                {
                    PetscReal dxdc = 0.25 *
                    (
                        coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                        coor[k-1][j  ][i  ].x + coor[k-1][j-1][i  ].x -
                        coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x -
                        coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x
                    );
                    PetscReal dydc = 0.25 *
                    (
                        coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                        coor[k-1][j  ][i  ].y + coor[k-1][j-1][i  ].y -
                        coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y -
                        coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y
                    );
                    PetscReal dzdc = 0.25 *
                    (
                        coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                        coor[k-1][j  ][i  ].z + coor[k-1][j-1][i  ].z -
                        coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z -
                        coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z
                    );
                    PetscReal dxde = 0.25 *
                    (
                        coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x +
                        coor[k-1][j  ][i  ].x + coor[k-1][j  ][i-1].x -
                        coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x -
                        coor[k-1][j-1][i  ].x - coor[k-1][j-1][i-1].x
                    );
                    PetscReal dyde = 0.25 *
                    (
                        coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y +
                        coor[k-1][j  ][i  ].y + coor[k-1][j  ][i-1].y -
                        coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y -
                        coor[k-1][j-1][i  ].y - coor[k-1][j-1][i-1].y
                    );
                    PetscReal dzde = 0.25 *
                    (
                        coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z +
                        coor[k-1][j  ][i  ].z + coor[k-1][j  ][i-1].z -
                        coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z -
                        coor[k-1][j-1][i  ].z - coor[k-1][j-1][i-1].z
                    );
                    PetscReal dxdz = 0.25 *
                    (
                        coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x +
                        coor[k  ][j  ][i-1].x + coor[k  ][j-1][i-1].x -
                        coor[k-1][j  ][i  ].x - coor[k-1][j-1][i  ].x -
                        coor[k-1][j  ][i-1].x - coor[k-1][j-1][i-1].x
                    );
                    PetscReal dydz = 0.25 *
                    (
                        coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y +
                        coor[k  ][j  ][i-1].y + coor[k  ][j-1][i-1].y -
                        coor[k-1][j  ][i  ].y - coor[k-1][j-1][i  ].y -
                        coor[k-1][j  ][i-1].y - coor[k-1][j-1][i-1].y
                    );
                    PetscReal dzdz = 0.25 *
                    (
                        coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z +
                        coor[k  ][j  ][i-1].z + coor[k  ][j-1][i-1].z -
                        coor[k-1][j  ][i  ].z - coor[k-1][j-1][i  ].z -
                        coor[k-1][j  ][i-1].z - coor[k-1][j-1][i-1].z
                    );

                    // rec. of the det. of the transformation matrix
                    aj[k][j][i] = dxdc * (dyde * dzdz - dzde * dydz) -
                                  dydc * (dxde * dzdz - dzde * dxdz) +
                                  dzdc * (dxde * dydz - dyde * dxdz);

                    aj[k][j][i] = 1./aj[k][j][i];
                }
            }
        }
    }

    // handle boundaries: mirror grid outside the boundaries
    if (xs==0)
    {
        i = xs;
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                csi[k][j][i] = csi[k][j][i+1];
                eta[k][j][i] = eta[k][j][i+1];
                zet[k][j][i] = zet[k][j][i+1];
                aj[k][j][i]  = aj[k][j][i+1];
            }
        }
    }
    if (xe==mx)
    {
        i = xe-1;
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                csi[k][j][i] = csi[k][j][i-1];
                eta[k][j][i] = eta[k][j][i-1];
                zet[k][j][i] = zet[k][j][i-1];
                aj[k][j][i]  = aj[k][j][i-1];
            }
        }
    }

    if (ys==0)
    {
        j = ys;
        for (k=zs; k<ze; k++)
        {
            for (i=xs; i<xe; i++)
            {
                eta[k][j][i] = eta[k][j+1][i];
                csi[k][j][i] = csi[k][j+1][i];
                zet[k][j][i] = zet[k][j+1][i];
                aj[k][j][i]  = aj[k][j+1][i];
            }
        }
    }

    if (ye==my)
    {
        j = ye-1;
        for (k=zs; k<ze; k++)
        {
            for (i=xs; i<xe; i++)
            {
                eta[k][j][i] = eta[k][j-1][i];
                csi[k][j][i] = csi[k][j-1][i];
                zet[k][j][i] = zet[k][j-1][i];
                aj[k][j][i]  = aj[k][j-1][i];
            }
        }
    }

    if (zs==0)
    {
        k = zs;
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                zet[k][j][i] = zet[k+1][j][i];
                eta[k][j][i] = eta[k+1][j][i];
                csi[k][j][i] = csi[k+1][j][i];
                aj[k][j][i]  = aj[k+1][j][i];
            }
        }
    }

    if (ze==mz)
    {
        k = ze-1;
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                zet[k][j][i] = zet[k-1][j][i];
                eta[k][j][i] = eta[k-1][j][i];
                csi[k][j][i] = csi[k-1][j][i];
                aj[k][j][i]  = aj[k-1][j][i];

            }
        }
    }

    DMDAVecRestoreArray(fda, Csi,   &csi);
    DMDAVecRestoreArray(fda, Eta,   &eta);
    DMDAVecRestoreArray(fda, Zet,   &zet);
    DMDAVecRestoreArray(da,  Aj,    &aj);

    // scatter csi, eta, zet, aj from global to local
    DMGlobalToLocalBegin(fda, Csi, INSERT_VALUES, mesh->lCsi);
    DMGlobalToLocalEnd  (fda, Csi, INSERT_VALUES, mesh->lCsi);
    DMGlobalToLocalBegin(fda, Eta, INSERT_VALUES, mesh->lEta);
    DMGlobalToLocalEnd  (fda, Eta, INSERT_VALUES, mesh->lEta);
    DMGlobalToLocalBegin(fda, Zet, INSERT_VALUES, mesh->lZet);
    DMGlobalToLocalEnd  (fda, Zet, INSERT_VALUES, mesh->lZet);
    DMGlobalToLocalBegin(da,  Aj,  INSERT_VALUES, mesh->lAj);
    DMGlobalToLocalEnd  (da,  Aj,  INSERT_VALUES, mesh->lAj);

    Cmpnts       ***lcsi, ***leta, ***lzet;
    PetscScalar  ***laj;

    // handle periodicity:
    // i,j,k    : ghost cell are the internal cells next to opp. boundary
    // ii,jj,kk : ghost cells are 2 external cells next to current boundary

    DMDAVecGetArray(fda, Csi, &csi);
    DMDAVecGetArray(fda, Eta, &eta);
    DMDAVecGetArray(fda, Zet, &zet);
    DMDAVecGetArray(da,  Aj,  &aj);

    DMDAVecGetArray(fda, mesh->lCsi, &lcsi);
    DMDAVecGetArray(fda, mesh->lEta, &leta);
    DMDAVecGetArray(fda, mesh->lZet, &lzet);
    DMDAVecGetArray(da,  mesh->lAj,  &laj);

    // loop over all cells (internal and physical ghosts)
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                PetscInt   flag = 0, i_flag = 0, j_flag = 0, k_flag = 0;
                PetscInt   a = i, b = j, c = k;

                if(mesh->i_periodic && i==0)
                {
                    a      = mx-2;
                    i_flag = 1;
                }
                else if(mesh->i_periodic && i==mx-1)
                {
                    a      = 1;
                    i_flag = 1;
                }

                if(mesh->j_periodic && j==0)
                {
                    b      = my-2;
                    j_flag = 1;
                }
                else if(mesh->j_periodic && j==my-1)
                {
                    b      = 1;
                    j_flag = 1;
                }

                if(mesh->k_periodic && k==0)
                {
                    c      = mz-2;
                    k_flag = 1;
                }
                else if(mesh->k_periodic && k==mz-1)
                {
                    c      = 1;
                    k_flag = 1;
                }

                if(mesh->ii_periodic && i==0)
                {
                    a      = -2;
                    i_flag = 1;
                }
                else if(mesh->ii_periodic && i==mx-1)
                {
                    a      = mx+1;
                    i_flag = 1;
                }

                if(mesh->jj_periodic && j==0)
                {
                    b      = -2;
                    j_flag = 1;
                }
                else if(mesh->jj_periodic && j==my-1)
                {
                    b      = my+1;
                    j_flag = 1;
                }

                if(mesh->kk_periodic && k==0)
                {
                    c      = -2;
                    k_flag = 1;
                }
                else if(mesh->kk_periodic && k==mz-1)
                {
                    c      = mz+1;
                    k_flag = 1;
                }

                flag = i_flag + j_flag + k_flag;

                if(flag)
                {
                    // modify local variables
                    lcsi[k][j][i] = lcsi[c][b][a];
                    leta[k][j][i] = leta[c][b][a];
                    lzet[k][j][i] = lzet[c][b][a];
                    laj[k][j][i]  = laj[c][b][a];

                    // modify global variables
                    csi[k][j][i]  = lcsi[k][j][i];
                    eta[k][j][i]  = leta[k][j][i];
                    zet[k][j][i]  = lzet[k][j][i];
                    aj[k][j][i]   = laj[k][j][i];

                }

            }
        }
    }

    DMDAVecRestoreArray(fda, Csi, &csi);
    DMDAVecRestoreArray(fda, Eta, &eta);
    DMDAVecRestoreArray(fda, Zet, &zet);
    DMDAVecRestoreArray(da,  Aj,  &aj);

    DMDAVecRestoreArray(fda, mesh->lCsi, &lcsi);
    DMDAVecRestoreArray(fda, mesh->lEta, &leta);
    DMDAVecRestoreArray(fda, mesh->lZet, &lzet);
    DMDAVecRestoreArray(da,  mesh->lAj,  &laj);

    // global to local scatter
    DMGlobalToLocalBegin(fda, Csi, INSERT_VALUES, mesh->lCsi);
    DMGlobalToLocalEnd  (fda, Csi, INSERT_VALUES, mesh->lCsi);
    DMGlobalToLocalBegin(fda, Eta, INSERT_VALUES, mesh->lEta);
    DMGlobalToLocalEnd  (fda, Eta, INSERT_VALUES, mesh->lEta);
    DMGlobalToLocalBegin(fda, Zet, INSERT_VALUES, mesh->lZet);
    DMGlobalToLocalEnd  (fda, Zet, INSERT_VALUES, mesh->lZet);
    DMGlobalToLocalBegin(da,  Aj,  INSERT_VALUES, mesh->lAj);
    DMGlobalToLocalEnd  (da,  Aj,  INSERT_VALUES, mesh->lAj);

    MPI_Barrier(mesh->MESH_COMM);

    //find the ghost node cell centers excluding the corner cells - only for cartesian mesh
    ghostnodesCellcenter(mesh);

    DMDAVecGetArray(fda, mesh->lCsi, &lcsi);
    DMDAVecGetArray(fda, mesh->lEta, &leta);
    DMDAVecGetArray(fda, mesh->lZet, &lzet);
    DMDAVecGetArray(da,  mesh->lAj,  &laj);

    // perform mesh checks
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                // exclude control at corners
                if(!isOnCornerCellCenters(i, j, k, info))
                {
                    if(nMag(lcsi[k][j][i]) < 1e-20)
                    {
                       char error[512];
                        sprintf(error, "negative or close-to-zero i-face area detected at i-face with [i][j][k] = [%ld][%ld][%ld], value = %lf\n", i, j, k, nMag(lcsi[k][j][i]));
                        fatalErrorInFunction("SetMeshMetrics",  error);
                    }

                    if(nMag(leta[k][j][i]) < 1e-20)
                    {
                       char error[512];
                        sprintf(error, "negative or close-to-zero j-face area detected at j-face with [i][j][k] = [%ld][%ld][%ld], value = %lf\n", i, j, k, nMag(leta[k][j][i]));
                        fatalErrorInFunction("SetMeshMetrics",  error);
                    }

                    if(nMag(lzet[k][j][i]) < 1e-20)
                    {
                       char error[512];
                        sprintf(error, "negative or close-to-zero k-face area detected at k-face with [i][j][k] = [%ld][%ld][%ld], value = %lf\n", i, j, k, nMag(lzet[k][j][i]));
                        fatalErrorInFunction("SetMeshMetrics",  error);
                    }

                    if(laj[k][j][i] < 1e-20)
                    {
                       char error[512];
                        sprintf(error, "negative or close-to-zero cell volume detected at cell[%ld][%ld][%ld], value = %lf\n", k, j, i, laj[k][j][i]);
                        fatalErrorInFunction("SetMeshMetrics",  error);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi, &lcsi);
    DMDAVecRestoreArray(fda, mesh->lEta, &leta);
    DMDAVecRestoreArray(fda, mesh->lZet, &lzet);
    DMDAVecRestoreArray(da,  mesh->lAj,  &laj);

    // ---------------------------------------------------------------------- //
    // create i-face centered contrav. basis and transformation jacobian
    // ---------------------------------------------------------------------- //

    DMDAVecGetArray(fda, ICsi, &icsi);
    DMDAVecGetArray(fda, IEta, &ieta);
    DMDAVecGetArray(fda, IZet, &izet);
    DMDAVecGetArray(da,  IAj,  &iaj);

    DMDAVecGetArray(fda, Centx, &centx);

    for(k=gzs+1; k<gze; k++)
    {
        for (j=gys+1; j<gye; j++)
        {
            for (i=gxs; i<gxe; i++)
            {
                centx[k][j][i].x = 0.25 *
                (
                    coor[k  ][j  ][i].x + coor[k-1][j  ][i].x +
                    coor[k  ][j-1][i].x + coor[k-1][j-1][i].x
                );
                centx[k][j][i].y = 0.25 *
                (
                    coor[k  ][j  ][i].y + coor[k-1][j  ][i].y +
                    coor[k  ][j-1][i].y + coor[k-1][j-1][i].y
                );
                centx[k][j][i].z = 0.25 *
                (
                    coor[k  ][j  ][i].z + coor[k-1][j  ][i].z +
                    coor[k  ][j-1][i].z + coor[k-1][j-1][i].z
                );
            }
        }
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                PetscScalar dxdc, dydc, dzdc,
                            dxde, dyde, dzde,
                            dxdz, dydz, dzdz;

                // d/dc - - - - - - - - - - - - - - - - - - - - - - - - - - - -

                // left boundary
                if (i==0)
                {
                    // forward differences
                    dxdc = centx[k][j][i+1].x - centx[k][j][i].x;
                    dydc = centx[k][j][i+1].y - centx[k][j][i].y;
                    dzdc = centx[k][j][i+1].z - centx[k][j][i].z;
                }

                // right boundary
                else if (i==mx-2)
                {
                    // non-ghosted periodicity
                    if(mesh->i_periodic)
                    {
                        // average derivatives at 0 and 1st right internal face
                        dxdc = 0.5 *
                        (
                            centx[k][j][1].x - centx[k][j][0].x +
                            centx[k][j][i].x - centx[k][j][i-1].x
                        );
                        dydc = 0.5 *
                        (
                            centx[k][j][1].y - centx[k][j][0].y +
                            centx[k][j][i].y - centx[k][j][i-1].y
                        );
                        dzdc = 0.5 *
                        (
                            centx[k][j][1].z - centx[k][j][0].z +
                            centx[k][j][i].z - centx[k][j][i-1].z
                        );
                    }
                    // ghosted periodicity
                    else if(mesh->ii_periodic)
                    {
                        // average derivatives at penultimate right ghost face and
                        // 1st right internal face
                        dxdc = 0.5 *
                        (
                            centx[k][j][mx+1].x - centx[k][j][mx+0].x +
                            centx[k][j][i].x - centx[k][j][i-1].x
                        );
                        dydc = 0.5 *
                        (
                            centx[k][j][mx+1].y - centx[k][j][mx+0].y +
                            centx[k][j][i].y - centx[k][j][i-1].y
                        );
                        dzdc = 0.5 *
                        (
                            centx[k][j][mx+1].z - centx[k][j][mx+0].z +
                            centx[k][j][i].z - centx[k][j][i-1].z
                        );
                    }
                    // backward differences
                    else
                    {
                        dxdc = centx[k][j][i].x - centx[k][j][i-1].x;
                        dydc = centx[k][j][i].y - centx[k][j][i-1].y;
                        dzdc = centx[k][j][i].z - centx[k][j][i-1].z;
                    }
                }
                // internal faces
                else
                {
                    // central differences
                    dxdc = 0.5 *
                    (
                        centx[k][j][i+1].x - centx[k][j][i-1].x
                    );
                    dydc = 0.5 *
                    (
                        centx[k][j][i+1].y - centx[k][j][i-1].y
                    );
                    dzdc = 0.5 *
                    (
                        centx[k][j][i+1].z - centx[k][j][i-1].z
                    );
                }

                // seba: use mesh points here instead of face centers for j and k directions
                dxde = 0.5 *
                (
                    coor[k][j  ][i].x + coor[k-1][j  ][i].x -
                    coor[k][j-1][i].x - coor[k-1][j-1][i].x
                );
                dyde = 0.5 *
                (
                    coor[k][j  ][i].y + coor[k-1][j  ][i].y -
                    coor[k][j-1][i].y - coor[k-1][j-1][i].y
                );
                dzde = 0.5 *
                (
                    coor[k][j  ][i].z + coor[k-1][j  ][i].z -
                    coor[k][j-1][i].z - coor[k-1][j-1][i].z
                );

                dxdz = 0.5 *
                (
                    coor[k  ][j-1][i].x + coor[k  ][j][i].x -
                    coor[k-1][j-1][i].x - coor[k-1][j][i].x
                );
                dydz = 0.5 *
                (
                    coor[k  ][j-1][i].y + coor[k  ][j][i].y -
                    coor[k-1][j-1][i].y - coor[k-1][j][i].y
                );
                dzdz = 0.5 *
                (
                    coor[k  ][j-1][i].z + coor[k  ][j][i].z -
                    coor[k-1][j-1][i].z - coor[k-1][j][i].z
                );

                // csi = eta X zet
                icsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                icsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
                icsi[k][j][i].z = dxde * dydz - dyde * dxdz;
                // eta = zet X csi
                ieta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                ieta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
                ieta[k][j][i].z = dxdz * dydc - dydz * dxdc;
                // zet = csi X eta
                izet[k][j][i].x = dydc * dzde - dzdc * dyde;
                izet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
                izet[k][j][i].z = dxdc * dyde - dydc * dxde;

                iaj[k][j][i]
                =
                dxdc * (dyde * dzdz - dzde * dydz) -
                dydc * (dxde * dzdz - dzde * dxdz) +
                dzdc * (dxde * dydz - dyde * dxdz);

                iaj[k][j][i] = 1./iaj[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(fda, ICsi, &icsi);
    DMDAVecRestoreArray(fda, IEta, &ieta);
    DMDAVecRestoreArray(fda, IZet, &izet);
    DMDAVecRestoreArray(da,  IAj,  &iaj);
    DMDAVecRestoreArray(fda, Centx, &centx);

    // ---------------------------------------------------------------------- //
    // create j-face centered contrav. basis and transformation jacobian
    // ---------------------------------------------------------------------- //

    DMDAVecGetArray(fda, JCsi, &jcsi);
    DMDAVecGetArray(fda, JEta, &jeta);
    DMDAVecGetArray(fda, JZet, &jzet);
    DMDAVecGetArray(da,  JAj,  &jaj);
    DMDAVecGetArray(fda, Centy, &centy);

    for(k=gzs+1; k<gze; k++)
    {
        for (j=gys; j<gye; j++)
        {
            for (i=gxs+1; i<gxe; i++)
            {
                centy[k][j][i].x = 0.25 *
                (
                    coor[k  ][j][i  ].x + coor[k-1][j][i  ].x +
                    coor[k  ][j][i-1].x + coor[k-1][j][i-1].x
                );
                centy[k][j][i].y = 0.25 *
                (
                    coor[k  ][j][i  ].y + coor[k-1][j][i  ].y +
                    coor[k  ][j][i-1].y + coor[k-1][j][i-1].y
                );
                centy[k][j][i].z = 0.25 *
                (
                    coor[k  ][j][i  ].z + coor[k-1][j][i  ].z +
                    coor[k  ][j][i-1].z + coor[k-1][j][i-1].z
                );
            }
        }
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscScalar dxdc, dydc, dzdc,
                            dxde, dyde, dzde,
                            dxdz, dydz, dzdz;

                // seba: use mesh points here instead of face centers for i and k directions
                dxdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k-1][j  ][i  ].x -
                    coor[k  ][j  ][i-1].x - coor[k-1][j  ][i-1].x
                );
                dydc = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k-1][j  ][i  ].y -
                    coor[k  ][j  ][i-1].y - coor[k-1][j  ][i-1].y
                );
                dzdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k-1][j  ][i  ].z -
                    coor[k  ][j  ][i-1].z - coor[k-1][j  ][i-1].z
                );

                if (j==0)
                {
                    dxde = centy[k][j+1][i].x - centy[k][j][i].x;
                    dyde = centy[k][j+1][i].y - centy[k][j][i].y;
                    dzde = centy[k][j+1][i].z - centy[k][j][i].z;
                }
                else if (j==my-2)
                {
                    if(mesh->j_periodic)
                    {
                        dxde = 0.5 *
                        (
                            centy[k][1][i].x - centy[k][0][i].x +
                            centy[k][j][i].x - centy[k][j-1][i].x
                        );
                        dyde = 0.5 *
                        (
                            centy[k][1][i].y - centy[k][0][i].y +
                            centy[k][j][i].y - centy[k][j-1][i].y
                        );
                        dzde = 0.5 *
                        (
                            centy[k][1][i].z - centy[k][0][i].z +
                            centy[k][j][i].z - centy[k][j-1][i].z
                        );
                    }
                    else if(mesh->jj_periodic)
                    {
                        dxde = 0.5 *
                        (
                            centy[k][my+1][i].x - centy[k][my+0][i].x +
                            centy[k][j][i].x - centy[k][j-1][i].x
                        );
                        dyde = 0.5 *
                        (
                            centy[k][my+1][i].y - centy[k][my+0][i].y +
                            centy[k][j][i].y - centy[k][j-1][i].y
                        );
                        dzde = 0.5 *
                        (
                            centy[k][my+1][i].z - centy[k][my+0][i].z +
                            centy[k][j][i].z - centy[k][j-1][i].z
                        );
                    }
                    else
                    {
                        dxde = centy[k][j][i].x - centy[k][j-1][i].x;
                        dyde = centy[k][j][i].y - centy[k][j-1][i].y;
                        dzde = centy[k][j][i].z - centy[k][j-1][i].z;
                    }
                }
                else
                {
                    dxde = 0.5 * (centy[k][j+1][i].x - centy[k][j-1][i].x);
                    dyde = 0.5 * (centy[k][j+1][i].y - centy[k][j-1][i].y);
                    dzde = 0.5 * (centy[k][j+1][i].z - centy[k][j-1][i].z);
                }

                dxdz = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                    coor[k-1][j  ][i  ].x - coor[k-1][j  ][i-1].x
                );
                dydz = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                    coor[k-1][j  ][i  ].y - coor[k-1][j  ][i-1].y
                );
                dzdz = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                    coor[k-1][j  ][i  ].z - coor[k-1][j  ][i-1].z
                );

                jcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                jcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
                jcsi[k][j][i].z = dxde * dydz - dyde * dxdz;

                jeta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                jeta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
                jeta[k][j][i].z = dxdz * dydc - dydz * dxdc;

                jzet[k][j][i].x = dydc * dzde - dzdc * dyde;
                jzet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
                jzet[k][j][i].z = dxdc * dyde - dydc * dxde;


                jaj[k][j][i]
                =
                dxdc * (dyde * dzdz - dzde * dydz) -
                dydc * (dxde * dzdz - dzde * dxdz) +
                dzdc * (dxde * dydz - dyde * dxdz);

                jaj[k][j][i] = 1./jaj[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(fda, JCsi, &jcsi);
    DMDAVecRestoreArray(fda, JEta, &jeta);
    DMDAVecRestoreArray(fda, JZet, &jzet);
    DMDAVecRestoreArray(da,  JAj,  &jaj);
    DMDAVecRestoreArray(fda, Centy, &centy);

    // ---------------------------------------------------------------------- //
    // create k-face centered contrav. basis and transformation jacobian
    // ---------------------------------------------------------------------- //

    DMDAVecGetArray(fda, KCsi, &kcsi);
    DMDAVecGetArray(fda, KEta, &keta);
    DMDAVecGetArray(fda, KZet, &kzet);
    DMDAVecGetArray(da,  KAj,  &kaj);
    DMDAVecGetArray(fda, Centz, &centz);

    for(k=gzs; k<gze; k++)
    {
        for (j=gys+1; j<gye; j++)
        {
            for (i=gxs+1; i<gxe; i++)
            {
                centz[k][j][i].x = 0.25 *
                (
                    coor[k  ][j][i  ].x + coor[k][j-1][i  ].x +
                    coor[k  ][j][i-1].x + coor[k][j-1][i-1].x
                );
                centz[k][j][i].y = 0.25 *
                (
                    coor[k  ][j][i  ].y + coor[k][j-1][i  ].y +
                    coor[k  ][j][i-1].y + coor[k][j-1][i-1].y
                );
                centz[k][j][i].z = 0.25 *
                (
                    coor[k  ][j][i  ].z + coor[k][j-1][i  ].z +
                    coor[k  ][j][i-1].z + coor[k][j-1][i-1].z
                );
            }
        }
    }

    for (k=zs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscScalar dxdc, dydc, dzdc,
                            dxde, dyde, dzde,
                            dxdz, dydz, dzdz;

                // seba: use mesh points here instead of face centers for i and j directions
                dxdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j-1][i  ].x -
                    coor[k  ][j  ][i-1].x - coor[k  ][j-1][i-1].x
                );
                dydc = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j-1][i  ].y -
                    coor[k  ][j  ][i-1].y - coor[k  ][j-1][i-1].y
                );
                dzdc = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j-1][i  ].z -
                    coor[k  ][j  ][i-1].z - coor[k  ][j-1][i-1].z
                );

                dxde = 0.5 *
                (
                    coor[k  ][j  ][i  ].x + coor[k  ][j  ][i-1].x -
                    coor[k  ][j-1][i  ].x - coor[k  ][j-1][i-1].x
                );
                dyde = 0.5 *
                (
                    coor[k  ][j  ][i  ].y + coor[k  ][j  ][i-1].y -
                    coor[k  ][j-1][i  ].y - coor[k  ][j-1][i-1].y
                );
                dzde = 0.5 *
                (
                    coor[k  ][j  ][i  ].z + coor[k  ][j  ][i-1].z -
                    coor[k  ][j-1][i  ].z - coor[k  ][j-1][i-1].z
                );

                if (k==0)
                {
                    dxdz = (centz[k+1][j][i].x - centz[k][j][i].x);
                    dydz = (centz[k+1][j][i].y - centz[k][j][i].y);
                    dzdz = (centz[k+1][j][i].z - centz[k][j][i].z);
                }
                else if (k==mz-2)
                {
                    if(mesh->k_periodic)
                    {
                        dxdz = 0.5 *
                        (
                            centz[1][j][i].x - centz[0][j][i].x +
                            centz[k][j][i].x - centz[k-1][j][i].x
                        );
                        dydz = 0.5 *
                        (
                            centz[1][j][i].y - centz[0][j][i].y +
                            centz[k][j][i].y - centz[k-1][j][i].y
                        );
                        dzdz = 0.5 *
                        (
                            centz[1][j][i].z - centz[0][j][i].z +
                            centz[k][j][i].z - centz[k-1][j][i].z
                        );
                    }
                    else if(mesh->kk_periodic)
                    {
                        dxdz = 0.5 *
                        (
                            centz[mz+1][j][i].x - centz[mz+0][j][i].x +
                            centz[k][j][i].x - centz[k-1][j][i].x
                        );
                        dydz = 0.5 *
                        (
                            centz[mz+1][j][i].y - centz[mz+0][j][i].y +
                            centz[k][j][i].y - centz[k-1][j][i].y
                        );
                        dzdz = 0.5 *
                        (
                            centz[mz+1][j][i].z - centz[mz+0][j][i].z +
                            centz[k][j][i].z - centz[k-1][j][i].z
                        );
                    }
                    else
                    {
                    dxdz = (centz[k][j][i].x - centz[k-1][j][i].x);
                    dydz = (centz[k][j][i].y - centz[k-1][j][i].y);
                    dzdz = (centz[k][j][i].z - centz[k-1][j][i].z);
                    }
                }
                else
                {
                    dxdz = 0.5 * (centz[k+1][j][i].x - centz[k-1][j][i].x);
                    dydz = 0.5 * (centz[k+1][j][i].y - centz[k-1][j][i].y);
                    dzdz = 0.5 * (centz[k+1][j][i].z - centz[k-1][j][i].z);
                }

                kcsi[k][j][i].x = dyde * dzdz - dzde * dydz;
                kcsi[k][j][i].y =-dxde * dzdz + dzde * dxdz;
                kcsi[k][j][i].z = dxde * dydz - dyde * dxdz;

                keta[k][j][i].x = dydz * dzdc - dzdz * dydc;
                keta[k][j][i].y =-dxdz * dzdc + dzdz * dxdc;
                keta[k][j][i].z = dxdz * dydc - dydz * dxdc;

                kzet[k][j][i].x = dydc * dzde - dzdc * dyde;
                kzet[k][j][i].y =-dxdc * dzde + dzdc * dxde;
                kzet[k][j][i].z = dxdc * dyde - dydc * dxde;

                kaj[k][j][i]
                =
                dxdc * (dyde * dzdz - dzde * dydz) -
                dydc * (dxde * dzdz - dzde * dxdz) +
                dzdc * (dxde * dydz - dyde * dxdz);

                kaj[k][j][i] = 1./kaj[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(fda, KCsi, &kcsi);
    DMDAVecRestoreArray(fda, KEta, &keta);
    DMDAVecRestoreArray(fda, KZet, &kzet);
    DMDAVecRestoreArray(da,  KAj,  &kaj);
    DMDAVecRestoreArray(fda, Centz, &centz);

    DMDAVecRestoreArray(fda, lCoor, &coor);

    VecAssemblyBegin(mesh->Cent); VecAssemblyEnd(mesh->Cent);

    VecAssemblyBegin(Csi);  VecAssemblyEnd(Csi);
    VecAssemblyBegin(Eta);  VecAssemblyEnd(Eta);
    VecAssemblyBegin(Zet);  VecAssemblyEnd(Zet);
    VecAssemblyBegin(Aj);   VecAssemblyEnd(Aj);
    VecAssemblyBegin(ICsi); VecAssemblyEnd(ICsi);
    VecAssemblyBegin(IEta); VecAssemblyEnd(IEta);
    VecAssemblyBegin(IZet); VecAssemblyEnd(IZet);
    VecAssemblyBegin(IAj);  VecAssemblyEnd(IAj);
    VecAssemblyBegin(JCsi); VecAssemblyEnd(JCsi);
    VecAssemblyBegin(JEta); VecAssemblyEnd(JEta);
    VecAssemblyBegin(JZet); VecAssemblyEnd(JZet);
    VecAssemblyBegin(JAj);  VecAssemblyEnd(JAj);
    VecAssemblyBegin(KCsi); VecAssemblyEnd(KCsi);
    VecAssemblyBegin(KEta); VecAssemblyEnd(KEta);
    VecAssemblyBegin(KZet); VecAssemblyEnd(KZet);
    VecAssemblyBegin(KAj);  VecAssemblyEnd(KAj);

    // scatter to local
    DMGlobalToLocalBegin(fda, Csi, INSERT_VALUES, mesh->lCsi);
    DMGlobalToLocalEnd  (fda, Csi, INSERT_VALUES, mesh->lCsi);

    DMGlobalToLocalBegin(fda, Eta, INSERT_VALUES, mesh->lEta);
    DMGlobalToLocalEnd  (fda, Eta, INSERT_VALUES, mesh->lEta);

    DMGlobalToLocalBegin(fda, Zet, INSERT_VALUES, mesh->lZet);
    DMGlobalToLocalEnd  (fda, Zet, INSERT_VALUES, mesh->lZet);

    DMGlobalToLocalBegin(fda, ICsi, INSERT_VALUES, mesh->lICsi);
    DMGlobalToLocalEnd  (fda, ICsi, INSERT_VALUES, mesh->lICsi);

    DMGlobalToLocalBegin(fda, IEta, INSERT_VALUES, mesh->lIEta);
    DMGlobalToLocalEnd  (fda, IEta, INSERT_VALUES, mesh->lIEta);

    DMGlobalToLocalBegin(fda, IZet, INSERT_VALUES, mesh->lIZet);
    DMGlobalToLocalEnd  (fda, IZet, INSERT_VALUES, mesh->lIZet);

    DMGlobalToLocalBegin(fda, JCsi, INSERT_VALUES, mesh->lJCsi);
    DMGlobalToLocalEnd  (fda, JCsi, INSERT_VALUES, mesh->lJCsi);

    DMGlobalToLocalBegin(fda, JEta, INSERT_VALUES, mesh->lJEta);
    DMGlobalToLocalEnd  (fda, JEta, INSERT_VALUES, mesh->lJEta);

    DMGlobalToLocalBegin(fda, JZet, INSERT_VALUES, mesh->lJZet);
    DMGlobalToLocalEnd  (fda, JZet, INSERT_VALUES, mesh->lJZet);

    DMGlobalToLocalBegin(fda, KCsi, INSERT_VALUES, mesh->lKCsi);
    DMGlobalToLocalEnd  (fda, KCsi, INSERT_VALUES, mesh->lKCsi);

    DMGlobalToLocalBegin(fda, KEta, INSERT_VALUES, mesh->lKEta);
    DMGlobalToLocalEnd  (fda, KEta, INSERT_VALUES, mesh->lKEta);

    DMGlobalToLocalBegin(fda, KZet, INSERT_VALUES, mesh->lKZet);
    DMGlobalToLocalEnd  (fda, KZet, INSERT_VALUES, mesh->lKZet);

    DMGlobalToLocalBegin(da, Aj, INSERT_VALUES, mesh->lAj);
    DMGlobalToLocalEnd  (da, Aj, INSERT_VALUES, mesh->lAj);

    DMGlobalToLocalBegin(da, IAj, INSERT_VALUES, mesh->lIAj);
    DMGlobalToLocalEnd  (da, IAj, INSERT_VALUES, mesh->lIAj);

    DMGlobalToLocalBegin(da, JAj, INSERT_VALUES, mesh->lJAj);
    DMGlobalToLocalEnd  (da, JAj, INSERT_VALUES, mesh->lJAj);

    DMGlobalToLocalBegin(da, KAj, INSERT_VALUES, mesh->lKAj);
    DMGlobalToLocalEnd  (da, KAj, INSERT_VALUES, mesh->lKAj);

    DMGlobalToLocalBegin(fda, mesh->Cent, INSERT_VALUES, mesh->lCent);
    DMGlobalToLocalEnd  (fda, mesh->Cent, INSERT_VALUES, mesh->lCent);

    // destroy global vectors
    VecDestroy(&Csi);
    VecDestroy(&Eta);
    VecDestroy(&Zet);

    VecDestroy(&ICsi);
    VecDestroy(&IEta);
    VecDestroy(&IZet);

    VecDestroy(&JCsi);
    VecDestroy(&JEta);
    VecDestroy(&JZet);

    VecDestroy(&KCsi);
    VecDestroy(&KEta);
    VecDestroy(&KZet);

    VecDestroy(&Aj);
    VecDestroy(&IAj);
    VecDestroy(&JAj);
    VecDestroy(&KAj);

    VecDestroy(&Centx);
    VecDestroy(&Centy);
    VecDestroy(&Centz);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetBoundingBox(mesh_ *mesh)
{
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs  = info.xs, xe = info.xs + info.xm;
    PetscInt      ys  = info.ys, ye = info.ys + info.ym;
    PetscInt      zs  = info.zs, ze = info.zs + info.zm;
    PetscInt      mx  = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    Vec           Coor;
    PetscReal     lxmin, lymin, lzmin,       // local min coordinates
                  lxmax, lymax, lzmax,       // local max coordinates
                  xmin, ymin, zmin,          // global min coordinates
                  xmax, ymax, zmax,          // global max coordinates
                  Lx, Ly, Lz;                // domain extension

    Cmpnts        ***coor;                   // point coordinates

    MPI_Comm      MESH_COMM = mesh->MESH_COMM;// this mesh communicator

    PetscMPIInt rank;     MPI_Comm_rank(MESH_COMM, &rank);

    DMGetCoordinatesLocal(da, &Coor);
    DMDAVecGetArray(fda, Coor, &coor);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // loop over cells and get domain extension, min coords and max coords
    lxmax = -1E30; lxmin = 1E30;
    lymax = -1E30; lymin = 1E30;
    lzmax = -1E30; lzmin = 1E30;

    for(k=zs; k<lze; k++)
    {
        for(j=ys; j<lye; j++)
        {
            for(i=xs; i<lxe; i++)
            {
                lxmin = PetscMin(lxmin, coor[k][j][i].x);
                lxmax = PetscMax(lxmax, coor[k][j][i].x);
                lymin = PetscMin(lymin, coor[k][j][i].y);
                lymax = PetscMax(lymax, coor[k][j][i].y);
                lzmin = PetscMin(lzmin, coor[k][j][i].z);
                lzmax = PetscMax(lzmax, coor[k][j][i].z);
            }
        }
    }

    // reduce the values
    MPI_Allreduce(&lxmin, &xmin, 1, MPIU_REAL, MPIU_MIN, MESH_COMM);
    MPI_Allreduce(&lymin, &ymin, 1, MPIU_REAL, MPIU_MIN, MESH_COMM);
    MPI_Allreduce(&lzmin, &zmin, 1, MPIU_REAL, MPIU_MIN, MESH_COMM);
    MPI_Allreduce(&lxmax, &xmax, 1, MPIU_REAL, MPIU_MAX, MESH_COMM);
    MPI_Allreduce(&lymax, &ymax, 1, MPIU_REAL, MPIU_MAX, MESH_COMM);
    MPI_Allreduce(&lzmax, &zmax, 1, MPIU_REAL, MPIU_MAX, MESH_COMM);

    Lx = xmax - xmin;
    Ly = ymax - ymin;
    Lz = zmax - zmin;

    // set domain bunding box information
    mesh->bounds.Lx   = Lx;
    mesh->bounds.Ly   = Ly;
    mesh->bounds.Lz   = Lz;

    mesh->bounds.xmin = xmin;
    mesh->bounds.xmax = xmax;
    mesh->bounds.ymin = ymin;
    mesh->bounds.ymax = ymax;
    mesh->bounds.zmin = zmin;
    mesh->bounds.zmax = zmax;

    DMDAVecRestoreArray(fda, Coor, &coor);

    PetscPrintf(mesh->MESH_COMM, "Domain bounding box : xmin = %lf, xmax = %lf\n", xmin, xmax);
    PetscPrintf(mesh->MESH_COMM, "                      ymin = %lf, ymax = %lf\n", ymin, ymax);
    PetscPrintf(mesh->MESH_COMM, "                      zmin = %lf, zmax = %lf\n", zmin, zmax);
    PetscPrintf(mesh->MESH_COMM, "Domain size         : Lx = %lf, Ly = %lf, Lz = %lf\n\n", Lx, Ly, Lz);

    //set ground level 
    mesh->grndLevel = mesh->bounds.zmin;
    
    return(0);
}

PetscErrorCode ghostnodesCellcenter(mesh_ *mesh){
    //cell center for the ghost nodes - as ghost node coordinate info is not available
    // assumed to be one cell length from the first internal cell, depending on the direction

    DM               da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;
    PetscReal        ***aj;
    Cmpnts           ***cent, ***lcent, ***csi, ***eta, ***zet;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->Cent, &cent);
    DMDAVecGetArray(fda, mesh->lCent, &lcent);
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    PetscReal areaI, areaJ, areaK;
    PetscReal lenI, lenJ, lenK;

    // currently only for cartesian mesh
    // can be generalized: find the normal and distance from cell center to face center of the first internal node
    // ghost node coordinate = 2*distance along the normal

    if (xs==0)
    {
        i = xs;
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                areaI = sqrt
                        (
                                csi[k][j][i+1].x*csi[k][j][i+1].x +
                                csi[k][j][i+1].y*csi[k][j][i+1].y +
                                csi[k][j][i+1].z*csi[k][j][i+1].z
                        );

                lenI = 1.0/aj[k][j][i+1]/areaI;

                cent[k][j][i].x = lcent[k][j][i+1].x;
                cent[k][j][i].y = lcent[k][j][i+1].y - lenI;
                cent[k][j][i].z = lcent[k][j][i+1].z;
            }
        }
    }

    if (xe==mx)
    {
        i = xe-1;
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                areaI = sqrt
                        (
                                csi[k][j][i-1].x*csi[k][j][i-1].x +
                                csi[k][j][i-1].y*csi[k][j][i-1].y +
                                csi[k][j][i-1].z*csi[k][j][i-1].z
                        );

                lenI = 1.0/aj[k][j][i-1]/areaI;

                cent[k][j][i].x = lcent[k][j][i-1].x;
                cent[k][j][i].y = lcent[k][j][i-1].y + lenI;
                cent[k][j][i].z = lcent[k][j][i-1].z;
            }
        }
    }

    if (ys==0)
    {
        j = ys;
        for (k=lzs; k<lze; k++)
        {
            for (i=lxs; i<lxe; i++)
            {

                areaJ = sqrt
                        (
                            eta[k][j+1][i].x*eta[k][j+1][i].x +
                            eta[k][j+1][i].y*eta[k][j+1][i].y +
                            eta[k][j+1][i].z*eta[k][j+1][i].z
                        );

                lenJ = 1.0/aj[k][j+1][i]/areaJ;

                cent[k][j][i].x = lcent[k][j+1][i].x;
                cent[k][j][i].y = lcent[k][j+1][i].y;
                cent[k][j][i].z = lcent[k][j+1][i].z  - lenJ;

            }
        }
    }

    if (ye==my)
    {
        j = ye-1;
        for (k=lzs; k<lze; k++)
        {
            for (i=lxs; i<lxe; i++)
            {

                areaJ = sqrt
                        (
                            eta[k][j-1][i].x*eta[k][j-1][i].x +
                            eta[k][j-1][i].y*eta[k][j-1][i].y +
                            eta[k][j-1][i].z*eta[k][j-1][i].z
                        );

                lenJ = 1.0/aj[k][j-1][i]/areaJ;

                cent[k][j][i].x = lcent[k][j-1][i].x;
                cent[k][j][i].y = lcent[k][j-1][i].y;
                cent[k][j][i].z = lcent[k][j-1][i].z  + lenJ;
            }
        }
    }

    if (zs==0)
    {
        k = zs;
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                areaK = sqrt
                        (
                            zet[k+1][j][i].x*zet[k+1][j][i].x +
                            zet[k+1][j][i].y*zet[k+1][j][i].y +
                            zet[k+1][j][i].z*zet[k+1][j][i].z
                        );

                lenK = 1.0/aj[k+1][j][i]/areaK;

                cent[k][j][i].x = lcent[k+1][j][i].x - lenK;
                cent[k][j][i].y = lcent[k+1][j][i].y;
                cent[k][j][i].z = lcent[k+1][j][i].z;

            }
        }
    }

    if (ze==mz)
    {
        k = ze-1;
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                areaK = sqrt
                        (
                            zet[k-1][j][i].x*zet[k-1][j][i].x +
                            zet[k-1][j][i].y*zet[k-1][j][i].y +
                            zet[k-1][j][i].z*zet[k-1][j][i].z
                        );

                lenK = 1.0/aj[k-1][j][i]/areaK;

                cent[k][j][i].x = lcent[k-1][j][i].x + lenK;
                cent[k][j][i].y = lcent[k-1][j][i].y;
                cent[k][j][i].z = lcent[k-1][j][i].z;
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->Cent, &cent);
    DMDAVecRestoreArray(fda, mesh->lCent, &lcent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    DMGlobalToLocalBegin(fda, mesh->Cent, INSERT_VALUES, mesh->lCent);
    DMGlobalToLocalEnd(fda, mesh->Cent, INSERT_VALUES, mesh->lCent);

    return 0;
}
