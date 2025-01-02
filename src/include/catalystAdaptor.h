#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.h>
#include <stdio.h>
#include <string.h>

//***************************************************************************************************************//

void catalystInitialize(domain_ *domain)
{
    // read input parameter
    io_   *io   = domain->io;
    mesh_ *mesh = domain->mesh;

    PetscPrintf(mesh->MESH_COMM, "Initializing TOSCA's catalyst-paraview interface...");

    readDictWord  ("sampling/catalystProperties", "ioType",      &(io->ioTypeCatalyst));
    readDictWord  ("sampling/catalystProperties", "outputType",  &(io->outputTypeCatalyst));
    readDictDouble("sampling/catalystProperties", "startTime",   &(io->startTimeCatalyst));
    readDictDouble("sampling/catalystProperties", "timeInterval",&(io->timeIntervalCatalyst));

    if(io->outputTypeCatalyst!="adjustableTime" && io->outputTypeCatalyst!="timeStep")
    {
        char error[512];
        sprintf(error, "unknown catalyst outputType. Available types are\n - adjustableTime\n - timeStep\n");
        fatalErrorInFunction("catalystInitialize",  error);
    }

    // create catalyst folder
    PetscMPIInt rank;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if(!rank)
    {
        errno = 0;
        PetscInt dirRes = mkdir("./catalyst", 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
            char error[512];
            sprintf(error, "could not create catalyst directory\n");
            fatalErrorInFunction("catalystInitialize",  error);
        }
    }

    if(io->ioTypeCatalyst == "general") {}
    else if(io->ioTypeCatalyst == "script")
    {
        readDictWord("sampling/catalystProperties", "scriptName", &(io->scriptNameCatalyst));
    }
    else
    {
        char error[512];
        sprintf(error, "unknown catalyst ioType. Available types are\n - general\n - script\n");
        fatalErrorInFunction("catalystInitialize",  error);
    }

    // initialize catalyst
    conduit_node* catalyst_init_params = conduit_node_create();
    if(io->ioTypeCatalyst == "general")
    {
        std::string path2save = "catalyst/dataset-%04t.vtpd";
        conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/type", "io");
        conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/filename", path2save.c_str());
        conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/channel", "grid");
    }
    else if(io->ioTypeCatalyst == "script")
    {
        // pass scripts on command line
        word scriptName = "./sampling/" + io->scriptNameCatalyst;
        conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/scripts/script0", scriptName.c_str());
    }

    conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/implementation", "paraview");
    conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/search_paths/paraview", PARAVIEW_IMPL_DIR);
    enum catalyst_status err = catalyst_initialize(catalyst_init_params);
    conduit_node_destroy(catalyst_init_params);

    if (err != catalyst_status_ok)
    {
        char error[512];
        sprintf(error, "\nTOSCA's catalyst-paraview interface initialization failed\n");
        fatalErrorInFunction("catalystInitialize",  error);
    }

    PetscPrintf(mesh->MESH_COMM, "done\n\n");
}

//***************************************************************************************************************//

void catalystExecute(domain_ *domain)
{
    PetscInt  executeCatalystOperations = 0;

    word      intervalType = domain->io->outputTypeCatalyst;
    PetscReal timeInterval = domain->io->timeIntervalCatalyst;

    // write every "timeInterval" seconds
    if
    (
        intervalType == "adjustableTime" &&
        mustWrite(domain->clock->time, domain->io->startTimeCatalyst, timeInterval)
    )
    {
        executeCatalystOperations = 1;
    }
    // write every "timeInterval" iterations
    else if
    (
        (intervalType == "timeStep") &&
        (
            domain->clock->it / timeInterval -
            std::floor
            (
                domain->clock->it / timeInterval
            ) < 1e-10
        )
    )
    {
        executeCatalystOperations = 1;
    }

    if(executeCatalystOperations)
    {
        // add time/cycle information
        conduit_node* catalyst_exec_params = conduit_node_create();

        conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", domain->clock->it);
        conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", domain->clock->time);

        //create the mesh in "grid" channel, using Conduit Mesh Blueprint
        conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh");
        conduit_node* condMesh = conduit_node_create();

        mesh_         *mesh  = domain->mesh;
        DM            fda = mesh->fda, da = mesh->da;
        DMDALocalInfo info = mesh->info;
        PetscInt      i, j, k;
        PetscInt      xs = info.xs, xe = info.xs + info.xm;
        PetscInt      ys = info.ys, ye = info.ys + info.ym;
        PetscInt      zs = info.zs, ze = info.zs + info.zm;
        PetscInt      mx = info.mx, my = info.my, mz = info.mz;

        PetscInt      lxs, lxe, lys, lye, lzs, lze;

        Vec           Coor;
        Cmpnts        ***coor, ***ucat;
        PetscReal     ***qcrit, ***p;

        lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
        lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
        lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

        PetscInt  numberOfCells = (lxe-lxs) * (lye-lys) * (lze-lzs);
        PetscReal *cellVar; cellVar = new PetscReal[numberOfCells];

        PetscReal ts, te;
        PetscTime(&ts);

        PetscPrintf(mesh->MESH_COMM, "Executed paraview-catalyst actions in ");

        DMGetCoordinatesLocal(da, &Coor);  // get the local vector Coor linked to the coordinates
        DMDAVecGetArray(fda, Coor, &coor); // access this vector using the 3D array coor

        // store Qcrit inside acquisition
        computeQCritIO(domain->acquisition);

        DMDAVecGetArray(fda, domain->ueqn->Ucat, &ucat);
        DMDAVecGetArray(da, domain->peqn->P, &p);
        DMDAVecGetArray(da,  domain->acquisition->fields->Q, &qcrit);

        PetscInt cartesian    = 1,
                 rectilinear  = 2,
                 unstructured = 3;
        PetscInt gridType     = rectilinear;

        if (gridType == cartesian)
        {
            // add coordsets
            conduit_node_set_path_char8_str(condMesh, "coordsets/coords/type", "uniform");

            // add topology
            conduit_node_set_path_char8_str(condMesh, "topologies/mesh/type", "uniform");
            conduit_node_set_path_char8_str(condMesh, "topologies/mesh/coordset", "coords");

            conduit_node_set_path_int64(condMesh, "coordsets/coords/dims/i", lze-lzs+1); // number of points in each dimension
            conduit_node_set_path_int64(condMesh, "coordsets/coords/dims/j", lxe-lxs+1);
            conduit_node_set_path_int64(condMesh, "coordsets/coords/dims/k", lye-lys+1);

            conduit_node_set_path_float64(condMesh, "coordsets/coords/origin/x", lzs);
            conduit_node_set_path_float64(condMesh, "coordsets/coords/origin/y", lxs);
            conduit_node_set_path_float64(condMesh, "coordsets/coords/origin/z", lys);

            conduit_node_set_path_float64(condMesh, "coordsets/coords/spacing/dx", 0.1);
            conduit_node_set_path_float64(condMesh, "coordsets/coords/spacing/dy", 0.1);
            conduit_node_set_path_float64(condMesh, "coordsets/coords/spacing/dz", 0.1);

            // add cell-based q-criterion
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = qcrit[k-lzs+zs][j-lys+ys][i-lxs+xs];
            }

            conduit_node_set_path_char8_str(condMesh, "fields/q/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/q/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/q/volume_dependent", "false");
            conduit_node_set_path_external_float64_ptr(condMesh, "fields/q/values", cellVar, numberOfCells);

            // add cell-based pressure
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = p[k-lzs+zs][j-lys+ys][i-lxs+xs];
            }

            conduit_node_set_path_char8_str(condMesh, "fields/p/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/p/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/p/volume_dependent", "false");
            conduit_node_set_path_external_float64_ptr(condMesh, "fields/p/values", cellVar, numberOfCells);

            // add cell-based velocity magnitude
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = nMag(ucat[k-lzs+zs][j-lys+ys][i-lxs+xs]);
            }

            conduit_node_set_path_char8_str(condMesh, "fields/uMag/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/uMag/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/uMag/volume_dependent", "false");
            conduit_node_set_path_external_float64_ptr(condMesh, "fields/uMag/values", cellVar, numberOfCells);

            // add the mesh info (condMesh) to catalyst_exec_params
            conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", condMesh);

            enum catalyst_status err = catalyst_execute(catalyst_exec_params);

            if (err != catalyst_status_ok)
            {
                char error[512];
                sprintf(error, "failed to execute catalyst-paraview interface\n");
                fatalErrorInFunction("catalystInitialize",  error);
            }
        }
        else if (gridType == rectilinear)
        {
            // add coordsets
            conduit_node_set_path_char8_str(condMesh, "coordsets/coords/type", "rectilinear");

            // add topology
            conduit_node_set_path_char8_str(condMesh, "topologies/mesh/type", "rectilinear");
            conduit_node_set_path_char8_str(condMesh, "topologies/mesh/coordset", "coords");

            PetscReal *xval; xval = new PetscReal[lze-lzs+1];
            PetscReal *yval; yval = new PetscReal[lxe-lxs+1];
            PetscReal *zval; zval = new PetscReal[lye-lys+1];

            for (i=lxs; i<lxe; i++) yval[i-lxs] = coor[lzs][lys][i].y;
            for (j=lys; j<lye; j++) zval[j-lys] = coor[lzs][j][lxs].z;
            for (k=lzs; k<lze; k++) xval[k-lzs] = coor[k][lys][lxs].x;

            // fill in the last point in each dimension
            PetscInt lastIndex;
            lastIndex = lxe-lxs; yval[lastIndex] = 2*yval[lastIndex-1] - yval[lastIndex-2];
            lastIndex = lye-lys; zval[lastIndex] = 2*zval[lastIndex-1] - zval[lastIndex-2];
            lastIndex = lze-lzs; xval[lastIndex] = 2*xval[lastIndex-1] - xval[lastIndex-2];

            conduit_node_set_path_float64_ptr(condMesh, "coordsets/coords/values/x", xval, lze-lzs+1);
            conduit_node_set_path_float64_ptr(condMesh, "coordsets/coords/values/y", yval, lxe-lxs+1);
            conduit_node_set_path_float64_ptr(condMesh, "coordsets/coords/values/z", zval, lye-lys+1);

            // add cell-based velocity mag
            /*
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = nMag(ucat[k-lzs+zs][j-lys+ys][i-lxs+xs]);
            }
            */

            // add cell-based velocity mag
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            for (k=lzs; k<lze; k++)
            {
                cellVar[(j-lys)*(lze-lzs)*(lxe-lxs) + (i-lxs)*(lze-lzs) + (k-lzs)] = nMag(ucat[k-lzs+zs][j-lys+ys][i-lxs+xs]);
            }

            conduit_node_set_path_char8_str(condMesh, "fields/uMag/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/uMag/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/uMag/volume_dependent", "false");
            conduit_node_set_path_float64_ptr(condMesh, "fields/uMag/values", cellVar, numberOfCells);

            // add cell-based pressure
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            for (k=lzs; k<lze; k++)
            {
                cellVar[(j-lys)*(lze-lzs)*(lxe-lxs) + (i-lxs)*(lze-lzs) + (k-lzs)] = p[k-lzs+zs][j-lys+ys][i-lxs+xs];
            }

            conduit_node_set_path_char8_str(condMesh, "fields/p/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/p/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/p/volume_dependent", "false");
            conduit_node_set_path_float64_ptr(condMesh, "fields/p/values", cellVar, numberOfCells);

            // add cell-based qcrit
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            for (k=lzs; k<lze; k++)
            {
                cellVar[(j-lys)*(lze-lzs)*(lxe-lxs) + (i-lxs)*(lze-lzs) + (k-lzs)] = qcrit[k-lzs+zs][j-lys+ys][i-lxs+xs];
            }

            conduit_node_set_path_char8_str(condMesh, "fields/q/association", "element");
            conduit_node_set_path_char8_str(condMesh, "fields/q/topology", "mesh");
            conduit_node_set_path_char8_str(condMesh, "fields/q/volume_dependent", "false");
            conduit_node_set_path_float64_ptr(condMesh, "fields/q/values", cellVar, numberOfCells);

            // for debugging
            // conduit_node_print(condMesh);

            // add the mesh info (condMesh) to catalyst_exec_params
            conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", condMesh);

            enum catalyst_status err = catalyst_execute(catalyst_exec_params);

            if (err != catalyst_status_ok)
            {
                char error[512];
                sprintf(error, "failed to execute catalyst-paraview interface\n");
                fatalErrorInFunction("catalystInitialize",  error);
            }

            delete [] xval, yval, zval;
        }
        else if (gridType == unstructured)
        {
          // not implemented
        }

        DMDAVecRestoreArray(fda, Coor, &coor);
        DMDAVecRestoreArray(fda, domain->ueqn->Ucat, &ucat);
        DMDAVecRestoreArray(da,  domain->peqn->P, &p);
        DMDAVecRestoreArray(da,  domain->acquisition->fields->Q, &qcrit);

        delete [] cellVar;

        conduit_node_destroy(catalyst_exec_params);
        conduit_node_destroy(condMesh);

        PetscTime(&te);
        PetscPrintf(mesh->MESH_COMM, "%lf s\n",te-ts);
    }
}

//***************************************************************************************************************//

void catalystFinalize()
{
    conduit_node* catalyst_fini_params = conduit_node_create();
    enum catalyst_status err = catalyst_finalize(catalyst_fini_params);

    if (err != catalyst_status_ok)
    {
        char error[512];
        sprintf(error, "failed to finalize catalyst-paraview interface\n");
        fatalErrorInFunction("catalystFinalize",  error);
    }

    conduit_node_destroy(catalyst_fini_params);
}

#endif
