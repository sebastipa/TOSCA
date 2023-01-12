#ifndef CatalystAdaptor_h
#define CatalystAdaptor_h

#include <catalyst.h>
#include <stdio.h>
#include <string.h>

// #include "include/base.h"
// #include "include/domain.h"
// #include "include/io.h"
// #include "include/inline.h"

//***************************************************************************************************************//

void catalystInitialize(int argc, char* argv[])
{
    conduit_node* catalyst_init_params = conduit_node_create();
    for (int cc = 1; cc < argc; ++cc)
    {
        if (strcmp(argv[cc], "--output") == 0 && (cc + 1) < argc)
        {
            std::string fileName(argv[cc + 1]);
            std::string path2save = "fields/" + fileName;
            conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/type", "io");
            conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/filename", path2save.c_str());
            conduit_node_set_path_char8_str(catalyst_init_params, "catalyst/pipelines/0/channel", "grid");
            ++cc;
        }
        else
        {
            // pass scripts on command line
            char buf[256];
            snprintf(buf, 256, "catalyst/scripts/script%d", (cc - 1));
            conduit_node_set_path_char8_str(catalyst_init_params, buf, argv[cc]);
        }
    }

    conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/implementation", "paraview");
    conduit_node_set_path_char8_str(catalyst_init_params, "catalyst_load/search_paths/paraview", PARAVIEW_IMPL_DIR);
    enum catalyst_status err = catalyst_initialize(catalyst_init_params);
    conduit_node_destroy(catalyst_init_params);

    if (err != catalyst_status_ok)
    {
        printf("Failed to initialize Catalyst: %d\n", err);
    }
}

//***************************************************************************************************************//

void do_catalyst_execute(int cycle, double time, /*mesh_ *domainMesh, */domain_ *domain) {
  // add time/cycle information
  conduit_node* catalyst_exec_params = conduit_node_create();
  conduit_node_set_path_int64(catalyst_exec_params, "catalyst/state/timestep", cycle);
  conduit_node_set_path_float64(catalyst_exec_params, "catalyst/state/time", time);

  //create the mesh in "grid" channel, using Conduit Mesh Blueprint
  conduit_node_set_path_char8_str(catalyst_exec_params, "catalyst/channels/grid/type", "mesh");
  conduit_node* mesh = conduit_node_create();

  mesh_ *domainMesh  = domain->mesh;
  int mpiRank = 0;
  MPI_Comm_rank(domainMesh->MESH_COMM, &mpiRank);
  DM fda = domainMesh->fda, da = domainMesh->da;
  Vec Coor;
  Cmpnts ***coor;
  DMGetCoordinatesLocal(da, &Coor);  // get the local vector Coor linked to the coordinates
  DMDAVecGetArray(fda, Coor, &coor); // access this vector using the 3D array coor
  DMDALocalInfo info = domainMesh->info;
  PetscInt i, j, k;
  PetscInt xs = info.xs, xe = info.xs + info.xm;
  PetscInt ys = info.ys, ye = info.ys + info.ym;
  PetscInt zs = info.zs, ze = info.zs + info.zm;
  PetscInt mx = info.mx, my = info.my, mz = info.mz;
  PetscInt           lxs, lxe, lys, lye, lzs, lze;
  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;
  int numberOfCells = (lxe-lxs) * (lye-lys) * (lze-lzs);

  // placeholder for point-based variables
  // int numberOfPoints = (lxe-lxs+1) * (lye-lys+1) * (lze-lzs+1);
  // double *pointVar; pointVar = new double[numberOfPoints];
  // for (k=lzs; k<lze; k++)
  //   for (j=lys; j<lye; j++)
  //     for (i=lxs; i<lxe; i++) {
  // 	//pointVar[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].x;
  // 	printf("%d %f %f %f\n", mpiRank, coor[k][j][i].x, coor[k][j][i].y, coor[k][j][i].z);
  //     }

  PetscReal ***pressure, ***dummy, ***qcriterion;
  DMDAVecGetArray(da, domain->peqn->P, &pressure);
  DMDAVecGetArray(da, domain->peqn->P, &dummy);
  DMDAVecGetArray(da, domain->acquisition->fields->Q, &qcriterion);
  double *cellVar; cellVar = new double[numberOfCells];

  int cartesian = 1, rectilinear = 2, unstructured = 3;
  int gridType = rectilinear;
  if (gridType == cartesian) {
    // add coordsets
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "uniform");

    conduit_node_set_path_int64(mesh, "coordsets/coords/dims/i", lxe-lxs+1); // number of points in each dimension
    conduit_node_set_path_int64(mesh, "coordsets/coords/dims/j", lye-lys+1);
    conduit_node_set_path_int64(mesh, "coordsets/coords/dims/k", lze-lzs+1);

    conduit_node_set_path_float64(mesh, "coordsets/coords/origin/x", lxs);
    conduit_node_set_path_float64(mesh, "coordsets/coords/origin/y", lys);
    conduit_node_set_path_float64(mesh, "coordsets/coords/origin/z", lzs);

    conduit_node_set_path_float64(mesh, "coordsets/coords/spacing/x", 0.1);
    conduit_node_set_path_float64(mesh, "coordsets/coords/spacing/y", 0.1);
    conduit_node_set_path_float64(mesh, "coordsets/coords/spacing/z", 0.1);

    // add topology
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "uniform");
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords");

    // add cell-based pressure
    for (k=lzs; k<lze; k++)
      for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = pressure[k-lzs+zs][j-lys+ys][i-lxs+xs];
	}
    conduit_node_set_path_char8_str(mesh, "fields/pressure/association", "element");
    conduit_node_set_path_char8_str(mesh, "fields/pressure/topology", "mesh");
    conduit_node_set_path_char8_str(mesh, "fields/pressure/volume_dependent", "false");
    conduit_node_set_path_external_float64_ptr(mesh, "fields/pressure/values", cellVar, numberOfCells);

    // add the mesh info (conduit mesh) to catalyst_exec_params
    conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh);

    enum catalyst_status err = catalyst_execute(catalyst_exec_params);
    if (err != catalyst_status_ok)
      printf("Failed to execute Catalyst: %d\n", err);
  }
  else if (gridType == rectilinear) {
    // add coordsets
    conduit_node_set_path_char8_str(mesh, "coordsets/coords/type", "rectilinear");

    // add topology
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/type", "rectilinear");
    conduit_node_set_path_char8_str(mesh, "topologies/mesh/coordset", "coords");

    double *xval; xval = new double[lxe-lxs+1];
    double *yval; yval = new double[lye-lys+1];
    double *zval; zval = new double[lze-lzs+1];
    for (i=lxs; i<lxe; i++) {
      xval[i-lxs] = coor[lzs][lys][i].y;
      // printf("%d > %f\n", mpiRank, xval[i-lxs]);
    }
    for (j=lys; j<lye; j++)
      yval[j-lys] = coor[lzs][j][lxs].z;
    for (k=lzs; k<lze; k++)
      zval[k-lzs] = coor[k][lys][lxs].x;
    // fill in the last point in each dimension
    int lastIndex = lxe-lxs; xval[lastIndex] = 2*xval[lastIndex-1] - xval[lastIndex-2];
    lastIndex = lye-lys; yval[lastIndex] = 2*yval[lastIndex-1] - yval[lastIndex-2];
    lastIndex = lze-lzs; zval[lastIndex] = 2*zval[lastIndex-1] - zval[lastIndex-2];
    //printf("%d > %d %d %d %d %d %d\n", mpiRank, lxs, lxe, lys, lye, lzs, lze);

    conduit_node_set_path_float64_ptr(mesh, "coordsets/coords/values/x", xval, lxe-lxs+1);
    conduit_node_set_path_float64_ptr(mesh, "coordsets/coords/values/y", yval, lye-lys+1);
    conduit_node_set_path_float64_ptr(mesh, "coordsets/coords/values/z", zval, lze-lzs+1);

    // add cell-based pressure
    for (k=lzs; k<lze; k++)
      for (j=lys; j<lye; j++)
	for (i=lxs; i<lxe; i++) {
	  cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = pressure[k-lzs+zs][j-lys+ys][i-lxs+xs];
	}
    conduit_node_set_path_char8_str(mesh, "fields/pressure/association", "element");
    conduit_node_set_path_char8_str(mesh, "fields/pressure/topology", "mesh");
    conduit_node_set_path_char8_str(mesh, "fields/pressure/volume_dependent", "false");
    conduit_node_set_path_float64_ptr(mesh, "fields/pressure/values", cellVar, numberOfCells);

    // add cell-based dummy
    for (k=lzs; k<lze; k++)
      for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++) {
          cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = 5. + abs(dummy[k-lzs+zs][j-lys+ys][i-lxs+xs])/1000.;
        }
    conduit_node_set_path_char8_str(mesh, "fields/dummy/association", "element");
    conduit_node_set_path_char8_str(mesh, "fields/dummy/topology", "mesh");
    conduit_node_set_path_char8_str(mesh, "fields/dummy/volume_dependent", "false");
    conduit_node_set_path_float64_ptr(mesh, "fields/dummy/values", cellVar, numberOfCells);

    // add cell-based qcriterion
    for (k=lzs; k<lze; k++)
      for (j=lys; j<lye; j++)
    	for (i=lxs; i<lxe; i++) {
    	  cellVar[(k-lzs)*(lxe-lxs)*(lye-lys) + (j-lys)*(lxe-lxs) + (i-lxs)] = qcriterion[k-lzs+zs][j-lys+ys][i-lxs+xs];
    	}
    conduit_node_set_path_char8_str(mesh, "fields/qcriterion/association", "element");
    conduit_node_set_path_char8_str(mesh, "fields/qcriterion/topology", "mesh");
    conduit_node_set_path_char8_str(mesh, "fields/qcriterion/volume_dependent", "false");
    conduit_node_set_path_float64_ptr(mesh, "fields/qcriterion/values", cellVar, numberOfCells);

    // conduit_node_print(mesh); // for debugging

    // add the mesh info (conduit mesh) to catalyst_exec_params
    conduit_node_set_path_external_node(catalyst_exec_params, "catalyst/channels/grid/data", mesh);

    enum catalyst_status err = catalyst_execute(catalyst_exec_params);
    if (err != catalyst_status_ok)
      printf("Failed to execute Catalyst: %d\n", err);

    delete [] xval, yval, zval;
  }
  else if (gridType == unstructured) {
    // placeholder for writing coordinates to an unstructured grid
    // conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/x",
    //   /*data=*/pointVar, /*num_elements=*/numberOfPoints, /*offset=*/0,
    //   /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    //   /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
    // for (k=zs; k<ze-1; k++)
    //   for (j=ys; j<ye-1; j++)
    //     for (i=xs; i<xe-1; i++) {
    //       pointVar[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].y;
    //     }
    // conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/y",
    //   /*data=*/pointVar, /*num_elements=*/numberOfPoints, /*offset=*/0,
    //   /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    //   /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
    // for (k=zs; k<ze-1; k++)
    //   for (j=ys; j<ye-1; j++)
    //     for (i=xs; i<xe-1; i++) {
    //       pointVar[(k-zs) * (mx-1)*(my-1) + (j-ys)*(mx-1) + (i-xs)] = coor[k][j][i].z;
    //     }
    // conduit_node_set_path_external_float64_ptr_detailed(mesh, "coordsets/coords/values/z",
    //   /*data=*/pointVar, /*num_elements=*/numberOfPoints, /*offset=*/0,
    //   /*stride=*/sizeof(double), /*element_bytes=*/sizeof(double),
    //   /*endianness=*/CONDUIT_ENDIANNESS_DEFAULT_ID);
  }

  DMDAVecRestoreArray(fda, Coor, &coor);
  DMDAVecRestoreArray(da, domain->peqn->P, &pressure);
  DMDAVecRestoreArray(da, domain->peqn->P, &dummy);
  DMDAVecRestoreArray(da, domain->acquisition->fields->Q, &qcriterion);

  // delete [] pointVar;
  delete [] cellVar;

  conduit_node_destroy(catalyst_exec_params);
  conduit_node_destroy(mesh);
}

//***************************************************************************************************************//

void catalystFinalize()
{
    conduit_node* catalyst_fini_params = conduit_node_create();
    enum catalyst_status err = catalyst_finalize(catalyst_fini_params);
    if (err != catalyst_status_ok)
    printf("Failed to execute Catalyst: %d\n", err);
    conduit_node_destroy(catalyst_fini_params);
}

#endif
