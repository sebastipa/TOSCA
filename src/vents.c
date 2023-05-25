//! \file  ibm.c
//! \brief Contains Immersed boundary method function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
//#include "include/inflow.h"
//#include "include/wallfunctions.h"


//***************************************************************************************************************//

PetscErrorCode initializeVents(domain_ *domain)
{
  PetscInt nDomains = domain[0].info.nDomains;
  flags_ flags = domain[0].flags;

  if(flags.isVentsActive)
  {

    // loop through the domains to set the domain vents pointer
    for(PetscInt d = 0; d < nDomains; d++)
    {
      readVentsProperties(domain[d].vents);

      MPI_Barrier(domain[d].mesh->MESH_COMM);

      ventSetAndPrint(domain[d].vents);

      MPI_Barrier(domain[d].mesh->MESH_COMM);


    }

  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode readVentsProperties(vents_ *vents)
{
  PetscMPIInt rank;

  mesh_ *mesh = vents->access->mesh;
  flags_ *flags = vents->access->flags;

  MPI_Comm_rank(mesh->MESH_COMM, &rank);

  PetscPrintf(PETSC_COMM_WORLD, "\nVents initialization...\n");

  // read number of vents total and room type. This info flags other functions later.
  readDictInt("./vents/ventsProperties.dat", "numberOfVents", &(vents->numberOfVents));
  //readDictWord("./vents/ventsProperties.dat", "roomPressure", &(vents->roomPressure)); addd IORatio later


  // allocate memory for each vent object
  PetscMalloc(vents->numberOfVents * sizeof(ventObject*), &(vents->vent));

  char ventName[256];

  for  (PetscInt  i=0; i < vents->numberOfVents; i++)
  {
      sprintf(ventName, "vent%ld", i);

    PetscMalloc(sizeof(ventObject), &(vents->vent[i]));

    //read vent info, start with vent location and direction of flow info
    readSubDictWord("./vents/ventsProperties.dat", ventName, "face", &(vents->vent[i]->face));
    readSubDictWord("./vents/ventsProperties.dat", ventName, "dir", &(vents->vent[i]->dir));

    if (vents->vent[i]->face  == "iLeft" || vents->vent[i]->face == "iRight")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "xBound1", &(vents->vent[i]->xBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "xBound2", &(vents->vent[i]->xBound2));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "zBound1", &(vents->vent[i]->zBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "zBound2", &(vents->vent[i]->zBound2));

    }

    if (vents->vent[i]->face  == "jLeft" || vents->vent[i]->face == "jRight")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "xBound1", &(vents->vent[i]->xBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "xBound2", &(vents->vent[i]->xBound2));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "yBound1", &(vents->vent[i]->yBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "yBound2", &(vents->vent[i]->yBound2));

    }

    if (vents->vent[i]->face  == "kLeft" || vents->vent[i]->face == "kRight")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "yBound1", &(vents->vent[i]->yBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "yBound2", &(vents->vent[i]->yBound2));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "zBound1", &(vents->vent[i]->zBound1));

        readSubDictDouble("./vents/ventsProperties.dat", ventName, "zBound2", &(vents->vent[i]->zBound2));

    }

    // read U and Nut BCs
    readSubDictWord("./vents/ventsProperties.dat", ventName, "ventBC", &(vents->vent[i]->ventBC));

    readSubDictWord("./vents/ventsProperties.dat", ventName, "ventBCNut", &(vents->vent[i]->ventBCNut));

    //read fixed value flux value if needed
    if (vents->vent[i]->ventBC == "fixedValue")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventDesiredFlux", &(vents->vent[i]->ventDesiredFlux));
    }

    //nut BC is usually ZG for vents. May change later.
    if (vents->vent[i]->ventBCNut == "fixedValue")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventBCNutReal", &(vents->vent[i]->ventBCNutReal));
    }

    // vent velocity, area, # of cells, all calculated in SetVentsAndPrint fnc using the info read above.
    vents->vent[i]->nCellsVent = 0;
    vents->vent[i]->lFluxVent = 0; //lFluxVent should converge to ventDesiredFlux.
    vents->vent[i]->ventArea = 0;
    vents->vent[i]->ventBCVec.x = 0;
    vents->vent[i]->ventBCVec.y = 0;
    vents->vent[i]->ventBCVec.z = 0;
    //vents->desiredLeakFlux = 0;

    /*if (vents->vent[i]->dir == "leak")
    {
        //printf("here...................");
        vents->desiredLeakFlux += vents->vent[i]->ventDesiredFlux;
    }*/

    // read temp BC info
    /*if(flags->isTeqnActive)
    {
        readSubDictWord("./vents/ventsProperties.dat", ventName, "ventTBC", &(vents->vent[i]->ventTBC));

        if (vents->vent[i]->ventTBC == "fixedValue" || vents->vent[i]->ventTBC == "fixedGradient")
        {
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventTBCVal", &(vents->vent[i]->ventTBCVal));
        }
    }*/

    // read temp BC info
    /*if(flags->isCeqnActive)
    {
        readSubDictWord("./vents/ventsProperties.dat", ventName, "ventCBC", &(vents->vent[i]->ventCBC));

        if (vents->vent[i]->ventCBC == "fixedValue" || vents->vent[i]->ventCBC == "fixedGradient")
        {
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventCBCVal", &(vents->vent[i]->ventCBCVal));
        }

        if (vents->vent[i]->ventCBC == "intermittent")
        {
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventCBCVal", &(vents->vent[i]->ventCBCVal));
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventCBCStart", &(vents->vent[i]->ventCBCStart));
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventCBCEnd", &(vents->vent[i]->ventCBCEnd));
        }
    }*/

  }

  vents->nCellsVentsOut = 0;
  vents->nCellsVentsIn = 0;
  //vents->nCellsLeak = 0;

  return 0;
}

 //***************************************************************************************************************//

 PetscErrorCode ventSetAndPrint(vents_ *vents)
 {
     mesh_        *mesh = vents->access->mesh;
     DM            da   = mesh->da, fda = mesh->fda;
     DMDALocalInfo info = mesh->info;
     PetscInt      xs   = info.xs, xe = info.xs + info.xm;
     PetscInt      ys   = info.ys, ye = info.ys + info.ym;
     PetscInt      zs   = info.zs, ze = info.zs + info.zm;
     PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

     Cmpnts        ***cent;

     PetscMPIInt   rank;

     PetscInt      lxs, lxe, lys, lye, lzs, lze;
     PetscInt      i, j, k;

     PetscInt      q;

     PetscInt     ***markVent;

     PetscInt     lNumCellsVent[vents->numberOfVents] = { };
     PetscInt     numCellsVent[vents->numberOfVents] = { };

     PetscInt     numCellsVentsOut = 0;
     PetscInt     numCellsVentsIn = 0;
     //PetscInt     numCellsLeak = 0;

     PetscScalar  cellArea = 0;

     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

     MPI_Comm_rank(mesh->MESH_COMM, &rank);

     DMDAVecGetArray(fda, mesh->lCent, &cent);
     DMDAVecGetArray(da, mesh->ventMarkers, &markVent);

     //printf("\n1. processor = %i, lzs = %li, lze = %li, lys = %li, lye = %li lxs = %li, lxe = %li\n", rank, lzs, lze, lys, lye, lxs, lxe);

      // Label each vent cell locations and track the number of cells in each vent.
     for (q=0; q < vents->numberOfVents; q++)
     {
          for (j=lys; j<lye; j++)
          {
           for (k=lzs; k<lze; k++)
           {
             for (i=lxs; i<lxe; i++)
             {
                 if (i == 1 || i == mx - 2 || j == 1 || j == my - 2 || k == 1 || k == mz - 2) // only completes loop if processor includes an edge of whole domain.
                 {
                     if (vents->vent[q]->face == "iLeft" && i == 1)
                     {

                         if (vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                         {

                          markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                          lNumCellsVent[q] = lNumCellsVent[q] + 1;
                         }

                     }

                     if (vents->vent[q]->face == "iRight" && i == mx - 2)
                     {
                         if (vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                         {

                          markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                          lNumCellsVent[q] = lNumCellsVent[q] + 1;
                         }

                     }

                     if (vents->vent[q]->face == "jLeft" && j == 1)
                     {
                         if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                         {

                          markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                          lNumCellsVent[q] = lNumCellsVent[q] + 1;
                         }

                     }

                     if (vents->vent[q]->face == "jRight" && j == my - 2)
                     {
                         if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                         {
                          markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                          lNumCellsVent[q] = lNumCellsVent[q] + 1;
                         }

                     }

                     if (vents->vent[q]->face == "kLeft" && k == 1)
                     {

                         if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z)
                         {

                          markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                          lNumCellsVent[q] = lNumCellsVent[q] + 1;
                         }

                      }

                     if (vents->vent[q]->face == "kRight" && k == mz - 2)
                     {
                          if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z)
                          {

                           markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                           lNumCellsVent[q] = lNumCellsVent[q] + 1;
                          }


                     }

                 }
             }
           }
          }

          //printf("\n1. processor = %i, num cells vent %i = %i\n", rank, q, lNumCellsVent[q]);
          MPI_Allreduce(&lNumCellsVent, &numCellsVent, vents->numberOfVents, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
          //printf("\n2. processor = %i, num cells vent %i = %i\n", rank, q, numCellsVent[q]);
          vents->vent[q]->nCellsVent = numCellsVent[q];
          //printf("\n3. processor = %i, num cells vent %li = %li\n", rank, q, vents->vent[q]->nCellsVent);
      }

     // find the area of each vent in m^2 and use that to set the velocity BC of each vent.
     for (q=0; q < vents->numberOfVents; q++)
     {
                if (vents->vent[q]->face == "iLeft")
                {
                    cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                    vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                    // determine the direction of the vent flow.
                    if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg")
                    {
                        vents->vent[q]->ventBCVec.x = 0;
                        vents->vent[q]->ventBCVec.y = vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea;
                        vents->vent[q]->ventBCVec.z = 0;

                    }

                    if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                    {
                        vents->vent[q]->ventBCVec.x = 0;
                        vents->vent[q]->ventBCVec.y = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                        vents->vent[q]->ventBCVec.z = 0;

                    }

                 }

                if (vents->vent[q]->face == "iRight")
                {
                     cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                     vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                     // determine the direction of the vent flow.
                     if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                     {
                         vents->vent[q]->ventBCVec.x = 0;
                         vents->vent[q]->ventBCVec.y = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                         vents->vent[q]->ventBCVec.z = 0;

                     }

                     if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                     {
                         vents->vent[q]->ventBCVec.x = 0;
                         vents->vent[q]->ventBCVec.y = (vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                         vents->vent[q]->ventBCVec.z = 0;

                     }

                  }

                if (vents->vent[q]->face == "jLeft")
                {
                     cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                     vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                     // determine the direction of the vent flow.
                     if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                     {
                         vents->vent[q]->ventBCVec.x = 0;
                         vents->vent[q]->ventBCVec.y = 0;
                         vents->vent[q]->ventBCVec.z = (vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);

                     }

                     if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                     {
                         vents->vent[q]->ventBCVec.x = 0;
                         vents->vent[q]->ventBCVec.y = 0;
                         vents->vent[q]->ventBCVec.z = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);

                     }

                 }

                if (vents->vent[q]->face == "jRight")
                {
                      cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                      vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                      // determine the direction of the vent flow.
                      if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                      {
                          vents->vent[q]->ventBCVec.x = 0;
                          vents->vent[q]->ventBCVec.y = 0;
                          vents->vent[q]->ventBCVec.z = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);

                      }

                      if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                      {
                          vents->vent[q]->ventBCVec.x = 0;
                          vents->vent[q]->ventBCVec.y = 0;
                          vents->vent[q]->ventBCVec.z = (vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);

                      }

                  }

                if (vents->vent[q]->face == "kLeft")
                {
                     cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
                     vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                     // determine the direction of the vent flow.
                     if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                     {
                         vents->vent[q]->ventBCVec.x = (vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                         vents->vent[q]->ventBCVec.y = 0;
                         vents->vent[q]->ventBCVec.z = 0;

                     }

                     if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                     {
                         vents->vent[q]->ventBCVec.x = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                         vents->vent[q]->ventBCVec.y = 0;
                         vents->vent[q]->ventBCVec.z = 0;

                     }

                 }

                if (vents->vent[q]->face == "kRight")
                {
                      cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
                      vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                      // determine the direction of the vent flow.
                      if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                      {
                          vents->vent[q]->ventBCVec.x = -(vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                          vents->vent[q]->ventBCVec.y = 0;
                          vents->vent[q]->ventBCVec.z = 0;

                      }

                      if (vents->vent[q]->dir == "outlet")// || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                      {
                          vents->vent[q]->ventBCVec.x = (vents->vent[q]->ventDesiredFlux / vents->vent[q]->ventArea);
                          vents->vent[q]->ventBCVec.y = 0;
                          vents->vent[q]->ventBCVec.z = 0;

                      }

                  }

           PetscPrintf(mesh->MESH_COMM, "vent%ld: vel x %lf, vel y %lf, vel z %lf\n", q, vents->vent[q]->ventBCVec.x, vents->vent[q]->ventBCVec.y, vents->vent[q]->ventBCVec.z);
      }

     //store the number of inlet and outlet cells
     for (q=0; q < vents->numberOfVents; q++)
     {
         if (vents->vent[q]->dir == "outlet")
         {
           numCellsVentsOut += vents->vent[q]->nCellsVent;
         }
         if (vents->vent[q]->dir == "inlet")
         {
            //printf("\nvent%li num cells = %li\n", q, vents->vent[q]->nCellsVent);
            numCellsVentsIn += vents->vent[q]->nCellsVent;
         }
         /*if (vents->vent[q]->dir == "leak")
         {
            numCellsLeak += vents->vent[q]->nCellsVent;
        }*/

      }

      vents->nCellsVentsOut = numCellsVentsOut;
      vents->nCellsVentsIn = numCellsVentsIn;
      //vents->nCellsLeak = numCellsLeak;

     DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
     DMDAVecRestoreArray(da, mesh->ventMarkers,  &markVent);
     printf("\nprocessor = %i, num cells out = %li\n", rank, vents->nCellsVentsOut);
     printf("\nprocessor = %i, num cells in = %li\n", rank, vents->nCellsVentsIn);
     //printf("\nprocessor = %i, num cells Leak = %li\n", rank, vents->nCellsLeak);


    return 0;
  }

 //***************************************************************************************************************//

 /*PetscErrorCode printVentTemp(vents_ *vents)
 {
     mesh_        *mesh = vents->access->mesh;
     DM            da   = mesh->da, fda = mesh->fda;
     teqn_        *teqn = vents->access->teqn;
     DMDALocalInfo info = mesh->info;
     PetscInt      xs   = info.xs, xe = info.xs + info.xm;
     PetscInt      ys   = info.ys, ye = info.ys + info.ym;
     PetscInt      zs   = info.zs, ze = info.zs + info.zm;
     PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

     PetscMPIInt   rank;

     PetscInt      lxs, lxe, lys, lye, lzs, lze;
     PetscInt      i, j, k;

     PetscInt      q;

     PetscReal     ***t, ***lt;
     PetscInt     ***markVent;

     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

     MPI_Comm_rank(mesh->MESH_COMM, &rank);

     DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
     DMDAVecGetArray(da, teqn->lTmprt, &lt);
     DMDAVecGetArray(da, teqn->Tmprt,  &t);

     //printf("\n1. processor = %i, lzs = %li, lze = %li, lys = %li, lye = %li lxs = %li, lxe = %li\n", rank, lzs, lze, lys, lye, lxs, lxe);

          for (j=lys; j<lye; j++)
          {
           for (k=lzs; k<lze; k++)
           {
             for (i=lxs; i<lxe; i++)
             {
                 if (markVent[k][j][i] > 0)
                 {
                     q = markVent[k][j][i] - 1;

                     if (vents->vent[q]->ventTBC == "zeroGradient")
                     {
                         printf(" ventTemp = %f   ", t[i][j][k]);
                         printf(" Local ventTemp = %f  ", lt[i][j][k]);
                     }
                 }
             }
           }
          }

     DMDAVecRestoreArray(da, mesh->ventMarkers,  &markVent);
     DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
     DMDAVecRestoreArray(da, teqn->Tmprt,  &t);
     //printf("\nprocessor = %i, num cells out = %li\n", rank, vents->nCellsVentsOut);
     //printf("\nprocessor = %i, num cells in = %li\n", rank, vents->nCellsVentsIn);
     //printf("\nprocessor = %i, num cells Leak = %li\n", rank, vents->nCellsLeak);
=======

     DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
     DMDAVecRestoreArray(da, mesh->ventMarkers,  &markVent);
>>>>>>> 2c4b51b0244c7b38e9bc9ced8e04c3606c221673


    return 0;
}*/

 //***************************************************************************************************************//
