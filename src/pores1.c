//! \file  pores.c
//! \brief Contains pore locations and constants

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"
#include "include/wallfunctions.h"


//***************************************************************************************************************//

PetscErrorCode InitializePores(domain_ *domain)
{
  PetscInt nDomains = domain[0].info.nDomains;
  flags_ flags = domain[0].flags;

  if(flags.isPoresActive)
  {

    // loop through the domains to set the domain pores pointer
    for(PetscInt d = 0; d < nDomains; d++)
    {
      readPoresProperties(domain[d].pores);

      MPI_Barrier(domain[d].mesh->MESH_COMM);

      poreSetAndPrint(domain[d].pores);

      MPI_Barrier(domain[d].mesh->MESH_COMM);


    }

  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode readPoresProperties(pores_ *pores)
{
  PetscMPIInt rank;

  mesh_ *mesh = pores->access->mesh;
  flags_ *flags = pores->access->flags;

  MPI_Comm_rank(mesh->MESH_COMM, &rank);

  PetscPrintf(PETSC_COMM_WORLD, "\nPores initialization...\n\n");

  // read number of poress total and momentum sink constants
  readDictInt("./porousZone/porousZone.dat", "numberOfPores", &(pores->numberOfPores));
  readDictDouble("./porousZone/porousZone.dat", "porosity", &(pores->porosity));
  readDictDouble("./porousZone/porousZone.dat", "alpha", &(pores->alpha));
  readDictDouble("./porousZone/porousZone.dat", "C2", &(pores->C2));
  readDictDouble("./porousZone/porousZone.dat", "filtrationEfficiency", &(pores->filtrationEfficiency));

  // allocate memory for each poreZone
  PetscMalloc(pores->numberOfPores * sizeof(porousZone*), &(pores->poreZone));

  char poreName[256];

  for  (PetscInt  i=0; i < pores->numberOfPores; i++)
  {
      sprintf(poreName, "poreZone%ld", i);

    PetscMalloc(sizeof(porousZone), &(pores->poreZone[i]));

    //read poreZone bounds

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "xBound1", &(pores->poreZone[i]->xBound1));

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "xBound2", &(pores->poreZone[i]->xBound2));

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "yBound1", &(pores->poreZone[i]->yBound1));

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "yBound2", &(pores->poreZone[i]->yBound2));

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "zBound1", &(pores->poreZone[i]->zBound1));

        readSubDictDouble("./porousZone/porousZone.dat", poreName, "zBound2", &(pores->poreZone[i]->zBound2));

  }

  return 0;
}

 //***************************************************************************************************************//

 PetscErrorCode poreSetAndPrint(pores_ *pores)
 {
     mesh_        *mesh = pores->access->mesh;
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

     PetscReal     ***markPore;


     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

     MPI_Comm_rank(mesh->MESH_COMM, &rank);

     DMDAVecGetArray(fda, mesh->lCent, &cent);
     DMDAVecGetArray(da, mesh->poreMarkers, &markPore);

     // Label each pore cell locations and track the number of cells in each vent.
     for (q=0; q < pores->numberOfPores; q++)
     {

         for (k=lzs; k<lze; k++)
         {
             for (j=lys; j<lye; j++)
             {
                 for (i=lxs; i<lxe; i++)
                 {

                         if (

                         (pores->poreZone[q]->xBound1 < cent[k][j][i].x && pores->poreZone[q]->xBound2 > cent[k][j][i].x) &&
                         (pores->poreZone[q]->yBound1 < cent[k][j][i].y && pores->poreZone[q]->yBound2 > cent[k][j][i].y) &&
                         (pores->poreZone[q]->zBound1 < cent[k][j][i].z && pores->poreZone[q]->zBound2 > cent[k][j][i].z)

                            )

                         {

                             markPore[k][j][i] = pores->porosity; //marks porezone with porosity value which will be less than 1. Greater than 0.

                         }
                 }
             }
         }
     }



     DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
     DMDAVecRestoreArray(da, mesh->poreMarkers,  &markPore);


    return 0;
  }

 //***************************************************************************************************************//
