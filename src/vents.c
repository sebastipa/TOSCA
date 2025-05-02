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

    // loop through the domains to set the vents labels
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
  //readDictWord("./vents/ventsProperties.dat", "roomPressure", &(vents->roomPressure)); add IORatio later

  // allocate memory for each vent object
  PetscMalloc(vents->numberOfVents * sizeof(ventObject*), &(vents->vent));

  char ventName[256];
  char smBC[256];
  char smBCVal[256];

  for  (PetscInt  i=0; i < vents->numberOfVents; i++)
  {
      sprintf(ventName, "vent%ld", i);

    PetscMalloc(sizeof(ventObject), &(vents->vent[i]));

    //read vent info, start with vent location and direction of flow info
    readSubDictWord("./vents/ventsProperties.dat", ventName, "face", &(vents->vent[i]->face));
    readSubDictWord("./vents/ventsProperties.dat", ventName, "dir", &(vents->vent[i]->dir));
    readSubDictWord("./vents/ventsProperties.dat", ventName, "shape", &(vents->vent[i]->shape));

    if (vents->vent[i]->shape == "rectangle")
    {
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
    }
    else if (vents->vent[i]->shape == "circle")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "diameter", &(vents->vent[i]->dia));
        readSubDictVector("./vents/ventsProperties.dat", ventName, "center", &(vents->vent[i]->center));
    }
    else
    {
        char error[512];
        sprintf(error, "unknown shape for vent %li U. Only circle and rectangle available", i);
        fatalErrorInFunction("readVentsProperties",  error);
    }

    // read U and Nut BCs
    readSubDictWord("./vents/ventsProperties.dat", ventName, "ventBC", &(vents->vent[i]->ventBC));

    readSubDictWord("./vents/ventsProperties.dat", ventName, "ventBCNut", &(vents->vent[i]->ventBCNut));

    //read fixed value flux value if needed
    if (vents->vent[i]->ventBC == "fixedValue")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventDesiredFlux", &(vents->vent[i]->ventDesiredFlux));
    }
    //read inletFunction values if needed
    else if (vents->vent[i]->ventBC == "inletFunction")
    {

            word           fileName, location, field;

            PetscPrintf(mesh->MESH_COMM, "Creating U inflow boundary data for vent%i...", i);

            // allocate memory for this patch inlet function data and init types
            PetscMalloc(sizeof(ventInletFunction), &vents->vent[i]->inletF);
            vents->vent[i]->inletF->typeU   = -1;

            // set local pointer to this inlet function type
            ventInletFunction *ifPtr = vents->vent[i]->inletF;

            location = "./vents/" + mesh->meshName + "/";
            field    = "ventsProperties.dat";
            fileName = location + field;

            readSubDictInt(fileName.c_str(), ventName, "type", &(ifPtr->typeU));

            // Inverse Fourier
            if (ifPtr->typeU == 7)
            {

                PetscPrintf(mesh->MESH_COMM, "... Reading Synthetic Turbulence Inflow for vent%i ...", i);

                //readSubDictDouble(fileName.c_str(), ventName, "Urms", &(ifPtr->Urms));
                //readSubDictVector(fileName.c_str(), ventName, "meanU", &(ifPtr->meanU));
                readSubDictDouble(fileName.c_str(), ventName, "ventDesiredFlux", &(ifPtr->desiredFlux));
                readSubDictDouble(fileName.c_str(), ventName, "TI", &(ifPtr->TI));
                readSubDictDouble(fileName.c_str(), ventName, "kolLScale", &(ifPtr->kolLScale));
                readSubDictDouble(fileName.c_str(), ventName, "intLScale", &(ifPtr->intLScale));
                readSubDictInt(fileName.c_str(), ventName, "fourierSumNum", &(ifPtr->FSumNum));
                readSubDictInt(fileName.c_str(), ventName, "iterTKE", &(ifPtr->iterTKE));
                readSubDictWord(fileName.c_str(), ventName, "genType", &(ifPtr->genType));

                ifPtr->kMax = 3.141592653589793238462643/(((mesh->bounds.Lx / mesh->KM) + (mesh->bounds.Lz / mesh->JM) + (mesh->bounds.Ly / mesh->IM))/3);
                ifPtr->kKol = 1/ifPtr->kolLScale;
                ifPtr->kEng = (9*3.141592653589793238462643*1.4256/(55*ifPtr->intLScale));
                ifPtr->kMin = ifPtr->kEng/2.75;
                ifPtr->dkn = (ifPtr->kMax - ifPtr->kMin)/(ifPtr->FSumNum-1);

                //allocate memory to vectors
                ifPtr->phaseN = (PetscReal *)malloc( sizeof(PetscReal) * ifPtr->FSumNum);
                ifPtr->uMagN = (PetscReal *)malloc( sizeof(PetscReal) * ifPtr->FSumNum);
                ifPtr->knMag = (PetscReal *)malloc( sizeof(PetscReal) * ifPtr->FSumNum);
                ifPtr->Ek = (PetscReal *)malloc( sizeof(PetscReal) * ifPtr->FSumNum);
                ifPtr->kn = (Cmpnts *)malloc( sizeof(Cmpnts) * ifPtr->FSumNum);
                ifPtr->Gn = (Cmpnts *)malloc( sizeof(Cmpnts) * ifPtr->FSumNum);

            }
            else
            {
                char error[512];
                sprintf(error, "unknown inflow profile on vent%li boundary, available profiles are:\n 7 : Sythetic turbulence (inverse fourier)", i);
                fatalErrorInFunction("SetInflowFunctions",  error);
            }

            PetscPrintf(mesh->MESH_COMM, "done\n\n");


    }
    else if (vents->vent[i]->ventBC == "zeroGradient")
    {}
    else
    {
        char error[512];
        sprintf(error, "unknown BC for vent %li U. Only fixedValue, inletFunction, and zeroGradient available", i);
        fatalErrorInFunction("readVentsProperties",  error);
    }

    //nut BC is usually ZG for vents. May change later.
    if (vents->vent[i]->ventBCNut == "fixedValue")
    {
        readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventBCNutReal", &(vents->vent[i]->ventBCNutReal));
    }
    else if (vents->vent[i]->ventBCNut == "zeroGradient")
    {

    }
    else
    {
        char error[512];
        sprintf(error, "unknown BC for vent %li nut. Only fixedValue and zeroGradient available", i);
        fatalErrorInFunction("readVentsProperties",  error);
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
    if(flags->isTeqnActive)
    {
        readSubDictWord("./vents/ventsProperties.dat", ventName, "ventTBC", &(vents->vent[i]->ventTBC));

        if (vents->vent[i]->ventTBC == "fixedValue")
        {
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventTBCVal", &(vents->vent[i]->ventTBCVal));
        }
        else if (vents->vent[i]->ventTBC == "zeroGradient")
        {

        }
        else
        {
            char error[512];
            sprintf(error, "unknown BC for vent %li temp. Only fixedValue and zeroGradient available", i);
            fatalErrorInFunction("readVentsProperties",  error);
        }
    }

    if(flags->isScalarMomentsActive)
    {
        // read sm BC info
        // allocate memory for smBC values for each vent object
        PetscMalloc(flags->isScalarMomentsActive * sizeof(ventSMObject*), &(vents->vent[i]->ventSMBC));

        readSubDictWord("./vents/ventsProperties.dat", ventName, "ventSMBC", &(vents->vent[i]->smBC));

    }

    for  (PetscInt  ii=0; ii < flags->isScalarMomentsActive; ii++)
    {
        sprintf(smBC, "ventSMBC%ld", ii);

        PetscMalloc(sizeof(ventSMObject), &(vents->vent[i]->ventSMBC[ii]));

        if (vents->vent[i]->smBC == "fixedValue")
        {
            sprintf(smBCVal, "ventSMBCVal%ld", ii);
            readSubDictDouble("./vents/ventsProperties.dat", ventName, smBCVal, &(vents->vent[i]->ventSMBC[ii]->smBCVal));
        }
        else if (vents->vent[i]->smBC == "zeroGradient")
        {

        }
        else
        {
            char error[512];
            sprintf(error, "unknown BC for vent %li SM%li. Only fixedValue and zeroGradient available", i, ii);
            fatalErrorInFunction("readVentsProperties",  error);
        }

        /*if (vents->vent[i]->ventSMBC == "intermittent")
        {
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventSMBCVal", &(vents->vent[i]->ventSMBCVal));
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventSMBCStart", &(vents->vent[i]->ventSMBCStart));
            readSubDictDouble("./vents/ventsProperties.dat", ventName, "ventSMBCEnd", &(vents->vent[i]->ventSMBCEnd));
        }*/
    }


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

     PetscScalar    lventCentX[vents->numberOfVents] = { };
     PetscScalar    lventCentY[vents->numberOfVents] = { };
     PetscScalar    lventCentZ[vents->numberOfVents] = { };

     PetscScalar    ventCentX[vents->numberOfVents] = { };
     PetscScalar    ventCentY[vents->numberOfVents] = { };
     PetscScalar    ventCentZ[vents->numberOfVents] = { };

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
         ventCentX[q] = 0.;
         ventCentY[q] = 0.;
         ventCentZ[q] = 0.;

          if (vents->vent[q]->shape == "rectangle")
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
                              ventCentX[q] += cent[k][j][i].x;
                              ventCentY[q] += cent[k][j][i].y;
                              ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if (vents->vent[q]->face == "iRight" && i == mx - 2)
                         {
                             if (vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                             {

                              markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                              lNumCellsVent[q] = lNumCellsVent[q] + 1;
                              ventCentX[q] += cent[k][j][i].x;
                              ventCentY[q] += cent[k][j][i].y;
                              ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if (vents->vent[q]->face == "jLeft" && j == 1)
                         {
                             if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                             {

                              markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                              lNumCellsVent[q] = lNumCellsVent[q] + 1;
                              ventCentX[q] += cent[k][j][i].x;
                              ventCentY[q] += cent[k][j][i].y;
                              ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if (vents->vent[q]->face == "jRight" && j == my - 2)
                         {
                             if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->xBound1 < cent[k][j][i].x && vents->vent[q]->xBound2 > cent[k][j][i].x)
                             {
                              markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                              lNumCellsVent[q] = lNumCellsVent[q] + 1;
                              ventCentX[q] += cent[k][j][i].x;
                              ventCentY[q] += cent[k][j][i].y;
                              ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if (vents->vent[q]->face == "kLeft" && k == 1)
                         {

                             if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z)
                             {

                              markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                              lNumCellsVent[q] = lNumCellsVent[q] + 1;
                              ventCentX[q] += cent[k][j][i].x;
                              ventCentY[q] += cent[k][j][i].y;
                              ventCentZ[q] += cent[k][j][i].z;

                             }

                          }

                         if (vents->vent[q]->face == "kRight" && k == mz - 2)
                         {
                              if (vents->vent[q]->yBound1 < cent[k][j][i].y && vents->vent[q]->yBound2 > cent[k][j][i].y && vents->vent[q]->zBound1 < cent[k][j][i].z && vents->vent[q]->zBound2 > cent[k][j][i].z)
                              {

                               markVent[k][j][i] = q+1; //marks vent with vent #q in ventMarkers array.
                               lNumCellsVent[q] = lNumCellsVent[q] + 1;
                               ventCentX[q] += cent[k][j][i].x;
                               ventCentY[q] += cent[k][j][i].y;
                               ventCentZ[q] += cent[k][j][i].z;

                              }


                         }

                     }
                 }
                }
               }

          }
          else if (vents->vent[q]->shape == "circle")
          {
              PetscScalar cX = vents->vent[q]->center.x;
              PetscScalar cY = vents->vent[q]->center.y;
              PetscScalar cZ = vents->vent[q]->center.z;
              PetscScalar rad = vents->vent[q]->dia/2.;

              for (j=lys; j<lye; j++)
              {
               for (k=lzs; k<lze; k++)
               {
                 for (i=lxs; i<lxe; i++)
                 {
                     if (i == 1 || i == mx - 2 || j == 1 || j == my - 2 || k == 1 || k == mz - 2) // only completes loop if processor includes an edge of whole domain.
                     {
                         if ((vents->vent[q]->face == "iLeft" && i == 1) || (vents->vent[q]->face == "iRight" && i == mx-2))
                         {
                             PetscReal DistX = cent[k][j][i].x - cX;
                             PetscReal DistZ = cent[k][j][i].z - cZ;
                             PetscReal Dist = sqrt(DistX*DistX + DistZ*DistZ);

                             if (Dist < rad)
                             {

                                 markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                                 lNumCellsVent[q] = lNumCellsVent[q] + 1;
                                 ventCentX[q] += cent[k][j][i].x;
                                 ventCentY[q] += cent[k][j][i].y;
                                 ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if ((vents->vent[q]->face == "jLeft" && j == 1) || (vents->vent[q]->face == "jRight" && j == my - 2))
                         {
                             PetscReal DistX = cent[k][j][i].x - cX;
                             PetscReal DistY = cent[k][j][i].y - cY;
                             PetscReal Dist = sqrt(DistX*DistX + DistY*DistY);

                             if (Dist < rad)
                             {

                                 markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                                 lNumCellsVent[q] = lNumCellsVent[q] + 1;
                                 ventCentX[q] += cent[k][j][i].x;
                                 ventCentY[q] += cent[k][j][i].y;
                                 ventCentZ[q] += cent[k][j][i].z;

                             }

                         }

                         if ((vents->vent[q]->face == "kLeft" && k == 1) || (vents->vent[q]->face == "kRight" && k == mz - 2))
                         {

                             PetscReal DistZ = cent[k][j][i].z - cZ;
                             PetscReal DistY = cent[k][j][i].y - cY;
                             PetscReal Dist = sqrt(DistZ*DistZ + DistY*DistY);

                             if (Dist < rad)
                             {

                                 markVent[k][j][i] = q+1; //marks vent with vent #q+1 in ventMarkers array.For vent0 q=1, Vent1 q=2.
                                 lNumCellsVent[q] = lNumCellsVent[q] + 1;
                                 ventCentX[q] += cent[k][j][i].x;
                                 ventCentY[q] += cent[k][j][i].y;
                                 ventCentZ[q] += cent[k][j][i].z;

                             }

                          }

                     }
                 }
               }
              }
          }
          else
          {
              char error[512];
              sprintf(error, "unknown shape for vent %li U. Only circle and rectangle available", i);
              fatalErrorInFunction("SetVentsAndPrint",  error);
          }

          //printf("\n1. processor = %i, num cells vent %i = %i\n", rank, q, lNumCellsVent[q]);
          MPI_Allreduce(&lNumCellsVent, &numCellsVent, vents->numberOfVents, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
          MPI_Allreduce(&ventCentX, &lventCentX, vents->numberOfVents, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
          MPI_Allreduce(&ventCentY, &lventCentY, vents->numberOfVents, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
          MPI_Allreduce(&ventCentZ, &lventCentZ, vents->numberOfVents, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
          //printf("\n2. processor = %i, num cells vent %i = %i\n", rank, q, numCellsVent[q]);
          vents->vent[q]->nCellsVent = numCellsVent[q];
          vents->vent[q]->ventCentX = lventCentX[q]/numCellsVent[q];
          vents->vent[q]->ventCentY = lventCentY[q]/numCellsVent[q];
          vents->vent[q]->ventCentZ = lventCentZ[q]/numCellsVent[q];

      }

     // find the area of each vent in m^2 and use that to set the velocity BC of each vent.
     for (q=0; q < vents->numberOfVents; q++)
     {

         if (vents->vent[q]->ventBC == "fixedValue")
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
         else if (vents->vent[q]->ventBC == "inletFunction")
         {
             ventInletFunction *ifPtr = vents->vent[q]->inletF;

             if (vents->vent[q]->face == "iLeft")
             {
                 cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                 vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                 // determine the direction of the vent flow.
                 if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg")
                 {
                     ifPtr->meanU.x = 0;
                     ifPtr->meanU.y = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                     ifPtr->meanU.z = 0;
                     ifPtr->Urms = ifPtr->meanU.y * ifPtr->TI;
                 }

                 if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                 {
                     ifPtr->meanU.x = 0;
                     ifPtr->meanU.y = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                     ifPtr->meanU.z = 0;
                     ifPtr->Urms = -ifPtr->meanU.y * ifPtr->TI;
                 }

              }

             if (vents->vent[q]->face == "iRight")
             {
                  cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                  // determine the direction of the vent flow.
                  if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                  {
                      ifPtr->meanU.x = 0;
                      ifPtr->meanU.y = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->meanU.z = 0;
                      ifPtr->Urms = -ifPtr->meanU.y * ifPtr->TI;
                  }

                  if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                  {
                      ifPtr->meanU.x = 0;
                      ifPtr->meanU.y = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->meanU.z = 0;
                      ifPtr->Urms = ifPtr->meanU.y * ifPtr->TI;
                  }

               }

             if (vents->vent[q]->face == "jLeft")
             {
                  cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                  // determine the direction of the vent flow.
                  if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                  {
                      ifPtr->meanU.x = 0;
                      ifPtr->meanU.y = 0;
                      ifPtr->meanU.z = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->Urms = ifPtr->meanU.z * ifPtr->TI;
                  }

                  if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                  {
                      ifPtr->meanU.x = 0;
                      ifPtr->meanU.y = 0;
                      ifPtr->meanU.z = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->Urms = -ifPtr->meanU.z * ifPtr->TI;
                  }

              }

             if (vents->vent[q]->face == "jRight")
             {
                   cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                   vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                   // determine the direction of the vent flow.
                   if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                   {
                       ifPtr->meanU.x = 0;
                       ifPtr->meanU.y = 0;
                       ifPtr->meanU.z = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                       ifPtr->Urms = -ifPtr->meanU.z * ifPtr->TI;
                   }

                   if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                   {
                       ifPtr->meanU.x = 0;
                       ifPtr->meanU.y = 0;
                       ifPtr->meanU.z = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                       ifPtr->Urms = ifPtr->meanU.z * ifPtr->TI;
                   }

               }

             if (vents->vent[q]->face == "kLeft")
             {
                  cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

                  // determine the direction of the vent flow.
                  if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
                  {
                      ifPtr->meanU.x = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->meanU.y = 0;
                      ifPtr->meanU.z = 0;
                      ifPtr->Urms = ifPtr->meanU.x * ifPtr->TI;
                  }

                  if (vents->vent[q]->dir == "outlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
                  {
                      ifPtr->meanU.x = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                      ifPtr->meanU.y = 0;
                      ifPtr->meanU.z = 0;
                      ifPtr->Urms = -ifPtr->meanU.x * ifPtr->TI;
                  }
         }

             if (vents->vent[q]->face == "kRight")
             {

               cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
               vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

               // determine the direction of the vent flow.
               if (vents->vent[q]->dir == "inlet") // || (vents->vent[q]->dir == "leak" && vents->roomPressure == "neg"))
               {
                   ifPtr->meanU.x = -ifPtr->desiredFlux / vents->vent[q]->ventArea;
                   ifPtr->meanU.y = 0;
                   ifPtr->meanU.z = 0;
                   ifPtr->Urms = -ifPtr->meanU.x * ifPtr->TI;
               }

               if (vents->vent[q]->dir == "outlet")// || (vents->vent[q]->dir == "leak" && vents->roomPressure == "pos"))
               {
                   ifPtr->meanU.x = ifPtr->desiredFlux / vents->vent[q]->ventArea;
                   ifPtr->meanU.y = 0;
                   ifPtr->meanU.z = 0;
                   ifPtr->Urms = ifPtr->meanU.x * ifPtr->TI;
               }

           }

             PetscPrintf(mesh->MESH_COMM, "vent%ld: vel x %lf, vel y %lf, vel z %lf\n", q, ifPtr->meanU.x, ifPtr->meanU.y, ifPtr->meanU.z);

         }
         else
         {
             if (vents->vent[q]->face == "iLeft")
             {
                 cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                 vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

             if (vents->vent[q]->face == "iRight")
             {
                  cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Lz / mesh->JM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

             if (vents->vent[q]->face == "jLeft")
             {
                  cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

             if (vents->vent[q]->face == "jRight")
             {
                   cellArea = (mesh->bounds.Lx / mesh->KM) * (mesh->bounds.Ly / mesh->IM);
                   vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

             if (vents->vent[q]->face == "kLeft")
             {
                  cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
                  vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

             if (vents->vent[q]->face == "kRight")
             {
                   cellArea = (mesh->bounds.Lz / mesh->JM) * (mesh->bounds.Ly / mesh->IM);
                   vents->vent[q]->ventArea = cellArea * vents->vent[q]->nCellsVent;

             }

         }

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
     //printf("\nprocessor = %i, num cells out = %li\n", rank, vents->nCellsVentsOut);
     //printf("\nprocessor = %i, num cells in = %li\n", rank, vents->nCellsVentsIn);
     //printf("\nprocessor = %i, num cells Leak = %li\n", rank, vents->nCellsLeak);


    return 0;
  }

 //***************************************************************************************************************//
