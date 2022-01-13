//! \file  initialField.c
//! \brief Contains initial field function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/initialField.h"

//***************************************************************************************************************//

PetscErrorCode SetInitialField(domain_ *domain)
{
  PetscInt    nDomains = domain[0].info.nDomains;
  flags_ *flags   = domain[0].access.flags;

  for(PetscInt d=0; d<nDomains; d++)
  {
    mesh_  *mesh       = domain[d].mesh;
    word   filenameU   = "./boundary/" + mesh->meshName + "/U";
    word   filenameT   = "./boundary/" + mesh->meshName + "/T";
    word   filenameNut = "./boundary/" + mesh->meshName + "/nut";

    // read the internal field U
    readDictWord(filenameU.c_str(), "internalField", &(domain[d].ueqn->initFieldType));

    // read the internal field T
    if(flags->isTeqnActive)
    {
      readDictWord(filenameT.c_str(), "internalField", &(domain[d].teqn->initFieldType));
    }

    // read the internal field nut
    if(flags->isLesActive)
    {
      readDictWord(filenameNut.c_str(), "internalField", &(domain[d].les->initFieldType));
    }

    SetInitialFieldU(domain[d].ueqn);

    if(flags->isTeqnActive)
    {
      SetInitialFieldT(domain[d].teqn);
    }

    SetInitialFieldP(domain[d].peqn);

    if(flags->isLesActive)
    {
      SetInitialFieldLES(domain[d].les);
    }

    // if readFields is on, read all the fields
    if(domain[d].ueqn->initFieldType == "readField")
    {
      PetscPrintf(mesh->MESH_COMM, "Setting initial field: %s\n\n", domain[d].ueqn->initFieldType.c_str());
      readFields(domain, domain[d].clock->startTime);
    }

    // save old fields
    VecCopy(domain[d].ueqn->Ucont, domain[d].ueqn->Ucont_o);

    if(flags->isTeqnActive)
    {
      VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
    }

  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInitialFieldU(ueqn_ *ueqn)
{
  clock_ *clock     = ueqn->access->clock;
  mesh_  *mesh      = ueqn->access->mesh;
  word   filename  = "./boundary/" + mesh->meshName + "/U";

  if(ueqn->initFieldType == "readField")
  {
    if(clock->time == 0)
    {
     char error[512];
      sprintf(error, "readField option not available at startTime 0. Use uniform field or spread the inflow\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    if(ueqn->access->flags->isTeqnActive)
    {
      if(ueqn->access->teqn->initFieldType != "readField")
      {
       char error[512];
        sprintf(error, "readField requires all fields to be set to this keyword. Temperature not set as readField\n");
        fatalErrorInFunction("SetInitialFieldU", error);
      }
    }

    if(ueqn->access->flags->isLesActive)
    {
      if(ueqn->access->les->initFieldType != "readField")
      {
       char error[512];
        sprintf(error, "readField requires all fields to be set to this keyword. nut not set as readField\n");
        fatalErrorInFunction("SetInitialFieldU", error);
      }
    }

    // all fields read together later
  }

  else if (ueqn->initFieldType == "uniform")
  {
    Cmpnts uRef;

    readSubDictVector(filename.c_str(), "uniform", "value", &(uRef));

    PetscPrintf(mesh->MESH_COMM, "Setting initial field for U: %s\n\n", ueqn->initFieldType.c_str());
    SetUniformFieldU(ueqn, uRef);
  }

  else if (ueqn->initFieldType == "ABLFlow")
  {
    if(!(ueqn->access->flags->isAblActive))
    {
     char error[512];
      sprintf(error, "activate ABL flag before setting the ABLFlow initial field\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    if(!(ueqn->access->flags->isTeqnActive))
    {
     char error[512];
      sprintf(error, "activate Teqn flag before setting the ABLFlow initial field\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    if(ueqn->access->teqn->initFieldType != "ABLFlow")
    {
     char error[512];
      sprintf(error, "Set initial field in /boundary/T to ABLFlow\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    PetscPrintf(mesh->MESH_COMM, "Setting initial field for U: %s\n\n", ueqn->initFieldType.c_str());
    SetABLInitialFlowU(ueqn);
  }

  else if (ueqn->initFieldType == "spreadInflow")
  {
    PetscPrintf(mesh->MESH_COMM, "Setting initial field for U: %s\n\n", ueqn->initFieldType.c_str());
    SpreadInletFlowU(ueqn);
  }

  else
  {
   char error[512];
    sprintf(error, "Invalid initial field keyword. Available initial fields are:\n\n        1. uniform\n        2. ABLFlow\n        3. spreadInflow\n        4. readField\n");
    fatalErrorInFunction("SetInitialFieldU", error);
  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInitialFieldT(teqn_ *teqn)
{
  clock_ *clock     = teqn->access->clock;
  mesh_  *mesh      = teqn->access->mesh;
  word   filename  = "./boundary/" + mesh->meshName + "/T";

  if(teqn->initFieldType == "readField")
  {
    if(clock->time == 0)
    {
     char error[512];
      sprintf(error, "readField option not available at startTime 0. Use uniform field or spread the inflow\n");
      fatalErrorInFunction("SetInitialFieldT", error);
    }

    if(teqn->access->flags->isTeqnActive)
    {
      if(teqn->access->ueqn->initFieldType != "readField")
      {
       char error[512];
        sprintf(error, "readField requires all fields to be set to this keyword. Velocity not set as readField\n");
        fatalErrorInFunction("SetInitialFieldT", error);
      }
    }

    // all fields read together later
  }

  else if (teqn->initFieldType == "uniform")
  {
    PetscReal tRef;

    readSubDictDouble(filename.c_str(), "uniform", "value", &(tRef));

    PetscPrintf(mesh->MESH_COMM, "Setting initial field for T: %s\n\n", teqn->initFieldType.c_str());
    SetUniformFieldT(teqn, tRef);
  }

  else if (teqn->initFieldType == "ABLFlow")
  {
    if(!(teqn->access->flags->isAblActive))
    {
     char error[512];
      sprintf(error, "activate ABL flag before setting the ABLFlow initial field\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    if(!(teqn->access->flags->isTeqnActive))
    {
     char error[512];
      sprintf(error, "activate Teqn flag before setting the ABLFlow initial field\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    if(teqn->access->ueqn->initFieldType != "ABLFlow")
    {
     char error[512];
      sprintf(error, "Set initial field in /boundary/U to ABLFlow\n");
      fatalErrorInFunction("SetInitialFieldU", error);
    }

    PetscPrintf(mesh->MESH_COMM, "Setting initial field for T: %s\n\n", teqn->initFieldType.c_str());
    SetABLInitialFlowT(teqn);
  }

  else if (teqn->initFieldType == "spreadInflow")
  {
    PetscPrintf(mesh->MESH_COMM, "Setting initial field for T: %s\n\n", teqn->initFieldType.c_str());
    SpreadInletFlowT(teqn);
  }

  else
  {
   char error[512];
    sprintf(error, "Invalid initial field keyword. Available initial fields at time=0 are:\n\n        1. uniform\n        2. ABLFlow\n        3. spreadInflow\n");
    fatalErrorInFunction("SetInitialFieldU", error);
  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInitialFieldP(peqn_ *peqn)
{
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetInitialFieldLES(les_ *les)
{
  // initialize fields to 0
  VecSet(les->lNu_t,0.);
  VecSet(les->lCs,0.);

  if(les->initFieldType == "readField")
  {
    // all fields read together later
  }

  return(0);
}

//***************************************************************************************************************//
PetscErrorCode SetUniformFieldU(ueqn_ *ueqn, Cmpnts &uRef)
{
  mesh_ *mesh = ueqn->access->mesh;
  DM            da   = mesh->da, fda = mesh->fda;
  DMDALocalInfo info = mesh->info;
  PetscInt      xs   = info.xs, xe = info.xs + info.xm;
  PetscInt      ys   = info.ys, ye = info.ys + info.ym;
  PetscInt      zs   = info.zs, ze = info.zs + info.zm;
  PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

  Cmpnts        ***ucont, ***ucat;
  Cmpnts        ***icsi, ***jeta, ***kzet;
  Cmpnts        uFaceI, uFaceJ, uFaceK;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  PetscInt      i, j, k;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);

  // loop on the internal cells and set the reference cartesian velocity
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
                  ucat[k][j][i].x = uRef.x;
                  ucat[k][j][i].y = uRef.y;
                  ucat[k][j][i].z = uRef.z;
          }
      }
  }

  DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
  DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

  UpdateCartesianBCs(ueqn);

  DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
  DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
  DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
  DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
  DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

  // interpolate contravariant velocity at internal faces
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=xs; i<lxe; i++)
          {
              // interpolate cartesian velocity at the i,j,k th face
              uFaceI.x = 0.5 * (ucat[k][j][i].x + ucat[k][j][i+1].x);
              uFaceI.y = 0.5 * (ucat[k][j][i].y + ucat[k][j][i+1].y);
              uFaceI.z = 0.5 * (ucat[k][j][i].z + ucat[k][j][i+1].z);

              ucont[k][j][i].x
              =
              (
                  uFaceI.x * icsi[k][j][i].x +
                  uFaceI.y * icsi[k][j][i].y +
                  uFaceI.z * icsi[k][j][i].z
              );

          }
      }
  }

  // loop over j-face centers
  for(k=lzs; k<lze; k++)
  {
      for(j=ys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              uFaceJ.x = 0.5 * (ucat[k][j][i].x + ucat[k][j+1][i].x);
              uFaceJ.y = 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
              uFaceJ.z = 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);

              ucont[k][j][i].y
              =
              (
                  uFaceJ.x * jeta[k][j][i].x +
                  uFaceJ.y * jeta[k][j][i].y +
                  uFaceJ.z * jeta[k][j][i].z
              );

          }
      }
  }

  // loop over k-face centers
  for(k=zs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              uFaceK.x = 0.5 * (ucat[k][j][i].x + ucat[k+1][j][i].x);
              uFaceK.y = 0.5 * (ucat[k][j][i].y + ucat[k+1][j][i].y);
              uFaceK.z = 0.5 * (ucat[k][j][i].z + ucat[k+1][j][i].z);

              ucont[k][j][i].z
              =
              (
                  uFaceK.x * kzet[k][j][i].x +
                  uFaceK.y * kzet[k][j][i].y +
                  uFaceK.z * kzet[k][j][i].z
              );

          }
      }
  }

  DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
  DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
  DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
  DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
  DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
  DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

  // update contravariant velocity at the boundaries
  UpdateContravariantBCs(ueqn);
  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetABLInitialFlowU(ueqn_ *ueqn)
{
  mesh_         *mesh = ueqn->access->mesh;
  abl_          *abl  = ueqn->access->abl;
  DM            da    = mesh->da, fda = mesh->fda;
  DMDALocalInfo info  = mesh->info;
  PetscInt      xs    = info.xs, xe = info.xs + info.xm;
  PetscInt      ys    = info.ys, ye = info.ys + info.ym;
  PetscInt      zs    = info.zs, ze = info.zs + info.zm;
  PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

  PetscReal     ***iaj, ***jaj, ***kaj;                       // cell face jacobians
  Cmpnts        ***cent;                                      // cell center coordinates
  Cmpnts        ***ucat, ***ucont;                            // cartesian and contravariant vel.
  Cmpnts        ***csi,  ***eta,  ***zet,                     // face area vectors at cell centers
                ***icsi, ***jeta, ***kzet;                    // face area vectors at cell faces

  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  PetscInt      i, j, k;

  PetscReal        uTau       = abl->uTau;
  PetscReal        hRough     = abl->hRough;
  PetscReal        uRef       = abl->uRef;
  PetscReal        hInversion = abl->hInv;
  PetscReal        vkConst    = abl->vkConst;

  PetscReal        Lx = mesh->bounds.Lx;
  PetscReal        Ly = mesh->bounds.Ly;
  PetscReal        Lz = mesh->bounds.Lz;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  DMDAVecGetArray(fda, mesh->lCsi,  &csi);
  DMDAVecGetArray(fda, mesh->lEta,  &eta);
  DMDAVecGetArray(fda, mesh->lZet,  &zet);
  DMDAVecGetArray(fda, mesh->lCent, &cent);

  DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);

  // loop over internal cells
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {

              PetscReal h = cent[k][j][i].z - mesh->bounds.zmin;
              PetscReal x = cent[k][j][i].x;
              PetscReal y = cent[k][j][i].y;
              PetscReal z = cent[k][j][i].z;

              // face area magnitude
              PetscReal faceArea
              =
              sqrt
              (
                  zet[k][j][i].x*zet[k][j][i].x +
                  zet[k][j][i].y*zet[k][j][i].y +
                  zet[k][j][i].z*zet[k][j][i].z
              );

              Cmpnts uCell;

              uCell.x = 0;
              uCell.y = 0;
              uCell.z = 0;

              if(h <= hInversion)
              {
                  uCell.z +=
                  (
                      PetscMax
                      (
                          (uTau/vkConst)*std::log(h/hRough),
                          1e-5
                      )
                  ) * faceArea;
              }
              else
              {
                  uCell.z += (uTau/vkConst)*std::log(hInversion/hRough) * faceArea;
              }

              ContravariantToCartesianPoint
              (
                  csi[k][j][i],
                  eta[k][j][i],
                  zet[k][j][i],
                  uCell,
                  &ucat[k][j][i]
              );
          }
      }
  }

  DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
  DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
  DMDAVecRestoreArray(fda, mesh->lZet,  &zet);
  DMDAVecRestoreArray(fda, mesh->lCent, &cent);

  DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
  DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

  // update cartesian BCs (sets cartesian velocity on ghost nodes)
  UpdateCartesianBCs(ueqn);

  DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
  DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
  DMDAVecGetArray(fda, mesh->lKZet,  &kzet);

  DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
  DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

  DMDAVecGetArray(da, mesh->lIAj,  &iaj);
  DMDAVecGetArray(da, mesh->lJAj,  &jaj);
  DMDAVecGetArray(da, mesh->lKAj,  &kaj);

  // loop over i-face centers
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=xs; i<lxe; i++)
          {
              // interpolate cartesian velocity at the i,j,k th face
              Cmpnts uFaceI;
              uFaceI.x = 0.5 * (ucat[k][j][i].x + ucat[k][j][i+1].x);
              uFaceI.y = 0.5 * (ucat[k][j][i].y + ucat[k][j][i+1].y);
              uFaceI.z = 0.5 * (ucat[k][j][i].z + ucat[k][j][i+1].z);

              ucont[k][j][i].x
              =
              (
                  uFaceI.x * icsi[k][j][i].x +
                  uFaceI.y * icsi[k][j][i].y +
                  uFaceI.z * icsi[k][j][i].z
              );
          }
      }
  }

  // loop over j-face centers
  for(k=lzs; k<lze; k++)
  {
      for(j=ys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              Cmpnts uFaceJ;
              uFaceJ.x = 0.5 * (ucat[k][j][i].x + ucat[k][j+1][i].x);
              uFaceJ.y = 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
              uFaceJ.z = 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);

              ucont[k][j][i].y
              =
              (
                  uFaceJ.x * jeta[k][j][i].x +
                  uFaceJ.y * jeta[k][j][i].y +
                  uFaceJ.z * jeta[k][j][i].z
              );

          }
      }
  }

  // loop over k-face centers
  for(k=zs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              Cmpnts uFaceK;
              uFaceK.x = 0.5 * (ucat[k][j][i].x + ucat[k+1][j][i].x);
              uFaceK.y = 0.5 * (ucat[k][j][i].y + ucat[k+1][j][i].y);
              uFaceK.z = 0.5 * (ucat[k][j][i].z + ucat[k+1][j][i].z);

              ucont[k][j][i].z
              =
              (
                  uFaceK.x * kzet[k][j][i].x +
                  uFaceK.y * kzet[k][j][i].y +
                  uFaceK.z * kzet[k][j][i].z
              );
          }
      }
  }

  DMDAVecRestoreArray(da, mesh->lIAj,  &iaj);
  DMDAVecRestoreArray(da, mesh->lJAj,  &jaj);
  DMDAVecRestoreArray(da, mesh->lKAj,  &kaj);

  DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
  DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
  DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

  DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
  DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
  DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

  // update contravariant BCs
  UpdateContravariantBCs(ueqn);

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SpreadInletFlowU(ueqn_ *ueqn)
{
  mesh_         *mesh = ueqn->access->mesh;
  DM            da    = mesh->da, fda = mesh->fda;
  DMDALocalInfo info  = mesh->info;
  PetscInt      xs    = info.xs, xe = info.xs + info.xm;
  PetscInt      ys    = info.ys, ye = info.ys + info.ym;
  PetscInt      zs    = info.zs, ze = info.zs + info.zm;
  PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

  Cmpnts        ***icsi, ***jeta, ***kzet; // face area vectors at cell faces
  Cmpnts        ***ucont, ***lucont;
  Cmpnts        ***ucat, ***lucat;

  PetscInt      i, j, k;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  // update cartesian BCs first (sets cartesian velocity on ghost nodes)
  UpdateCartesianBCs(ueqn);

  DMDAVecGetArray(fda, ueqn->lUcat,  &lucat);
  DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);

  if
  (
      mesh->boundaryU.kLeft == "fixedValue" ||
      mesh->boundaryU.kLeft == "inletFunction" ||
      mesh->boundaryU.kLeft == "oversetInterpolate"
  )
  {
      // allocate patch field
      std::vector<std::vector<Cmpnts>> lpatchField(my);
      std::vector<std::vector<Cmpnts>> gpatchField(my);

      for( j=0; j<my; j++)
      {
          lpatchField[j].resize(mx);
          gpatchField[j].resize(mx);
      }

      // set to zero in every processor
      for( j=0; j<my; j++)
      {
          for( i=0; i<mx; i++)
          {
              lpatchField[j][i].x = 0.0;
              lpatchField[j][i].y = 0.0;
              lpatchField[j][i].z = 0.0;

              gpatchField[j][i].x = 0.0;
              gpatchField[j][i].y = 0.0;
              gpatchField[j][i].z = 0.0;
          }
      }

      // store the cartesian velocity field
      if(zs==0)
      {
          for( j=lys; j<lye; j++)
          {
              for( i=lxs; i<lxe; i++)
              {
                  lpatchField[j][i].x = lucat[0][j][i].x;
                  lpatchField[j][i].y = lucat[0][j][i].y;
                  lpatchField[j][i].z = lucat[0][j][i].z;
              }
          }
      }

      // store the inflow data by scattering the information on all nodes
      for(j=0; j<my; j++)
      {
          MPI_Allreduce(&lpatchField[j][0], &gpatchField[j][0], mx*3, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
      }

      // loop on the internal cells and set the cartesian velocity
      for(k=lzs; k<lze; k++)
      {
          for(j=lys; j<lye; j++)
          {
              for(i=lxs; i<lxe; i++)
              {
                      ucat[k][j][i].x = gpatchField[j][i].x;
                      ucat[k][j][i].y = gpatchField[j][i].y;
                      ucat[k][j][i].z = gpatchField[j][i].z;
              }
          }
      }

      // clear the vector indices
      std::vector<std::vector<Cmpnts>> ().swap(lpatchField);
      std::vector<std::vector<Cmpnts>> ().swap(gpatchField);

  }
  else if
  (
      mesh->boundaryU.kRight == "fixedValue" ||
      mesh->boundaryU.kRight == "inletFunction"
  )
  {
      // allocate patch field
      std::vector<std::vector<Cmpnts>> lpatchField(my);
      std::vector<std::vector<Cmpnts>> gpatchField(my);

      for( j=0; j<my; j++)
      {
          lpatchField[j].resize(mx);
          gpatchField[j].resize(mx);
      }

      // set to zero in every processor
      for( j=0; j<my; j++)
      {
          for( i=0; i<mx; i++)
          {
              lpatchField[j][i].x = 0.0;
              lpatchField[j][i].y = 0.0;
              lpatchField[j][i].z = 0.0;

              gpatchField[j][i].x = 0.0;
              gpatchField[j][i].y = 0.0;
              gpatchField[j][i].z = 0.0;
          }
      }

      // store the cartesian velocity field
      if(ze==mz)
      {
          for( j=lys; j<lye; j++)
          {
              for( i=lxs; i<lxe; i++)
              {
                  lpatchField[j][i].x = lucat[mz-1][j][i].x;
                  lpatchField[j][i].y = lucat[mz-1][j][i].y;
                  lpatchField[j][i].z = lucat[mz-1][j][i].z;
              }
          }
      }

      // store the inflow data by scattering the information on all nodes
      for(j=0; j<my; j++)
      {
          MPI_Allreduce(&lpatchField[j][0], &gpatchField[j][0], mx*3, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
      }

      // loop on the internal cells and set the cartesian velocity
      for(k=lzs; k<lze; k++)
      {
          for(j=lys; j<lye; j++)
          {
              for(i=lxs; i<lxe; i++)
              {
                      ucat[k][j][i].x = gpatchField[j][i].x;
                      ucat[k][j][i].y = gpatchField[j][i].y;
                      ucat[k][j][i].z = gpatchField[j][i].z;
              }
          }
      }

      // clear the vector indices
      std::vector<std::vector<Cmpnts>> ().swap(lpatchField);
      std::vector<std::vector<Cmpnts>> ().swap(gpatchField);
  }
  else if
  (
      mesh->boundaryU.jLeft == "fixedValue"
  )
  {
      // allocate patch field
      std::vector<std::vector<Cmpnts>> lpatchField(mz);
      std::vector<std::vector<Cmpnts>> gpatchField(mz);

      for( k=0; k<mz; k++)
      {
          lpatchField[k].resize(mx);
          gpatchField[k].resize(mx);
      }

      // set to zero in every processor
      for( k=0; k<mz; k++)
      {
          for( i=0; i<mx; i++)
          {
              lpatchField[k][i].x = 0.0;
              lpatchField[k][i].y = 0.0;
              lpatchField[k][i].z = 0.0;

              gpatchField[k][i].x = 0.0;
              gpatchField[k][i].y = 0.0;
              gpatchField[k][i].z = 0.0;
          }
      }

      // store the cartesian velocity field
      if(ys==0)
      {
          for( k=lzs; k<lze; k++)
          {
              for( i=lxs; i<lxe; i++)
              {
                  lpatchField[k][i].x = lucat[k][0][i].x;
                  lpatchField[k][i].y = lucat[k][0][i].y;
                  lpatchField[k][i].z = lucat[k][0][i].z;
              }
          }
      }

      // store the inflow data by scattering the information on all nodes
      for(k=0; k<mz; k++)
      {
          MPI_Allreduce(&lpatchField[k][0], &gpatchField[k][0], mx*3, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
      }

      // loop on the internal cells and set the cartesian velocity
      for(k=lzs; k<lze; k++)
      {
          for(j=lys; j<lye; j++)
          {
              for(i=lxs; i<lxe; i++)
              {
                      ucat[k][j][i].x = gpatchField[k][i].x;
                      ucat[k][j][i].y = gpatchField[k][i].y;
                      ucat[k][j][i].z = gpatchField[k][i].z;
              }
          }
      }

      // clear the vector indices
      std::vector<std::vector<Cmpnts>> ().swap(lpatchField);
      std::vector<std::vector<Cmpnts>> ().swap(gpatchField);

  }
  else if
  (
      mesh->boundaryU.jRight == "fixedValue"
  )
  {
      // allocate patch field
      std::vector<std::vector<Cmpnts>> lpatchField(mz);
      std::vector<std::vector<Cmpnts>> gpatchField(mz);

      for( k=0; k<mz; k++)
      {
          lpatchField[k].resize(mx);
          gpatchField[k].resize(mx);
      }

      // set to zero in every processor
      for( k=0; k<mz; k++)
      {
          for( i=0; i<mx; i++)
          {
              lpatchField[k][i].x = 0.0;
              lpatchField[k][i].y = 0.0;
              lpatchField[k][i].z = 0.0;

              gpatchField[k][i].x = 0.0;
              gpatchField[k][i].y = 0.0;
              gpatchField[k][i].z = 0.0;
          }
      }

      // store the cartesian velocity field
      if(ye==my)
      {
          for( k=lzs; k<lze; k++)
          {
              for( i=lxs; i<lxe; i++)
              {
                  lpatchField[k][i].x = lucat[k][my-1][i].x;
                  lpatchField[k][i].y = lucat[k][my-1][i].y;
                  lpatchField[k][i].z = lucat[k][my-1][i].z;
              }
          }
      }

      // store the inflow data by scattering the information on all nodes
      for(k=0; k<mz; k++)
      {
          MPI_Allreduce(&lpatchField[k][0], &gpatchField[k][0], mx*3, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
      }

      // loop on the internal cells and set the cartesian velocity
      for(k=lzs; k<lze; k++)
      {
          for(j=lys; j<lye; j++)
          {
              for(i=lxs; i<lxe; i++)
              {
                      ucat[k][j][i].x = gpatchField[k][i].x;
                      ucat[k][j][i].y = gpatchField[k][i].y;
                      ucat[k][j][i].z = gpatchField[k][i].z;
              }
          }
      }

      // clear the vector indices
      std::vector<std::vector<Cmpnts>> ().swap(lpatchField);
      std::vector<std::vector<Cmpnts>> ().swap(gpatchField);

  }

  DMDAVecRestoreArray(fda, ueqn->lUcat,  &lucat);
  DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
  DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

  // update cartesian BCs again (sets cartesian velocity on ghost nodes)
  UpdateCartesianBCs(ueqn);

  DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
  DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
  DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
  DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
  DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

  // loop over i-face centers
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=xs; i<lxe; i++)
          {
              // interpolate cartesian velocity at the i,j,k th face
              Cmpnts uFaceI;
              uFaceI.x = 0.5 * (ucat[k][j][i].x + ucat[k][j][i+1].x);
              uFaceI.y = 0.5 * (ucat[k][j][i].y + ucat[k][j][i+1].y);
              uFaceI.z = 0.5 * (ucat[k][j][i].z + ucat[k][j][i+1].z);

              ucont[k][j][i].x
              =
              (
                  uFaceI.x * icsi[k][j][i].x +
                  uFaceI.y * icsi[k][j][i].y +
                  uFaceI.z * icsi[k][j][i].z
              );

          }
      }
  }

  // loop over j-face centers
  for(k=lzs; k<lze; k++)
  {
      for(j=ys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              Cmpnts uFaceJ;
              uFaceJ.x = 0.5 * (ucat[k][j][i].x + ucat[k][j+1][i].x);
              uFaceJ.y = 0.5 * (ucat[k][j][i].y + ucat[k][j+1][i].y);
              uFaceJ.z = 0.5 * (ucat[k][j][i].z + ucat[k][j+1][i].z);

              ucont[k][j][i].y
              =
              (
                  uFaceJ.x * jeta[k][j][i].x +
                  uFaceJ.y * jeta[k][j][i].y +
                  uFaceJ.z * jeta[k][j][i].z
              );

          }
      }
  }

  // loop over k-face centers
  for(k=zs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              Cmpnts uFaceK;
              uFaceK.x = 0.5 * (ucat[k][j][i].x + ucat[k+1][j][i].x);
              uFaceK.y = 0.5 * (ucat[k][j][i].y + ucat[k+1][j][i].y);
              uFaceK.z = 0.5 * (ucat[k][j][i].z + ucat[k+1][j][i].z);

              ucont[k][j][i].z
              =
              (
                  uFaceK.x * kzet[k][j][i].x +
                  uFaceK.y * kzet[k][j][i].y +
                  uFaceK.z * kzet[k][j][i].z
              );

          }
      }
  }

  DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
  DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
  DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
  DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
  DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);

  // scatter data to local values
  DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
  DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

  // update contravariant BCs
  UpdateContravariantBCs(ueqn);

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetUniformFieldT(teqn_ *teqn, PetscReal &tRef)
{
  mesh_         *mesh = teqn->access->mesh;
  DM            da    = mesh->da, fda = mesh->fda;
  DMDALocalInfo info  = mesh->info;
  PetscInt      xs    = info.xs, xe = info.xs + info.xm;
  PetscInt      ys    = info.ys, ye = info.ys + info.ym;
  PetscInt      zs    = info.zs, ze = info.zs + info.zm;
  PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

  PetscReal     ***temp;

  PetscInt      lxs, lxe, lys, lye, lzs, lze;
  PetscInt      i, j, k;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  DMDAVecGetArray(da, teqn->Tmprt,  &temp);

  // loop on the internal cells and set the reference cartesian velocity
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              temp[k][j][i] = tRef;
          }
      }
  }

  DMDAVecRestoreArray(da, teqn->Tmprt,  &temp);

  // scatter data to local values
  DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
  DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

  UpdateTemperatureBCs(teqn);

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetABLInitialFlowT(teqn_ *teqn)
{
  mesh_         *mesh = teqn->access->mesh;
  abl_          *abl  = teqn->access->abl;
  DM            da    = mesh->da, fda = mesh->fda;
  DMDALocalInfo info  = mesh->info;
  PetscInt      xs    = info.xs, xe = info.xs + info.xm;
  PetscInt      ys    = info.ys, ye = info.ys + info.ym;
  PetscInt      zs    = info.zs, ze = info.zs + info.zm;
  PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

  PetscInt           i, j, k;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  Cmpnts        ***cent;                   // cell center coordinates
  PetscReal     ***tmprt;                  // potential temperature

  PetscReal        thetaRef   = abl->tRef;
  PetscReal        gamma      = abl->gTop;
  PetscReal        deltaInv   = abl->dInv;
  PetscReal        hInv       = abl->hInv;
  PetscReal        gradInv    = abl->gInv;

  PetscReal        Lx = mesh->bounds.Lx;
  PetscReal        Ly = mesh->bounds.Ly;
  PetscReal        Lz = mesh->bounds.Lz;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  DMDAVecGetArray(da, teqn->Tmprt,  &tmprt);
  DMDAVecGetArray(fda, mesh->Cent, &cent);

  // Rampanelli and Zardi model parameters
  PetscReal smearing = abl->smear;
  PetscReal b  = smearing * gamma * deltaInv;
  PetscReal a  = gradInv - b;
  PetscReal h0 = hInv - deltaInv/2;

  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              PetscReal h = cent[k][j][i].z - mesh->bounds.zmin;

              // non dimensional height eta
              PetscReal eta = (h - hInv) / smearing / deltaInv;

              // non dimensional functions
              PetscReal f_eta = (std::tanh(eta) + 1) / 2;
              PetscReal g_eta = (std::log(2 * std::cosh(eta)) + eta) / 2;

              // potential temperature
              tmprt[k][j][i] = thetaRef + a * f_eta + b * g_eta;
          }
      }
  }

  DMDAVecRestoreArray(fda, mesh->Cent, &cent);
  DMDAVecRestoreArray(da, teqn->Tmprt,  &tmprt);

  // scatter data to local values
  DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
  DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

  UpdateTemperatureBCs(teqn);

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode SpreadInletFlowT(teqn_ *teqn)
{
  mesh_         *mesh = teqn->access->mesh;
  DM            da    = mesh->da, fda = mesh->fda;
  DMDALocalInfo info  = mesh->info;
  PetscInt      xs    = info.xs, xe = info.xs + info.xm;
  PetscInt      ys    = info.ys, ye = info.ys + info.ym;
  PetscInt      zs    = info.zs, ze = info.zs + info.zm;
  PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

  PetscReal     ***t;

  PetscInt      i, j, k;
  PetscInt      lxs, lxe, lys, lye, lzs, lze;

  lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
  lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
  lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

  // first update BC so that have field on ghost nodes
  UpdateTemperatureBCs(teqn);

  DMDAVecGetArray(da, teqn->Tmprt, &t);

  // allocate patch field
  std::vector<std::vector<PetscReal>> lpatchField(my);
  std::vector<std::vector<PetscReal>> gpatchField(my);

  for( j=0; j<my; j++)
  {
      lpatchField[j].resize(mx);
      gpatchField[j].resize(mx);
  }

  // set to zero in every processor
  for( j=0; j<my; j++)
  {
      for( i=0; i<mx; i++)
      {
          lpatchField[j][i] = 0.0;

          gpatchField[j][i] = 0.0;
      }
  }

  // store the temperature field
  if(zs==0)
  {
      for( j=lys; j<lye; j++)
      {
          for( i=lxs; i<lxe; i++)
          {
              lpatchField[j][i] = t[0][j][i];
          }
      }
  }

  // store the inflow data by scattering the information on all nodes
  for(j=0; j<my; j++)
  {
      MPI_Allreduce(&lpatchField[j][0], &gpatchField[j][0], mx, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);
  }

  // loop on the internal cells and set the temperature
  for(k=lzs; k<lze; k++)
  {
      for(j=lys; j<lye; j++)
      {
          for(i=lxs; i<lxe; i++)
          {
              t[k][j][i] = gpatchField[j][i];
          }
      }
  }

  // clear the vector indices
  std::vector<std::vector<PetscReal>> ().swap(lpatchField);
  std::vector<std::vector<PetscReal>> ().swap(gpatchField);

  DMDAVecRestoreArray(da, teqn->Tmprt,  &t);

  // scatter data to local values
  DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
  DMGlobalToLocalEnd(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

  // update cartesian BCs again (sets cartesian velocity on ghost nodes)
  UpdateTemperatureBCs(teqn);

  return(0);
}
