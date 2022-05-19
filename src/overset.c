//! \file  overset.c
//! \brief overset function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/initialField.h"
#include "include/overset.h"

//***************************************************************************************************************//

PetscErrorCode InitializeOverset(domain_ *domain)
{

  PetscInt nDomains = domain[0].info.nDomains;
  flags_ flags = domain[0].flags;

  for(PetscInt d = 0; d < nDomains; d++)
  {
    if(domain[d].os!=NULL)
    {
      overset_ *os   = domain[d].os;
      mesh_    *mesh = domain[d].mesh;

      // only for overset meshes. Base mesh is excluded
      if(os->parentMeshId[0] != -1)
      {
        PetscPrintf(mesh->MESH_COMM, "\nOverset initialization...");

        // read the overset parameters for the current domain
        readOversetProperties(os);

        // Acceptor cells of the background and overset mesh for current user
        createAcceptorCell(os);

        // Find donors and setup acceptor donor connectivity
        for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
        {

            findClosestDonor(domain[os->parentMeshId[pi]].mesh, mesh);

            if(os->interpolationType != "trilinear")
            {
                acellDcellConnectivity(domain[os->parentMeshId[pi]].mesh, mesh);

                if(os->interpolationType == "LS1")
                {
                  getLSWeights(domain[os->parentMeshId[pi]].mesh, mesh);
                }
                else if (os->interpolationType == "LS2")
                {
                  getLSWeights_2nd(domain[os->parentMeshId[pi]].mesh, mesh);
                }
                else if (os->interpolationType == "LS3")
                {
                  getLSWeights_3rd(domain[os->parentMeshId[pi]].mesh, mesh);
                }
            }

            // initial interpolation
            if(os->interpolationType == "inverseDist")
            {
              interpolateACellInvD(domain[os->parentMeshId[pi]].mesh, mesh);
            }
            else if(os->interpolationType == "LS1" || os->interpolationType == "LS2" || os->interpolationType == "LS3")
            {
              interpolateACellLS(domain[os->parentMeshId[pi]].mesh, mesh);
            }
            else  //trilinear interpolation by default
            {
              interpolateACellTrilinear(domain[os->parentMeshId[pi]].mesh, mesh);
            }

        }

        PetscPrintf(mesh->MESH_COMM, "done.\n\n");

        // initialize the field after interpolation
        SetInitialFieldU(domain[d].ueqn);

        if(flags.isTeqnActive)
        {
          SetInitialFieldT(domain[d].teqn);
        }

        SetInitialFieldP(domain[d].peqn);

        if(flags.isLesActive)
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

        if(flags.isTeqnActive)
        {
          VecCopy(domain[d].teqn->Tmprt, domain[d].teqn->Tmprt_o);
        }

      }

    }

  }

  return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateOversetInterpolation(domain_ *domain)
{
  PetscInt nDomains = domain[0].info.nDomains;
  PetscReal     ts,te;

  PetscTime(&ts);

  for(PetscInt d = 0; d < nDomains; d++)
  {
    overset_ *os   = domain[d].os;
    mesh_    *mesh = domain[d].mesh;

    // only for overset meshes. Base mesh is excluded
    if(os->parentMeshId[0] != -1)
    {
      for (PetscInt pi = 0; pi < os->parentMeshId.size(); pi++)
      {
        if(os->interpolationType == "inverseDist")
        {
          interpolateACellInvD(domain[os->parentMeshId[pi]].mesh, mesh);
        }
        else if(os->interpolationType == "LS1" || os->interpolationType == "LS2" || os->interpolationType == "LS3")
        {
          interpolateACellLS(domain[os->parentMeshId[pi]].mesh, mesh);
        }
        else  //trilinear interpolation by default
        {
          interpolateACellTrilinear(domain[os->parentMeshId[pi]].mesh, mesh);
        }

        MPI_Barrier(mesh->MESH_COMM);
      }

      //Update mesh if dynamic overset
      if(os->dynamicOverset)
      {
        // update the acceptor cell co-ordinates for the next time step
            updateAcceptorCoordinates(os);

        // if background mesh processors intersecting with the overset mesh have changed
        if(os->procChange)
        {
          // update the intersecting processors
            // updateIntersectingProcessors(user);

            // updateDonorCells(user);
        }
        else
        {
            // updateDonorCells(user);
        }
      }

      // update bc at boundaries not interpolated
      // update cartesian BC
      UpdateCartesianBCs(domain[d].ueqn);

      // update contravariant BC
      UpdateContravariantBCs(domain[d].ueqn);
    }

  }

  PetscTime(&te);
  PetscPrintf(PETSC_COMM_WORLD,"Overset Interpolation Elapsed Time = %lf\n", te-ts);

    return 0;
}

//***************************************************************************************************************//

// overset simulation properties
PetscErrorCode readOversetProperties(overset_ *os)
{
  // to set dynamic overset on
  readDictInt("Overset/OversetInput.dat", "dynamicOverset", &(os->dynamicOverset));

  // read the interpolation type
  new(&(os->interpolationType)) word{};
  readDictWord("Overset/OversetInput.dat", "interpolationType", &(os->interpolationType));

  // search radius size = cell size X cell factors
  if( (os->interpolationType == "LS1") || (os->interpolationType == "LS2") || (os->interpolationType == "LS3"))
  {
    readDictDouble("Overset/OversetInput.dat", "cellFactor", &(os->cellFactor));
  }

  //allocate memory for the overset motion struc if dynamicOverset is on
  if(os->dynamicOverset)
  {
    PetscMalloc(sizeof(oversetMotion), &(os->oMotion));

    oversetMotion *osetMotion = os->oMotion;

    // read the prescribed motion switch
    readSubDictInt("Overset/OversetInput.dat", "oversetMotion", "setMotion", &(osetMotion->setMotion));

    if(osetMotion->setMotion)
    {
      // read the interpolation type
      readSubDictWord("Overset/OversetInput.dat", "oversetMotion", "motionType", &(osetMotion->motionType));

      if(osetMotion->motionType == "Translation")
      {
        readSubDictVector("Overset/OversetInput.dat", "oversetMotion", "prescribedVel", &(osetMotion->prescribedVel));
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

PetscErrorCode interpolateACellTrilinear(mesh_ *meshP, mesh_ *mesh){

    overset_         *os    = mesh->access->os;
    ueqn_            *ueqn  = mesh->access->ueqn;
    ueqn_            *ueqnP = meshP->access->ueqn;
    teqn_            *teqn  = mesh->access->teqn;
    teqn_            *teqnP = meshP->access->teqn;
    flags_           *flags = mesh->access->flags;
    DM               da1    = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1  = mesh->info;
    DM               da0    = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0  = meshP->info;

    PetscInt         xs = info1.xs, xe = info1.xs + info1.xm;
    PetscInt         ys = info1.ys, ye = info1.ys + info1.ym;
    PetscInt         zs = info1.zs, ze = info1.zs + info1.zm;
    PetscInt         mx = info1.mx, my = info1.my, mz = info1.mz;

    PetscInt         i, j, k, ic, kc, jc, b, n, m;

    Cmpnts           ***lucat0, ***ucat1, ***cent, ucart;
    PetscReal        ***ltemp0, ***temp1, Temp;

    Cmpnts           pCoor;

    PetscMPIInt      rank, size, rankP, sizeP;
    PetscInt         sum_ind1 = 0;

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    MPI_Comm_size(meshP->MESH_COMM, &sizeP);
    MPI_Comm_rank(meshP->MESH_COMM, &rankP);

    std::vector<Acell> aCell = os->aCell;
    std::vector<Dcell> dCell = os->closestDonor;
    std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMat;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;

    DMDAVecGetArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecGetArray(fda1, ueqn->Ucat, &ucat1);
    DMDAVecGetArray(fda0, meshP->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecGetArray(da1, teqn->Tmprt, &temp1);
    }

    // loop through the ranks
    for(n = 0; n < size; n++)
    {

        if(NumAcellPerProc[n]!=0)
        {

            if(AcellProcMat[n][rankP] !=MPI_UNDEFINED)
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

                    if (rankP == dCell[b].rank)
                    {
                        ic = dCell[b].indi;
                        jc = dCell[b].indj;
                        kc = dCell[b].indk;

                        vectorPointLocalVolumeInterpolation
                        (
                                meshP,
                                pCoor.x, pCoor.y, pCoor.z,
                                ic, jc, kc,
                                cent,
                                lucat0,
                                ucart
                        );

                        if (flags->isTeqnActive)
                        {
                            scalarPointLocalVolumeInterpolation
                            (
                                    meshP,
                                    pCoor.x, pCoor.y, pCoor.z,
                                    ic, jc, kc,
                                    cent,
                                    ltemp0,
                                    Temp
                            );
                        }

                        MPI_Send(&ucart, 3, MPIU_REAL, aCell[b].rank, 0, meshP->MESH_COMM);
                        MPI_Send(&Temp, 1, MPIU_REAL, aCell[b].rank, 1, meshP->MESH_COMM);

                    }

                    if (rank == aCell[b].rank)
                    {

                      MPI_Recv(&ucart, 3, MPIU_REAL, dCell[b].rank, 0, meshP->MESH_COMM, MPI_STATUS_IGNORE);
                      MPI_Recv(&Temp, 1, MPIU_REAL, dCell[b].rank, 1, meshP->MESH_COMM, MPI_STATUS_IGNORE);

                      if(aCell[b].face == 0)
                      {
                          ucat1[k][j][i].x = ucart.x;
                          ucat1[k][j][i].y = ucart.y;
                          ucat1[k][j][i].z = ucart.z;

                          if (flags->isTeqnActive)
                          {
                              temp1[k][j][i] = Temp;
                          }
                      }
                      else
                      {
                          oversetContravariantBC(mesh, i, j, k, ucart, aCell[b].face);
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

    DMDAVecRestoreArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecRestoreArray(fda1, ueqn->Ucat, &ucat1);
    DMDAVecRestoreArray(fda0, meshP->lCent, &cent);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecRestoreArray(da1, teqn->Tmprt, &temp1);

        DMGlobalToLocalBegin(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    }

    DMGlobalToLocalBegin(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    DMGlobalToLocalBegin(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return 0;
}

//***************************************************************************************************************//

// split the global aCell vector based on processors. For the cells within each processor
// interpolate from the donor cell connectivity and using the processor matrix (to know which processors have the donor cells)

PetscErrorCode interpolateACellLS(mesh_ *meshP, mesh_ *mesh)
{
    overset_         *os    = mesh->access->os;
    ueqn_            *ueqn  = mesh->access->ueqn;
    ueqn_            *ueqnP = meshP->access->ueqn;
    teqn_            *teqn  = mesh->access->teqn;
    teqn_            *teqnP = meshP->access->teqn;
    flags_           *flags = mesh->access->flags;
    DM               da1    = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1  = mesh->info;
    DM               da0    = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0  = meshP->info;

    PetscInt         xs = info1.xs, xe = info1.xs + info1.xm;
    PetscInt         ys = info1.ys, ye = info1.ys + info1.ym;
    PetscInt         zs = info1.zs, ze = info1.zs + info1.zm;
    PetscInt         mx = info1.mx, my = info1.my, mz = info1.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, b, m, n;
    PetscInt         ii, jj, kk;

    Cmpnts           ***lucat0, ***ucat1, lucart, gucart;

    PetscReal        dist, ds;
    PetscReal        ***nvert, ***ltemp0, ***temp1, lT, gT;
    PetscMPIInt      rank, size, rankP, sizeP;
    PetscInt         sum_ind1 = 0;

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    MPI_Comm_size(meshP->MESH_COMM, &sizeP);
    MPI_Comm_rank(meshP->MESH_COMM, &rankP);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    std::vector<Acell> aCell = os->aCell;
    std::vector<std::vector<Dcell>> dCell = os->dCell;
    std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMat;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;
    std::vector<std::vector<PetscReal>> DWeights = os->DWeights;

    DMDAVecGetArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecGetArray(fda1, ueqn->Ucat, &ucat1);

    if (flags->isTeqnActive)
    {
        DMDAVecGetArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecGetArray(da1, teqn->Tmprt, &temp1);
    }

    // to set the local rank within the new communicator
    std::vector<PetscInt> local_rank_Vec;
    PetscInt aCellLocalrank = 0;

    // loop through the ranks
    for(n = 0; n < size; n++){

        for(m = 0; m < sizeP; m++){
            if(AcellProcMat[n][m] == 1)
                local_rank_Vec.push_back(m);
        }

        // find the rank of aCell cells within the new communicator
        for(m = 0; m < local_rank_Vec.size(); m++){
            if(n == local_rank_Vec[m])
                aCellLocalrank = m;
        }

        std::vector<PetscInt> ().swap(local_rank_Vec);

        // proceed only if there are non zero acceptor cells in processor n
        if(NumAcellPerProc[n]!=0){

            // proceed only if the current processor has donor cells to processor n
            if(AcellProcMat[n][rankP] !=MPI_UNDEFINED){

                // now data is exchanged only between the acceptor donor processors through the new communicator
                // loop through the aCell cells of a given processor n
                for(b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++){

                    lucart.x = 0.; gucart.x = 0.;
                    lucart.y = 0.; gucart.y = 0.;
                    lucart.z = 0.; gucart.z = 0.;
                    lT = 0;        gT = 0;

                    // aCell cell index
                    i = aCell[b].indi;
                    j = aCell[b].indj;
                    k = aCell[b].indk;

                    // loop through the donor cells of the current acceptor
                    for(m = 0; m < dCell[b].size(); m++){

                        kk = dCell[b][m].indk;
                        jj = dCell[b][m].indj;
                        ii = dCell[b][m].indi;

                        lucart.x += DWeights[b][m] * lucat0[kk][jj][ii].x;
                        lucart.y += DWeights[b][m] * lucat0[kk][jj][ii].y;
                        lucart.z += DWeights[b][m] * lucat0[kk][jj][ii].z;

                        if (flags->isTeqnActive)
                        {
                            lT += DWeights[b][m] * ltemp0[kk][jj][ii];
                        }

                   }

                    // reduce the contribution of all valid processors to the local rank aCellLocalrank within the new communicator
                    MPI_Reduce(&lucart, &gucart, 3, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
                    MPI_Reduce(&lT, &gT, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);

                    if (rank == aCell[b].rank){

                        if(aCell[b].face == 0)
                        {
                            ucat1[k][j][i].x = gucart.x;
                            ucat1[k][j][i].y = gucart.y;
                            ucat1[k][j][i].z = gucart.z;

                            if (flags->isTeqnActive)
                            {
                                temp1[k][j][i] = gT;
                            }
                        }
                        else
                        {
                            oversetContravariantBC(mesh, i, j, k, gucart, aCell[b].face);
                        }

                    }

                }
            }

            sum_ind1 +=NumAcellPerProc[n];

        }

    }

    std::vector<Acell> ().swap(aCell);
    std::vector<std::vector<Dcell>> ().swap(dCell);
    std::vector<std::vector<PetscReal>> ().swap(DWeights);
    std::vector<std::vector<PetscInt>> ().swap(AcellProcMat);
    std::vector<PetscInt> ().swap(NumAcellPerProc);

    DMDAVecRestoreArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecRestoreArray(fda1, ueqn->Ucat, &ucat1);

    if (flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecRestoreArray(da1, teqn->Tmprt, &temp1);

        DMGlobalToLocalBegin(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    }

    DMGlobalToLocalBegin(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    DMGlobalToLocalBegin(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode interpolateACellInvD(mesh_ *meshP, mesh_ *mesh)
{
    overset_         *os    = mesh->access->os;
    ueqn_            *ueqn  = mesh->access->ueqn;
    ueqn_            *ueqnP = meshP->access->ueqn;
    teqn_            *teqn  = mesh->access->teqn;
    teqn_            *teqnP = meshP->access->teqn;

    DM               da1    = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1  = mesh->info;
    DM               da0    = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0  = meshP->info;

    PetscInt         xs = info1.xs, xe = info1.xs + info1.xm;
    PetscInt         ys = info1.ys, ye = info1.ys + info1.ym;
    PetscInt         zs = info1.zs, ze = info1.zs + info1.zm;
    PetscInt         mx = info1.mx, my = info1.my, mz = info1.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, b, m, n;
    PetscInt         ii, jj, kk;
    PetscInt         ip, im, jp, jm, kp, km;

    Cmpnts           ***lucat0, ***ucat1, lucart, gucart;

    PetscReal        dist;
    PetscReal        ***ltemp0, ***temp1, lT, gT;

    PetscMPIInt      rank, size, rankP, sizeP;
    PetscInt         sum_ind1 = 0;

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    MPI_Comm_size(meshP->MESH_COMM, &sizeP);
    MPI_Comm_rank(meshP->MESH_COMM, &rankP);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    std::vector<Acell> aCell = os->aCell;
    std::vector<std::vector<Dcell>> dCell = os->dCell;
    std::vector<std::vector<PetscInt>> AcellProcMat = os->AcellProcMat;
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;
    std::vector<std::vector<PetscReal>> DWeights = os->DWeights;

    DMDAVecGetArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecGetArray(fda1, ueqn->Ucat, &ucat1);

    if (mesh->access->flags->isTeqnActive)
    {
        DMDAVecGetArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecGetArray(da1, teqn->Tmprt, &temp1);
    }

    // to set the local rank within the new communicator
    std::vector<PetscInt> local_rank_Vec;
    PetscInt aCellLocalrank = 0;

    // loop through the ranks
    for(n = 0; n < size; n++){

        for(m = 0; m < sizeP; m++){
            if(AcellProcMat[n][m] == 1)
                local_rank_Vec.push_back(m);
        }

        // find the rank of aCell cells within the new communicator
        for(m = 0; m < local_rank_Vec.size(); m++){
            if(n == local_rank_Vec[m])
                aCellLocalrank = m;
        }

        std::vector<PetscInt> ().swap(local_rank_Vec);

        if(NumAcellPerProc[n]!=0){

            if(AcellProcMat[n][rankP] !=MPI_UNDEFINED){

                // loop through the aCell cells of a given processor n
                for(b = sum_ind1; b < sum_ind1 + NumAcellPerProc[n]; b++){

                    PetscReal   lsumwt = 0.;
                    PetscReal   gsumwt = 0.;

                    lucart.x = 0.; gucart.x = 0.;
                    lucart.y = 0.; gucart.y = 0.;
                    lucart.z = 0.; gucart.z = 0.;
                    lT = 0;        gT = 0;

                    // aCell cell index
                    i = aCell[b].indi;
                    j = aCell[b].indj;
                    k = aCell[b].indk;

                    // loop through the donor cells within the current processor - gnumdCell[b][rank]
                    for(m = 0; m < dCell[b].size(); m++){

                        kk = dCell[b][m].indk;
                        jj = dCell[b][m].indj;
                        ii = dCell[b][m].indi;

                        dist = dCell[b][m].dist2p;

                        lucart.x += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].x;
                        lucart.y += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].y;
                        lucart.z += (1.0/(PetscMax(dist, 1e-10))) * lucat0[kk][jj][ii].z;

                        if (mesh->access->flags->isTeqnActive)
                        {
                            lT += (1.0/(PetscMax(dist, 1e-10))) * ltemp0[kk][jj][ii];
                        }

                        lsumwt += 1.0/(PetscMax(dist, 1e-10));

                    }

                    // reduce the contribution of all valid processors to the local rank aCellLocalrank within the new communicator
                    MPI_Reduce(&lucart, &gucart, 3, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
                    MPI_Reduce(&lT, &gT, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);
                    MPI_Reduce(&lsumwt, &gsumwt, 1, MPIU_REAL, MPIU_SUM, aCellLocalrank, os->oset_comm[n]);

                    gucart.x /= gsumwt;
                    gucart.y /= gsumwt;
                    gucart.z /= gsumwt;
                    gT /= gsumwt;

                    if (rank == aCell[b].rank){

                       if(aCell[b].face == 0)
                       {
                         ucat1[k][j][i].x = gucart.x;
                         ucat1[k][j][i].y = gucart.y;
                         ucat1[k][j][i].z = gucart.z;

                         if (mesh->access->flags->isTeqnActive)
                         {
                             temp1[k][j][i] = gT;
                         }
                       }
                       else
                       {
                         oversetContravariantBC(mesh, i, j, k, gucart, aCell[b].face);
                       }

                    }

                }
            }

            sum_ind1 +=NumAcellPerProc[n];

        }

    }

    std::vector<Acell> ().swap(aCell);
    std::vector<std::vector<Dcell>> ().swap(dCell);
    std::vector<std::vector<PetscReal>> ().swap(DWeights);
    std::vector<std::vector<PetscInt>> ().swap(AcellProcMat);
    std::vector<PetscInt> ().swap(NumAcellPerProc);

    DMDAVecRestoreArray(fda0, ueqnP->lUcat, &lucat0);
    DMDAVecRestoreArray(fda1, ueqn->Ucat, &ucat1);

    if (mesh->access->flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da0, teqnP->lTmprt, &ltemp0);
        DMDAVecRestoreArray(da1, teqn->Tmprt, &temp1);

        DMGlobalToLocalBegin(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        DMGlobalToLocalEnd(da1, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    }

    DMGlobalToLocalBegin(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda1, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    DMGlobalToLocalBegin(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda1, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return 0;
}

//***************************************************************************************************************//

//! \brief Create the list of acceptor cells which will be interpolated
PetscErrorCode createAcceptorCell(overset_ *os){

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

    PetscMPIInt   rank, size;
    PetscInt      count = 0, sum_ind = 0, gcount = 0;

    Cmpnts        ***cent, ***csi, ***eta, ***zet;

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::vector <PetscInt>   lnum(size);                                                // local number of Acceptor cells in each processor
    std::vector <PetscInt>   gnum(size);                                                // scattered to find the Acceptor cells in all the processors
    std::vector <Acell> laCell;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da,  mesh->lAj,   &aj);
    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);

    for (b=0; b < size; b++){
        lnum[b] = 0;
        gnum[b] = 0;
    }

    //Calculate the local number of aCell cells in each processor. to scatter the local vectors correctly

    for (k=zs; k<ze; k++)
        for (j=ys; j<ye; j++)
            for (i=xs; i<xe; i++){

                // exclude corner ghost cells
                if(isOnCornerCellCenters(i, j, k, info)) continue;

                // count the ghost cells as the acceptor cells at the boundary
                if((k == 0) || (k == mz-1) || (j==0) || (j==my-1) || (i == 0) || (i == mx-1))
                {
                    count++;
                }

                // count the faces for the acceptor cells at the boundary of first internal cells

                if(k == 1 && (j!=0 && i!=0 && j!=my-1 && i!=mx-1)){
                    count++;
                }

                if(k == mz-2 && (j!=0 && i!=0 && j!=my-1 && i!=mx-1)){
                    count++;
                }

                if(j==1 && (k!=0 && i!=0 && k!=mz-1 && i!=mx-1)){
                    count++;
                }

                if(j==my-2 && (k!=0 && i!=0 && k!=mz-1 && i!=mx-1)){
                    count++;
                }

                if(i == 1 && (k!=0 && j!=0 && k!=mz-1 && j!=my-1)){
                    count++;
                }

                if(i == mx-2 && (k!=0 && j!=0 && k!=mz-1 && j!=my-1)){
                    count++;
                }

            }

    lnum[rank] = count;

    MPI_Allreduce(&lnum[0], &gnum[0], size, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    os->NumAcellPerProc = gnum; // store the number of acceptor cells in each processor for interpolation later

    // sum_ind is the starting cellIds of acceptor cells in each processor
    for (b=0; b< rank; b++){
        sum_ind +=gnum[b];
    }

    for (b=0; b< size; b++){
        gcount +=gnum[b];   // total number of acceptor cells
    }

    laCell.resize(gcount);
    os->aCell.resize(gcount);

    for (b=0; b< gcount; b++){
        laCell[b].indi = 0; laCell[b].indj = 0; laCell[b].indk = 0;
        laCell[b].coorx = 0.; laCell[b].coory = 0.; laCell[b].coorz = 0.;
        laCell[b].rank = 0;
        laCell[b].cell_size = 0;
        laCell[b].face = 0;

        os->aCell[b].indi = 0; os->aCell[b].indj = 0; os->aCell[b].indk = 0;
        os->aCell[b].coorx = 0.; os->aCell[b].coory = 0.; os->aCell[b].coorz = 0.;
        os->aCell[b].rank = 0;
        os->aCell[b].cell_size = 0;
        os->aCell[b].face = 0;
    }

    PetscReal areaI, areaJ, areaK ;
    PetscReal lenI, lenJ, lenK;

    for (k=zs; k<ze; k++)
        for (j=ys; j<ye; j++)
            for (i=xs; i<xe; i++)
            {

                // exclude corner ghost cells as cell centers not defined
                if(isOnCornerCellCenters(i, j, k, info)) continue;

                if((k == 0) || (k == mz-1) || (j==0) || (j==my-1) || (i == 0) || (i == mx-1))
                {
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x;
                    laCell[sum_ind].coory = cent[k][j][i].y;
                    laCell[sum_ind].coorz = cent[k][j][i].z;

                    laCell[sum_ind].rank = rank;

                    sum_ind++;
                }

                areaI = sqrt
                        (
                                csi[k][j][i].x*csi[k][j][i].x +
                                csi[k][j][i].y*csi[k][j][i].y +
                                csi[k][j][i].z*csi[k][j][i].z
                        );

                areaJ = sqrt
                        (
                            eta[k][j][i].x*eta[k][j][i].x +
                            eta[k][j][i].y*eta[k][j][i].y +
                            eta[k][j][i].z*eta[k][j][i].z
                        );

                areaK = sqrt
                        (
                            zet[k][j][i].x*zet[k][j][i].x +
                            zet[k][j][i].y*zet[k][j][i].y +
                            zet[k][j][i].z*zet[k][j][i].z
                        );

                lenI = 1.0/aj[k][j][i]/areaI;
                lenJ = 1.0/aj[k][j][i]/areaJ;
                lenK = 1.0/aj[k][j][i]/areaK;

                if(k == 1 && (j!=0 && i!=0 && j!=my-1 && i!=mx-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x - 0.5*lenK;
                    laCell[sum_ind].coory = cent[k][j][i].y;
                    laCell[sum_ind].coorz = cent[k][j][i].z;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 5;
                    sum_ind++;
                }

                if(k == mz-2 && (j!=0 && i!=0 && j!=my-1 && i!=mx-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x + 0.5*lenK;
                    laCell[sum_ind].coory = cent[k][j][i].y;
                    laCell[sum_ind].coorz = cent[k][j][i].z;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 6;
                    sum_ind++;
                }

                if(j==1 && (k!=0 && i!=0 && k!=mz-1 && i!=mx-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x;
                    laCell[sum_ind].coory = cent[k][j][i].y;
                    laCell[sum_ind].coorz = cent[k][j][i].z - 0.5*lenJ;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 3;
                    sum_ind++;
                }

                if(j==my-2 && (k!=0 && i!=0 && k!=mz-1 && i!=mx-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x;
                    laCell[sum_ind].coory = cent[k][j][i].y;
                    laCell[sum_ind].coorz = cent[k][j][i].z + 0.5*lenJ;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 4;
                    sum_ind++;
                }

                if(i == 1 && (k!=0 && j!=0 && k!=mz-1 && j!=my-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x;
                    laCell[sum_ind].coory = cent[k][j][i].y - 0.5*lenI;
                    laCell[sum_ind].coorz = cent[k][j][i].z;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 1;
                    sum_ind++;
                }

                if(i == mx-2 && (k!=0 && j!=0 && k!=mz-1 && j!=my-1)){
                    laCell[sum_ind].indk = k;
                    laCell[sum_ind].indj = j;
                    laCell[sum_ind].indi = i;

                    laCell[sum_ind].coorx = cent[k][j][i].x;
                    laCell[sum_ind].coory = cent[k][j][i].y + 0.5*lenI;
                    laCell[sum_ind].coorz = cent[k][j][i].z;

                    laCell[sum_ind].rank = rank;

                    laCell[sum_ind].face = 2;
                    sum_ind++;
                }

            }

    // define the datatype and operation for the mpi reduce operation
    MPI_Datatype mpi_Acell;
    defineStruct_Acell(&mpi_Acell);

    MPI_Op       sumstruct;
    MPI_Op_create(sum_struct_Acell, gcount, &sumstruct);

    // all reduce so that the acceptor cell list is available to all the processors
    MPI_Allreduce(&laCell[0], &(os->aCell[0]), gcount, mpi_Acell, sumstruct, mesh->MESH_COMM);

   // for (b = 0; b < gcount; b++){
   //     PetscPrintf(mesh->MESH_COMM, "rank = %ld aCell[%ld] %ld %ld %ld %lf %lf %lf\n", rank, b, os->aCell[b].indi, os->aCell[b].indj, os->aCell[b].indk, os->aCell[b].coorx, os->aCell[b].coory, os->aCell[b].coorz);
   // }

    std::vector<Acell> ().swap(laCell);
    std::vector<PetscInt> ().swap(lnum);
    std::vector<PetscInt> ().swap(gnum);

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode acellDcellConnectivity(mesh_ *meshP, mesh_ *mesh)
{

    overset_         *os  = mesh->access->os;
    DM               da   = meshP->da, fda = meshP->fda;
    DMDALocalInfo    info = meshP->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, b, m;

    PetscMPIInt      rankP, sizeP, rank, size;
    PetscReal        dist = 0., ds;

    Cmpnts           ***cent;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    MPI_Comm_size(meshP->MESH_COMM, &sizeP);
    MPI_Comm_rank(meshP->MESH_COMM, &rankP);

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, meshP->lCent, &cent);

    std::vector<Acell> aCell = os->aCell;

    // vector that stores the number of aCell cells in each processor
    std::vector<PetscInt> NumAcellPerProc = os->NumAcellPerProc;

    std::vector<Dcell> dCellVec;
    Dcell              dCell;

    std::vector<std::vector<PetscInt>> lAcellProcMat(size);
    os->AcellProcMat.resize(size);

    os->oset_comm.resize(size);

    for(b = 0; b < size; b++){
        lAcellProcMat[b].resize(sizeP);
        os->AcellProcMat[b].resize(sizeP);
    }

    os->dCell.resize( aCell.size() );   // stores the dCell elements within a given radius of each element of aCell vector

    for(b = 0; b < aCell.size(); b++){

        dist = aCell[b].cell_size * os->cellFactor;

        for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++){

                    ds = sqrt(pow(cent[k][j][i].x - aCell[b].coorx,2.)
                            +pow(cent[k][j][i].y - aCell[b].coory,2.)
                            +pow(cent[k][j][i].z - aCell[b].coorz,2.));

                    if (ds < dist){

                        lAcellProcMat[aCell[b].rank][rankP] = 1;

                        dCell.indi = i;
                        dCell.indj = j;
                        dCell.indk = k;

                        dCell.rank = rank;
                        dCell.dist2p = ds;

                        dCellVec.push_back(dCell);

                    }

                }

        // store the connectivity for all the acceptor cells
        os->dCell[b] = dCellVec;
        std::vector<Dcell> ().swap(dCellVec);

    }

    // reduce the local processor matrix to get the global processor matrix
    // this is a matrix of dim - size x size - rows - processors of aCell cells, columns- processors of dCell cells
    for(b = 0; b < size; b++){
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMat[b][0], sizeP, MPIU_INT, MPI_SUM, meshP->MESH_COMM);
    }

    // ensure the communicator colours are set right - include the diagonal element and set all that are 0 - MPI_UNDEFINED
    // communicators will be created for each row of the processor matrix

    if(size != sizeP)
    {
      char error[512];
      sprintf(error, "current implementation requires the 2 meshes to have same number of processors. Recheck acellDcellConnectivity function\n");
      fatalErrorInFunction("acellDcellConnectivity", error);
    }

    for(b = 0; b < size; b++){
        if(NumAcellPerProc[b] != 0)
            os->AcellProcMat[b][b] = 1; // set to 1  if the aCell cells are within the given processor
        for(m = 0; m < sizeP; m++){
            if (os->AcellProcMat[b][m] == 0)
                os->AcellProcMat[b][m] = MPI_UNDEFINED;
        }
    }

    //create communicator for each row of the AcellProcMat that are non 0
    for(b = 0; b < size; b++)
    {
      if(NumAcellPerProc[b]!=0)
      {
        MPI_Comm_split(PETSC_COMM_WORLD, os->AcellProcMat[b][rankP], rankP, &(os->oset_comm[b]));
      }
    }

    std::vector<Acell> ().swap(aCell);

    // vector that stores the number of aCell cells in each processor
    std::vector<PetscInt> ().swap(NumAcellPerProc);

    std::vector<std::vector<PetscInt>> ().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, meshP->lCent, &cent);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode findClosestDonor(mesh_ *meshP, mesh_ *mesh)
{
    overset_         *os  = mesh->access->os;
    DM               da   = meshP->da, fda = meshP->fda;
    DMDALocalInfo    info = meshP->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, b, m;

    PetscMPIInt      rankP, sizeP, rank, size;
    PetscReal        ds, ***aj;

    Dcell            dCell;

    Cmpnts           ***cent;

    MPI_Comm_size(meshP->MESH_COMM, &sizeP);
    MPI_Comm_rank(meshP->MESH_COMM, &rankP);

    MPI_Comm_size(mesh->MESH_COMM, &size);
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    std::vector <Acell> aCell = os->aCell;

    std::vector <PetscInt> NumAcellPerProc = os->NumAcellPerProc;                  // vector that stores the number of aCell cells in each processor

    std::vector <std::vector<PetscInt>> lAcellProcMat(size);                         //Acell processor matrix - shows the processor connectivity between the acceptor and donor cell processors

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, meshP->lCent, &cent);
    DMDAVecGetArray(da, meshP->lAj, &aj);

    //vector that stores the closest donor cell of all acceptor cells
    os->closestDonor.resize( aCell.size() );

    os->AcellProcMat.resize(size);

    for(b = 0; b < size; b++){
        lAcellProcMat[b].resize(sizeP);
        os->AcellProcMat[b].resize(sizeP);
    }

    // loop through the query acceptor cells
    for(b = 0; b < aCell.size(); b++)
    {

        // max perturbation amplitude
        PetscReal maxPerturb   = 1e-10;
        PetscReal lClosestSize = 0.0;

        // initialize donor rank to -1
        dCell.rank = -1;

        // processor perturbation (changes between processors)
        PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)sizeP;

        // test for points accounted for twice or none
        PetscReal lminDist = 1e20;
        PetscReal gminDist = 1e20;
        std::vector<PetscInt> indices {0, 0, 0} ;

        for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {

                  ds = sqrt
                  (
                      (cent[k][j][i].x - aCell[b].coorx - procContrib) * (cent[k][j][i].x - aCell[b].coorx - procContrib) +
                      (cent[k][j][i].y - aCell[b].coory - procContrib) * (cent[k][j][i].y - aCell[b].coory - procContrib) +
                      (cent[k][j][i].z - aCell[b].coorz - procContrib) * (cent[k][j][i].z - aCell[b].coorz - procContrib)
                  );

                    if (ds < lminDist)
                    {

                        // save distance value
                        lminDist = ds + procContrib;
                        // save indices
                        indices[0] = k;
                        indices[1] = j;
                        indices[2] = i;

                    }

                }

        // scatter the local distance to global using MIN operator
        MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, meshP->MESH_COMM);

        // now compare the distances: where they agree, this processor controls the probe
        if(lminDist == gminDist)
        {
            // store the closest donor cell for all the acceptor cells
            dCell.indi   = indices[2];
            dCell.indj   = indices[1];
            dCell.indk   = indices[0];
            dCell.dist2p = gminDist;
            dCell.rank   = rankP;

            lClosestSize = pow( 1./aj[dCell.indk][dCell.indj][dCell.indi], 1./3.);
            os->closestDonor[b] = dCell;

            //set the procmatrix for this rank to 1 as the closest donor has been found
            lAcellProcMat[aCell[b].rank][rankP] = 1;
        }

        //scatter the cell size of the closest donor cell, to be used as search radius for Least square method
        MPI_Allreduce(&lClosestSize, &os->aCell[b].cell_size, 1, MPIU_REAL, MPI_SUM, meshP->MESH_COMM);

        // scatter the donor rank for acceptor cell b to all processes
        MPI_Allreduce(&dCell.rank, &os->closestDonor[b].rank, 1, MPI_INT, MPI_MAX, meshP->MESH_COMM);

    }

    // reduce the local processor matrix to get the global processor matrix
    // this is a matrix of dim - size x size - rows - processors of aCell cells, columns- processors of dCell cells
    for(b = 0; b < size; b++){
        MPI_Allreduce(&lAcellProcMat[b][0], &os->AcellProcMat[b][0], sizeP, MPIU_INT, MPI_SUM, meshP->MESH_COMM);
    }

    // ensure the communicator colours are set right - include the diagonal element and set all that are 0 - MPI_UNDEFINED
    // communicators will be created for each row of the processor matrix

    if(size != sizeP)
    {
     char error[512];
      sprintf(error, "current implementation requires the 2 meshes to have same number of processors. Recheck findClosestDonor function\n");
      fatalErrorInFunction("findClosestDonor", error);
    }

    for(b = 0; b < size; b++){
        if(NumAcellPerProc[b] != 0)
            os->AcellProcMat[b][b] = 1; // set to 1  if the aCell cells are within the given processor
        for(m = 0; m < sizeP; m++){
            if (os->AcellProcMat[b][m] == 0)
                os->AcellProcMat[b][m] = MPI_UNDEFINED;
        }
    }

    std::vector <Acell> ().swap(aCell);
    std::vector <PetscInt>   ().swap(NumAcellPerProc);
    std::vector<std::vector<PetscInt>> ().swap(lAcellProcMat);

    DMDAVecRestoreArray(fda, meshP->lCent, &cent);
    DMDAVecRestoreArray(da, meshP->lAj, &aj);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode getLSWeights(mesh_ *meshP, mesh_ *mesh){

    overset_         *os  = mesh->access->os;

    DM               da1   = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1 = mesh->info;
    DM               da0   = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0 = meshP->info;

    PetscInt         i, j, k, ii, jj, kk, b, n, m;
    PetscInt         nsupport = 0;
    PetscReal        *W, *PHI, **P, **B;
    PetscReal        A[4][4]={0.}, inv_A[4][4]={0.}, lA1[16] = {0.}, A1[16] = {0.};

    Cmpnts           ***cent1, ***cent0;
    PetscReal        ***nvert, ds, dist;
    cellIds          supportNode;

    std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
    std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

    os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

    DMDAVecGetArray(fda1, mesh->lCent, &cent1);
    DMDAVecGetArray(fda0, meshP->lCent, &cent0);

    for(b = 0; b < aCell.size(); b++){

        os->DWeights[b].resize(dCell[b].size());

        // aCell cell index
        i = aCell[b].indi;
        j = aCell[b].indj;
        k = aCell[b].indk;

        // support radius size
        ds = aCell[b].cell_size * os->cellFactor;

        nsupport = dCell[b].size();

        // allocate local variables for MLS interpolation
        B = (PetscReal**) malloc(4*sizeof(PetscReal*));

        for (n=0;n<4;n++)
        {
            B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
        }

        P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
        for (n=0;n<nsupport;n++)
        {
            P[n] = (PetscReal*) malloc(4*sizeof(PetscReal));
        }

        W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        //initialise matrix and vectors
        for (ii=0; ii<4; ii++)
        {
            for (jj=0; jj<4; jj++)
            {
                A[ii][jj] = 0.;
            }
        }
        for (ii=0; ii<16; ii++){
            A1[ii] = 0.;
            lA1[ii] = 0.;
        }

        for (n=0; n<nsupport; n++){
            PHI[n] = 0.;
            os->DWeights[b][n] = 0.;
        }

        // loop through the donor cells within the current processor - gnumdCell[b][rank]
        for(m = 0; m < nsupport; m++){

            kk = dCell[b][m].indk;
            jj = dCell[b][m].indj;
            ii = dCell[b][m].indi;

            P[m][0] = 1.0;
            P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
            P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
            P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;

            // get normalized distance
            dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
                    +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
                    +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

            // get interpolation weights
            W[m]
              =
                      (dist < 0.5)
                      ?
                              2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
            :
                              4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

            // set B
            for (ii=0; ii<4; ii++){
                B[ii][m] = W[m] * P[m][ii];
            }


            // set A
            for (ii=0; ii<4; ii++)
            {
                for (jj=0; jj<4; jj++)
                {

                    A[ii][jj] += B[ii][m] * P[m][jj];
                }
            }

        }

        // set A1 elements to scatter A matrix elements
        for (ii=0; ii<4; ii++)
        {
            for (jj=0; jj<4; jj++)
            {
                lA1[4*ii + jj] = A[ii][jj];
            }
        }

        // reduce and scatter A matrix to all processors
        MPI_Allreduce(&lA1, &A1, 16, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        for (ii=0; ii<4; ii++)
        {
            for (jj=0; jj<4; jj++)
            {
                A[ii][jj] = A1[4*ii + jj];
            }
        }

        //invert matrix A
        inv_4by4(A, inv_A);

        for(m = 0; m < nsupport; m++){
            for (ii=0; ii<4; ii++)
            {
                PHI[m] += inv_A[0][ii]*B[ii][m];
            }

            os->DWeights[b][m] = PHI[m];
        }


        // free all the vectors
        for(PetscInt ind=0; ind<4; ind++)
        {
            free(B[ind]);
        }

        for(PetscInt ind=0; ind<nsupport; ind++)
        {
            free(P[ind]);
        }

        free(B);
        free(P);
        free(W);
        free(PHI);
    }


    DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
    DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode getLSWeights_2nd(mesh_ *meshP, mesh_ *mesh)
{
    overset_         *os  = mesh->access->os;

    DM               da1   = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1 = mesh->info;
    DM               da0   = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0 = meshP->info;

    PetscInt         i, j, k, ii, jj, kk, b, n, m;
    PetscInt         nsupport = 0;
    PetscReal        *W, *PHI, **P, **B;
    PetscReal        A[10][10]={0.}, inv_A[10][10]={0.}, lA1[100] = {0.}, A1[100] = {0.};

    Cmpnts           ***cent1, ***cent0;
    PetscReal        ***nvert, ds, dist;
    cellIds          supportNode;

    std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
    std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

    os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

    DMDAVecGetArray(fda1, mesh->lCent, &cent1);
    DMDAVecGetArray(fda0, meshP->lCent, &cent0);

    for(b = 0; b < aCell.size(); b++){

        os->DWeights[b].resize(dCell[b].size());

        // aCell cell index
        i = aCell[b].indi;
        j = aCell[b].indj;
        k = aCell[b].indk;

        // support radius size
        ds = aCell[b].cell_size * os->cellFactor;

        nsupport = dCell[b].size();

        // allocate local variables for MLS interpolation
        B = (PetscReal**) malloc(10*sizeof(PetscReal*));
        for (n=0;n<10;n++)
        {
            B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
        }

        P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
        for (n=0;n<nsupport;n++)
        {
            P[n] = (PetscReal*) malloc(10*sizeof(PetscReal));
        }

        W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        //initialise matrix and vectors
        for (ii=0; ii<10; ii++)
        {
            for (jj=0; jj<10; jj++)
            {
                A[ii][jj] = 0.;
            }
        }
        for (ii=0; ii<100; ii++){
            A1[ii] = 0.;
            lA1[ii] = 0.;
        }

        for (n=0; n<nsupport; n++){
            PHI[n] = 0.;
            os->DWeights[b][n] = 0.;
        }

        // loop through the donor cells within the current processor - gnumdCell[b][rank]
        for(m = 0; m < nsupport; m++){

            kk = dCell[b][m].indk;
            jj = dCell[b][m].indj;
            ii = dCell[b][m].indi;

            P[m][0] = 1.0;
            P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
            P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
            P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;
            P[m][4]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
            P[m][5]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
            P[m][6]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
            P[m][7]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
            P[m][8]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
            P[m][9]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;

            // get normalized distance
            dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
                    +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
                    +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

            // get interpolation weights
            W[m]
              =
                      (dist < 0.5)
                      ?
                              2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
            :
                              4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

            // set B
            for (ii=0; ii<10; ii++){
                B[ii][m] = W[m] * P[m][ii];
            }


            // set A
            for (ii=0; ii<10; ii++)
            {
                for (jj=0; jj<10; jj++)
                {

                    A[ii][jj] += B[ii][m] * P[m][jj];
                }
            }

        }

        // set A1 elements to scatter A matrix elements
        for (ii=0; ii<10; ii++)
        {
            for (jj=0; jj<10; jj++)
            {
                lA1[10*ii + jj] = A[ii][jj];
            }
        }

        // reduce and scatter A matrix to all processors
        MPI_Allreduce(&lA1, &A1, 100, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

        for (ii=0; ii<10; ii++)
        {
            for (jj=0; jj<10; jj++)
            {
                A[ii][jj] = A1[10*ii + jj];
            }
        }

        //invert matrix A
        inv(A,inv_A,10);

        for(m = 0; m < nsupport; m++){
            for (ii=0; ii<10; ii++)
            {
                PHI[m] += inv_A[0][ii]*B[ii][m];
            }

            os->DWeights[b][m] = PHI[m];
        }


        // free all the vectors
        for(PetscInt ind=0; ind<10; ind++)
        {
            free(B[ind]);
        }

        for(PetscInt ind=0; ind<nsupport; ind++)
        {
            free(P[ind]);
        }

        free(B);
        free(P);
        free(W);
        free(PHI);
    }


    DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
    DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
    return 0;
}

//***************************************************************************************************************//

// third order least square method for finding the weights for the donor cells of the background mesh to interpolate the overset mesh acceptor cells
PetscErrorCode getLSWeights_3rd(mesh_ *meshP, mesh_ *mesh)
{
    overset_         *os  = mesh->access->os;

    DM               da1   = mesh->da, fda1 = mesh->fda;
    DMDALocalInfo    info1 = mesh->info;
    DM               da0   = meshP->da, fda0 = meshP->fda;
    DMDALocalInfo    info0 = meshP->info;

    PetscInt         ii, jj, kk, b, n, m;
    PetscInt         nsupport = 0;
    PetscReal        *W, *PHI, **P, **B;
    PetscReal        A[20][20]={0.}, inv_A[20][20]={0.}, lA1[400] = {0.}, A1[400] = {0.};

    Cmpnts           ***cent1, ***cent0;
    PetscReal        ***nvert, ds, dist;
    cellIds          supportNode;

    std::vector<Acell> aCell = os->aCell;   // acceptor cells of overset mesh
    std::vector<std::vector<Dcell>> dCell = os->dCell;  // donor cells of background mesh

    os->DWeights.resize(aCell.size());      // vector to store the LS weights for each acceptor cell

    DMDAVecGetArray(fda1, mesh->lCent, &cent1);
    DMDAVecGetArray(fda0, meshP->lCent, &cent0);

    for(b = 0; b < aCell.size(); b++){

        os->DWeights[b].resize(dCell[b].size());

        // support radius size
        ds = aCell[b].cell_size * os->cellFactor;

        nsupport = dCell[b].size();

        // allocate local variables for MLS interpolation
        B = (PetscReal**) malloc(20*sizeof(PetscReal*));
        for (n=0;n<20;n++)
        {
            B[n] = (PetscReal*) malloc(nsupport*sizeof(PetscReal));
        }

        P = (PetscReal**) malloc(nsupport*sizeof(PetscReal*));
        for (n=0;n<nsupport;n++)
        {
            P[n] = (PetscReal*) malloc(20*sizeof(PetscReal));
        }

        W =  (PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        PHI =(PetscReal* ) malloc(nsupport*sizeof(PetscReal));

        //initialise matrix and vectors
        for (ii=0; ii<20; ii++)
        {
            for (jj=0; jj<20; jj++)
            {
                A[ii][jj] = 0.;
            }
        }
        for (ii=0; ii<400; ii++){
            A1[ii] = 0.;
            lA1[ii] = 0.;
        }

        for (n=0; n<nsupport; n++){
            PHI[n] = 0.;
            os->DWeights[b][n] = 0.;
        }

        // loop through the donor cells within the current processor - gnumdCell[b][rank]
        for(m = 0; m < nsupport; m++){

            kk = dCell[b][m].indk;
            jj = dCell[b][m].indj;
            ii = dCell[b][m].indi;

            P[m][0] = 1.0;
            P[m][1] = (cent0[kk][jj][ii].x - aCell[b].coorx) / ds;
            P[m][2] = (cent0[kk][jj][ii].y - aCell[b].coory) / ds;
            P[m][3] = (cent0[kk][jj][ii].z - aCell[b].coorz) / ds;
            P[m][4]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
            P[m][5]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
            P[m][6]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
            P[m][7]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds;
            P[m][8]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds;
            P[m][9]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds;
            P[m][10]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
            P[m][11]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
            P[m][12]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
            P[m][13]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
            P[m][14]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
            P[m][15]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
            P[m][16]=(cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;
            P[m][17]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].x - aCell[b].coorx)/ds/ds/ds;
            P[m][18]=(cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].z - aCell[b].coorz) * (cent0[kk][jj][ii].y - aCell[b].coory)/ds/ds/ds;
            P[m][19]=(cent0[kk][jj][ii].x - aCell[b].coorx) * (cent0[kk][jj][ii].y - aCell[b].coory) * (cent0[kk][jj][ii].z - aCell[b].coorz)/ds/ds/ds;


            // get normalized distance
            dist=sqrt(pow(cent0[kk][jj][ii].x - aCell[b].coorx,2.)
                    +pow(cent0[kk][jj][ii].y - aCell[b].coory,2.)
                    +pow(cent0[kk][jj][ii].z - aCell[b].coorz,2.)) / ds;

            // get interpolation weights
            W[m]
              =
                      (dist < 0.5)
                      ?
                              2.0 / 3.0 - 4.0 * dist * dist + 4.0 * pow(dist, 3.0)
            :
                              4.0 / 3.0 - 4.0 * dist + 4.0 * pow(dist, 2.0) - 4.0 / 3.0 * pow(dist, 3.0);

            // set B
            for (ii=0; ii<20; ii++){
                B[ii][m] = W[m] * P[m][ii];
            }


            // set A
            for (ii=0; ii<20; ii++)
            {
                for (jj=0; jj<20; jj++)
                {

                    A[ii][jj] += B[ii][m] * P[m][jj];
                }
            }

        }

        // set A1 elements to scatter A matrix elements
        for (ii=0; ii<20; ii++)
        {
            for (jj=0; jj<20; jj++)
            {
                lA1[20*ii + jj] = A[ii][jj];
            }
        }

        // reduce and scatter A matrix to all processors
        MPI_Allreduce(&lA1, &A1, 400, MPIU_REAL, MPIU_SUM, PETSC_COMM_WORLD);

        for (ii=0; ii<20; ii++)
        {
            for (jj=0; jj<20; jj++)
            {
                A[ii][jj] = A1[20*ii + jj];
            }
        }

        //invert matrix A
        inv_20(A,inv_A,20);

        for(m = 0; m < nsupport; m++){
            for (ii=0; ii<20; ii++)
            {
                PHI[m] += inv_A[0][ii]*B[ii][m];
            }

            os->DWeights[b][m] = PHI[m];
        }


        // free all the vectors
        for(PetscInt ind=0; ind<20; ind++)
        {
            free(B[ind]);
        }

        for(PetscInt ind=0; ind<nsupport; ind++)
        {
            free(P[ind]);
        }

        free(B);
        free(P);
        free(W);
        free(PHI);
    }


    DMDAVecRestoreArray(fda1, mesh->lCent, &cent1);
    DMDAVecRestoreArray(fda0, meshP->lCent, &cent0);
    return 0;
}

//***************************************************************************************************************//

PetscErrorCode oversetContravariantBC(mesh_ *mesh, PetscInt i, PetscInt j, PetscInt k, Cmpnts ucart, PetscInt face)
{

    ueqn_         *ueqn = mesh->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts           ***ucont, ***icsi, ***jeta, ***kzet;

    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    // flux bc based on the face where ucart is interpolated
    if(face == 1 && !(mesh->i_periodic) && !(mesh->ii_periodic))
    {
        ucont[k][j][i-1].x = (ucart.x * icsi[k][j][i-1].x + ucart.y * icsi[k][j][i-1].y + ucart.z * icsi[k][j][i-1].z);
    }

    if(face == 2 && !(mesh->i_periodic) && !(mesh->ii_periodic))
    {
        ucont[k][j][i].x = (ucart.x * icsi[k][j][i].x + ucart.y * icsi[k][j][i].y + ucart.z * icsi[k][j][i].z);
    }

    if(face == 3 && !(mesh->j_periodic) && !(mesh->jj_periodic))
    {
        ucont[k][j-1][i].y = (ucart.x * jeta[k][j-1][i].x + ucart.y * jeta[k][j-1][i].y + ucart.z * jeta[k][j-1][i].z);
    }

    if(face == 4 && !(mesh->j_periodic) && !(mesh->jj_periodic))
    {
        ucont[k][j][i].y = (ucart.x * jeta[k][j][i].x + ucart.y * jeta[k][j][i].y + ucart.z * jeta[k][j][i].z);
    }

    if(face == 5 && !(mesh->k_periodic) && !(mesh->kk_periodic))
    {
        ucont[k-1][j][i].z = (ucart.x * kzet[k-1][j][i].x + ucart.y * kzet[k-1][j][i].y + ucart.z * kzet[k-1][j][i].z );
    }

    if(face == 6 && !(mesh->k_periodic) && !(mesh->kk_periodic))
    {
        ucont[k][j][i].z = (ucart.x * kzet[k][j][i].x + ucart.y * kzet[k][j][i].y + ucart.z * kzet[k][j][i].z );
    }

    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

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

void defineStruct_Acell(MPI_Datatype *tstype) {
    const PetscInt    count = 9;
    int               blocklens[count];
    MPI_Aint          disps[count];

    for (PetscInt i=0; i < count; i++) {
        blocklens[i] = 1;
    }

    MPI_Datatype types[count] = {MPIU_INT, MPIU_INT, MPIU_INT, MPIU_REAL, MPIU_REAL, MPIU_REAL, MPIU_INT, MPIU_REAL, MPIU_INT};

    disps[0] = offsetof(Acell,indi);
    disps[1] = offsetof(Acell,indj);
    disps[2] = offsetof(Acell,indk);
    disps[3] = offsetof(Acell,coorx);
    disps[4] = offsetof(Acell,coory);
    disps[5] = offsetof(Acell,coorz);
    disps[6] = offsetof(Acell,rank);
    disps[7] = offsetof(Acell,cell_size);
    disps[8] = offsetof(Acell,face);

    MPI_Type_create_struct(count, blocklens, disps, types, tstype);
    MPI_Type_commit(tstype);
    return;
}

//***************************************************************************************************************//

//MPI operation function to find the sum the elements of the vector of structs
void sum_struct_Acell(void *in, void *inout, int *len, MPI_Datatype *type){
    /* ignore type, just trust that it's our struct type */

    Acell *invals    = (Acell*)in;
    Acell *inoutvals = (Acell*)inout;

    for (PetscInt i=0; i<*len; i++) {
        inoutvals[i].indi  += invals[i].indi;
        inoutvals[i].indj  += invals[i].indj;
        inoutvals[i].indk  += invals[i].indk;

        inoutvals[i].coorx  += invals[i].coorx;
        inoutvals[i].coory  += invals[i].coory;
        inoutvals[i].coorz  += invals[i].coorz;

        inoutvals[i].rank  += invals[i].rank;
        inoutvals[i].cell_size  += invals[i].cell_size;
        inoutvals[i].face  += invals[i].face;

    }

    return;
}
