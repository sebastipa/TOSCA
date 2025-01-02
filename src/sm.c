//! \brief Contains SM equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/boundary.h"
#include <Eigen/Dense>
#include <iostream>


//***************************************************************************************************************//

PetscErrorCode SNESMonitorSM(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeSMObject(SMObj_ *smObject)
{
    if(smObject != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = smObject->sm[0]->access->mesh;
        flags_ *flags = smObject->sm[0]->access->flags;

        for (PetscInt eig = 0; eig < 3; ++eig)
        {
            VecDuplicate(mesh->Nvert,   &(smObject->weightAbsc[eig]->weight));      VecSet(smObject->weightAbsc[eig]->weight,   0.0);
            VecDuplicate(mesh->Nvert,   &(smObject->weightAbsc[eig]->absc));      VecSet(smObject->weightAbsc[eig]->absc,   0.0);
            VecDuplicate(mesh->Nvert,   &(smObject->weightAbsc[eig]->tauP));      VecSet(smObject->weightAbsc[eig]->tauP,   0.0);

        }

        VecDuplicate(mesh->Nvert,   &(smObject->Eps));   VecSet(smObject->Eps,0.);
        VecDuplicate(mesh->Nvert,   &(smObject->quant));   VecSet(smObject->quant,0.);
        VecDuplicate(mesh->Nvert,   &(smObject->Dq));   VecSet(smObject->Dq,0.);
        VecDuplicate(mesh->Nvert,   &(smObject->probI));   VecSet(smObject->probI,0.);

        VecDuplicate(mesh->Nvert, &(smObject->ExCount));      VecSet(smObject->ExCount,   0.0);
        VecDuplicate(mesh->Nvert, &(smObject->DepCount));      VecSet(smObject->DepCount,   0.0);

        smObject->dissWeight = 0;
        smObject->OGConc = 0;
        smObject->monoDFrac = 0;
        smObject->rhoPart = 0;

    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeSM(sm_ *sm)
{
    if(sm != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = sm->access->mesh;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        VecDuplicate(mesh->Nvert, &(sm->smTmp));    VecSet(sm->smTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(sm->smVal));       VecSet(sm->smVal,    0.0);
        VecDuplicate(mesh->Nvert, &(sm->sm_o));     VecSet(sm->sm_o,  0.0);
        VecDuplicate(mesh->Nvert, &(sm->Rhs));         VecSet(sm->Rhs,      0.0);
        VecDuplicate(mesh->Nvert, &(sm->Rhs_o));       VecSet(sm->Rhs_o,    0.0);
        VecDuplicate(mesh->lAj,   &(sm->lsmVal));      VecSet(sm->lsmVal,   0.0);
        VecDuplicate(mesh->lAj,   &(sm->lsm_o));    VecSet(sm->lsm_o, 0.0);

        VecDuplicate(mesh->lCent, &(sm->lDivSM));       VecSet(sm->lDivSM,    0.0);
        VecDuplicate(mesh->lCent, &(sm->lViscSM));      VecSet(sm->lViscSM,   0.0);

        VecDuplicate(mesh->Cent, &(sm->Sed));      VecSet(sm->Sed,   0.0);
        VecDuplicate(mesh->lCent, &(sm->lSed));      VecSet(sm->lSed,   0.0);

        VecDuplicate(mesh->Nvert, &(sm->sedCent));       VecSet(sm->sedCent,    0.0);
        VecDuplicate(mesh->lAj,   &(sm->lSedCent));      VecSet(sm->lSedCent,   0.0);

        VecDuplicate(mesh->Cent, &(sm->Dev));      VecSet(sm->Dev,   0.0);
        VecDuplicate(mesh->lCent, &(sm->lDev));      VecSet(sm->lDev,   0.0);

        VecDuplicate(mesh->Nvert, &(sm->devCent));       VecSet(sm->devCent,    0.0);
        VecDuplicate(mesh->lAj,   &(sm->lDevCent));      VecSet(sm->lDevCent,   0.0);

        VecDuplicate(mesh->Cent, &(sm->Dep));      VecSet(sm->Dep,   0.0);
        VecDuplicate(mesh->lCent, &(sm->lDep));      VecSet(sm->lDep,   0.0);

        VecDuplicate(mesh->Nvert, &(sm->coagSource));      VecSet(sm->coagSource,   0.0);

        // read time discretization scheme
        readDictWord("control.dat", "-dSMdtScheme", &(sm->ddtScheme));

        if (sm->ddtScheme=="rungeKutta4")
        {

        }
        else
        {
            char error[512];
            sprintf(error, "unknown ddtScheme %s for SM equation, available schemes are\n    1. rungeKutta4\n", sm->ddtScheme.c_str());
            fatalErrorInFunction("InitializeSM", error);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii)
{
    // In this function the viscous + divergence term of the scalarMoment equation are
    // discretized at cell centers.
    // First the divergence and viscous fluxes are evaluated at cell faces, then
    // their budget is evaluated at the internal cells, forming the Rhs.

    mesh_         *mesh  = sm->access->mesh;
    ueqn_         *ueqn  = sm->access->ueqn;
    les_          *les   = sm->access->les;
    constants_    *cst   = sm->access->constants;

    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***ucont;

    Cmpnts        ***csi, ***eta, ***zet, ***cent;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;

    PetscReal     ***lsmVal, ***rhs, ***nvert;

    Cmpnts        ***div, ***visc;                                                // divergence and viscous terms
    Cmpnts        ***limiter;                                                     // flux limiter
    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;                                  // cell and face jacobians
    PetscReal     ***lnu_t, ***Sabs, ***lch, ***lcs;

    PetscReal     dsdc, dsde, dsdz;                                              // velocity der. w.r.t. curvil. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components


    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  Rhs,               &rhs);
    DMDAVecGetArray(da,  sm->lsmVal,      &lsmVal);
    DMDAVecGetArray(fda, ueqn->lUcont,      &ucont);

    DMDAVecGetArray(fda, mesh->lCsi,        &csi);
    DMDAVecGetArray(fda, mesh->lEta,        &eta);
    DMDAVecGetArray(fda, mesh->lZet,        &zet);
    DMDAVecGetArray(fda, mesh->lICsi,       &icsi);
    DMDAVecGetArray(fda, mesh->lIEta,       &ieta);
    DMDAVecGetArray(fda, mesh->lIZet,       &izet);
    DMDAVecGetArray(fda, mesh->lJCsi,       &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta,       &jeta);
    DMDAVecGetArray(fda, mesh->lJZet,       &jzet);
    DMDAVecGetArray(fda, mesh->lKCsi,       &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta,       &keta);
    DMDAVecGetArray(fda, mesh->lKZet,       &kzet);
    DMDAVecGetArray(da,  mesh->lAj,         &aj);
    DMDAVecGetArray(da,  mesh->lIAj,        &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,        &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,        &kaj);
    DMDAVecGetArray(fda, mesh->lCent,       &cent);
    DMDAVecGetArray(da,  mesh->lNvert,      &nvert);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    VecSet(sm->lDivSM,  0.0);
    VecSet(sm->lViscSM, 0.0);

    DMDAVecGetArray(fda, sm->lDivSM, &div);
    DMDAVecGetArray(fda, sm->lViscSM, &visc);

    if (sm->access->flags->isLesActive)
    {
        if(sm->access->flags->isLesActive != 2)
        {
            DMDAVecGetArray(da, les->lNu_t, &lnu_t);
        }
        else
        {
            DMDAVecGetArray(da,  les->lS,  &Sabs);
            DMDAVecGetArray(da,  les->lCh, &lch);
            DMDAVecGetArray(da,  les->lCs, &lcs);
        }
    }

    // ---------------------------------------------------------------------- //
    //                FORM DIVERGENCE AND VISCOUS CONTRIBUTIONS               //
    // ---------------------------------------------------------------------- //

    // i direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                //skip non-internal faces
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(j==0 || k==0) continue;

                //skip all solid cells and fluid cells with solid on right face
                if (isIBMSolidCell(k, j, i, nvert) || isIBMSolidCell(k, j, i+1, nvert)) continue;

                //divergence term using central4 scheme
                if (i == 0 || isIBMSolidCell(k, j, i-1, nvert))
                {
                    div[k][j][i].x =
                    - ucont[k][j][i].x
                    * central4
                    (
                        lsmVal[k][j][i],
                        lsmVal[k][j][i],
                        lsmVal[k][j][i+1],
                        lsmVal[k][j][i+2]
                    );
                }
                else if (i == mx-2 || isIBMSolidCell(k, j, i+2, nvert))
                {
                    div[k][j][i].x =
                    - ucont[k][j][i].x
                    * central4
                    (
                        lsmVal[k][j][i-1],
                        lsmVal[k][j][i],
                        lsmVal[k][j][i+1],
                        lsmVal[k][j][i+1]
                    );
                }
                else
                {
                    div[k][j][i].x =
                    - ucont[k][j][i].x
                    * central4
                    (
                        lsmVal[k][j][i-1],
                        lsmVal[k][j][i],
                        lsmVal[k][j][i+1],
                        lsmVal[k][j][i+2]
                    );
                }


                //move on to visc term
                // get 1/V at the i-face
                PetscReal ajc = iaj[k][j][i];

                // get face normals
                csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;

                // compute metric tensor - WARNING: there is a factor of 1/J^2 if using face area vectors
                //                                  must multiply for ajc in viscous term!!!
                g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                Compute_dscalar_i(mesh, i, j, k, mx, my, mz, lsmVal, nvert, &dsdc, &dsde, &dsdz);

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;

                // viscous terms
                if
                (
                    sm->access->flags->isLesActive
                )
                {
                    if(sm->access->flags->isLesActive != 2)
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            nut = lnu_t[k][j][i+1];
                        }
                        else if (isIBMCell(k, j, i+1, nvert))
                        {
                            nut = lnu_t[k][j][i];
                        }
                        else
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);
                        }

                        gammaEff = (nu / cst->Sc) + (nut / cst->ScT);
                    }
                    else
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k][j][i+1] / 0.1 * Sabs[k][j][i+1] * lcs[k][j][i+1] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else if (isIBMCell(k, j, i+1, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k][j][i] / 0.1 * Sabs[k][j][i] * lcs[k][j][i] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else
                        {
                            gammaEff = (nu / cst->Sc) + 0.5 *
                            (
                                lch[k][j][i+1] / 0.1 * Sabs[k][j][i+1] * lcs[k][j][i+1] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }
                    }

                }
                else
                {
                    gammaEff = nu / cst->Sc;
                }

                visc[k][j][i].x += (g11 * dsdc + g21 * dsde + g31 * dsdz) * ajc * (gammaEff);

            }
        }
    }

    // j direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                //skip non-internal faces
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || k==0) continue;

                //skip all solid cells and fluid cells with solid on right face
                if (isIBMSolidCell(k, j, i, nvert) || isIBMSolidCell(k, j+1, i, nvert)) continue;

                //divergence term using central4 scheme
                if (j == 0 || isIBMSolidCell(k, j-1, i, nvert))
                {
                    div[k][j][i].y =
                    - ucont[k][j][i].y
                    * central4
                    (
                        lsmVal[k][j][i],
                        lsmVal[k][j][i],
                        lsmVal[k][j+1][i],
                        lsmVal[k][j+2][i]
                    );
                }
                else if (j == my-2 || isIBMSolidCell(k, j+2, i, nvert))
                {
                    div[k][j][i].y =
                    - ucont[k][j][i].y
                    * central4
                    (
                        lsmVal[k][j-1][i],
                        lsmVal[k][j][i],
                        lsmVal[k][j+1][i],
                        lsmVal[k][j+1][i]
                    );
                }
                else
                {
                    div[k][j][i].y =
                    - ucont[k][j][i].y
                    * central4
                    (
                        lsmVal[k][j-1][i],
                        lsmVal[k][j][i],
                        lsmVal[k][j+1][i],
                        lsmVal[k][j+2][i]
                    );
                }

                // get 1/V at the j-face
                PetscReal ajc = jaj[k][j][i];

                // get face normals
                csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;

                // compute metric tensor
                g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
                g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
                g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

                Compute_dscalar_j(mesh, i, j, k, mx, my, mz, lsmVal, nvert, &dsdc, &dsde, &dsdz);

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;

                if
                (
                    sm->access->flags->isLesActive
                )
                {
                    if(sm->access->flags->isLesActive != 2)
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            nut = lnu_t[k][j+1][i];
                        }
                        else if (isIBMCell(k, j+1, i, nvert))
                        {
                            nut = lnu_t[k][j][i];
                        }
                        else
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);
                        }

                        gammaEff = (nu / cst->Sc) + (nut / cst->ScT);
                    }
                    else
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k][j+1][i] / 0.1 * Sabs[k][j+1][i] * lcs[k][j+1][i] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else if (isIBMCell(k, j+1, i, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k][j][i] / 0.1 * Sabs[k][j][i] * lcs[k][j][i] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else
                        {
                            gammaEff = (nu / cst->Sc) + 0.5 *
                            (
                                lch[k][j+1][i] / 0.1 * Sabs[k][j+1][i] * lcs[k][j+1][i] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }
                    }

                }
                else
                {
                    gammaEff = nu / cst->Sc;
                }

                visc[k][j][i].y += (g11 * dsdc + g21 * dsde + g31 * dsdz) * ajc * (gammaEff);

            }
        }
    }

    // k direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                //skip non-internal faces
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || j==0) continue;

                //skip all solid cells and fluid cells with solid on right face
                if (isIBMSolidCell(k, j, i, nvert) || isIBMSolidCell(k+1, j, i, nvert)) continue;

                //divergence term using central4 scheme
                if (k == 0 || isIBMSolidCell(k-1, j, i, nvert))
                {
                    div[k][j][i].z =
                    - ucont[k][j][i].z
                    * central4
                    (
                        lsmVal[k][j][i],
                        lsmVal[k][j][i],
                        lsmVal[k+1][j][i],
                        lsmVal[k+2][j][i]
                    );
                }
                else if (k == mz-2 || isIBMSolidCell(k+2, j, i, nvert))
                {
                    div[k][j][i].z =
                    - ucont[k][j][i].z
                    * central4
                    (
                        lsmVal[k-1][j][i],
                        lsmVal[k][j][i],
                        lsmVal[k+1][j][i],
                        lsmVal[k+1][j][i]
                    );
                }
                else
                {
                    div[k][j][i].z =
                    - ucont[k][j][i].z
                    * central4
                    (
                        lsmVal[k-1][j][i],
                        lsmVal[k][j][i],
                        lsmVal[k+1][j][i],
                        lsmVal[k+2][j][i]
                    );
                }

                // get 1/V at the j-face
                PetscReal ajc = kaj[k][j][i];

                // get face normals
                csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                // compute metric tensor
                g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                Compute_dscalar_k(mesh, i, j, k, mx, my, mz, lsmVal, nvert, &dsdc, &dsde, &dsdz);

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;

                if
                (
                    sm->access->flags->isLesActive
                )
                {
                    if(sm->access->flags->isLesActive != 2)
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            nut = lnu_t[k+1][j][i];
                        }
                        else if (isIBMCell(k+1, j, i, nvert))
                        {
                            nut = lnu_t[k][j][i];
                        }
                        else
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);
                        }

                        gammaEff = (nu / cst->Sc) + (nut / cst->ScT);
                    }
                    else
                    {
                        if(isIBMCell(k, j, i, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k+1][j][i] / 0.1 * Sabs[k+1][j][i] * lcs[k+1][j][i] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else if (isIBMCell(k+1, j, i, nvert))
                        {
                            gammaEff = (nu / cst->Sc) + lch[k][j][i] / 0.1 * Sabs[k][j][i] * lcs[k][j][i] * pow(1.0 / ajc, 1.0/3.0);
                        }
                        else
                        {
                            gammaEff = (nu / cst->Sc) + 0.5 *
                            (
                                lch[k+1][j][i] / 0.1 * Sabs[k+1][j][i] * lcs[k+1][j][i] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }
                    }
                }
                else
                {
                    gammaEff = nu / cst->Sc;
                }

                visc[k][j][i].z += (g11 * dsdc + g21 * dsde + g31 * dsdz) * ajc * (gammaEff);

            }
        }
    }

    DMDAVecRestoreArray(fda, sm->lDivSM, &div);
    DMDAVecRestoreArray(fda, sm->lViscSM, &visc);

    DMLocalToLocalBegin(fda, sm->lDivSM,  INSERT_VALUES, sm->lDivSM);
    DMLocalToLocalEnd  (fda, sm->lDivSM,  INSERT_VALUES, sm->lDivSM);
    DMLocalToLocalBegin(fda, sm->lViscSM, INSERT_VALUES, sm->lViscSM);
    DMLocalToLocalEnd  (fda, sm->lViscSM, INSERT_VALUES, sm->lViscSM);

    resetFacePeriodicFluxesVector(mesh, sm->lDivSM,   sm->lDivSM, "localToLocal");
    resetFacePeriodicFluxesVector(mesh, sm->lViscSM, sm->lViscSM, "localToLocal");

    DMDAVecGetArray(fda, sm->lDivSM, &div);
    DMDAVecGetArray(fda, sm->lViscSM, &visc);

    // ---------------------------------------------------------------------- //
    //                      FORM THE CUMULATIVE FLUXES                        //
    // ---------------------------------------------------------------------- //

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                //don't skip ibfluid in this caseto allow for SM build up
                if (isIBMSolidCell(k, j, i, nvert))
                {
                    rhs[k][j][i] = 0;
                }
                else
                {
                    rhs[k][j][i]
                    =
                    scale *
                    (
                         div[k][j][i].x  - div[k  ][j  ][i-1].x       +
                         div[k][j][i].y  - div[k  ][j-1][i  ].y       +
                         div[k][j][i].z  - div[k-1][j  ][i  ].z       +
                         visc[k][j][i].x - visc[k  ][j  ][i-1].x      +
                         visc[k][j][i].y - visc[k  ][j-1][i  ].y      +
                         visc[k][j][i].z - visc[k-1][j  ][i  ].z
                    ) * aj[k][j][i];

                }

            }
        }
    }

    if(sm->access->flags->isLesActive)
    {
        if(sm->access->flags->isLesActive != 2)
        {
            DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
        }
        else
        {
            DMDAVecRestoreArray(da,  les->lS,  &Sabs);
            DMDAVecRestoreArray(da,  les->lCh, &lch);
            DMDAVecRestoreArray(da,  les->lCs, &lcs);
        }
    }


    DMDAVecRestoreArray(fda, sm->lDivSM,       &div);
    DMDAVecRestoreArray(fda, sm->lViscSM,      &visc);

    DMDAVecRestoreArray(da,  Rhs,              &rhs);
    DMDAVecRestoreArray(da,  sm->lsmVal,      &lsmVal);
    DMDAVecRestoreArray(fda, ueqn->lUcont,     &ucont);

    DMDAVecRestoreArray(fda, mesh->lCsi,        &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,        &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,        &zet);
    DMDAVecRestoreArray(fda, mesh->lICsi,       &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta,       &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet,       &izet);
    DMDAVecRestoreArray(fda, mesh->lJCsi,       &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,       &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet,       &jzet);
    DMDAVecRestoreArray(fda, mesh->lKCsi,       &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta,       &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet,       &kzet);
    DMDAVecRestoreArray(da,  mesh->lAj,         &aj);
    DMDAVecRestoreArray(da,  mesh->lIAj,        &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,        &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,        &kaj);
    DMDAVecRestoreArray(fda, mesh->lCent,       &cent);
    DMDAVecRestoreArray(da,  mesh->lNvert,      &nvert);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    return(0);
}


//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsSM(sm_ *sm, PetscInt ii)
{
    mesh_ *mesh   = sm->access->mesh;
    PetscReal ts1,te1;

    // reset scalarMoment periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, sm->smVal, sm->lsmVal, "scalar", "globalToLocal");

    // initialize the rhs vector
    VecSet(sm->Rhs, 0.0);

    // add viscous and transport terms
    //PetscTime(&ts1);
    FormSM(sm, sm->Rhs, 1.0, ii);
    //PetscTime(&te1);
    //PetscPrintf(mesh->MESH_COMM,"%li Elapsed Time VISC = %f\n", ii, te1-ts1);

    // add coag source terms
    if(sm->access->flags->isCoagSourceActive)
    {
        //PetscTime(&ts1);
        formCoagSourceExp(sm, ii);
        sourceSMCoag(sm, sm->Rhs, -1.0);
        //PetscTime(&te1);
        //PetscPrintf(mesh->MESH_COMM,"%li Elapsed Time COAG = %f\n", ii, te1-ts1);

    }

    // add deposition source terms
    if(sm->access->flags->isDepoSourceActive)
    {
        //PetscTime(&ts);
        //PetscTime(&ts1);
        formDepSourceExp(sm, ii);
        sourceSMDep(sm, sm->Rhs, 1.0);
        //PetscTime(&te1);
        //PetscPrintf(mesh->MESH_COMM,"%li Elapsed Time DEP = %f\n", ii, te1-ts1);

    }

    // add deposition source terms
    if(sm->access->flags->isSediFluxActive)
    {
        //PetscTime(&ts1);
        VecSet(sm->sedCent, 0.0);
        sedFluxSM(sm, sm->Rhs, -1.0, ii);
        //PetscTime(&te1);
        //PetscPrintf(mesh->MESH_COMM,"%li Elapsed Time SED = %f\n", ii, te1-ts1);

    }

    // add deposition source terms
    if(sm->access->flags->isDeviFluxActive)
    {
        //PetscTime(&ts1);
        VecSet(sm->devCent, 0.0);
        devFluxSM(sm, sm->Rhs, -1.0, ii);
        //PetscTime(&te1);
        //PetscPrintf(mesh->MESH_COMM,"%li Elapsed Time DEVV = %f\n", ii, te1-ts1);

    }

    // set to zero at non-resolved cell faces, this deos not includeIBFluid for SM values!
    resetNonResolvedCellCentersScalarIBMSolid(mesh, sm->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SMRK4(SMObj_ *smObject)
{
    mesh_  *mesh  = smObject->sm[0]->access->mesh;
    clock_ *clock = smObject->sm[0]->access->clock;
    flags_ *flags = smObject->sm[0]->access->flags;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for SM, Stage ");

    PetscInt  s = 4;
    PetscReal b[4];
    PetscReal a[4];

    b[0] = 1.0 / 6.0;
    b[1] = 1.0 / 3.0;
    b[2] = 1.0 / 3.0;
    b[3] = 1.0 / 6.0;

    a[0] = 0.0;
    a[1] = 0.5;
    a[2] = 0.5;
    a[3] = 1.0;

    PetscReal dt = clock->dt;

    for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
    {
        // Ucont_o contribution
        VecCopy(smObject->sm[ii]->sm_o, smObject->sm[ii]->smTmp);
    }


    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate U guess and evaluate RHS
        if(i!=0)
        {
            for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
            {
                VecWAXPY(smObject->sm[ii]->smVal, a[i] * dt, smObject->sm[ii]->Rhs, smObject->sm[ii]->sm_o);
                DMGlobalToLocalBegin(mesh->da, smObject->sm[ii]->smVal, INSERT_VALUES, smObject->sm[ii]->lsmVal);
                DMGlobalToLocalEnd(mesh->da, smObject->sm[ii]->smVal, INSERT_VALUES, smObject->sm[ii]->lsmVal);

                UpdateScalarMomentBCs(smObject->sm[ii], ii);
            }

        }

        if (flags->isCoagSourceActive || flags->isDepoSourceActive || flags->isSediFluxActive || flags->isDeviFluxActive)
        {
            quickUpdateWeightsAndAbscissi(smObject);
        }


        for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
        {
            // compute function guess for div and visc terms. Also finds sm's at faces
            FormExplicitRhsSM(smObject->sm[ii], ii);
        }

        for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
        {
            // add contribution from K1, K2, K3, K4
            VecAXPY(smObject->sm[ii]->smTmp, dt * b[i], smObject->sm[ii]->Rhs);
        }


    }

    for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
    {
        VecCopy(smObject->sm[ii]->smTmp, smObject->sm[ii]->smVal);
        DMGlobalToLocalBegin(mesh->da, smObject->sm[ii]->smVal, INSERT_VALUES, smObject->sm[ii]->lsmVal);
        DMGlobalToLocalEnd(mesh->da, smObject->sm[ii]->smVal, INSERT_VALUES, smObject->sm[ii]->lsmVal);
    }

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveSM(SMObj_ *smObject)
{
    mesh_          *mesh = smObject->sm[0]->access->mesh;
    flags_         *flags = smObject->sm[0]->access->flags;

    // set the right hand side to zero
    for (PetscInt ii = 0; ii < flags->isScalarMomentsActive; ii++)
    {
       VecSet (smObject->sm[ii]->Rhs, 0.0);
    }

    if (smObject->sm[0]->ddtScheme=="rungeKutta4")
    {
        SMRK4(smObject);
    }
    else
    {
        char error[512];
        sprintf(error, "SM solver is only compatible with rungeKutta4");
        fatalErrorInFunction("SolveSM", error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode quickUpdateWeightsAndAbscissi(SMObj_ *smObject)
{

    mesh_          *mesh = smObject->sm[0]->access->mesh; //just needs to use the access from any sm, sm[0] will always be available if sm flag is 1 or greater.
    DM             da = mesh->da, fda = mesh->fda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, cols, rows, idxm, idxn, n, eig, nev;

    PetscReal     ***nvert;
    PetscReal     value, *a, *An, *Bn, *Vec0, *Vec1, *Vec2, *Vec3, *Vec4, *Vec5;

    PetscReal     ***weight0, ***absc0, ***weight1, ***absc1, ***weight2, ***absc2; // **Pmat;
    PetscReal     ***sm0, ***sm1, ***sm2, ***sm3, ***sm4, ***sm5;

    PetscScalar   p1, p2, p3, p4, eigVal, eigVecVal, *array;

    Vec an, bn, alpha, vec0, vec1, vec2, vec3, vec4, vec5;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecCreateSeq(MPI_COMM_SELF, 6, &alpha);
    VecCreateSeq(MPI_COMM_SELF, 3, &an);
    VecCreateSeq(MPI_COMM_SELF, 3, &bn);

    VecCreateSeq(MPI_COMM_SELF, 7, &vec0);
    VecCreateSeq(MPI_COMM_SELF, 7, &vec1);
    VecCreateSeq(MPI_COMM_SELF, 7, &vec2);
    VecCreateSeq(MPI_COMM_SELF, 7, &vec3);
    VecCreateSeq(MPI_COMM_SELF, 7, &vec4);
    VecCreateSeq(MPI_COMM_SELF, 7, &vec5);

    VecGetArray(vec0, &Vec0);
    VecGetArray(vec1, &Vec1);
    VecGetArray(vec2, &Vec2);
    VecGetArray(vec3, &Vec3);
    VecGetArray(vec4, &Vec4);
    VecGetArray(vec5, &Vec5);

    VecGetArray(alpha, &a);
    VecGetArray(an, &An);
    VecGetArray(bn, &Bn);

    DMDAVecGetArray(da, smObject->sm[0]->smVal, &sm0);
    DMDAVecGetArray(da, smObject->sm[1]->smVal, &sm1);
    DMDAVecGetArray(da, smObject->sm[2]->smVal, &sm2);
    DMDAVecGetArray(da, smObject->sm[3]->smVal, &sm3);
    DMDAVecGetArray(da, smObject->sm[4]->smVal, &sm4);
    DMDAVecGetArray(da, smObject->sm[5]->smVal, &sm5);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    //loop to all cells and set weight and abscissi
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {

                if (sm0[k][j][i] <= 1e-3 || isIBMSolidCell(k, j, i, nvert))
                {
                    weight0[k][j][i] = 0.;
                    weight1[k][j][i] = 0.;
                    weight2[k][j][i] = 0.;

                    absc0[k][j][i] = 0.;
                    absc1[k][j][i] = 0.;
                    absc2[k][j][i] = 0.;

                    continue;
                }

                //loop to set P(rows,cols) at each cell
                for (cols=0; cols<7; cols++)
                {
                    if (cols == 0)
                    {
                        Vec0[cols] = 1.0;
                        Vec1[cols] = 0.0;
                        Vec2[cols] = 0.0;
                        Vec3[cols] = 0.0;
                        Vec4[cols] = 0.0;
                        Vec5[cols] = 0.0;
                    }

                    if (cols == 1)
                    {
                        Vec0[cols] = sm0[k][j][i];
                        Vec1[cols] = -sm1[k][j][i];
                        Vec2[cols] = sm2[k][j][i];
                        Vec3[cols] = -sm3[k][j][i];
                        Vec4[cols] = sm4[k][j][i];
                        Vec5[cols] = -sm5[k][j][i];
                    }

                    if (cols > 1)
                    {
                        Vec0[cols] = Vec0[cols-1]*Vec1[cols-2] - Vec0[cols-2]*Vec1[cols-1];
                        Vec1[cols] = Vec0[cols-1]*Vec2[cols-2] - Vec0[cols-2]*Vec2[cols-1];
                        Vec2[cols] = Vec0[cols-1]*Vec3[cols-2] - Vec0[cols-2]*Vec3[cols-1];
                        Vec3[cols] = Vec0[cols-1]*Vec4[cols-2] - Vec0[cols-2]*Vec4[cols-1];
                        Vec4[cols] = Vec0[cols-1]*Vec5[cols-2] - Vec0[cols-2]*Vec5[cols-1];
                        Vec5[cols] = 0.0;
                    }

                }

                //next, find alphas, An, Bn, and then set Jacobian at each cell.

                //set alpha
                for (rows = 0; rows < 6; ++rows)
                {
                    if (rows == 0)
                    {
                        a[rows] = 0;
                        //PetscPrintf(PETSC_COMM_WORLD, "\nalpha values: %f\n", a[rows]);
                    }
                    else
                    {
                        //MatGetValue(Pmat, 0, rows+1, &p1);
                        //p1 = Pmat[0][rows+1];
                        p1 = Vec0[rows + 1];
                        //MatGetValue(Pmat, 0, rows, &p2);
                        //p2 = Pmat[0][rows];
                        p2 = Vec0[rows];
                        //MatGetValue(Pmat, 0, rows-1, &p3);
                        //p3 = Pmat[0][rows-1];
                        p3 = Vec0[rows-1];

                        if (p3 == 0)
                        {
                            a[rows] == 0;
                        }
                        else
                        {
                            a[rows] = p1/(p2*p3);
                            //PetscPrintf(PETSC_COMM_WORLD, "\nalpha values: %f\n", a[rows]);
                        }


                    }

                }

                //set An and Bn
                for (n = 0; n < 3; ++n)
                {
                    An[n] = a[2*n+1] + a[2*n];

                    if (n < 2)
                    {
                        Bn[n] = sqrt(abs(a[2*n+2]*a[2*n+1]));
                    }


                }

                //PetscPrintf(PETSC_COMM_WORLD, "\nJ values: %f, %f, %f, %f, %f\n", An[0], An[1], An[2], Bn[0], Bn[1]);

                Eigen::Matrix3d A;
                A << An[0], Bn[0], 0,
                Bn[0], An[1], Bn[1],
                0, Bn[1], An[2];

                Eigen::EigenSolver<Eigen::Matrix3d> solver(A);
                Eigen::Vector3d realEigenvalues = solver.eigenvalues().real();       // Real parts of eigenvalues
                Eigen::Matrix3d realEigenvectors = solver.eigenvectors().real();     // Real parts of eigenvectors

                for (PetscInt eig = 0; eig < realEigenvectors.cols(); ++eig)
                {

                    eigVecVal = realEigenvectors(0, eig);
                    eigVal = realEigenvalues(eig);

                    if (PetscIsInfOrNanReal(eigVecVal))
                    {
                        char error[512];
                        sprintf(error, "EigenVecvalue %li is NaN or Inf: %f at %li %li %li", eig, eigVecVal, k, j, i);
                        fatalErrorInFunction("UpdateWeightsAndAbscissi", error);
                    }

                    if (PetscIsInfOrNanReal(eigVal))
                    {
                        char error[512];
                        sprintf(error, "EigenValue %li is NaN or Inf: %f at %li %li %li", eig, eigVal, k, j, i);
                        fatalErrorInFunction("UpdateWeightsAndAbscissi", error);
                    }

                    if (eig == 0)
                    {
                        absc0[k][j][i] = eigVal;
                        weight0[k][j][i] = pow(eigVecVal, 2) * Vec0[1];
                        //printf("A0=%f, w0=%f\n",  absc0[k][j][i], weight0[k][j][i]);
                    }
                    else if (eig ==1)
                    {
                        absc1[k][j][i] = eigVal;
                        weight1[k][j][i] = pow(eigVecVal, 2) * Vec0[1];
                        //printf("A1=%f, w1=%f\n",  absc1[k][j][i], weight1[k][j][i]);
                    }
                    else if (eig == 2)
                    {
                        absc2[k][j][i] = eigVal;
                        weight2[k][j][i] = pow(eigVecVal, 2) * Vec0[1];
                        //printf("A2=%f, w2=%f\n",  absc2[k][j][i], weight2[k][j][i]);
                    }
                }


            }
        }
    }

    VecRestoreArray(vec0, &Vec0);
    VecRestoreArray(vec1, &Vec1);
    VecRestoreArray(vec2, &Vec2);
    VecRestoreArray(vec3, &Vec3);
    VecRestoreArray(vec4, &Vec4);
    VecRestoreArray(vec5, &Vec5);

    VecRestoreArray(alpha, &a);
    VecRestoreArray(an, &An);
    VecRestoreArray(bn, &Bn);

    DMDAVecRestoreArray(da, smObject->sm[0]->smVal, &sm0);
    DMDAVecRestoreArray(da, smObject->sm[1]->smVal, &sm1);
    DMDAVecRestoreArray(da, smObject->sm[2]->smVal, &sm2);
    DMDAVecRestoreArray(da, smObject->sm[3]->smVal, &sm3);
    DMDAVecRestoreArray(da, smObject->sm[4]->smVal, &sm4);
    DMDAVecRestoreArray(da, smObject->sm[5]->smVal, &sm5);
    //PetscPrintf(PETSC_COMM_WORLD, "mesh POST LOOP\n");

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    VecDestroy(&alpha);
    VecDestroy(&an);
    VecDestroy(&bn);

    VecDestroy(&vec0);
    VecDestroy(&vec1);
    VecDestroy(&vec2);
    VecDestroy(&vec3);
    VecDestroy(&vec4);
    VecDestroy(&vec5);


    return 0;
}

//***************************************************************************************************************//

PetscErrorCode DissipationCalc(SMObj_ *smObject)
{
    mesh_  *mesh   = smObject->sm[0]->access->mesh;
    flags_ *flags  = smObject->sm[0]->access->flags;

    ueqn_  *ueqn   = smObject->sm[0]->access->ueqn;
    les_   *les    = NULL;

    DMDALocalInfo info = mesh->info;
    DM            da = mesh->da, fda = mesh->fda; //sda = mesh->sda;

    PetscInt       xs = info.xs, xe = info.xs + info.xm;
    PetscInt       ys = info.ys, ye = info.ys + info.ym;
    PetscInt       zs = info.zs, ze = info.zs + info.zm;
    PetscInt       mx = info.mx, my = info.my, mz = info.mz;

    PetscInt       i, j, k;
    PetscInt       lxs, lxe, lys, lye, lzs, lze;

    PetscReal      ***nut, ***nvert, ***aj;

    PetscReal      TauijSij;

    PetscReal      tau11_SGS, tau12_SGS, tau13_SGS,
                   tau21_SGS, tau22_SGS, tau23_SGS,
                   tau31_SGS, tau32_SGS, tau33_SGS;

    PetscReal      dudc, dvdc, dwdc,
                   dude, dvde, dwde,
                   dudz, dvdz, dwdz;

    PetscReal      du_dx, du_dy, du_dz,
                   dv_dx, dv_dy, dv_dz,
                   dw_dx, dw_dy, dw_dz;

    PetscReal      csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc,
                   aN, m1, m2;

    Cmpnts         ***ucat, ***csi,  ***eta,   ***zet;

    PetscReal      ***keeps;

    PetscReal      nu = smObject->sm[0]->access->constants->nu;
    PetscReal      nuEff = nu;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    // get solution arrays
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);

    DMDAVecGetArray(da,  smObject->Eps,     &keeps);

    if(flags->isLesActive)
    {
        les = smObject->sm[0]->access->les;

        DMDAVecGetArray(da, les->lNu_t, &nut);
    }


    // compute averaging weights
    aN = (PetscReal)smObject->dissWeight; //update to smObject version
    m1 = aN  / (aN + 1.0);
    m2 = 1.0 / (aN + 1.0);

    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {

                csi0 = csi[k][j][i].x,
                csi1 = csi[k][j][i].y,
                csi2 = csi[k][j][i].z;
                eta0 = eta[k][j][i].x,
                eta1 = eta[k][j][i].y,
                eta2 = eta[k][j][i].z;
                zet0 = zet[k][j][i].x,
                zet1 = zet[k][j][i].y,
                zet2 = zet[k][j][i].z;
                ajc  = aj[k][j][i];

                Compute_du_center
                (
                    mesh,
                    i, j, k, mx, my, mz, ucat, nvert, &dudc,
                    &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz
                );

                Compute_du_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                    zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                    dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                    &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                // Tau_ij
                if(flags->isLesActive) {nuEff += nut[k][j][i];}

                tau11_SGS = - nuEff*(2.0*du_dx);
                tau12_SGS = - nuEff*(du_dy + dv_dx);
                tau13_SGS = - nuEff*(du_dz + dw_dx);
                tau21_SGS =   tau12_SGS;
                tau22_SGS = - nuEff*(2.0*dv_dy);
                tau23_SGS = - nuEff*(dv_dz + dw_dy);
                tau31_SGS =   tau13_SGS;
                tau32_SGS =   tau23_SGS;
                tau33_SGS = - nuEff*(2.0*dw_dz);

                // dissipation calculation
                TauijSij  = 0.5 *
                (
                    tau11_SGS*(2.0*du_dx) + tau12_SGS*(du_dy + dv_dx) + tau13_SGS*(du_dz + dw_dx) +
                    tau21_SGS*(du_dy + dv_dx) + tau22_SGS*(2.0*dv_dy) + tau23_SGS*(dv_dz + dw_dy) +
                    tau31_SGS*(du_dz + dw_dx) + tau32_SGS*(dv_dz + dw_dy) + tau33_SGS*(2.0*dw_dz)
                );

                // save dissipation after time averaging
                keeps[k][j][i] = m1 * keeps[k][j][i] - m2 * TauijSij;

                if (i==3 && j==3 && k ==3)
                {
                    PetscPrintf(PETSC_COMM_WORLD, "Keps: %f\n", keeps[k][j][i]);
                }

            }
        }
    }

    smObject->dissWeight++;

    // restore solution arrays
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);

    if(flags->isLesActive)
    {
        les = smObject->sm[0]->access->les;
        DMDAVecRestoreArray(da, les->lNu_t, &nut);
    }


    DMDAVecRestoreArray(da,  smObject->Eps,     &keeps);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode infectProb(SMObj_ *smObject)
{
    mesh_          *mesh = smObject->sm[0]->access->mesh; //just needs to use the access from any sm, sm[0] will always be available if sm flag is 1 or greater.
    //flags_         *flags = smObject->sm[0]->access->flags;
    ibm_           *ibm = smObject->sm[0]->access->ibm;
    clock_         *clock = smObject->sm[0]->access->clock;
    vents_         *vents = smObject->sm[0]->access->vents;

    DM             da = mesh->da, fda = mesh->fda;// daCoag = smObject->daCoag;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***sm, ***sm_1, ***sm_2, ***sm_3, ***sm_4, ***sm_5;

    PetscReal     ***Quant, ***dq, ***ProbI, ***nvert, meanD;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscInt      ***markVent;

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    DMDAVecGetArray(da, smObject->quant, &Quant);
    DMDAVecGetArray(da, smObject->Dq, &dq);
    DMDAVecGetArray(da, smObject->probI, &ProbI);

    DMDAVecGetArray(da, smObject->sm[0]->smVal, &sm);
    DMDAVecGetArray(da, smObject->sm[1]->smVal, &sm_1);
    DMDAVecGetArray(da, smObject->sm[2]->smVal, &sm_2);
    DMDAVecGetArray(da, smObject->sm[3]->smVal, &sm_3);
    DMDAVecGetArray(da, smObject->sm[4]->smVal, &sm_4);
    DMDAVecGetArray(da, smObject->sm[5]->smVal, &sm_5);


    //loop to all cells
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                //don't apply any sources at vents
               if (sm[k][j][i] <= 1e-6)
               {
                    continue;
               }

               meanD = sm_1[k][j][i]/sm[k][j][i];
               Quant[k][j][i] = sm[k][j][i]*smObject->OGConc*meanD*meanD*meanD*(M_PI/6)*pow(10, 7)*pow(10, -12)/130/210;
               dq[k][j][i] += 0.0004166667*Quant[k][j][i]*clock->dt;
               ProbI[k][j][i] = 1 - exp(-dq[k][j][i]);

               if (sm[k][j][i] > 1e-6)
               {
                    //printf("%f, %f, %f, %f, %f, %f\n", meanD, sm[k][j][i], smObject->OGConc, 10000000000*Quant[k][j][i], dq[k][j][i], ProbI[k][j][i]);
               }

            }
        }
    }

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    DMDAVecRestoreArray(da, smObject->sm[0]->smVal, &sm);
    DMDAVecRestoreArray(da, smObject->sm[1]->smVal, &sm_1);
    DMDAVecRestoreArray(da, smObject->sm[2]->smVal, &sm_2);
    DMDAVecRestoreArray(da, smObject->sm[3]->smVal, &sm_3);
    DMDAVecRestoreArray(da, smObject->sm[4]->smVal, &sm_4);
    DMDAVecRestoreArray(da, smObject->sm[5]->smVal, &sm_5);

    DMDAVecRestoreArray(da, smObject->quant, &Quant);
    DMDAVecRestoreArray(da, smObject->Dq, &dq);
    DMDAVecRestoreArray(da, smObject->probI, &ProbI);


    return 0;
}

//***************************************************************************************************************//

PetscErrorCode partCount(SMObj_ *smObject)
{
    mesh_          *mesh = smObject->sm[0]->access->mesh; //just needs to use the access from any sm, sm[0] will always be available if sm flag is 1 or greater.
    //flags_         *flags = smObject->sm[0]->access->flags;
    ibm_           *ibm = smObject->sm[0]->access->ibm;
    clock_         *clock = smObject->sm[0]->access->clock;
    vents_         *vents = smObject->sm[0]->access->vents;
    ueqn_          *ueqn = smObject->sm[0]->access->ueqn;

    DM             da = mesh->da, fda = mesh->fda;// daCoag = smObject->daCoag;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***sm, ***nvert, ***depCount, ***exCount, ***aj;

    PetscInt      ***markVent;

    Cmpnts        ***lsed, ***dep, ***ucont;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    DMDAVecGetArray(da, smObject->DepCount, &depCount);
    DMDAVecGetArray(da, smObject->ExCount, &exCount);

    DMDAVecGetArray(fda, smObject->sm[0]->lSed, &lsed);
    DMDAVecGetArray(fda, smObject->sm[0]->Dep, &dep);
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);

    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    DMDAVecGetArray(da, smObject->sm[0]->smVal, &sm);

    //loop to all cells
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                //particles deposited due to sedimatation (sed) + difusion and impaction (dep)
                if (i==1 || isIBMSolidCell(k, j, i-1, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].y * clock->dt);
                    // -                  (- lsed[k  ][j  ][i-1].y * aj[k][j][i] * clock->dt);

                }

                if (i==mx-2 || isIBMSolidCell(k, j, i+1, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].y * clock->dt);
                    //-                    (lsed[k  ][j  ][i].y * aj[k][j][i] * clock->dt);

                }

                if (j==1 || isIBMSolidCell(k, j-1, i, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].z * clock->dt);
                    //-                  (- lsed[k  ][j-1][i].z * aj[k][j][i] * clock->dt);

                }

                if (j==my-2 || isIBMSolidCell(k, j+1, i, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].z * clock->dt) -
                    (lsed[k  ][j][i].z * aj[k][j][i] * clock->dt);

                }

                if (k==1 || isIBMSolidCell(k-1, j, i, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].x * clock->dt);
                    //-                  (- lsed[k-1][j][i].x * aj[k][j][i] * clock->dt);

                }

                if (k==mz-2 || isIBMSolidCell(k+1, j, i, nvert))
                {
                    depCount[k][j][i]
                    +=
                    (dep[k][j][i].x * clock->dt);
                    //-                (- lsed[k][j][i].x * aj[k][j][i] * clock->dt);

                }

                if (markVent[k][j][i] > 0)
                {
                    PetscInt q = markVent[k][j][i] - 1;

                    if (vents->vent[q]->face == "iLeft" || vents->vent[q]->face == "iRight")
                    {
                        exCount[k][j][i] += -sm[k][j][i] * ucont[k][j][i].x * aj[k][j][i] * clock->dt;
                    }

                    if (vents->vent[q]->face == "jLeft" || vents->vent[q]->face == "jRight")
                    {
                        exCount[k][j][i] += -sm[k][j][i] * ucont[k][j][i].y * aj[k][j][i] * clock->dt;
                    }

                    if (vents->vent[q]->face == "kLeft" || vents->vent[q]->face == "kRight")
                    {
                        exCount[k][j][i] += -sm[k][j][i] * ucont[k][j][i].z * aj[k][j][i] * clock->dt;
                    }

                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    DMDAVecRestoreArray(da, smObject->DepCount, &depCount);
    DMDAVecRestoreArray(da, smObject->ExCount, &exCount);

    DMDAVecRestoreArray(fda, smObject->sm[0]->lSed, &lsed);
    DMDAVecRestoreArray(fda, smObject->sm[0]->Dep, &dep);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    DMDAVecRestoreArray(da, smObject->sm[0]->smVal, &sm);


    return 0;
}

//***************************************************************************************************************//

PetscErrorCode formCoagSourceExp(sm_ *sm, PetscInt ii)
{
    mesh_          *mesh = sm->access->mesh; //just needs to use the access from any sm, sm[0] will always be available if sm flag is 1 or greater.
    flags_         *flags = sm->access->flags;
    peqn_          *peqn = sm->access->peqn;
    ueqn_          *ueqn = sm->access->ueqn;
    ibm_           *ibm = sm->access->ibm;
    teqn_          *teqn = sm->access->teqn;
    constants_     *cst   = sm->access->constants;
    clock_         *clock = sm->access->clock;
    SMObj_         *smObject = sm->access->smObject;

    vents_         *vents = sm->access->vents;

    DM             da = mesh->da, fda = mesh->fda;// daCoag = smObject->daCoag;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, cols, rows, idxm, idxn, n, eig, nev;

    PetscReal     ***sm_o4, ***sm_o3, ***sm_o;
    PetscReal     ***weight, ***absc, ***weight1, ***absc1, ***weight2, ***absc2;
    PetscReal     ***pres, ***keeps, ***tmprt;

    PetscScalar   a1, a2, a3, w1, w2, k_b, mfp, meanDia,
                  sigma, t_kol, Kn, Cu, D_b, y1, m1, m2, m, D_f,
                  term1, term2, term3, term4, power, temp, kk, a1f, a2f;

    PetscScalar   c0, c1, c2, c3, c4, c5, c6, c7, m3Dia;
    PetscScalar   ap0, ap1, ap2;
    PetscScalar   cap0, cap1, cap2;
    PetscScalar   c6p, expMfp;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscReal     ***coag;
    Cmpnts	      ***ucat, ***cent;
    PetscReal     ***nvert;
    PetscInt      ***markVent;

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    DMDAVecGetArray(da, peqn->P, &pres); //just using 101kPA for now, figure out units before implemnting pres.
    DMDAVecGetArray(da, smObject->Eps, &keeps);
    DMDAVecGetArray(da, smObject->sm[0]->lsmVal, &sm_o);
    DMDAVecGetArray(da, smObject->sm[3]->lsmVal, &sm_o3);
    DMDAVecGetArray(da, smObject->sm[4]->lsmVal, &sm_o4);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->weight, &weight);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->absc, &absc);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->absc, &absc2);

    VecSet(sm->coagSource, 0.0);
    DMDAVecGetArray(da, sm->coagSource, &coag);

    k_b = (1.380649E-23); //*1E18; //boltzman constant in kg*nm^2/(K*s^2)
    D_f = 2.2;

    c0 = smObject->monoDFrac/2;
    c1 = (0.414*(D_f)-0.211)/smObject->monoDFrac;
    c2 = 0;
    c6 = 0;
    c4 = sqrt(8*M_PI)*smObject->OGConc*pow(10, -18);
    c7 = M_PI/2;

    //loop to all cells
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                //don't apply any sources at vents
                if (markVent[k][j][i] > 0 || isIBMCell(k, j, i, nvert))
                {
                    continue;
                }

                if (weight[k][j][i] <=0.0 || weight1[k][j][i] <= 0.0 || weight2[k][j][i] <= 0.0 || absc[k][j][i] <= 0.0 || absc1[k][j][i] <= 0.0 || absc2[k][j][i] <= 0.0)
                {
                    continue;
                }

                if (flags->isTeqnActive)
                {
                    DMDAVecGetArray(da, teqn->Tmprt, &tmprt);
                    temp = tmprt[k][j][i];
                    DMDAVecRestoreArray(da, teqn->Tmprt, &tmprt);
                }
                else
                {
                    temp = cst->tRef;
                }

                mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                if (keeps[k][j][i] != 0)
                {
                    t_kol = sqrt(cst->nu/abs(keeps[k][j][i]));
                    c2 = c0*c0/(15*t_kol*t_kol);
                    c6 = (0.15 * 3.0 * M_PI * cst->nu * cst->rho)/(t_kol*k_b*temp*pow(10, 18));
                }

                c3 = 6*k_b*temp*pow(10, 30)/(smObject->rhoPart*M_PI);

                c5 = c4*c0*3*M_PI*cst->nu*cst->rho/(4*M_PI*k_b*temp*pow(10, 18));

                m3Dia = sm_o4[k][j][i]/sm_o3[k][j][i]; //volume mean diameter of previous it.

                ap0 = PetscPowScalar(absc[k][j][i], 3.);
                ap1 = PetscPowScalar(absc1[k][j][i], 3.);
                ap2 = PetscPowScalar(absc2[k][j][i], 3.);
                cap0 = PetscPowScalar(c1*absc[k][j][i]/(m3Dia), (3.0/D_f));
                cap1 = PetscPowScalar(c1*absc1[k][j][i]/(m3Dia), (3.0/D_f));
                cap2 = PetscPowScalar(c1*absc2[k][j][i]/(m3Dia), (3.0/D_f));
                expMfp = exp(-0.88/(2*mfp/(m3Dia)));
                c6p = PetscPowScalar(c6/((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia)), 1./2.);


                //set coag source value
                coag[k][j][i] =

                (0.5*PetscPowScalar((ap0 + ap0), ((PetscScalar)ii)/3.) - PetscPowScalar(absc[k][j][i], ((PetscScalar)ii))) * weight[k][j][i]*weight[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap0, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap0, 2.0)) + c3*(ap0 + ap0)/(ap0*ap0), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap0 + cap0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap0, 2.0)) + c3*(ap0 + ap0)/(ap0*ap0), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap0 + cap0)) * c6p + c0 * ((m3Dia)*(cap0 + cap0)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap0 + cap0)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap0 + ap1), ((PetscScalar)ii)/3.) - PetscPowScalar(absc[k][j][i], ((PetscScalar)ii)))*weight[k][j][i]*weight1[k][j][i]* c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap1, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap1, 2.0)) + c3*(ap0 + ap1)/(ap0*ap1), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap0 + cap1)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap1, 2.0)) + c3*(ap0 + ap1)/(ap0*ap1), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap0 + cap1)) * c6p + c0 * ((m3Dia)*(cap0 + cap1)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap0 + cap1)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap0 + ap2), ((PetscScalar)ii)/3.) - PetscPowScalar(absc[k][j][i], ((PetscScalar)ii)))*weight[k][j][i]*weight2[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap2, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap2, 2.0)) + c3*(ap0 + ap2)/(ap0*ap2), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap0 + cap2)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap0 + cap2, 2.0)) + c3*(ap0 + ap2)/(ap0*ap2), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap0 + cap2)) * c6p + c0 * ((m3Dia)*(cap0 + cap2)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap0 + cap2)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap1 + ap0), ((PetscScalar)ii)/3.) - PetscPowScalar(absc1[k][j][i], ((PetscScalar)ii)))*weight1[k][j][i]*weight[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap0, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap0, 2.0)) + c3*(ap1 + ap0)/(ap1*ap0), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap1 + cap0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap0, 2.0)) + c3*(ap1 + ap0)/(ap1*ap0), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap1 + cap0)) * c6p + c0 * ((m3Dia)*(cap1 + cap0)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap1 + cap0)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap1 + ap1), ((PetscScalar)ii)/3.) - PetscPowScalar(absc1[k][j][i], ((PetscScalar)ii)))*weight1[k][j][i]*weight1[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap1, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap1, 2.0)) + c3*(ap1 + ap1)/(ap1*ap1), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap1 + cap1)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap1, 2.0)) + c3*(ap1 + ap1)/(ap1*ap1), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap1 + cap1)) * c6p + c0 * ((m3Dia)*(cap1 + cap1)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap1 + cap1)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap1 + ap2), ((PetscScalar)ii)/3.) - PetscPowScalar(absc1[k][j][i], ((PetscScalar)ii)))*weight1[k][j][i]*weight2[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap2, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap2, 2.0)) + c3*(ap1 + ap2)/(ap1*ap2), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap1 + cap2)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap1 + cap2, 2.0)) + c3*(ap1 + ap2)/(ap1*ap2), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap1 + cap2)) * c6p + c0 * ((m3Dia)*(cap1 + cap2)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap1 + cap2)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap2 + ap0), ((PetscScalar)ii)/3.) - PetscPowScalar(absc2[k][j][i], ((PetscScalar)ii)))*weight2[k][j][i]*weight[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap0, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap0, 2.0)) + c3*(ap2 + ap0)/(ap2*ap0), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap2 + cap0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap0, 2.0)) + c3*(ap2 + ap0)/(ap2*ap0), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap2 + cap0)) * c6p + c0 * ((m3Dia)*(cap2 + cap0)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap2 + cap0)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap2 + ap1), ((PetscScalar)ii)/3.) - PetscPowScalar(absc2[k][j][i], ((PetscScalar)ii)))*weight2[k][j][i]*weight1[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap1, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap1, 2.0)) + c3*(ap2 + ap1)/(ap2*ap1), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap2 + cap1)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap1, 2.0)) + c3*(ap2 + ap1)/(ap2*ap1), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap2 + cap1)) * c6p + c0 * ((m3Dia)*(cap2 + cap1)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap2 + cap1)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  ) -

                (0.5*PetscPowScalar((ap2 + ap2), ((PetscScalar)ii)/3.) - PetscPowScalar(absc2[k][j][i], ((PetscScalar)ii)))*weight2[k][j][i]*weight2[k][j][i] * c4*c0*c0 * ((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap2, 2.0)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap2, 2.0)) + c3*(ap2 + ap2)/(ap2*ap2), 1./2.) /
                (1 + c5 * ((m3Dia)*(cap2 + cap2)) *
                PetscPowScalar(c2*((m3Dia)*(m3Dia)*PetscPowScalar(cap2 + cap2, 2.0)) + c3*(ap2 + ap2)/(ap2*ap2), 1./2.) *
                (1 - c7*c0 * ((m3Dia)*(cap2 + cap2)) * c6p + c0 * ((m3Dia)*(cap2 + cap2)) *
                c6p *
                atan(c0 * ((m3Dia)*(cap2 + cap2)) * c6p)) /
                ((1 + 2*mfp*(1.2+0.4*expMfp)/(m3Dia))/(m3Dia))  );


                if (coag[k][j][i] != 0)
                {
                   //printf("coag = %f\n", coag[k][j][i]);
                }
            }
        }
    }


    DMDAVecRestoreArray(da, sm->coagSource, &coag);


    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);

    DMDAVecRestoreArray(da, peqn->P, &pres);
    DMDAVecRestoreArray(da, smObject->Eps, &keeps);
    DMDAVecRestoreArray(da, smObject->sm[0]->lsmVal, &sm_o);
    DMDAVecRestoreArray(da, smObject->sm[3]->lsmVal, &sm_o3);
    DMDAVecRestoreArray(da, smObject->sm[4]->lsmVal, &sm_o4);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->weight, &weight);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->absc, &absc);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode formDepSourceExp(sm_ *sm, PetscInt ii)
{
    mesh_          *mesh = sm->access->mesh; //just needs to use the access from any sm, sm[0] will always be available if sm flag is 1 or greater.
    flags_         *flags = sm->access->flags;
    peqn_          *peqn = sm->access->peqn;
    ueqn_          *ueqn = sm->access->ueqn;
    ibm_           *ibm = sm->access->ibm;
    teqn_          *teqn = sm->access->teqn;
    constants_     *cst   = sm->access->constants;
    clock_         *clock = sm->access->clock;
    SMObj_         *smObject = sm->access->smObject;

    vents_         *vents = sm->access->vents;

    DM             da = mesh->da, fda = mesh->fda;// daCoag = smObject->daCoag;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    //depostion variables
    PetscReal     du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

    PetscReal     dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz;
    PetscReal     dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords

    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc;      // surface area vectors components

    Cmpnts	      ***icsi, ***ieta, ***izet;
    Cmpnts	      ***jcsi, ***jeta, ***jzet;
    Cmpnts	      ***kcsi, ***keta, ***kzet;

    Cmpnts	      ***ucat, ***cent, ustar, ***dep;

    PetscReal     ***iaj, ***jaj, ***kaj, ***aj, ***nvert;

    PetscReal     tauWall, gradNorm, yPlus, temp;
    PetscReal     Ip0, Ip1, Ip2, uDep0, uDep1, uDep2;
    PetscReal     TauP0, TauP1, TauP2, TauP0Plus, TauP1Plus, TauP2Plus;
    PetscReal     Cu0, Cu1, Cu2, expMfp0, expMfp1, expMfp2;
    PetscReal     gPlus0, gPlus1, gPlus2, ustarCont, ustarMag;

    PetscReal     ***sm_o;
    PetscReal     ***weight0, ***absc0, ***weight1, ***absc1, ***weight2, ***absc2;
    PetscReal     ***tmprt;

    PetscScalar   mfp;
    PetscReal     k_b = 1.380649E-23; //*1E18; //boltzman constant in kg*nm^2/(K*s^2)

    PetscInt      i1, j1, k1;

    PetscReal Karman = 0.41; PetscReal lamda0 = 13.7; PetscReal lamda1 = 2./3.; //empircal constants for ventilated rooms.
    PetscReal ome = 1700.;

    PetscInt      ***markVent;

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    DMDAVecGetArray(da, sm->lsmVal, &sm_o);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->absc, &absc2);

    VecSet(smObject->sm[ii]->Dep, 0.0);
    DMDAVecGetArray(fda, smObject->sm[ii]->Dep, &dep);

    DMDAVecGetArray(da, mesh->lAj, &aj);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);

    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);

    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    //dep at boundary and IBMSolid Face
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {

                //don't apply any depsition at vents or solid IBm cells
               if (markVent[k][j][i] > 0 || isIBMSolidCell(k, j, i, nvert))
               {
                    continue;
               }

               if (weight0[k][j][i] <=0.0 || weight1[k][j][i] <= 0.0 || weight2[k][j][i] <= 0.0 || absc0[k][j][i] <= 0.0 || absc1[k][j][i] <= 0.0 || absc2[k][j][i] <= 0.0)
               {
                   continue;
               }

               //i-face dep.y (Y-face in paraview)
               if ((i==1 && mesh->boundaryU.iLeft == "noSlip") || (i==mx-2 && mesh->boundaryU.iRight == "noSlip")) //formula for deposition and wall shear varies for vertical vs. horizontal walls.
               {

                   csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                   eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                   zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
                   ajc  = iaj[k][j][i];

                   if (i==1) // to ensure derivitative it between boundary and first fluid cell.
                   {
                       i1 = i-1; j1 = j; k1 = k;
                   }
                   else
                   {
                       i1 = i; j1 = j; k1 = k;
                   }

                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_i
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );

                   // compute cartesian velocity derivatives w.r.t cartesian coords
                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates i-face is y-cart _dy used
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dy), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dy), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dy), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at i-face
                   =
                   (
                       ustar.x * icsi[k][j][i].x +
                       ustar.y * icsi[k][j][i].y +
                       ustar.z * icsi[k][j][i].z
                   );

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   //note all WaA value are at center, but since the boundary condition at all locations where Deposition is applicable is ZG, this is equivalent to face values/
                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   //depositon regimes from nerisson et al.
                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].y - cent[k][j][i+1].y)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   uDep0 = 1.0/Ip0;
                   uDep1 = 1.0/Ip1;
                   uDep2 = 1.0/Ip2;

                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].y +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }
               //deposition from cell beside IBM solid i-face
               if (isIBMSolidCell(k, j, i+1, nvert) || isIBMSolidCell(k, j, i-1, nvert))
               {
                   csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                   eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                   zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
                   ajc  = iaj[k][j][i];

                   if (isIBMSolidCell(k, j, i+1, nvert)) // to ensure derivitative it between ibmSolid and fluid.
                   {
                       i1 = i; j1 = j; k1 = k;
                   }
                   else
                   {
                       i1 = i-1; j1 = j; k1 = k;
                   }

                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_i
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );

                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates iface is _dy
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dy), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dy), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dy), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at i-face
                   =
                   (
                       ustar.x * icsi[k][j][i].x +
                       ustar.y * icsi[k][j][i].y +
                       ustar.z * icsi[k][j][i].z
                   );

                   //printf("\n%f, %li, %li, %li, %f, %f, %f\n", ustarCont, k, j, i, ustar.x, ustar.y, ustar.z);

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   //note all WaA value are at center, but since the boundary condition at all locations where Deposition is applicable is ZG, this is equivalent to face values/
                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].y - cent[k][j][i+1].y)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   uDep0 = 1.0/Ip0;
                   uDep1 = 1.0/Ip1;
                   uDep2 = 1.0/Ip2;

                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].y +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }

               //j-face dep.y (Zface in Paraview)
               if ((j==1 && mesh->boundaryU.jLeft == "noSlip") || (j==my-2 && mesh->boundaryU.jRight == "noSlip")) //formula for deposition and wall shear varies for vertical vs. horizontal walls.
               {
                   csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                   eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                   zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
                   ajc  = jaj[k][j][i];

                   if (j==1) // to ensure derivitative it between boundary and first fluid cell.
                   {
                       i1 = i; j1 = j-1; k1 = k;
                   }
                   else
                   {
                       i1 = i; j1 = j; k1 = k;
                   }

                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_j
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );

                   // compute cartesian velocity derivatives w.r.t cartesian coords
                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates at j-face which is _dz
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dz), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dz), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dz), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at j-face
                   =
                   (
                       ustar.x * jeta[k][j][i].x +
                       ustar.y * jeta[k][j][i].y +
                       ustar.z * jeta[k][j][i].z
                   );

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].z - cent[k][j+1][i].z)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   if (ustarMag > 0)
                   {
                       gPlus0 = 9.81*TauP0/ustarMag;
                       gPlus1 = 9.81*TauP1/ustarMag;
                       gPlus2 = 9.81*TauP2/ustarMag;
                   }
                   else
                   {
                       gPlus0 = 0.;
                       gPlus1 = 0.;
                       gPlus2 = 0.;
                   }

                   if (gPlus0 <= 0 || gPlus1 <= 0 || gPlus2 <= 0)
                   {
                       uDep0 = 1/Ip0;
                       uDep1 = 1/Ip1;
                       uDep2 = 1/Ip2;
                   }
                   else
                   {
                       if (j==1)
                       {
                           uDep0 = gPlus0/(1-exp(-gPlus0*Ip0));
                           uDep1 = gPlus1/(1-exp(-gPlus1*Ip1));
                           uDep2 = gPlus2/(1-exp(-gPlus2*Ip2));
                       }
                       else
                       {
                           uDep0 = -gPlus0/(1-exp(gPlus0*Ip0));
                           uDep1 = -gPlus1/(1-exp(gPlus1*Ip1));
                           uDep2 = -gPlus2/(1-exp(gPlus2*Ip2));
                       }
                   }


                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].z +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }
               //deposition from cell beside IBM solid j-face
               if (isIBMSolidCell(k, j+1, i, nvert) || isIBMSolidCell(k, j-1, i, nvert))
               {
                   csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                   eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                   zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
                   ajc  = jaj[k][j][i];

                   if (isIBMSolidCell(k, j+1, i, nvert)) // to ensure derivitative it between boundary and first fluid cell.
                   {
                       i1 = i; j1 = j; k1 = k;
                   }
                   else
                   {
                       i1 = i; j1 = j-1; k1 = k;
                   }
                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_j
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );
                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates at j-face is _dz
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dz), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dz), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dz), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at j-face
                   =
                   (
                       ustar.x * jeta[k][j][i].x +
                       ustar.y * jeta[k][j][i].y +
                       ustar.z * jeta[k][j][i].z
                   );

                   //printf("\n%f, %li, %li, %li, %f, %f, %f\n", ustarCont, k, j, i, ustar.x, ustar.y, ustar.z);

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].z - cent[k][j+1][i].z)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   if (ustarMag > 0)
                   {
                       gPlus0 = 9.81*TauP0/ustarMag;
                       gPlus1 = 9.81*TauP1/ustarMag;
                       gPlus2 = 9.81*TauP2/ustarMag;
                   }
                   else
                   {
                       gPlus0 = 0.;
                       gPlus1 = 0.;
                       gPlus2 = 0.;
                   }

                   if (gPlus0 <= 0 || gPlus1 <= 0 || gPlus2 <= 0)
                   {
                       uDep0 = 1/Ip0;
                       uDep1 = 1/Ip1;
                       uDep2 = 1/Ip2;
                   }
                   else
                   {
                       if (isIBMSolidCell(k, j-1, i, nvert))
                       {
                           uDep0 = gPlus0/(1-exp(-gPlus0*Ip0));
                           uDep1 = gPlus1/(1-exp(-gPlus1*Ip1));
                           uDep2 = gPlus2/(1-exp(-gPlus2*Ip2));
                       }
                       else
                       {
                           uDep0 = -gPlus0/(1-exp(gPlus0*Ip0));
                           uDep1 = -gPlus1/(1-exp(gPlus1*Ip1));
                           uDep2 = -gPlus2/(1-exp(gPlus2*Ip2));
                       }
                   }

                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].z +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }

               //k-face dep.x (XFace in paraview)
               if ((k==1 && mesh->boundaryU.kLeft == "noSlip") || (k==mz-2 && mesh->boundaryU.kRight == "noSlip")) //formula for deposition and wall shear varies for vertical vs. horizontal walls.
               {

                   csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                   eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                   zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
                   ajc  = kaj[k][j][i];

                   if (k==1) // to ensure derivitative it between boundary and first fluid cell.
                   {
                       i1 = i; j1 = j; k1 = k-1;
                   }
                   else
                   {
                       i1 = i; j1 = j; k1 = k;
                   }

                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_k
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );

                   // compute cartesian velocity derivatives w.r.t cartesian coords
                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates at k-face is _dx
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dx), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dx), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dx), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at k-face
                   =
                   (
                       ustar.x * kzet[k][j][i].x +
                       ustar.y * kzet[k][j][i].y +
                       ustar.z * kzet[k][j][i].z
                   );

                   //printf("\n%f, %li, %li, %li, %f, %f, %f\n", ustarCont, k, j, i, ustar.x, ustar.y, ustar.z);

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0* PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].x - cent[k+1][j][i].x)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   uDep0 = 1.0/Ip0;
                   uDep1 = 1.0/Ip1;
                   uDep2 = 1.0/Ip2;

                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].x +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }
               //deposition from cell beside IBM k-face
               else if (isIBMSolidCell(k+1, j, i, nvert) || isIBMSolidCell(k-1, j, i, nvert))
               {
                   csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                   eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                   zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
                   ajc  = kaj[k][j][i];

                   if (isIBMSolidCell(k+1, j, i, nvert)) // to ensure derivitative it between boundary and first fluid cell.
                   {
                       i1 = i; j1 = j; k1 = k;
                   }
                   else
                   {
                       i1 = i; j1 = j; k1 = k-1;
                   }

                   // compute cartesian velocity derivatives w.r.t. curvilinear coords
                   Compute_du_k
                   (   mesh, i1, j1, k1, mx, my, mz, ucat, nvert,
                       &dudc, &dvdc, &dwdc,
                       &dude, &dvde, &dwde,
                       &dudz, &dvdz, &dwdz
                   );

                   // compute cartesian velocity derivatives w.r.t cartesian coords
                   Compute_du_dxyz
                   (
                       mesh,
                       csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                       zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                       dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                       &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                   );

                   //friction velocity in cartesian coordinates at k-face is _dx
                   ustar.x = PetscPowScalar(cst->nu*fabs(du_dx), 1./2.);
                   ustar.y = PetscPowScalar(cst->nu*fabs(dv_dx), 1./2.);
                   ustar.z = PetscPowScalar(cst->nu*fabs(dw_dx), 1./2.);

                   ustarMag = PetscPowScalar(PetscPowScalar(ustar.x, 2.) + PetscPowScalar(ustar.y, 2.) + PetscPowScalar(ustar.z, 2.), 1./2.);

                   ustarCont // only need ustar at k-face
                   =
                   (
                       ustar.x * kzet[k][j][i].x +
                       ustar.y * kzet[k][j][i].y +
                       ustar.z * kzet[k][j][i].z
                   );

                   //printf("\n%f, %li, %li, %li, %f, %f, %f\n", ustarCont, k, j, i, ustar.x, ustar.y, ustar.z);

                   if (sm->access->flags->isTeqnActive)
                   {
                       DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                       temp = tmprt[k][j][i]; // need to get at face value??

                       DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                   }
                   else
                   {
                       temp = cst->tRef; //need at face??
                   }

                   mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                   expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                   expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                   expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                   Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                   Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                   Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                   TauP0 = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*cst->nu*cst->rho);
                   TauP1 = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*cst->nu*cst->rho);
                   TauP2 = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*cst->nu*cst->rho);

                   TauP0Plus = (TauP0* PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP1Plus = (TauP1 * PetscPowScalar(ustarMag, 2))/cst->nu;
                   TauP2Plus = (TauP0 * PetscPowScalar(ustarMag, 2))/cst->nu;

                   if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 0.1)
                   {
                       Ip0 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip1 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                       Ip2 = 1./(PetscPowScalar(cst->Sc, -lamda1)/lamda0);
                   }
                   else if ((TauP0Plus + TauP1Plus + TauP2Plus)/3. < 10)
                   {
                       Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                   }
                   else
                   {
                       yPlus = ustarMag*0.5*fabs(cent[k][j][i].x - cent[k+1][j][i].x)/cst->nu;

                       if (yPlus > 1)
                       {
                           Ip0 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = (cst->ScT/Karman)*log(yPlus) + PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                       else
                       {
                           Ip0 = PetscPowScalar(((PetscPowScalar(TauP0Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip1 = PetscPowScalar(((PetscPowScalar(TauP1Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                           Ip2 = PetscPowScalar(((PetscPowScalar(TauP2Plus, 2)/ome)+(PetscPowScalar(cst->Sc, -lamda1)/lamda0)), -1);
                       }
                   }

                   uDep0 = 1.0/Ip0;
                   uDep1 = 1.0/Ip1;
                   uDep2 = 1.0/Ip2;

                   // depositon is only at 1 side of cell, so no subtraction is present.
                   dep[k][j][i].x +=
                   - aj[k][j][i] * ustarCont *
                   ((weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)*uDep0) +
                    (weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii)*uDep1) +
                    (weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii)*uDep2));

               }

            }
        }
    }

    DMDAVecRestoreArray(fda, sm->Dep, &dep);

    DMGlobalToLocalBegin(fda, sm->Dep,  INSERT_VALUES, sm->lDep);
    DMGlobalToLocalEnd  (fda, sm->Dep,  INSERT_VALUES, sm->lDep);

    resetFacePeriodicFluxesVector(mesh, sm->Dep,   sm->lDep, "globalToLocal");

    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);


    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);

    DMDAVecRestoreArray(da, sm->lsmVal, &sm_o);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode sourceSMCoag(sm_ *sm, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = sm->access->mesh;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***coag, ***rhs;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  sm->coagSource, &coag);
    DMDAVecGetArray(da,  Rhs,  &rhs);

    // loop over internal cell faces
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhs[k][j][i] += scale * coag[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(da,  sm->coagSource, &coag);
    DMDAVecRestoreArray(da,  Rhs,  &rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode sourceSMDep(sm_ *sm, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh = sm->access->mesh;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***rhs;
    Cmpnts        ***dep;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda,  sm->lDep, &dep);
    DMDAVecGetArray(da,  Rhs,  &rhs);

    // loop over internal cell faces
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhs[k][j][i] += scale * (dep[k][j][i].x + dep[k][j][i].y + dep[k][j][i].z);
            }
        }
    }

    DMDAVecRestoreArray(fda,  sm->lDep, &dep);
    DMDAVecRestoreArray(da,  Rhs,  &rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode sedFluxSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii)
{
    mesh_         *mesh  = sm->access->mesh;
    SMObj_        *smObject = sm->access->smObject;
    constants_    *Cst   = sm->access->constants;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***gcont, ***sed, ***lsed, ***limiter;
    Cmpnts        ***icsi, ***jeta, ***kzet, ***cent, ***ucont;
    PetscReal     ***nvert, ***aj, ***tmprt, ***rhs, ***sedCent, ***lSedCent;
    PetscReal     ***weight0, ***absc0, ***weight1, ***absc1, ***weight2, ***absc2, ***TauP0, ***TauP1, ***TauP2, ***sm0;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     temp;
    PetscReal     mfp;
    PetscReal     m3Dia;
    PetscReal     expMfp0, expMfp1, expMfp2;
    PetscReal     Cu0, Cu1, Cu2;

    PetscInt    iL, iR, jL, jR, kL, kR;
    PetscReal   denom;

    PetscInt      ***markVent;

    PetscReal     k_b = 1.380649E-23; //*1E18; //boltzman constant in kg*nm^2/(K*s^2)

    Cmpnts        gravity;
                  gravity.x = 0;
                  gravity.y = 0;
                  gravity.z = -9.81;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(fda, sm->access->ueqn->gCont,  &gcont);
    DMDAVecGetArray(fda, sm->access->ueqn->lUcont,      &ucont);

    DMDAVecGetArray(da, sm->smVal, &sm0);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    // set contravariant gravity (may turn off when teqn is on) and sed values at cell centers
    for(k=zs; k<ze; k++)
    {
        for(j=ys; j<ye; j++)
        {
            for(i=xs; i<xe; i++)
            {
                if(j > 0 && k > 0 && j < my-1 && k < mz-1)
                {
                    gcont[k][j][i].x
                    =
                    (
                        gravity.x * icsi[k][j][i].x +
                        gravity.y * icsi[k][j][i].y +
                        gravity.z * icsi[k][j][i].z
                    );
                }

                if(i > 0 && k > 0 && i < mx-1 && k < mz-1)
                {
                    gcont[k][j][i].y
                    =
                    (
                        gravity.x * jeta[k][j][i].x +
                        gravity.y * jeta[k][j][i].y +
                        gravity.z * jeta[k][j][i].z
                    );
                }

                if(i > 0 && j > 0 && j < my-1 && i < mx-1)
                {
                    gcont[k][j][i].z
                    =
                    (
                        gravity.x * kzet[k][j][i].x +
                        gravity.y * kzet[k][j][i].y +
                        gravity.z * kzet[k][j][i].z
                    );
                }

            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

    VecSet(sm->sedCent,  0.0);
    DMDAVecGetArray(da, sm->sedCent, &sedCent);

    // set sed value at cell centers
    for(k=lzs; k<lze; k++)
    {
        for(j=lys; j<lye; j++)
        {
            for(i=lxs; i<lxe; i++)
            {

                //skip edges
                /*if(i==0 && j==0) continue;
                if(j==0 && k==0) continue;
                if(i==0 && k==0) continue;
                if(i==0 && j==my-1) continue;
                if(j==0 && k==mz-1) continue;
                if(i==0 && k==mz-1) continue;
                if(i==mx-1 && j==my-1) continue;
                if(j==my-1 && k==mz-1) continue;
                if(i==mx-1 && k==mz-1) continue;
                if(i==mx-1 && j==0) continue;
                if(j==my-1 && k==0) continue;
                if(i==mx-1 && k==0) continue;*/

                //no need for SM value at solid cells. Will be tracked later on.
                if (isIBMSolidCell(k, j, i, nvert)) continue;

                if (weight0[k][j][i] <=0.0 || weight1[k][j][i] <= 0.0 || weight2[k][j][i] <= 0.0 || absc0[k][j][i] <= 0.0 || absc1[k][j][i] <= 0.0 || absc2[k][j][i] <= 0.0)
                {
                    continue;
                }

                if (sm->access->flags->isTeqnActive)
                {
                    DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                    temp = tmprt[k][j][i]; // need to get at face value??

                    DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                }
                else
                {
                    temp = Cst->tRef; //need at face??
                }

                mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                //m3Dia = sm4[k][j][i]/sm3[k][j][i];

                expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                TauP0[k][j][i] = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*Cst->nu*Cst->rho);
                TauP1[k][j][i] = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*Cst->nu*Cst->rho);
                TauP2[k][j][i] = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*Cst->nu*Cst->rho);

                sedCent[k][j][i]
                =
                (TauP0[k][j][i]*weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)  +  TauP1[k][j][i]*weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii) +  TauP2[k][j][i]*weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii));


            }
        }
    }

    DMDAVecRestoreArray(da, sm->sedCent, &sedCent);

    DMGlobalToLocalBegin(da, sm->sedCent,  INSERT_VALUES, sm->lSedCent);
    DMGlobalToLocalEnd  (da, sm->sedCent,  INSERT_VALUES, sm->lSedCent);

    DMDAVecGetArray(da, sm->lSedCent, &lSedCent);
    VecSet(sm->Sed, 0.0);
    DMDAVecGetArray(fda, sm->Sed, &sed);

    //set sed flux at faces, skip edges since dep takes sed into account at wall.
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                /*//skip left edges
                if(i==0 && j==0) continue;
                if(j==0 && k==0) continue;
                if(i==0 && k==0) continue;

                //skip right faces
                if(i==mx-1 || j==my-1 || k==mz-1) continue;*/

                //i-face Sed.y
                //skip all fluid cells with solid on right face
                if (isIBMSolidCell(k, j, i+1, nvert) || isIBMSolidCell(k, j, i, nvert))
                {
                    sed[k][j][i].y = 0.;
                }
                //no sedimentation from upper walls where vents are not located
                else if (!mesh->ii_periodic && (i == mx-2) && markVent[k][j][i] == 0)
                {
                    sed[k][j][i].y = 0.;
                }
                else
                {
                    sed[k][j][i].y
                    =
                    (
                        gcont[k][j][i].x * lSedCent[k][j][i+1]
                    );
                }

                //j-face Sed.z
                //skip all solid cells and fluid cells with solid on right face
                if (isIBMSolidCell(k, j+1, i, nvert) || isIBMSolidCell(k, j, i, nvert))
                {
                    sed[k][j][i].z = 0.;
                }
                //no sedimentation through upper walls where vents are not located
                else if (!mesh->jj_periodic && (j == my-2) && markVent[k][j][i] == 0)
                {
                    sed[k][j][i].z = 0.;
                }
                else
                {
                    sed[k][j][i].z
                    =
                    (
                        gcont[k][j][i].y * lSedCent[k][j+1][i]
                    );

                }

                //k-face Sed.x
                //skip all solid cells and fluid cells with solid on right face
                if (isIBMSolidCell(k+1, j, i, nvert) || isIBMSolidCell(k, j, i, nvert))
                {
                    sed[k][j][i].x = 0.;
                }
                //no sedimentation through upper walls where vents are not located
                else if (!mesh->kk_periodic && (k == mz-2) && markVent[k][j][i] == 0)
                {
                    sed[k][j][i].x = 0.;
                }
                else
                {
                    sed[k][j][i].x
                    =
                    (
                        gcont[k][j][i].z * lSedCent[k+1][j][i]
                    );
                }

            }
        }
    }


    DMDAVecRestoreArray(da, sm->lSedCent, &lSedCent);
    DMDAVecRestoreArray(fda, sm->Sed, &sed);

    DMGlobalToLocalBegin(fda, sm->Sed,  INSERT_VALUES, sm->lSed);
    DMGlobalToLocalEnd  (fda, sm->Sed,  INSERT_VALUES, sm->lSed);

    resetFacePeriodicFluxesVector(mesh, sm->Sed,   sm->lSed, "globalToLocal");

    DMDAVecGetArray(fda, sm->lSed, &lsed);
    DMDAVecGetArray(da, Rhs, &rhs);

    // ---------------------------------------------------------------------- //
    //                      FORM THE CUMULATIVE FLUXES                        //
    // ---------------------------------------------------------------------- //

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if (isIBMSolidCell(k, j, i, nvert))
                {
                    rhs[k][j][i] = 0;
                }
                else
                {
                    rhs[k][j][i]
                    +=
                    scale *
                    (
                         lsed[k][j][i].y  - lsed[k  ][j  ][i-1].y       +
                         lsed[k][j][i].z  - lsed[k  ][j-1][i  ].z       +
                         lsed[k][j][i].x  - lsed[k-1][j  ][i  ].x
                    ) * aj[k][j][i];

                }

            }
        }
    }

    DMDAVecRestoreArray(fda, sm->lSed, &lsed);
    DMDAVecRestoreArray(da, Rhs, &rhs);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
    DMDAVecRestoreArray(fda, sm->access->ueqn->gCont,  &gcont);

    DMDAVecRestoreArray(da, sm->smVal, &sm0);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);
    DMDAVecRestoreArray(fda, sm->access->ueqn->lUcont,      &ucont);

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode devFluxSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii)
{
    mesh_         *mesh  = sm->access->mesh;
    SMObj_        *smObject = sm->access->smObject;
    constants_    *Cst   = sm->access->constants;
    clock_         *clock = sm->access->clock;

    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***gcont, ***dev, ***ldev, ***limiter;
    Cmpnts        ***cent, ***ucont;
    PetscReal     ***nvert, ***tmprt, ***rhs, ***devCent, ***lDevCent;
    PetscReal     ***weight0, ***absc0, ***weight1, ***absc1, ***weight2, ***absc2, ***TauP0, ***TauP1, ***TauP2, ***sm0;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     temp;
    PetscReal     mfp;
    PetscReal     m3Dia;
    PetscReal     expMfp0, expMfp1, expMfp2;
    PetscReal     Cu0, Cu1, Cu2;

    PetscInt    iL, iR, jL, jR, kL, kR;
    PetscReal   denom;

    PetscReal     du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

    PetscReal     dxdc, dxde, dxdz, dydc, dyde, dydz, dzdc, dzde, dzdz;
    PetscReal     dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords

    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2, ajc;      // surface area vectors components

    Cmpnts	      ***icsi, ***ieta, ***izet;
    Cmpnts	      ***jcsi, ***jeta, ***jzet;
    Cmpnts	      ***kcsi, ***keta, ***kzet;

    PetscReal     ***iaj, ***jaj, ***kaj, ***aj;

    PetscInt      ***markVent;

    PetscReal     k_b = 1.380649E-23; //*1E18; //boltzman constant in kg*nm^2/(K*s^2)

    Cmpnts        ***ucat, ***ucat_o, matDer, matDerCont;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->lAj, &aj);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);

    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);

    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(fda, sm->access->ueqn->gCont,  &gcont);
    DMDAVecGetArray(fda, sm->access->ueqn->lUcont,      &ucont);
    DMDAVecGetArray(fda, sm->access->ueqn->lUcat,      &ucat);
    DMDAVecGetArray(fda, sm->access->ueqn->Ucat_o,      &ucat_o);

    DMDAVecGetArray(da, sm->smVal, &sm0);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecGetArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecGetArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecGetArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecGetArray(da, smObject->weightAbsc[2]->tauP, &TauP2);

    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    VecSet(sm->devCent, 0.0);
    DMDAVecGetArray(da, sm->devCent, &devCent);

    // set dev value at cell centers
    for(k=lzs; k<lze; k++)
    {
        for(j=lys; j<lye; j++)
        {
            for(i=lxs; i<lxe; i++)
            {

                if(isIBMSolidCell(k, j, i, nvert)) continue;

                if (weight0[k][j][i] <= 0.0 || weight1[k][j][i] <= 0.0 || weight2[k][j][i] <= 0.0 || absc0[k][j][i] <= 0.0 || absc1[k][j][i] <= 0.0 || absc2[k][j][i] <= 0.0)
                {
                    continue;
                }

                if (sm->access->flags->isTeqnActive)
                {
                    DMDAVecGetArray(da, sm->access->teqn->lTmprt, &tmprt);

                    temp = tmprt[k][j][i]; // need to get at face value??

                    DMDAVecRestoreArray(da, sm->access->teqn->lTmprt, &tmprt);
                }
                else
                {
                    temp = Cst->tRef; //need at face??
                }

                mfp = pow(10, 6)*((k_b)*temp/(sqrt(2)*M_PI*3.46E-10*3.46E-10*101325)); // mfp in um

                //m3Dia = sm4[k][j][i]/sm3[k][j][i];

                expMfp0 = exp(-0.88/(2*mfp/(absc0[k][j][i])));
                expMfp1 = exp(-0.88/(2*mfp/(absc1[k][j][i])));
                expMfp2 = exp(-0.88/(2*mfp/(absc2[k][j][i])));

                Cu0 = (1 + 2*mfp*(1.2+0.4*expMfp0)/(absc0[k][j][i]));
                Cu1 = (1 + 2*mfp*(1.2+0.4*expMfp1)/(absc1[k][j][i]));
                Cu2 = (1 + 2*mfp*(1.2+0.4*expMfp2)/(absc2[k][j][i]));

                TauP0[k][j][i] = pow(10, -12)*absc0[k][j][i]*absc0[k][j][i]*smObject->rhoPart*Cu0/(18*Cst->nu*Cst->rho);
                TauP1[k][j][i] = pow(10, -12)*absc1[k][j][i]*absc1[k][j][i]*smObject->rhoPart*Cu1/(18*Cst->nu*Cst->rho);
                TauP2[k][j][i] = pow(10, -12)*absc2[k][j][i]*absc2[k][j][i]*smObject->rhoPart*Cu2/(18*Cst->nu*Cst->rho);

                devCent[k][j][i]
                =
                (TauP0[k][j][i]*weight0[k][j][i]*PetscPowScalar(absc0[k][j][i], ii)  +  TauP1[k][j][i]*weight1[k][j][i]*PetscPowScalar(absc1[k][j][i], ii) +  TauP2[k][j][i]*weight2[k][j][i]*PetscPowScalar(absc2[k][j][i], ii));


            }
        }
    }

    DMDAVecRestoreArray(da, sm->devCent, &devCent);

    DMGlobalToLocalBegin(da, sm->devCent,  INSERT_VALUES, sm->lDevCent);
    DMGlobalToLocalEnd  (da, sm->devCent,  INSERT_VALUES, sm->lDevCent);

    DMDAVecGetArray(da, sm->lDevCent, &lDevCent);
    VecSet(sm->Dev, 0.0);
    DMDAVecGetArray(fda, sm->Dev, &dev);

    //set dev flux at faces, but no dev through boundaries
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if (isIBMSolidCell(k, j, i, nvert))
                {
                    continue;
                }

                //i-faces dev.y
                csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;
                ajc  = iaj[k][j][i];

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_i
                (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                    &dudc, &dvdc, &dwdc,
                    &dude, &dvde, &dwde,
                    &dudz, &dvdz, &dwdz
                );

                // compute cartesian velocity derivatives w.r.t cartesian coords
                Compute_du_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                    zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                    dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                    &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                //skip if at left IBFluid cell (IB solid face)
                if (isIBMSolidCell(k, j, i+1, nvert) || (i == mx-2 && !mesh->ii_periodic))
                {
                    dev[k][j][i].y = 0.0;
                }
                else
                {
                    matDer.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDer.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDer.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    matDerCont.x
                    =
                    (
                        matDer.x * icsi[k][j][i].x +
                        matDer.y * icsi[k][j][i].y +
                        matDer.z * icsi[k][j][i].z
                    );

                    dev[k][j][i].y
                    =
                    -
                    (
                        matDerCont.x * upwind(lDevCent[k][j][i], lDevCent[k][j][i+1], k, j, i, ucont[k][j][i].x)
                    );
                }

                //j-face dev.z
                csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;
                ajc  = jaj[k][j][i];

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_j
                (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                    &dudc, &dvdc, &dwdc,
                    &dude, &dvde, &dwde,
                    &dudz, &dvdz, &dwdz
                );

                // compute cartesian velocity derivatives w.r.t cartesian coords
                Compute_du_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                    zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                    dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                    &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                //skip if at left IBFluid cell (IB solid face)
                if (isIBMSolidCell(k, j+1, i, nvert) || (j == my-2 && !mesh->jj_periodic))
                {
                    continue;
                }
                else
                {
                    matDer.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDer.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDer.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    matDerCont.y
                    =
                    (
                        matDer.x * jeta[k][j][i].x +
                        matDer.y * jeta[k][j][i].y +
                        matDer.z * jeta[k][j][i].z
                    );

                    dev[k][j][i].z
                    =
                    -
                    (
                        matDerCont.y * upwind(lDevCent[k][j][i], lDevCent[k][j+1][i], k, j, i, ucont[k][j][i].y)
                    );
                }

                //k-face sed.x
                csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;
                ajc  = kaj[k][j][i];

                // compute cartesian velocity derivatives w.r.t. curvilinear coords
                Compute_du_k
                (   mesh, i, j, k, mx, my, mz, ucat, nvert,
                    &dudc, &dvdc, &dwdc,
                    &dude, &dvde, &dwde,
                    &dudz, &dvdz, &dwdz
                );

                // compute cartesian velocity derivatives w.r.t cartesian coords
                Compute_du_dxyz
                (
                    mesh,
                    csi0, csi1, csi2, eta0, eta1, eta2, zet0,
                    zet1, zet2, ajc, dudc, dvdc, dwdc, dude, dvde, dwde,
                    dudz, dvdz, dwdz, &du_dx, &dv_dx, &dw_dx, &du_dy,
                    &dv_dy, &dw_dy, &du_dz, &dv_dz, &dw_dz
                );

                //skip if at left IBFluid cell (IB solid face)
                if (isIBMSolidCell(k+1, j, i, nvert) || (k == mz-2 && !mesh->kk_periodic))
                {
                    continue;
                }
                else
                {
                    matDer.x = (ucat[k][j][i].x - ucat_o[k][j][i].x)/clock->dt + ucat[k][j][i].x * du_dx + ucat[k][j][i].y * du_dy + ucat[k][j][i].z * du_dz;
                    matDer.y = (ucat[k][j][i].y - ucat_o[k][j][i].y)/clock->dt + ucat[k][j][i].x * dv_dx + ucat[k][j][i].y * dv_dy + ucat[k][j][i].z * dv_dz;
                    matDer.z = (ucat[k][j][i].z - ucat_o[k][j][i].z)/clock->dt + ucat[k][j][i].x * dw_dx + ucat[k][j][i].y * dw_dy + ucat[k][j][i].z * dw_dz;

                    matDerCont.z
                    =
                    (
                        matDer.x * kzet[k][j][i].x +
                        matDer.y * kzet[k][j][i].y +
                        matDer.z * kzet[k][j][i].z
                    );

                    dev[k][j][i].x
                    =
                    -
                    (
                        matDerCont.z * upwind(lDevCent[k][j][i], lDevCent[k+1][j][i], k, j, i, ucont[k][j][i].z)
                    );
                }

            }
        }
    }

    DMDAVecRestoreArray(da, sm->lDevCent, &lDevCent);
    DMDAVecRestoreArray(fda, sm->Dev, &dev);

    DMGlobalToLocalBegin(fda, sm->Dev,  INSERT_VALUES, sm->lDev);
    DMGlobalToLocalEnd  (fda, sm->Dev,  INSERT_VALUES, sm->lDev);

    resetFacePeriodicFluxesVector(mesh, sm->Dev,   sm->lDev, "globalToLocal");

    DMDAVecGetArray(fda, sm->lDev, &ldev);
    DMDAVecGetArray(da, Rhs, &rhs);

    // ---------------------------------------------------------------------- //
    //                      FORM THE CUMULATIVE FLUXES                        //
    // ---------------------------------------------------------------------- //

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if (isIBMSolidCell(k, j, i, nvert))
                {
                    rhs[k][j][i] = 0;
                }
                else
                {
                    rhs[k][j][i]
                    +=
                    scale *
                    (
                         ldev[k][j][i].y  - ldev[k  ][j  ][i-1].y       +
                         ldev[k][j][i].z  - ldev[k  ][j-1][i  ].z       +
                         ldev[k][j][i].x  - ldev[k-1][j  ][i  ].x
                    ) * aj[k][j][i];

                }

            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lAj, &aj);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, sm->lDev, &ldev);
    DMDAVecRestoreArray(da, Rhs, &rhs);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
    DMDAVecRestoreArray(fda, sm->access->ueqn->gCont,  &gcont);

    DMDAVecRestoreArray(da, sm->smVal, &sm0);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->weight, &weight0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->absc, &absc0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->weight, &weight1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->absc, &absc1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->weight, &weight2);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->absc, &absc2);

    DMDAVecRestoreArray(da, smObject->weightAbsc[0]->tauP, &TauP0);
    DMDAVecRestoreArray(da, smObject->weightAbsc[1]->tauP, &TauP1);
    DMDAVecRestoreArray(da, smObject->weightAbsc[2]->tauP, &TauP2);


    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);
    DMDAVecRestoreArray(fda, sm->access->ueqn->lUcont,      &ucont);
    DMDAVecRestoreArray(fda, sm->access->ueqn->lUcat,      &ucat);
    DMDAVecRestoreArray(fda, sm->access->ueqn->Ucat_o,      &ucat_o);

    DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    return(0);
}
