//! \file  aeqn_solvers.c
//! \brief Alpha-Water-equation solver function implementations.

#include "aeqn_solvers.h"

PetscErrorCode FormALowOrder(aeqn_ *aeqn)
{
    // compute first-order upwind face fluxes and store in lDivA

    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    PetscReal     ***alpha;
    Cmpnts        ***ucont, ***div;

    DMDAVecGetArray(da,  aeqn->lAlpha,  &alpha);
    DMDAVecGetArray(fda, ueqn->lUcont,  &ucont);

    VecSet(aeqn->lDivALO, 0.0);
    DMDAVecGetArray(fda, aeqn->lDivALO, &div);

    // ---------------------------------------------------------------------- //
    //                  upwind face flux assembly (i, j, k)                   //
    // ---------------------------------------------------------------------- //

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                if(j!=0 && k!=0)
                    div[k][j][i].x = -ucont[k][j][i].x * upwind(alpha[k][j][i], alpha[k][j][i+1], ucont[k][j][i].x);

                if(i!=0 && k!=0)
                    div[k][j][i].y = -ucont[k][j][i].y * upwind(alpha[k][j][i], alpha[k][j+1][i], ucont[k][j][i].y);

                if(i!=0 && j!=0)
                    div[k][j][i].z = -ucont[k][j][i].z * upwind(alpha[k][j][i], alpha[k+1][j][i], ucont[k][j][i].z);
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivALO, &div);
    DMDAVecRestoreArray(da,  aeqn->lAlpha, &alpha);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);

    DMLocalToLocalBegin(fda, aeqn->lDivALO, INSERT_VALUES, aeqn->lDivALO);
    DMLocalToLocalEnd  (fda, aeqn->lDivALO, INSERT_VALUES, aeqn->lDivALO);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivALO, aeqn->lDivALO, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormAHighOrder(aeqn_ *aeqn)
{
    // compute fourth-order central face fluxes and store in lDivAHO

    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    PetscReal     ***alpha, ***nvert, ***meshTag;
    Cmpnts        ***ucont, ***div;

    DMDAVecGetArray(da,  aeqn->lAlpha,    &alpha);
    DMDAVecGetArray(fda, ueqn->lUcont,    &ucont);
    DMDAVecGetArray(da,  mesh->lNvert,    &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,  &meshTag);

    VecSet(aeqn->lDivAHO, 0.0);
    DMDAVecGetArray(fda, aeqn->lDivAHO, &div);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                // i-face
                if(j!=0 && k!=0)
                {
                    if(isOversetIFace(k, j, i, i+1, meshTag))
                    {
                        div[k][j][i].x = -ucont[k][j][i].x * central(alpha[k][j][i], alpha[k][j][i+1]);
                    }
                    else
                    {
                        PetscInt  iL, iR;
                        PetscReal denom;
                        getFace2Cell4StencilCsi(mesh, k, j, i, mx, &iL, &iR, &denom, nvert, meshTag);
                        div[k][j][i].x = -ucont[k][j][i].x * central4(alpha[k][j][iL], alpha[k][j][i], alpha[k][j][i+1], alpha[k][j][iR]);
                        //div[k][j][i].x = -ucont[k][j][i].x * central(alpha[k][j][i], alpha[k][j][i+1]);
                    }
                }

                // j-face
                if(i!=0 && k!=0)
                {
                    if(isOversetJFace(k, j, i, j+1, meshTag))
                    {
                        div[k][j][i].y = -ucont[k][j][i].y * central(alpha[k][j][i], alpha[k][j+1][i]);
                    }
                    else
                    {
                        PetscInt  jL, jR;
                        PetscReal denom;
                        getFace2Cell4StencilEta(mesh, k, j, i, my, &jL, &jR, &denom, nvert, meshTag);
                        div[k][j][i].y = -ucont[k][j][i].y * central4(alpha[k][jL][i], alpha[k][j][i], alpha[k][j+1][i], alpha[k][jR][i]);
                        //div[k][j][i].y = -ucont[k][j][i].y * central(alpha[k][j][i], alpha[k][j+1][i]);
                    }
                }

                // k-face
                if(i!=0 && j!=0)
                {
                    if(isOversetKFace(k, j, i, k+1, meshTag))
                    {
                        div[k][j][i].z = -ucont[k][j][i].z * central(alpha[k][j][i], alpha[k+1][j][i]);
                    }
                    else
                    {
                        PetscInt  kL, kR;
                        PetscReal denom;
                        getFace2Cell4StencilZet(mesh, k, j, i, mz, &kL, &kR, &denom, nvert, meshTag);
                        div[k][j][i].z = -ucont[k][j][i].z * central4(alpha[kL][j][i], alpha[k][j][i], alpha[k+1][j][i], alpha[kR][j][i]);
                        //div[k][j][i].z = -ucont[k][j][i].z * central(alpha[k][j][i], alpha[k+1][j][i]);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivAHO,   &div);
    DMDAVecRestoreArray(da,  aeqn->lAlpha,    &alpha);
    DMDAVecRestoreArray(fda, ueqn->lUcont,    &ucont);
    DMDAVecRestoreArray(da,  mesh->lNvert,    &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,  &meshTag);

    DMLocalToLocalBegin(fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);
    DMLocalToLocalEnd  (fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivAHO, aeqn->lDivAHO, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ComputeLowOrderUpdate(aeqn_ *aeqn)
{
    // provisional low-order cell update: alpha_LO = alpha + dt * aj * divBudget_LO
    // divBudget is net inflow (positive) because FormALowOrder stores -u*alpha at each face
    // clips to [0,1] and syncs ghosts for use in the mules limiter R+/R- bounds

    mesh_         *mesh  = aeqn->access->mesh;
    clock_        *clock = aeqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***alpha, ***alphaLO, ***aj;
    Cmpnts        ***div;

    PetscReal     dt = clock->dt;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  aeqn->lAlpha,   &alpha);
    DMDAVecGetArray(da,  aeqn->lAlphaLO, &alphaLO);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(fda, aeqn->lDivALO,  &div);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal divBudget =
                    div[k][j][i].x - div[k][j][i-1].x +
                    div[k][j][i].y - div[k][j-1][i].y +
                    div[k][j][i].z - div[k-1][j][i].z;

                // provisional low-order update: alpha + dt * inflow_rate
                alphaLO[k][j][i] = PetscMax(0.0, PetscMin(1.0, alpha[k][j][i] + dt * aj[k][j][i] * divBudget));
            }
        }
    }

    DMDAVecRestoreArray(da,  aeqn->lAlpha,   &alpha);
    DMDAVecRestoreArray(da,  aeqn->lAlphaLO, &alphaLO);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivALO,  &div);

    // sync ghost cells so neighboring processors can read alpha_LO
    DMLocalToLocalBegin(da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);
    DMLocalToLocalEnd  (da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AddCompressionFlux(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    Vec           GradA, AlphaSmooth;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     dadcsi, dadeta, dadzet;                                        // alpha der. w.r.t. curvil. coords, 
    PetscReal     dadx, dady, dadz;                                              // alpha der. w.r.t. cart. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components
    
    PetscReal     ***aj  , ***iaj , ***jaj , ***kaj;
    Cmpnts        ***csi , ***eta , ***zet ;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;

    PetscReal     ***alpha, ***alpha_smooth, ***nvert, ***meshTag;
    Cmpnts        ***ucont, ***ucatc, ***div, ***gradAlpha;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // get numerical compression coefficient
    PetscReal compcoeff = aeqn->compCoeff; 

    VecDuplicate(mesh->lCent,  &GradA);
    VecDuplicate(aeqn->lAlpha, &AlphaSmooth);
    VecSet(GradA, 0.0);
    VecSet(AlphaSmooth, 0.0);

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

    DMDAVecGetArray(fda, ueqn->lUcont,     &ucont);
    DMDAVecGetArray(da,  aeqn->lAlpha,     &alpha);
    DMDAVecGetArray(da,  mesh->lNvert,     &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,   &meshTag);
    DMDAVecGetArray(fda, aeqn->lDivAHO,    &div);

    DMDAVecGetArray(da,  AlphaSmooth,      &alpha_smooth);

    // step 0: smooth alpha by averaging over a 3x3x3 stencil
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal sum = 0.0;
                for (PetscInt kk=-1; kk<=1; kk++)
                {
                    for (PetscInt jj=-1; jj<=1; jj++)
                    {
                        for (PetscInt ii=-1; ii<=1; ii++)
                        {
                            sum += alpha[k+kk][j+jj][i+ii];
                        }
                    }
                }
                alpha_smooth[k][j][i] = sum / 27.0;
            }
        }
    }

    DMDAVecRestoreArray(da,  AlphaSmooth,      &alpha_smooth);

    // sync ghost cells so neighboring processors can read GradA
    DMLocalToLocalBegin(da, AlphaSmooth, INSERT_VALUES, AlphaSmooth);
    DMLocalToLocalEnd  (da, AlphaSmooth, INSERT_VALUES, AlphaSmooth);


    DMDAVecGetArray(fda, GradA,       &gradAlpha);
    DMDAVecGetArray(da,  AlphaSmooth, &alpha_smooth);

    // step 1: compute alpha gradient at cell centers 
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // get 1/V at the i-face
                PetscReal ajc = aj[k][j][i];

                // get face normals
                csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
                eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
                zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

                // compute metric tensor
                g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                Compute_dscalar_center
                (
                    mesh,
                    i, j, k,  mx, my, mz,
                    alpha_smooth, nvert, meshTag, 
                    &dadcsi, &dadeta, &dadzet
                );

                Compute_dscalar_dxyz
                (
                    mesh,
                    csi0, csi1, csi2,
                    eta0, eta1, eta2,
                    zet0, zet1, zet2,
                    ajc,
                    dadcsi, dadeta, dadzet,
                    &dadx, &dady, &dadz
                );

                gradAlpha[k][j][i].x = dadx;
                gradAlpha[k][j][i].y = dady;
                gradAlpha[k][j][i].z = dadz;
            }
        }
    }

    DMDAVecRestoreArray(fda, GradA, &gradAlpha);
    DMDAVecRestoreArray(da,  AlphaSmooth, &alpha_smooth);

    // sync ghost cells so neighboring processors can read GradA
    DMLocalToLocalBegin(fda, GradA, INSERT_VALUES, GradA);
    DMLocalToLocalEnd  (fda, GradA, INSERT_VALUES, GradA);

    // reset periodic fluxes in GradA so that they are correct at periodic boundaries
    resetCellPeriodicFluxes(mesh, GradA, GradA, "vector", "localToLocal");

    DMDAVecGetArray(fda, GradA, &gradAlpha);

    // compute interface velocity at cell faces as uc_face = c * mag_u_face * grad_alpha_face / mag_grad_alpha_face
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                // i-faces
                if(j!=0 && k!=0)
                {
                    PetscReal alpha_face;
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j][i+1]);

                    Cmpnts gradAlpha_face;
                    gradAlpha_face.x = 0.5 * (gradAlpha[k][j][i].x + gradAlpha[k][j][i+1].x);
                    gradAlpha_face.y = 0.5 * (gradAlpha[k][j][i].y + gradAlpha[k][j][i+1].y);
                    gradAlpha_face.z = 0.5 * (gradAlpha[k][j][i].z + gradAlpha[k][j][i+1].z);

                    PetscReal ucat_face_mag      = fabs(ucont[k][j][i].x) / nMag(icsi[k][j][i]); 
                    PetscReal gradAlpha_face_mag = nMag(gradAlpha_face);

                    Cmpnts ucatc_face;
                    ucatc_face.x = compcoeff * ucat_face_mag * gradAlpha_face.x / (gradAlpha_face_mag + 1e-10);
                    ucatc_face.y = compcoeff * ucat_face_mag * gradAlpha_face.y / (gradAlpha_face_mag + 1e-10);
                    ucatc_face.z = compcoeff * ucat_face_mag * gradAlpha_face.z / (gradAlpha_face_mag + 1e-10);

                    PetscReal coeff = alpha_face*(1.0 - alpha_face);
                    div[k][j][i].x -= coeff * nDot(ucatc_face, icsi[k][j][i]);
                }

                // j-faces
                if(i!=0 && k!=0)
                {
                    PetscReal alpha_face;
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j+1][i]);

                    Cmpnts gradAlpha_face;
                    gradAlpha_face.x = 0.5 * (gradAlpha[k][j][i].x + gradAlpha[k][j+1][i].x);
                    gradAlpha_face.y = 0.5 * (gradAlpha[k][j][i].y + gradAlpha[k][j+1][i].y);
                    gradAlpha_face.z = 0.5 * (gradAlpha[k][j][i].z + gradAlpha[k][j+1][i].z);

                    PetscReal ucat_face_mag      = fabs(ucont[k][j][i].y) / nMag(jeta[k][j][i]);
                    PetscReal gradAlpha_face_mag = nMag(gradAlpha_face);

                    Cmpnts ucatc_face;
                    ucatc_face.x = compcoeff * ucat_face_mag * gradAlpha_face.x / (gradAlpha_face_mag + 1e-10);
                    ucatc_face.y = compcoeff * ucat_face_mag * gradAlpha_face.y / (gradAlpha_face_mag + 1e-10);
                    ucatc_face.z = compcoeff * ucat_face_mag * gradAlpha_face.z / (gradAlpha_face_mag + 1e-10);

                    PetscReal coeff = alpha_face*(1.0 - alpha_face);
                    div[k][j][i].y -= coeff * nDot(ucatc_face, jeta[k][j][i]);
                }

                // k-faces
                if(i!=0 && j!=0)
                {
                    PetscReal alpha_face;
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k+1][j][i]);

                    Cmpnts gradAlpha_face;
                    gradAlpha_face.x = 0.5 * (gradAlpha[k][j][i].x + gradAlpha[k+1][j][i].x);
                    gradAlpha_face.y = 0.5 * (gradAlpha[k][j][i].y + gradAlpha[k+1][j][i].y);
                    gradAlpha_face.z = 0.5 * (gradAlpha[k][j][i].z + gradAlpha[k+1][j][i].z);

                    PetscReal ucat_face_mag      = fabs(ucont[k][j][i].z) / nMag(kzet[k][j][i]);
                    PetscReal gradAlpha_face_mag = nMag(gradAlpha_face);

                    Cmpnts ucatc_face;
                    ucatc_face.x = compcoeff * ucat_face_mag * gradAlpha_face.x / (gradAlpha_face_mag + 1e-15);
                    ucatc_face.y = compcoeff * ucat_face_mag * gradAlpha_face.y / (gradAlpha_face_mag + 1e-15);
                    ucatc_face.z = compcoeff * ucat_face_mag * gradAlpha_face.z / (gradAlpha_face_mag + 1e-15);

                    PetscReal coeff = alpha_face*(1.0 - alpha_face);
                    div[k][j][i].z -= coeff * nDot(ucatc_face, kzet[k][j][i]);
                }
            }
        }
    }

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

    DMDAVecRestoreArray(fda, ueqn->lUcont,     &ucont);
    DMDAVecRestoreArray(da,  aeqn->lAlpha,     &alpha);
    DMDAVecRestoreArray(da,  mesh->lNvert,     &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,   &meshTag);
    DMDAVecRestoreArray(fda, aeqn->lDivAHO,    &div);

    DMDAVecRestoreArray(fda, GradA, &gradAlpha);

    DMLocalToLocalBegin(fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);
    DMLocalToLocalEnd  (fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);

    VecDestroy(&GradA);
    VecDestroy(&AlphaSmooth);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivAHO, aeqn->lDivAHO, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ApplyMULESLimiter(aeqn_ *aeqn)
{
    // apply mules flux limiter to produce a bounded, conservative Rhs
    // step 1: correction = lDivAHO - lDivA
    // step 2: provisional low-order cell update into lAlphaLO
    // step 3: per-cell R_plus, R_minus (how much correction each cell can absorb)
    // step 4: per-face lambda_f = min(R of owner, R of neighbor)
    // step 5: final Rhs = low-order budget + lambda-limited correction budget

    mesh_         *mesh  = aeqn->access->mesh;
    clock_        *clock = aeqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***alphaLO, ***aj;
    PetscReal     ***Rplus, ***Rminus;
    Cmpnts        ***div, ***divLO, ***divCor, ***lambda;

    PetscReal     dt  = clock->dt;
    PetscReal     eps = 1.0e-12;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // correction flux lDivACor = lDivAHO - lDivALO               //
    VecCopy (aeqn->lDivAHO,  aeqn->lDivACor);
    VecAXPY (aeqn->lDivACor, -1.0, aeqn->lDivALO);

    // initialize to 1.0: physical boundary ghost cells allow full correction
    VecSet(aeqn->lRplusA,  1.0);
    VecSet(aeqn->lRminusA, 1.0);

    // alphaLO must be up to date
    DMDAVecGetArray(da,  aeqn->lAlphaLO,  &alphaLO);
    DMDAVecGetArray(da,  aeqn->lRplusA,   &Rplus);
    DMDAVecGetArray(da,  aeqn->lRminusA,  &Rminus);
    DMDAVecGetArray(da,  mesh->lAj,       &aj);
    DMDAVecGetArray(fda, aeqn->lDivACor,  &divCor);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // room for alpha to increase or decrease within [0,1]
                PetscReal P_plus  = (1.0 / (aj[k][j][i] * dt)) * (1.0 - alphaLO[k][j][i]);
                PetscReal P_minus = (1.0 / (aj[k][j][i] * dt)) *  alphaLO[k][j][i];

                // sum correction contributions that tend to increase (C_plus) or decrease (C_minus) alpha
                // right faces enter divBudget with + sign; left faces with - sign
                PetscReal C_plus = 0.0, C_minus = 0.0;
                PetscReal c;

                c =  divCor[k][j][i].x;    // right i-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -divCor[k][j][i-1].x;  // left i-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c =  divCor[k][j][i].y;    // right j-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -divCor[k][j-1][i].y;  // left j-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c =  divCor[k][j][i].z;    // right k-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -divCor[k-1][j][i].z;  // left k-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                Rplus [k][j][i] = (C_plus  > eps) ? PetscMin(1.0, P_plus  / C_plus)  : 1.0;
                Rminus[k][j][i] = (C_minus > eps) ? PetscMin(1.0, P_minus / C_minus) : 1.0;
            }
        }
    }

    DMDAVecRestoreArray(da,  aeqn->lAlphaLO,  &alphaLO);
    DMDAVecRestoreArray(da,  aeqn->lRplusA,   &Rplus);
    DMDAVecRestoreArray(da,  aeqn->lRminusA,  &Rminus);
    DMDAVecRestoreArray(da,  mesh->lAj,       &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivACor,  &divCor);

    // sync R values so each processor can read neighbor R for face lambda
    DMLocalToLocalBegin(da, aeqn->lRplusA,  INSERT_VALUES, aeqn->lRplusA);
    DMLocalToLocalEnd  (da, aeqn->lRplusA,  INSERT_VALUES, aeqn->lRplusA);
    DMLocalToLocalBegin(da, aeqn->lRminusA, INSERT_VALUES, aeqn->lRminusA);
    DMLocalToLocalEnd  (da, aeqn->lRminusA, INSERT_VALUES, aeqn->lRminusA);

    // compute per-face limiter coefficient lambda_f in [0,1]
    VecSet(aeqn->lLambdaA, 1.0);
    DMDAVecGetArray(fda, aeqn->lDivACor,  &divCor);
    DMDAVecGetArray(fda, aeqn->lLambdaA,  &lambda);
    DMDAVecGetArray(da,  aeqn->lRplusA,   &Rplus);
    DMDAVecGetArray(da,  aeqn->lRminusA,  &Rminus);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                // i-faces: owner=(k,j,i), neighbor=(k,j,i+1)
                if(j!=0 && k!=0)
                {
                    PetscReal dF = divCor[k][j][i].x;
                    if(dF < 0.0)
                        lambda[k][j][i].x = PetscMin(Rminus[k][j][i], Rplus [k][j][i+1]);
                    else
                        lambda[k][j][i].x = PetscMin(Rplus [k][j][i], Rminus[k][j][i+1]);
                }

                // j-face: owner=(k,j,i), neighbor=(k,j+1,i)
                if(i!=0 && k!=0)
                {
                    PetscReal dF = divCor[k][j][i].y;
                    if(dF < 0.0)
                        lambda[k][j][i].y = PetscMin(Rminus[k][j][i], Rplus [k][j+1][i]);
                    else
                        lambda[k][j][i].y = PetscMin(Rplus [k][j][i], Rminus[k][j+1][i]);
                }

                // k-face: owner=(k,j,i), neighbor=(k+1,j,i)
                if(i!=0 && j!=0)
                {
                    PetscReal dF = divCor[k][j][i].z;
                    if(dF < 0.0)
                        lambda[k][j][i].z = PetscMin(Rminus[k][j][i], Rplus [k+1][j][i]);
                    else
                        lambda[k][j][i].z = PetscMin(Rplus [k][j][i], Rminus[k+1][j][i]);
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivACor,  &divCor);
    DMDAVecRestoreArray(fda, aeqn->lLambdaA,  &lambda);
    DMDAVecRestoreArray(da,  aeqn->lRplusA,   &Rplus);
    DMDAVecRestoreArray(da,  aeqn->lRminusA,  &Rminus);

    DMLocalToLocalBegin(fda, aeqn->lLambdaA, INSERT_VALUES, aeqn->lLambdaA);
    DMLocalToLocalEnd  (fda, aeqn->lLambdaA, INSERT_VALUES, aeqn->lLambdaA);

    resetFacePeriodicFluxesVector(mesh, aeqn->lLambdaA, aeqn->lLambdaA, "localToLocal");

    // compute final flux  = low-order + lambda-weighted correction  
    DMDAVecGetArray(fda, aeqn->lDivA,    &div);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(fda, aeqn->lDivALO,  &divLO);
    DMDAVecGetArray(fda, aeqn->lDivACor, &divCor);
    DMDAVecGetArray(fda, aeqn->lLambdaA, &lambda);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                if(j!=0 && k!=0)
                    div[k][j][i].x = divLO [k][j][i].x + lambda[k][j][i].x   * divCor[k][j][i].x;

                if(i!=0 && k!=0)
                    div[k][j][i].y = divLO [k][j][i].y + lambda[k][j][i].y   * divCor[k][j][i].y;

                if(i!=0 && j!=0)
                    div[k][j][i].z = divLO [k][j][i].z + lambda[k][j][i].z   * divCor[k][j][i].z;
            }
        }
    }

    DMDAVecRestoreArray(fda,  aeqn->lDivA,    &div);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivALO,  &divLO);
    DMDAVecRestoreArray(fda, aeqn->lDivACor, &divCor);
    DMDAVecRestoreArray(fda, aeqn->lLambdaA, &lambda);

    DMLocalToLocalBegin(fda, aeqn->lDivA, INSERT_VALUES, aeqn->lDivA);
    DMLocalToLocalEnd  (fda, aeqn->lDivA, INSERT_VALUES, aeqn->lDivA);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivA, aeqn->lDivA, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormA(aeqn_ *aeqn, Vec Div, Vec Rhs)
{
    mesh_         *mesh  = aeqn->access->mesh;
    clock_        *clock = aeqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***rhs, ***aj;
    Cmpnts        ***div;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, Div,            &div);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(da,  Rhs,            &rhs);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhs[k][j][i] += aj[k][j][i] * 
                (
                    div [k][j][i].x - div [k][j][i-1].x +
                    div [k][j][i].y - div [k][j-1][i].y +
                    div [k][j][i].z - div [k-1][j][i].z
                );
            }
        }
    }

    DMDAVecRestoreArray(fda, Div,            &div);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(da,  Rhs,            &rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AeqnSNES(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    clock_        *clock = aeqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     norm;
    PetscInt      iter;
    PetscReal     ts, te;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "%s%s: Solving for Alpha-Water, Initial residual = ", aeqn->ddtScheme.c_str(), aeqn->solverType.c_str());

    // Step 1: implicit upwind BE/BDF2 solve → alpha* (smooth residual, Newton converges)
    VecCopy(aeqn->Alpha_o, aeqn->AlphaTmp);
    SNESSolve(aeqn->snesA, PETSC_NULL, aeqn->AlphaTmp);
    SNESGetFunctionNorm(aeqn->snesA, &norm);
    SNESGetIterationNumber(aeqn->snesA, &iter);

    PetscInt linIter = 0;
    SNESGetLinearSolveIterations(aeqn->snesA, &linIter);

    // save solution and sync ghost cells
    VecCopy(aeqn->AlphaTmp, aeqn->Alpha);
    DMGlobalToLocalBegin(da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM, "Final residual = %e, Iterations = %ld (linear = %ld), Elapsed Time = %lf\n", norm, iter, linIter, te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SNESFuncEvalA(SNES snes, Vec Alpha, Vec Rhs, void *ptr)
{
    aeqn_         *aeqn  = (aeqn_*)ptr;
    mesh_         *mesh  = aeqn->access->mesh;
    clock_        *clock = aeqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // scatter trial vector to local ghost-extended array
    VecCopy(Alpha, aeqn->Alpha);
    DMGlobalToLocalBegin(da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // compute low-order divergence fluxes
    FormALowOrder (aeqn);
    FormAHighOrder(aeqn);
    AddCompressionFlux(aeqn);
    ComputeLowOrderUpdate(aeqn);   
    ApplyMULESLimiter(aeqn);

    VecSet(Rhs, 0.0);
    FormA(aeqn, aeqn->lDivA, Rhs);

    // scale by dt
    VecScale(Rhs, clock->dt);

    // zero non resolved cells
    resetNonResolvedCellCentersScalar(mesh, Rhs);

    // add time-derivative residual terms
    if(aeqn->ddtScheme == "BDF2" && clock->it > clock->itStart)
    {
        VecAXPY(Rhs, -3.0/2.0, Alpha);
        VecAXPY(Rhs,  2.0,     aeqn->Alpha_o);
        VecAXPY(Rhs, -1.0/2.0, aeqn->Alpha_oo);
    }
    else
    {
        VecAXPY(Rhs, -1.0, Alpha);
        VecAXPY(Rhs,  1.0, aeqn->Alpha_o);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsA(aeqn_ *aeqn)
{
    mesh_ *mesh = aeqn->access->mesh;

    // ensure alpha periodic fluxes are consistent
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // build low-order and high-order face fluxes, apply mules limiter -> result in aeqn->Rhs
    FormALowOrder (aeqn);
    FormAHighOrder(aeqn);
    AddCompressionFlux(aeqn);
    ComputeLowOrderUpdate(aeqn);   
    ApplyMULESLimiter(aeqn);

    VecSet(aeqn->Rhs, 0.0);
    FormA(aeqn, aeqn->lDivA, aeqn->Rhs);

    // zero out non-resolved cells
    resetNonResolvedCellCentersScalar(mesh, aeqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AeqnRK4(aeqn_ *aeqn)
{
    mesh_  *mesh  = aeqn->access->mesh;
    clock_ *clock = aeqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for Alpha-Water, Stage ");

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

    // Alpha_o contribution
    VecCopy(aeqn->Alpha_o, aeqn->AlphaTmp);

    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate guess and evaluate RHS
        if(i!=0)
        {
            VecWAXPY(aeqn->Alpha, a[i] * dt, aeqn->Rhs, aeqn->Alpha_o);
            DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
            DMGlobalToLocalEnd(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
        }

        // compute function guess
        FormExplicitRhsA(aeqn);

        // add contribution from K1, K2, K3, K4
        VecAXPY(aeqn->AlphaTmp, dt * b[i], aeqn->Rhs);
    }

    VecCopy(aeqn->AlphaTmp, aeqn->Alpha);
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateRho(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    clock_        *clock = aeqn->access->clock;
    constants_    *constants = aeqn->access->constants;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     rhoAir   = constants->rho;
    PetscReal     rhoWater = constants->rhoWater;

    Cmpnts        ***ucont, ***rhoFace, 
                  ***div;

    PetscReal     ***alpha, ***rhoCell;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray( da, aeqn->lAlpha,   &alpha);
    DMDAVecGetArray(fda, aeqn->lDivA,    &div);
    DMDAVecGetArray(fda, aeqn->lRhoFace, &rhoFace);
    DMDAVecGetArray(fda, ueqn->lUcont,   &ucont);
    

    PetscScalar tol = 1.0e-12;

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                if(j!=0 && k!=0)
                {
                    PetscReal alpha_face;

                    // like OpenFOAM does it
                    //alpha_face =  - div[k][j][i].x  / ucont[k][j][i].x;

                    // harmonic mean
                    //alpha_face = 2.0 * alpha[k][j][i] * alpha[k][j][i+1] / (alpha[k][j][i] + alpha[k][j][i+1]);

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j][i+1]);
                    
                    // bound alpha
                    alpha_face = PetscMax(0.0, PetscMin(1.0, alpha_face));

                    // compute rho at cell faces
                    rhoFace[k][j][i].x  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;

                    if(rhoFace[k][j][i].x < tol)
                    {
                        printf("Warning: near-zero/negative density at csi face: rho(%ld,%ld,%ld) = %e\n", i, j, k, rhoFace[k][j][i].x);
                    }
                }

                if(i!=0 && k!=0)
                {
                    PetscReal alpha_face;

                    // like OpenFOAM does it
                    //alpha_face =  - div[k][j][i].y  / ucont[k][j][i].y;

                    // harmonic mean
                    //alpha_face = 2.0 * alpha[k][j][i] * alpha[k][j+1][i] / (alpha[k][j][i] + alpha[k][j+1][i]);

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j+1][i]);   

                    // bound alpha
                    alpha_face = PetscMax(0.0, PetscMin(1.0, alpha_face));

                    // compute rho at cell faces
                    rhoFace[k][j][i].y  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;

                    if(rhoFace[k][j][i].y < tol)
                    {
                        printf("Warning: near-zero/negative density at eta face: rho(%ld,%ld,%ld) = %e\n", i, j, k, rhoFace[k][j][i].y);
                    }
                }

                if(i!=0 && j!=0)
                {
                    PetscReal alpha_face;

                    // like OpenFOAM does it
                    // alpha_face =  - div[k][j][i].z  / ucont[k][j][i].z;

                    // harmonic mean
                    //alpha_face = 2.0 * alpha[k][j][i] * alpha[k+1][j][i] / (alpha[k][j][i] + alpha[k+1][j][i]);

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k+1][j][i]);
                    
                    // bound alpha
                    alpha_face = PetscMax(0.0, PetscMin(1.0, alpha_face));

                    // compute rho at cell faces
                    rhoFace[k][j][i].z  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;

                    if(rhoFace[k][j][i].z < tol)          
                    {
                        printf("Warning: near-zero/negative density at zeta face: rho(%ld,%ld,%ld) = %e\n", i, j, k, rhoFace[k][j][i].z);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray( da, aeqn->lAlpha,   &alpha);
    DMDAVecRestoreArray(fda, aeqn->lDivA,    &div);
    DMDAVecRestoreArray(fda, aeqn->lRhoFace, &rhoFace);
    DMDAVecRestoreArray(fda, ueqn->lUcont,   &ucont);

    // local to local scatter 
    DMLocalToLocalBegin(fda, aeqn->lRhoFace,  INSERT_VALUES, aeqn->lRhoFace);
    DMLocalToLocalEnd  (fda, aeqn->lRhoFace,  INSERT_VALUES, aeqn->lRhoFace);

    // handle periodicity
    resetFacePeriodicFluxesVector(mesh, aeqn->lRhoFace,  aeqn->lRhoFace,  "localToLocal");

    // interpolate rho at cell centers by averaging face values
    DMDAVecGetArray( da, aeqn->lRho,     &rhoCell);
    DMDAVecGetArray(fda, aeqn->lRhoFace, &rhoFace);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                rhoCell[k][j][i] = (rhoFace[k][j][i].x + rhoFace[k][j][i-1].x +
                                    rhoFace[k][j][i].y + rhoFace[k][j-1][i].y +
                                    rhoFace[k][j][i].z + rhoFace[k-1][j][i].z) / 6.0;

            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lRhoFace, &rhoFace);
    DMDAVecRestoreArray( da, aeqn->lRho,     &rhoCell);

    DMLocalToLocalBegin(da, aeqn->lRho,  INSERT_VALUES, aeqn->lRho);
    DMLocalToLocalEnd  (da, aeqn->lRho,  INSERT_VALUES, aeqn->lRho);

    // update boundary conditions (first cells and then faces)
    UpdateRhoBCs(aeqn);

    return(0);
}