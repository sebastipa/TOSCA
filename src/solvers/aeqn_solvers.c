//! \file  aeqn_solvers.c
//! \brief Alpha-Water-equation solver function implementations.

#include "aeqn_solvers.h"

PetscErrorCode FormALowOrder(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i,  j,  k,
                  ii, jj, kk;

    PetscReal     ***alpha, ***alphaSmooth;
    Cmpnts        ***ucont, ***div;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // first copy alpha into alpha smooth
    VecCopy(aeqn->lAlpha, aeqn->lAlphaSmooth);

    DMDAVecGetArray(da,  aeqn->lAlpha,  &alpha);
    DMDAVecGetArray(fda, ueqn->lUcont,  &ucont);
    DMDAVecGetArray(da,  aeqn->lAlphaSmooth, &alphaSmooth);

    VecSet(aeqn->lDivALO, 0.0);
    DMDAVecGetArray(fda, aeqn->lDivALO, &div);

    // compute first-order upwind face fluxes and store in lDivALO
    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(j!=0 && k!=0)
                    div[k][j][i].x = -ucont[k][j][i].x * upwind(alpha[k][j][i], alpha[k][j][i+1], ucont[k][j][i].x);

                if(i!=0 && k!=0)
                    div[k][j][i].y = -ucont[k][j][i].y * upwind(alpha[k][j][i], alpha[k][j+1][i], ucont[k][j][i].y);

                if(i!=0 && j!=0)
                    div[k][j][i].z = -ucont[k][j][i].z * upwind(alpha[k][j][i], alpha[k+1][j][i], ucont[k][j][i].z);

                // compute alpha smooth for subsequent compression term
                if(i!=0 && j!=0 && k!=0)
                {
                    PetscReal sum = 0.0;
                    for (kk=-1; kk<=1; kk++)
                    {
                        for (jj=-1; jj<=1; jj++)
                        {
                            for (ii=-1; ii<=1; ii++)
                            {
                                sum += alpha[k+kk][j+jj][i+ii];
                            }
                        }
                    }
                    alphaSmooth[k][j][i] = sum / 27.0;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivALO, &div);
    DMDAVecRestoreArray(da,  aeqn->lAlpha, &alpha);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(da,  aeqn->lAlphaSmooth, &alphaSmooth);

    // sync lDivALO across processors
    DMLocalToLocalBegin(fda, aeqn->lDivALO, INSERT_VALUES, aeqn->lDivALO);
    DMLocalToLocalEnd  (fda, aeqn->lDivALO, INSERT_VALUES, aeqn->lDivALO);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivALO, aeqn->lDivALO, "localToLocal");

    // sync alpha_smooth across processors
    DMLocalToLocalBegin(da, aeqn->lAlphaSmooth, INSERT_VALUES, aeqn->lAlphaSmooth);
    DMLocalToLocalEnd  (da, aeqn->lAlphaSmooth, INSERT_VALUES, aeqn->lAlphaSmooth);

    resetCellPeriodicFluxes(mesh, aeqn->lAlphaSmooth, aeqn->lAlphaSmooth, "scalar", "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormAHighOrder(aeqn_ *aeqn)
{
    mesh_         *mesh  = aeqn->access->mesh;
    ueqn_         *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscInt      i,  j,  k;

    PetscReal     dadcsi, dadeta, dadzet;                                        // alpha der. w.r.t. curvil. coords, 
    PetscReal     dadx, dady, dadz;                                              // alpha der. w.r.t. cart. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;     

    PetscReal     ***alpha, ***alphaSmooth, 
                  ***nvert, ***meshTag, ***aj;
    Cmpnts        ***ucont, ***div, ***gradAlpha,
                  ***csi , ***eta , ***zet ;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(aeqn->lDivAHO, 0.0);

    DMDAVecGetArray(fda, aeqn->lDivAHO,      &div);
    DMDAVecGetArray(da,  aeqn->lAlpha,       &alpha);
    DMDAVecGetArray(fda, ueqn->lUcont,       &ucont);
    DMDAVecGetArray(da,  mesh->lNvert,       &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,     &meshTag);
    DMDAVecGetArray(da,  mesh->lAj,          &aj);
    DMDAVecGetArray(fda, mesh->lCsi,         &csi);
    DMDAVecGetArray(fda, mesh->lEta,         &eta);
    DMDAVecGetArray(fda, mesh->lZet,         &zet);
    DMDAVecGetArray(fda, aeqn->lGradAlpha,   &gradAlpha);
    DMDAVecGetArray(da,  aeqn->lAlphaSmooth, &alphaSmooth);
    
    // compute high-order face fluxes and store in lDivAHO
    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
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
                        //div[k][j][i].x = -ucont[k][j][i].x * central4(alpha[k][j][iL], alpha[k][j][i], alpha[k][j][i+1], alpha[k][j][iR]);
                        div[k][j][i].x = -ucont[k][j][i].x * central(alpha[k][j][i], alpha[k][j][i+1]);
                        //div[k][j][i].x = -ucont[k][j][i].x * tvd12(alpha[k][j][iL], alpha[k][j][i], alpha[k][j][i+1], alpha[k][j][iR], ucont[k][j][i].x);
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
                        //div[k][j][i].y = -ucont[k][j][i].y * central4(alpha[k][jL][i], alpha[k][j][i], alpha[k][j+1][i], alpha[k][jR][i]);
                        div[k][j][i].y = -ucont[k][j][i].y * central(alpha[k][j][i], alpha[k][j+1][i]);
                        //div[k][j][i].y = -ucont[k][j][i].y * tvd12(alpha[k][jL][i], alpha[k][j][i], alpha[k][j+1][i], alpha[k][jR][i], ucont[k][j][i].y);
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
                        //div[k][j][i].z = -ucont[k][j][i].z * central4(alpha[kL][j][i], alpha[k][j][i], alpha[k+1][j][i], alpha[kR][j][i]);
                        div[k][j][i].z = -ucont[k][j][i].z * central(alpha[k][j][i], alpha[k+1][j][i]);
                        //div[k][j][i].z = -ucont[k][j][i].z * tvd12(alpha[kL][j][i], alpha[k][j][i], alpha[k+1][j][i], alpha[kR][j][i], ucont[k][j][i].z);
                    }
                }

                if(i!=0 && j!=0 && k!=0)
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
                        alphaSmooth, nvert, meshTag, 
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
    }

    DMDAVecRestoreArray(fda, aeqn->lDivAHO,      &div);
    DMDAVecRestoreArray(da,  aeqn->lAlpha,       &alpha);
    DMDAVecRestoreArray(fda, ueqn->lUcont,       &ucont);
    DMDAVecRestoreArray(da,  mesh->lNvert,       &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,     &meshTag);
    DMDAVecRestoreArray(fda, aeqn->lGradAlpha,   &gradAlpha);
    DMDAVecRestoreArray(da,  aeqn->lAlphaSmooth, &alphaSmooth);

    // sync lDivAHO across processors
    DMLocalToLocalBegin(fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);
    DMLocalToLocalEnd  (fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivAHO, aeqn->lDivAHO, "localToLocal");

    // sync ghost cells so neighboring processors can read GradA
    DMLocalToLocalBegin(fda, aeqn->lGradAlpha, INSERT_VALUES, aeqn->lGradAlpha);
    DMLocalToLocalEnd  (fda, aeqn->lGradAlpha, INSERT_VALUES, aeqn->lGradAlpha);

    // reset BCs in GradA
    UpdateGradABCs(aeqn);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ComputeLowOrderUpdate(aeqn_ *aeqn)
{
    // provisional low-order cell update: alpha_LO = alpha + dt * aj * divBudget_LO
    // Note: do NOT clip alphaLO as it must allow values outside [0,1] so that P_plus/P_minus
    // can go negative, enabling MULES to self-correct when low-order fluxes violate bounds

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
                alphaLO[k][j][i] = alpha[k][j][i] + dt * aj[k][j][i] * divBudget;
            }
        }
    }

    DMDAVecRestoreArray(da,  aeqn->lAlpha,   &alpha);
    DMDAVecRestoreArray(da,  aeqn->lAlphaLO, &alphaLO);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivALO,  &div);

    // scatter before periodicity reset
    DMLocalToLocalBegin(da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);
    DMLocalToLocalEnd  (da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);

    resetCellPeriodicFluxes(mesh, aeqn->lAlphaLO, aeqn->lAlphaLO, "scalar", "localToLocal");

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

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;                                           
    
    Cmpnts        ***icsi, ***jeta, ***kzet;

    PetscReal     ***alpha, ***nvert, ***meshTag;
    Cmpnts        ***ucont, ***div, ***gradAlpha;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // get numerical compression coefficient
    PetscReal compcoeff = aeqn->compCoeff; 

    DMDAVecGetArray(fda, mesh->lICsi,       &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,       &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,       &kzet);
    DMDAVecGetArray(da,  mesh->lNvert,      &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,    &meshTag);
    
    DMDAVecGetArray(fda, ueqn->lUcont,      &ucont);
    DMDAVecGetArray(da,  aeqn->lAlpha,      &alpha);
    DMDAVecGetArray(fda, aeqn->lDivAHO,     &div);
    DMDAVecGetArray(fda, aeqn->lGradAlpha,  &gradAlpha);

    // compute interface velocity at cell faces as uc_face = c * mag_u_face * grad_alpha_face / mag_grad_alpha_face
    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                // i-faces (exclude IBM faces)
                if(j!=0 && k!=0 && !isOversetIFace(k, j, i, i+1, meshTag))
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

                    PetscReal raw_flux = nDot(ucatc_face, icsi[k][j][i]);
                    PetscReal max_flux = fabs(ucont[k][j][i].x);

                    if (fabs(raw_flux) > max_flux)
                    {
                        raw_flux = (raw_flux > 0.0 ? max_flux : -max_flux);
                    }

                    div[k][j][i].x -= alpha_face*(1.0 - alpha_face) * raw_flux;
                }

                // j-faces (exclude IBM faces)
                if(i!=0 && k!=0 && !isOversetJFace(k, j, i, j+1, meshTag))
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

                    PetscReal raw_flux = nDot(ucatc_face, jeta[k][j][i]);
                    PetscReal max_flux = fabs(ucont[k][j][i].y);

                    if (fabs(raw_flux) > max_flux)
                    {
                        raw_flux = (raw_flux > 0.0 ? max_flux : -max_flux);
                    }

                    div[k][j][i].y -= alpha_face*(1.0 - alpha_face) * raw_flux;
                }

                // k-faces (exclude IBM faces)
                if (i!=0 && j!=0 && !isOversetKFace(k, j, i, k+1, meshTag))
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

                    PetscReal raw_flux = nDot(ucatc_face, kzet[k][j][i]);
                    PetscReal max_flux = fabs(ucont[k][j][i].z);

                    if (fabs(raw_flux) > max_flux)
                    {
                        raw_flux = (raw_flux > 0.0 ? max_flux : -max_flux);
                    }

                    div[k][j][i].z -= alpha_face*(1.0 - alpha_face) * raw_flux;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,       &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,       &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,       &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert,      &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,    &meshTag);
    
    DMDAVecRestoreArray(fda, ueqn->lUcont,      &ucont);
    DMDAVecRestoreArray(da,  aeqn->lAlpha,      &alpha);
    DMDAVecRestoreArray(fda, aeqn->lDivAHO,     &div);
    DMDAVecRestoreArray(fda, aeqn->lGradAlpha,  &gradAlpha);

    DMLocalToLocalBegin(fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);
    DMLocalToLocalEnd  (fda, aeqn->lDivAHO, INSERT_VALUES, aeqn->lDivAHO);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivAHO, aeqn->lDivAHO, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ApplyMULESLimiter(aeqn_ *aeqn)
{
    // apply mules flux limiter to produce a bounded, conservative fluxes
    // step 1: correction = lDivAHO - lDivA
    // step 2: provisional low-order cell update into lAlphaLO
    // step 3: computation of R_plus, R_minus (how much correction each cell can absorb)
    // step 4: computation of lambda_f = min(R of owner, R of neighbor)
    // step 5: final fluxes = low-order fluxes + lambda-limited high-order fluxes

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
    Cmpnts        ***div, ***divLO, ***divCor, 
                  ***lambda;

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

    // enforce periodic consistency for R values at periodic boundaries
    resetCellPeriodicFluxes(mesh, aeqn->lRplusA,  aeqn->lRplusA,  "scalar", "localToLocal");
    resetCellPeriodicFluxes(mesh, aeqn->lRminusA, aeqn->lRminusA, "scalar", "localToLocal");

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

    // compute final flux = low-order + lambda-weighted correction  
    DMDAVecGetArray(fda, aeqn->lDivA,    &div);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(fda, aeqn->lDivALO,  &divLO);
    DMDAVecGetArray(fda, aeqn->lDivACor, &divCor);
    DMDAVecGetArray(fda, aeqn->lLambdaA, &lambda);

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(j!=0 && k!=0)
                    div[k][j][i].x = divLO [k][j][i].x + lambda[k][j][i].x   * divCor[k][j][i].x;

                if(i!=0 && k!=0)
                    div[k][j][i].y = divLO [k][j][i].y + lambda[k][j][i].y   * divCor[k][j][i].y;

                if(i!=0 && j!=0)
                    div[k][j][i].z = divLO [k][j][i].z + lambda[k][j][i].z   * divCor[k][j][i].z;
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivA,    &div);
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
    ueqn_        *ueqn  = aeqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***rhs, ***aj, ***lalpha;
    Cmpnts        ***div, ***lucont;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, Div,            &div);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(da,  Rhs,            &rhs);
    DMDAVecGetArray(fda, ueqn->lUcont,     &lucont);
    DMDAVecGetArray(da,  aeqn->lAlpha,     &lalpha);

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

                // correction for non-divergence free cells (this should be all zero)
                rhs[k][j][i] += aj[k][j][i] * lalpha[k][j][i] *
                (
                    lucont [k][j][i].x - lucont [k][j][i-1].x +
                    lucont [k][j][i].y - lucont [k][j-1][i].y +
                    lucont [k][j][i].z - lucont [k-1][j][i].z
                );
            }
        }
    }

    DMDAVecRestoreArray(fda, Div,            &div);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(da,  Rhs,            &rhs);
    DMDAVecRestoreArray(fda, ueqn->lUcont,   &lucont);
    DMDAVecRestoreArray(da,  aeqn->lAlpha,   &lalpha);

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

    // cumulate fluxes
    VecAXPY(aeqn->lDivMean, 1.0/aeqn->nAlphaSubCycles, aeqn->lDivA); 

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
    
    // ensure alpha periodic fluxes are consistent
    enforceInteriorCellPeriodicity(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar");
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // compute low-order divergence fluxes
    FormALowOrder (aeqn);
    FormAHighOrder(aeqn);
    AddCompressionFlux(aeqn);
    ComputeLowOrderUpdate(aeqn);   
    ApplyMULESLimiter(aeqn);

    VecSet(Rhs, 0.0);
    FormA(aeqn, aeqn->lDivA, Rhs);

    if(aeqn->access->flags->isKLeftAlphaDampingActive || aeqn->access->flags->isKRightAlphaDampingActive)
    {
        dampingSourceA(aeqn, Rhs, 1.0);
    }

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

    // scatter alpha to ensure ghost cells are updated
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    // ensure alpha periodic fluxes are consistent
    // enforceInteriorCellPeriodicity(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar");
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // build low-order and high-order face fluxes, apply mules limiter -> result in aeqn->Rhs
    FormALowOrder (aeqn);
    FormAHighOrder(aeqn);
    AddCompressionFlux(aeqn);
    ComputeLowOrderUpdate(aeqn);   
    ApplyMULESLimiter(aeqn);

    VecSet(aeqn->Rhs, 0.0);
    FormA(aeqn, aeqn->lDivA, aeqn->Rhs);

    if
    (
        aeqn->access->flags->isKLeftAlphaDampingActive || 
        aeqn->access->flags->isKRightAlphaDampingActive
    )
    {
        dampingSourceA(aeqn, aeqn->Rhs, 1.0);
    }
    
    // zero out non-resolved cells
    resetNonResolvedCellCentersScalar(mesh, aeqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AeqnRK3(aeqn_ *aeqn)
{
    // SSP-RK3 (Strong Stability Preserving Runge-Kutta 3rd order)
    // Third-order accurate, guarantees boundedness when MULES limits each stage
    // 
    // Mathematical form:
    //   Stage 1: u1 = u0 + dt*L(u0)
    //   Stage 2: u2 = 3/4*u0 + 1/4*(u1 + dt*L(u1))
    //   Stage 3: u_new = 1/3*u0 + 2/3*(u2 + dt*L(u2))
    //
    // Each stage: compute full forward Euler step (dtCoeff=1.0), 
    // then take convex combination with u0 to preserve boundedness

    mesh_  *mesh  = aeqn->access->mesh;
    clock_ *clock = aeqn->access->clock;
    PetscReal ts, te, dt = clock->dt;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "SSP-RK3: Solving for Alpha-Water, Stage ");

    // Stage 1: Alpha_1 = Alpha_o + dt*L(Alpha_o)
    PetscPrintf(mesh->MESH_COMM, "1, ");
    
    FormExplicitRhsA(aeqn);
    VecWAXPY(aeqn->Alpha, dt, aeqn->Rhs, aeqn->Alpha_o);  // Alpha = Alpha_o + dt*L
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    // cumulate fluxes
    VecAXPY(aeqn->lDivMean, 1.0/6.0/aeqn->nAlphaSubCycles, aeqn->lDivA); 

    // Stage 2: Alpha_2 = 3/4*Alpha_o + 1/4*(Alpha_1 + dt*L(Alpha_1))
    PetscPrintf(mesh->MESH_COMM, "2, ");
    
    FormExplicitRhsA(aeqn);  
    VecAXPY (aeqn->Alpha, dt, aeqn->Rhs);               // Alpha = Alpha_1 + dt*L
    VecAXPBY(aeqn->Alpha, 0.75, 0.25, aeqn->Alpha_o);  // Alpha = 0.75*Alpha_o + 0.25*Alpha
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    // cumulate fluxes
    VecAXPY(aeqn->lDivMean, 1.0/6.0/aeqn->nAlphaSubCycles, aeqn->lDivA);
    
    // Stage 3: Alpha_new = 1/3*Alpha_o + 2/3*(Alpha_2 + dt*L(Alpha_2)) 
    PetscPrintf(mesh->MESH_COMM, "3, ");

    FormExplicitRhsA(aeqn);
    VecAXPY (aeqn->Alpha, dt, aeqn->Rhs);                     // Alpha = Alpha_2 + dt*L
    VecAXPBY(aeqn->Alpha, 1.0/3.0, 2.0/3.0, aeqn->Alpha_o);  // Alpha = (1/3)*Alpha_o + (2/3)*Alpha
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    // cumulate fluxes
    VecAXPY(aeqn->lDivMean, 2.0/3.0/aeqn->nAlphaSubCycles, aeqn->lDivA);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM, "Elapsed Time = %f\n", te-ts);

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

    Cmpnts        ***rhoFace;

    PetscReal     ***alpha, ***rhoCell;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray( da, aeqn->lAlpha,   &alpha);
    DMDAVecGetArray(fda, aeqn->lRhoFace, &rhoFace);
    
    PetscScalar tol = 1.0e-12;

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(j!=0 && k!=0)
                {
                    PetscReal alpha_face;

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j][i+1]);

                    rhoFace[k][j][i].x  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;
                }

                if(i!=0 && k!=0)
                {
                    PetscReal alpha_face;

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k][j+1][i]);  

                    rhoFace[k][j][i].y  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;
                }

                if(i!=0 && j!=0)
                {
                    PetscReal alpha_face;

                    // algebraic mean
                    alpha_face = 0.5 * (alpha[k][j][i] + alpha[k+1][j][i]);

                    rhoFace[k][j][i].z  = alpha_face * rhoWater + (1.0 - alpha_face) * rhoAir;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lRhoFace, &rhoFace);

    // local to local scatter 
    DMLocalToLocalBegin(fda, aeqn->lRhoFace,  INSERT_VALUES, aeqn->lRhoFace);
    DMLocalToLocalEnd  (fda, aeqn->lRhoFace,  INSERT_VALUES, aeqn->lRhoFace);

    // handle periodicity
    resetFacePeriodicFluxesVector(mesh, aeqn->lRhoFace,  aeqn->lRhoFace,  "localToLocal");

    // compute rho at cell centers 
    DMDAVecGetArray( da, aeqn->lRho,     &rhoCell);
    DMDAVecGetArray(fda, aeqn->lRhoFace, &rhoFace);

    for (k=zs; k<ze; k++)
    {
        for (j=lys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                // this is a smoothening: more stable by less accurate
                //rhoCell[k][j][i] = (rhoFace[k][j][i].x + rhoFace[k][j][i-1].x +
                //                    rhoFace[k][j][i].y + rhoFace[k][j-1][i].y +
                //                    rhoFace[k][j][i].z + rhoFace[k-1][j][i].z) / 6.0;
                
                // this is the exact density arising from the alpha solution
                rhoCell[k][j][i] = alpha[k][j][i] * rhoWater + (1.0 - alpha[k][j][i]) * rhoAir;
            }
        }
    }

    DMDAVecRestoreArray( da, aeqn->lAlpha,   &alpha);
    DMDAVecRestoreArray(fda, aeqn->lRhoFace, &rhoFace);
    DMDAVecRestoreArray( da, aeqn->lRho,     &rhoCell);

    DMLocalToLocalBegin(da, aeqn->lRho,  INSERT_VALUES, aeqn->lRho);
    DMLocalToLocalEnd  (da, aeqn->lRho,  INSERT_VALUES, aeqn->lRho);

    // rho BCs update is not required as this field is derived from alpha which is alredy in sync

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode dampingSourceA(aeqn_ *aeqn, Vec &Rhs, PetscReal scale)
{
    mesh_         *mesh  = aeqn->access->mesh;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***rhs, ***alpha;
    Cmpnts        ***cent;

    PetscReal     z, y, x, xmin, xmax, nud, alphaBar;
    PetscReal     start, end, delta, coeff;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    if(aeqn->access->flags->isKLeftAlphaDampingActive)
    {
        start = aeqn->kLeftAlphaDampingStart;
        end   = aeqn->kLeftAlphaDampingEnd;
        delta = aeqn->kLeftAlphaDampingDelta;
        coeff = aeqn->kLeftAlphaDampingCoeff;
    }

    DMDAVecGetArray(da,  Rhs,   &rhs);
    DMDAVecGetArray(da,  aeqn->lAlpha, &alpha);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // compute cell center x at i,j,k
                x     = cent[k][j][i].x;
                y     = cent[k][j][i].y;
                z     = cent[k][j][i].z;

                if(aeqn->access->flags->isKLeftAlphaDampingActive)
                {
                    // here we have to prescribe the wave in half a meter 
                    PetscReal t      = mesh->access->clock->time;
                    PetscReal d      = aeqn->waveLevel - mesh->bounds.zmin;
                    PetscReal theta  = aeqn->waveDirection;
                    PetscReal z_surface, eta_instant;

                    if (aeqn->waveType == "linear" || aeqn->waveType == "stokes2")
                    {
                        // wave parameters
                        PetscReal H      = aeqn->waveHeight;
                        PetscReal k_wave = aeqn->waveNumber;
                        PetscReal omega  = aeqn->waveOmega;
                        PetscReal phi    = aeqn->wavePhase;
                        
                        // wave phase argument
                        PetscReal arg = k_wave * (x * std::cos(theta) + y * std::sin(theta)) - omega * t + phi;
                        
                        // compute instantaneous wave elevation at this (x,y,t)
                        if (aeqn->waveType == "linear")
                        {
                            eta_instant = (H/2.0) * std::cos(arg);
                        }
                        else  if (aeqn->waveType == "stokes2")
                        {
                            PetscReal sinh_kd  = std::sinh(      k_wave * d);
                            PetscReal cosh_kd  = std::cosh(      k_wave * d);
                            PetscReal cosh_2kd = std::cosh(2.0 * k_wave * d);
                            PetscReal eta1     = (H/2.0) * std::cos(arg);
                            PetscReal eta2     = (k_wave * H * H / 8.0) * cosh_kd / (sinh_kd * sinh_kd * sinh_kd) 
                                                * (2.0 + cosh_2kd) * std::cos(2.0 * arg);
                            eta_instant        = eta1 + eta2;
                        }
                        
                        // compute the instantaneous surface elevation
                        z_surface = aeqn->waveLevel + eta_instant;
                    }
                    else if (aeqn->waveType == "dispersiveFocusing")
                    {
                        eta_instant = 0.0;

                        for(PetscInt g=0; g<aeqn->nGroups; g++)
                        {
                            PetscReal A      = aeqn->groupAmplitudes[g];
                            PetscReal k_wave = aeqn->groupWaveNumbers[g];
                            PetscReal omega  = aeqn->groupOmegas[g];
                            PetscReal phi    = aeqn->groupPhases[g];

                            // wave phase argument
                            PetscReal arg = omega * t - k_wave * (x * std::cos(theta) + y * std::sin(theta)) + phi;

                            // compute instantaneous wave elevation at this (x,y,t)
                            eta_instant += A * std::cos(arg);
                        }

                        // compute the instantaneous surface elevation
                        z_surface = aeqn->waveLevel + eta_instant;
                    }

                    alphaBar = (z < z_surface) ? 1.0 : 0.0;
                    nud      = viscNordstrom(coeff, start, end, delta, x);

                    rhs[k][j][i] += scale * nud * (alphaBar - alpha[k][j][i]);
                }
                else if(aeqn->access->flags->isKRightAlphaDampingActive)
                {
                    xmin  = mesh->bounds.xmax - aeqn->kRightAlphaDampingDelta;
                    xmax  = mesh->bounds.xmax;

                    // bar state 
                    if(z < aeqn->kRightAlphaDampingWaterLevel)
                    {
                        alphaBar = 1.0;
                    }
                    else
                    {
                        alphaBar = 0.0;
                    }

                    // compute Cosine viscosity 
                    nud   = viscCosAscending(aeqn->kRightAlphaDampingCoeff, xmin, xmax, x);

                    rhs[k][j][i] += scale * nud * (alphaBar - alpha[k][j][i]);
                }
            }
        }
    }

    DMDAVecRestoreArray(da,  Rhs,   &rhs);
    DMDAVecRestoreArray(da,  aeqn->lAlpha, &alpha);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    return(0);
}