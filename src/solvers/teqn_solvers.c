//! \file  teqn_solvers.c
//! \brief Temperature-equation solver function implementations.

#include "teqn_solvers.h"

PetscErrorCode FormT(teqn_ *teqn, Vec &Rhs, PetscReal scale, PetscInt formMode)
{
    // In this function the viscous + divergence term of the temperature equation are
    // discretized at cell centers.
    // First the divergence and viscous fluxes are evaluated at cell faces, then
    // their budget is evaluated at the internal cells, forming the Rhs.

    mesh_         *mesh  = teqn->access->mesh;
    ueqn_         *ueqn  = teqn->access->ueqn;
    les_          *les   = teqn->access->les;
    constants_    *cst   = teqn->access->constants;
    flags_        *flags = teqn->access->flags;
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

    PetscReal     ***tmprt, ***rhs, ***nvert, ***meshTag;

    Cmpnts        ***div, ***visc, ***viscIBM;                                                // divergence and viscous terms
    Cmpnts        ***limiter;                                                     // flux limiter
    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;                                  // cell and face jacobians
    PetscReal     ***lnu_t, ***Sabs, ***lch, ***lcs, ***lkt;

    PetscReal     dtdc, dtde, dtdz;                                              // velocity der. w.r.t. curvil. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components

    PetscReal     tRef;
    if(ueqn->access->flags->isAblActive) tRef = teqn->access->abl->tRef;
    else                                 tRef = teqn->access->constants->tRef;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  Rhs,               &rhs);
    DMDAVecGetArray(da,  teqn->lTmprt,      &tmprt);
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
    DMDAVecGetArray(da,  mesh->lmeshTag,    &meshTag);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    if(flags->isIBMActive)
    {
        DMDAVecGetArray(fda, teqn->lViscIBMT, &viscIBM);
    }

    VecSet(teqn->lDivT,  0.0);
    VecSet(teqn->lViscT, 0.0);;

    DMDAVecGetArray(fda, teqn->lDivT, &div);
    DMDAVecGetArray(fda, teqn->lViscT, &visc);

    if (teqn->access->flags->isLesActive)
    {
        if(les->model != STABILITY_BASED)
        {
            DMDAVecGetArray(da, les->lNu_t, &lnu_t);

            if(les->model == AMD   || 
               les->model == DSM   || 
               les->model == DLASI || 
               les->model == DLASD || 
               les->model == DPASD ||
               les->model == BAMD)            
            {
                DMDAVecGetArray(da, les->lk_t, &lkt);
            }
        }
        else
        {
            DMDAVecGetArray(da,  les->lS,  &Sabs);
            DMDAVecGetArray(da,  les->lCh, &lch);
            DMDAVecGetArray(da,  les->lCs, &lcs);
        }
    }

    // ---------------------------------------------------------------------- //
    //         FORM DIVERGENCE AND VISCOUS CONTRIBUTIONS (fused i/j/k)       //
    // ---------------------------------------------------------------------- //

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                // ---- i-face (skip when j==0 or k==0) ----
                if(j!=0 && k!=0)
                {
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

                    Compute_dscalar_i(mesh, i, j, k, mx, my, mz, tmprt, nvert, meshTag, &dtdc, &dtde, &dtdz);

                    PetscInt  iL, iR;
                    PetscReal denom;
                    getFace2Cell4StencilCsi(mesh, k, j, i, mx, &iL, &iR, &denom, nvert, meshTag);

                    div[k][j][i].x =
                    - ucont[k][j][i].x
                    * centralUpwind
                    (
                        tmprt[k][j][iL],
                        tmprt[k][j][i],
                        tmprt[k][j][i+1],
                        tmprt[k][j][iR],
                        ucont[k][j][i].x
                        ,limiter[k][j][i].x
                    );

                    PetscReal nu = cst->nu, nut;
                    PetscReal kappaEff;

                    if(teqn->access->flags->isLesActive)
                    {
                        nu = 0;

                        if(les->model != STABILITY_BASED)
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                            if( les->model == AMD   ||
                                les->model == DSM   ||
                                les->model == DLASI ||
                                les->model == DLASD ||
                                les->model == DPASD ||
                                les->model == BAMD)
                            {
                                PetscReal diff = 0.5 * (lkt[k][j][i] + lkt[k][j][i+1]);
                                kappaEff = (nu / cst->Pr) + diff;
                            }
                            else
                            {
                                PetscReal gradTdotG = dtde*(-9.81);
                                PetscReal l, delta = pow(1./ajc, 1./3.);
                                if(gradTdotG < 0.)
                                    l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                                else
                                    l = delta;
                                PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));
                                kappaEff = (nu / cst->Pr) + (nut / Prt);
                            }
                        }
                        else
                        {
                            kappaEff = (nu / cst->Pr) + 0.5 *
                            (
                                lch[k][j][i+1] / 0.1 * Sabs[k][j][i+1] * lcs[k][j][i+1] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }

                        // wall model i-left/right patch
                        if
                        (
                            (mesh->boundaryT.iLeft=="thetaWallFunction" && i==0) ||
                            (mesh->boundaryT.iRight=="thetaWallFunction" && i==mx-2)
                        )
                        {
                            PetscReal signQ = 1.0;
                            if(i==0) signQ = -1.0;
                            visc[k][j][i].x =
                            signQ *
                            (
                                teqn->iRWM->qWall.x[k-zs][j-ys] * icsi[k][j][i].x +
                                teqn->iRWM->qWall.y[k-zs][j-ys] * icsi[k][j][i].y +
                                teqn->iRWM->qWall.z[k-zs][j-ys] * icsi[k][j][i].z
                            );
                            kappaEff = 0.0;
                        }

                        // IBM wall model
                        if(isIBMFluidIFace(k, j, i, i+1, nvert))
                        {
                            if(teqn->access->ibm->wallShearOn && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                            {
                                if(isIBMFluidCell(k, j, i, nvert))
                                    visc[k][j][i].x = viscIBM[k][j][i].x;
                                else if(isIBMFluidCell(k, j, i+1, nvert))
                                    visc[k][j][i].x = viscIBM[k][j][i+1].x;
                                kappaEff = 0;
                            }
                        }
                    }
                    else
                    {
                        kappaEff = nu / cst->Pr;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc[k][j][i].x += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * kappaEff;
                }

                // ---- j-face (skip when i==0 or k==0) ----
                if(i!=0 && k!=0)
                {
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

                    Compute_dscalar_j(mesh, i, j, k, mx, my, mz, tmprt, nvert, meshTag, &dtdc, &dtde, &dtdz);

                    PetscInt  jL, jR;
                    PetscReal denom;
                    getFace2Cell4StencilEta(mesh, k, j, i, my, &jL, &jR, &denom, nvert, meshTag);

                    div[k][j][i].y =
                    - ucont[k][j][i].y
                    * centralUpwind
                    (
                        tmprt[k][jL][i],
                        tmprt[k][j][i],
                        tmprt[k][j+1][i],
                        tmprt[k][jR][i],
                        ucont[k][j][i].y
                        ,limiter[k][j][i].y
                    );

                    PetscReal nu = cst->nu, nut;
                    PetscReal kappaEff;

                    if(teqn->access->flags->isLesActive)
                    {
                        nu = 0;

                        if(les->model != STABILITY_BASED)
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                            if( les->model == AMD   ||
                                les->model == DSM   ||
                                les->model == DLASI ||
                                les->model == DLASD ||
                                les->model == DPASD ||
                                les->model == BAMD)
                            {
                                PetscReal diff = 0.5 * (lkt[k][j][i] + lkt[k][j+1][i]);
                                kappaEff = (nu / cst->Pr) + diff;
                            }
                            else
                            {
                                PetscReal gradTdotG = dtde*(-9.81);
                                PetscReal l, delta = pow(1./ajc, 1./3.);
                                if(gradTdotG < 0.)
                                    l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                                else
                                    l = delta;
                                PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));
                                kappaEff = (nu / cst->Pr) + (nut / Prt);
                            }
                        }
                        else
                        {
                            kappaEff = (nu / cst->Pr) + 0.5 *
                            (
                                lch[k][j+1][i] / 0.1 * Sabs[k][j+1][i] * lcs[k][j+1][i] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }

                        // wall model j-left/right patch
                        if
                        (
                            (mesh->boundaryT.jLeft=="thetaWallFunction"  && j==0) ||
                            (mesh->boundaryT.jRight=="thetaWallFunction" && j==my-2)
                        )
                        {
                            PetscReal signQ = 1.0;
                            if(j==0) signQ = -1.0;
                            visc[k][j][i].y =
                            signQ *
                            (
                                teqn->jLWM->qWall.x[k-zs][i-xs] * jeta[k][j][i].x +
                                teqn->jLWM->qWall.y[k-zs][i-xs] * jeta[k][j][i].y +
                                teqn->jLWM->qWall.z[k-zs][i-xs] * jeta[k][j][i].z
                            );
                            kappaEff = 0.0;
                        }

                        if(isIBMFluidJFace(k, j, i, j+1, nvert) && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                        {
                            if(teqn->access->ibm->wallShearOn)
                            {
                                if(isIBMFluidCell(k, j, i, nvert))
                                    visc[k][j][i].y = viscIBM[k][j][i].y;
                                else if(isIBMFluidCell(k, j+1, i, nvert))
                                    visc[k][j][i].y = viscIBM[k][j+1][i].y;
                                kappaEff = 0;
                            }
                        }
                    }
                    else
                    {
                        kappaEff = nu / cst->Pr;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc[k][j][i].y += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * kappaEff;
                }

                // ---- k-face (skip when i==0 or j==0) ----
                if(i!=0 && j!=0)
                {
                    // get 1/V at the k-face
                    PetscReal ajc = kaj[k][j][i];

                    // get face normals
                    csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                    eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                    zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                    // compute metric tensor
                    g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                    g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                    g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                    Compute_dscalar_k(mesh, i, j, k, mx, my, mz, tmprt, nvert, meshTag, &dtdc, &dtde, &dtdz);

                    PetscInt  kL, kR;
                    PetscReal denom;
                    getFace2Cell4StencilZet(mesh, k, j, i, mz, &kL, &kR, &denom, nvert, meshTag);

                    div[k][j][i].z =
                    - ucont[k][j][i].z
                    * centralUpwind
                    (
                        tmprt[kL][j][i],
                        tmprt[k][j][i],
                        tmprt[k+1][j][i],
                        tmprt[kR][j][i],
                        ucont[k][j][i].z
                        ,limiter[k][j][i].z
                    );

                    PetscReal nu = cst->nu, nut;
                    PetscReal kappaEff;

                    if(teqn->access->flags->isLesActive)
                    {
                        nu = 0;

                        if(les->model != STABILITY_BASED)
                        {
                            nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

                            if( les->model == AMD   ||
                                les->model == DSM   ||
                                les->model == DLASI ||
                                les->model == DLASD ||
                                les->model == DPASD ||
                                les->model == BAMD)
                            {
                                PetscReal diff = 0.5 * (lkt[k][j][i] + lkt[k+1][j][i]);
                                kappaEff = (nu / cst->Pr) + diff;
                            }
                            else
                            {
                                PetscReal gradTdotG = dtde*(-9.81);
                                PetscReal l, delta = pow(1./ajc, 1./3.);
                                if(gradTdotG < 0.)
                                    l = PetscMin(delta, 7.6*nut/delta*std::sqrt(tRef / std::fabs(gradTdotG)));
                                else
                                    l = delta;
                                PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));
                                kappaEff = (nu / cst->Pr) + (nut / Prt);
                            }
                        }
                        else
                        {
                            kappaEff = (nu / cst->Pr) + 0.5 *
                            (
                                lch[k+1][j][i] / 0.1 * Sabs[k+1][j][i] * lcs[k+1][j][i] +
                                lch[k][j][i]   / 0.1 * Sabs[k][j][i]   * lcs[k][j][i]
                            ) * pow(1.0 / ajc, 1.0/3.0);
                        }

                        if(isIBMFluidKFace(k, j, i, k+1, nvert) && teqn->access->ibm->ibmBody[0]->tempBC == "thetaWallFunction")
                        {
                            if(teqn->access->ibm->wallShearOn)
                            {
                                if(isIBMFluidCell(k, j, i, nvert))
                                    visc[k][j][i].z = viscIBM[k][j][i].z;
                                else if(isIBMFluidCell(k+1, j, i, nvert))
                                    visc[k][j][i].z = viscIBM[k+1][j][i].z;
                                kappaEff = 0;
                            }
                        }
                    }
                    else
                    {
                        kappaEff = nu / cst->Pr;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc[k][j][i].z += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * kappaEff;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, teqn->lDivT, &div);
    DMDAVecRestoreArray(fda, teqn->lViscT, &visc);

    if(flags->isIBMActive)
    {
        DMDAVecRestoreArray(fda, teqn->lViscIBMT, &viscIBM);
    }

    DMLocalToLocalBegin(fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalEnd  (fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalBegin(fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);
    DMLocalToLocalEnd  (fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);

    if(teqn->access->flags->isLesActive && (les->model == BAMD || les->model == BV))
    {
        updateLESScalarStructuralModel(les); 
    }

    resetFacePeriodicFluxesVector(mesh, teqn->lDivT,   teqn->lDivT, "localToLocal");
    resetFacePeriodicFluxesVector(mesh, teqn->lViscT, teqn->lViscT, "localToLocal");

    DMDAVecGetArray(fda, teqn->lDivT, &div);
    DMDAVecGetArray(fda, teqn->lViscT, &visc);

    // ---------------------------------------------------------------------- //
    //                      FORM THE CUMULATIVE FLUXES                        //
    // ---------------------------------------------------------------------- //

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal divBudget =
                    div[k][j][i].x  - div[k  ][j  ][i-1].x       +
                    div[k][j][i].y  - div[k  ][j-1][i  ].y       +
                    div[k][j][i].z  - div[k-1][j  ][i  ].z;

                PetscReal viscBudget =
                    visc[k][j][i].x - visc[k  ][j  ][i-1].x      +
                    visc[k][j][i].y - visc[k  ][j-1][i  ].y      +
                    visc[k][j][i].z - visc[k-1][j  ][i  ].z;

                rhs[k][j][i]
                =
                scale *
                (
                    (formMode != 2 ? divBudget  : 0.0) +
                    (formMode != 1 ? viscBudget : 0.0)
                ) * aj[k][j][i];
            }
        }
    }

    if(teqn->access->flags->isLesActive)
    {
        if(les->model != STABILITY_BASED)
        {
            DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);

            if( les->model == AMD   || 
                les->model == DSM   || 
                les->model == DLASI || 
                les->model == DLASD || 
                les->model == DPASD ||
                les->model == BAMD)             
            {
                DMDAVecRestoreArray(da, les->lk_t, &lkt);
            }
        }
        else
        {
            DMDAVecRestoreArray(da,  les->lS,  &Sabs);
            DMDAVecRestoreArray(da,  les->lCh, &lch);
            DMDAVecRestoreArray(da,  les->lCs, &lcs);
        }
    }

    DMDAVecRestoreArray(fda, teqn->lDivT,       &div);
    DMDAVecRestoreArray(fda, teqn->lViscT,      &visc);

    DMDAVecRestoreArray(da,  Rhs,              &rhs);
    DMDAVecRestoreArray(da,  teqn->lTmprt,      &tmprt);
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
    DMDAVecRestoreArray(da,  mesh->lmeshTag,    &meshTag);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode IMEXTMatVec(Mat A, Vec v, Vec Av)
{
    // IMEX MatShell operator: A*v = v - dt*D(v)
    // D(v) is the diffusion-only operator evaluated at the current iterate v.
    // Used by kspIMEX to solve the linear IMEX system A*T^{n+1} = T^n + dt*bT.
    // Called by PETSc KSP for every matrix-vector product during the iterative solve.
    // Note: backward Euler for diffusion (scale=1.0 always) avoids velocity-temperature coupling issues of CN.

    teqn_  *teqn;
    MatShellGetContext(A, (void**)&teqn);

    mesh_  *mesh  = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;
    PetscReal dt  = clock->dt;

    // sync v → Tmprt/lTmprt so FormT reads the correct state
    VecCopy(v, teqn->Tmprt);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // reset periodic BCs
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // compute D(v) into scratch buffer Rhs (diffusion-only, formMode=2, scale=1.0 backward Euler)
    VecSet(teqn->Rhs, 0.0);
    FormT(teqn, teqn->Rhs, 1.0, 2);
    VecScale(teqn->Rhs, dt);                // Rhs = dt*D(v)
    resetNonResolvedCellCentersScalar(mesh, teqn->Rhs);

    // Av = v - dt*D(v)
    VecCopy(v, Av);
    VecAXPY(Av, -1.0, teqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TeqnIMEX(teqn_ *teqn)
{
    mesh_  *mesh  = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;

    PetscReal ts, te;
    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "IMEX-CNAB: Solving for T, ");

    // Step 1: ensure T^n state is synchronised (lTmprt current)
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // Step 2: convective RHS at T^n → RhsConv (conv-only, formMode=1)
    VecSet(teqn->RhsConv, 0.0);
    FormT(teqn, teqn->RhsConv, 1.0, 1);

    // Step 3: build bT = AB2 convection + dampingSourceT + sourceT/dt
    //         (no explicit diffusion: full backward-Euler diffusion is handled by the linear operator A)
    VecSet(teqn->bT, 0.0);

    // AB2 convection extrapolation: c1=3/2, c2=-1/2 (Forward Euler on first step)
    if(clock->it > clock->itStart)
    {
        VecAXPY(teqn->bT,  1.5, teqn->RhsConv);     // (3/2) * Conv^n
        VecAXPY(teqn->bT, -0.5, teqn->RhsConv_o);   // (-1/2) * Conv^{n-1}
    }
    else
    {
        VecAXPY(teqn->bT, 1.0, teqn->RhsConv);      // Forward Euler on first step
    }

    // sponge/Rayleigh damping contribution (T/time units, added to bT)
    if(teqn->access->flags->isXDampingActive)
        dampingSourceT(teqn, teqn->bT, 1.0);

    // Step 4: build RHS b = T^n + dt*bT
    //         sourceT is added WITHOUT dt scaling, matching TeqnSNES convention
    VecCopy(teqn->Tmprt_o, teqn->TmprtTmp);
    VecAXPY(teqn->TmprtTmp, clock->dt, teqn->bT);

    if(teqn->access->flags->isAblActive)
        if(teqn->access->abl->controllerActiveT)
            sourceT(teqn, teqn->TmprtTmp, 1.0);

    resetNonResolvedCellCentersScalar(mesh, teqn->TmprtTmp);

    // Step 5: KSP direct solve: A*T^{n+1} = TmprtTmp
    // RhsConv_o is used as the solution buffer because:
    //   - its AB2 history value was already consumed into bT above,
    //   - IMEXTMatVec writes to Tmprt/lTmprt as workspace (different from the solution buffer),
    //   - Step 7 will overwrite RhsConv_o with RhsConv, so order matters.
    VecCopy(teqn->TmprtTmp, teqn->RhsConv_o);         // non-zero initial guess
    KSPSolve(teqn->kspIMEX, teqn->TmprtTmp, teqn->RhsConv_o);

    PetscReal norm;  PetscInt iter;
    KSPGetResidualNorm(teqn->kspIMEX, &norm);
    KSPGetIterationNumber(teqn->kspIMEX, &iter);

    // Step 6: copy solution back and scatter
    VecCopy(teqn->RhsConv_o, teqn->Tmprt);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // Step 7: rotate Conv buffer for AB2 at the next step (must follow Step 6)
    VecCopy(teqn->RhsConv, teqn->RhsConv_o);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,
        "Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n",
        norm, iter, te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TeqnSNES(SNES snes, Vec T, Vec Rhs, void *ptr)
{
    teqn_ *teqn   = (teqn_*)ptr;
    mesh_ *mesh   = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;
    VecCopy(T, teqn->Tmprt);

    // scatter temperature from global to local
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // reset temperature periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // update wall model (optional)
    // UpdateWallModelsT(teqn);

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    const PetscReal dt = clock->dt;

    // add viscous and transport terms.
    // implicit-half scale: 0.5 for Crank-Nicolson (explicit half prebuilt in bT), 1.0 for backward Euler
    FormT(teqn, Rhs, 1.0);

    if(teqn->access->flags->isXDampingActive)
    {
        // this is causing spurious oscillation: deactivate
        dampingSourceT(teqn, Rhs, 1.0);
    }

    // multiply for dt: 0.5*dt for CN (implicit half only), dt for backward Euler
    PetscReal dtScale = (teqn->ddtScheme == "crankNicholson") ? 0.5 : 1.0;
    VecScale(Rhs, dtScale * dt);

    // for CN: add prebuilt explicit half bT = 0.5*dt*(FormT + damp)(T^n)
    if(teqn->ddtScheme == "crankNicholson")
        VecAXPY(Rhs, 1.0, teqn->bT);

    // add driving source terms after as it is not scaled by 1/dt
    if(teqn->access->flags->isAblActive)
    {
        if(teqn->access->abl->controllerActiveT)
        {
            sourceT(teqn, Rhs, 1.0);
        }
    }

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, Rhs);

    VecAXPY(Rhs, -1.0, T);
    VecAXPY(Rhs,  1.0, teqn->Tmprt_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsT(teqn_ *teqn)
{
    mesh_ *mesh   = teqn->access->mesh;

    // reset temperature periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, teqn->Tmprt, teqn->lTmprt, "scalar", "globalToLocal");

    // update wall model (optional)
    // UpdateWallModelsT(teqn);

    // initialize the rhs vector
    VecSet(teqn->Rhs, 0.0);

    // add viscous and transport terms
    FormT(teqn, teqn->Rhs, 1.0);

    if(teqn->access->flags->isXDampingActive)
    {
        // this is causing spurious oscillation: deactivate
        dampingSourceT(teqn, teqn->Rhs, 1.0);
    }

    // add driving source terms after as it is not scaled by 1/dt
    if(teqn->access->flags->isAblActive)
    {
        if(teqn->access->abl->controllerActiveT)
        {
            sourceT(teqn, teqn->Rhs, 1.0 / teqn->access->clock->dt);
        }
    }

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, teqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode TeqnRK4(teqn_ *teqn)
{
    mesh_  *mesh  = teqn->access->mesh;
    clock_ *clock = teqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for T, Stage ");

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

    // Tmprt_o contribution
    VecCopy(teqn->Tmprt_o, teqn->TmprtTmp);

    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate U guess and evaluate RHS
        if(i!=0)
        {
            VecWAXPY(teqn->Tmprt, a[i] * dt, teqn->Rhs, teqn->Tmprt_o);
            DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
            DMGlobalToLocalEnd(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
        }

        // compute function guess
        FormExplicitRhsT(teqn);

        // add contribution from K1, K2, K3, K4
        VecAXPY(teqn->TmprtTmp, dt * b[i], teqn->Rhs);
    }

    VecCopy(teqn->TmprtTmp, teqn->Tmprt);
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

