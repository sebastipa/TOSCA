//! \file  ueqn_solvers.c
//! \brief Momentum-equation solver function implementations.

#include "ueqn_solvers.h"

PetscErrorCode FormU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale, PetscInt formMode)
{
    // In this function the viscous + divergence term of the momentum equation are
    // discretized at cell faces. The disposition is staggered, meaning that the
    // component of the face-centered quantities are not located at the same point
    // but rather at the corresponding faces.
    // First at every cell face (in the 1st 2nd and 3rd curvilinear direction) the
    // fluxes are evaluated in cartesian coordinates. Secondly the net flux is evaluated
    // for every cell in cartesian coordinates, summing up the previously calculated fluxes.
    // Third, this cell centered flux is interpolated at the faces and dotted with that face's
    // area vector, meaning that it will point in that face's curiviliear coordinate direction.
    // Note: if the flow is non-periodic in a given direction, only the contravariant velocity at
    // the internal faces is solved for. If the flow is periodic in that direction also the
    // contravariant velocity at the right boundary faces is solved for and fluxes at the left
    // boundary faces are calculated as if that face was the right boundary face (as if left and
    // right boundaries were the same boundary).

    // formMode: 0=full (conv+visc, default), 1=conv-only, 2=visc-only 

    mesh_            *mesh   = ueqn->access->mesh;
    les_             *les    = ueqn->access->les;
    constants_       *cst    = ueqn->access->constants;
    flags_           *flags  = ueqn->access->flags;
    DM               da      = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info    = mesh->info;
    PetscInt         xs      = info.xs, xe = info.xs + info.xm;
    PetscInt         ys      = info.ys, ye = info.ys + info.ym;
    PetscInt         zs      = info.zs, ze = info.zs + info.zm;
    PetscInt         mx      = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont, ***ucat;

    Cmpnts           ***csi,  ***eta,  ***zet;
    Cmpnts           ***icsi, ***ieta, ***izet;
    Cmpnts           ***jcsi, ***jeta, ***jzet;
    Cmpnts           ***kcsi, ***keta, ***kzet;
    Cmpnts           ***cent;

    PetscReal        ***nvert, ***lnu_t, ***meshTag;

    Cmpnts           ***div1,  ***div2,  ***div3;                               // divergence & cumulative fluxes
    Cmpnts           ***visc1, ***visc2, ***visc3;
    Cmpnts           ***viscIBM1, ***viscIBM2, ***viscIBM3;                     // viscous terms                             // viscous terms
    Cmpnts           ***rhs,   ***fp,    ***limiter;                            // right hand side
    PetscReal        ***aj,    ***iaj,   ***jaj, ***kaj;                        // cell and face jacobians

    PetscReal        dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;      // velocity der. w.r.t. curvil. coords
    PetscReal        csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;      // surface area vectors components
    PetscReal        g11, g21, g31;                                             // metric tensor components
    PetscReal        r11, r21, r31,
                     r12, r22, r32,
                     r13, r23, r33;

    PetscScalar      solid = 0.5;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, Rhs,  &rhs);

    // get fundamental distributed arrays
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);

    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);

    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lIAj, &iaj);
    DMDAVecGetArray(da,  mesh->lJAj, &jaj);
    DMDAVecGetArray(da,  mesh->lKAj, &kaj);

    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    // ---------------------------------------------------------------------- //
    // FORM DIVERGENCE AND VISCOUS CONTRIBUTIONS                              //
    // ---------------------------------------------------------------------- //

    VecSet(ueqn->lFp,   0.0);
    if(formMode != 2) { VecSet(ueqn->lDiv1,  0.0); VecSet(ueqn->lDiv2,  0.0); VecSet(ueqn->lDiv3,  0.0); }
    if(formMode != 1) { VecSet(ueqn->lVisc1, 0.0); VecSet(ueqn->lVisc2, 0.0); VecSet(ueqn->lVisc3, 0.0); }

    // get distributed arrays
    DMDAVecGetArray(fda, ueqn->lDiv1, &div1);
    DMDAVecGetArray(fda, ueqn->lDiv2, &div2);
    DMDAVecGetArray(fda, ueqn->lDiv3, &div3);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);

    if(flags->isIBMActive)
    {
        DMDAVecGetArray(fda, ueqn->lViscIBM1, &viscIBM1);
        DMDAVecGetArray(fda, ueqn->lViscIBM2, &viscIBM2);
        DMDAVecGetArray(fda, ueqn->lViscIBM3, &viscIBM3);
    }

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
    }

    // damping viscosity for fringe region advection damping
    PetscReal nuD, nuDY;

    // fringe region parameters (set only if active)
    PetscReal xS, yS;
    PetscReal xE, yE;
    PetscReal xDS, yDS;
    PetscReal xDE, yDE;

    PetscInt  advectionDamping = 0, advectionDampingY = 0;
	PetscReal advDampH = 0, advDampYH = 0;

    if(ueqn->access->flags->isAdvectionDampingXActive)
    {
        xS     = ueqn->access->abl->advDampingXStart;
        xE     = ueqn->access->abl->advDampingXEnd;
        xDE    = ueqn->access->abl->advDampingXDeltaEnd;
        xDS    = ueqn->access->abl->advDampingXDeltaStart;

        advectionDamping = 1;

        advDampH = ueqn->access->abl->hInv - 0.5*ueqn->access->abl->dInv;
    }
    else
    {
        nuD = 1.0;
    }

    if(ueqn->access->flags->isAdvectionDampingYActive)
    {
        yS     = ueqn->access->abl->advDampingYStart;
        yE     = ueqn->access->abl->advDampingYEnd;
        yDE    = ueqn->access->abl->advDampingYDeltaEnd;
        yDS    = ueqn->access->abl->advDampingYDeltaStart;

        advectionDampingY = 1;

        advDampYH = ueqn->access->abl->hInv - 0.5*ueqn->access->abl->dInv;
    }
    else
    {
        nuDY = 1.0;
    }

    // FUSED LOOP: i/j/k direction faces computed in a single pass for improved cache reuse.
    // Each cell's ucat, ucont, nvert etc. are loaded once and all 3 face contributions
    // (div1/visc1, div2/visc2, div3/visc3) are computed before moving to the next cell.
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;

                // i faces 
                if(j!=0 && k!=0)
                {
                    // viscous setup (skip when conv-only: fields only needed for visc)
                    PetscReal ajc = 0.0;
                    if(formMode != 1)
                    {
                        ajc = iaj[k][j][i];

                        // get face normals
                        csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;
                        eta0 = ieta[k][j][i].x, eta1 = ieta[k][j][i].y, eta2 = ieta[k][j][i].z;
                        zet0 = izet[k][j][i].x, zet1 = izet[k][j][i].y, zet2 = izet[k][j][i].z;

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_i (mesh, i, j, k, mx, my, mz, ucat, nvert, meshTag, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                        // compute metric tensor - WARNING: there is a factor of 1/J^2 if using face area vectors
                        //                                  must multiply for ajc in viscous term!!!
                        g11 = csi0 * csi0 + csi1 * csi1 + csi2 * csi2;
                        g21 = eta0 * csi0 + eta1 * csi1 + eta2 * csi2;
                        g31 = zet0 * csi0 + zet1 * csi1 + zet2 * csi2;

                        // compute cartesian velocity derivatives w.r.t. cartesian coords
                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;    //du_dx / J -> another factor of / J is added in the viscous term
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;    //dv_dx / J -> another factor of / J is added in the viscous term
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;    //dw_dx / J -> another factor of / J is added in the viscous term

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
                    }

                    // divergence (skip when visc-only)
                    if(formMode != 2)
                    {
                    PetscInt    iL, iR;
                    PetscReal denom;
                    getFace2Cell4StencilCsi(mesh, k, j, i, mx, &iL, &iR, &denom, nvert, meshTag);

                    // test: inviscid flow or weno3Div
                    if(ueqn->inviscid || ueqn->weno3Div)
                    {
                        div1[k][j][i] = nScale
                        (
                            - ucont[k][j][i].x,
                            weno3Vec
                            (
                                ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                ucont[k][j][i].x
                            )
                        );
                    }
                    else
                    {
                        // central scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k][j][i+1]
                                )
                            );
                        }
                        else if(ueqn->central4Div)
                        {
                            
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                centralVec4thCsi(mesh, k, j, i, mx, nvert, meshTag, ucat, ucont[k][j][i].x, ueqn->hyperVisc)
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                centralUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    ucont[k][j][i].x, limiter[k][j][i].x
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][j][iL ] * nMag(csi[k][j][iL ]));
                            PetscReal d1 = 1.0 / (aj[k][j][i  ] * nMag(csi[k][j][i  ]));
                            PetscReal d2 = 1.0 / (aj[k][j][i+1] * nMag(csi[k][j][i+1]));
                            PetscReal d3 = 1.0 / (aj[k][j][iR ] * nMag(csi[k][j][iR ]));

                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                wCentralUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].x, limiter[k][j][i].x
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div1[k][j][i] = nScale
                            (
                                - ucont[k][j][i].x,
                                quadraticUpwindVec
                                (
                                    ucat[k][j][iL], ucat[k][j][i], ucat[k][j][i+1], ucat[k][j][iR],
                                    ucont[k][j][i].x
                                )
                            );
                        }
                    }
                    } // end formMode != 2

                    // viscous accumulation (skip when conv-only)
                    if(formMode != 1)
                    {
                    PetscReal nuEff, nu = cst->nu, nut;

                    // viscous terms
                    if
                    (
                        ueqn->access->flags->isLesActive
                    )
                    {

                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                        // wall model i-left patch
                        if
                        (
                                (mesh->boundaryU.iLeft=="velocityWallFunction"  && i==0) ||
                                (mesh->boundaryU.iRight=="velocityWallFunction" && i==mx-2)
                        )
                        {
                            visc1[k][j][i].x = - ueqn->iLWM->tauWall.x[k-zs][j-ys];
                            visc1[k][j][i].y = - ueqn->iLWM->tauWall.y[k-zs][j-ys];
                            visc1[k][j][i].z = - ueqn->iLWM->tauWall.z[k-zs][j-ys];

                            nuEff = nu;
                        }
                        // slip boundary condition on U (set nuEff to 0)
                        else if
                        (
                            (mesh->boundaryU.iLeft =="slip" && i==0   ) ||
                            (mesh->boundaryU.iRight=="slip" && i==mx-2)
                        )
                        {
                            nuEff = 0.0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

                        if(isIBMFluidIFace(k, j, i, i+1, nvert))
                        {
                            if(ueqn->access->ibm->wallShearOn)
                            {
                                if(isIBMFluidCell(k, j, i, nvert))
                                {
                                    visc1[k][j][i] = nSet(viscIBM1[k][j][i]);
                                }
                                else if(isIBMFluidCell(k, j, i+1, nvert))
                                {
                                    visc1[k][j][i] = nSet(viscIBM1[k][j][i+1]);
                                }

                                nuEff = 0;
                            }
                            else
                            {
                                nuEff = nu + nut;
                            }

                        }
                    }
                    else
                    {
                        nuEff = nu;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc1[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * csi0 + r21 * csi1 + r31 * csi2) * ajc * (nuEff);
                    visc1[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * csi0 + r22 * csi1 + r32 * csi2) * ajc * (nuEff);
                    visc1[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * csi0 + r23 * csi1 + r33 * csi2) * ajc * (nuEff);
                    } // end formMode != 1
                }

                // j faces 
                if(i!=0 && k!=0)
                {
                    // viscous setup (skip when conv-only)
                    PetscReal ajc = 0.0;
                    if(formMode != 1)
                    {
                        ajc = jaj[k][j][i];

                        // get face normals
                        csi0 = jcsi[k][j][i].x, csi1 = jcsi[k][j][i].y, csi2 = jcsi[k][j][i].z;
                        eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;
                        zet0 = jzet[k][j][i].x, zet1 = jzet[k][j][i].y, zet2 = jzet[k][j][i].z;

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_j (mesh, i, j, k, mx, my, mz, ucat, nvert, meshTag, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                        // compute metric tensor
                        g11 = csi0 * eta0 + csi1 * eta1 + csi2 * eta2;
                        g21 = eta0 * eta0 + eta1 * eta1 + eta2 * eta2;
                        g31 = zet0 * eta0 + zet1 * eta1 + zet2 * eta2;

                        // compute cartesian velocity derivatives w.r.t. cartesian coords
                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
                    }

                    // divergence (skip when visc-only)
                    if(formMode != 2)
                    {
                    PetscInt    jL, jR;
                    PetscReal denom;
                    getFace2Cell4StencilEta(mesh, k, j, i, my, &jL, &jR, &denom, nvert, meshTag);

                    // test: inviscid flow or weno3Div
                    if( ueqn->inviscid || ueqn->weno3Div)
                    {
                        div2[k][j][i] = nScale
                        (
                            - ucont[k][j][i].y,
                            weno3Vec
                            (
                                ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                ucont[k][j][i].y
                            )
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k][j+1][i]
                                )
                            );
                        }
                        else if(ueqn->central4Div)
                        {
                            
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                centralVec4thEta(mesh, k, j, i, my, nvert, meshTag, ucat, ucont[k][j][i].y, ueqn->hyperVisc)
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                centralUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    ucont[k][j][i].y, limiter[k][j][i].y
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[k][jL ][i] * nMag(eta[k][jL ][i]));
                            PetscReal d1 = 1.0 / (aj[k][j  ][i] * nMag(eta[k][j  ][i]));
                            PetscReal d2 = 1.0 / (aj[k][j+1][i] * nMag(eta[k][j+1][i]));
                            PetscReal d3 = 1.0 / (aj[k][jR ][i] * nMag(eta[k][jR ][i]));

                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                wCentralUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].y, limiter[k][j][i].y
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div2[k][j][i] = nScale
                            (
                                - ucont[k][j][i].y,
                                quadraticUpwindVec
                                (
                                    ucat[k][jL][i], ucat[k][j][i], ucat[k][j+1][i], ucat[k][jR][i],
                                    ucont[k][j][i].y
                                )
                            );
                        }
                    }
                    } // end formMode != 2

                    // viscous accumulation (skip when conv-only)
                    if(formMode != 1)
                    {
                    PetscReal nuEff, nu = cst->nu, nut;

                    if
                    (
                        ueqn->access->flags->isLesActive
                    )
                    {

                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                        // wall model j-left patch
                        if
                        (
                                (mesh->boundaryU.jLeft=="velocityWallFunction" && j==0) ||
                                (mesh->boundaryU.jRight=="velocityWallFunction" && j==my-2)
                        )
                        {
                            visc2[k][j][i].x = - ueqn->jLWM->tauWall.x[k-zs][i-xs];
                            visc2[k][j][i].y = - ueqn->jLWM->tauWall.y[k-zs][i-xs];
                            visc2[k][j][i].z = - ueqn->jLWM->tauWall.z[k-zs][i-xs];

                            nuEff = nu;
                        }
                        // slip boundary condition on U (set nuEff to 0)
                        else if
                        (
                            (mesh->boundaryU.jLeft =="slip" && j==0   ) ||
                            (mesh->boundaryU.jRight=="slip" && j==my-2)
                        )
                        {
                            nuEff = 0.0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

                        if(isIBMFluidJFace(k, j, i, j+1, nvert))
                        {
                            if(ueqn->access->ibm->wallShearOn)
                            {

                                if(isIBMFluidCell(k, j, i, nvert))
                                {
                                    visc2[k][j][i] = nSet(viscIBM2[k][j][i]);
                                }
                                else if(isIBMFluidCell(k, j+1, i, nvert))
                                {
                                    visc2[k][j][i] = nSet(viscIBM2[k][j+1][i]);
                                }
                               
                                nuEff = 0;
                            }
                            else
                            {
                                nuEff = nu + nut;
                            }

                        }
                    }
                    else
                    {
                        nuEff = nu;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc2[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * eta0 + r21 * eta1 + r31 * eta2) * ajc * (nuEff);
                    visc2[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * eta0 + r22 * eta1 + r32 * eta2) * ajc * (nuEff);
                    visc2[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * eta0 + r23 * eta1 + r33 * eta2) * ajc * (nuEff);
                    } // end formMode != 1
                }

                // k faces 
                if(i!=0 && j!=0)
                {
                    // viscous setup (skip when conv-only)
                    PetscReal ajc = 0.0;
                    if(formMode != 1)
                    {
                        ajc = kaj[k][j][i];

                        // get face normals
                        csi0 = kcsi[k][j][i].x, csi1 = kcsi[k][j][i].y, csi2 = kcsi[k][j][i].z;
                        eta0 = keta[k][j][i].x, eta1 = keta[k][j][i].y, eta2 = keta[k][j][i].z;
                        zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                        // compute cartesian velocity derivatives w.r.t. curvilinear coords
                        Compute_du_k (mesh, i, j, k, mx, my, mz, ucat, nvert, meshTag, &dudc, &dvdc, &dwdc, &dude, &dvde, &dwde, &dudz, &dvdz, &dwdz);

                        // compute metric tensor
                        g11 = csi0 * zet0 + csi1 * zet1 + csi2 * zet2;
                        g21 = eta0 * zet0 + eta1 * zet1 + eta2 * zet2;
                        g31 = zet0 * zet0 + zet1 * zet1 + zet2 * zet2;

                        // compute cartesian velocity derivatives w.r.t. cartesian coords
                        r11 = dudc * csi0 + dude * eta0 + dudz * zet0;
                        r21 = dvdc * csi0 + dvde * eta0 + dvdz * zet0;
                        r31 = dwdc * csi0 + dwde * eta0 + dwdz * zet0;

                        r12 = dudc * csi1 + dude * eta1 + dudz * zet1;
                        r22 = dvdc * csi1 + dvde * eta1 + dvdz * zet1;
                        r32 = dwdc * csi1 + dwde * eta1 + dwdz * zet1;

                        r13 = dudc * csi2 + dude * eta2 + dudz * zet2;
                        r23 = dvdc * csi2 + dvde * eta2 + dvdz * zet2;
                        r33 = dwdc * csi2 + dwde * eta2 + dwdz * zet2;
                    }

                    // divergence (skip when visc-only)
                    if(formMode != 2)
                    {
                    PetscInt    kL, kR;
                    PetscReal denom;
                    getFace2Cell4StencilZet(mesh, k, j, i, mz, &kL, &kR, &denom, nvert, meshTag);

                    // inviscid flow or weno3Div
                    if(ueqn->inviscid || ueqn->weno3Div)
                    {
                        div3[k][j][i] = nScale
                        (
                            - ucont[k][j][i].z,
                            weno3Vec
                            (
                                ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                ucont[k][j][i].z
                            )
                        );
                    }
                    else
                    {
                        // second order divergence scheme
                        if(ueqn->centralDiv)
                        {
                            // ucat is interpolated at the face
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                centralVec
                                (
                                    ucat[k][j][i], ucat[k+1][j][i]
                                )
                            );
                        }
                        else if(ueqn->central4Div)
                        {
                            
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                centralVec4thZet(mesh, k, j, i, mz, nvert, meshTag, ucat, ucont[k][j][i].z, ueqn->hyperVisc)
                            );
                        }
                        else if(ueqn->centralUpwindDiv)
                        {
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                centralUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    ucont[k][j][i].z, limiter[k][j][i].z
                                )
                            );
                        }
                        else if(ueqn->centralUpwindWDiv)
                        {
                            // compute cell widths
                            PetscReal d0 = 1.0 / (aj[kL ][j][i] * nMag(zet[kL ][j][i]));
                            PetscReal d1 = 1.0 / (aj[k  ][j][i] * nMag(zet[k  ][j][i]));
                            PetscReal d2 = 1.0 / (aj[k+1][j][i] * nMag(zet[k+1][j][i]));
                            PetscReal d3 = 1.0 / (aj[kR ][j][i] * nMag(zet[kR ][j][i]));

                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                wCentralUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    d0, d1, d2, d3,
                                    ucont[k][j][i].z, limiter[k][j][i].z
                                )
                            );
                        }
                        // quickDiv scheme (3rd order upwind)
                        else if(ueqn->quickDiv)
                        {
                            div3[k][j][i] = nScale
                            (
                                - ucont[k][j][i].z,
                                quadraticUpwindVec
                                (
                                    ucat[kL][j][i], ucat[k][j][i], ucat[k+1][j][i], ucat[kR][j][i],
                                    ucont[k][j][i].z
                                )
                            );
                        }
                    }
                    } // end formMode != 2

                    // viscous accumulation (skip when conv-only)
                    if(formMode != 1)
                    {
                    PetscReal nuEff, nu = cst->nu, nut;

                    if
                    (
                        ueqn->access->flags->isLesActive
                    )
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

                        // slip boundary condition on U (set nuEff to 0)
                        if
                        (
                            (mesh->boundaryU.kLeft =="slip" && k==0   ) ||
                            (mesh->boundaryU.kRight=="slip" && k==mz-2)
                        )
                        {
                            nuEff = 0.0;
                        }
                        else
                        {
                            nuEff = nu + nut;
                        }

                        if(isIBMFluidKFace(k, j, i, k+1, nvert))
                        {
                            if(ueqn->access->ibm->wallShearOn)
                            {
                                if(isIBMFluidCell(k, j, i, nvert))
                                {
                                    visc3[k][j][i] = nSet(viscIBM3[k][j][i]);
                                }
                                else if(isIBMFluidCell(k+1, j, i, nvert))
                                {
                                    visc3[k][j][i] = nSet(viscIBM3[k+1][j][i]);
                                }

                                nuEff = 0;
                            }
                            else
                            {
                                nuEff = nu + nut;
                            }

                        }
                    }
                    else
                    {
                        nuEff = nu;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc3[k][j][i].x += (g11 * dudc + g21 * dude + g31 * dudz + r11 * zet0 + r21 * zet1 + r31 * zet2) * ajc * (nuEff);
                    visc3[k][j][i].y += (g11 * dvdc + g21 * dvde + g31 * dvdz + r12 * zet0 + r22 * zet1 + r32 * zet2) * ajc * (nuEff);
                    visc3[k][j][i].z += (g11 * dwdc + g21 * dwde + g31 * dwdz + r13 * zet0 + r23 * zet1 + r33 * zet2) * ajc * (nuEff);
                    } // end formMode != 1
                }
            }
        }
    }

    // restore distributed arrays and scatter local to local
    DMDAVecRestoreArray(fda, ueqn->lDiv1, &div1);
    DMDAVecRestoreArray(fda, ueqn->lDiv2, &div2);
    DMDAVecRestoreArray(fda, ueqn->lDiv3, &div3);

    // restore all arrays before posting scatters so all 6 can be pipelined
    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);

    if(flags->isIBMActive)
    {
        DMDAVecRestoreArray(fda, ueqn->lViscIBM1, &viscIBM1);
        DMDAVecRestoreArray(fda, ueqn->lViscIBM2, &viscIBM2);
        DMDAVecRestoreArray(fda, ueqn->lViscIBM3, &viscIBM3);
    }

    // update inter-processor ghost values for divergence and viscous contributions 
    DMLocalToLocalBegin(fda, ueqn->lDiv1,  INSERT_VALUES, ueqn->lDiv1);
    DMLocalToLocalEnd  (fda, ueqn->lDiv1,  INSERT_VALUES, ueqn->lDiv1);
    DMLocalToLocalBegin(fda, ueqn->lDiv2,  INSERT_VALUES, ueqn->lDiv2);
    DMLocalToLocalEnd  (fda, ueqn->lDiv2,  INSERT_VALUES, ueqn->lDiv2);
    DMLocalToLocalBegin(fda, ueqn->lDiv3,  INSERT_VALUES, ueqn->lDiv3);
    DMLocalToLocalEnd  (fda, ueqn->lDiv3,  INSERT_VALUES, ueqn->lDiv3);
    DMLocalToLocalBegin(fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalEnd  (fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalBegin(fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalEnd  (fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalBegin(fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);
    DMLocalToLocalEnd  (fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);

    // update LES strucural model in contravariant form
    if(ueqn->access->flags->isLesActive && (les->model == BAMD || les->model == BV))
    {
        updateLESStructuralModelContravariantForm(les); 
    }

    // ---------------------------------------------------------------------- //
    // FORM THE RIGHT HAND SIDE
    // ---------------------------------------------------------------------- //

    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv1,   ueqn->lDiv1,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv2,   ueqn->lDiv2,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lDiv3,   ueqn->lDiv3,   "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc1,  ueqn->lVisc1,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc2,  ueqn->lVisc2,  "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ueqn->lVisc3,  ueqn->lVisc3,  "localToLocal");

    // get the arrays
    DMDAVecGetArray(fda, ueqn->lDiv1, &div1);
    DMDAVecGetArray(fda, ueqn->lDiv2, &div2);
    DMDAVecGetArray(fda, ueqn->lDiv3, &div3);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);

    DMDAVecGetArray(fda, ueqn->lFp, &fp);

    // loop over cells and compute fp (viscous + divergence contributions)
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // divergence contribution (formMode == 2: visc-only, skip convection)
                if(formMode != 2)
                {
                fp[k][j][i].x
                +=
                (
                    div1[k][j][i].x - div1[k][j][i-1].x +
                    div2[k][j][i].x - div2[k][j-1][i].x +
                    div3[k][j][i].x - div3[k-1][j][i].x
                );

                fp[k][j][i].y
                +=
                (
                    div1[k][j][i].y - div1[k][j][i-1].y +
                    div2[k][j][i].y - div2[k][j-1][i].y +
                    div3[k][j][i].y - div3[k-1][j][i].y
                );

                fp[k][j][i].z
                +=
                (
                    div1[k][j][i].z - div1[k][j][i-1].z +
                    div2[k][j][i].z - div2[k][j-1][i].z +
                    div3[k][j][i].z - div3[k-1][j][i].z
                );

                // damp convective term in vertical momentum equation for streamwise fringe
				if(advectionDamping)
				{
					PetscReal height = cent[k][j][i].z - mesh->grndLevel;
					nuD              = viscStipaDelta(xS, xE, xDS, xDE, cent[k][j][i].x, height, advDampH);
					fp[k][j][i].z    = nuD * fp[k][j][i].z;
				}

                // damp convective term in vertical momentum equation for lateral fringe
                if(advectionDampingY)
				{
					PetscReal height = cent[k][j][i].z - mesh->grndLevel;
					nuDY             = viscStipaDelta(yS, yE, yDS, yDE, cent[k][j][i].y, height, advDampYH);
                    fp[k][j][i].z    = nuDY * fp[k][j][i].z;
				}
                } // end formMode != 2

                // viscous contribution (formMode == 1: conv-only, skip viscosity)
                if(formMode != 1 && !ueqn->inviscid)
                {
                    fp[k][j][i].x
                    +=
                    (
                        visc1[k][j][i].x - visc1[k][j][i-1].x +
                        visc2[k][j][i].x - visc2[k][j-1][i].x +
                        visc3[k][j][i].x - visc3[k-1][j][i].x
                    );

                    fp[k][j][i].y
                    +=
                    (
                        visc1[k][j][i].y - visc1[k][j][i-1].y +
                        visc2[k][j][i].y - visc2[k][j-1][i].y +
                        visc3[k][j][i].y - visc3[k-1][j][i].y
                    );

                    fp[k][j][i].z
                    +=
                    (
                        visc1[k][j][i].z - visc1[k][j][i-1].z +
                        visc2[k][j][i].z - visc2[k][j-1][i].z +
                        visc3[k][j][i].z - visc3[k-1][j][i].z
                    );
                }

                // handle periodicity regardless of formMode
                if( (i == 1) && !(mesh->i_periodic) && !(mesh->ii_periodic))
                {
                    fp[k][j][i-1].x = fp[k][j][i].x;
                    fp[k][j][i-1].y = fp[k][j][i].y;
                    fp[k][j][i-1].z = fp[k][j][i].z;
                }

                if( (j == 1) && !(mesh->j_periodic) && !(mesh->jj_periodic))
                {
                    fp[k][j-1][i].x = fp[k][j][i].x;
                    fp[k][j-1][i].y = fp[k][j][i].y;
                    fp[k][j-1][i].z = fp[k][j][i].z;
                }

                if( (k == 1) && !(mesh->k_periodic) && !(mesh->kk_periodic))
                {
                    fp[k-1][j][i].x = fp[k][j][i].x;
                    fp[k-1][j][i].y = fp[k][j][i].y;
                    fp[k-1][j][i].z = fp[k][j][i].z;
                }

                if( (i == mx-2) && !(mesh->i_periodic) && !(mesh->ii_periodic))
                {
                    fp[k][j][i+1].x = fp[k][j][i].x;
                    fp[k][j][i+1].y = fp[k][j][i].y;
                    fp[k][j][i+1].z = fp[k][j][i].z;
                }

                if( (j == my-2) && !(mesh->j_periodic) && !(mesh->jj_periodic))
                {
                    fp[k][j+1][i].x = fp[k][j][i].x;
                    fp[k][j+1][i].y = fp[k][j][i].y;
                    fp[k][j+1][i].z = fp[k][j][i].z;
                }

                if((k == mz-2) && !(mesh->k_periodic) && !(mesh->kk_periodic))
                {
                    fp[k+1][j][i].x = fp[k][j][i].x;
                    fp[k+1][j][i].y = fp[k][j][i].y;
                    fp[k+1][j][i].z = fp[k][j][i].z;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lFp, &fp);

    DMLocalToLocalBegin(fda, ueqn->lFp, INSERT_VALUES, ueqn->lFp);
    DMLocalToLocalEnd  (fda, ueqn->lFp, INSERT_VALUES, ueqn->lFp);

    resetCellPeriodicFluxes(mesh, ueqn->lFp, ueqn->lFp, "vector", "localToLocal");

    DMDAVecGetArray(fda, ueqn->lFp, &fp);

    // projection loop: fp (cell-centred combined flux) → face contravariant RHS.
    // Dead code removed: v0/v3 and d0-d3 cell-width computations were computed
    // but unused (only central(v1, v2) appears in the rhs assignment).
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal v1, v2;

                v1 = csi[k][j][i  ].x * fp[k][j][i  ].x + csi[k][j][i  ].y * fp[k][j][i  ].y + csi[k][j][i  ].z * fp[k][j][i  ].z;
                v2 = csi[k][j][i+1].x * fp[k][j][i+1].x + csi[k][j][i+1].y * fp[k][j][i+1].y + csi[k][j][i+1].z * fp[k][j][i+1].z;

                rhs[k][j][i].x += scale * iaj[k][j][i] * central(v1, v2);

                v1 = eta[k][j  ][i].x * fp[k][j  ][i].x + eta[k][j  ][i].y * fp[k][j  ][i].y + eta[k][j  ][i].z * fp[k][j  ][i].z;
                v2 = eta[k][j+1][i].x * fp[k][j+1][i].x + eta[k][j+1][i].y * fp[k][j+1][i].y + eta[k][j+1][i].z * fp[k][j+1][i].z;

                rhs[k][j][i].y += scale * jaj[k][j][i] * central(v1, v2);

                v1 = zet[k  ][j][i].x * fp[k  ][j][i].x + zet[k  ][j][i].y * fp[k  ][j][i].y + zet[k  ][j][i].z * fp[k  ][j][i].z;
                v2 = zet[k+1][j][i].x * fp[k+1][j][i].x + zet[k+1][j][i].y * fp[k+1][j][i].y + zet[k+1][j][i].z * fp[k+1][j][i].z;

                rhs[k][j][i].z += scale * kaj[k][j][i] * central(v1, v2);
            }
        }
    }

    if(ueqn->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);

    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);

    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    DMDAVecRestoreArray(fda, ueqn->lFp, &fp);
    DMDAVecRestoreArray(fda, ueqn->lDiv1, &div1);
    DMDAVecRestoreArray(fda, ueqn->lDiv2, &div2);
    DMDAVecRestoreArray(fda, ueqn->lDiv3, &div3);

    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, Rhs,  &rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode IMEXMatVec(Mat A, Vec v, Vec Av)
{

    // IMEX MatShell operator: A*v = v - dt*scale*Visc(v)
    // Used by kspIMEX to solve the linear IMEX system A*U^{n+1} = U^n + dt*bU.
    // Called by PETSc KSP for every matrix-vector product during the iterative solve

    ueqn_  *ueqn;
    MatShellGetContext(A, (void**)&ueqn);

    mesh_  *mesh    = ueqn->access->mesh;
    clock_ *clock   = ueqn->access->clock;
    PetscReal dt    = clock->dt;
    PetscReal scale;

    if (clock->it > clock->itStart)
    {
        scale = 0.5;   // Crank-Nicolson implicit half for subsequent steps
    }
    else
    {
        scale = 1.0;   // full backward-Euler implicit viscosity for first step to initialize the IMEX scheme
    }

    // sync v → Ucont/lUcont so FormU has the correct state
    VecCopy(v, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    
    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset contravariant periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");
    
    // transform to cartesian
    contravariantToCartesian(ueqn);

    // if(ueqn->access->flags->isIBMActive)
    // {
    //     UpdateImmersedBCs(ueqn->access->ibm);
    // }

    if(ueqn->access->flags->isOversetActive)
    {
        setBackgroundBC(mesh);
    }

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // compute scale*Visc(v) into scratch buffer ueqn->Rhs
    VecSet(ueqn->Rhs, 0.0);
    FormU(ueqn, ueqn->Rhs, scale, 2);   // visc-only: Rhs = scale*Visc(v)
    VecScale(ueqn->Rhs, dt);            // Rhs = dt*scale*Visc(v)
    resetNonResolvedCellFaces(mesh, ueqn->Rhs);

    // Av = v - dt*scale*Visc(v)
    VecCopy(v, Av);
    VecAXPY(Av, -1.0, ueqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnSNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
{
    ueqn_  *ueqn  = (ueqn_*)ptr;
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    VecCopy(Ucont, ueqn->Ucont);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset contravariant periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // transform to cartesian
    contravariantToCartesian(ueqn);

    // if(ueqn->access->flags->isIBMActive)
    // {
    //     UpdateImmersedBCs(ueqn->access->ibm);
    // }

    if(ueqn->access->flags->isOversetActive)
    {
        setBackgroundBC(mesh);
    }

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // initialize Rhs from the prebuilt constant-RHS buffer (bU).
    // bU = -dP [± buoyancy] [+ 0.5*Rhs_o (it>itStart)] [+ farm source]
    // Built once per time step in SolveUEqn, avoiding the equivalent VecSet + multiple
    // VecAXPY calls on every SNES function evaluation (called O(iter × ksp) times).
    VecCopy(ueqn->bU, Rhs);

    // get time step
    PetscReal dt    = clock->dt;
    PetscReal scale = (clock->it > clock->itStart) ? 0.5 : 1.0;

    // convection and diffusion terms 
    FormU(ueqn, Rhs, scale);

    // add coriolis term
    if(ueqn->access->flags->isAblActive && ueqn->access->abl->coriolisActive)
        Coriolis(ueqn, Rhs, 1.0);

    // add canopy term
    if(ueqn->access->flags->isCanopyActive)
        CanopyForce(ueqn, Rhs, 1.0);

    // add source terms 
    if(ueqn->access->flags->isAblActive && ueqn->access->abl->controllerActive)
        sourceU(ueqn, Rhs, 1.0);

    // add damping terms 
    if
    (
        ueqn->access->flags->isXDampingActive ||
        ueqn->access->flags->isYDampingActive ||
        ueqn->access->flags->isZDampingActive ||
        ueqn->access->flags->isKLeftRayleighDampingActive ||
        ueqn->access->flags->isKRightRayleighDampingActive
    )
    {
        dampingSourceU(ueqn, Rhs, 1.0);
    }

    // add mean pressure gradient forcing term
    if(ueqn->access->flags->isMeangradPForcingActive)
    {
        meanGradPForcing(ueqn, Rhs, 1.0);
    }

    // multiply for dt
    VecScale(Rhs, dt);
    
    resetNonResolvedCellFaces(mesh, Rhs);

    // add time derivative term
    VecAXPY(Rhs, -1.0, Ucont);
    VecAXPY(Rhs,  1.0, ueqn->Ucont_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsU(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset contravariant periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // transform to cartesian
    contravariantToCartesian(ueqn);

    // if(ueqn->access->flags->isIBMActive)
    // {
    //     UpdateImmersedBCs(ueqn->access->ibm);
    // }

    if(ueqn->access->flags->isOversetActive)
    {
        setBackgroundBC(mesh);;
    }

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // initialize the rhs vector
    VecSet(ueqn->Rhs, 0.0);

    // add pressure gradient term
    VecAXPY(ueqn->Rhs, -1.0, ueqn->dP);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->coriolisActive)
        {
            Coriolis(ueqn, ueqn->Rhs, 1.0);
        }
    }

    // add side force term
    if(ueqn->access->flags->isCanopyActive)
    {
        CanopyForce(ueqn, ueqn->Rhs, 1.0);
    }

    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;

        // add buoyancy gradient term
        if(teqn->pTildeFormulation)
        {
            VecAXPY(ueqn->Rhs, -1.0, teqn->ghGradRhok);
        }
        else
        {
            VecAXPY(ueqn->Rhs,  1.0, ueqn->bTheta);
        }
    }

    // add convection and diffusion terms
    FormU(ueqn, ueqn->Rhs, 1.0);

    // add wind farm model
    if(ueqn->access->flags->isWindFarmActive)
    {
        VecAXPY(ueqn->Rhs, 1.0, ueqn->access->farm->sourceFarmCont);
    }

    // add damping terms
    if
    (
        ueqn->access->flags->isXDampingActive ||
        ueqn->access->flags->isYDampingActive ||
        ueqn->access->flags->isZDampingActive ||
        ueqn->access->flags->isKLeftRayleighDampingActive ||
        ueqn->access->flags->isKRightRayleighDampingActive
    )
    {
        dampingSourceU(ueqn, ueqn->Rhs, 1.0);
    }

    // add driving source terms
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->controllerActive)
        {
            sourceU(ueqn, ueqn->Rhs, 1.0);
        }
    }

    if(ueqn->access->flags->isMeangradPForcingActive)
    {
        meanGradPForcing(ueqn, ueqn->Rhs, 1.0);
    }

    resetNonResolvedCellFaces(mesh, ueqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnEuler(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    // compute the right hand side
    FormExplicitRhsU(ueqn);

    // multiply for dt
    VecScale(ueqn->Rhs, clock->dt);

    // add time derivative term
    VecAXPY(ueqn->Rhs, 1, ueqn->Ucont);

    // store the solution
    VecCopy(ueqn->Rhs,  ueqn->Ucont);

    // scatter to local values
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnRK4(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RungeKutta-4: Solving for Ucont, Stage ");

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

    // Ucont_o contribution
    VecCopy(ueqn->Ucont_o, ueqn->Utmp);

    // contribution from K2, K3, K4
    for (PetscInt i=0; i<s; i++)
    {
        PetscPrintf(mesh->MESH_COMM, "%ld, ", i+1);

        // compute intermediate U guess and evaluate RHS
        if(i!=0)
        {
            VecWAXPY(ueqn->Ucont, a[i] * dt, ueqn->Rhs, ueqn->Ucont_o);
            DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
            DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
        }

        // compute function guess
        FormExplicitRhsU(ueqn);

        // add contribution from K1, K2, K3, K4
        VecAXPY(ueqn->Utmp, dt * b[i], ueqn->Rhs);
    }

    VecCopy(ueqn->Utmp, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode hyperViscosityU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale)
{
    if(ueqn->hyperVisc4i == 0.0 && ueqn->hyperVisc4j == 0.0 && ueqn->hyperVisc4k == 0.0) return(0);

    // Explicit biharmonic (4th-order, index-space) hyperviscosity for the IMEX scheme.
    // Each index direction has an independent coefficient, allowing e.g. horizontal-only
    // damping (i+j) without touching vertical (k) turbulence transport.

    mesh_        *mesh  = ueqn->access->mesh;
    clock_       *clock = ueqn->access->clock;
    DM            fda   = mesh->fda;
    DM            da    = mesh->da;
    DMDALocalInfo info  = mesh->info;
    PetscInt xs = info.xs, xe = info.xs + info.xm;
    PetscInt ys = info.ys, ye = info.ys + info.ym;
    PetscInt zs = info.zs, ze = info.zs + info.zm;
    PetscInt mx = info.mx, my = info.my, mz = info.mz;

    PetscInt lxs = xs; if(lxs == 0) lxs=lxs+1; PetscInt lxe = xe; if(lxe == mx) lxe=lxe-1;
    PetscInt lys = ys; if(lys == 0) lys=lys+1; PetscInt lye = ye; if(lye == my) lye=lye-1;
    PetscInt lzs = zs; if(lzs == 0) lzs=lzs+1; PetscInt lze = ze; if(lze == mz) lze=lze-1;

    PetscInt    i, j, k;
    Cmpnts    ***ucont, ***rhs;
    PetscReal ***nvert;
    PetscScalar solid = 0.5;

    // eps_{i,j,k}/dt: after VecScale(Rhs, dt) the net effect is
    //   Rhs += scale*dt * ( eps_i*δ⁴_i(U) + eps_j*δ⁴_j(U) + eps_k*δ⁴_k(U) )
    PetscReal eps_i = scale * ueqn->hyperVisc4i / clock->dt;
    PetscReal eps_j = scale * ueqn->hyperVisc4j / clock->dt;
    PetscReal eps_k = scale * ueqn->hyperVisc4k / clock->dt;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, Rhs,          &rhs);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    for(k = lzs; k < lze; k++)
    for(j = lys; j < lye; j++)
    for(i = lxs; i < lxe; i++)
    {
        PetscInt im2 = i-2, im1 = i-1, ip1 = i+1, ip2 = i+2;
        if(mesh->i_periodic)
        {
            if     (i == 1)    { im2 = mx-3; }               // im1=0 is already correct
            else if(i == mx-2) { ip2 = 2; }                  // ip1=mx-1 (ghost) is fine
        }
        else if(mesh->ii_periodic)
        {
            if     (i == 1)    { im1 = -2;   im2 = -3; }
            else if(i == 2)    { im2 = -2; }                  // 0 → -2 for ii_periodic
            else if(i == mx-3) { ip1 = mx-2; ip2 = mx+1; }
            else if(i == mx-2) { ip1 = mx+1; ip2 = mx+2; }   // mx-1→mx+1, mx→mx+2
        }
        else
        {
            // non-periodic: pad with the nearest interior value (zero-gradient / slip)
            if(im2 < 1)    im2 = 1;
            if(im1 < 1)    im1 = 1;
            if(ip1 > mx-2) ip1 = mx-2;
            if(ip2 > mx-2) ip2 = mx-2;
        }

        PetscInt jm2 = j-2, jm1 = j-1, jp1 = j+1, jp2 = j+2;
        if(mesh->j_periodic)
        {
            if     (j == 1)    { jm2 = my-3; }
            else if(j == my-2) { jp2 = 2; }
        }
        else if(mesh->jj_periodic)
        {
            if     (j == 1)    { jm1 = -2;   jm2 = -3; }
            else if(j == 2)    { jm2 = -2; }
            else if(j == my-3) { jp1 = my-2; jp2 = my+1; }
            else if(j == my-2) { jp1 = my+1; jp2 = my+2; }
        }
        else
        {
            if(jm2 < 1)    jm2 = 1;
            if(jm1 < 1)    jm1 = 1;
            if(jp1 > my-2) jp1 = my-2;
            if(jp2 > my-2) jp2 = my-2;
        }

        PetscInt km2 = k-2, km1 = k-1, kp1 = k+1, kp2 = k+2;
        if(mesh->k_periodic)
        {
            if     (k == 1)    { km2 = mz-3; }
            else if(k == mz-2) { kp2 = 2; }
        }
        else if(mesh->kk_periodic)
        {
            if     (k == 1)    { km1 = -2;   km2 = -3; }
            else if(k == 2)    { km2 = -2; }
            else if(k == mz-3) { kp1 = mz-2; kp2 = mz+1; }
            else if(k == mz-2) { kp1 = mz+1; kp2 = mz+2; }
        }
        else
        {
            if(km2 < 1)    km2 = 1;
            if(km1 < 1)    km1 = 1;
            if(kp1 > mz-2) kp1 = mz-2;
            if(kp2 > mz-2) kp2 = mz-2;
        }

        // i-face flux (ucont.x)
        if(nvert[k][j][i] + nvert[k][j][i+1] < 2.0 * solid)
        {
            rhs[k][j][i].x -=
                  eps_i * (ucont[k  ][j  ][ip2].x - 4.0*ucont[k  ][j  ][ip1].x + 6.0*ucont[k][j][i].x - 4.0*ucont[k  ][j  ][im1].x + ucont[k  ][j  ][im2].x)
                + eps_j * (ucont[k  ][jp2][i  ].x - 4.0*ucont[k  ][jp1][i  ].x + 6.0*ucont[k][j][i].x - 4.0*ucont[k  ][jm1][i  ].x + ucont[k  ][jm2][i  ].x)
                + eps_k * (ucont[kp2][j  ][i  ].x - 4.0*ucont[kp1][j  ][i  ].x + 6.0*ucont[k][j][i].x - 4.0*ucont[km1][j  ][i  ].x + ucont[km2][j  ][i  ].x);
        }

        // j-face flux (ucont.y)
        if(nvert[k][j][i] + nvert[k][j+1][i] < 2.0 * solid)
        {
            rhs[k][j][i].y -=
                  eps_i * (ucont[k  ][j  ][ip2].y - 4.0*ucont[k  ][j  ][ip1].y + 6.0*ucont[k][j][i].y - 4.0*ucont[k  ][j  ][im1].y + ucont[k  ][j  ][im2].y)
                + eps_j * (ucont[k  ][jp2][i  ].y - 4.0*ucont[k  ][jp1][i  ].y + 6.0*ucont[k][j][i].y - 4.0*ucont[k  ][jm1][i  ].y + ucont[k  ][jm2][i  ].y)
                + eps_k * (ucont[kp2][j  ][i  ].y - 4.0*ucont[kp1][j  ][i  ].y + 6.0*ucont[k][j][i].y - 4.0*ucont[km1][j  ][i  ].y + ucont[km2][j  ][i  ].y);
        }

        // k-face flux
        if(nvert[k][j][i] + nvert[k+1][j][i] < 2.0 * solid)
        {
            rhs[k][j][i].z -=
                  eps_i * (ucont[k  ][j  ][ip2].z - 4.0*ucont[k  ][j  ][ip1].z + 6.0*ucont[k][j][i].z - 4.0*ucont[k  ][j  ][im1].z + ucont[k  ][j  ][im2].z)
                + eps_j * (ucont[k  ][jp2][i  ].z - 4.0*ucont[k  ][jp1][i  ].z + 6.0*ucont[k][j][i].z - 4.0*ucont[k  ][jm1][i  ].z + ucont[k  ][jm2][i  ].z)
                + eps_k * (ucont[kp2][j  ][i  ].z - 4.0*ucont[kp1][j  ][i  ].z + 6.0*ucont[k][j][i].z - 4.0*ucont[km1][j  ][i  ].z + ucont[km2][j  ][i  ].z);
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, Rhs,          &rhs);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);

    return(0);
}

//***************************************************************************************************************//
// IMEX-CNAB time-step routine.  Call order:
//   1. Evaluate Conv^n and Visc^n from U^n (current Ucont at step start).
//   2. Build the explicit RHS buffer bU once.  bU holds ALL terms that do NOT
//      depend on U^{n+1}: pressure, buoyancy, AB2 convection, CN explicit
//      viscous half, farm source, Coriolis, sourceU, damping, meanGradP.
//   3. Build the KSP RHS: Utmp = U^n + dt * bU (sanitized at boundary faces).
//   4. KSPSolve: A*U^{n+1} = Utmp, where A*v = v - dt*scale*Visc(v).
//      A is linear; no SNES needed.  Converges in O(1) KSP iterations.
//   5. Rotate RhsConv storage for the next step's AB2 stencil.
//***************************************************************************************************************//

PetscErrorCode UeqnIMEX(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts, te;
    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "IMEX-CNAB: Solving for Ucont, ");

    // Step 1: prepare U^n state (synchronise lUcont / lUcat) 

    resetNoPenetrationFluxes(ueqn);
    resetFacePeriodicFluxesVector(mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");
    contravariantToCartesian(ueqn);
    if(ueqn->access->flags->isOversetActive) setBackgroundBC(mesh);
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // Step 2a: convective RHS at U^n → RhsConv  
    VecSet(ueqn->RhsConv, 0.0);
    FormU(ueqn, ueqn->RhsConv, 1.0, 1);   // conv-only

    // Step 2b: viscous RHS at U^n (scratch in Rhs Vec)  
    VecSet(ueqn->Rhs, 0.0);
    FormU(ueqn, ueqn->Rhs, 1.0, 2);        // visc-only

    // Step 3: build bU (explicit part, evaluated once at U^n) 
    //
    //   bU = -dP ± buoy  +  AB2_conv  +  CN_visc_explicit  +  farm
    //         + Coriolis  +  CanopyForce  +  sourceU  +  damping  +  meanGradP
    //

    // pressure gradient
    VecCopy(ueqn->dP, ueqn->bU);
    VecScale(ueqn->bU, -1.0);                           

    // buoyancy — AB2: b^{n+1/2} ≈ 1.5*b^n - 0.5*b^{n-1} (Forward Euler on first step)
    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;
        if(teqn->pTildeFormulation)
        {
            if(clock->it > clock->itStart)
            {
                VecAXPY(ueqn->bU, -1.5, teqn->ghGradRhok);
                VecAXPY(ueqn->bU,  0.5, teqn->ghGradRhok_o);
            }
            else
            {
                VecAXPY(ueqn->bU, -1.0, teqn->ghGradRhok);
            }
        }
        else
        {
            if(clock->it > clock->itStart)
            {
                // AB2 extrapolation: 1.5*b(T^n) - 0.5*b(T^{n-1})
                VecAXPY(ueqn->bU,  1.5, ueqn->bTheta);
                VecAXPY(ueqn->bU, -0.5, ueqn->bTheta_o);
            }
            else
            {
                VecAXPY(ueqn->bU, 1.0, ueqn->bTheta);
            }
        }
    }

    // wind farm source
    if(ueqn->access->flags->isWindFarmActive)
        VecAXPY(ueqn->bU, 1.0, ueqn->access->farm->sourceFarmCont);

    // AB2 convection extrapolation: c1=3/2, c2=-1/2 (Forward Euler on first step)
    if(clock->it > clock->itStart)
    {
        VecAXPY(ueqn->bU,  1.5, ueqn->RhsConv);   // (3/2) * Conv^n
        VecAXPY(ueqn->bU, -0.5, ueqn->RhsConv_o); // (-1/2) * Conv^{n-1}
    }
    else
    {
        VecAXPY(ueqn->bU, 1.0, ueqn->RhsConv);    // Forward Euler on first step
    }

    // Crank-Nicolson viscous: explicit half 0.5*Visc^n
    // first step: visc is fully implicit (backward Euler), nothing added here (IMEXMatVec uses scale=1.0 on that step)
    if(clock->it > clock->itStart)
        VecAXPY(ueqn->bU, 0.5, ueqn->Rhs);            // +0.5 * Visc^n

    // Coriolis
    if(ueqn->access->flags->isAblActive)
        if(ueqn->access->abl->coriolisActive)
            Coriolis(ueqn, ueqn->bU, 1.0);

    // canopy force
    if(ueqn->access->flags->isCanopyActive)
        CanopyForce(ueqn, ueqn->bU, 1.0);

    // ABL driving source (PI controller)
    if(ueqn->access->flags->isAblActive)
        if(ueqn->access->abl->controllerActive)
            sourceU(ueqn, ueqn->bU, 1.0);

    // sponge/Rayleigh damping
    if(ueqn->access->flags->isXDampingActive ||
       ueqn->access->flags->isYDampingActive ||
       ueqn->access->flags->isZDampingActive ||
       ueqn->access->flags->isKLeftRayleighDampingActive ||
       ueqn->access->flags->isKRightRayleighDampingActive)
    {
        dampingSourceU(ueqn, ueqn->bU, 1.0);
    }

    // mean-gradient pressure forcing
    if(ueqn->access->flags->isMeangradPForcingActive)
        meanGradPForcing(ueqn, ueqn->bU, 1.0);

    // biharmonic hyperviscosity: damp checkerboard / Nyquist modes.
    if(ueqn->hyperVisc4i > 0.0 || ueqn->hyperVisc4j > 0.0 || ueqn->hyperVisc4k > 0.0)
        hyperViscosityU(ueqn, ueqn->bU, 1.0);

    // Step 4: initial guess = explicit Euler estimate of U^{n+1} 
    //         b = Utmp = U^n + dt * bU  (also serves as non-zero initial guess)
    VecCopy(ueqn->Ucont_o, ueqn->Utmp);
    VecAXPY(ueqn->Utmp, clock->dt, ueqn->bU);
    
    // zero boundary/IBM faces before KSPSolve
    resetNonResolvedCellFaces(mesh, ueqn->Utmp);

    // Step 5: KSP direct solve: A*Ucont = Utmp  (A is linear, no SNES needed) ──
    //   RhsConv_o is used as the KSP solution buffer to avoid aliasing:
    //   IMEXMatVec writes to ueqn->Ucont/lUcont as FormU workspace, so the KSP
    //   solution vector must be a *different* Vec to prevent in-place update
    //   (VecAXPY on x) from interfering with the MatVec workspace.
    //   RhsConv_o is free at this point (its value was already consumed into bU).
    //   Step 7 will overwrite RhsConv_o with RhsConv, so the order matters.
    
    // initial guess 
    VecCopy(ueqn->Utmp, ueqn->RhsConv_o);    

    // solve 
    KSPSolve(ueqn->kspIMEX, ueqn->Utmp, ueqn->RhsConv_o);

    PetscReal norm;  PetscInt iter;
    KSPGetResidualNorm(ueqn->kspIMEX, &norm);
    KSPGetIterationNumber(ueqn->kspIMEX, &iter);

    // Step 6: copy solution and scatter 
    VecCopy(ueqn->RhsConv_o, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // Step 7: rotate Conv buffer for AB2 at next step 
    VecCopy(ueqn->RhsConv, ueqn->RhsConv_o);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,
        "Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n",
        norm, iter, te-ts);

    return(0);
}

//***************************************************************************************************************//

