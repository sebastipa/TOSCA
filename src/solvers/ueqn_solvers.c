//! \file  ueqn_solvers.c
//! \brief Momentum-equation solver function implementations.

#include "ueqn_solvers.h"  

PetscErrorCode FormU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale, PetscInt formMode)
{
    // Diffusion + convection terms of the momentum right hand side at cell faces. 
    // The Rhs is staggered, i.e. components of each ikj are not at the same location.
    // 1. Diffusion and convection fluxes in cartesian coordinates are evaluated at 
    //    every cell face along each curvilinear directions.  
    // 2. divergence in cartesian coordinates is calculated at cell centers by cumulating fluxes on each cell face. 
    // 3. The cell centered flux in cartesian coordinates is interpolated at the faces and dotted with the face area vectors,
    //    for each curvilinear direction.

    // Notes: 
    // A: formMode: 0 = conv+visc (default), 1 = conv-only, 2 = visc-only
    // B: non-periodic boundaries: only internal faces are solved (0 and m-2 excluded)
    // C: periodic boundaries: contravariant velocity at the right boundary faces is solved for (m-2). 
    //    Since 0 is coincident with m-2 in this case, velocity at m-2 is copied to 0 after the solve. 
    //    right boundaries were the same boundary) 

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

    Cmpnts           ***div1,  ***div2,  ***div3;
    Cmpnts           ***visc1, ***visc2, ***visc3;
    Cmpnts           ***viscIBM1, ***viscIBM2, ***viscIBM3;
    Cmpnts           ***rhs,   ***fp,    ***limiter;
    PetscReal        ***aj,    ***iaj,   ***jaj, ***kaj;

    PetscReal        dudc, dude, dudz, dvdc, dvde, dvdz, dwdc, dwde, dwdz;
    PetscReal        csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;
    PetscReal        g11, g21, g31;
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
    // FORM DIVERGENCE AND VISCOUS FLUXES AT CELL FACES                       //
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
                    } 

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
                    }
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
                    } 

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
                    } 
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
                    } 

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
                    } 
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

    // add LES strucural model contribution in contravariant form
    if(formMode != 1 && ueqn->access->flags->isLesActive && (les->model == BAMD || les->model == BV))
    {
        updateLESStructuralModelContravariantForm(les); 
    }

    // ---------------------------------------------------------------------- //
    // FORM THE RIGHT HAND SIDE AT CELL CENTERS                               //
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
                } 

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

                // handle non periodic boundaries
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

    // ---------------------------------------------------------------------- //
    // INTERPOLATE AT CELL FACES AND DOT WITH FACE AREA VECTORS               //
    // ---------------------------------------------------------------------- //

    DMDAVecGetArray(fda, ueqn->lFp, &fp);

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

PetscErrorCode UeqnSNES(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal     norm;
    PetscInt      iter;
    PetscReal     ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "CN%s: Solving for Ucont, Initial residual = ", ueqn->solverType.c_str());

    // build the constant-RHS buffer for this time step

    // 1. pressure gradient
    VecCopy(ueqn->dP, ueqn->bU);
    VecScale(ueqn->bU, -1.0);

    // 2. buoyancy (explicit for first step, AB2 for subsequent steps)
    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;
        
        if(teqn->pTildeFormulation)
        {
            if(clock->it > clock->itStart)
            {
                VecAXPY(ueqn->bU, -3.0/2.0, teqn->ghGradRhok);
                VecAXPY(ueqn->bU,  1.0/2.0, teqn->ghGradRhok_o);
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
                VecAXPY(ueqn->bU, 3.0/2.0, ueqn->bTheta);
                VecAXPY(ueqn->bU, -1.0/2.0, ueqn->bTheta_o);
            }
            else
            {
                VecAXPY(ueqn->bU, 1.0, ueqn->bTheta);
            }
        }
    }

    // 3. explicit filtering 
    if
    (
        ueqn->hyperVisc4i > 0.0 || 
        ueqn->hyperVisc4j > 0.0 || 
        ueqn->hyperVisc4k > 0.0
    )
    {
        hyperViscosityU(ueqn, ueqn->bU, 1.0);
    }
    
    // 4. half Rhs for Crank-Nicolson time stepping 
    if(clock->it > clock->itStart)
        VecAXPY(ueqn->bU, 0.5, ueqn->Rhs_o);

    // 5. wind farm source
    if(ueqn->access->flags->isWindFarmActive)
        VecAXPY(ueqn->bU, 1.0, ueqn->access->farm->sourceFarmCont);

    // 6. source terms 
    if(ueqn->access->flags->isAblActive)
        if(ueqn->access->abl->controllerActive)
            sourceU(ueqn, ueqn->bU, 1.0);
    
    // initial guess is old velocity (might change that in the future but it doesnt make a big difference)
    VecCopy(ueqn->Ucont_o, ueqn->Utmp);

    // solve momentum equation 
    SNESSolve(ueqn->snesU, PETSC_NULL, ueqn->Utmp);
    SNESGetFunctionNorm(ueqn->snesU, &norm);
    SNESGetIterationNumber(ueqn->snesU, &iter);

    // report total inner linear iterations (= total MatMFFD FormU calls for J*v)
    PetscInt linIter = 0;
    SNESGetLinearSolveIterations(ueqn->snesU, &linIter);

    // store the solution and global to local scatter
    VecCopy(ueqn->Utmp, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Final residual = %e, Iterations = %ld (linear = %ld), Elapsed Time = %lf\n", norm, iter, linIter, te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SNESFuncEval(SNES snes, Vec Ucont, Vec Rhs, void *ptr)
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
    VecCopy(ueqn->bU, Rhs);

    // scaling factotor (scale is 1 for the first time step (BDF1), 0.5 for subsequent steps CN)
    PetscReal scale = (clock->it > clock->itStart) ? 0.5 : 1.0;

    // implicit convection and diffusion terms 
    FormU(ueqn, Rhs, scale);

    // implicit Coriolis term
    if(ueqn->access->flags->isAblActive && ueqn->access->abl->coriolisActive)
        Coriolis(ueqn, Rhs, 1.0);

    // implicit canopy term
    if(ueqn->access->flags->isCanopyActive)
        CanopyForce(ueqn, Rhs, 1.0);

    // implicit damping terms 
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

    // implicit mean pressure gradient forcing term
    if(ueqn->access->flags->isMeangradPForcingActive)
    {
        meanGradPForcing(ueqn, Rhs, 1.0);
    } 

    // multiply for dt
    VecScale(Rhs, clock->dt);
    
    // set the Rhs to zero at non-resolved cell faces (override numerical errors)
    resetNonResolvedCellFaces(mesh, Rhs);

    // add time derivative term
    VecAXPY(Rhs, -1.0, Ucont);
    VecAXPY(Rhs,  1.0, ueqn->Ucont_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ExplicitRhsU(ueqn_ *ueqn)
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

    // add buoyancy gradient term
    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;
        if(teqn->pTildeFormulation)
        {
            if(clock->it > clock->itStart)
            {
                VecAXPY(ueqn->Rhs, -3.0/2.0, teqn->ghGradRhok);
                VecAXPY(ueqn->Rhs,  1.0/2.0, teqn->ghGradRhok_o);
            }
            else
            {
                VecAXPY(ueqn->Rhs, -1.0, teqn->ghGradRhok);
            }
        }
        else
        {
            if(clock->it > clock->itStart)
            {
                VecAXPY(ueqn->Rhs, 3.0/2.0, ueqn->bTheta);
                VecAXPY(ueqn->Rhs, -1.0/2.0, ueqn->bTheta_o);
            }
            else
            {
                VecAXPY(ueqn->Rhs, 1.0, ueqn->bTheta);
            }
        }
    }

    // add wind farm model
    if(ueqn->access->flags->isWindFarmActive)
    {
        VecAXPY(ueqn->Rhs, 1.0, ueqn->access->farm->sourceFarmCont);
    }

    // add convection and diffusion terms
    FormU(ueqn, ueqn->Rhs, 1.0);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->coriolisActive)
        {
            Coriolis(ueqn, ueqn->Rhs, 1.0);
        }
    }

    // canopy force
    if(ueqn->access->flags->isCanopyActive)
    {
        CanopyForce(ueqn, ueqn->Rhs, 1.0);
    }

    // add driving source terms
    if(ueqn->access->flags->isAblActive)
    {
        if(ueqn->access->abl->controllerActive)
        {
            sourceU(ueqn, ueqn->Rhs, 1.0);
        }
    }

    // sponge/Rayleigh damping
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

    // mean-gradient pressure forcing
    if(ueqn->access->flags->isMeangradPForcingActive)
    {
        meanGradPForcing(ueqn, ueqn->Rhs, 1.0);
    }

    // biharmonic hyperviscosity: damp checkerboard / Nyquist modes.
    if(ueqn->hyperVisc4i > 0.0 || ueqn->hyperVisc4j > 0.0 || ueqn->hyperVisc4k > 0.0)
        hyperViscosityU(ueqn, ueqn->Rhs, 1.0);

    // set the Rhs to zero at non-resolved cell faces (override numerical errors)
    resetNonResolvedCellFaces(mesh, ueqn->Rhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnFE(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "FE%s: Solving for Ucont, ", ueqn->solverType.c_str());

    // compute the right hand side
    ExplicitRhsU(ueqn);

    // multiply by dt
    VecScale(ueqn->Rhs, clock->dt);

    // add time derivative term
    VecAXPY(ueqn->Rhs, 1, ueqn->Ucont_o);

    // store the solution
    VecCopy(ueqn->Rhs,  ueqn->Ucont);

    // scatter to local values
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %f\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UeqnRK4(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts,te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "RK4%s: Solving for Ucont, Stage ",ueqn->solverType.c_str());

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
        ExplicitRhsU(ueqn);

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

PetscErrorCode UeqnCNAB(ueqn_ *ueqn)
{
    mesh_  *mesh  = ueqn->access->mesh;
    clock_ *clock = ueqn->access->clock;

    PetscReal ts, te;
    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "CNAB%s: Solving for Ucont, ",ueqn->solverType.c_str());

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

    // set convective and viscous RHS at U^n then build the right and side
    VecSet(ueqn->RhsConv, 0.0); FormU(ueqn, ueqn->RhsConv, 1.0, 1);
    VecSet(ueqn->RhsVisc, 0.0); FormU(ueqn, ueqn->RhsVisc, 1.0, 2);

    // pressure gradient
    VecCopy(ueqn->dP, ueqn->bU);
    VecScale(ueqn->bU, -1.0);

    // add buoyancy gradient term
    if(ueqn->access->flags->isTeqnActive)
    {
        teqn_ *teqn = ueqn->access->teqn;
        if(teqn->pTildeFormulation)
        {
            if(clock->it > clock->itStart)
            {
                VecAXPY(ueqn->bU, -3.0/2.0, teqn->ghGradRhok);
                VecAXPY(ueqn->bU,  1.0/2.0, teqn->ghGradRhok_o);
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
                VecAXPY(ueqn->bU, 3.0/2.0, ueqn->bTheta);
                VecAXPY(ueqn->bU, -1.0/2.0, ueqn->bTheta_o);
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

    
    // add convection
    if(clock->it > clock->itStart)
    {
        // AB2: c1=3/2, c2=-1/2 
        VecAXPY(ueqn->bU,  1.5, ueqn->RhsConv);   
        VecAXPY(ueqn->bU, -0.5, ueqn->RhsConv_o);
    }
    else
    {
        // forward Euler on first step 
        VecAXPY(ueqn->bU, 1.0, ueqn->RhsConv);   
    }

    // add Crank-Nicolson viscous half (at the start this is fully implicit)
    if(clock->it > clock->itStart)
        VecAXPY(ueqn->bU, 0.5, ueqn->RhsVisc);

    // add coriolis term
    if(ueqn->access->flags->isAblActive)
        if(ueqn->access->abl->coriolisActive)
            Coriolis(ueqn, ueqn->bU, 1.0);

    // add canopy force
    if(ueqn->access->flags->isCanopyActive)
        CanopyForce(ueqn, ueqn->bU, 1.0);

    // ABL driving source 
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
    if
    (
        ueqn->hyperVisc4i > 0.0 || 
        ueqn->hyperVisc4j > 0.0 || 
        ueqn->hyperVisc4k > 0.0
    )
    {
        hyperViscosityU(ueqn, ueqn->bU, 1.0);
    }

    // set the Rhs to zero at non-resolved cell faces (override numerical errors)
    resetNonResolvedCellFaces(mesh, ueqn->bU);

    // multiply by dt
    VecScale(ueqn->bU, clock->dt);

    // add U^n contribution: bU = U^n + dt * bU
    VecAXPY(ueqn->bU, 1.0, ueqn->Ucont_o);

    // form initial guess as explicit Euler estimate
    VecCopy(ueqn->bU, ueqn->Utmp);

    // add missing viscous hal:  Utmp = U^n + dt*bU + 0.5*dt*Visc^n 
    // note: this is essential to get good convergence from the KSP solver
    VecAXPY(ueqn->Utmp, 0.5, ueqn->RhsVisc);

    // solve linear system
    KSPSolve(ueqn->kspIMEX, ueqn->bU, ueqn->Utmp);

    PetscReal norm;  PetscInt iter;
    KSPGetResidualNorm(ueqn->kspIMEX, &norm);
    KSPGetIterationNumber(ueqn->kspIMEX, &iter);

    // Step 6: copy solution and scatter 
    VecCopy(ueqn->Utmp, ueqn->Ucont);
    DMGlobalToLocalBegin(mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (mesh->fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // update the cartesian velocity (this is actually required, as there is no guarantee that MatMult is called on the last Ucont)
    contravariantToCartesian(ueqn);

    // reset cartesian periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ueqn->Ucat, ueqn->lUcat, "vector", "globalToLocal");

    // rotate Conv buffer for AB2 at next step.
    VecCopy(ueqn->RhsConv, ueqn->RhsConv_o);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM, "Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n", norm, iter, te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CNABMatVec(Mat A, Vec v, Vec Av)
{

    // IMEX MatShell operator: A*v = v - dt*scale*Visc(v)
    // Called by kspIMEX to solve the linear IMEX system A*U^{n+1} = U^n + dt*bU for every matrix-vector product during the iterative solve

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

    VecCopy(v, ueqn->Ucont);

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

    // visc-only: Rhs = scale*Visc(v)
    FormU(ueqn, Av, -1.0*scale, 2);

    // Rhs = dt*scale*Visc(v)
    VecScale(Av, dt);   

    resetNonResolvedCellFaces(mesh, Av);

    // Av = v - dt*scale*Visc(v)
    VecAXPY(Av, 1.0, v);

    return(0);
}

//***************************************************************************************************************//

