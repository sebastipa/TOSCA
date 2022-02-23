//! \file  teqn.c
//! \brief Contains T equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorT(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeTEqn(teqn_ *teqn)
{
    if(teqn != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = teqn->access->mesh;

        //default parameters
        teqn->absExitTol        = 1e-5;
        teqn->relExitTol        = 1e-30;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolT",  &(teqn->absExitTol), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolT",  &(teqn->relExitTol), PETSC_NULL);

        VecDuplicate(mesh->Nvert, &(teqn->TmprtTmp));    VecSet(teqn->TmprtTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt));       VecSet(teqn->Tmprt,    0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Tmprt_o));     VecSet(teqn->Tmprt_o,  0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Rhs));         VecSet(teqn->Rhs,      0.0);
        VecDuplicate(mesh->Nvert, &(teqn->Rhs_o));       VecSet(teqn->Rhs_o,    0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt));      VecSet(teqn->lTmprt,   0.0);
        VecDuplicate(mesh->lAj,   &(teqn->lTmprt_o));    VecSet(teqn->lTmprt_o, 0.0);

        VecDuplicate(mesh->lCent, &(teqn->lDivT));       VecSet(teqn->lDivT,    0.0);
        VecDuplicate(mesh->lCent, &(teqn->lViscT));      VecSet(teqn->lViscT,   0.0);

        // create the SNES solver
        SNESCreate(mesh->MESH_COMM,&(teqn->snesT));
        SNESMonitorSet(teqn->snesT, SNESMonitorT, (void*)&(mesh->MESH_COMM), PETSC_NULL);

        // set the SNES evaluating function
        SNESSetFunction(teqn->snesT, teqn->Rhs, TeqnSNES, (void *)teqn);

        // create jacobian matrix
        MatCreateSNESMF(teqn->snesT, &(teqn->JT));
        SNESSetJacobian(teqn->snesT, teqn->JT, teqn->JT, MatMFFDComputeJacobian, (void *)teqn);

        // set SNES solver type
        //SNESSetType(teqn->snesT, SNESNEWTONTR);           //SNESTR
        SNESSetType(teqn->snesT, SNESNEWTONLS);        //SNESLS is better for stiff PDEs such as the one including IB but slower

        // set SNES solve and step failures
        SNESSetMaxLinearSolveFailures(teqn->snesT,10000);
        SNESSetMaxNonlinearStepFailures(teqn->snesT,10000);
        SNESKSPSetUseEW(teqn->snesT, PETSC_TRUE);

        // set SNES Krylov Sub-Space parameters
        SNESKSPSetParametersEW(teqn->snesT,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

        // SNES tolerances
        // 2nd arg: absolute tolerance
        // 3rd arg: relative tolerance
        // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
        // 5th arg: maximum number of iterations
        // 6th arg: maximum function evaluations
        SNESSetTolerances(teqn->snesT, teqn->absExitTol, 1e-30, 1e-30, 20, 1000);

        SNESGetKSP(teqn->snesT, &(teqn->ksp));
        KSPGetPC(teqn->ksp,&(teqn->pc));

        // set KSP solver type
        KSPSetType(teqn->ksp, KSPGMRES);

        //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);    //2009.09.22 poor performance
        //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);      //2009.09.22

        //KSPFischerGuess itg;
        //KSPFischerGuessCreate(ksp,1,100,&itg);
        //KSPSetFischerGuess(ksp, itg);                  //2009.09.22

        //KSPGMRESSetPreAllocateVectors(ksp);            --> crazy thing consumes memory

        PCSetType(teqn->pc, PCNONE);
        PetscReal rtol=teqn->relExitTol, atol=teqn->absExitTol, dtol=PETSC_DEFAULT;
        KSPSetTolerances(teqn->ksp, rtol, atol, dtol, 1000);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode dampingSourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
{
    abl_          *abl  = teqn->access->abl;
    mesh_         *mesh = teqn->access->mesh;
    ueqn_         *ueqn = teqn->access->ueqn;
    les_          *les  = teqn->access->les;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***cent;
    PetscReal     ***rhs, ***t, ***tP;

    precursor_    *precursor;
    domain_       *pdomain;
    PetscInt      kStart;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent,  &cent);
    DMDAVecGetArray(da,  teqn->lTmprt, &t);
    DMDAVecGetArray(da,  Rhs,  &rhs);

    if(teqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            precursor = abl->precursor;
            pdomain   = precursor->domain;

            if(precursor->thisProcessorInFringe)
            {
                DMDAVecGetArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
                kStart = precursor->map.kStart;
            }
        }
    }

    // x damping layer
    PetscReal alphaX = abl->xDampingAlpha;
    PetscReal xS     = abl->xDampingStart;
    PetscReal xE     = abl->xDampingEnd;
    PetscReal xD     = abl->xDampingDelta;

    // loop over internal cell faces
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(teqn->access->flags->isXDampingActive)
                {
                    // compute cell center x at i,j,k
                    PetscReal x     = cent[k][j][i].x;

                    // compute Nordstrom viscosity at i,j,k
                    PetscReal nud_x   = viscNordstrom(alphaX, xS, xE, xD, x);
                    PetscReal tBarInstX;

                    if(abl->xFringeUBarSelectionType == 1 || abl->xFringeUBarSelectionType == 2)
                    {
                        // set desired temperature
                        tBarInstX  = abl->tBarInstX[j][i];
                    }
                    else if(abl->xFringeUBarSelectionType == 3)
                    {
                        if(precursor->thisProcessorInFringe)
                        {
                            // set desired temperature
                            tBarInstX  = tP[k+kStart][j][i];
                        }
                        else
                        {
                            tBarInstX  = t[k][j][i];
                        }
                    }

                    rhs[k][j][i] += scale * nud_x * (tBarInstX - t[k][j][i]);
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  teqn->lTmprt, &t);
    DMDAVecRestoreArray(da,  Rhs,  &rhs);

    if(teqn->access->flags->isXDampingActive)
    {
        if(abl->xFringeUBarSelectionType == 3)
        {
            if(precursor->thisProcessorInFringe)
            {
                DMDAVecRestoreArray(pdomain->mesh->da, pdomain->teqn->lTmprt,  &tP);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
{
    // In this function the viscous + divergence term of the temperature equation are
    // discretized at cell centers.
    // First the divergence and viscous fluxes are evaluated at cell faces, then
    // their budget is evaluated at the internal cells, forming the Rhs.

    mesh_         *mesh  = teqn->access->mesh;
    ueqn_         *ueqn  = teqn->access->ueqn;
    les_          *les   = teqn->access->les;
    constants_    *cst   = teqn->access->constants;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***ucont;

    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;

    PetscReal     ***tmprt, ***rhs, ***nvert;

    Cmpnts        ***div, ***visc;                                                // divergence and viscous terms
    Cmpnts        ***limiter;                                                     // flux limiter
    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;                                  // cell and face jacobians
    PetscReal     ***lnu_t;

    PetscReal     dtdc, dtde, dtdz;                                              // velocity der. w.r.t. curvil. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  Rhs,              &rhs);
    DMDAVecGetArray(da,  teqn->lTmprt,      &tmprt);
    DMDAVecGetArray(fda, ueqn->lUcont,     &ucont);

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
    DMDAVecGetArray(da,  mesh->lNvert,      &nvert);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    VecSet(teqn->lDivT,  0.0);
    VecSet(teqn->lViscT, 0.0);;

    DMDAVecGetArray(fda, teqn->lDivT, &div);
    DMDAVecGetArray(fda, teqn->lViscT, &visc);

    if(teqn->access->flags->isLesActive)
    {
        DMDAVecGetArray(da, les->lNu_t, &lnu_t);
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
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(j==0 || k==0) continue;

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

                Compute_dscalar_i(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    iL, iR;
                PetscReal denom;
                getFace2Cell4StencilCsi(mesh, i, mx, &iL, &iR, &denom);

                if( isIBMIFace(k, j, iL, iR, nvert) )
                {
                    iL = i;
                    iR = i+1;
                }

                PetscReal ucon = ucont[k][j][i].x;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                //IBM active and on the right face of the IBM fluid cell
                // Use second order if positive flux and quick if negative flux
                if(teqn->access->flags->isIBMActive && i!=mx-2 && isIBMCell(k, j, i, nvert))

                {
                    div[k][j][i].x =
                            -um * (0.125 * (-    tmprt[k][j][i+2] - 2. * tmprt[k][j][i+1] + 3. * tmprt[k][j][i]) + tmprt[k][j][i+1]) +
                            -up * (0.125 * (-    tmprt[k][j][i] -  2. * tmprt[k][j][i] +  3. * tmprt[k][j][i+1]) +  tmprt[k][j][i]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if (teqn->access->flags->isIBMActive && i!=0 && isIBMCell(k, j, i+1, nvert))
                {
                    div[k][j][i].x =
                            -um * (0.125 * (-    tmprt[k][j][i+1] -  2. * tmprt[k][j][i+1] +  3. * tmprt[k][j][i]) + tmprt[k][j][i+1]) +
                            -up * (0.125 * (-    tmprt[k][j][i-1] -  2. * tmprt[k][j][i] +  3. * tmprt[k][j][i+1]) + tmprt[k][j][i]);

                }
                else
                {
                    div[k][j][i].x =
                    - ucon
                    * centralUpwind
                    (
                        tmprt[k][j][iL],
                        tmprt[k][j][i],
                        tmprt[k][j][i+1],
                        tmprt[k][j][iR],
                        ucont[k][j][i].x
                        ,limiter[k][j][i].x
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                // viscous terms
                if
                (
                    teqn->access->flags->isLesActive
                )
                {
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                    // compute stability dependent turbulent Prandtl number
                    PetscReal gradTdotG = dtde*(-9.81);
                    PetscReal l, delta = pow( 1./ajc, 1./3. );
                    if(gradTdotG < 0.)
                    {
                        l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                    }
                    else
                    {
                        l = delta;
                    }

                    PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));

                    kappaEff = (nu / cst->Pr) + (nut / Prt);
                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].x += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
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
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || k==0) continue;

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

                Compute_dscalar_j(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    jL, jR;
                PetscReal denom;
                getFace2Cell4StencilEta(mesh, j, my, &jL, &jR, &denom);

                if( isIBMJFace(k, jL, i, jR, nvert) )
                {
                    jL = j;
                    jR = j+1;
                }

                PetscReal ucon = ucont[k][j][i].y;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                //IBM active and on the right face of the IBM fluid cell
                // Use second order if positive flux and quick if negative flux
                if ( teqn->access->flags->isIBMActive &&  j!=my-2 && isIBMCell(k, j, i, nvert) )
                {
                    div[k][j][i].y =
                            -um * (0.125 * (-    tmprt[k][j+2][i] -  2. * tmprt[k][j+1][i] +  3. * tmprt[k][j  ][i]) +  tmprt[k][j+1][i]) +
                            -up * (0.125 * (-    tmprt[k][j  ][i] -  2. * tmprt[k][j  ][i] +  3. * tmprt[k][j+1][i]) +  tmprt[k][j][i  ]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if ( teqn->access->flags->isIBMActive &&  j!=0 && isIBMCell(k, j+1, i, nvert) )
                {
                    div[k][j][i].y =
                            -um * (0.125 * (-    tmprt[k][j+1][i] -  2. * tmprt[k][j+1][i] +  3. * tmprt[k][j  ][i]) + tmprt[k][j+1][i]) +
                            -up * (0.125 * (-    tmprt[k][j-1][i] -  2. * tmprt[k][j  ][i] +  3. * tmprt[k][j+1][i]) + tmprt[k][j][i  ]);

                }

                else
                {
                    div[k][j][i].y =
                    - ucon
                    * centralUpwind
                    (
                        tmprt[k][jL][i],
                        tmprt[k][j][i],
                        tmprt[k][j+1][i],
                        tmprt[k][jR][i],
                        ucont[k][j][i].y
                        ,limiter[k][j][i].y
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                if
                (
                    teqn->access->flags->isLesActive
                )
                {
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                    // compute stability depentend turbulent Prandtl number
                    PetscReal gradTdotG = dtde*(-9.81);
                    PetscReal l, delta = pow( 1./ajc, 1./3. );
                    if(gradTdotG < 0.)
                    {
                        l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                    }
                    else
                    {
                        l = delta;
                    }

                    PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));

                    kappaEff = (nu / cst->Pr) + (nut / Prt);
                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].y += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
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
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || j==0) continue;

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

                Compute_dscalar_k(mesh, i, j, k, mx, my, mz, tmprt, nvert, &dtdc, &dtde, &dtdz);

                PetscInt    kL, kR;
                PetscReal denom;
                getFace2Cell4StencilZet(mesh, k, mz, &kL, &kR, &denom);

                if( isIBMKFace(kL, j, i, kR, nvert) )
                {
                    kL = k;
                    kR=k+1;
                    denom=1.;
                }

                PetscReal ucon = ucont[k][j][i].z;

                PetscReal up = 0.5 * ( ucon + fabs(ucon) );
                PetscReal um = 0.5 * ( ucon - fabs(ucon) );

                //IBM active and on the right face of the IBM fluid cell
                // Use second order if positive flux and quick if negative flux
                if ( teqn->access->flags->isIBMActive && k!=mz-2 && isIBMCell(k, j, i, nvert) )
                {
                    div[k][j][i].z =
                            -um * (0.125 * (-    tmprt[k+2][j][i] -  2. * tmprt[k+1][j][i] +  3. * tmprt[k  ][j][i]) + tmprt[k+1][j][i])   +
                            -up * (0.125 * (-    tmprt[k  ][j][i] -  2. * tmprt[k  ][j][i] +  3. * tmprt[k+1][j][i]) +  tmprt[k][j][i  ]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if ( teqn->access->flags->isIBMActive &&  k!=0 && isIBMCell(k+1, j, i, nvert) )
                {
                    div[k][j][i].z =
                            -um * (0.125 * (-    tmprt[k+1][j][i] -  2. * tmprt[k+1][j][i] +  3. * tmprt[k  ][j][i]) + tmprt[k+1][j][i])   +
                            -up * (0.125 * (-    tmprt[k-1][j][i] -  2. * tmprt[k  ][j][i] +  3. * tmprt[k+1][j][i]) +  tmprt[k][j][i  ]);

                }

                else
                {
                    div[k][j][i].z =
                    - ucon
                    * centralUpwind
                    (
                        tmprt[kL][j][i],
                        tmprt[k][j][i],
                        tmprt[k+1][j][i],
                        tmprt[kR][j][i],
                        ucont[k][j][i].z
                        ,limiter[k][j][i].z
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal kappaEff;

                if
                (
                    teqn->access->flags->isLesActive
                )
                {
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

                    // compute stability depentend turbulent Prandtl number
                    PetscReal gradTdotG = dtde*(-9.81);
                    PetscReal l, delta = pow( 1./ajc, 1./3. );
                    if(gradTdotG < 0.)
                    {
                        l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                    }
                    else
                    {
                        l = delta;
                    }

                    PetscReal Prt = 1.0 / (1.0 + (2.0 * l / delta));

                    kappaEff = (nu / cst->Pr) + (nut / Prt);
                }
                else
                {
                    kappaEff = nu / cst->Pr;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].z += (g11 * dtdc + g21 * dtde + g31 * dtdz) * ajc * (kappaEff);
            }
        }
    }

    DMDAVecRestoreArray(fda, teqn->lDivT, &div);
    DMDAVecRestoreArray(fda, teqn->lViscT, &visc);

    DMLocalToLocalBegin(fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalEnd  (fda, teqn->lDivT,  INSERT_VALUES, teqn->lDivT);
    DMLocalToLocalBegin(fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);
    DMLocalToLocalEnd  (fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);

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
                if (isIBMCell(k, j, i, nvert))
                {
                    rhs[k][j][i] = 0;
                }
                else
                {
                    rhs[k][j][i]
                    =
                    scale *
                    (
                         div[k][j][i].x  - div[k  ][j  ][i-1].x +
                         div[k][j][i].y  - div[k  ][j-1][i  ].y +
                         div[k][j][i].z  - div[k-1][j  ][i  ].z +
                        visc[k][j][i].x - visc[k  ][j  ][i-1].x +
                        visc[k][j][i].y - visc[k  ][j-1][i  ].y +
                        visc[k][j][i].z - visc[k-1][j  ][i  ].z
                    ) * aj[k][j][i];
                }

            }
        }
    }

    if(teqn->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
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
    DMDAVecRestoreArray(da,  mesh->lNvert,      &nvert);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

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

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    const PetscReal dt = clock->dt;

    FormT(teqn, Rhs, 1.0);

    // add viscous and transport terms
    if(clock->it > clock->itStart)
    {
        VecAXPY(Rhs, 0.5, teqn->Rhs_o);
        FormT(teqn, Rhs, 0.5);
    }
    else
    {
        FormT(teqn, Rhs, 1.0);
    }

    if(teqn->access->flags->isXDampingActive)
    {
        dampingSourceT(teqn, Rhs, 1.0);
    }

    // multiply for dt
    VecScale(Rhs, dt);

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, Rhs);

    VecAXPY(Rhs, -1.0, T);
    VecAXPY(Rhs,  1.0, teqn->Tmprt_o);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveTEqn(teqn_ *teqn)
{
    PetscReal     norm;
    PetscInt      iter;
    PetscReal     ts,te;

    mesh_          *mesh = teqn->access->mesh;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "TRSNES: Solving for T, Initial residual = ");

    // set the right hand side to zero and copy the solution in the tmp variable
    VecSet (teqn->Rhs, 0.0);
    VecCopy(teqn->Tmprt, teqn->TmprtTmp);

    // solve temperature equation and compute solution norm and iteration number
    SNESSolve(teqn->snesT, PETSC_NULL, teqn->TmprtTmp);
    SNESGetFunctionNorm(teqn->snesT, &norm);
    SNESGetIterationNumber(teqn->snesT, &iter);

    // store the solution and global to local scatter
    VecCopy(teqn->TmprtTmp, teqn->Tmprt);

    // scatter
    DMGlobalToLocalBegin(mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (mesh->da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    // compute elapsed time
    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM,"Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n", norm, iter, te-ts);

    return(0);
}
