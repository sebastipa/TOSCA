//! \file  ceqn.c
//! \brief Contains C equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorC(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeCEqn(ceqn_ *ceqn)
{
    if(ceqn != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = ceqn->access->mesh;

        VecDuplicate(mesh->Nvert, &(ceqn->ConcTmp));    VecSet(ceqn->ConcTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(ceqn->Conc));       VecSet(ceqn->Conc,    0.0);
        VecDuplicate(mesh->Nvert, &(ceqn->Conc_o));     VecSet(ceqn->Conc_o,  0.0);
        VecDuplicate(mesh->Nvert, &(ceqn->Rhs));         VecSet(ceqn->Rhs,      0.0);
        VecDuplicate(mesh->Nvert, &(ceqn->Rhs_o));       VecSet(ceqn->Rhs_o,    0.0);
        VecDuplicate(mesh->lAj,   &(ceqn->lConc));      VecSet(ceqn->lConc,   0.0);
        VecDuplicate(mesh->lAj,   &(ceqn->lConc_o));    VecSet(ceqn->lConc_o, 0.0);

        VecDuplicate(mesh->lCent, &(ceqn->lDivC));       VecSet(ceqn->lDivC,    0.0);
        VecDuplicate(mesh->lCent, &(ceqn->lViscC));      VecSet(ceqn->lViscC,   0.0);

        ceqn->cOutTotal = 0;
        ceqn->cInTotal = 0;
        ceqn->cOutProc = 0;
        ceqn->cInProc = 0;

        // read time discretization scheme
        readDictWord("control.dat", "-dCdtScheme", &(ceqn->ddtScheme));

        //need to set cRef wen abl is not active still
        if (ceqn->access->abl == NULL)
        {
            PetscPrintf(mesh->MESH_COMM, "Reading NON ABL CREF ...\n");
            readDictDouble("./boundary/C", "cRefNoAbl", &(ceqn->cRefNoAbl));
            //printf("tref = %lf", teqn->tRefNoAbl);
        }

        // create the SNES solver
        if(ceqn->ddtScheme=="backwardEuler")
        {
            // readDictDouble("control.dat", "-absTolT", &(teqn->absExitTol));
            // readDictDouble("control.dat", "-relTolT", &(teqn->relExitTol));

            // default parameters
            ceqn->absExitTol        = 1e-5;
            ceqn->relExitTol        = 1e-30;

            // input file
            PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolC",  &(ceqn->absExitTol), PETSC_NULL);
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolC",  &(ceqn->relExitTol), PETSC_NULL);

            SNESCreate(mesh->MESH_COMM,&(ceqn->snesC));
            SNESMonitorSet(ceqn->snesC, SNESMonitorC, (void*)&(mesh->MESH_COMM), PETSC_NULL);

            // set the SNES evaluating function
            SNESSetFunction(ceqn->snesC, ceqn->Rhs, CeqnSNES, (void *)ceqn);

            // create jacobian matrix
            MatCreateSNESMF(ceqn->snesC, &(ceqn->JC));
            SNESSetJacobian(ceqn->snesC, ceqn->JC, ceqn->JC, MatMFFDComputeJacobian, (void *)ceqn);

            // set SNES solver type
            //SNESSetType(teqn->snesT, SNESNEWTONTR);           //SNESTR
            SNESSetType(ceqn->snesC, SNESNEWTONLS);        //SNESLS is better for stiff PDEs such as the one including IB but slower

            // set SNES solve and step failures
            SNESSetMaxLinearSolveFailures(ceqn->snesC,10000);
            SNESSetMaxNonlinearStepFailures(ceqn->snesC,10000);
            SNESKSPSetUseEW(ceqn->snesC, PETSC_TRUE);

            // set SNES Krylov Sub-Space parameters
            SNESKSPSetParametersEW(ceqn->snesC,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

            // SNES tolerances
            // 2nd arg: absolute tolerance
            // 3rd arg: relative tolerance
            // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
            // 5th arg: maximum number of iterations
            // 6th arg: maximum function evaluations
            SNESSetTolerances(ceqn->snesC, ceqn->absExitTol, 1e-30, 1e-30, 20, 1000);

            SNESGetKSP(ceqn->snesC, &(ceqn->ksp));
            KSPGetPC(ceqn->ksp,&(ceqn->pc));

            // set KSP solver type
            KSPSetType(ceqn->ksp, KSPGMRES);

            //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);    //2009.09.22 poor performance
            //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);      //2009.09.22

            //KSPFischerGuess itg;
            //KSPFischerGuessCreate(ksp,1,100,&itg);
            //KSPSetFischerGuess(ksp, itg);                  //2009.09.22

            //KSPGMRESSetPreAllocateVectors(ksp);            --> crazy thing consumes memory

            PCSetType(ceqn->pc, PCNONE);
            PetscReal rtol=ceqn->relExitTol, atol=ceqn->absExitTol, dtol=PETSC_DEFAULT;
            KSPSetTolerances(ceqn->ksp, rtol, atol, dtol, 1000);
        }
        else if (ceqn->ddtScheme=="rungeKutta4")
        {

        }
        else
        {
            char error[512];
            sprintf(error, "unknown ddtScheme %s for C equation, available schemes are\n    1. backwardEuler\n", ceqn->ddtScheme.c_str());
            fatalErrorInFunction("InitializeTEqn", error);
        }
    }

    return(0);
}

//***************************************************************************************************************//
//not used in ceqn yet.
/*PetscErrorCode dampingSourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale)
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
}*/

//***************************************************************************************************************//

PetscErrorCode FormC(ceqn_ *ceqn, Vec &Rhs, PetscReal scale)
{
    // In this function the viscous + divergence term of the concentration equation are
    // discretized at cell centers.
    // First the divergence and viscous fluxes are evaluated at cell faces, then
    // their budget is evaluated at the internal cells, forming the Rhs.

    mesh_         *mesh  = ceqn->access->mesh;
    ueqn_         *ueqn  = ceqn->access->ueqn;
    les_          *les   = ceqn->access->les;
    constants_    *cst   = ceqn->access->constants;
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

    PetscReal     ***conc, ***conc_o, ***rhs, ***nvert;

    PetscReal     ***markPore;

    Cmpnts        ***div, ***visc;                                                // divergence and viscous terms
    Cmpnts        ***limiter;                                                     // flux limiter
    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;                                  // cell and face jacobians
    PetscReal     ***lnu_t;

    PetscReal     dcdc, dcde, dcdz;                                              // concentr der. w.r.t. curvil. coords
    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;          // surface area vectors components
    PetscReal     g11, g21, g31;                                                 // metric tensor components

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  Rhs,              &rhs);
    DMDAVecGetArray(da,  ceqn->lConc,      &conc);
    //DMDAVecGetArray(da,  ceqn->lConc_o,    &conc_o);
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

    VecSet(ceqn->lDivC,  0.0);
    VecSet(ceqn->lViscC, 0.0);;

    DMDAVecGetArray(fda, ceqn->lDivC, &div);
    DMDAVecGetArray(fda, ceqn->lViscC, &visc);
    DMDAVecGetArray(da, mesh->poreMarkers,  &markPore);

    if(ceqn->access->flags->isLesActive)
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

                Compute_dscalar_i(mesh, i, j, k, mx, my, mz, conc, nvert, &dcdc, &dcde, &dcdz);

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
                if(ceqn->access->flags->isIBMActive && i!=mx-2 && isIBMCell(k, j, i, nvert))

                {
                    div[k][j][i].x =
                            -um * (0.125 * (-    conc[k][j][i+2] - 2. * conc[k][j][i+1] + 3. * conc[k][j][i]) + conc[k][j][i+1]) +
                            -up * (0.125 * (-    conc[k][j][i] -  2. * conc[k][j][i] +  3. * conc[k][j][i+1]) +  conc[k][j][i]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if (ceqn->access->flags->isIBMActive && i!=0 && isIBMCell(k, j, i+1, nvert))
                {
                    div[k][j][i].x =
                            -um * (0.125 * (-    conc[k][j][i+1] -  2. * conc[k][j][i+1] +  3. * conc[k][j][i]) + conc[k][j][i+1]) +
                            -up * (0.125 * (-    conc[k][j][i-1] -  2. * conc[k][j][i] +  3. * conc[k][j][i+1]) + conc[k][j][i]);

                }
                else
                {
                    div[k][j][i].x =
                    - ucon
                    * centralUpwind
                    (
                        conc[k][j][iL],
                        conc[k][j][i],
                        conc[k][j][i+1],
                        conc[k][j][iR],
                        ucont[k][j][i].x
                        ,limiter[k][j][i].x
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;


                    // viscous terms
                    if
                    (
                        ceqn->access->flags->isLesActive
                    )
                    {
                        nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j][i+1]);

                        /*// compute stability dependent turbulent Prandtl number
                        PetscReal gradTdotG = dcde*(-9.81);
                        PetscReal l, delta = pow( 1./ajc, 1./3. );
                        if(gradTdotG < 0.)
                        {
                            l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                        }
                        else
                        {
                            l = delta;
                        }*/

                        PetscReal Scht = 0.7;

                        gammaEff = (nu / cst->Sch) + (nut / Scht);

                        // wall model i-left/right patch. Will always be skipped for concentration
                        /*if
                        (
                            (mesh->boundaryC.iLeft=="thetaWallFunction" && i==0) ||
                            (mesh->boundaryC.iRight=="thetaWallFunction" && i==mx-2)
                        )
                        {
                            PetscReal signQ =  1.0;
                            if(i==0)  signQ = -1.0;

                            visc[k][j][i].y
                            =
                            signQ *
                            (
                                teqn->iRWM->qWall.x[k][j] * icsi[k][j][i].x +
                                teqn->iRWM->qWall.y[k][j] * icsi[k][j][i].y +
                                teqn->iRWM->qWall.z[k][j] * icsi[k][j][i].z
                            );

                            kappaEff = 0.0;
                        }*/
                    }
                    else
                    {
                        gammaEff = nu / cst->Sch;
                    }

                    // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                    visc[k][j][i].x += (g11 * dcdc + g21 * dcde + g31 * dcdz) * ajc * (gammaEff);



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

                Compute_dscalar_j(mesh, i, j, k, mx, my, mz, conc, nvert, &dcdc, &dcde, &dcdz);

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
                if ( ceqn->access->flags->isIBMActive &&  j!=my-2 && isIBMCell(k, j, i, nvert) )
                {
                    div[k][j][i].y =
                            -um * (0.125 * (-    conc[k][j+2][i] -  2. * conc[k][j+1][i] +  3. * conc[k][j  ][i]) +  conc[k][j+1][i]) +
                            -up * (0.125 * (-    conc[k][j  ][i] -  2. * conc[k][j  ][i] +  3. * conc[k][j+1][i]) +  conc[k][j][i  ]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if ( ceqn->access->flags->isIBMActive &&  j!=0 && isIBMCell(k, j+1, i, nvert) )
                {
                    div[k][j][i].y =
                            -um * (0.125 * (-    conc[k][j+1][i] -  2. * conc[k][j+1][i] +  3. * conc[k][j  ][i]) + conc[k][j+1][i]) +
                            -up * (0.125 * (-    conc[k][j-1][i] -  2. * conc[k][j  ][i] +  3. * conc[k][j+1][i]) + conc[k][j][i  ]);

                }

                else
                {
                    div[k][j][i].y =
                    - ucon
                    * centralUpwind
                    (
                        conc[k][jL][i],
                        conc[k][j][i],
                        conc[k][j+1][i],
                        conc[k][jR][i],
                        ucont[k][j][i].y
                        ,limiter[k][j][i].y
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;


                if
                (
                    ceqn->access->flags->isLesActive
                )
                {
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k][j+1][i]);

                    /*// compute stability depentend turbulent Prandtl number
                    PetscReal gradTdotG = dcde*(-9.81);
                    PetscReal l, delta = pow( 1./ajc, 1./3. );
                    if(gradTdotG < 0.)
                    {
                        l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                    }
                    else
                    {
                        l = delta;
                    }*/

                    PetscReal Scht = 0.7;

                    gammaEff = (nu / cst->Sch) + (nut / Scht);

                    // wall model j-left patch
                    /*if
                    (
                        (mesh->boundaryT.jLeft=="thetaWallFunction"  && j==0) ||
                        (mesh->boundaryT.jRight=="thetaWallFunction" && j==my-2)
                    )
                    {
                        PetscReal signQ =  1.0;
                        if(j==0)  signQ = -1.0;

                        visc[k][j][i].y
                        =
                        signQ *
                        (
                            teqn->jLWM->qWall.x[k][i] * jeta[k][j][i].x +
                            teqn->jLWM->qWall.y[k][i] * jeta[k][j][i].y +
                            teqn->jLWM->qWall.z[k][i] * jeta[k][j][i].z
                        );

                        kappaEff = 0.0;
                    }/*/
                }
                else
                {
                    gammaEff = nu / cst->Sch;
                }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].y += (g11 * dcdc + g21 * dcde + g31 * dcdz) * ajc * (gammaEff);


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

                Compute_dscalar_k(mesh, i, j, k, mx, my, mz, conc, nvert, &dcdc, &dcde, &dcdz);

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
                if ( ceqn->access->flags->isIBMActive && k!=mz-2 && isIBMCell(k, j, i, nvert) )
                {
                    div[k][j][i].z =
                            -um * (0.125 * (-    conc[k+2][j][i] -  2. * conc[k+1][j][i] +  3. * conc[k  ][j][i]) + conc[k+1][j][i])   +
                            -up * (0.125 * (-    conc[k  ][j][i] -  2. * conc[k  ][j][i] +  3. * conc[k+1][j][i]) +  conc[k][j][i  ]);

                }

                //IBM active and on the left face of the IBM fluid cell
                // Use second order if negative flux and quickDiv if positive flux

                else if ( ceqn->access->flags->isIBMActive &&  k!=0 && isIBMCell(k+1, j, i, nvert) )
                {
                    div[k][j][i].z =
                            -um * (0.125 * (-    conc[k+1][j][i] -  2. * conc[k+1][j][i] +  3. * conc[k  ][j][i]) + conc[k+1][j][i])   +
                            -up * (0.125 * (-    conc[k-1][j][i] -  2. * conc[k  ][j][i] +  3. * conc[k+1][j][i]) +  conc[k][j][i  ]);

                }

                else
                {
                    div[k][j][i].z =
                    - ucon
                    * centralUpwind
                    (
                        conc[kL][j][i],
                        conc[k][j][i],
                        conc[k+1][j][i],
                        conc[kR][j][i],
                        ucont[k][j][i].z
                        ,limiter[k][j][i].z
                    );
                }

                PetscReal nu = cst->nu, nut;
                PetscReal gammaEff;


                 if
                 (
                     ceqn->access->flags->isLesActive
                 )
                 {
                    nut = 0.5 * (lnu_t[k][j][i] + lnu_t[k+1][j][i]);

                    /*// compute stability depentend turbulent Prandtl number
                    PetscReal gradTdotG = dcde*(-9.81);
                    PetscReal l, delta = pow( 1./ajc, 1./3. );
                    if(gradTdotG < 0.)
                    {
                        l = PetscMin(delta, 0.76*std::sqrt(300 / std::fabs(gradTdotG)));
                    }
                    else
                    {
                        l = delta;
                    }*/

                    PetscReal Scht = 0.7;

                    gammaEff = (nu / cst->Sch) + (nut / Scht);
                 }
                 else
                 {
                    gammaEff = nu / cst->Sch;
                 }

                // note: 1/J is the original term, here terms arrive already with a factor of 1/J^2 so actually we multiply for J (ajc)
                visc[k][j][i].z += (g11 * dcdc + g21 * dcde + g31 * dcdz) * ajc * (gammaEff);


            }
        }
    }

    DMDAVecRestoreArray(fda, ceqn->lDivC, &div);
    DMDAVecRestoreArray(fda, ceqn->lViscC, &visc);

    DMLocalToLocalBegin(fda, ceqn->lDivC,  INSERT_VALUES, ceqn->lDivC);
    DMLocalToLocalEnd  (fda, ceqn->lDivC,  INSERT_VALUES, ceqn->lDivC);
    DMLocalToLocalBegin(fda, ceqn->lViscC, INSERT_VALUES, ceqn->lViscC);
    DMLocalToLocalEnd  (fda, ceqn->lViscC, INSERT_VALUES, ceqn->lViscC);

    resetFacePeriodicFluxesVector(mesh, ceqn->lDivC,   ceqn->lDivC, "localToLocal");
    resetFacePeriodicFluxesVector(mesh, ceqn->lViscC, ceqn->lViscC, "localToLocal");

    DMDAVecGetArray(fda, ceqn->lDivC, &div);
    DMDAVecGetArray(fda, ceqn->lViscC, &visc);
    DMDAVecGetArray(da,  ceqn->Conc_o,      &conc_o);


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
                /*else if (markPore[k][j][i] < 1) //if in porous zone apply filtrationEfficiency
                {

                    rhs[k][j][i]
                    =
                    scale * (1 - ceqn->access->pores->filtrationEfficiency) *
                    (
                         div[k][j][i].x  - div[k  ][j  ][i-1].x +
                         div[k][j][i].y  - div[k  ][j-1][i  ].y +
                         div[k][j][i].z  - div[k-1][j  ][i  ].z +
                        visc[k][j][i].x - visc[k  ][j  ][i-1].x +
                        visc[k][j][i].y - visc[k  ][j-1][i  ].y +
                        visc[k][j][i].z - visc[k-1][j  ][i  ].z
                    ) * aj[k][j][i];

                    //conc_o[k][j][i] *= (1 - ceqn->access->pores->filtrationEfficiency);

                }*/
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

    if(ceqn->access->flags->isLesActive)
    {
        DMDAVecRestoreArray(da, les->lNu_t, &lnu_t);
    }

    DMDAVecRestoreArray(fda, ceqn->lDivC,       &div);
    DMDAVecRestoreArray(fda, ceqn->lViscC,      &visc);

    DMDAVecRestoreArray(da,  Rhs,              &rhs);
    DMDAVecRestoreArray(da,  ceqn->lConc,      &conc);
    DMDAVecRestoreArray(da,  ceqn->Conc_o,    &conc_o);
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
    DMDAVecRestoreArray(da, mesh->poreMarkers,  &markPore);

    return(0);
}


//***************************************************************************************************************//


PetscErrorCode CeqnSNES(SNES snes, Vec C, Vec Rhs, void *ptr)
{
    ceqn_ *ceqn   = (ceqn_*)ptr;
    mesh_ *mesh   = ceqn->access->mesh;
    clock_ *clock = ceqn->access->clock;
    VecCopy(C, ceqn->Conc);

    // scatter concentration from global to local
    DMGlobalToLocalBegin(mesh->da, ceqn->Conc, INSERT_VALUES, ceqn->lConc);
    DMGlobalToLocalEnd  (mesh->da, ceqn->Conc, INSERT_VALUES, ceqn->lConc);

    // reset concentration periodic fluxes to be consistent if the flow is periodic
    resetCellPeriodicFluxes(mesh, ceqn->Conc, ceqn->lConc, "scalar", "globalToLocal");

    // update wall model (optional)
    //UpdateWallModelsT(ceqn); no wall models for C yet.

    // initialize the rhs vector
    VecSet(Rhs, 0.0);

    // get time step
    const PetscReal dt = clock->dt;

    // add viscous and transport terms
    if(clock->it > clock->itStart)
    {
        VecAXPY(Rhs, 0.5, ceqn->Rhs_o);
        FormC(ceqn, Rhs, 0.5);
    }
    else
    {
        FormC(ceqn, Rhs, 1.0);
    }

    if(ceqn->access->flags->isXDampingActive)
    {
        // this is causing spurious oscillation: deactivate
        // dampingSourceT(ceqn, Rhs, 1.0);
    }

    // multiply for dt
    VecScale(Rhs, dt);

    // set to zero at non-resolved cell faces
    resetNonResolvedCellCentersScalar(mesh, Rhs);

    VecAXPY(Rhs, -1.0, C);

    //filterConcOld(ceqn); //apply filtration to Conc_o term in porous zones.
    VecAXPY(Rhs,  1.0, ceqn->Conc_o);

    return(0);
}


//***************************************************************************************************************//

PetscErrorCode SolveCEqn(ceqn_ *ceqn)
{
    mesh_          *mesh = ceqn->access->mesh;

    // set the right hand side to zero
    VecSet (ceqn->Rhs, 0.0);

    if(ceqn->ddtScheme=="backwardEuler")
    {
        PetscReal     norm;
        PetscInt      iter;
        PetscReal     ts,te;

        PetscTime(&ts);
        PetscPrintf(mesh->MESH_COMM, "TRSNES: Solving for C, Initial residual = ");

        // copy the solution in the tmp variable
        VecCopy(ceqn->Conc, ceqn->ConcTmp);

        // solve temperature equation and compute solution norm and iteration number
        SNESSolve(ceqn->snesC, PETSC_NULL, ceqn->ConcTmp);
        SNESGetFunctionNorm(ceqn->snesC, &norm);
        SNESGetIterationNumber(ceqn->snesC, &iter);

        // store the solution and global to local scatter
        VecCopy(ceqn->ConcTmp, ceqn->Conc);

        resetNoNegativeConc(ceqn);
        filterConc(ceqn);

        // scatter
        DMGlobalToLocalBegin(mesh->da, ceqn->Conc, INSERT_VALUES, ceqn->lConc);
        DMGlobalToLocalEnd  (mesh->da, ceqn->Conc, INSERT_VALUES, ceqn->lConc);

        // compute elapsed time
        PetscTime(&te);
        PetscPrintf(mesh->MESH_COMM,"Final residual = %e, Iterations = %ld, Elapsed Time = %lf\n", norm, iter, te-ts);
    }
    /*else if (ceqn->ddtScheme=="rungeKutta4")
    {
        TeqnRK4(teqn);
    }*/ //only backwardeuler available for C.


    concSummary(ceqn);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode concSummary(ceqn_ *ceqn)
{
    mesh_           *mesh = ceqn->access->mesh;
    ueqn_           *ueqn  = ceqn->access->ueqn;
    vents_          *vents = ceqn->access->vents;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscMPIInt   rank;

    Cmpnts        ***ucont;

    PetscReal     ***conc;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscInt      q;

    PetscReal     ***markPore;
    PetscInt      ***markVent;

    PetscScalar      lConcIn = 0.0, lConcOut = 0.0, lConcPore = 0.0, lConcAir = 0.0,
                     concIn = 0.0, concOut = 0.0, concPore = 0.0, concAir = 0.0;


    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, mesh->poreMarkers, &markPore);
    DMDAVecGetArray(da, mesh->ventMarkers, &markVent);
    DMDAVecGetArray(da,  ceqn->lConc,      &conc);
    DMDAVecGetArray(fda, ueqn->lUcont,     &ucont);

    //sum Conc in differetn areas of domain.

      for (k=lzs; k<lze; k++)
      {
         for (j=lys; j<lye; j++)
         {
            for (i=lxs; i<lxe; i++)
            {

                if (markPore[k][j][i] < 1) //sum local particle concentration in and out of porous zones.
                {
                    lConcPore += conc[k][j][i];
                }
                else
                {
                    lConcAir += conc[k][j][i];
                }

                if (markVent[k][j][i] > 0) //sum local concentration fluxes in and out of domain. assumes all inlets and outlets are vents. add leaks laters on.
                {
                    q=markVent[k][j][i]-1;

                    if (vents->vent[q]->dir == "outlet")
                    {

                        if (vents->vent[q]->face == "iLeft")
                        {

                            if (ucont[k][j][i-1].x < 0) //only assume Conc is leavin domain if te associated Ucont is out of domain too.
                            {
                                lConcOut   +=  conc[k][j][i];
                            }

                        }

                        if (vents->vent[q]->face == "iRight")
                        {

                            if (ucont[k][j][i].x > 0)
                            {
                                lConcOut   +=  conc[k][j][i+1];
                            }

                        }

                        if (vents->vent[q]->face == "jLeft")
                        {

                            if (ucont[k][j-1][i].y < 0)
                            {
                                lConcOut   +=  conc[k][j-1][i];

                            }

                        }

                        if (vents->vent[q]->face == "jRight")
                        {

                            if (ucont[k][j][i].y > 0)
                            {
                                lConcOut   +=  conc[k][j+1][i];
                            }

                        }

                        if (vents->vent[q]->face == "kLeft")
                        {

                            if (ucont[k-1][j][i].z < 0)
                            {
                                lConcOut   +=  conc[k-1][j][i];
                            }

                        }

                        if (vents->vent[q]->face == "kRight")
                        {

                            if (ucont[k][j][i].z > 0)
                            {
                                lConcOut   +=  conc[k+1][j][i];
                            }

                        }

                    }
                    if (vents->vent[q]->dir == "inlet")
                    {

                        if (vents->vent[q]->face == "iLeft")
                        {

                            if (ucont[k][j][i-1].x > 0)
                            {
                                lConcIn   +=  conc[k][j][i-1];
                            }

                        }

                        if (vents->vent[q]->face == "iRight")
                        {

                            if (ucont[k][j][i].x < 0)
                            {
                                lConcIn   +=  conc[k][j][i+1];
                            }

                        }

                        if (vents->vent[q]->face == "jLeft")
                        {

                            if (ucont[k][j-1][i].y > 0)
                            {
                                lConcIn   +=  conc[k][j-1][i];
                            }

                        }

                        if (vents->vent[q]->face == "jRight")
                        {

                            if (ucont[k][j][i].y < 0)
                            {
                                lConcIn   +=  conc[k][j+1][i];
                            }

                        }

                        if (vents->vent[q]->face == "kLeft")
                        {

                            if (ucont[k-1][j][i].z > 0)
                            {
                                lConcIn   +=  conc[k-1][j][i];
                            }

                        }

                        if (vents->vent[q]->face == "kRight")
                        {

                            if (ucont[k][j][i].z < 0)
                            {
                               lConcIn   +=  conc[k+1][j][i];
                            }

                        }

                    }
                }


            }
          }
       }

      //MPI_Allreduce(&lConcPore, &concPore, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
      //MPI_Allreduce(&lConcAir, &concAir, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
      //MPI_Allreduce(&lConcOut, &concOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
      //MPI_Allreduce(&lConcIn, &concIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

      //PetscMPIInt nProcs;

      //MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);

      //MPI_Barrier(mesh->MESH_COMM);

      ceqn->cOutProc += lConcOut;
      ceqn->cInProc += lConcIn;

      MPI_Barrier(mesh->MESH_COMM);

      MPI_Allreduce(&ceqn->cOutProc, &ceqn->cOutTotal, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
      MPI_Allreduce(&ceqn->cInProc, &ceqn->cInTotal, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);


      PetscPrintf(mesh->MESH_COMM, "\n concOutTotal = %lf, concInTotal= % lf\n", ceqn->cOutTotal, ceqn->cInTotal);
      //printf("\n concOutTotal = %lf, concInTotal= % lf\n", ceqn->cOutTotal, ceqn->cInTotal);

      DMDAVecRestoreArray(da, mesh->poreMarkers, &markPore);
      DMDAVecRestoreArray(da, mesh->ventMarkers, &markVent);
      DMDAVecRestoreArray(da,  ceqn->lConc,      &conc);
      DMDAVecRestoreArray(fda, ueqn->lUcont,     &ucont);

    return(0);
 }

 //***************************************************************************************************************//

 PetscErrorCode filterConc(ceqn_ *ceqn)
 {
     mesh_           *mesh = ceqn->access->mesh;
     DM               da = mesh->da, fda = mesh->fda;
     DMDALocalInfo    info = mesh->info;
     PetscInt         xs = info.xs, xe = info.xs + info.xm;
     PetscInt         ys = info.ys, ye = info.ys + info.ym;
     PetscInt         zs = info.zs, ze = info.zs + info.zm;
     PetscInt         mx = info.mx, my = info.my, mz = info.mz;

     PetscMPIInt   rank;

     PetscReal     ***conc;

     PetscInt      lxs, lxe, lys, lye, lzs, lze;
     PetscInt      i, j, k;

     PetscReal     ***markPore;

     lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
     lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
     lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

     DMDAVecGetArray(da, mesh->poreMarkers, &markPore);
     DMDAVecGetArray(da,  ceqn->Conc,      &conc);


     //change Conc_o to (1-filtrationEfficiency)*Conc_o in porous zone

       for (k=lzs; k<lze; k++)
       {
          for (j=lys; j<lye; j++)
          {
             for (i=lxs; i<lxe; i++)
             {

                 if (markPore[k][j][i] < 1)
                 {
                     conc[k][j][i] *= (1-ceqn->access->pores->filtrationEfficiency);
                 }


             }
           }
        }

       DMDAVecRestoreArray(da, mesh->poreMarkers, &markPore);
       DMDAVecRestoreArray(da,  ceqn->Conc,      &conc);

     return(0);
  }
