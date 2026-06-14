//! \file  aeqn_solvers.c
//! \brief Alpha-Water-equation solver function implementations.

#include "aeqn_solvers.h"

PetscErrorCode FormALowOrder(aeqn_ *aeqn)
{
    // compute first-order upwind face fluxes and store in lDivA
    // divergence accumulation is done separately in ComputeLowOrderUpdate / ApplyMULESLimiter

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

    VecSet(aeqn->lDivA, 0.0);
    DMDAVecGetArray(fda, aeqn->lDivA, &div);

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

    DMDAVecRestoreArray(fda, aeqn->lDivA, &div);
    DMDAVecRestoreArray(da,  aeqn->lAlpha, &alpha);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);

    DMLocalToLocalBegin(fda, aeqn->lDivA, INSERT_VALUES, aeqn->lDivA);
    DMLocalToLocalEnd  (fda, aeqn->lDivA, INSERT_VALUES, aeqn->lDivA);

    resetFacePeriodicFluxesVector(mesh, aeqn->lDivA, aeqn->lDivA, "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormAHighOrder(aeqn_ *aeqn)
{
    // compute fourth-order central face fluxes and store in lDivAHO
    // falls back to second-order central at overset interpolated faces

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

    // ---------------------------------------------------------------------- //
    //               central4 face flux assembly (i, j, k)                   //
    // ---------------------------------------------------------------------- //

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
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, aeqn->lDivAHO,  &div);
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
    // provisional low-order cell update: alpha_LO = alpha - dt * aj * divBudget_LO
    // clips to [0,1] and syncs ghosts for use in the mules limiter

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
    DMDAVecGetArray(fda, aeqn->lDivA,    &div);

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

                // clip to [0,1] for the provisional estimate used by the limiter budgets
                alphaLO[k][j][i] = PetscMax(0.0, PetscMin(1.0, alpha[k][j][i] - dt * aj[k][j][i] * divBudget));
            }
        }
    }

    DMDAVecRestoreArray(da,  aeqn->lAlpha,   &alpha);
    DMDAVecRestoreArray(da,  aeqn->lAlphaLO, &alphaLO);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivA,    &div);

    // sync ghost cells so neighboring processors can read alpha_LO
    DMLocalToLocalBegin(da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);
    DMLocalToLocalEnd  (da, aeqn->lAlphaLO, INSERT_VALUES, aeqn->lAlphaLO);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AddCompressionFlux(aeqn_ *aeqn)
{
    // compression flux is 
    // div ( - uc_face dot face_normal) times alpha_face * (1 - alpha_face)

    // where uc_face = c * mag_u_face * grad_alpha_face / mag_grad_alpha_face

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

    PetscReal     ***alphaLO, ***aj, ***rhs;
    PetscReal     ***Rplus, ***Rminus;
    Cmpnts        ***divLO, ***divCor, ***lambda;

    PetscReal     dt  = clock->dt;
    PetscReal     eps = 1.0e-12;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // ---------------------------------------------------------------------- //
    //    step 1: correction flux = high-order minus low-order               //
    // ---------------------------------------------------------------------- //

    VecCopy (aeqn->lDivAHO,  aeqn->lDivACor);
    VecAXPY (aeqn->lDivACor, -1.0, aeqn->lDivA);

    // ---------------------------------------------------------------------- //
    //    step 2: provisional low-order cell update                           //
    // ---------------------------------------------------------------------- //

    ComputeLowOrderUpdate(aeqn);

    // ---------------------------------------------------------------------- //
    //    step 3: per-cell R_plus and R_minus                                 //
    // ---------------------------------------------------------------------- //

    // initialize to 1.0: physical boundary ghost cells allow full correction
    VecSet(aeqn->lRplusA,  1.0);
    VecSet(aeqn->lRminusA, 1.0);

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

                c =  dt * aj[k][j][i] * divCor[k][j][i].x;    // right i-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -dt * aj[k][j][i] * divCor[k][j][i-1].x;  // left i-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c =  dt * aj[k][j][i] * divCor[k][j][i].y;    // right j-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -dt * aj[k][j][i] * divCor[k][j-1][i].y;  // left j-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c =  dt * aj[k][j][i] * divCor[k][j][i].z;    // right k-face
                if(c > 0.0) C_plus += c; else C_minus -= c;

                c = -dt * aj[k][j][i] * divCor[k-1][j][i].z;  // left k-face
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

    // ---------------------------------------------------------------------- //
    //    step 4: per-face limiter coefficient lambda_f in [0,1]              //
    // ---------------------------------------------------------------------- //

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

                // i-face: owner=(k,j,i), neighbor=(k,j,i+1)
                // negative dF -> correction decreases owner alpha, increases neighbor
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

    // ---------------------------------------------------------------------- //
    //    step 5: final limited Rhs = low-order + lambda-weighted correction  //
    // ---------------------------------------------------------------------- //

    DMDAVecGetArray(da,  aeqn->Rhs,      &rhs);
    DMDAVecGetArray(da,  mesh->lAj,      &aj);
    DMDAVecGetArray(fda, aeqn->lDivA,    &divLO);
    DMDAVecGetArray(fda, aeqn->lDivACor, &divCor);
    DMDAVecGetArray(fda, aeqn->lLambdaA, &lambda);

    VecSet(aeqn->Rhs, 0.0);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal loBudget =
                    divLO [k][j][i].x - divLO [k][j][i-1].x +
                    divLO [k][j][i].y - divLO [k][j-1][i].y +
                    divLO [k][j][i].z - divLO [k-1][j][i].z;

                PetscReal corBudget =
                    lambda[k][j][i].x   * divCor[k][j][i].x   - lambda[k][j][i-1].x * divCor[k][j][i-1].x +
                    lambda[k][j][i].y   * divCor[k][j][i].y   - lambda[k][j-1][i].y * divCor[k][j-1][i].y +
                    lambda[k][j][i].z   * divCor[k][j][i].z   - lambda[k-1][j][i].z * divCor[k-1][j][i].z;

                rhs[k][j][i] = (loBudget + corBudget) * aj[k][j][i];
            }
        }
    }

    DMDAVecRestoreArray(da,  aeqn->Rhs,      &rhs);
    DMDAVecRestoreArray(da,  mesh->lAj,      &aj);
    DMDAVecRestoreArray(fda, aeqn->lDivA,    &divLO);
    DMDAVecRestoreArray(fda, aeqn->lDivACor, &divCor);
    DMDAVecRestoreArray(fda, aeqn->lLambdaA, &lambda);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AeqnSNES(aeqn_ *aeqn)
{
    mesh_  *mesh  = aeqn->access->mesh;
    clock_ *clock = aeqn->access->clock;

    PetscReal     norm;
    PetscInt      iter;
    PetscReal     ts, te;

    PetscTime(&ts);
    PetscPrintf(mesh->MESH_COMM, "%s%s: Solving for Alpha-Water, Initial residual = ", aeqn->ddtScheme.c_str() ,aeqn->solverType.c_str());
    
    // compute initial guess
    VecCopy(aeqn->Alpha_o, aeqn->AlphaTmp);    

    SNESSolve(aeqn->snesA, PETSC_NULL, aeqn->AlphaTmp);
    SNESGetFunctionNorm(aeqn->snesA, &norm);
    SNESGetIterationNumber(aeqn->snesA, &iter);

    // report total inner linear iterations (= total MatMFFD FormA calls for J*v)
    PetscInt linIter = 0;
    SNESGetLinearSolveIterations(aeqn->snesA, &linIter);

    VecCopy(aeqn->AlphaTmp, aeqn->Alpha);

    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);

    PetscTime(&te);
    PetscPrintf(mesh->MESH_COMM, "Final residual = %e, Iterations = %ld (linear = %ld), Elapsed Time = %lf\n", norm, iter, linIter, te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SNESFuncEvalA(SNES snes, Vec Alpha, Vec Rhs, void *ptr)
{
    aeqn_ *aeqn   = (aeqn_*)ptr;
    mesh_ *mesh   = aeqn->access->mesh;
    clock_ *clock = aeqn->access->clock;

    // scatter trial solution from snes into aeqn->Alpha and lAlpha
    VecCopy(Alpha, aeqn->Alpha);
    DMGlobalToLocalBegin(mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    DMGlobalToLocalEnd  (mesh->da, aeqn->Alpha, INSERT_VALUES, aeqn->lAlpha);
    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // build low-order and high-order face fluxes, apply mules limiter
    FormALowOrder (aeqn);
    FormAHighOrder(aeqn);
    AddCompressionFlux(aeqn);
    ApplyMULESLimiter(aeqn);
    resetNonResolvedCellCentersScalar(mesh, aeqn->Rhs);

    // scale by dt and add time-derivative residual terms
    PetscReal dt = clock->dt;
    VecScale(aeqn->Rhs, dt);

    if(aeqn->ddtScheme == "BDF2" && clock->it > clock->itStart)
    {
        VecAXPY(aeqn->Rhs, -3.0/2.0, Alpha);
        VecAXPY(aeqn->Rhs,  2.0,     aeqn->Alpha_o);
        VecAXPY(aeqn->Rhs, -1.0/2.0, aeqn->Alpha_oo);
    }
    else
    {
        VecAXPY(aeqn->Rhs, -1.0, Alpha);
        VecAXPY(aeqn->Rhs,  1.0, aeqn->Alpha_o);
    }

    // copy residual to output
    VecCopy(aeqn->Rhs, Rhs);

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
    ApplyMULESLimiter(aeqn);

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