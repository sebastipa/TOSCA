//! \file  les.c
//! \brief Contains LES model function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inline.h"

const PetscReal wall_cs   = 0.001;
const PetscReal std_cs    = 0.0289;
const PetscReal amd_cs    = 0.1;
const PetscReal vreman_cs = 0.07225;
const PetscReal sa_cs     = 0.325; 
//***************************************************************************************************************//

PetscErrorCode InitializeLESScalar(les_ *les)
{
    // TOSCA LES models
    // 1. standard Smagorinsky
    // 2. stability dependent
    // 3. dynamic Smagorinsky with box averaging (DSM)
    // 4. dynamic Smagorinsky scale invariant with lagrangian averaging (LASI)
    // 5. dynamic smagorinsky scale dependent with lagrangian averaging (LASD)
    // 6. dynamic smagorinsky scale dependent with plane averaging (PASD)
    // 7. Anisotropic minimum dissipation model (AMD)
    // 8. Vreman model //needs verification
    // 9. Bardina vreman model 
    //10. Bardina AMD model

    if(les != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = les->access->mesh;

        if(les->access->flags->isLesActive)
        {
            VecDuplicate(mesh->lAj, &(les->lk_t)); VecSet(les->lk_t,0.);
            VecDuplicate(mesh->lAj, &(les->lCsk)); VecSet(les->lCsk,0.);

            if( les->model == DSM || 
                les->model == DLASI || 
                les->model == DLASD || 
                les->model == DPASD
            )
            {
                VecDuplicate(mesh->lAj, &(les->lKX));       VecSet(les->lKX, 0.);
                VecDuplicate(mesh->lAj, &(les->lXX));       VecSet(les->lXX, 0.);
                VecDuplicate(mesh->lAj, &(les->ldTheta));   VecSet(les->ldTheta, 0.);

                if (les->model == DLASI || les->model == DLASD)
                {
                    VecDuplicate(mesh->lAj, &(les->lKX_old));   VecSet(les->lKX_old, 0.);
                    VecDuplicate(mesh->lAj, &(les->lXX_old));   VecSet(les->lXX_old, 0.);
                }

                if (les->model == DLASD || les->model == DPASD)
                {
                    VecDuplicate(mesh->lAj, &(les->lPY));  VecSet(les->lPY, 0.);
                    VecDuplicate(mesh->lAj, &(les->lYY));  VecSet(les->lYY, 0.);

                    if(les->model == DLASD)
                    {
                        VecDuplicate(mesh->lAj, &(les->lPY_old));   VecSet(les->lPY_old, 0.);
                        VecDuplicate(mesh->lAj, &(les->lYY_old));   VecSet(les->lYY_old, 0.);
                    }
                }
            }

            if (les->model == BAMD || les->model == BV)
            {
                VecDuplicate(mesh->lCent, &(les->lQ)); VecSet(les->lQ, 0.);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateCsk (les_ *les)
{
    mesh_          *mesh  = les->access->mesh;
    ueqn_          *ueqn  = les->access->ueqn;
    clock_         *clock = les->access->clock;
    teqn_          *teqn  = les->access->teqn;
    DM             da     = mesh->da, fda = mesh->fda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    Cmpnts        ***ucat, ***ucont, ***dtheta;
    Cmpnts        ***csi, ***eta, ***zet,
                  ***icsi, ***ieta, ***izet,
                  ***jcsi, ***jeta, ***jzet,
                  ***kcsi, ***keta, ***kzet;

    PetscReal     ***iaj, ***jaj, ***kaj, ***aj, ***lt;

    PetscReal     ***nvert, ***csk, ***meshTag;

    Cmpnts        ***Ax, ***Ay, ***Az, ***cent;
    PetscReal     ***KX, ***XX, ***PY, ***YY;
    PetscReal     ***Sabs;

    PetscReal     ajc;

    PetscReal     dudc, dude, dudz,
                  dvdc, dvde, dvdz,
                  dwdc, dwde, dwdz;

    const char*   lagrangianAverageCs = "trilinear";    // nearestPoint
                                                        // trilinear

    // standard smagorinsky model
    if(les->model == SMAGORINSKY)
    {
        VecSet(les->lCsk, std_cs);
        DMLocalToLocalBegin(da, les->lCsk, INSERT_VALUES, les->lCsk);
        DMLocalToLocalEnd  (da, les->lCsk, INSERT_VALUES, les->lCsk);
        return(0);
    }
    else if(les->model == AMD || les->model == BAMD)
    {
        VecSet(les->lCsk, les->amdCs);
        DMLocalToLocalBegin(da, les->lCsk, INSERT_VALUES, les->lCsk);
        DMLocalToLocalEnd  (da, les->lCsk, INSERT_VALUES, les->lCsk);
        return(0);
    }
    else if(les->model == VREMAN || les->model == BV)
    {
        VecSet(les->lCsk, vreman_cs);
        DMLocalToLocalBegin(da, les->lCsk, INSERT_VALUES, les->lCsk);
        DMLocalToLocalEnd  (da, les->lCsk, INSERT_VALUES, les->lCsk);
        return(0);
    }
    else if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD)
    {
        //compute temperature gradients 
        // get jacobians
        DMDAVecGetArray(da, mesh->lAj,  &aj);

        // get face area vectors
        DMDAVecGetArray(fda, mesh->lCsi, &csi);
        DMDAVecGetArray(fda, mesh->lEta, &eta);
        DMDAVecGetArray(fda, mesh->lZet, &zet);
        DMDAVecGetArray(da,  mesh->lNvert, &nvert);
        DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
        DMDAVecGetArray(da, teqn->lTmprt, &lt);

        DMDAVecGetArray(da,  les->ldTheta, &dtheta);

        DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
        DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
        DMDAVecGetArray(fda, mesh->lCent,  &cent);

        // 1 - loop over internal cells and compute temperature gradient
        for (k=lzs; k<lze; k++)
        {
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // if on body skip
                    if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                    {
                        continue;
                    }

                    // get cell metrics
                    ajc = aj[k][j][i];
                    PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
                    PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
                    PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

                    PetscReal dtdc, dtde, dtdz;
                    PetscReal dt_dx, dt_dy, dt_dz;

                    Compute_dscalar_center
                    (
                        mesh,
                        i, j, k,
                        mx, my, mz,
                        lt, nvert, meshTag,
                        &dtdc, &dtde, &dtdz
                    );
        
                    // compute temperature derivative w.r.t. the cartesian coordinates
                    Compute_dscalar_dxyz
                    (
                        mesh,
                        csi0, csi1, csi2,
                        eta0, eta1, eta2,
                        zet0, zet1, zet2,
                        ajc,
                        dtdc, dtde, dtdz,
                        &dt_dx, &dt_dy, &dt_dz
                    );

                    dtheta[k][j][i].x = dt_dx; dtheta[k][j][i].y = dt_dy; dtheta[k][j][i].z = dt_dz;
                }
            }
        }

        DMDAVecRestoreArray(da,  les->ldTheta, &dtheta);
        DMLocalToLocalBegin(da, les->ldTheta, INSERT_VALUES, les->ldTheta);
        DMLocalToLocalEnd  (da, les->ldTheta, INSERT_VALUES, les->ldTheta);

        // get arrays
        DMDAVecGetArray(da,  les->lS,  &Sabs);
        DMDAVecGetArray(fda, les->lSx, &Ax);
        DMDAVecGetArray(fda, les->lSy, &Ay);
        DMDAVecGetArray(fda, les->lSz, &Az);
        DMDAVecGetArray(da,  les->ldTheta, &dtheta);
        DMDAVecGetArray(da,  les->lCsk, &csk);

        VecSet(les->lKX, 0.0);
        VecSet(les->lXX, 0.0);
        DMDAVecGetArray(da,  les->lKX, &KX);
        DMDAVecGetArray(da,  les->lXX, &XX);

        //first test filter - compute Ki and Xi
        for (k=lzs; k<lze; k++)
        for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++)
        {
            // set to zero if solid and skip
            if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                KX[k][j][i] = XX[k][j][i] = 0;
                continue;
            }

            // 2.1 - Build the variables to be test-filtered

            // get cell metrics
            ajc = aj[k][j][i];
            PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
            PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
            PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

            // indices for 2-nd rank tensor labeling
            PetscInt    a, b;

            // local tensors after the test-filtering
            PetscReal   Ki[3],                                   
                        S_hat,
                        dtheta_hat[3],
                        Sij_hat[3][3],                                  
                        Sdt_hat[3],                              
                        Xi[3],                                   
                        Ki_cat[3],                               
                        Xi_cat[3];                               

            PetscReal filter, test_filter;                         // grid filter and dynamic model filter

            // local tensors before the test-filtering (3x3x3 cell stencil)
            PetscReal S[3][3][3];                                  
            PetscReal u[3][3][3] ,   v[3][3][3] ,   w[3][3][3]   ;
            PetscReal U[3][3][3]   , V[3][3][3]   , W[3][3][3]   ; 
            PetscReal Sdt1[3][3][3], Sdt2[3][3][3], Sdt3[3][3][3];
            PetscReal S11[3][3][3] , S12[3][3][3] , S13[3][3][3] , // strain rate tensor
                      S21[3][3][3] , S22[3][3][3] , S23[3][3][3] ,
                      S31[3][3][3] , S32[3][3][3] , S33[3][3][3] ;
            PetscReal Utheta[3][3][3]  , Vtheta[3][3][3]  , Wtheta[3][3][3]  ;
            PetscReal theta[3][3][3]  ;
            PetscReal dthetax[3][3][3] ,   dthetay[3][3][3] ,   dthetaz[3][3][3]   ;    


            // creating variables to be test-filtered
            PetscInt p,q,r;
            for(p=-1; p<=1; p++)
            for(q=-1; q<=1; q++)
            for(r=-1; r<=1; r++)
            {
                PetscInt R = r+1, Q = q+1, P = p+1;                  // local stencil labeling: putting data
                PetscInt K = k+r, J = j+q, I = i+p;                  // mesh cell labeling    : getting data

                // //mirroring ghost node values
                // if(K == 0) K = 1;
                // if(K == mz-1) K = mz-2;
                // if(J == 0) J = 1;
                // if(J == my-1) J = my-2;
                // if(I == 0) I = 1;
                // if(I == mx-1) I = mx-2;

                // cartesian velocity
                u[R][Q][P] = ucat[K][J][I].x;
                v[R][Q][P] = ucat[K][J][I].y;
                w[R][Q][P] = ucat[K][J][I].z;

                // contravariant fluxes
                U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
                V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
                W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;

                theta[R][Q][P] = lt[K][J][I];

                Utheta[R][Q][P] = U[R][Q][P] * theta[R][Q][P];
                Vtheta[R][Q][P] = V[R][Q][P] * theta[R][Q][P];
                Wtheta[R][Q][P] = W[R][Q][P] * theta[R][Q][P];

                dthetax[R][Q][P] = dtheta[K][J][I].x;
                dthetay[R][Q][P] = dtheta[K][J][I].y;
                dthetaz[R][Q][P] = dtheta[K][J][I].z;

                // strain rate tensor S
                const PetscReal du_dx = Ax[K][J][I].x,
                                du_dy = Ax[K][J][I].y,
                                du_dz = Ax[K][J][I].z;
                const PetscReal dv_dx = Ay[K][J][I].x,
                                dv_dy = Ay[K][J][I].y,
                                dv_dz = Ay[K][J][I].z;
                const PetscReal dw_dx = Az[K][J][I].x,
                                dw_dy = Az[K][J][I].y,
                                dw_dz = Az[K][J][I].z;

                PetscReal Sxx = 0.5*(du_dx + du_dx),
                          Sxy = 0.5*(du_dy + dv_dx),
                          Sxz = 0.5*(du_dz + dw_dx);

                PetscReal Syx = Sxy,
                          Syy = 0.5*(dv_dy + dv_dy),
                          Syz = 0.5*(dv_dz + dw_dy);
                PetscReal Szx = Sxz,
                          Szy=Syz,
                          Szz = 0.5*(dw_dz + dw_dz);

                S11[R][Q][P]  = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
                S21[R][Q][P]  = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
                S31[R][Q][P]  = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;

                // |S|
                S[R][Q][P]    = Sabs[K][J][I];

                // |S|*dtheta
                Sdt1[R][Q][P] = dthetax[R][Q][P]*S[R][Q][P], Sdt2[R][Q][P] = dthetay[R][Q][P]*S[R][Q][P], Sdt3[R][Q][P] = dthetaz[R][Q][P]*S[R][Q][P];
            }

            // 2.2 - Build the filters
            // binomial filter coefficient n=2, m = 1/4, d= [1 2 1] (from pascal triangle)
            // binomial filter of order n is given by f = m[d]

            // filter kernel
            PetscReal coef[3][3][3] = 
            {
                0.015625, 0.03125, 0.015625,
                0.03125, 0.0625, 0.03125,
                0.015625, 0.03125, 0.015625,
            
            
                0.03125, 0.0625, 0.03125,
                0.0625, 0.125, 0.0625,
                0.03125, 0.0625, 0.03125,
            
            
                0.015625, 0.03125, 0.015625,
                0.03125, 0.0625, 0.03125,
                0.015625, 0.03125, 0.015625
            };

            // cell volume times stencil integration coeff. and fluid/body weight
            PetscReal weight[3][3][3];
            PetscReal sum_vol    = 0;

            for(p=-1; p<=1; p++)
            for(q=-1; q<=1; q++)
            for(r=-1; r<=1; r++)
            {
                PetscInt R = r+1, Q = q+1, P = p+1;                  // local stencil labeling: putting data
                PetscInt K = k+r, J = j+q, I = i+p;                  // mesh cell labeling    : getting data

                // //mirroring ghost node values
                // if(K == 0) K = 1;
                // if(K == mz-1) K = mz-2;
                // if(J == 0) J = 1;
                // if(J == my-1) J = my-2;
                // if(I == 0) I = 1;
                // if(I == mx-1) I = mx-2;

                // fluid
                if
                (
                    ((isFluidCell(K, J, I, nvert) || isIBMFluidCell(K, J, I, nvert)) && isCalculatedCell(K, J, I, meshTag)) &&
                    //!isOnCornerCellCenters(I, J, K, mesh->info) <- use also ghost? Not for now
                    (I!=0 && I!=mx-1 && J!=0 && J!=my-1 && K!=0 && K!=mz-1)
                )
                {
                    sum_vol += 1./aj[K][J][I] * 8.0 * coef[R][Q][P];
                    weight[R][Q][P] = 1;
                }
                // body or physical ghost nodes
                else
                {
                    // note: we could use physical ghost, but not certanely the corner points,
                    // where variables is undefined
                    weight[R][Q][P] = 0;
                }
            }

            // filter dimension is cell volume cubic root
            filter      = pow( 1./aj[k][j][i], 1./3.);

            // use the cubic root of total weighted volume
            test_filter = pow( sum_vol, 1./3. );

            // apply test-filtering
            PetscReal _U = integrateTestfilterSimpson(U, weight);
            PetscReal _V = integrateTestfilterSimpson(V, weight);
            PetscReal _W = integrateTestfilterSimpson(W, weight);

            PetscReal _theta = integrateTestfilterSimpson(theta, weight);
            
            dtheta_hat[0] = integrateTestfilterSimpson(dthetax, weight);
            dtheta_hat[1] = integrateTestfilterSimpson(dthetay, weight);
            dtheta_hat[2] = integrateTestfilterSimpson(dthetaz, weight);

            Ki[0] = integrateTestfilterSimpson(Utheta, weight) - _U*_theta;
            Ki[1] = integrateTestfilterSimpson(Vtheta, weight) - _V*_theta;
            Ki[2] = integrateTestfilterSimpson(Wtheta, weight) - _W*_theta;

            Sij_hat[0][0] = integrateTestfilterSimpson(S11, weight);
            Sij_hat[0][1] = integrateTestfilterSimpson(S12, weight);
            Sij_hat[0][2] = integrateTestfilterSimpson(S13, weight);
            Sij_hat[1][0] = integrateTestfilterSimpson(S21, weight);
            Sij_hat[1][1] = integrateTestfilterSimpson(S22, weight);
            Sij_hat[1][2] = integrateTestfilterSimpson(S23, weight);
            Sij_hat[2][0] = integrateTestfilterSimpson(S31, weight);
            Sij_hat[2][1] = integrateTestfilterSimpson(S32, weight);
            Sij_hat[2][2] = integrateTestfilterSimpson(S33, weight);

            S_hat = 0;

            for(a=0; a<3; a++)
            {
                for(b=0; b<3; b++)
                {
                    S_hat += pow( Sij_hat[a][b], 2. );
                }
            }

            S_hat = sqrt ( 2 * S_hat );

            Sdt_hat[0] = integrateTestfilterSimpson(Sdt1, weight);
            Sdt_hat[1] = integrateTestfilterSimpson(Sdt2, weight);
            Sdt_hat[2] = integrateTestfilterSimpson(Sdt3, weight);

            PetscReal gg[3][3], ggc[3][3], G[3][3];
            PetscReal   xcsi, xeta, xzet,
                        ycsi, yeta, yzet,
                        zcsi, zeta, zzet;

            // contravariant base matrix
            gg[0][0] = csi0, gg[0][1] = csi1, gg[0][2] = csi2;
            gg[1][0] = eta0, gg[1][1] = eta1, gg[1][2] = eta2;
            gg[2][0] = zet0, gg[2][1] = zet1, gg[2][2] = zet2;

            // covariant base matrix
            Calculate_Covariant_metrics(gg, ggc);

            // covariant base vectors
            xcsi = ggc[0][0], xeta = ggc[0][1], xzet = ggc[0][2];
            ycsi = ggc[1][0], yeta = ggc[1][1], yzet = ggc[1][2];
            zcsi = ggc[2][0], zeta = ggc[2][1], zzet = ggc[2][2];

            // covariant metric tensor
            G[0][0] = xcsi * xcsi + ycsi * ycsi + zcsi * zcsi;
            G[1][1] = xeta * xeta + yeta * yeta + zeta * zeta;
            G[2][2] = xzet * xzet + yzet * yzet + zzet * zzet;
            G[0][1] = G[1][0] = xeta * xcsi + yeta * ycsi + zeta * zcsi;
            G[0][2] = G[2][0] = xzet * xcsi + yzet * ycsi + zzet * zcsi;
            G[1][2] = G[2][1] = xeta * xzet + yeta * yzet + zeta * zzet;

            // compute cartesian Xi vector
            for(a=0; a<3; a++)
            {
                // Xi first defined in cartesian and then transformed
                Xi_cat[a]
                =
                (
                    - pow( test_filter, 2. ) * S_hat * dtheta_hat[a]
                    + pow( filter, 2. ) * Sdt_hat[a]
                );
            }

            // project the Xi
            Xi[0] = Xi_cat[0] * csi0 + Xi_cat[1] * csi1 + Xi_cat[2] * csi2;
            Xi[1] = Xi_cat[0] * eta0 + Xi_cat[1] * eta1 + Xi_cat[2] * eta2;
            Xi[2] = Xi_cat[0] * zet0 + Xi_cat[1] * zet1 + Xi_cat[2] * zet2;

            PetscReal num = 0, num1 = 0, denom = 0;
            PetscInt  m, n, l;

            // compute numerator
            for(b=0; b<3; b++)
            for(q=0; q<3; q++)
            {

                num += Ki[b] * Xi[q] * G[b][q];
            }

            // compute denominator
            for(m=0; m<3; m++)
            for(l=0; l<3; l++)
            {

                denom += Xi[m] * Xi[l] * G[m][l];
            }

            // store in KX and XX
            KX[k][j][i] = num;
            XX[k][j][i] = denom;

        }

        DMDAVecRestoreArray(da, les->lKX, &KX);
        DMDAVecRestoreArray(da, les->lXX, &XX);

        // update ghost nodes
        DMLocalToLocalBegin(da, les->lKX, INSERT_VALUES, les->lKX);
        DMLocalToLocalEnd  (da, les->lKX, INSERT_VALUES, les->lKX);
        DMLocalToLocalBegin(da, les->lXX, INSERT_VALUES, les->lXX);
        DMLocalToLocalEnd  (da, les->lXX, INSERT_VALUES, les->lXX);

        if(les->model == DLASD || les->model == DPASD)
        {
            VecSet(les->lPY, 0.0);
            VecSet(les->lYY, 0.0);
            DMDAVecGetArray(da,  les->lPY, &PY);
            DMDAVecGetArray(da,  les->lYY, &YY);

            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                // set to zero if solid and skip
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    PY[k][j][i] = YY[k][j][i] = 0;
                    continue;
                }

                // 2.1 - Build the variables to be test-filtered

                // get cell metrics
                ajc = aj[k][j][i];
                PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
                PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
                PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

                // indices for 2-nd rank tensor labeling
                PetscInt    a, b;

                // local tensors after the test-filtering
                PetscReal   Pi[3],                                   
                            S_hat,
                            dtheta_hat[3],
                            Sij_hat[3][3],                                  
                            Sdt_hat[3],                              
                            Yi[3],                                   
                            Pi_cat[3],                               
                            Yi_cat[3];                               

                PetscReal filter, test_filter;                         // grid filter and dynamic model filter

                // local tensors before the test-filtering (5x5x5 cell stencil)
                PetscReal S[5][5][5];                                  // scalar |S| = sqrt(2*Sij*Sij)
                PetscReal u[5][5][5] ,   v[5][5][5] ,   w[5][5][5]   ; // cartesian velocity   
                PetscReal U[5][5][5]   , V[5][5][5]   , W[5][5][5]   ; // contravariant fluxes
                PetscReal Sdt1[5][5][5], Sdt2[5][5][5], Sdt3[5][5][5];
                PetscReal S11[5][5][5] , S12[5][5][5] , S13[5][5][5] , // strain rate tensor
                          S21[5][5][5] , S22[5][5][5] , S23[5][5][5],
                          S31[5][5][5] , S32[5][5][5] , S33[5][5][5] ;
                PetscReal Utheta[5][5][5]  , Vtheta[5][5][5]  , Wtheta[5][5][5]  ;
                PetscReal theta[5][5][5]  ;
                PetscReal dthetax[5][5][5] ,   dthetay[5][5][5] ,   dthetaz[5][5][5]   ;  

                // creating variables to be test-filtered
                PetscInt p,q,r;
                for(p=-2; p<=2; p++)
                for(q=-2; q<=2; q++)
                for(r=-2; r<=2; r++)
                {
                    PetscInt R = r+2, Q = q+2, P = p+2;                  // local stencil labeling: putting data
                    PetscInt K = k+r, J = j+q, I = i+p;                  // mesh cell labeling    : getting data

                    // // ghost node values are mirrored
                    // if(K <= 0) K = 1;
                    // if(K >= mz-1) K = mz-2;
                    // if(J <= 0) J = 1;
                    // if(J >= my-1) J = my-2;
                    // if(I <= 0) I = 1;
                    // if(I >= mx-1) I = mx-2; 

                    // cartesian velocity
                    u[R][Q][P] = ucat[K][J][I].x;
                    v[R][Q][P] = ucat[K][J][I].y;
                    w[R][Q][P] = ucat[K][J][I].z;

                    // contravariant fluxes
                    U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
                    V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
                    W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;

                    theta[R][Q][P] = lt[K][J][I];

                    Utheta[R][Q][P] = U[R][Q][P] * theta[R][Q][P];
                    Vtheta[R][Q][P] = V[R][Q][P] * theta[R][Q][P];
                    Wtheta[R][Q][P] = W[R][Q][P] * theta[R][Q][P];

                    dthetax[R][Q][P] = dtheta[K][J][I].x;
                    dthetay[R][Q][P] = dtheta[K][J][I].y;
                    dthetaz[R][Q][P] = dtheta[K][J][I].z;

                    // strain rate tensor S
                    const PetscReal du_dx = Ax[K][J][I].x,
                                    du_dy = Ax[K][J][I].y,
                                    du_dz = Ax[K][J][I].z;
                    const PetscReal dv_dx = Ay[K][J][I].x,
                                    dv_dy = Ay[K][J][I].y,
                                    dv_dz = Ay[K][J][I].z;
                    const PetscReal dw_dx = Az[K][J][I].x,
                                    dw_dy = Az[K][J][I].y,
                                    dw_dz = Az[K][J][I].z;

                    PetscReal Sxx = 0.5*(du_dx + du_dx),
                            Sxy = 0.5*(du_dy + dv_dx),
                            Sxz = 0.5*(du_dz + dw_dx);

                    PetscReal Syx = Sxy,
                            Syy = 0.5*(dv_dy + dv_dy),
                            Syz = 0.5*(dv_dz + dw_dy);
                    PetscReal Szx = Sxz,
                            Szy=Syz,
                            Szz = 0.5*(dw_dz + dw_dz);

                    S11[R][Q][P]  = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
                    S21[R][Q][P]  = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
                    S31[R][Q][P]  = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;

                    // |S|
                    S[R][Q][P]    = Sabs[K][J][I];

                    // |S|*dtheta
                    Sdt1[R][Q][P] = dthetax[R][Q][P]*S[R][Q][P], Sdt2[R][Q][P] = dthetay[R][Q][P]*S[R][Q][P], Sdt3[R][Q][P] = dthetaz[R][Q][P]*S[R][Q][P];
                }

                // 2.2 - Build the filters
                // binomial filter coefficient n=4, m = 1/16, d= [1 4 6 4 1] (from pascal triangle)
                // binomial filter of order n is given by f = m[d]

                // filter kernel
                PetscReal coef[5][5][5]=
                {
                    0.000244140625000, 0.000976562500000, 0.001464843750000, 0.000976562500000, 0.000244140625000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.001464843750000, 0.005859375000000, 0.008789062500000, 0.005859375000000, 0.001464843750000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.000244140625000, 0.000976562500000, 0.001464843750000, 0.000976562500000, 0.000244140625000,

                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.003906250000000, 0.015625000000000, 0.023437500000000, 0.015625000000000, 0.003906250000000,
                    0.005859375000000, 0.023437500000000, 0.035156250000000, 0.023437500000000, 0.005859375000000,
                    0.003906250000000, 0.015625000000000, 0.023437500000000, 0.015625000000000, 0.003906250000000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,

                    0.001464843750000, 0.005859375000000, 0.008789062500000, 0.005859375000000, 0.001464843750000,
                    0.005859375000000, 0.023437500000000, 0.035156250000000, 0.023437500000000, 0.005859375000000,
                    0.008789062500000, 0.035156250000000, 0.052734375000000, 0.035156250000000, 0.008789062500000,
                    0.005859375000000, 0.023437500000000, 0.035156250000000, 0.023437500000000, 0.005859375000000,
                    0.001464843750000, 0.005859375000000, 0.008789062500000, 0.005859375000000, 0.001464843750000,

                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.003906250000000, 0.015625000000000, 0.023437500000000, 0.015625000000000, 0.003906250000000,
                    0.005859375000000, 0.023437500000000, 0.035156250000000, 0.023437500000000, 0.005859375000000,
                    0.003906250000000, 0.015625000000000, 0.023437500000000, 0.015625000000000, 0.003906250000000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,

                    0.000244140625000, 0.000976562500000, 0.001464843750000, 0.000976562500000, 0.000244140625000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.001464843750000, 0.005859375000000, 0.008789062500000, 0.005859375000000, 0.001464843750000,
                    0.000976562500000, 0.003906250000000, 0.005859375000000, 0.003906250000000, 0.000976562500000,
                    0.000244140625000, 0.000976562500000, 0.001464843750000, 0.000976562500000, 0.000244140625000
                };

                // cell volume times stencil integration coeff. and fluid/body weight
                PetscReal weight[5][5][5];
                PetscReal sum_vol    = 0;

                for(p=-2; p<=2; p++)
                for(q=-2; q<=2; q++)
                for(r=-2; r<=2; r++)
                {
                    PetscInt R = r+2, Q = q+2, P = p+2;                  // local stencil labeling: putting data
                    PetscInt K = k+r, J = j+q, I = i+p;                  // mesh cell labeling    : getting data

                    // //ghost node values are mirrored
                    // if(K <= 0) K = 1;
                    // if(K >= mz-1) K = mz-2;
                    // if(J <= 0) J = 1;
                    // if(J >= my-1) J = my-2;
                    // if(I <= 0) I = 1;
                    // if(I >= mx-1) I = mx-2; 

                    // fluid
                    if
                    (
                        ((isFluidCell(K, J, I, nvert) || isIBMFluidCell(K, J, I, nvert)) && isCalculatedCell(K, J, I, meshTag)) &&
                        //!isOnCornerCellCenters(I, J, K, mesh->info) <- use also ghost? Not for now
                        (I>0 && I<mx-1 && J>0 && J<my-1 && K>0 && K<mz-1)
                    )
                    {
                        sum_vol += 1./aj[K][J][I] * 64.0 * coef[R][Q][P];
                        weight[R][Q][P] = 1;
                    }
                    // body or physical ghost nodes
                    else
                    {
                        // note: we could use physical ghost, but not certanely the corner points,
                        // where variables is undefined
                        weight[R][Q][P] = 0;
                    }
                }

                // filter dimension is cell volume cubic root
                filter      = pow( 1./aj[k][j][i], 1./3.);

                // use the cubic root of total weighted volume
                test_filter = pow( sum_vol, 1./3. );

                // apply test-filtering
                PetscReal _U = integrateTestfilterSimpson5x5(U, weight);
                PetscReal _V = integrateTestfilterSimpson5x5(V, weight);
                PetscReal _W = integrateTestfilterSimpson5x5(W, weight);

                PetscReal _theta = integrateTestfilterSimpson5x5(theta, weight);

                Pi[0] = integrateTestfilterSimpson5x5(Utheta, weight) - _U*_theta;
                Pi[1] = integrateTestfilterSimpson5x5(Vtheta, weight) - _V*_theta;
                Pi[2] = integrateTestfilterSimpson5x5(Wtheta, weight) - _W*_theta;

                Sij_hat[0][0] = integrateTestfilterSimpson5x5(S11, weight);
                Sij_hat[0][1] = integrateTestfilterSimpson5x5(S12, weight);
                Sij_hat[0][2] = integrateTestfilterSimpson5x5(S13, weight);
                Sij_hat[1][0] = integrateTestfilterSimpson5x5(S21, weight);
                Sij_hat[1][1] = integrateTestfilterSimpson5x5(S22, weight);
                Sij_hat[1][2] = integrateTestfilterSimpson5x5(S23, weight);
                Sij_hat[2][0] = integrateTestfilterSimpson5x5(S31, weight);
                Sij_hat[2][1] = integrateTestfilterSimpson5x5(S32, weight);
                Sij_hat[2][2] = integrateTestfilterSimpson5x5(S33, weight);

                S_hat = 0;

                for(a=0; a<3; a++)
                {
                    for(b=0; b<3; b++)
                    {
                        S_hat += pow( Sij_hat[a][b], 2. );
                    }
                }

                S_hat = sqrt ( 2 * S_hat );

                Sdt_hat[0] = integrateTestfilterSimpson5x5(Sdt1, weight);
                Sdt_hat[1] = integrateTestfilterSimpson5x5(Sdt2, weight);
                Sdt_hat[2] = integrateTestfilterSimpson5x5(Sdt3, weight);

                PetscReal gg[3][3], ggc[3][3], G[3][3];
                PetscReal xcsi, xeta, xzet,
                        ycsi, yeta, yzet,
                        zcsi, zeta, zzet;

                // contravariant base matrix
                gg[0][0] = csi0, gg[0][1] = csi1, gg[0][2] = csi2;
                gg[1][0] = eta0, gg[1][1] = eta1, gg[1][2] = eta2;
                gg[2][0] = zet0, gg[2][1] = zet1, gg[2][2] = zet2;

                // covariant base matrix
                Calculate_Covariant_metrics(gg, ggc);

                // covariant base vectors
                xcsi = ggc[0][0], xeta = ggc[0][1], xzet = ggc[0][2];
                ycsi = ggc[1][0], yeta = ggc[1][1], yzet = ggc[1][2];
                zcsi = ggc[2][0], zeta = ggc[2][1], zzet = ggc[2][2];

                // covariant metric tensor
                G[0][0] = xcsi * xcsi + ycsi * ycsi + zcsi * zcsi;
                G[1][1] = xeta * xeta + yeta * yeta + zeta * zeta;
                G[2][2] = xzet * xzet + yzet * yzet + zzet * zzet;
                G[0][1] = G[1][0] = xeta * xcsi + yeta * ycsi + zeta * zcsi;
                G[0][2] = G[2][0] = xzet * xcsi + yzet * ycsi + zzet * zcsi;
                G[1][2] = G[2][1] = xeta * xzet + yeta * yzet + zeta * zzet;

                // compute cartesian Yi tensor
                for(a=0; a<3; a++)
                {
                    // Xi first defined in cartesian and then transformed
                    Yi_cat[a]
                    =
                    (
                        - pow( test_filter, 2. ) * S_hat * dtheta_hat[a]
                        + pow( filter, 2. ) * Sdt_hat[a]
                    );
                }

                // project the Xi
                Yi[0] = Yi_cat[0] * csi0 + Yi_cat[1] * csi1 + Yi_cat[2] * csi2;
                Yi[1] = Yi_cat[0] * eta0 + Yi_cat[1] * eta1 + Yi_cat[2] * eta2;
                Yi[2] = Yi_cat[0] * zet0 + Yi_cat[1] * zet1 + Yi_cat[2] * zet2;

                PetscReal num = 0, num1 = 0, denom = 0;
                PetscInt  m, n, l;

                // compute numerator
                for(b=0; b<3; b++)
                for(q=0; q<3; q++)
                {

                    num += Pi[b] * Yi[q] * G[b][q];
                }

                // compute denominator
                for(m=0; m<3; m++)
                for(l=0; l<3; l++)
                {

                    denom += Yi[m] * Yi[l] * G[m][l];
                }

                // store in KX and XX
                PY[k][j][i] = num;
                YY[k][j][i] = denom;
            }

            DMDAVecRestoreArray(da,  les->lPY, &PY);
            DMDAVecRestoreArray(da,  les->lYY, &YY);

            DMLocalToLocalBegin(da, les->lPY, INSERT_VALUES, les->lPY);
            DMLocalToLocalEnd  (da, les->lPY, INSERT_VALUES, les->lPY);
            DMLocalToLocalBegin(da, les->lYY, INSERT_VALUES, les->lYY);
            DMLocalToLocalEnd  (da, les->lYY, INSERT_VALUES, les->lYY);
        }

        PetscReal ***KX_old, ***XX_old, ***PY_old, ***YY_old;

        // langrangian - local convective weighting
        if(les->model == DLASI || les->model == DLASD)
        {
            PetscInt initializeLes4;
            if (clock->it <= clock->itStart + 1) initializeLes4 = 1;
            else                                 initializeLes4 = 0;

            DMDAVecGetArray(da, les->lKX, &KX);
            DMDAVecGetArray(da, les->lXX, &XX);

            DMDAVecGetArray(da, les->lKX_old, &KX_old);
            DMDAVecGetArray(da, les->lXX_old, &XX_old);

            // use standard Smagorinsky for first iteration (since need old values of KX and XX)
            if(initializeLes4)
            {
                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            KX_old[k][j][i] = std_cs*XX[k][j][i];
                            XX_old[k][j][i] = XX[k][j][i];
                        }
                    }
                }

                VecSet(les->lCsk, std_cs);
            }
            else
            {
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    // filter dimension
                    PetscReal filter = pow( 1./aj[k][j][i], 1./3. );

                    // The higher the time scale, the strongest the lagrangian scmoothing.
                    // To low time scales could result in weighting too much the new
                    // value, which in turn could be clipped at zero. Hence we have noisier fields.

                    // Bou-Zaid scale dependent model time scale
                    PetscReal T_scale
                    =
                    1.5*filter*
                    pow
                    (
                        fabs(KX_old[k][j][i]*XX_old[k][j][i]) +
                        1.e-19,
                        -1.0/8.0
                    ) +
                    1.e-19;
                    

                    PetscReal KX_new, XX_new;

                    KX_new = KX[k][j][i];
                    XX_new = XX[k][j][i];

                    // find the coordinate where old values will be gahtered
                    Cmpnts X_new, X_old;

                    X_new = cent[k][j][i];

                    X_old.x = X_new.x - ucat[k][j][i].x * clock->dt;
                    X_old.y = X_new.y - ucat[k][j][i].y * clock->dt;
                    X_old.z = X_new.z - ucat[k][j][i].z * clock->dt;

                    PetscInt    i1, j1, k1;
                    PetscReal dmin = 10e6, d;
                    PetscInt    i_old, j_old, k_old;

                    i_old = i; j_old = j; k_old = k;

                    // must be close due to CFL so don't loop over all cells
                    for (k1=k-2; k1<k+3; k1++)
                    for (j1=j-2; j1<j+3; j1++)
                    for (i1=i-2; i1<i+3; i1++)
                    {
                        // clip the search on internal cells. Must update values if
                        // periodic after that otherwise inlet values are not advected
                        // from outlet
                        if
                        ((
                            k1>=1 && k1<mz-1 &&
                            j1>=1 && j1<my-1 &&
                            i1>=1 && i1<mx-1
                        ) && ( !(isIBMSolidCell(k1, j1, i1, nvert) || isZeroedCell(k1, j1, i1, meshTag)) ) )
                        {
                            d = pow((X_old.x - cent[k1][j1][i1].x), 2) +
                                pow((X_old.y - cent[k1][j1][i1].y), 2) +
                                pow((X_old.z - cent[k1][j1][i1].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                i_old = i1;
                                j_old = j1;
                                k_old = k1;
                            }
                        }
                    }

                    PetscReal _KX_old, _XX_old;

                    // interpolate the value at the exat old point (slower)
                    // if there is an IBM or overset cell around the given cell, use nearest cell
                    if(isBoxIBMCell(k_old, j_old, i_old, nvert))
                    {
                        PetscInt intId[6];
                        PetscInt ibmCellCtr = 0;

                        // get the trilinear interpolation cells
                        PointInterpolationCells
                        (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                intId
                        );

                        //see if any one of the trilinear interpolation cells is ibmSolid
                        for (PetscInt kk = 0; kk<2; kk++)
                        for (PetscInt jj = 0; jj<2; jj++)
                        for (PetscInt ii = 0; ii<2; ii++)
                        {
                            if(isIBMSolidCell(intId[kk], intId[jj+2], intId[ii+4], nvert) || isZeroedCell(intId[kk], intId[jj+2], intId[ii+4], meshTag))
                            {
                                ibmCellCtr ++;
                            }

                        }

                        if(ibmCellCtr)
                        {
                            //use the value at the closest cell
                            _KX_old = KX_old[k_old][j_old][i_old];
                            _XX_old = XX_old[k_old][j_old][i_old];
                        }
                        else
                        {
                            //safe to trilinear interpolate
                            scalarPointLocalVolumeInterpolation
                            (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                KX_old,
                                _KX_old
                            );

                            scalarPointLocalVolumeInterpolation
                            (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                XX_old,
                                _XX_old
                            );
                        }
                    }
                    else
                    {
                      if(lagrangianAverageCs == "trilinear")
                      {
                          scalarPointLocalVolumeInterpolation
                          (
                              mesh,
                              X_old.x, X_old.y, X_old.z,
                              i_old, j_old, k_old,
                              cent,
                              KX_old,
                              _KX_old
                          );

                          scalarPointLocalVolumeInterpolation
                          (
                              mesh,
                              X_old.x, X_old.y, X_old.z,
                              i_old, j_old, k_old,
                              cent,
                              XX_old,
                              _XX_old
                          );
                      }
                      else if(lagrangianAverageCs == "nearestPoint")
                      {
                          // take the closest cell value (faster, but a lot of cs clipping from below)
                          _KX_old = KX_old[k_old][j_old][i_old];
                          _XX_old = XX_old[k_old][j_old][i_old];
                      }

                    }

                    // interpolation weights depend on the time scale.
                    // if dt >> T_scale new value is used
                    // if dt << T_scale old value is used
                    PetscReal eps = (clock->dt/T_scale) / (1.0 + clock->dt/T_scale);

                    // apply lagrangian weighting (relaxation)
                    KX[k][j][i]
                    =
                    PetscMax
                    (
                        (
                            (1-eps) * _KX_old +
                            eps * KX_new
                        ),
                        0.0
                    );

                    XX[k][j][i]
                    =
                    (1-eps) * _XX_old +
                    eps * XX_new;
                }

                // store KX and XX old values
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    KX_old[k][j][i] = KX[k][j][i];
                    XX_old[k][j][i] = XX[k][j][i];
                }
            }

            DMDAVecRestoreArray(da, les->lKX, &KX);
            DMDAVecRestoreArray(da, les->lXX, &XX);

            DMDAVecRestoreArray(da, les->lKX_old, &KX_old);
            DMDAVecRestoreArray(da, les->lXX_old, &XX_old);

            DMLocalToLocalBegin(da, les->lKX_old, INSERT_VALUES, les->lKX_old);
            DMLocalToLocalEnd  (da, les->lKX_old, INSERT_VALUES, les->lKX_old);

            DMLocalToLocalBegin(da, les->lXX_old, INSERT_VALUES, les->lXX_old);
            DMLocalToLocalEnd  (da, les->lXX_old, INSERT_VALUES, les->lXX_old);

        }

        if(les->model == DLASD)
        {
            PetscInt initializeLes4;
            if (clock->it <= clock->itStart + 1) initializeLes4 = 1;
            else                                 initializeLes4 = 0;

            DMDAVecGetArray(da, les->lPY, &PY);
            DMDAVecGetArray(da, les->lYY, &YY);

            DMDAVecGetArray(da, les->lPY_old, &PY_old);
            DMDAVecGetArray(da, les->lYY_old, &YY_old);

            // use standard Smagorinsky for first iteration (since need old values of PY and YY)
            if(initializeLes4)
            {
                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            PY_old[k][j][i] = std_cs*YY[k][j][i];
                            YY_old[k][j][i] = YY[k][j][i];
                        }
                    }
                }

                VecSet(les->lCsk, std_cs);
            }
            else
            {
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    // filter dimension
                    PetscReal filter = pow( 1./aj[k][j][i], 1./3. );

                    // The higher the time scale, the strongest the lagrangian scmoothing.
                    // To low time scales could result in weighting too much the new
                    // value, which in turn could be clipped at zero. Hence we have noisier fields.

                    // Bou-Zaid scale dependent model time scale
                    PetscReal T_scale
                    =
                    1.5*filter*
                    pow
                    (
                        fabs(PY_old[k][j][i]*YY_old[k][j][i]) +
                        1.e-19,
                        -1.0/8.0
                    ) +
                    1.e-19;
                    

                    PetscReal PY_new, YY_new;

                    PY_new = PY[k][j][i];
                    YY_new = YY[k][j][i];

                    // find the coordinate where old values will be gahtered
                    Cmpnts X_new, X_old;

                    X_new = cent[k][j][i];

                    X_old.x = X_new.x - ucat[k][j][i].x * clock->dt;
                    X_old.y = X_new.y - ucat[k][j][i].y * clock->dt;
                    X_old.z = X_new.z - ucat[k][j][i].z * clock->dt;

                    PetscInt    i1, j1, k1;
                    PetscReal dmin = 10e6, d;
                    PetscInt    i_old, j_old, k_old;

                    i_old = i; j_old = j; k_old = k;

                    // must be close due to CFL so don't loop over all cells
                    for (k1=k-2; k1<k+3; k1++)
                    for (j1=j-2; j1<j+3; j1++)
                    for (i1=i-2; i1<i+3; i1++)
                    {
                        // clip the search on internal cells. Must update values if
                        // periodic after that otherwise inlet values are not advected
                        // from outlet
                        if
                        ((
                            k1>=1 && k1<mz-1 &&
                            j1>=1 && j1<my-1 &&
                            i1>=1 && i1<mx-1
                        ) && ( !(isIBMSolidCell(k1, j1, i1, nvert) || isZeroedCell(k1, j1, i1, meshTag)) ) )
                        {
                            d = pow((X_old.x - cent[k1][j1][i1].x), 2) +
                                pow((X_old.y - cent[k1][j1][i1].y), 2) +
                                pow((X_old.z - cent[k1][j1][i1].z), 2);

                            if
                            (
                                d < dmin
                            )
                            {
                                dmin  = d;
                                i_old = i1;
                                j_old = j1;
                                k_old = k1;
                            }
                        }
                    }

                    PetscReal _PY_old, _YY_old;

                    // interpolate the value at the exact old point (slower)
                    // if there is an IBM or overset cell around the given cell, use nearest cell
                    if(isBoxIBMCell(k_old, j_old, i_old, nvert))
                    {
                        PetscInt intId[6];
                        PetscInt ibmCellCtr = 0;

                        // get the trilinear interpolation cells
                        PointInterpolationCells
                        (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                intId
                        );

                        //see if any one of the trilinear interpolation cells is ibmSolid
                        for (PetscInt kk = 0; kk<2; kk++)
                        for (PetscInt jj = 0; jj<2; jj++)
                        for (PetscInt ii = 0; ii<2; ii++)
                        {
                            if(isIBMSolidCell(intId[kk], intId[jj+2], intId[ii+4], nvert) || isZeroedCell(intId[kk], intId[jj+2], intId[ii+4], meshTag))
                            {
                                ibmCellCtr ++;
                            }

                        }

                        if(ibmCellCtr)
                        {
                            //use the value at the closest cell
                            _PY_old = PY_old[k_old][j_old][i_old];
                            _YY_old = YY_old[k_old][j_old][i_old];
                        }
                        else
                        {
                            //safe to trilinear interpolate
                            scalarPointLocalVolumeInterpolation
                            (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                PY_old,
                                _PY_old
                            );

                            scalarPointLocalVolumeInterpolation
                            (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                YY_old,
                                _YY_old
                            );
                        }
                    }
                    else
                    {
                      if(lagrangianAverageCs == "trilinear")
                      {
                          scalarPointLocalVolumeInterpolation
                          (
                              mesh,
                              X_old.x, X_old.y, X_old.z,
                              i_old, j_old, k_old,
                              cent,
                              PY_old,
                              _PY_old
                          );

                          scalarPointLocalVolumeInterpolation
                          (
                              mesh,
                              X_old.x, X_old.y, X_old.z,
                              i_old, j_old, k_old,
                              cent,
                              YY_old,
                              _YY_old
                          );
                      }
                      else if(lagrangianAverageCs == "nearestPoint")
                      {
                          // take the closest cell value (faster, but a lot of cs clipping from below)
                          _PY_old = PY_old[k_old][j_old][i_old];
                          _YY_old = YY_old[k_old][j_old][i_old];
                      }

                    }

                    // interpolation weights depend on the time scale.
                    // if dt >> T_scale new value is used
                    // if dt << T_scale old value is used
                    PetscReal eps = (clock->dt/T_scale) / (1.0 + clock->dt/T_scale);

                    // apply lagrangian weighting (relaxation)
                    PY[k][j][i]
                    =
                    PetscMax
                    (
                        (
                            (1-eps) * _PY_old +
                            eps * PY_new
                        ),
                        0.0
                    );

                    YY[k][j][i]
                    =
                    (1-eps) * _YY_old +
                    eps * YY_new;
                }

                // store PY and YY old values
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    PY_old[k][j][i] = PY[k][j][i];
                    YY_old[k][j][i] = YY[k][j][i];
                }
            }

            DMDAVecRestoreArray(da, les->lPY, &PY);
            DMDAVecRestoreArray(da, les->lYY, &YY);
            DMDAVecRestoreArray(da, les->lPY_old, &PY_old);
            DMDAVecRestoreArray(da, les->lYY_old, &YY_old);

            DMLocalToLocalBegin(da, les->lPY_old, INSERT_VALUES, les->lPY_old);
            DMLocalToLocalEnd  (da, les->lPY_old, INSERT_VALUES, les->lPY_old);
            DMLocalToLocalBegin(da, les->lYY_old, INSERT_VALUES, les->lYY_old);
            DMLocalToLocalEnd  (da, les->lYY_old, INSERT_VALUES, les->lYY_old);

            // update ghost values
            DMLocalToLocalBegin(da, les->lPY, INSERT_VALUES, les->lPY);
            DMLocalToLocalEnd  (da, les->lPY, INSERT_VALUES, les->lPY);
            DMLocalToLocalBegin(da, les->lYY, INSERT_VALUES, les->lYY);
            DMLocalToLocalEnd  (da, les->lYY, INSERT_VALUES, les->lYY);

            resetCellPeriodicFluxes(mesh, les->lPY, les->lPY, "scalar", "localToLocal");
            resetCellPeriodicFluxes(mesh, les->lYY, les->lYY, "scalar", "localToLocal");
        }

        // update ghost values
        DMLocalToLocalBegin(da, les->lKX, INSERT_VALUES, les->lKX);
        DMLocalToLocalEnd  (da, les->lKX, INSERT_VALUES, les->lKX);
        DMLocalToLocalBegin(da, les->lXX, INSERT_VALUES, les->lXX);
        DMLocalToLocalEnd  (da, les->lXX, INSERT_VALUES, les->lXX);

        resetCellPeriodicFluxes(mesh, les->lKX, les->lKX, "scalar", "localToLocal");
        resetCellPeriodicFluxes(mesh, les->lXX, les->lXX, "scalar", "localToLocal");

        DMDAVecRestoreArray(da,  les->lS,  &Sabs);
        DMDAVecRestoreArray(fda, les->lSx, &Ax);
        DMDAVecRestoreArray(fda, les->lSy, &Ay);
        DMDAVecRestoreArray(fda, les->lSz, &Az);
        
        // 3 - Compute Csk

        if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD)
        {
            DMDAVecGetArray(da, les->lKX, &KX);
            DMDAVecGetArray(da, les->lXX, &XX);
    
            if(les->model == DLASD || les->model == DPASD)
            {
                DMDAVecGetArray(da,  les->lPY, &PY);
                DMDAVecGetArray(da,  les->lYY, &YY);
            }
        }

        // loop over internal cells
        for (k=lzs; k<lze; k++)
        for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++)
        {
            // body
            if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                csk[k][j][i] = 0.0;
                continue;
            }

            // Csk value at current point
            PetscReal C = 0.0;

            // dynamic Smagorinsky models
            if(les->model == DSM || les->model == DLASI)
            {
                PetscReal KX_avg, XX_avg;

                // box averaging
                if(les->model == DSM)
                {
                    PetscReal weight[3][3][3];

                    PetscReal KX0[3][3][3],
                              XX0[3][3][3];

                    PetscInt  a, b, c;

                    for(a=-1; a<=1; a++)
                    for(b=-1; b<=1; b++)
                    for(c=-1; c<=1; c++)
                    {
                        PetscInt R = c+1, Q = b+1, P = a+1;
                        PetscInt K = k+c, J = j+b, I = i+a;

                        // volume weighting the mean
                        weight[R][Q][P] = 1./aj[K][J][I];

                        if(isIBMSolidCell(K, J, I, nvert) || isZeroedCell(K, J, I, meshTag))
                        {
                            weight[R][Q][P] = 0;
                        }

                        if(mesh->i_periodic)
                        {
                            if(I==0) I=mx-2;
                            else if(I==mx-1) I = 1;
                        }
                        else if(mesh->ii_periodic)
                        {
                            if( I==0 ) I = -2;
                            else if( I==mx-1 ) I = mx+1;
                        }
                        else if( I==0 || I==mx-1)
                        {
                            // set to zero the weighting on physical ghost nodes
                            weight[R][Q][P] = 0;
                        }

                        if(mesh->j_periodic)
                        {
                            if( J==0 ) J = my-2;
                            else if( J==my-1 ) J = 1;
                        }
                        else if(mesh->jj_periodic)
                        {
                            if( J==0 ) J=-2;
                            else if( J==my-1 ) J = my+1;
                        }
                        else if( J==0 || j==my-1)
                        {
                            // set to zero the weighting on physical ghost nodes
                            weight[R][Q][P] = 0;
                        }

                        if(mesh->k_periodic)
                        {
                            if( K==0 ) K = mz-2;
                            else if( K==mz-1 ) K = 1;
                        }
                        else if(mesh->kk_periodic)
                        {
                            if( K==0 ) K=-2;
                            else if( K==mz-1 ) K=mz+1;
                        }
                        else if( K==0 || K==mz-1)
                        {
                            // set to zero the weighting on physical ghost nodes
                            weight[R][Q][P] = 0;
                        }

                        KX0[R][Q][P] = KX[K][J][I];
                        XX0[R][Q][P] = XX[K][J][I];
                    }

                    // integrate over the testfilter volume locally
                    KX_avg = integrateTestfilterSimpson(KX0, weight);
                    XX_avg = integrateTestfilterSimpson(XX0, weight);
                }
                // lagrangian averaging
                else if(les->model == DLASI)
                {
                    // no average procedure
                    KX_avg = KX[k][j][i];
                    XX_avg = XX[k][j][i];
                }

                // set Smagorinsky coefficient
                C = 0.5 * KX_avg / (XX_avg + 1e-10);
            }
            else if(les->model == DLASD)
            {
                PetscReal num, denom, beta;

                num   = KX[k][j][i]/(XX[k][j][i] + 1e-12);

                //ratio of Csk^2 at different scale
                beta  = (PY[k][j][i] * XX[k][j][i])/(YY[k][j][i] * KX[k][j][i] + 1e-12);

                denom = PetscMin(PetscMax(0.125, beta), 100.0);

                C = 0.5 * num/denom;
            }

            // set csk to zero if solid or on box corners
            if( isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                csk[k][j][i] = 0;
            }
            else
            {
                // clip csk between upper and lower bound
                csk[k][j][i] = PetscMin( PetscMax(C, 0.0), les->maxCs);
            }
        }

        if(les->model == DPASD)
        {
            std::vector<PetscInt> count, total_count;
            std::vector<PetscReal> J_KX(my), J_XX(my), J_PY(my), J_YY(my), KX_tmp(my), XX_tmp(my), PY_tmp(my), YY_tmp(my);
            
            count.resize(my);
            total_count.resize(my);
            
            for(j=0; j<my; j++) 
            {
                KX_tmp[j] = 0;
                XX_tmp[j] = 0;
                PY_tmp[j] = 0;
                YY_tmp[j] = 0;
                count[j]  = 0;
                total_count[j] = 0;
            }

            // loop over internal cells
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                if(isFluidCell(k,j,i, nvert)) 
                {
                    KX_tmp[j] += KX[k][j][i];
                    XX_tmp[j] += XX[k][j][i];
                    PY_tmp[j] += PY[k][j][i];
                    YY_tmp[j] += YY[k][j][i];            
                    count[j] ++;
                }
            }
            
            MPI_Allreduce( &KX_tmp[0], &J_KX[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
            MPI_Allreduce( &XX_tmp[0], &J_XX[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
            MPI_Allreduce( &PY_tmp[0], &J_PY[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
            MPI_Allreduce( &YY_tmp[0], &J_YY[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
            MPI_Allreduce( &count[0], &total_count[0], my, MPIU_INT, MPI_SUM, mesh->MESH_COMM);	
            
            for(j=0; j<my; j++) 
            {
                if( total_count[j]>0) 
                {
                    J_KX[j] = J_KX[j]/total_count[j];
                    J_XX[j] = J_XX[j]/total_count[j];
                    J_PY[j] = J_PY[j]/total_count[j];
                    J_YY[j] = J_YY[j]/total_count[j];
                }
            }
            
            PetscReal beta;
            // loop over internal cells
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                // body
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    csk[k][j][i] = 0.0;
                    continue;
                }

                beta = (J_PY[j] * J_XX[j])/(J_YY[j] * J_KX[j] +  1e-10);
                csk[k][j][i] = 0.5 * (J_KX[j] / ( J_XX[j] + 1e-10))/beta;

                // clip csk between upper and lower bound
                csk[k][j][i] = PetscMin( PetscMax(csk[k][j][i], 0.0), les->maxCs);
            }
        }

        if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD)
        {
            DMDAVecRestoreArray(da, les->lKX, &KX);
            DMDAVecRestoreArray(da, les->lXX, &XX);
    
            if(les->model == DLASD || les->model == DPASD)
            {
                DMDAVecRestoreArray(da,  les->lPY, &PY);
                DMDAVecRestoreArray(da,  les->lYY, &YY);
            }
        }

        // get face area vectors
        DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
        DMDAVecRestoreArray(fda, mesh->lEta, &eta);
        DMDAVecRestoreArray(fda, mesh->lZet, &zet);
        DMDAVecRestoreArray(da, mesh->lAj,  &aj);
        DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
        DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
        DMDAVecRestoreArray(da, teqn->lTmprt, &lt);

        DMDAVecRestoreArray(da,  les->ldTheta, &dtheta);
        DMDAVecRestoreArray(da,  les->lCsk, &csk);

        DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
        DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
        DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

        // scatter the model coefficient
        DMLocalToLocalBegin(da, les->lCsk, INSERT_VALUES, les->lCsk);
        DMLocalToLocalEnd  (da, les->lCsk, INSERT_VALUES, les->lCsk);

        // handle periodicity and scatter
        resetCellPeriodicFluxes(mesh, les->lCsk, les->lCsk, "scalar", "localToLocal");
    }

    return(0);
}

PetscErrorCode UpdatekT(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    ueqn_         *ueqn = les->access->ueqn;
    teqn_         *teqn = les->access->teqn;
    constants_    *cst  = teqn->access->constants;
    DM            da    = mesh->da, fda = mesh->fda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***csk, ***nvert, ***meshTag, ***aj, ***lkt, ***lt, ***nut;
    Cmpnts        ***csi, ***eta, ***zet, ***ucat;

    PetscReal     ajc;

    PetscReal     dudc, dude, dudz,
                  dvdc, dvde, dvdz,
                  dwdc, dwde, dwdz;

    const char*   lagrangianAverageCs = "trilinear";    // nearestPoint
                                                        // trilinear

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(les->lk_t, 0.);

    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da,  mesh->lAj, &aj);
    
    DMDAVecGetArray(da,  les->lCsk, &csk);
    DMDAVecGetArray(da,  les->lNu_t, &nut);
    DMDAVecGetArray(da, teqn->lTmprt, &lt);
    DMDAVecGetArray(da,  les->lk_t, &lkt);

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
        {
            lkt[k][j][i]=0;
            continue;
        }

        PetscReal ajc  = aj[k][j][i];
        PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

        PetscReal dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;
        PetscReal dtdc, dtde, dtdz;
        PetscReal dt_dx, dt_dy, dt_dz;
        PetscReal gradU[3][3], gradT[3];

        // compute velocity derivatives w.r.t the curvilinear coordinates
        if(ueqn->central4Div)
        {
            Compute_du_center4th
            (
                mesh,
                i, j, k,
                mx, my, mz,
                ucat, nvert, meshTag,
                &dudc, &dvdc, &dwdc,
                &dude, &dvde, &dwde,
                &dudz, &dvdz, &dwdz
            );
        }
        else
        { 
            Compute_du_center
            (
                mesh,
                i, j, k,
                mx, my, mz,
                ucat, nvert, meshTag,
                &dudc, &dvdc, &dwdc,
                &dude, &dvde, &dwde,
                &dudz, &dvdz, &dwdz
            );
        }

        // compute velocity derivative w.r.t. the cartesian coordinates
        Compute_du_dxyz
        (
            mesh,
            csi0, csi1, csi2,
            eta0, eta1, eta2,
            zet0, zet1, zet2,
            ajc,
            dudc, dvdc, dwdc,
            dude, dvde, dwde,
            dudz, dvdz, dwdz,
            &du_dx, &dv_dx, &dw_dx,
            &du_dy, &dv_dy, &dw_dy,
            &du_dz, &dv_dz, &dw_dz
        );

        PetscReal S[3][3] = {
            {0.5 * (du_dx + du_dx), 0.5 * (du_dy + dv_dx), 0.5 * (du_dz + dw_dx)},
            {0.5 * (dv_dx + du_dy), 0.5 * (dv_dy + dv_dy), 0.5 * (dv_dz + dw_dy)},
            {0.5 * (dw_dx + du_dz), 0.5 * (dw_dy + dv_dz), 0.5 * (dw_dz + dw_dz)}
        };

        // compute rate of strain tensor norm
        PetscReal Sabs = sqrt(2.0 * (
            S[0][0]*S[0][0] + S[0][1]*S[0][1] + S[0][2]*S[0][2] +
            S[1][0]*S[1][0] + S[1][1]*S[1][1] + S[1][2]*S[1][2] +
            S[2][0]*S[2][0] + S[2][1]*S[2][1] + S[2][2]*S[2][2]
        ));

        if (les->model == AMD || les->model == BAMD)
        {

            // Compute volume
            PetscReal V = 1.0 / aj[k][j][i];
            
            // Compute face areas
            PetscReal A_x = nMag(zet[k][j][i]);
            PetscReal A_y = nMag(csi[k][j][i]);
            PetscReal A_z = nMag(eta[k][j][i]); 
                                                    
            // Compute delta_x, delta_y, delta_z
            PetscReal delta_x = V / A_x;
            PetscReal delta_y = V / A_y;
            PetscReal delta_z = V / A_z;

            PetscReal invDelta2Sum =
                  1.0 / (delta_x * delta_x)
                + 1.0 / (delta_y * delta_y)
                + 1.0 / (delta_z * delta_z);

            PetscReal delta2_harm = 3.0 / invDelta2Sum; 

            PetscReal num = 0.0, denom = 0.0;

            gradU[0][0] = du_dx; gradU[0][1] = dv_dx; gradU[0][2] = dw_dx;
            gradU[1][0] = du_dy; gradU[1][1] = dv_dy; gradU[1][2] = dw_dy;
            gradU[2][0] = du_dz; gradU[2][1] = dv_dz; gradU[2][2] = dw_dz;
            
            Compute_dscalar_center
            (
                mesh,
                i, j, k,
                mx, my, mz,
                lt, nvert, meshTag,
                &dtdc, &dtde, &dtdz
            );

            // compute temperature derivative w.r.t. the cartesian coordinates
            Compute_dscalar_dxyz
            (
                mesh,
                csi0, csi1, csi2,
                eta0, eta1, eta2,
                zet0, zet1, zet2,
                ajc,
                dtdc, dtde, dtdz,
                &dt_dx, &dt_dy, &dt_dz
            );

            gradT[0] = dt_dx;
            gradT[1] = dt_dy;
            gradT[2] = dt_dz;

            for (PetscInt a = 0; a < 3; a++)
            {
                for (PetscInt c = 0; c < 3; c++)
                {
                    num += (-gradU[c][a] * gradT[c] * gradT[a]);
                }
            }

            denom = gradT[0]*gradT[0] + gradT[1]*gradT[1] + gradT[2]*gradT[2];

            PetscReal diff_e = num / (denom + 1e-10);

            diff_e = csk[k][j][i] * delta2_harm * PetscMax(diff_e, 0.0);

            lkt[k][j][i] = diff_e;
        }
        else if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD)
        {
            PetscReal filter = pow( 1./aj[k][j][i], 1./3.);
            lkt[k][j][i] = csk[k][j][i] * pow(filter, 2.0) * Sabs;
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);
    DMDAVecRestoreArray(da,  les->lCsk, &csk);
    DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
    DMDAVecRestoreArray(da,  les->lNu_t, &nut);
    DMDAVecRestoreArray(da,  les->lk_t, &lkt);

    DMLocalToLocalBegin(da,  les->lk_t, INSERT_VALUES, les->lk_t);
    DMLocalToLocalEnd  (da,  les->lk_t, INSERT_VALUES, les->lk_t);    

    UpdatektBCs(les);
    return(0);
}

PetscErrorCode updateLESScalarStructuralModel(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    ueqn_         *ueqn = les->access->ueqn;
    teqn_         *teqn = les->access->teqn;

    DM            da    = mesh->da, fda = mesh->fda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***nvert, ***meshTag;
    Cmpnts        ***csi, ***eta, ***zet, ***ucat;
    Cmpnts        ***visc;
    PetscReal     ***lt;
    
    Cmpnts        ***qSGS, qSGSLocal;                                      //subgrid scale heat flux 
    
    lxs = xs; lxe = xe; if (xs == 0) lxs = xs+1; if (xe == mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys == 0) lys = ys+1; if (ye == my) lye = ye-1;
    lzs = zs; lze = ze; if (zs == 0) lzs = zs+1; if (ze == mz) lze = mz-1;

    VecSet(les->lQ, 0.0);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(da, teqn->lTmprt, &lt);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);

    DMDAVecGetArray(fda, les->lQ, &qSGS);

    //bardina model
    if (les->model == BAMD || les->model == BV)
    {

        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {
            /* Skip solid cells */
            if( isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                qSGS[k][j][i] = nSetZero();
                continue;
            }

            PetscReal u[3][3][3] ,   v[3][3][3] ,   w[3][3][3];    // cartesian velocity
            
            PetscReal U[3][3][3]   , V[3][3][3]   , W[3][3][3];

            PetscReal Utheta[3][3][3]  , Vtheta[3][3][3]  , Wtheta[3][3][3]  ;
            
            PetscReal theta[3][3][3]  ;
                        
            PetscInt p, q, r;

            for (p = -1; p <= 1; p++)
            for (q = -1; q <= 1; q++)
            for (r = -1; r <= 1; r++)
            {
                PetscInt R = r + 1, Q = q + 1, P = p + 1;
                PetscInt K = k + r, J = j + q, I = i + p;

                //cartesian velocity
                u[R][Q][P] = ucat[K][J][I].x;
                v[R][Q][P] = ucat[K][J][I].y;
                w[R][Q][P] = ucat[K][J][I].z;
            
                // contravariant fluxes
                U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
                V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
                W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;

                theta[R][Q][P] = lt[K][J][I];

                Utheta[R][Q][P] = U[R][Q][P] * theta[R][Q][P];
                Vtheta[R][Q][P] = V[R][Q][P] * theta[R][Q][P];
                Wtheta[R][Q][P] = W[R][Q][P] * theta[R][Q][P];
            }

            PetscReal weight[3][3][3];

            for (p = -1; p <= 1; p++)
            for (q = -1; q <= 1; q++)
            for (r = -1; r <= 1; r++)
            {
                PetscInt R = r + 1, Q = q + 1, P = p + 1;
      
                PetscInt K = k + r, J = j + q, I = i + p;
                
                if 
                (
                    (isFluidCell(K, J, I, nvert) || isIBMFluidCell(K, J, I, nvert)) &&
                    
                    (I != 0 && I != mx-1 && J != 0 && J != my-1 && K != 0 && K != mz-1))

                {
                    weight[R][Q][P] = 1;                
                } 
                
                else 
                {
                    weight[R][Q][P] = 0;
                }
            
            }
            
            //Apply test filter to local velocity fields
            
            PetscReal _U = integrateTestfilterSimpson(U, weight);
            PetscReal _V = integrateTestfilterSimpson(V, weight);
            PetscReal _W = integrateTestfilterSimpson(W, weight);

            PetscReal _theta = integrateTestfilterSimpson(theta, weight);

            //Apply test filter to product
            PetscReal _Utheta = integrateTestfilterSimpson(Utheta, weight);        
            PetscReal _Vtheta = integrateTestfilterSimpson(Vtheta, weight);
            PetscReal _Wtheta = integrateTestfilterSimpson(Wtheta, weight);

            //Compute the Bardina SGS stress ---
            //q_i = test-filtered(U_i*theta) - (test-filtered(U_i))*(test-filtered(theta))
            
            qSGS[k][j][i].x = _Utheta - _U*_theta;   
            qSGS[k][j][i].y = _Vtheta - _V*_theta;   
            qSGS[k][j][i].z = _Wtheta - _W*_theta;   
        } 
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    
    DMDAVecRestoreArray(fda, les->lQ, &qSGS);

    DMLocalToLocalBegin(fda, les->lQ, INSERT_VALUES, les->lQ);
    DMLocalToLocalEnd  (fda, les->lQ, INSERT_VALUES, les->lQ);

    //set tau boundary conditions
    DMDAVecGetArray(fda, les->lQ, &qSGS);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                PetscInt a=i, b=j, c=k, flag=0;

                if(i==0)
                {
                    if(mesh->i_periodic)       a=mx-2, flag=1;
                    else if(mesh->ii_periodic) a=-2, flag=1;
                    else                      a=1, flag=1;
                }
                if(i==mx-1)
                {
                    if(mesh->i_periodic)       a=1, flag=1;
                    else if(mesh->ii_periodic) a=mx+1, flag=1;
                    else                      a=mx-2, flag=1;
                }
                if(j==0)
                {
                    if(mesh->j_periodic)       b=my-2, flag=1;
                    else if(mesh->jj_periodic) b=-2, flag=1;
                    else                      b=1, flag=1;
                }
                if(j==my-1)
                {
                    if(mesh->j_periodic)       b=1, flag=1;
                    else if(mesh->jj_periodic) b=my+1, flag=1;
                    else                      b=my-2, flag=1;
                }
                if(k==0)
                {
                    if(mesh->k_periodic)       c=mz-2, flag=1;
                    else if(mesh->kk_periodic) c=-2, flag=1;
                    else                      c=1, flag=1;
                }
                if(k==mz-1)
                {
                    if(mesh->k_periodic)       c=1, flag=1;
                    else if(mesh->kk_periodic) c=mz+1, flag=1;
                    else                      c=mz-2, flag=1;
                }

                if(flag)
                {
                    qSGS[k][j][i].x = qSGS[c][b][a].x;
                    qSGS[k][j][i].y = qSGS[c][b][a].y;
                    qSGS[k][j][i].z = qSGS[c][b][a].z;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, les->lQ, &qSGS);
    DMLocalToLocalBegin(fda, les->lQ, INSERT_VALUES, les->lQ);
    DMLocalToLocalEnd  (fda, les->lQ, INSERT_VALUES, les->lQ);

    DMDAVecGetArray(fda, teqn->lViscT, &visc);
    DMDAVecGetArray(fda, les->lQ, &qSGS);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);

    //i direction
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(j==0 || k==0) continue;

                qSGSLocal.x = 0.5 * (qSGS[k][j][i].x + qSGS[k][j][i+1].x);

                if
                (
                    (mesh->boundaryT.iLeft=="thetaWallFunction" && i==0) ||
                    (mesh->boundaryT.iRight=="thetaWallFunction" && i==mx-2)
                )
                {

                }
                else if(isIBMFluidIFace(k, j, i, i+1, nvert))
                {

                }
                else 
                {
                    visc[k][j][i].x -= qSGSLocal.x;
                }
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

                qSGSLocal.y = 0.5 * (qSGS[k][j][i].y + qSGS[k][j+1][i].y);

                if
                (
                    (mesh->boundaryT.jLeft=="thetaWallFunction"  && j==0) ||
                    (mesh->boundaryT.jRight=="thetaWallFunction" && j==my-2)
                )
                {

                }
                else if(isIBMFluidJFace(k, j, i, j+1, nvert))
                {

                }
                else
                {
                    // if(k ==15 && j==2 && i==2)	
                    // {
                    //     PetscPrintf(PETSC_COMM_SELF,"before visc %.10lf qSGSLocal %.10lf\n", visc[k][j][i].y, qSGSLocal.y);
                    // }
                    visc[k][j][i].y -= qSGSLocal.y;

                    // if(k ==15 && j==2 && i==2)
                    // PetscPrintf(PETSC_COMM_SELF,"after visc %.10lf \n", visc[k][j][i].y);
                } 
            }
        }
    }

    // k direction faces are from 0 to mz-2
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if(i==mx-1 || j==my-1 || k==mz-1) continue;
                if(i==0 || j==0) continue;

                qSGSLocal.z = 0.5 * (qSGS[k][j][i].z + qSGS[k+1][j][i].z);

                if(isIBMFluidKFace(k, j, i, k+1, nvert))
                {
                    
                }
                else
                {
                    visc[k][j][i].z -= qSGSLocal.z;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, teqn->lViscT, &visc);
    DMDAVecRestoreArray(fda, les->lQ, &qSGS);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    DMLocalToLocalBegin(fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);
    DMLocalToLocalEnd  (fda, teqn->lViscT, INSERT_VALUES, teqn->lViscT);

    return 0;
}
