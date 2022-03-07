//! \file  les.c
//! \brief Contains LES model function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

const PetscReal wall_cs = 0.001 ;
const PetscReal std_cs  = 0.0289; // standard Cs value for HIT

//***************************************************************************************************************//

PetscErrorCode InitializeLES(les_ *les)
{
    if(les != NULL)
    {
        // set pointer to mesh
        mesh_ *mesh = les->access->mesh;

        les->maxCs  = 0.5;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-max_cs",  &(les->maxCs), PETSC_NULL);

        if(les->access->flags->isLesActive)
        {
            VecDuplicate(mesh->lAj, &(les->lCs));     VecSet(les->lCs,0.);
            VecDuplicate(mesh->lAj, &(les->lNu_t));   VecSet(les->lNu_t,0.);

            if (les->access->flags->isLesActive > 1)
            {
                VecDuplicate(mesh->lAj, &(les->lLM));  VecSet(les->lLM,0.);
                VecDuplicate(mesh->lAj, &(les->lMM));  VecSet(les->lMM,0.);

                VecDuplicate(mesh->lCsi, &les->lSx);   VecSet(les->lSx,0.);
                VecDuplicate(mesh->lCsi, &les->lSy);   VecSet(les->lSy,0.);
                VecDuplicate(mesh->lCsi, &les->lSz);   VecSet(les->lSz,0.);
                VecDuplicate(mesh->lAj, &(les->lS));   VecSet(les->lS,0.);

                if (les->access->flags->isLesActive == 4)
                {
                    VecDuplicate(mesh->lAj, &(les->lLM_old));   VecSet(les->lLM_old,0.);
                    VecDuplicate(mesh->lAj, &(les->lMM_old));   VecSet(les->lMM_old,0.);
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateCs (les_ *les)
{
    mesh_          *mesh  = les->access->mesh;
    ueqn_          *eqn   = les->access->ueqn;
    clock_         *clock = les->access->clock;
    DM             da     = mesh->da, fda = mesh->fda;

    // standard smagorinsky model
    if(les->access->flags->isLesActive==1)
    {
        VecSet(les->lCs, std_cs);
        DMLocalToLocalBegin(da, les->lCs, INSERT_VALUES, les->lCs);
        DMLocalToLocalEnd  (da, les->lCs, INSERT_VALUES, les->lCs);
        return(0);
    }

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    Cmpnts        ***ucont, ***ucat;
    Cmpnts        ***csi, ***eta, ***zet,
                  ***icsi, ***ieta, ***izet,
                  ***jcsi, ***jeta, ***jzet,
                  ***kcsi, ***keta, ***kzet;

    PetscReal     ***iaj, ***jaj, ***kaj, ***aj;

    PetscReal     ***nvert, ***Cs;

    Cmpnts        ***Ax, ***Ay, ***Az, ***cent;
    PetscReal     ***LM, ***MM;
    PetscReal     ***Sabs;

    PetscReal     ajc;

    PetscReal     dudc, dude, dudz,
                  dvdc, dvde, dvdz,
                  dwdc, dwde, dwdz;

    const char*   lagrangianAverageCs = "trilinear";    // nearestPoint
                                                        // trilinear

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecSet(les->lSx, 0.0);
    VecSet(les->lSy, 0.0);
    VecSet(les->lSz, 0.0);
    VecSet(les->lS,  0.0);
    VecSet(les->lLM, 0.0);
    VecSet(les->lMM, 0.0);

    // get model parameters
    DMDAVecGetArray(fda, les->lSx, &Ax);
    DMDAVecGetArray(fda, les->lSy, &Ay);
    DMDAVecGetArray(fda, les->lSz, &Az);
    DMDAVecGetArray(da,  les->lS, &Sabs);
    DMDAVecGetArray(da,  les->lLM, &LM);
    DMDAVecGetArray(da,  les->lMM, &MM);
    DMDAVecGetArray(da,  les->lCs, &Cs);

    // get fields
    DMDAVecGetArray(fda, eqn->lUcont, &ucont);
    DMDAVecGetArray(fda, eqn->lUcat,  &ucat);

    // get jacobians
    DMDAVecGetArray(da, mesh->lAj, &aj);
    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    // get face area vectors
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

    // get IBM markup and cell centers
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // 1 - loop over internal cells and compute strain rate tensor S and |S|
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // if on body skip
                if( isIBMSolidCell(k, j, i, nvert))
                {
                    continue;
                }

                // get cell metrics
                ajc = aj[k][j][i];
                PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
                PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
                PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

                PetscReal dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
                PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

                // compute velocity derivatices
                Compute_du_center
                (
                    mesh,
                    i, j, k,
                    mx, my, mz,
                    ucat, nvert,
                    &dudc, &dvdc, &dwdc,
                    &dude, &dvde, &dwde,
                    &dudz, &dvdz, &dwdz
                );

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
                    &du_dx, &dv_dx,
                    &dw_dx, &du_dy,
                    &dv_dy, &dw_dy,
                    &du_dz, &dv_dz,
                    &dw_dz
                );

                PetscReal Sxx = 0.5*(du_dx + du_dx),
                          Sxy = 0.5*(du_dy + dv_dx),
                          Sxz = 0.5*(du_dz + dw_dx);

                PetscReal Syx = Sxy,
                          Syy = 0.5*(dv_dy + dv_dy),
                          Syz = 0.5*(dv_dz + dw_dy);
                PetscReal Szx = Sxz,
                          Szy = Syz,
                          Szz = 0.5*(dw_dz + dw_dz);

                // |S| = sqrt(2*SijSij)
                Sabs[k][j][i]
                =
                sqrt
                (
                    2.0*
                    (
                        Sxx*Sxx + Sxy*Sxy + Sxz*Sxz +
                        Syx*Syx + Syy*Syy + Syz*Syz +
                        Szx*Szx + Szy*Szy + Szz*Szz
                    )
                );

                // Sij
                Ax[k][j][i].x = du_dx; Ax[k][j][i].y = du_dy; Ax[k][j][i].z = du_dz;
                Ay[k][j][i].x = dv_dx; Ay[k][j][i].y = dv_dy; Ay[k][j][i].z = dv_dz;
                Az[k][j][i].x = dw_dx; Az[k][j][i].y = dw_dy; Az[k][j][i].z = dw_dz;
            }
        }
    }

    // restore arrays
    DMDAVecRestoreArray(fda, les->lSx, &Ax);
    DMDAVecRestoreArray(fda, les->lSy, &Ay);
    DMDAVecRestoreArray(fda, les->lSz, &Az);
    DMDAVecRestoreArray(da,  les->lS,  &Sabs);

    // update ghost values
    DMLocalToLocalBegin(fda, les->lSx, INSERT_VALUES, les->lSx);
    DMLocalToLocalEnd  (fda, les->lSx, INSERT_VALUES, les->lSx);
    DMLocalToLocalBegin(fda, les->lSy, INSERT_VALUES, les->lSy);
    DMLocalToLocalEnd  (fda, les->lSy, INSERT_VALUES, les->lSy);
    DMLocalToLocalBegin(fda, les->lSz, INSERT_VALUES, les->lSz);
    DMLocalToLocalEnd  (fda, les->lSz, INSERT_VALUES, les->lSz);
    DMLocalToLocalBegin(da,  les->lS,  INSERT_VALUES, les->lS);
    DMLocalToLocalEnd  (da,  les->lS,  INSERT_VALUES, les->lS);

    // handle periodicity
    resetCellPeriodicFluxes(mesh, les->lSx, les->lSx, "vector", "localToLocal");
    resetCellPeriodicFluxes(mesh, les->lSy, les->lSy, "vector", "localToLocal");
    resetCellPeriodicFluxes(mesh, les->lSz, les->lSz, "vector", "localToLocal");
    resetCellPeriodicFluxes(mesh, les->lS,  les->lS,  "scalar", "localToLocal");

    // get arrays
    DMDAVecGetArray(fda, les->lSx, &Ax);
    DMDAVecGetArray(fda, les->lSy, &Ay);
    DMDAVecGetArray(fda, les->lSz, &Az);
    DMDAVecGetArray(da,  les->lS,  &Sabs);

    // 2 - loop over cells and compute LijMij (LM) and MijMij (MM)
    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        // set to zero if solid and skip
        if(isIBMSolidCell(k, j, i, nvert))
        {
            LM[k][j][i] = MM[k][j][i] = 0;
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
        PetscReal Lij[3][3],                                   // Germano tensor projected on generalized coordinates
                  Sij_hat[3][3],                               // filtered strain rate tensor S projected on generalized coordinates
                  S_hat,                                       // scalar |S| = sqrt(2*Sij*Sij)
                  SSij_hat[3][3],                              // filtered |S|*S tensor
                  Mij[3][3],                                   // Mij tensor
                  Nij[3][3],                                   // Nij tensor
                  Lij_cat[3][3],                               // cartesian Germano tensor
                  Mij_cat[3][3];                               // cartesian Mij tensor

        PetscReal filter, test_filter;                         // grid filter and dynamic model filter

        // local tensors before the test-filtering (3x3x3 cell stencil)
        PetscReal S[3][3][3];                                  // scalar |S| = sqrt(2*Sij*Sij)
        PetscReal u[3][3][3] ,   v[3][3][3] ,   w[3][3][3] ,   // cartesian velocity
                  d[3][3][3];
        PetscReal U[3][3][3]   , V[3][3][3]   , W[3][3][3]   ; // contravariant fluxes
        PetscReal Uu[3][3][3]  , Uv[3][3][3]  , Uw[3][3][3]  ; //         |Uu Uv Uw|
        PetscReal Vu[3][3][3]  , Vv[3][3][3]  , Vw[3][3][3]  ; // tensor: |Vu Vv Vw|
        PetscReal Wu[3][3][3]  , Wv[3][3][3]  , Ww[3][3][3]  ; //         |Wu Wv Ww|
        PetscReal S11[3][3][3] , S12[3][3][3] , S13[3][3][3] , // strain rate tensor
                  S21[3][3][3] , S22[3][3][3] , S23[3][3][3] ,
                  S31[3][3][3] , S32[3][3][3] , S33[3][3][3] ;
        PetscReal SS11[3][3][3], SS12[3][3][3], SS13[3][3][3], // |S|*S tensor
                  SS21[3][3][3], SS22[3][3][3], SS23[3][3][3],
                  SS31[3][3][3], SS32[3][3][3], SS33[3][3][3];

        // creating variables to be test-filtered
        PetscInt p,q,r;
        for(p=-1; p<=1; p++)
        for(q=-1; q<=1; q++)
        for(r=-1; r<=1; r++)
        {
            PetscInt R = r+1, Q = q+1, P = p+1;                  // local stencil labeling: putting data
            PetscInt K = k+r, J = j+q, I = i+p;                  // mesh cell labeling    : getting data

            // cartesian velocity
            u[R][Q][P] = ucat[K][J][I].x;
            v[R][Q][P] = ucat[K][J][I].y;
            w[R][Q][P] = ucat[K][J][I].z;

            // contravariant fluxes
            U[R][Q][P] = u[R][Q][P]*csi[K][J][I].x + v[R][Q][P]*csi[K][J][I].y + w[R][Q][P]*csi[K][J][I].z;
            V[R][Q][P] = u[R][Q][P]*eta[K][J][I].x + v[R][Q][P]*eta[K][J][I].y + w[R][Q][P]*eta[K][J][I].z;
            W[R][Q][P] = u[R][Q][P]*zet[K][J][I].x + v[R][Q][P]*zet[K][J][I].y + w[R][Q][P]*zet[K][J][I].z;

            // Uu tensor components
            Uu[R][Q][P] = U[R][Q][P] * u[R][Q][P];
            Uv[R][Q][P] = U[R][Q][P] * v[R][Q][P];
            Uw[R][Q][P] = U[R][Q][P] * w[R][Q][P];

            Vu[R][Q][P] = V[R][Q][P] * u[R][Q][P];
            Vv[R][Q][P] = V[R][Q][P] * v[R][Q][P];
            Vw[R][Q][P] = V[R][Q][P] * w[R][Q][P];

            Wu[R][Q][P] = W[R][Q][P] * u[R][Q][P];
            Wv[R][Q][P] = W[R][Q][P] * v[R][Q][P];
            Ww[R][Q][P] = W[R][Q][P] * w[R][Q][P];

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


            // Mohammad: the model is expressed in cartesian coordinates
            S11[R][Q][P]  = Sxx, S12[R][Q][P] = Sxy, S13[R][Q][P] = Sxz;
            S21[R][Q][P]  = Syx, S22[R][Q][P] = Syy, S23[R][Q][P] = Syz;
            S31[R][Q][P]  = Szx, S32[R][Q][P] = Szy, S33[R][Q][P] = Szz;

            /*
            // seba: the strain rate tensor is transformed in curvilinear coordinates
            //       Need more work since the value of nut is too small. Probably need
            //       to use also the transformed Sabs instead of the cartesian one.
            S11[R][Q][P]  = Sxx*csi[K][J][I].x + Sxy*csi[K][J][I].y + Sxz*csi[K][J][I].z;
            S12[R][Q][P]  = Sxx*eta[K][J][I].x + Sxy*eta[K][J][I].y + Sxz*eta[K][J][I].z;
            S13[R][Q][P]  = Sxx*zet[K][J][I].x + Sxy*zet[K][J][I].y + Sxz*zet[K][J][I].z;
            S21[R][Q][P]  = Syx*csi[K][J][I].x + Syy*csi[K][J][I].y + Syz*csi[K][J][I].z;
            S22[R][Q][P]  = Syx*eta[K][J][I].x + Syy*eta[K][J][I].y + Syz*eta[K][J][I].z;
            S23[R][Q][P]  = Syx*zet[K][J][I].x + Syy*zet[K][J][I].y + Syz*zet[K][J][I].z;
            S31[R][Q][P]  = Szx*csi[K][J][I].x + Szy*csi[K][J][I].y + Szz*csi[K][J][I].z;
            S32[R][Q][P]  = Szx*eta[K][J][I].x + Szy*eta[K][J][I].y + Szz*eta[K][J][I].z;
            S33[R][Q][P]  = Szx*zet[K][J][I].x + Szy*zet[K][J][I].y + Szz*zet[K][J][I].z;
            */

            // |S|
            S[R][Q][P]    = Sabs[K][J][I];

            // |S|*S
            SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
            SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
            SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];
        }

        // 2.2 - Build the filters

        // filter kernel
        PetscReal coef[3][3][3]=
        {
            0.125, 0.250, 0.125,
            0.250, 0.500, 0.250,
            0.125, 0.250, 0.125,

            0.250, 0.500, 0.250,
            0.500, 1.000, 0.500,
            0.250, 0.500, 0.250,

            0.125, 0.250, 0.125,
            0.250, 0.500, 0.250,
            0.125, 0.250, 0.125
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

            // fluid
            if
            (
                (isFluidCell(K, J, I, nvert)) &&
                //!isOnCornerCellCenters(I, J, K, mesh->info) <- use also ghost? Not for now
                (I!=0 && I!=mx-1 && J!=0 && J!=my-1 && K!=0 && K!=mz-1)
            )
            {
                sum_vol += 1./aj[K][J][I] * coef[R][Q][P];
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
        filter = pow( 1./aj[k][j][i], 1./3.);

        // use the cubic root of total weighted volume
        test_filter = pow( sum_vol, 1./3. );

        // apply test-filtering
        PetscReal _U = integrateTestfilterSimpson(U, weight);
        PetscReal _V = integrateTestfilterSimpson(V, weight);
        PetscReal _W = integrateTestfilterSimpson(W, weight);

        PetscReal _u = integrateTestfilterSimpson(u, weight);
        PetscReal _v = integrateTestfilterSimpson(v, weight);
        PetscReal _w = integrateTestfilterSimpson(w, weight);
        PetscReal _d = 1;

        Lij[0][0] = integrateTestfilterSimpson(Uu, weight) - _U*_u;
        Lij[0][1] = integrateTestfilterSimpson(Uv, weight) - _U*_v;
        Lij[0][2] = integrateTestfilterSimpson(Uw, weight) - _U*_w;
        Lij[1][0] = integrateTestfilterSimpson(Vu, weight) - _V*_u;
        Lij[1][1] = integrateTestfilterSimpson(Vv, weight) - _V*_v;
        Lij[1][2] = integrateTestfilterSimpson(Vw, weight) - _V*_w;
        Lij[2][0] = integrateTestfilterSimpson(Wu, weight) - _W*_u;
        Lij[2][1] = integrateTestfilterSimpson(Wv, weight) - _W*_v;
        Lij[2][2] = integrateTestfilterSimpson(Ww, weight) - _W*_w;

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

        SSij_hat[0][0] = integrateTestfilterSimpson(SS11, weight);
        SSij_hat[0][1] = integrateTestfilterSimpson(SS12, weight);
        SSij_hat[0][2] = integrateTestfilterSimpson(SS13, weight);
        SSij_hat[1][0] = integrateTestfilterSimpson(SS21, weight);
        SSij_hat[1][1] = integrateTestfilterSimpson(SS22, weight);
        SSij_hat[1][2] = integrateTestfilterSimpson(SS23, weight);
        SSij_hat[2][0] = integrateTestfilterSimpson(SS31, weight);
        SSij_hat[2][1] = integrateTestfilterSimpson(SS32, weight);
        SSij_hat[2][2] = integrateTestfilterSimpson(SS33, weight);


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

        // compute cartesian Mij tensor
        for(a=0; a<3; a++)
        for(b=0; b<3; b++)
        {
            // Mohammad: Mij first defined in cartesian and then transformed
            Mij_cat[a][b]

            // Seba: Mij readily defined in curvilinear coordinates
            //Mij[a][b]
            =
            (
                - pow( test_filter, 2. ) * S_hat * Sij_hat[a][b]
                + pow( filter, 2. ) * SSij_hat[a][b]
            );

        }
        // Mohammad: project the Mij
        Mij[0][0] = Mij_cat[0][0] * csi0 + Mij_cat[0][1] * csi1 + Mij_cat[0][2] * csi2;
        Mij[0][1] = Mij_cat[0][0] * eta0 + Mij_cat[0][1] * eta1 + Mij_cat[0][2] * eta2;
        Mij[0][2] = Mij_cat[0][0] * zet0 + Mij_cat[0][1] * zet1 + Mij_cat[0][2] * zet2;
        Mij[1][0] = Mij_cat[1][0] * csi0 + Mij_cat[1][1] * csi1 + Mij_cat[1][2] * csi2;
        Mij[1][1] = Mij_cat[1][0] * eta0 + Mij_cat[1][1] * eta1 + Mij_cat[1][2] * eta2;
        Mij[1][2] = Mij_cat[1][0] * zet0 + Mij_cat[1][1] * zet1 + Mij_cat[1][2] * zet2;
        Mij[2][0] = Mij_cat[2][0] * csi0 + Mij_cat[2][1] * csi1 + Mij_cat[2][2] * csi2;
        Mij[2][1] = Mij_cat[2][0] * eta0 + Mij_cat[2][1] * eta1 + Mij_cat[2][2] * eta2;
        Mij[2][2] = Mij_cat[2][0] * zet0 + Mij_cat[2][1] * zet1 + Mij_cat[2][2] * zet2;

        PetscReal num = 0, num1 = 0, denom = 0;
        PetscInt  m, n, l;

        // compute numerator
        for(q=0; q<3; q++)
        for(a=0; a<3; a++)
        for(b=0; b<3; b++)
        {
            // this contraction is given by
            // L G M = Lba * (Maq * Gbq)
            // Note: for how tensors are build the matrix producti is row by row and
            //       not row by colums as usual. In both products.
            num += Lij[b][a] * Mij[a][q] * G[b][q];

        }

        // compute denominator
        for(m=0; m<3; m++)
        for(n=0; n<3; n++)
        for(l=0; l<3; l++)
        {
            // this contraction is given by
            // M G M = Mnm * (Mnl * Gml)
            // Note: matrix product is row by row in the product between parenthesis
            //       and row by column between the first matrix and the result of previous one.
            denom += Mij[n][m] * Mij[n][l] * G[m][l];
        }

        // store in LM and MM
        LM[k][j][i] = num;
        MM[k][j][i] = denom;
    }

    DMDAVecRestoreArray(da, les->lLM, &LM);
    DMDAVecRestoreArray(da, les->lMM, &MM);

    // update ghost nodes
    DMLocalToLocalBegin(da, les->lLM, INSERT_VALUES, les->lLM);
    DMLocalToLocalEnd  (da, les->lLM, INSERT_VALUES, les->lLM);
    DMLocalToLocalBegin(da, les->lMM, INSERT_VALUES, les->lMM);
    DMLocalToLocalEnd  (da, les->lMM, INSERT_VALUES, les->lMM);

    PetscReal ***LM_old, ***MM_old;

    // local convective weighting
    if(les->access->flags->isLesActive==4)
    {
        // use standard Smagorinsky for first iteration (since need old values of LM and MM)
        if ( clock->it < clock->itStart + 2)
        {
            DMDAVecGetArray(da, les->lLM, &LM);
            DMDAVecGetArray(da, les->lMM, &MM);

            DMDAVecGetArray(da, les->lLM_old, &LM_old);
            DMDAVecGetArray(da, les->lMM_old, &MM_old);

            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        LM_old[k][j][i] = std_cs*MM[k][j][i];
                        MM_old[k][j][i] = MM[k][j][i];
                    }
                }
            }

            DMDAVecRestoreArray(da, les->lLM, &LM);
            DMDAVecRestoreArray(da, les->lMM, &MM);

            DMDAVecRestoreArray(da, les->lLM_old, &LM_old);
            DMDAVecRestoreArray(da, les->lMM_old, &MM_old);

            DMLocalToLocalBegin(da, les->lLM_old, INSERT_VALUES, les->lLM_old);
            DMLocalToLocalEnd  (da, les->lLM_old, INSERT_VALUES, les->lLM_old);

            DMLocalToLocalBegin(da, les->lMM_old, INSERT_VALUES, les->lMM_old);
            DMLocalToLocalEnd  (da, les->lMM_old, INSERT_VALUES, les->lMM_old);

            VecSet(les->lCs, std_cs);

            return(0);

        }
        else
        {

            DMDAVecGetArray(da, les->lLM, &LM);
            DMDAVecGetArray(da, les->lMM, &MM);

            DMDAVecGetArray(da, les->lLM_old, &LM_old);
            DMDAVecGetArray(da, les->lMM_old, &MM_old);

            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                // filter dimension
                PetscReal filter = pow( 1./aj[k][j][i], 1./3. );

                // The higher the time scale, the strongest the lagrangian scmoothing.
                // To low time scales could result in weighting too much the new
                // value, which in turn could be clipped at zero. Hence we have noisier fields.

                // Meneveau time scale
                PetscReal T_scale
                =
                PetscMax
                (
                    2.0*filter*
                    pow
                    (
                        fabs(LM_old[k][j][i]),
                        -1./4.
                    ),
                    1e-20
                );

                // Piomelli time scale
                /*
                PetscReal T_scale
                =
                PetscMax
                (
                    1.5*filter*
                    pow
                    (
                        8.*fabs(LM_old[k][j][i]*MM_old[k][j][i]),
                        -1./8.
                    ),
                    1e-15
                );
                */

                // Mohammad time scale
                /*
                PetscReal T_scale
                =
                1.5*filter*
                pow
                (
                    fabs(LM_old[k][j][i]*MM_old[k][j][i]) +
                    1.e-19,
                    -1.0/8.0
                ) +
                1.e-19;
                */

                PetscReal LM_new, MM_new;

                LM_new = LM[k][j][i];
                MM_new = MM[k][j][i];

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
                    ) && isFluidCell(k1, j1, i1, nvert))
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

                PetscReal _LM_old, _MM_old;

                // interpolate the value at the exat old point (slower)
                // if there is an IBM or overset cell around the given cell, use nearest cell
                if(isBoxIBMCell(k_old, j_old, i_old, nvert))
                {
                    _LM_old = LM_old[k_old][j_old][i_old];
                    _MM_old = MM_old[k_old][j_old][i_old];
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
                          LM_old,
                          _LM_old
                      );

                      scalarPointLocalVolumeInterpolation
                      (
                          mesh,
                          X_old.x, X_old.y, X_old.z,
                          i_old, j_old, k_old,
                          cent,
                          MM_old,
                          _MM_old
                      );
                  }
                  else if(lagrangianAverageCs == "nearestPoint")
                  {
                      // take the closest cell value (faster, but a lot of cs clipping from below)
                      _LM_old = LM_old[k_old][j_old][i_old];
                      _MM_old = MM_old[k_old][j_old][i_old];
                  }

                }

                // interpolation weights depend on the time scale.
                // if dt >> T_scale new value is used
                // if dt << T_scale old value is used
                PetscReal eps = (clock->dt/T_scale) / (1.0 + clock->dt/T_scale);

                // apply lagrangian weighting (relaxation)
                LM[k][j][i]
                =
                PetscMax
                (
                    (
                        (1-eps) * _LM_old +
                        eps * LM_new
                    ),
                    0.0
                );

                MM[k][j][i]
                =
                (1-eps) * _MM_old +
                eps * MM_new;
            }

            // store LM and MM old values
            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                LM_old[k][j][i] = LM[k][j][i];
                MM_old[k][j][i] = MM[k][j][i];
            }

            DMDAVecRestoreArray(da, les->lLM, &LM);
            DMDAVecRestoreArray(da, les->lMM, &MM);

            DMDAVecRestoreArray(da, les->lLM_old, &LM_old);
            DMDAVecRestoreArray(da, les->lMM_old, &MM_old);

            DMLocalToLocalBegin(da, les->lLM_old, INSERT_VALUES, les->lLM_old);
            DMLocalToLocalEnd  (da, les->lLM_old, INSERT_VALUES, les->lLM_old);

            DMLocalToLocalBegin(da, les->lMM_old, INSERT_VALUES, les->lMM_old);
            DMLocalToLocalEnd  (da, les->lMM_old, INSERT_VALUES, les->lMM_old);
        }
    }

    // update ghost values
    DMLocalToLocalBegin(da, les->lLM, INSERT_VALUES, les->lLM);
    DMLocalToLocalEnd  (da, les->lLM, INSERT_VALUES, les->lLM);
    DMLocalToLocalBegin(da, les->lMM, INSERT_VALUES, les->lMM);
    DMLocalToLocalEnd  (da, les->lMM, INSERT_VALUES, les->lMM);

    resetCellPeriodicFluxes(mesh, les->lLM, les->lLM, "scalar", "localToLocal");
    resetCellPeriodicFluxes(mesh, les->lMM, les->lMM, "scalar", "localToLocal");

    DMDAVecGetArray(da, les->lLM, &LM);
    DMDAVecGetArray(da, les->lMM, &MM);

    // 3 - Compute Cs

    // loop over internal cells
    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        // body
        if(isIBMSolidCell(k, j, i, nvert))
        {
            Cs[k][j][i] = 0;
            continue;
        }

        PetscReal weight[3][3][3];

        PetscReal LM0[3][3][3],
                  MM0[3][3][3];

        PetscInt  a, b, c;

        for(a=-1; a<=1; a++)
        for(b=-1; b<=1; b++)
        for(c=-1; c<=1; c++)
        {
            PetscInt R = c+1, Q = b+1, P = a+1;
            PetscInt K = k+c, J = j+b, I = i+a;

            // volume weighting the mean
            weight[R][Q][P] = 1./aj[K][J][I];

            if( isIBMSolidCell(K, J, I, nvert) )
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
                weight[R][Q][P]=0;
            }

            LM0[R][Q][P] = LM[K][J][I];
            MM0[R][Q][P] = MM[K][J][I];
        }

        PetscReal C = 0;

        PetscReal LM_avg, MM_avg;

        if
        (
            les->access->flags->isLesActive == 3 ||                                // damped
            les->access->flags->isLesActive == 4                                   // lagrangian relaxation
        )
        {
            // no average procedure
            LM_avg = LM[k][j][i];
            MM_avg = MM[k][j][i];
        }
        else
        {
            // integrate over the testfilter volume locally
            LM_avg = integrateTestfilterSimpson(LM0, weight);
            MM_avg = integrateTestfilterSimpson(MM0, weight);
        }

        // set Smagorinsky constant
        C = 0.5 * LM_avg / (MM_avg + 1.e-7);

        if ( les->access->flags->isLesActive==3 )
        {
            if(clock->it < 100)
            {

            }
            else
            {
                PetscReal rlx = 0.001;
                C = (1.0 - rlx) * Cs[k][j][i] + rlx * C;
            }
        }

        if(les->access->flags->isLesActive==1)
        {
            Cs[k][j][i] = std_cs;
        }
        else
        {
            Cs[k][j][i] = PetscMax(C, 0);
        }
    }

    for (k=zs; k<ze; k++)
    for (j=ys; j<ye; j++)
    for (i=xs; i<xe; i++)
    {
        // set Cs to zero if solid or on box corners
        if(isIBMSolidCell(k, j, i, nvert) || k==0 || k==mz-1 || j==0 || j==my-1 || i==0 || i==mx-1)
        {
            Cs[k][j][i] = 0;
        }
        else
        {
            // add viscosity near the wall with IBM to stabilize the solution
            if(isIBMFluidCell(k, j, i, nvert))
            {
                Cs[k][j][i] = PetscMax(wall_cs, Cs[k][j][i]);    // stabilize at high Re, osl 0.005
            }

            // clip Cs between upper and lower bound
            Cs[k][j][i] = PetscMin( PetscMax(Cs[k][j][i], 0.0), les->maxCs);
        }
    }

    // restore numerator and denominator of the model coefficient
    DMDAVecRestoreArray(fda, les->lSx, &Ax);
    DMDAVecRestoreArray(fda, les->lSy, &Ay);
    DMDAVecRestoreArray(fda, les->lSz, &Az);
    DMDAVecRestoreArray(da,  les->lS,  &Sabs);
    DMDAVecRestoreArray(da,  les->lLM, &LM);
    DMDAVecRestoreArray(da,  les->lMM, &MM);
    DMDAVecRestoreArray(da,  les->lCs, &Cs);

    // restore fields
    DMDAVecRestoreArray(fda, eqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, eqn->lUcat,  &ucat);

    // restore jacobians
    DMDAVecRestoreArray(da, mesh->lAj, &aj);
    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    // restore face area vectors
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

    // restore IB markup and cell centers
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    // scatter the model coefficient
    DMLocalToLocalBegin(da, les->lCs, INSERT_VALUES, les->lCs);
    DMLocalToLocalEnd  (da, les->lCs, INSERT_VALUES, les->lCs);

    // handle periodicity and scatter
    resetCellPeriodicFluxes(mesh, les->lCs, les->lCs, "scalar", "localToLocal");

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateNut(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    ueqn_         *eqn  = les->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***Cs, ***lnu_t, ***nvert, ***aj, ***ustar;
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

    VecSet(les->lNu_t, 0.);

    DMDAVecGetArray(fda, eqn->lUcat,  &ucat);
    DMDAVecGetArray(da,  eqn->lUstar, &ustar);

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lAj, &aj);

    DMDAVecGetArray(da,  les->lNu_t, &lnu_t);
    DMDAVecGetArray(da,  les->lCs, &Cs);


    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        if(isIBMSolidCell(k, j, i, nvert))
        {
            lnu_t[k][j][i]=0;
            continue;
        }

        PetscReal ajc  = aj[k][j][i];
        PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

        PetscReal dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

        // compute velocity derivatives w.r.t the curvilinear coordinates
        Compute_du_center
        (
            mesh,
            i, j, k,
            mx, my, mz,
            ucat, nvert,
            &dudc, &dvdc, &dwdc,
            &dude, &dvde, &dwde,
            &dudz, &dvdz, &dwdz
        );

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

        // compute rate of strain tensor
        PetscReal Sxx = 0.5*(du_dx + du_dx),
               Sxy = 0.5*(du_dy + dv_dx),
               Sxz = 0.5*(du_dz + dw_dx);
        PetscReal Syx = Sxy,
               Syy = 0.5*(dv_dy + dv_dy),
               Syz = 0.5*(dv_dz + dw_dy);
        PetscReal Szx = Sxz,
               Szy = Syz,
               Szz = 0.5*(dw_dz + dw_dz);

        // compute rate of strain tensor norm
        PetscReal Sabs
        =
        sqrt
        (
            2.0*
            (
                Sxx*Sxx +
                Sxy*Sxy +
                Sxz*Sxz +
                Syx*Syx +
                Syy*Syy +
                Syz*Syz +
                Szx*Szx +
                Szy*Szy +
                Szz*Szz
            )
        );

        PetscReal filter;

        filter = pow( 1./aj[k][j][i], 1./3.);

        lnu_t[k][j][i] = Cs[k][j][i] * pow ( filter, 2.0 ) * Sabs;

    }

    DMDAVecRestoreArray(fda, eqn->lUcat,  &ucat);
    DMDAVecRestoreArray(da,  eqn->lUstar, &ustar);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);

    DMDAVecRestoreArray(da,  les->lNu_t, &lnu_t);
    DMDAVecRestoreArray(da,  les->lCs, &Cs);

    DMLocalToLocalBegin(da,  les->lNu_t, INSERT_VALUES, les->lNu_t);
    DMLocalToLocalEnd  (da,  les->lNu_t, INSERT_VALUES, les->lNu_t);

    // correct boundary conditions
    UpdateNutBCs(les);

    return(0);
}
