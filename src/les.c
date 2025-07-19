//! \file  les.c
//! \brief Contains LES model function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

const PetscReal wall_cs = 0.001;
const PetscReal std_cs  = 0.0289;
const PetscReal amd_cs  =  0.1;
const PetscReal vreman_cs = 0.07225;
const PetscReal sa_cs     = 0.325; 
//***************************************************************************************************************//

static const char *LesModelNames[] = {
    "smagorinsky",
    "stabilityDependent",
    "dynamicSmagorinsky",
    "dynamicLASI",
    "dynamicLASD",
    "dynamicPASD",
    "amd",
    "vreman",
    "bardinaVreman",
    "bardinaAMD",
    "bardinaDSM",
    "scaleAdaptive",
    "m43",
    "LesModel", "", PETSC_NULL
};

PetscErrorCode InitializeLES(les_ *les)
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
    // 9. Bardina model 
    //10. Scale adaptive model

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
            PetscBool modelSet;
            PetscOptionsGetEnum(PETSC_NULL, PETSC_NULL, "-lesModel", LesModelNames, (PetscEnum*)&(les->model), &modelSet);

            // Check if a valid model was set
            if (!modelSet)
            {
                PetscPrintf(PETSC_COMM_WORLD, "Error: No valid LES model specified with -lesModel in control.dat.\n");
                PetscPrintf(PETSC_COMM_WORLD, "Available models:\n");
                for (int i = 0; LesModelNames[i] != PETSC_NULL && strcmp(LesModelNames[i], "LesModel") != 0; i++)
                {
                    PetscPrintf(PETSC_COMM_WORLD, "  %s\n", LesModelNames[i]);
                }
                char error[512];
                sprintf(error, "\nUse one of the available models.\n");
                fatalErrorInFunction("InitializeLES",  error);
            }    

            if(les->model == AMD)
            {
                les->amdCs = amd_cs;
                PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-amd_cs", &(les->amdCs), PETSC_NULL);
                PetscPrintf(PETSC_COMM_WORLD, "AMD Constant is set to: %.5lf \n",les->amdCs);
            }

            VecDuplicate(mesh->lAj, &(les->lCs));     VecSet(les->lCs, 0.);
            VecDuplicate(mesh->lAj, &(les->lNu_t));   VecSet(les->lNu_t, 0.);
            VecDuplicate(mesh->lAj, &(les->lKsgs));   VecSet(les->lKsgs, 0.);

            if (les->model == SMAGORINSKY ||
                les->model == STABILITY_BASED ||
                les->model == DSM    ||
                les->model == DLASI  ||
                les->model == DLASD  ||
                les->model == DPASD  ||
                les->model == AMD    ||
                les->model == VREMAN ||
                les->model == BV     ||
                les->model == BAMD   ||
                les->model == BDS)
            {
                VecDuplicate(mesh->lCsi, &(les->lSx)); VecSet(les->lSx, 0.);
                VecDuplicate(mesh->lCsi, &(les->lSy)); VecSet(les->lSy, 0.);
                VecDuplicate(mesh->lCsi, &(les->lSz)); VecSet(les->lSz, 0.);
                VecDuplicate(mesh->lAj, &(les->lS));   VecSet(les->lS, 0.);
            }

            // STABILITY_BASED specific
            if (les->model == STABILITY_BASED)
            {
                VecDuplicate(mesh->lAj, &(les->lN));   VecSet(les->lN, 0.);
                VecDuplicate(mesh->lAj, &(les->lCh));  VecSet(les->lCh, 0.);
                VecDuplicate(mesh->Nvert, &(les->L));  VecSet(les->L, 0.);
            }

            // Dynamic models (DSM, DLASI, DLASD, DPASD)
            if (les->model == DSM ||
                les->model == DLASI ||
                les->model == DLASD ||
                les->model == DPASD || 
                les->model == BDS)
            {
                VecDuplicate(mesh->lAj, &(les->lLM)); VecSet(les->lLM, 0.);
                VecDuplicate(mesh->lAj, &(les->lMM)); VecSet(les->lMM, 0.);
            }

            // Lagrangian averaging (DLASI, DLASD)
            if (les->model == DLASI || les->model == DLASD)
            {
                VecDuplicate(mesh->lAj, &(les->lLM_old)); VecSet(les->lLM_old, 0.);
                VecDuplicate(mesh->lAj, &(les->lMM_old)); VecSet(les->lMM_old, 0.);
            }

            // Scale-dependent dynamic models (DLASD, DPASD)
            if (les->model == DLASD || les->model == DPASD)
            {
                VecDuplicate(mesh->lAj, &(les->lQN)); VecSet(les->lQN, 0.);
                VecDuplicate(mesh->lAj, &(les->lNN)); VecSet(les->lNN, 0.);
            }

            // DLASD-specific
            if (les->model == DLASD)
            {
                VecDuplicate(mesh->lAj, &(les->lQN_old)); VecSet(les->lQN_old, 0.);
                VecDuplicate(mesh->lAj, &(les->lNN_old)); VecSet(les->lNN_old, 0.);
            }

            if (les->model == BAMD || les->model == BV || les->model == BDS)
            {
                DMCreateLocalVector(mesh->tda, &(les->lTau)); VecSet(les->lTau,0.);
            }
        }
        
        if (les->access->flags->isLesActive)
        {
            PetscPrintf(PETSC_COMM_WORLD, "\nSelected LES Model: %s\n", LesModelNames[les->model]);
        }
        else
        {
            PetscPrintf(PETSC_COMM_WORLD, "Selected LES Model: none (LES inactive)\n");
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateCs (les_ *les)
{
    mesh_          *mesh  = les->access->mesh;
    ueqn_          *ueqn   = les->access->ueqn;
    clock_         *clock = les->access->clock;
    DM             da     = mesh->da, fda = mesh->fda;

    // standard smagorinsky model
    if(les->model == SMAGORINSKY)
    {
        VecSet(les->lCs, std_cs);
        DMLocalToLocalBegin(da, les->lCs, INSERT_VALUES, les->lCs);
        DMLocalToLocalEnd  (da, les->lCs, INSERT_VALUES, les->lCs);
        return(0);
    }

    if(les->model == AMD || les->model == BAMD)
    {
        VecSet(les->lCs, les->amdCs);
        DMLocalToLocalBegin(da, les->lCs, INSERT_VALUES, les->lCs);
        DMLocalToLocalEnd  (da, les->lCs, INSERT_VALUES, les->lCs);
        return(0);
    }

    if(les->model == VREMAN || les->model == BV)
    {
        VecSet(les->lCs, vreman_cs);
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

    PetscReal     ***nvert, ***Cs, ***meshTag;

    Cmpnts        ***Ax, ***Ay, ***Az, ***cent;
    PetscReal     ***LM, ***MM, ***QN, ***NN;
    PetscReal     ***Sabs, ***n, ***lch, ***scale;

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

    // get model parameters
    DMDAVecGetArray(fda, les->lSx, &Ax);
    DMDAVecGetArray(fda, les->lSy, &Ay);
    DMDAVecGetArray(fda, les->lSz, &Az);
    DMDAVecGetArray(da,  les->lS,  &Sabs);
    DMDAVecGetArray(da,  les->lCs, &Cs);

    // get fields
    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    // get jacobians
    DMDAVecGetArray(da, mesh->lAj,  &aj);
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
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // 1 - loop over internal cells and compute strain rate tensor S and |S|
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // if on body skip
                if( isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
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
                    ucat, nvert, meshTag, 
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

    // 2 - loop over cells and compute specific model parameters (LijMij and MijMij for dynamic model and N for stability dep. model)
    if(les->model == STABILITY_BASED)
    {
        DMDAVecGetArray(da,  les->lN, &n);

        PetscReal g     = 9.81;
        PetscReal tRef  = les->access->abl->tRef;
        teqn_     *teqn = les->access->teqn;
        PetscReal ***tmprt;

        DMDAVecGetArray(da,  teqn->lTmprt,      &tmprt);

        for (k=lzs; k<lze; k++)
        for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++)
        {
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
                tmprt, nvert, meshTag,
                &dtdc, &dtde, &dtdz
            );

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

            n[k][j][i] = std::max(sqrt(std::max(g / tRef * dt_dz, 0.0)), 1e-10);
        }

        DMDAVecRestoreArray(da,  teqn->lTmprt,      &tmprt);
        DMDAVecRestoreArray(da,  les->lN, &n);

        DMLocalToLocalBegin(da, les->lN, INSERT_VALUES, les->lN);
        DMLocalToLocalEnd  (da, les->lN, INSERT_VALUES, les->lN);
    }
    else if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD || les->model == BDS)
    {
        // get arrays
        DMDAVecGetArray(fda, les->lSx, &Ax);
        DMDAVecGetArray(fda, les->lSy, &Ay);
        DMDAVecGetArray(fda, les->lSz, &Az);
        DMDAVecGetArray(da,  les->lS,  &Sabs);

        VecSet(les->lLM, 0.0);
        VecSet(les->lMM, 0.0);
        DMDAVecGetArray(da,  les->lLM, &LM);
        DMDAVecGetArray(da,  les->lMM, &MM);

        //first test filter - compute Lij and Mij
        for (k=lzs; k<lze; k++)
        for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++)
        {
            // set to zero if solid and skip
            if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
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

                // |S|
                S[R][Q][P]    = Sabs[K][J][I];

                // |S|*S
                SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
                SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
                SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];

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

        //second test filter - compute Qij and Nij
        if(les->model == DLASD || les->model == DPASD)
        {
            VecSet(les->lQN, 0.0);
            VecSet(les->lNN, 0.0);
            DMDAVecGetArray(da,  les->lQN, &QN);
            DMDAVecGetArray(da,  les->lNN, &NN);

            for (k=lzs; k<lze; k++)
            for (j=lys; j<lye; j++)
            for (i=lxs; i<lxe; i++)
            {
                // set to zero if solid and skip
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    QN[k][j][i] = NN[k][j][i] = 0;
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
                PetscReal Qij[3][3],                                   // Germano tensor projected on generalized coordinates
                        Sij_hat[3][3],                               // filtered strain rate tensor S projected on generalized coordinates
                        S_hat,                                       // scalar |S| = sqrt(2*Sij*Sij)
                        SSij_hat[3][3],                              // filtered |S|*S tensor
                        Nij[3][3],                                   // Nij tensor
                        Qij_cat[3][3],                               // cartesian Germano tensor
                        Nij_cat[3][3];                               // cartesian Nij tensor

                PetscReal filter, test_filter;                         // grid filter and dynamic model filter

                // local tensors before the test-filtering (5x5x5 cell stencil)
                PetscReal S[5][5][5];                                  // scalar |S| = sqrt(2*Sij*Sij)
                PetscReal u[5][5][5] ,   v[5][5][5] ,   w[5][5][5] ,   // cartesian velocity
                        d[5][5][5];
                PetscReal U[5][5][5]   , V[5][5][5]   , W[5][5][5]   ; // contravariant fluxes
                PetscReal Uu[5][5][5]  , Uv[5][5][5]  , Uw[5][5][5]  ; //         |Uu Uv Uw|
                PetscReal Vu[5][5][5]  , Vv[5][5][5]  , Vw[5][5][5]  ; // tensor: |Vu Vv Vw|
                PetscReal Wu[5][5][5]  , Wv[5][5][5]  , Ww[5][5][5]  ; //         |Wu Wv Ww|
                PetscReal S11[5][5][5] , S12[5][5][5] , S13[5][5][5] , // strain rate tensor
                        S21[5][5][5] , S22[5][5][5] , S23[5][5][5],
                        S31[5][5][5] , S32[5][5][5] , S33[5][5][5] ;
                PetscReal SS11[5][5][5], SS12[5][5][5], SS13[5][5][5], // |S|*S tensor
                        SS21[5][5][5], SS22[5][5][5], SS23[5][5][5],
                        SS31[5][5][5], SS32[5][5][5], SS33[5][5][5];

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

                    // |S|
                    S[R][Q][P]    = Sabs[K][J][I];

                    // |S|*S
                    SS11[R][Q][P] = S11[R][Q][P]*S[R][Q][P], SS12[R][Q][P] = S12[R][Q][P]*S[R][Q][P], SS13[R][Q][P] = S13[R][Q][P]*S[R][Q][P];
                    SS21[R][Q][P] = S21[R][Q][P]*S[R][Q][P], SS22[R][Q][P] = S22[R][Q][P]*S[R][Q][P], SS23[R][Q][P] = S23[R][Q][P]*S[R][Q][P];
                    SS31[R][Q][P] = S31[R][Q][P]*S[R][Q][P], SS32[R][Q][P] = S32[R][Q][P]*S[R][Q][P], SS33[R][Q][P] = S33[R][Q][P]*S[R][Q][P];
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

                PetscReal _u = integrateTestfilterSimpson5x5(u, weight);
                PetscReal _v = integrateTestfilterSimpson5x5(v, weight);
                PetscReal _w = integrateTestfilterSimpson5x5(w, weight);
                PetscReal _d = 1;

                Qij[0][0] = integrateTestfilterSimpson5x5(Uu, weight) - _U*_u;
                Qij[0][1] = integrateTestfilterSimpson5x5(Uv, weight) - _U*_v;
                Qij[0][2] = integrateTestfilterSimpson5x5(Uw, weight) - _U*_w;
                Qij[1][0] = integrateTestfilterSimpson5x5(Vu, weight) - _V*_u;
                Qij[1][1] = integrateTestfilterSimpson5x5(Vv, weight) - _V*_v;
                Qij[1][2] = integrateTestfilterSimpson5x5(Vw, weight) - _V*_w;
                Qij[2][0] = integrateTestfilterSimpson5x5(Wu, weight) - _W*_u;
                Qij[2][1] = integrateTestfilterSimpson5x5(Wv, weight) - _W*_v;
                Qij[2][2] = integrateTestfilterSimpson5x5(Ww, weight) - _W*_w;

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

                SSij_hat[0][0] = integrateTestfilterSimpson5x5(SS11, weight);
                SSij_hat[0][1] = integrateTestfilterSimpson5x5(SS12, weight);
                SSij_hat[0][2] = integrateTestfilterSimpson5x5(SS13, weight);
                SSij_hat[1][0] = integrateTestfilterSimpson5x5(SS21, weight);
                SSij_hat[1][1] = integrateTestfilterSimpson5x5(SS22, weight);
                SSij_hat[1][2] = integrateTestfilterSimpson5x5(SS23, weight);
                SSij_hat[2][0] = integrateTestfilterSimpson5x5(SS31, weight);
                SSij_hat[2][1] = integrateTestfilterSimpson5x5(SS32, weight);
                SSij_hat[2][2] = integrateTestfilterSimpson5x5(SS33, weight);


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
                    // Mohammad: Nij first defined in cartesian and then transformed
                    Nij_cat[a][b]

                    // Seba: Mij readily defined in curvilinear coordinates
                    //Mij[a][b]
                    =
                    (
                        - pow( test_filter, 2. ) * S_hat * Sij_hat[a][b]
                        + pow( filter, 2. ) * SSij_hat[a][b]
                    );

                }
                // Mohammad: project the Mij
                Nij[0][0] = Nij_cat[0][0] * csi0 + Nij_cat[0][1] * csi1 + Nij_cat[0][2] * csi2;
                Nij[0][1] = Nij_cat[0][0] * eta0 + Nij_cat[0][1] * eta1 + Nij_cat[0][2] * eta2;
                Nij[0][2] = Nij_cat[0][0] * zet0 + Nij_cat[0][1] * zet1 + Nij_cat[0][2] * zet2;
                Nij[1][0] = Nij_cat[1][0] * csi0 + Nij_cat[1][1] * csi1 + Nij_cat[1][2] * csi2;
                Nij[1][1] = Nij_cat[1][0] * eta0 + Nij_cat[1][1] * eta1 + Nij_cat[1][2] * eta2;
                Nij[1][2] = Nij_cat[1][0] * zet0 + Nij_cat[1][1] * zet1 + Nij_cat[1][2] * zet2;
                Nij[2][0] = Nij_cat[2][0] * csi0 + Nij_cat[2][1] * csi1 + Nij_cat[2][2] * csi2;
                Nij[2][1] = Nij_cat[2][0] * eta0 + Nij_cat[2][1] * eta1 + Nij_cat[2][2] * eta2;
                Nij[2][2] = Nij_cat[2][0] * zet0 + Nij_cat[2][1] * zet1 + Nij_cat[2][2] * zet2;

                PetscReal num = 0, num1 = 0, denom = 0;
                PetscInt  m, n, l;

                // compute numerator
                for(q=0; q<3; q++)
                for(a=0; a<3; a++)
                for(b=0; b<3; b++)
                {
                    // this contraction is given by
                    // Q G N = Qba * (Naq * Gbq)
                    // Note: for how tensors are build the matrix producti is row by row and
                    //       not row by colums as usual. In both products.
                    num += Qij[b][a] * Nij[a][q] * G[b][q];
                }

                // compute denominator
                for(m=0; m<3; m++)
                for(n=0; n<3; n++)
                for(l=0; l<3; l++)
                {
                    // this contraction is given by
                    // N G N = Nnm * (Nnl * Gml)
                    // Note: matrix product is row by row in the product between parenthesis
                    //       and row by column between the first matrix and the result of previous one.
                    denom += Nij[n][m] * Nij[n][l] * G[m][l];

                }

                // store in LM and MM
                QN[k][j][i] = num;
                NN[k][j][i] = denom;
            }

            DMDAVecRestoreArray(da,  les->lQN, &QN);
            DMDAVecRestoreArray(da,  les->lNN, &NN);

            DMLocalToLocalBegin(da, les->lQN, INSERT_VALUES, les->lQN);
            DMLocalToLocalEnd  (da, les->lQN, INSERT_VALUES, les->lQN);
            DMLocalToLocalBegin(da, les->lNN, INSERT_VALUES, les->lNN);
            DMLocalToLocalEnd  (da, les->lNN, INSERT_VALUES, les->lNN);
        }

        PetscReal ***LM_old, ***MM_old, ***QN_old, ***NN_old;

        // lagrangian - local convective weighting
        if(les->model == DLASI || les->model == DLASD)
        {
            PetscInt initializeLes4;
            if (clock->it <= clock->itStart + 1) initializeLes4 = 1;
            else                                 initializeLes4 = 0;

            DMDAVecGetArray(da, les->lLM, &LM);
            DMDAVecGetArray(da, les->lMM, &MM);

            DMDAVecGetArray(da, les->lLM_old, &LM_old);
            DMDAVecGetArray(da, les->lMM_old, &MM_old);

            // use standard Smagorinsky for first iteration (since need old values of LM and MM)
            if(initializeLes4)
            {
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

                VecSet(les->lCs, std_cs);
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

                    // Meneveau time scale
                    // PetscReal T_scale
                    // =
                    // PetscMax
                    // (
                    //     2.0*filter*
                    //     pow
                    //     (
                    //         fabs(LM_old[k][j][i]),
                    //         -1./4.
                    //     ),
                    //     1e-20
                    // );

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

                    // Bou-Zaid scale dependent model time scale
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

                    PetscReal _LM_old, _MM_old;

                    // interpolate the value at the exat old point (slower)
                    // if there is an IBM or overset cell around the given cell, use nearest cell
                    if(isBoxIBMCell(k_old, j_old, i_old, nvert)|| isBoxOversetCell(k_old, j_old, i_old, meshTag))
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
                            _LM_old = LM_old[k_old][j_old][i_old];
                            _MM_old = MM_old[k_old][j_old][i_old];
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

        if(les->model == DLASD)
        {
            PetscInt initializeLes4;
            if (clock->it <= clock->itStart + 1) initializeLes4 = 1;
            else                                 initializeLes4 = 0;

            DMDAVecGetArray(da, les->lQN, &QN);
            DMDAVecGetArray(da, les->lNN, &NN);

            DMDAVecGetArray(da, les->lQN_old, &QN_old);
            DMDAVecGetArray(da, les->lNN_old, &NN_old);

            // use standard Smagorinsky for first iteration (since need old values of QN and NN)
            if(initializeLes4)
            {
                for (k=lzs; k<lze; k++)
                {
                    for (j=lys; j<lye; j++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            QN_old[k][j][i] = std_cs*NN[k][j][i];
                            NN_old[k][j][i] = NN[k][j][i];
                        }
                    }
                }

                VecSet(les->lCs, std_cs);
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

                    // Meneveau time scale
                    // PetscReal T_scale
                    // =
                    // PetscMax
                    // (
                    //     2.0*filter*
                    //     pow
                    //     (
                    //         fabs(QN_old[k][j][i]),
                    //         -1./4.
                    //     ),
                    //     1e-20
                    // );

                    // Piomelli time scale
                    /*
                    PetscReal T_scale
                    =
                    PetscMax
                    (
                        1.5*filter*
                        pow
                        (
                            8.*fabs(QN_old[k][j][i]*NN_old[k][j][i]),
                            -1./8.
                        ),
                        1e-15
                    );
                    */

                    // Bou-Zaid scale dependent model time scale
                    PetscReal T_scale
                    =
                    1.5*filter*
                    pow
                    (
                        fabs(QN_old[k][j][i]*NN_old[k][j][i]) +
                        1.e-19,
                        -1.0/8.0
                    ) +
                    1.e-19;
                    

                    PetscReal QN_new, NN_new;

                    QN_new = QN[k][j][i];
                    NN_new = NN[k][j][i];

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

                    PetscReal _QN_old, _NN_old;

                    // interpolate the value at the exact old point (slower)
                    // if there is an IBM or overset cell around the given cell, use nearest cell
                    if(isBoxIBMCell(k_old, j_old, i_old, nvert)|| isBoxOversetCell(k_old, j_old, i_old, meshTag))
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
                            _QN_old = QN_old[k_old][j_old][i_old];
                            _NN_old = NN_old[k_old][j_old][i_old];
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
                                QN_old,
                                _QN_old
                            );

                            scalarPointLocalVolumeInterpolation
                            (
                                mesh,
                                X_old.x, X_old.y, X_old.z,
                                i_old, j_old, k_old,
                                cent,
                                NN_old,
                                _NN_old
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
                              QN_old,
                              _QN_old
                          );

                          scalarPointLocalVolumeInterpolation
                          (
                              mesh,
                              X_old.x, X_old.y, X_old.z,
                              i_old, j_old, k_old,
                              cent,
                              NN_old,
                              _NN_old
                          );
                      }
                      else if(lagrangianAverageCs == "nearestPoint")
                      {
                          // take the closest cell value (faster, but a lot of cs clipping from below)
                          _QN_old = QN_old[k_old][j_old][i_old];
                          _NN_old = NN_old[k_old][j_old][i_old];
                      }

                    }

                    // interpolation weights depend on the time scale.
                    // if dt >> T_scale new value is used
                    // if dt << T_scale old value is used
                    PetscReal eps = (clock->dt/T_scale) / (1.0 + clock->dt/T_scale);

                    // apply lagrangian weighting (relaxation)
                    QN[k][j][i]
                    =
                    PetscMax
                    (
                        (
                            (1-eps) * _QN_old +
                            eps * QN_new
                        ),
                        0.0
                    );

                    NN[k][j][i]
                    =
                    (1-eps) * _NN_old +
                    eps * NN_new;
                }

                // store QN and NN old values
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    QN_old[k][j][i] = QN[k][j][i];
                    NN_old[k][j][i] = NN[k][j][i];
                }
            }

            DMDAVecRestoreArray(da, les->lQN, &QN);
            DMDAVecRestoreArray(da, les->lNN, &NN);
            DMDAVecRestoreArray(da, les->lQN_old, &QN_old);
            DMDAVecRestoreArray(da, les->lNN_old, &NN_old);

            DMLocalToLocalBegin(da, les->lQN_old, INSERT_VALUES, les->lQN_old);
            DMLocalToLocalEnd  (da, les->lQN_old, INSERT_VALUES, les->lQN_old);
            DMLocalToLocalBegin(da, les->lNN_old, INSERT_VALUES, les->lNN_old);
            DMLocalToLocalEnd  (da, les->lNN_old, INSERT_VALUES, les->lNN_old);

            // update ghost values
            DMLocalToLocalBegin(da, les->lQN, INSERT_VALUES, les->lQN);
            DMLocalToLocalEnd  (da, les->lQN, INSERT_VALUES, les->lQN);
            DMLocalToLocalBegin(da, les->lNN, INSERT_VALUES, les->lNN);
            DMLocalToLocalEnd  (da, les->lNN, INSERT_VALUES, les->lNN);

            resetCellPeriodicFluxes(mesh, les->lQN, les->lQN, "scalar", "localToLocal");
            resetCellPeriodicFluxes(mesh, les->lNN, les->lNN, "scalar", "localToLocal");
        }

        // update ghost values
        DMLocalToLocalBegin(da, les->lLM, INSERT_VALUES, les->lLM);
        DMLocalToLocalEnd  (da, les->lLM, INSERT_VALUES, les->lLM);
        DMLocalToLocalBegin(da, les->lMM, INSERT_VALUES, les->lMM);
        DMLocalToLocalEnd  (da, les->lMM, INSERT_VALUES, les->lMM);

        resetCellPeriodicFluxes(mesh, les->lLM, les->lLM, "scalar", "localToLocal");
        resetCellPeriodicFluxes(mesh, les->lMM, les->lMM, "scalar", "localToLocal");

        DMDAVecRestoreArray(fda, les->lSx, &Ax);
        DMDAVecRestoreArray(fda, les->lSy, &Ay);
        DMDAVecRestoreArray(fda, les->lSz, &Az);
        DMDAVecRestoreArray(da,  les->lS,  &Sabs);
    }

    // 3 - Compute Cs
    if(les->model == STABILITY_BASED)
    {
        // reset to zero for eventual IBM
        VecSet(les->lCh, 0.0);

        DMDAVecGetArray(da,  les->lS,  &Sabs);
        DMDAVecGetArray(da,  les->lN,  &n);
        DMDAVecGetArray(da,  les->lCh, &lch);
        DMDAVecGetArray(da,  les->L,   &scale);
    }
    else if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD || les->model == BDS)
    {
        DMDAVecGetArray(da, les->lLM, &LM);
        DMDAVecGetArray(da, les->lMM, &MM);

        if(les->model == DLASD || les->model == DPASD)
        {
            DMDAVecGetArray(da,  les->lQN, &QN);
            DMDAVecGetArray(da,  les->lNN, &NN);
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
            Cs[k][j][i] = 0.0;
            continue;
        }

        // Cs value at current point
        PetscReal C = 0.0;

        // stability dependent model
        if(les->model == STABILITY_BASED)
        {
            // model coefficients
            PetscReal cm_les2  = 0.1,
                      ch1_les2 = 0.1,
                      ch2_les2 = 0.2,
                      ce1_les2 = 0.225,
                      ce2_les2 = 0.705,
                      cl_les2  = 0.82,
                      Ric_les2 = 0.23,
                      Ris_les2 = 0.059,
                      cet_les2 = ce2_les2 + ch2_les2*pow(cl_les2,2.0),
                      cht_les2 = ce1_les2 / pow(cl_les2,2.0) + ch1_les2;

            // filter delta
            PetscReal delta = pow( 1.0 / aj[k][j][i], 1.0/3.0);

            // Richardson number
            PetscReal Ri = pow(n[k][j][i],2.0) / pow(Sabs[k][j][i],2.0);

            // length scale based on richardson number
            PetscReal l, ch, ce, cs;

            if(Ri >= Ric_les2)
            {
                // above the BL: Ri can be huge (inf)
                l  = 1e-10;
            }
            else if (Ri <= Ris_les2)
            {
                // neutral or high stress regions
                l  = delta;
            }
            else
            {
                // stability regions where length scale decreases
                PetscReal e  = pow(cm_les2 * cl_les2 * delta * Sabs[k][j][i] * (1.0 - Ri / Ric_les2)  / cet_les2, 2.0) / Ri ;
                PetscReal l  = cl_les2 * sqrt(e) / n[k][j][i];
            }

            PetscReal walldist = cent[k][j][i].z - mesh->bounds.zmin + les->access->abl->hRough;
            l  = sqrt(1.0 / ( 1.0/pow(delta,2.0) + 1.0/pow(l,2.0) + 1.0/pow((0.41*walldist),2.0) ));
            ch = ch1_les2 + ch2_les2 * l / delta;
            ce = ce1_les2 + ce2_les2 * l / delta;
            cs = pow(pow(cm_les2, 3.0) / ce, 0.25);

            // save ch value for temperature equation
            lch[k][j][i] = ch;

            // save scale for visualization
            scale[k][j][i] = l;

            // compute C coefficient
            if(Ri >= Ric_les2) C = 0.0;
            else               C = pow(cs*l,2.0)*sqrt(1.0-ch/cm_les2*Ri);

            // Note: l would have to be saved. To save memory we incorporate it into Cs and divide by delta,
            //       so that delta will be simplified later and l will be taken as the length scale.
            C = C / pow (delta, 2.0);
        }
        // dynamic Smagorinsky models
        else if(les->model == DSM || les->model == DLASI || les->model == BDS)
        {
            PetscReal LM_avg, MM_avg;

            // box averaging
            if(les->model == DSM || les->model == BDS)
            {
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

                    LM0[R][Q][P] = LM[K][J][I];
                    MM0[R][Q][P] = MM[K][J][I];
                }

                // integrate over the testfilter volume locally
                LM_avg = integrateTestfilterSimpson(LM0, weight);
                MM_avg = integrateTestfilterSimpson(MM0, weight);
            }
            // lagrangian averaging
            else if(les->model == DLASI)
            {
                // no average procedure
                LM_avg = LM[k][j][i];
                MM_avg = MM[k][j][i];
            }

            // set Smagorinsky coefficient
            C = 0.5 * LM_avg / (MM_avg + 1e-10);
        }
        else if(les->model == DLASD)
        {
            PetscReal num, denom, beta;

            num   = LM[k][j][i]/(MM[k][j][i] + 1e-12);

            //ratio of Cs^2 at different scale
            beta  = (QN[k][j][i] * MM[k][j][i])/(NN[k][j][i] * LM[k][j][i] + 1e-12);

            denom = PetscMin(PetscMax(0.125, beta), 100.0);

            C = 0.5 * num/denom;
        }

        // set Cs to zero if solid or on box corners
        if( isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
        {
            Cs[k][j][i] = 0;
        }
        else
        {

            // clip Cs between upper and lower bound
            Cs[k][j][i] = PetscMin( PetscMax(C, 0.0), les->maxCs);
        }
    }

    if(les->model == DPASD)
    {
        std::vector<PetscInt> count, total_count;
        std::vector<PetscReal> J_LM(my), J_MM(my), J_QN(my), J_NN(my), LM_tmp(my), MM_tmp(my), QN_tmp(my), NN_tmp(my);
        
        count.resize(my);
        total_count.resize(my);
        
        for(j=0; j<my; j++) 
        {
            LM_tmp[j] = 0;
            MM_tmp[j] = 0;
            QN_tmp[j] = 0;
            NN_tmp[j] = 0;
            count[j]  = 0;
            total_count[j] = 0;
        }

        // loop over internal cells
        for (k=lzs; k<lze; k++)
        for (j=lys; j<lye; j++)
        for (i=lxs; i<lxe; i++)
        {
            if(isFluidCell(k,j,i, nvert) && isCalculatedCell(k, j, i, meshTag)) 
            {
                LM_tmp[j] += LM[k][j][i];
                MM_tmp[j] += MM[k][j][i];
                QN_tmp[j] += QN[k][j][i];
                NN_tmp[j] += NN[k][j][i];            
                count[j] ++;
            }
        }
        
        MPI_Allreduce( &LM_tmp[0], &J_LM[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
        MPI_Allreduce( &MM_tmp[0], &J_MM[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
        MPI_Allreduce( &QN_tmp[0], &J_QN[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
        MPI_Allreduce( &NN_tmp[0], &J_NN[0], my, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
        MPI_Allreduce( &count[0], &total_count[0], my, MPIU_INT, MPI_SUM, mesh->MESH_COMM);    
        
        for(j=0; j<my; j++) 
        {
            if( total_count[j]>0) 
            {
                J_LM[j] = J_LM[j]/total_count[j];
                J_MM[j] = J_MM[j]/total_count[j];
                J_QN[j] = J_QN[j]/total_count[j];
                J_NN[j] = J_NN[j]/total_count[j];
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
                Cs[k][j][i] = 0.0;
                continue;
            }

            beta = (J_QN[j] * J_MM[j])/(J_NN[j] * J_LM[j] +  1e-10);
            Cs[k][j][i] = 0.5 * (J_LM[j] / ( J_MM[j] + 1e-10))/beta;

            // clip Cs between upper and lower bound
            Cs[k][j][i] = PetscMin( PetscMax(Cs[k][j][i], 0.0), les->maxCs);
        }
    }

    if(les->model == STABILITY_BASED)
    {
        DMDAVecRestoreArray(da,  les->lS,  &Sabs);
        DMDAVecRestoreArray(da,  les->lN,  &n);
        DMDAVecRestoreArray(da,  les->lCh, &lch);
        DMDAVecRestoreArray(da,  les->L,   &scale);

        DMLocalToLocalBegin(da,  les->lCh,  INSERT_VALUES, les->lCh);
        DMLocalToLocalEnd  (da,  les->lCh,  INSERT_VALUES, les->lCh);
    }
    else if(les->model == DSM || les->model == DLASI || les->model == DLASD || les->model == DPASD || les->model == BDS)
    {
        DMDAVecRestoreArray(da, les->lLM, &LM);
        DMDAVecRestoreArray(da, les->lMM, &MM);

        if(les->model == DLASD || les->model == DPASD)
        {
            DMDAVecRestoreArray(da,  les->lQN, &QN);
            DMDAVecRestoreArray(da,  les->lNN, &NN);
        }
    }

    // restore numerator and denominator of the model coefficient
    DMDAVecRestoreArray(da,  les->lCs, &Cs);

    // restore fields
    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

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
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
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
    teqn_         *teqn = les->access->teqn;
    ueqn_         *ueqn = les->access->ueqn;
    constants_    *cst  = ueqn->access->constants;
    DM             da    = mesh->da, fda = mesh->fda, tda=mesh->tda;
                                                   
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***Cs, ***lnu_t, ***nvert, ***meshTag, ***aj, ***lt, ***ksg, ***ltempDiff;
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
    VecSet(les->lKsgs, 0.);

    std::vector<PetscReal> tempJAvg;
    Vec tempMAvg; 

    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da,  mesh->lAj, &aj);

    DMDAVecGetArray(da,  les->lNu_t, &lnu_t);
    DMDAVecGetArray(da,  les->lKsgs, &ksg);

    DMDAVecGetArray(da,  les->lCs, &Cs);

    if (les->model == AMD || les->model == BAMD)
    {

        if (les->access->flags->isTeqnActive)
        {
            DMDAVecGetArray(da, teqn->lTmprt, &lt);
            
            if(les->access->flags->isAblActive)
            {
                //j plane temperature average 
                tempJAvg = jPlaneScalarMean(mesh, lt, my-2);
                VecDuplicate(teqn->lTmprt, &tempMAvg);
                VecSet(tempMAvg, 0.);  

                DMDAVecGetArray(da, tempMAvg, &ltempDiff);
    
                for (k=lzs; k<lze; k++)
                for (j=lys; j<lye; j++)
                for (i=lxs; i<lxe; i++)
                {
                    ltempDiff[k][j][i] = lt[k][j][i] - tempJAvg[j-1];
                }
    
                DMDAVecRestoreArray(da, tempMAvg, &ltempDiff);
                DMLocalToLocalBegin(da,  tempMAvg, INSERT_VALUES, tempMAvg);
                DMLocalToLocalEnd  (da,  tempMAvg, INSERT_VALUES, tempMAvg);

                DMDAVecGetArray(da, tempMAvg, &ltempDiff);
            }
        }
    }

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
        {
            lnu_t[k][j][i]=0;
            ksg[k][j][i]=0;
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
            ucat, nvert, meshTag,
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
            
            PetscReal gradU[3][3] = 
            {
                {du_dx, dv_dx, dw_dx}, 
                {du_dy, dv_dy, dw_dy}, 
                {du_dz, dv_dz, dw_dz}  
            };

            PetscReal gradU_hat[3][3];
            {
                PetscReal deltaCell[3] = {delta_x, delta_y, delta_z};
                for (PetscInt c = 0; c < 3; c++) 
                {
                    for (PetscInt a = 0; a < 3; a++) 
                    {
                        gradU_hat[c][a] = (deltaCell[c] / deltaCell[a]) * gradU[c][a];
                    }
                }
            }

            PetscReal S_hat[3][3];
            for (PetscInt c = 0; c < 3; c++) 
            {
                for (PetscInt a = 0; a < 3; a++) 
                {
                    S_hat[c][a] = 0.5 * (gradU_hat[c][a] + gradU_hat[a][c]);
                }
            }

            PetscReal num = 0.0, denom = 0.0;

            for (PetscInt a = 0; a < 3; a++)
            {
                for (PetscInt b = 0; b < 3; b++) 
                {
                    for (PetscInt c = 0; c < 3; c++)
                    {
                        num += (-gradU_hat[c][a] * gradU_hat[c][b] * S_hat[a][b]);
                    }
                }
            }

            for (PetscInt a = 0; a < 3; a++)
            {
                for (PetscInt c = 0; c < 3; c++)
                {
                    denom += gradU_hat[c][a] * gradU_hat[c][a];
                }
            }

            // add abl stratification contribution - not used - this term seems to affect the ABL velocity near the inversion height.
            if (les->access->flags->isTeqnActive)
            {
                if(les->access->flags->isAblActive)
                {
                    //compute eddy diffusivity 
                    PetscReal dtdc, dtde, dtdz;
                    PetscReal dt_dx, dt_dy, dt_dz;
                    PetscReal gradT_hat[3];

                    PetscReal gravity = 9.81;
                    PetscReal tRef = les->access->abl->tRef;
                    PetscReal beta = 0; //gravity/tRef;

                    Compute_dscalar_center
                    (
                        mesh,
                        i, j, k,
                        mx, my, mz,
                        ltempDiff, nvert, meshTag,
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

                    gradT_hat[0] = delta_x * dt_dx;
                    gradT_hat[1] = delta_y * dt_dy;
                    gradT_hat[2] = delta_z * dt_dz;
                    
                    PetscInt a = 2; //only the z component
                    for (PetscInt c = 0; c < 3; c++)
                    {
                        num += (beta/delta_z) * gradU_hat[c][a] * gradT_hat[c];
                    }
                }
                
            }

            PetscReal nu_e = num / (denom + 1e-10);

            nu_e = Cs[k][j][i] * delta2_harm * PetscMax(nu_e, 0.0);
            
            lnu_t[k][j][i] = nu_e;
        }
        else if (les->model == VREMAN || les->model == BV)
        {
            
            PetscReal V = 1.0 / aj[k][j][i];

            PetscReal delta  = pow(V, 1.0/3.0);

            PetscReal delta2 = delta * delta;
            
        
            PetscReal gradU[3][3] 
                                = 
                                {
                                    {du_dx, dv_dx, dw_dx},
                                    {du_dy, dv_dy, dw_dy},
                                    {du_dz, dv_dz, dw_dz}
                                };

            //Compute G_ij = delta^2 * beta_ij
           
            PetscReal G[3][3] 
                            = {
                                {0,0,0},
                                {0,0,0},
                                {0,0,0}
                              };

            for (PetscInt i_comp = 0; i_comp < 3; i_comp++)
            {
                for (PetscInt j_comp = 0; j_comp < 3; j_comp++)
                {
                    for (PetscInt k_dir = 0; k_dir < 3; k_dir++)
                    {
                        G[i_comp][j_comp] += gradU[k_dir][i_comp] * gradU[k_dir][j_comp];
                    }
                    // Multiply by delta^2
                    G[i_comp][j_comp] *= delta2;
                }
            }

            PetscReal Gnorm2 = 0.0;
            
            for (PetscInt p = 0; p < 3; p++)
            {
                for (PetscInt q = 0; q < 3; q++)
                {
                    Gnorm2 += G[p][q] * G[p][q];
                }
            }

            PetscReal Bbeta = 0.0;
            
            if (delta2 > 1e-20)  // safety check
            {
                Bbeta = Gnorm2 / (delta2*delta2);  // i.e. / delta^4
            }

            PetscReal traceG = (G[0][0] + G[1][1] + G[2][2]);
            
            PetscReal A_beta = 0.0;
            
            if (delta2 > 1e-20)
            {
                A_beta = traceG / delta2;
            }

 
            PetscReal eps  = 1.0e-12;  // small positive regularization
            
            PetscReal nu_t = 0.0;

            if ((A_beta + eps) != 0.0)
            {
                nu_t = Cs[k][j][i] * delta2 * PetscSqrtReal(Bbeta) / (A_beta + eps);
            }
            
            else
            {
                nu_t = 0.0;
            }

            // Ensure non-negative
            nu_t = PetscMax(nu_t, 0.0);

            
            lnu_t[k][j][i] = nu_t;


        }
        else
        {
            PetscReal filter = pow( 1./aj[k][j][i], 1./3.);
            lnu_t[k][j][i] = Cs[k][j][i] * pow(filter, 2.0) * Sabs;
        }

        //compute the subgrid scale kinetic energy - from local equilibrium balance (https://caefn.com/openfoam/smagorinsky-sgs-model)
        PetscReal Ck = 0.094;
        PetscReal V = 1.0 / aj[k][j][i]; 

        PetscReal A_x = nMag(zet[k][j][i]);
        PetscReal A_y = nMag(csi[k][j][i]);
        PetscReal A_z = nMag(eta[k][j][i]);

        PetscReal delta_x = V / A_x;
        PetscReal delta_y = V / A_y;
        PetscReal delta_z = V / A_z;

        PetscReal invDelta2Sum =
              1.0 / (delta_x * delta_x)
            + 1.0 / (delta_y * delta_y)
            + 1.0 / (delta_z * delta_z);

        PetscReal delta2_harm = 3.0 / invDelta2Sum; 

        ksg[k][j][i] = PetscPowReal(lnu_t[k][j][i], 2.0)/(PetscPowReal(Ck, 2.0) * delta2_harm);
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);

    DMDAVecRestoreArray(da,  les->lNu_t, &lnu_t);
    DMDAVecRestoreArray(da,  les->lKsgs, &ksg);
    DMDAVecRestoreArray(da,  les->lCs, &Cs);

    if (les->model == AMD || les->model == BAMD)
    {
        if (les->access->flags->isTeqnActive)
        {
            DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
            if(les->access->flags->isAblActive)
            {
                DMDAVecRestoreArray(da, tempMAvg, &ltempDiff);
                VecDestroy(&tempMAvg);
            }

            UpdateDiffBCs(les);
        }
    }

    DMLocalToLocalBegin(da,  les->lNu_t, INSERT_VALUES, les->lNu_t);
    DMLocalToLocalEnd  (da,  les->lNu_t, INSERT_VALUES, les->lNu_t);

    DMLocalToLocalBegin(da,  les->lKsgs, INSERT_VALUES, les->lKsgs);
    DMLocalToLocalEnd  (da,  les->lKsgs, INSERT_VALUES, les->lKsgs);
    // correct boundary conditions
    UpdateNutBCs(les);

    if(les->turbStat)
    {
        ComputeTurbulenceStatistics(les);
    }

    return (0);
}

//*******************************************************************************************
PetscErrorCode ComputeTurbulenceStatistics(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    teqn_         *teqn = les->access->teqn;
    ueqn_         *ueqn = les->access->ueqn;
    constants_    *cst  = ueqn->access->constants;
    DM             da    = mesh->da, fda = mesh->fda, tda=mesh->tda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***ksg, ***lnu_t, ***nvert, ***meshTag, ***aj;
    Cmpnts        ***csi, ***eta, ***zet, ***ucat;

    PetscMPIInt rank;

    PetscReal lEns, gEns, lVol1, gVol1,lVol2, gVol2, lSS, gSS, lW, gW, lVS, gVS, lSS_disp, gSS_disp;
    PetscReal lE, gE;
    PetscReal lUx, lUy, lUz, gUx, gUy, gUz;

    PetscReal luu, lvv, lww, luv,luw,lvw;
    
    PetscReal guu, gvv, gww, guv,guw,gvw;
    
    PetscReal lTauS, gTauS;
    PetscReal lSgs_disp, gSgs_disp;

    DMDAVecGetArray(da,  les->lNu_t, &lnu_t);
    DMDAVecGetArray(da,  les->lKsgs, &ksg);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da,  mesh->lAj, &aj);

    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    lEns = 0.0; gEns = 0.0;
    lSS  = 0.0; gSS  = 0.0;
    lSS_disp=0.0; gSS_disp = 0.0;
    lW   = 0.0; gW   = 0.0;
    lVS  = 0.0; gVS  = 0.0;
    lVol1 = 0.0; gVol1 = 0.0;
    lVol2 = 0.0; gVol2 = 0.0;
    lE   = 0.0; gE   =  0.0;
    lUx  = 0.0; gUx  = 0.0;
    lUy  = 0.0; gUy  = 0.0;
    lUz  = 0.0; gUz  = 0.0;
    luu = 0.0; guu = 0.0;
    lvv = 0.0; gvv = 0.0;
    lww = 0.0; gww = 0.0;
    luv  = 0.0; guv  = 0.0;
    luw  = 0.0; guw  = 0.0;
    lvw  = 0.0; gvw  = 0.0;
    lTauS = 0.0; gTauS = 0.0;
    lSgs_disp = 0.0, gSgs_disp = 0.0;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        PetscReal volCell = 1.0 / aj[k][j][i]; 

    
        lUx  += (ucat[k][j][i].x) * volCell;
        lUy  += (ucat[k][j][i].y) * volCell;
        lUz  += (ucat[k][j][i].z) * volCell;

        lVol1 += volCell;
    }

    MPI_Allreduce(&lVol1, &gVol1, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lUx,  &gUx,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lUy,  &gUy,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lUz,  &gUz,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    Cmpnts gVelocity;
    
    gVelocity.x = gUx / gVol1;
    gVelocity.y = gUy / gVol1;
    gVelocity.z = gUz / gVol1;

    for (k=lzs; k<lze; k++)
    for (j=lys; j<lye; j++)
    for (i=lxs; i<lxe; i++)
    {
        //storing time series data of volume averaged Enstrophy, kinetic energy, taylor micro-scale, and skewness
        PetscReal ajc = aj[k][j][i];
        PetscReal csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
        PetscReal eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
        PetscReal zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

        PetscReal dudc, dvdc, dwdc, dude, dvde, dwde, dudz, dvdz, dwdz;
        PetscReal du_dx, du_dy, du_dz, dv_dx, dv_dy, dv_dz, dw_dx, dw_dy, dw_dz;

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

        PetscReal volCell = 1.0/ajc;
        
        Tensor VGT;

        //intialize to 0
        zeroTensor(&VGT);

        VGT.xx = du_dx;  VGT.xy = dv_dx;  VGT.xz = dw_dx;
        VGT.yx = du_dy;  VGT.yy = dv_dy;  VGT.yz = dw_dy;
        VGT.zx = du_dz;  VGT.zy = dv_dz;  VGT.zz = dw_dz;

        //strain rate
        Tensor SS= symm(VGT);     
      
        Cmpnts omega_curl = vorticity(VGT);

        PetscReal omegaMagSqr = (omega_curl.x * omega_curl.x + omega_curl.y * omega_curl.y + omega_curl.z * omega_curl.z);

        // Compute enstrophy: 0.5 * (omega_x^2 + omega_y^2 + omega_z^2)
        PetscReal enstrophy = 0.5 * omegaMagSqr;

        PetscReal SS_mag = PetscSqrtReal(tensorInnerProduct(SS));

        PetscReal SS_disp = tensorInnerProduct(SS);

        PetscReal vs = omega_curl.x * (SS.xx*omega_curl.x + SS.xy*omega_curl.y + SS.xz*omega_curl.z)
                     + omega_curl.y * (SS.xy*omega_curl.x + SS.yy*omega_curl.y + SS.yz*omega_curl.z)
                     + omega_curl.z * (SS.xz*omega_curl.x + SS.yz*omega_curl.y + SS.zz*omega_curl.z);
        
        PetscReal u2 = 0.5 * (ucat[k][j][i].x * ucat[k][j][i].x 
                     + ucat[k][j][i].y * ucat[k][j][i].y 
                     + ucat[k][j][i].z * ucat[k][j][i].z);

        lEns += enstrophy    * volCell;
        lSS  += SS_mag       * volCell;
        lW   += omegaMagSqr  * volCell;
        lVS  += vs           * volCell;
        lE += u2             * volCell;
        lSS_disp += SS_disp  * volCell;
        lSgs_disp += 2*lnu_t[k][j][i] * SS_disp * volCell;

        PetscReal up = ucat[k][j][i].x - gVelocity.x;  
        PetscReal vp = ucat[k][j][i].y - gVelocity.y; 
        PetscReal wp = ucat[k][j][i].z - gVelocity.z; 
        
        luu += (up * up) * volCell; 
        lvv += (vp * vp) * volCell;
        lww += (wp * wp) * volCell;
        luv  += (up * vp) * volCell;
        luw  += (up * wp) * volCell;
        lvw  += (vp * wp) * volCell;

        symmTensor Rij;

        PetscReal nusgs = lnu_t[k][j][i];
        Rij.xx = (2.0/3.0 * ksg[k][j][i] - 2.0 * (nusgs) * du_dx); 
        Rij.yy = (2.0/3.0 * ksg[k][j][i] - 2.0 * (nusgs) * dv_dy);
        Rij.zz = (2.0/3.0 * ksg[k][j][i] - 2.0 * (nusgs) * dw_dz);
        Rij.xy = (- (nusgs) * (dv_dx + du_dy));
        Rij.xz = (- (nusgs) * (dw_dx + du_dz));
        Rij.yz = (- (nusgs) * (dv_dz + dw_dy));
        
        PetscReal tauSij = (Rij.xx) * SS.xx + (Rij.yy) * SS.yy + (Rij.zz) * SS.zz
        + 2.0 * ((Rij.xy) * SS.xy + (Rij.xz) * SS.xz + (Rij.yz) * SS.yz);

        lTauS += tauSij * volCell;

        lVol2 += volCell;
    }

    MPI_Allreduce(&lTauS, &gTauS, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&luu, &guu, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lvv, &gvv, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lww, &gww, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    
    MPI_Allreduce(&luv,  &guv,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&luw,  &guw,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lvw,  &gvw,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    
    MPI_Allreduce(&lE,   &gE, 1, MPIU_REAL, MPI_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lEns, &gEns, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lSS_disp,  &gSS_disp,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lSS,  &gSS,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lW,   &gW,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lVS,  &gVS,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lVol2, &gVol2, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lSgs_disp, &gSgs_disp, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    guu = guu / gVol2; 
    gvv = gvv / gVol2; 
    gww = gww / gVol2; 
    guv =  guv / gVol2; 
    guw = guw  / gVol2; 
    gvw = gvw  / gVol2; 
    gTauS = gTauS / gVol2;

    gEns = gEns / gVol2;  
    
    gSS  = gSS  / gVol2;
    gSS_disp = gSS_disp / gVol2; 
    gW   = gW   / gVol2; 
    gVS  = gVS  / gVol2; 
    gE    = gE / gVol2;
    gSgs_disp = gSgs_disp / gVol2;

    PetscReal urms = PetscSqrtReal(guu);  
    PetscReal vrms = PetscSqrtReal(gvv);
    PetscReal wrms = PetscSqrtReal(gww);

    //viscous dissipation
    PetscReal dispS = 2 * cst->nu * gSS_disp;
    
    PetscReal dispW  = cst->nu * gW;

    //kolmogorov scale
    PetscReal neta = PetscPowReal(cst->nu, 3.0/4.0)/PetscPowReal(dispS, 1.0/4.0);

    //Taylor-Microscale(davidson)
    // PetscReal TayL  = sqrt(15.0 * PetscPowReal(urms, 2.0)/ dispS);
        PetscReal TayL = sqrt(15.0 * gE / gSS_disp);
    
    //Taylor-Reynolds number
    PetscReal ReTayL = (urms * TayL) / cst->nu;
    
    //skewness
    PetscReal Sk = - ((6.0/7.0) * PetscSqrtReal(15.0) * gVS) / pow(gW, 3.0/2.0);

    //volume averaged resolved energy
    PetscReal k_res = 0.5 * (guu + gvv + gww);
    
    PetscPrintf(PETSC_COMM_WORLD, "TayL: %.5lf\n", TayL);
    PetscPrintf(PETSC_COMM_WORLD, "Sk:   %.5lf\n", Sk);

    // Write the source terms (time and volume-averaged enstrophy) to file.
    if (!rank) 
    {
        // Create the postProcessing directory on the first iteration.
        if (les->access->clock->it == les->access->clock->itStart) 
        {
            errno = 0;
            PetscInt dirRes = mkdir("./postProcessing", 0777);
            if (dirRes != 0 && errno != EEXIST) 
            {
                char error[512];
                sprintf(error, "could not create postProcessing directory\n");
                fatalErrorInFunction("turbStat", error);
            }
        }

        // Construct the file name.
        word fileName = "postProcessing/turbstat_" + getStartTimeName(les->access->clock);
        FILE *fp = fopen(fileName.c_str(), "a");

        if (!fp) 
        {
            char error[512];
            sprintf(error, "cannot open file postProcessing/turbstat file\n");
            fatalErrorInFunction("turbStat", error);
        } 
        else 
        {
            PetscInt width = -15;

            // Write the current time.
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, les->access->clock->time);
            // Write the volume-averaged enstrophy.
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, neta);   //kolmogorov scale
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, TayL);  //taylor microscale
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, ReTayL); //taylor reynolds number
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, dispS);  //dissipation with SijSij
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, dispW);  //dissipation with enstrophy
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, k_res);  //resolved kinetic energy
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gE);     //energy
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gTauS);  //Energy flux <tauijSij>
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, Sk);     //Skewness S0
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gEns);   //Enstrophy
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gSS_disp);    //SijSij
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gW);     //||omega^2||
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gVS);    //Vortex-stretching term
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.6e\t", width, gSgs_disp);    //sgs dissipation
            PetscFPrintf(mesh->MESH_COMM, fp, "\n");
            
            fclose(fp);
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);
    DMDAVecRestoreArray(da,  les->lNu_t, &lnu_t);
    DMDAVecRestoreArray(da,  les->lKsgs, &ksg);

    return (0);   
}
//*******************************************************************************************
PetscErrorCode updateLESStructuralModelCartesianForm(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    ueqn_         *ueqn = les->access->ueqn;
    teqn_         *teqn = les->access->teqn;

    DM            da    = mesh->da, fda = mesh->fda, tda = mesh->tda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***aj, ***Sabs, ***nut,  ***ksg, ***nvert, ***meshTag;
    PetscReal     ***iaj, ***jaj, ***kaj;
    Cmpnts        ***csi, ***eta, ***zet, ***ucat;
    Cmpnts        ***visc1, ***visc2, ***visc3;
    Cmpnts        ***icsi, ***jeta, ***kzet;

    PetscReal     ajc;
    PetscMPIInt   rank;

    PetscReal     dudc, dude, dudz,
                  dvdc, dvde, dvdz,
                  dwdc, dwde, dwdz;

    Tensor        ***tauSGSCat, tauLocalCat;                                      //subgrid scale stress tensor in cartesian form
    symmTensor    SS;

    PetscReal     du_dx, du_dy, du_dz, 
                  dv_dx, dv_dy, dv_dz, 
                  dw_dx, dw_dy, dw_dz;

    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;      // surface area vectors components

    lxs = xs; lxe = xe; if (xs == 0) lxs = xs+1; if (xe == mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys == 0) lys = ys+1; if (ye == my) lye = ye-1;
    lzs = zs; lze = ze; if (zs == 0) lzs = zs+1; if (ze == mz) lze = mz-1;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    VecSet(les->lTau, 0.0);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da,  les->lNu_t, &nut);
    DMDAVecGetArray(da,  les->lKsgs, &ksg);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da, mesh->lAj,  &aj);

    DMDAVecGetArray(tda, les->lTau, &tauSGSCat);

    PetscReal lVol2 = 0., gVol2 = 0., lTauS = 0., gTauS = 0.;

    //bardina model
    if (les->model == BAMD || les->model == BV || les->model == BDS)
    {

        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {
            /* Skip solid cells */
            if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                zeroTensor(&tauSGSCat[k][j][i]); //inline function each element 0
                continue;
            }

            // get 1/V at the i-face
            ajc = aj[k][j][i];

            csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
            eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
            zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

            PetscReal filter, test_filter;                         // grid filter and dynamic model filter

            PetscReal u[3][3][3] ,   v[3][3][3] ,   w[3][3][3];    // cartesian velocity
            
            PetscReal uu[3][3][3]  , uv[3][3][3]  , uw[3][3][3];   //         |uu uv uw|
            
            PetscReal vu[3][3][3]  , vv[3][3][3]  , vw[3][3][3];   // tensor: |vu vv vw|
            
            PetscReal wu[3][3][3]  , wv[3][3][3]  , ww[3][3][3];   //         |wu wv ww|
            
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
            
                //uu tensor components
                
                uu[R][Q][P] = u[R][Q][P] * u[R][Q][P];
                uv[R][Q][P] = u[R][Q][P] * v[R][Q][P];
                uw[R][Q][P] = u[R][Q][P] * w[R][Q][P];

                vu[R][Q][P] = v[R][Q][P] * u[R][Q][P];
                vv[R][Q][P] = v[R][Q][P] * v[R][Q][P];
                vw[R][Q][P] = v[R][Q][P] * w[R][Q][P];

                wu[R][Q][P] = w[R][Q][P] * u[R][Q][P];
                wv[R][Q][P] = w[R][Q][P] * v[R][Q][P];
                ww[R][Q][P] = w[R][Q][P] * w[R][Q][P];
            }

            //Build filter weights (using a binomial 3x3x3 kernel)
            
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

            PetscReal weight[3][3][3];

            PetscReal sum_vol = 0;

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

                    sum_vol += (1.0/aj[K][J][I]) * 8.0 * coef[R][Q][P];
                    weight[R][Q][P] = 1;
                
                } 
                
                else 
                {
                    weight[R][Q][P] = 0;
                }
            
            }
            
            //Apply test filter to local velocity fields
            
            PetscReal _u = integrateTestfilterSimpson(u, weight);
            PetscReal _v = integrateTestfilterSimpson(v, weight);
            PetscReal _w = integrateTestfilterSimpson(w, weight);
            
            //Apply test filter to product
            PetscReal _uu = integrateTestfilterSimpson(uu, weight);
            PetscReal _uv = integrateTestfilterSimpson(uv, weight);
            PetscReal _uw = integrateTestfilterSimpson(uw, weight);
            
            PetscReal _vu = integrateTestfilterSimpson(vu, weight);
            PetscReal _vv = integrateTestfilterSimpson(vv, weight);
            PetscReal _vw = integrateTestfilterSimpson(vw, weight);

            PetscReal _wu = integrateTestfilterSimpson(wu, weight);
            PetscReal _wv = integrateTestfilterSimpson(wv, weight);
            PetscReal _ww = integrateTestfilterSimpson(ww, weight);

            //Compute the Bardina SGS stress ---
            //tau_ij = test-filtered(u_i*u_j) - (test-filtered(u_i))*(test-filtered(u_j))
            //and store the full 3x3 tensor in tauSGSCat.
            
            tauSGSCat[k][j][i].xx = _uu - _u*_u;   // tau_xx
            tauSGSCat[k][j][i].xy = _uv - _u*_v;   // tau_xy
            tauSGSCat[k][j][i].xz = _uw - _u*_w;   // tau_xz
            
            tauSGSCat[k][j][i].yx = _vu - _v*_u;   // tau_yx
            tauSGSCat[k][j][i].yy = _vv - _v*_v;   // tau_yy
            tauSGSCat[k][j][i].yz = _vw - _v*_w;   // tau_yz
            
            tauSGSCat[k][j][i].zx = _wu - _w*_u;   // tau_zx
            tauSGSCat[k][j][i].zy = _wv - _w*_v;   // tau_zy
            tauSGSCat[k][j][i].zz = _ww - _w*_w;   // tau_zz
            
           /*( if(ueqn->central4Div)
            {
                Compute_du_center4th
                (
                    mesh,
                    i, j, k,
                    mx, my, mz,
                    ucat, nvert,
                    &dudc, &dvdc, &dwdc,
                    &dude, &dvde, &dwde,
                    &dudz, &dvdz, &dwdz
                );
            }
            else*/
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

            Compute_du_dxyz
            (
                mesh, csi0, csi1, csi2,
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

            //compute strain rate tensor
            SS.xx = 0.5 * (du_dx + du_dx);   
            SS.xy = 0.5 * (du_dy + dv_dx);
            SS.xz = 0.5 * (du_dz + dw_dx);

            SS.yy = 0.5 * (dv_dy + dv_dy);   
            SS.yz = 0.5 * (dv_dz + dw_dy);

            SS.zz = 0.5 * (dw_dz + dw_dz);  

            //kinetic energy: trace of tauSGSCat
            ksg[k][j][i] = 0.5 * trace(tauSGSCat[k][j][i]);   

            // provide dissipation to bardina model
            tauSGSCat[k][j][i].xx -= 2.0 * nut[k][j][i] * SS.xx;
            tauSGSCat[k][j][i].xy -= 2.0 * nut[k][j][i] * SS.xy;
            tauSGSCat[k][j][i].xz -= 2.0 * nut[k][j][i] * SS.xz;

            tauSGSCat[k][j][i].yx -= 2.0 * nut[k][j][i] * SS.xy; 
            tauSGSCat[k][j][i].yy -= 2.0 * nut[k][j][i] * SS.yy;
            tauSGSCat[k][j][i].yz -= 2.0 * nut[k][j][i] * SS.yz;

            tauSGSCat[k][j][i].zx -= 2.0 * nut[k][j][i] * SS.xz; 
            tauSGSCat[k][j][i].zy -= 2.0 * nut[k][j][i] * SS.yz; 
            tauSGSCat[k][j][i].zz -= 2.0 * nut[k][j][i] * SS.zz;
            
            PetscReal tauSij = (tauSGSCat[k][j][i].xx) * SS.xx + (tauSGSCat[k][j][i].yy) * SS.yy + (tauSGSCat[k][j][i].zz) * SS.zz
            + 2.0 * ((tauSGSCat[k][j][i].xy) * SS.xy + (tauSGSCat[k][j][i].xz) * SS.xz + (tauSGSCat[k][j][i].yz) * SS.yz);
        
            lTauS += tauSij * (1.0/aj[k][j][i]);
            lVol2 += (1.0/aj[k][j][i]);
        } 
    }

    DMDAVecRestoreArray(tda, les->lTau, &tauSGSCat);
    DMDAVecRestoreArray(da,  les->lKsgs, &ksg);
    DMDAVecRestoreArray(da,  les->lNu_t, &nut);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da, mesh->lAj,  &aj);

    DMLocalToLocalBegin(tda, les->lTau, INSERT_VALUES, les->lTau);
    DMLocalToLocalEnd  (tda, les->lTau, INSERT_VALUES, les->lTau);

    //set tau boundary conditions
    DMDAVecGetArray(tda, les->lTau, &tauSGSCat);

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
                    tauSGSCat[k][j][i].xx = tauSGSCat[c][b][a].xx;
                    tauSGSCat[k][j][i].xy = tauSGSCat[c][b][a].xy;
                    tauSGSCat[k][j][i].xz = tauSGSCat[c][b][a].xz;
                    tauSGSCat[k][j][i].yx = tauSGSCat[c][b][a].yx;
                    tauSGSCat[k][j][i].yy = tauSGSCat[c][b][a].yy;
                    tauSGSCat[k][j][i].yz = tauSGSCat[c][b][a].yz;
                    tauSGSCat[k][j][i].zx = tauSGSCat[c][b][a].zx;
                    tauSGSCat[k][j][i].zy = tauSGSCat[c][b][a].zy;
                    tauSGSCat[k][j][i].zz = tauSGSCat[c][b][a].zz;
                }
            }
        }
    }

    DMDAVecRestoreArray(tda, les->lTau, &tauSGSCat);
    DMLocalToLocalBegin(tda, les->lTau, INSERT_VALUES, les->lTau);
    DMLocalToLocalEnd  (tda, les->lTau, INSERT_VALUES, les->lTau);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);
    DMDAVecGetArray(tda, les->lTau, &tauSGSCat);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    //i direction
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

                tauLocalCat.xx = 0.5 * (tauSGSCat[k][j][i].xx + tauSGSCat[k][j][i+1].xx);
                tauLocalCat.xy = 0.5 * (tauSGSCat[k][j][i].xy + tauSGSCat[k][j][i+1].xy);
                tauLocalCat.xz = 0.5 * (tauSGSCat[k][j][i].xz + tauSGSCat[k][j][i+1].xz);
                tauLocalCat.yx = 0.5 * (tauSGSCat[k][j][i].yx + tauSGSCat[k][j][i+1].yx);
                tauLocalCat.yy = 0.5 * (tauSGSCat[k][j][i].yy + tauSGSCat[k][j][i+1].yy);
                tauLocalCat.yz = 0.5 * (tauSGSCat[k][j][i].yz + tauSGSCat[k][j][i+1].yz);
                tauLocalCat.zx = 0.5 * (tauSGSCat[k][j][i].zx + tauSGSCat[k][j][i+1].zx);
                tauLocalCat.zy = 0.5 * (tauSGSCat[k][j][i].zy + tauSGSCat[k][j][i+1].zy);
                tauLocalCat.zz = 0.5 * (tauSGSCat[k][j][i].zz + tauSGSCat[k][j][i+1].zz);

                if
                (
                    (mesh->boundaryU.iLeft=="velocityWallFunction"  && i==0) ||
                    (mesh->boundaryU.iRight=="velocityWallFunction" && i==mx-2)
                )
                {
                    visc1[k][j][i].x += (- ueqn->iLWM->tauWall.x[k-zs][j-ys]);
                    visc1[k][j][i].y += (- ueqn->iLWM->tauWall.y[k-zs][j-ys]);
                    visc1[k][j][i].z += (- ueqn->iLWM->tauWall.z[k-zs][j-ys]);

                }
                else if
                (
                    (mesh->boundaryU.iLeft =="slip" && i==0   ) ||
                    (mesh->boundaryU.iRight=="slip" && i==mx-2)
                )
                {
                    visc1[k][j][i].x = 0.0;
                    visc1[k][j][i].y = 0.0;
                    visc1[k][j][i].z = 0.0;
                }
                else 
                {
                    visc1[k][j][i].x -= (tauLocalCat.xx * csi0 + tauLocalCat.xy * csi1 + tauLocalCat.xz * csi2);
                    visc1[k][j][i].y -= (tauLocalCat.yx * csi0 + tauLocalCat.yy * csi1 + tauLocalCat.yz * csi2);
                    visc1[k][j][i].z -= (tauLocalCat.zx * csi0 + tauLocalCat.zy * csi1 + tauLocalCat.zz * csi2);
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

                // get 1/V at the j-face
                PetscReal ajc = jaj[k][j][i];

                // get face normals
                eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;

                tauLocalCat.xx = 0.5 * (tauSGSCat[k][j][i].xx + tauSGSCat[k][j+1][i].xx);
                tauLocalCat.xy = 0.5 * (tauSGSCat[k][j][i].xy + tauSGSCat[k][j+1][i].xy);
                tauLocalCat.xz = 0.5 * (tauSGSCat[k][j][i].xz + tauSGSCat[k][j+1][i].xz);
                tauLocalCat.yx = 0.5 * (tauSGSCat[k][j][i].yx + tauSGSCat[k][j+1][i].yx);
                tauLocalCat.yy = 0.5 * (tauSGSCat[k][j][i].yy + tauSGSCat[k][j+1][i].yy);
                tauLocalCat.yz = 0.5 * (tauSGSCat[k][j][i].yz + tauSGSCat[k][j+1][i].yz);
                tauLocalCat.zx = 0.5 * (tauSGSCat[k][j][i].zx + tauSGSCat[k][j+1][i].zx);
                tauLocalCat.zy = 0.5 * (tauSGSCat[k][j][i].zy + tauSGSCat[k][j+1][i].zy);
                tauLocalCat.zz = 0.5 * (tauSGSCat[k][j][i].zz + tauSGSCat[k][j+1][i].zz);

                if
                (
                    (mesh->boundaryU.jLeft=="velocityWallFunction" && j==0) ||
                    (mesh->boundaryU.jRight=="velocityWallFunction" && j==my-2)
                )
                {
                    visc2[k][j][i].x += (- ueqn->jLWM->tauWall.x[k-zs][i-xs]);
                    visc2[k][j][i].y += (- ueqn->jLWM->tauWall.y[k-zs][i-xs]);
                    visc2[k][j][i].z += (- ueqn->jLWM->tauWall.z[k-zs][i-xs]);
                }
                else if
                (
                    (mesh->boundaryU.jLeft =="slip" && j==0   ) ||
                    (mesh->boundaryU.jRight=="slip" && j==my-2)
                )
                {
                    visc2[k][j][i].x = 0.0;
                    visc2[k][j][i].y = 0.0;
                    visc2[k][j][i].z = 0.0;
                }
                else
                {
                    // if(k ==15 && j==2 && i==2)	
                    // {
                    //     PetscPrintf(PETSC_COMM_WORLD,"before %.10lf %.10lf %.10lf\n", visc2[k][j][i].x, visc2[k][j][i].y, visc2[k][j][i].z);
                    //     PetscPrintf(PETSC_COMM_WORLD,"taulocal %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", tauLocalCat.xx, tauLocalCat.yy, tauLocalCat.zz, tauLocalCat.xy, tauLocalCat.xz, tauLocalCat.yz);
                    //     PetscPrintf(PETSC_COMM_WORLD,"eta %.10lf %.10lf %.10lf\n", eta0, eta1, eta2);
                    // }
                    visc2[k][j][i].x -= (tauLocalCat.xx * eta0 + tauLocalCat.xy * eta1 + tauLocalCat.xz * eta2);
                    visc2[k][j][i].y -= (tauLocalCat.yx * eta0 + tauLocalCat.yy * eta1 + tauLocalCat.yz * eta2);
                    visc2[k][j][i].z -= (tauLocalCat.zx * eta0 + tauLocalCat.zy * eta1 + tauLocalCat.zz * eta2);
                    // if(k ==15 && j==2 && i==2)
                    // PetscPrintf(PETSC_COMM_WORLD,"after %.10lf %.10lf %.10lf\n", visc2[k][j][i].x, visc2[k][j][i].y, visc2[k][j][i].z);
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

                // get 1/V at the j-face
                PetscReal ajc = kaj[k][j][i];

                // get face normals
                zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                tauLocalCat.xx = 0.5 * (tauSGSCat[k][j][i].xx + tauSGSCat[k+1][j][i].xx);
                tauLocalCat.xy = 0.5 * (tauSGSCat[k][j][i].xy + tauSGSCat[k+1][j][i].xy);
                tauLocalCat.xz = 0.5 * (tauSGSCat[k][j][i].xz + tauSGSCat[k+1][j][i].xz);
                tauLocalCat.yx = 0.5 * (tauSGSCat[k][j][i].yx + tauSGSCat[k+1][j][i].yx);
                tauLocalCat.yy = 0.5 * (tauSGSCat[k][j][i].yy + tauSGSCat[k+1][j][i].yy);
                tauLocalCat.yz = 0.5 * (tauSGSCat[k][j][i].yz + tauSGSCat[k+1][j][i].yz);
                tauLocalCat.zx = 0.5 * (tauSGSCat[k][j][i].zx + tauSGSCat[k+1][j][i].zx);
                tauLocalCat.zy = 0.5 * (tauSGSCat[k][j][i].zy + tauSGSCat[k+1][j][i].zy);
                tauLocalCat.zz = 0.5 * (tauSGSCat[k][j][i].zz + tauSGSCat[k+1][j][i].zz);

                if
                (
                    (mesh->boundaryU.kLeft =="slip" && k==0   ) ||
                    (mesh->boundaryU.kRight=="slip" && k==mz-2)
                )
                {
                    visc3[k][j][i].x = 0.0;
                    visc3[k][j][i].y = 0.0;
                    visc3[k][j][i].z = 0.0;
                }
                else
                {
                    visc3[k][j][i].x -= (tauLocalCat.xx * zet0 + tauLocalCat.xy * zet1 + tauLocalCat.xz * zet2);
                    visc3[k][j][i].y -= (tauLocalCat.yx * zet0 + tauLocalCat.yy * zet1 + tauLocalCat.yz * zet2);
                    visc3[k][j][i].z -= (tauLocalCat.zx * zet0 + tauLocalCat.zy * zet1 + tauLocalCat.zz * zet2);
                }
            }
        }
    }

    // compute postprocess data
    MPI_Allreduce(&lTauS, &gTauS, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lVol2, &gVol2, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    gTauS = gTauS / gVol2;
    
    if (!rank) 
    {
        // Create the postProcessing directory on the first iteration.
        if (les->access->clock->it == les->access->clock->itStart) 
        {
            errno = 0;
            PetscInt dirRes = mkdir("./postProcessing", 0777);
            if (dirRes != 0 && errno != EEXIST) 
            {
                char error[512];
                sprintf(error, "could not create postProcessing directory\n");
                fatalErrorInFunction("turbStat", error);
            }
        }

        // Construct the file name.
        word fileName = "postProcessing/turbstat_Bardina_" + getStartTimeName(les->access->clock);
        FILE *fp = fopen(fileName.c_str(), "a");

        if (!fp) 
        {
            char error[512];
            sprintf(error, "cannot open file postProcessing/turbstat file\n");
            fatalErrorInFunction("turbStat", error);
        } 
        else 
        {
            PetscInt width = -15;

            // Write the current time.
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.5f\t", width, les->access->clock->time);
            // Write the volume-averaged enstrophy.
            PetscFPrintf(mesh->MESH_COMM, fp, "%*.12e\t", width, gTauS);  //Energy flux <tauijSij>
            PetscFPrintf(mesh->MESH_COMM, fp, "\n");
            

            fclose(fp);
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);
    DMDAVecRestoreArray(tda, les->lTau, &tauSGSCat);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);



    DMLocalToLocalBegin(fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalEnd  (fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalBegin(fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalEnd  (fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalBegin(fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);
    DMLocalToLocalEnd  (fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);

    return 0;
}

PetscErrorCode updateLESStructuralModelContravariantForm(les_ *les)
{
    mesh_         *mesh = les->access->mesh;
    ueqn_         *ueqn = les->access->ueqn;
    teqn_         *teqn = les->access->teqn;

    DM            da    = mesh->da, fda = mesh->fda, tda = mesh->tda;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***aj, ***Sabs, ***nut,  ***ksg, ***nvert, ***meshTag;
    PetscReal     ***iaj, ***jaj, ***kaj;
    Cmpnts        ***csi, ***eta, ***zet, ***ucat, ***ucont;
    Cmpnts        ***visc1, ***visc2, ***visc3;
    Cmpnts        ***icsi, ***jeta, ***kzet;

    PetscReal     ajc;
    PetscMPIInt   rank;

    PetscReal     dudc, dude, dudz,
                  dvdc, dvde, dvdz,
                  dwdc, dwde, dwdz;

    Tensor        ***tauSGS, tauLocal;                                      //subgrid scale stress tensor in cartesian form
    symmTensor    SS;

    PetscReal     du_dx, du_dy, du_dz, 
                  dv_dx, dv_dy, dv_dz, 
                  dw_dx, dw_dy, dw_dz;

    PetscReal     csi0, csi1, csi2, eta0, eta1, eta2, zet0, zet1, zet2;      // surface area vectors components

    lxs = xs; lxe = xe; if (xs == 0) lxs = xs+1; if (xe == mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys == 0) lys = ys+1; if (ye == my) lye = ye-1;
    lzs = zs; lze = ze; if (zs == 0) lzs = zs+1; if (ze == mz) lze = mz-1;

    MPI_Comm_rank(mesh->MESH_COMM, &rank);
    VecSet(les->lTau, 0.0);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(da,  les->lNu_t, &nut);
    DMDAVecGetArray(da,  les->lKsgs, &ksg);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da, mesh->lAj,  &aj);

    DMDAVecGetArray(tda, les->lTau, &tauSGS);

    //bardina model
    if (les->model == BAMD || les->model == BV || les->model == BDS)
    {

        for (k = lzs; k < lze; k++)
        for (j = lys; j < lye; j++)
        for (i = lxs; i < lxe; i++)
        {
            /* Skip solid cells */
            if( isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
            {
                zeroTensor(&tauSGS[k][j][i]); //inline function each element 0
                continue;
            }

            // get 1/V at the i-face
            ajc = aj[k][j][i];

            csi0 = csi[k][j][i].x, csi1 = csi[k][j][i].y, csi2 = csi[k][j][i].z;
            eta0 = eta[k][j][i].x, eta1 = eta[k][j][i].y, eta2 = eta[k][j][i].z;
            zet0 = zet[k][j][i].x, zet1 = zet[k][j][i].y, zet2 = zet[k][j][i].z;

            PetscReal filter, test_filter;                         // grid filter and dynamic model filter

            PetscReal u[3][3][3] ,   v[3][3][3] ,   w[3][3][3];    // cartesian velocity
            
            PetscReal U[3][3][3]   , V[3][3][3]   , W[3][3][3];

            PetscReal Uu[3][3][3]  , Uv[3][3][3]  , Uw[3][3][3]; //         |Uu Uv Uw|
            
            PetscReal Vu[3][3][3]  , Vv[3][3][3]  , Vw[3][3][3]; // tensor: |Vu Vv Vw|
            
            PetscReal Wu[3][3][3]  , Wv[3][3][3]  , Ww[3][3][3]; //        |wu wv ww|
            
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

                //uu tensor components
                
                Uu[R][Q][P] = U[R][Q][P] * u[R][Q][P];
                Uv[R][Q][P] = U[R][Q][P] * v[R][Q][P];
                Uw[R][Q][P] = U[R][Q][P] * w[R][Q][P];

                Vu[R][Q][P] = V[R][Q][P] * u[R][Q][P];
                Vv[R][Q][P] = V[R][Q][P] * v[R][Q][P];
                Vw[R][Q][P] = V[R][Q][P] * w[R][Q][P];

                Wu[R][Q][P] = W[R][Q][P] * u[R][Q][P];
                Wv[R][Q][P] = W[R][Q][P] * v[R][Q][P];
                Ww[R][Q][P] = W[R][Q][P] * w[R][Q][P];
            }

            //Build filter weights (using a binomial 3x3x3 kernel)
            
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

            PetscReal weight[3][3][3];

            PetscReal sum_vol = 0;

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

                    sum_vol += (1.0/aj[K][J][I]) * 8.0 * coef[R][Q][P];
                    weight[R][Q][P] = 1;
                
                } 
                
                else 
                {
                    weight[R][Q][P] = 0;
                }
            
            }
            
            //Apply test filter to local velocity fields
            
            PetscReal _u = integrateTestfilterSimpson(u, weight);
            PetscReal _v = integrateTestfilterSimpson(v, weight);
            PetscReal _w = integrateTestfilterSimpson(w, weight);
            
            PetscReal _U = integrateTestfilterSimpson(U, weight);
            PetscReal _V = integrateTestfilterSimpson(V, weight);
            PetscReal _W = integrateTestfilterSimpson(W, weight);

            //Apply test filter to product
            PetscReal _Uu = integrateTestfilterSimpson(Uu, weight);
            PetscReal _Uv = integrateTestfilterSimpson(Uv, weight);
            PetscReal _Uw = integrateTestfilterSimpson(Uw, weight);
            
            PetscReal _Vu = integrateTestfilterSimpson(Vu, weight);
            PetscReal _Vv = integrateTestfilterSimpson(Vv, weight);
            PetscReal _Vw = integrateTestfilterSimpson(Vw, weight);

            PetscReal _Wu = integrateTestfilterSimpson(Wu, weight);
            PetscReal _Wv = integrateTestfilterSimpson(Wv, weight);
            PetscReal _Ww = integrateTestfilterSimpson(Ww, weight);

            //Compute the Bardina SGS stress ---
            //tau_ij = test-filtered(u_i*u_j) - (test-filtered(u_i))*(test-filtered(u_j))
            //and store the full 3x3 tensor in tauSGS.
            
            tauSGS[k][j][i].xx = _Uu - _U*_u;   // tau_xx
            tauSGS[k][j][i].xy = _Uv - _U*_v;   // tau_xy
            tauSGS[k][j][i].xz = _Uw - _U*_w;   // tau_xz
            
            tauSGS[k][j][i].yx = _Vu - _V*_u;   // tau_yx
            tauSGS[k][j][i].yy = _Vv - _V*_v;   // tau_yy
            tauSGS[k][j][i].yz = _Vw - _V*_w;   // tau_yz
            
            tauSGS[k][j][i].zx = _Wu - _W*_u;   // tau_zx
            tauSGS[k][j][i].zy = _Wv - _W*_v;   // tau_zy
            tauSGS[k][j][i].zz = _Ww - _W*_w;   // tau_zz
        } 
    }

    DMDAVecRestoreArray(tda, les->lTau, &tauSGS);
    DMDAVecRestoreArray(da,  les->lKsgs, &ksg);
    DMDAVecRestoreArray(da,  les->lNu_t, &nut);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da, mesh->lAj,  &aj);

    DMLocalToLocalBegin(tda, les->lTau, INSERT_VALUES, les->lTau);
    DMLocalToLocalEnd  (tda, les->lTau, INSERT_VALUES, les->lTau);

    //set tau boundary conditions
    DMDAVecGetArray(tda, les->lTau, &tauSGS);

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
                    tauSGS[k][j][i].xx = tauSGS[c][b][a].xx;
                    tauSGS[k][j][i].xy = tauSGS[c][b][a].xy;
                    tauSGS[k][j][i].xz = tauSGS[c][b][a].xz;
                    tauSGS[k][j][i].yx = tauSGS[c][b][a].yx;
                    tauSGS[k][j][i].yy = tauSGS[c][b][a].yy;
                    tauSGS[k][j][i].yz = tauSGS[c][b][a].yz;
                    tauSGS[k][j][i].zx = tauSGS[c][b][a].zx;
                    tauSGS[k][j][i].zy = tauSGS[c][b][a].zy;
                    tauSGS[k][j][i].zz = tauSGS[c][b][a].zz;
                }
            }
        }
    }

    DMDAVecRestoreArray(tda, les->lTau, &tauSGS);
    DMLocalToLocalBegin(tda, les->lTau, INSERT_VALUES, les->lTau);
    DMLocalToLocalEnd  (tda, les->lTau, INSERT_VALUES, les->lTau);

    DMDAVecGetArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecGetArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecGetArray(fda, ueqn->lVisc3, &visc3);
    DMDAVecGetArray(tda, les->lTau, &tauSGS);

    DMDAVecGetArray(da, mesh->lIAj, &iaj);
    DMDAVecGetArray(da, mesh->lJAj, &jaj);
    DMDAVecGetArray(da, mesh->lKAj, &kaj);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

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

                // get 1/V at the i-face
                PetscReal ajc = iaj[k][j][i];

                // get face normals
                csi0 = icsi[k][j][i].x, csi1 = icsi[k][j][i].y, csi2 = icsi[k][j][i].z;

                tauLocal.xx = 0.5 * (tauSGS[k][j][i].xx + tauSGS[k][j][i+1].xx);
                tauLocal.xy = 0.5 * (tauSGS[k][j][i].xy + tauSGS[k][j][i+1].xy);
                tauLocal.xz = 0.5 * (tauSGS[k][j][i].xz + tauSGS[k][j][i+1].xz);


                if
                (
                    (mesh->boundaryU.iLeft=="velocityWallFunction"  && i==0) ||
                    (mesh->boundaryU.iRight=="velocityWallFunction" && i==mx-2)
                )
                {

                }
                else if
                (
                    (mesh->boundaryU.iLeft =="slip" && i==0   ) ||
                    (mesh->boundaryU.iRight=="slip" && i==mx-2)
                )
                {

                }
                else if(isIBMFluidIFace(k, j, i, i+1, nvert))
                {

                }
                else 
                {
                    visc1[k][j][i].x -= tauLocal.xx;
                    visc1[k][j][i].y -= tauLocal.xy;
                    visc1[k][j][i].z -= tauLocal.xz;
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

                // get 1/V at the j-face
                PetscReal ajc = jaj[k][j][i];

                // get face normals
                eta0 = jeta[k][j][i].x, eta1 = jeta[k][j][i].y, eta2 = jeta[k][j][i].z;

                tauLocal.yx = 0.5 * (tauSGS[k][j][i].yx + tauSGS[k][j+1][i].yx);
                tauLocal.yy = 0.5 * (tauSGS[k][j][i].yy + tauSGS[k][j+1][i].yy);
                tauLocal.yz = 0.5 * (tauSGS[k][j][i].yz + tauSGS[k][j+1][i].yz);

                if
                (
                    (mesh->boundaryU.jLeft=="velocityWallFunction" && j==0) ||
                    (mesh->boundaryU.jRight=="velocityWallFunction" && j==my-2)
                )
                {

                }
                else if
                (
                    (mesh->boundaryU.jLeft =="slip" && j==0   ) ||
                    (mesh->boundaryU.jRight=="slip" && j==my-2)
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
                    //     PetscPrintf(PETSC_COMM_WORLD,"before %.10lf %.10lf %.10lf\n", visc2[k][j][i].x, visc2[k][j][i].y, visc2[k][j][i].z);
                    //     PetscPrintf(PETSC_COMM_WORLD,"taulocal %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", tauLocal.xx, tauLocal.yy, tauLocal.zz, tauLocal.xy, tauLocal.xz, tauLocal.yz);
                    //     PetscPrintf(PETSC_COMM_WORLD,"eta %.10lf %.10lf %.10lf\n", eta0, eta1, eta2);
                    // }
                    visc2[k][j][i].x -= tauLocal.yx;
                    visc2[k][j][i].y -= tauLocal.yy;
                    visc2[k][j][i].z -= tauLocal.yz;
                    // if(k ==15 && j==2 && i==2)
                    // PetscPrintf(PETSC_COMM_WORLD,"after %.10lf %.10lf %.10lf\n", visc2[k][j][i].x, visc2[k][j][i].y, visc2[k][j][i].z);
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

                // get 1/V at the j-face
                PetscReal ajc = kaj[k][j][i];

                // get face normals
                zet0 = kzet[k][j][i].x, zet1 = kzet[k][j][i].y, zet2 = kzet[k][j][i].z;

                tauLocal.zx = 0.5 * (tauSGS[k][j][i].zx + tauSGS[k+1][j][i].zx);
                tauLocal.zy = 0.5 * (tauSGS[k][j][i].zy + tauSGS[k+1][j][i].zy);
                tauLocal.zz = 0.5 * (tauSGS[k][j][i].zz + tauSGS[k+1][j][i].zz);

                if
                (
                    (mesh->boundaryU.kLeft =="slip" && k==0   ) ||
                    (mesh->boundaryU.kRight=="slip" && k==mz-2)
                )
                {

                }
                else if(isIBMFluidKFace(k, j, i, k+1, nvert))
                {
                    
                }
                else
                {
                    visc3[k][j][i].x -= tauLocal.zx;
                    visc3[k][j][i].y -= tauLocal.zy;
                    visc3[k][j][i].z -= tauLocal.zz;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lVisc1, &visc1);
    DMDAVecRestoreArray(fda, ueqn->lVisc2, &visc2);
    DMDAVecRestoreArray(fda, ueqn->lVisc3, &visc3);
    DMDAVecRestoreArray(tda, les->lTau, &tauSGS);

    DMDAVecRestoreArray(da, mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da, mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da, mesh->lKAj, &kaj);

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    DMLocalToLocalBegin(fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalEnd  (fda, ueqn->lVisc1, INSERT_VALUES, ueqn->lVisc1);
    DMLocalToLocalBegin(fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalEnd  (fda, ueqn->lVisc2, INSERT_VALUES, ueqn->lVisc2);
    DMLocalToLocalBegin(fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);
    DMLocalToLocalEnd  (fda, ueqn->lVisc3, INSERT_VALUES, ueqn->lVisc3);

    return 0;
}
