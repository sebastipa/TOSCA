//! \file  ueqn.c
//! \brief Contains U equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"

#include "sources/ueqn_sources.h"
#include "solvers/ueqn_solvers.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorU(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==0)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeUEqn(ueqn_ *ueqn)
{
    // set pointer to mesh
    mesh_  *mesh  = ueqn->access->mesh;
    flags_ *flags = ueqn->access->flags;

    VecDuplicate(mesh->Cent, &(ueqn->Utmp));      VecSet(ueqn->Utmp,    0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Rhs));       VecSet(ueqn->Rhs,     0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Rhs_o));     VecSet(ueqn->Rhs_o,   0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucont));     VecSet(ueqn->Ucont,   0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucont_o));   VecSet(ueqn->Ucont_o, 0.0);
    VecDuplicate(mesh->Cent, &(ueqn->Ucat));      VecSet(ueqn->Ucat,    0.0);
    VecDuplicate(mesh->Cent, &(ueqn->dP));        VecSet(ueqn->dP,      0.0);
    VecDuplicate(mesh->Cent, &(ueqn->gCont));     VecSet(ueqn->gCont,   0.0);

    VecDuplicate(mesh->lCent, &(ueqn->lUcont));   VecSet(ueqn->lUcont,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lUcat));    VecSet(ueqn->lUcat,   0.0);

    VecDuplicate(mesh->lCent, &(ueqn->lFp));      VecSet(ueqn->lFp,     0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv1));    VecSet(ueqn->lDiv1,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv2));    VecSet(ueqn->lDiv2,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lDiv3));    VecSet(ueqn->lDiv3,   0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc1));   VecSet(ueqn->lVisc1,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc2));   VecSet(ueqn->lVisc2,  0.0);
    VecDuplicate(mesh->lCent, &(ueqn->lVisc3));   VecSet(ueqn->lVisc3,  0.0);
	VecDuplicate(mesh->lCent, &(ueqn->sourceU));  VecSet(ueqn->sourceU, 0.0);

    VecDuplicate(mesh->lAj, &(ueqn->lUstar));     VecSet(ueqn->lUstar,  0.0);

    if(ueqn->access->flags->isTeqnActive)
    {
        VecDuplicate(mesh->Cent, &(ueqn->bTheta));   VecSet(ueqn->bTheta,   0.0);
        VecDuplicate(mesh->Cent, &(ueqn->bTheta_o)); VecSet(ueqn->bTheta_o, 0.0);
    }

    if(flags->isIBMActive)
    {
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM1));   VecSet(ueqn->lViscIBM1,  0.0);
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM2));   VecSet(ueqn->lViscIBM2,  0.0);
        VecDuplicate(mesh->lCent, &(ueqn->lViscIBM3));   VecSet(ueqn->lViscIBM3,  0.0);
    }

	if(ueqn->access->flags->isGravityWaveModelingActive)
    {
		VecDuplicate(mesh->Cent, &(ueqn->dPAGW));        VecSet(ueqn->dPAGW,      0.0);
	}

    // read additional input parameters 
    PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

    // default parameters
    ueqn->inviscid          = 0;
    PetscOptionsGetInt(PETSC_NULL, PETSC_NULL,  "-inviscid", &(ueqn->inviscid),   PETSC_NULL);

    // read divergence scheme
    readDictWord("control.dat", "-divScheme", &(ueqn->divScheme));

    ueqn->centralDiv        = 0;
    ueqn->central4Div       = 0;
    ueqn->centralUpwindDiv  = 0;
    ueqn->centralUpwindWDiv = 0;
    ueqn->weno3Div          = 0;
    ueqn->quickDiv          = 0;

    if      (ueqn->divScheme == "centralUpwind")   ueqn->centralUpwindDiv  = 1;
    else if (ueqn->divScheme == "centralUpwindW")  ueqn->centralUpwindWDiv = 1;
    else if (ueqn->divScheme == "central")         ueqn->centralDiv        = 1;
    else if (ueqn->divScheme == "central4")        ueqn->central4Div       = 1;
    else if (ueqn->divScheme == "weno3")           ueqn->weno3Div          = 1;
    else if (ueqn->divScheme == "quickDiv")        ueqn->quickDiv          = 1;
    else
    {
        char error[512];
        sprintf(error,
            "unknown divScheme %s for U equation, available schemes are\n"
            "    - centralUpwind\n"
            "    - centralUpwindW\n"
            "    - central\n"
            "    - central4\n"
            "    - weno3\n"
            "    - quickDiv",
            ueqn->divScheme.c_str());
        fatalErrorInFunction("InitializeUEqn", error);
    }

    // read hyperviscosity parameter if applicable
    if(ueqn->central4Div)
    {
        ueqn->hyperVisc = 1.0;
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL,  "-hyperVisc", &(ueqn->hyperVisc),   PETSC_NULL);
        PetscPrintf(PETSC_COMM_WORLD, " > hyperviscosity parameter = %lf\n", ueqn->hyperVisc);
    }

    PetscPrintf(mesh->MESH_COMM, "selected %s divergence scheme for U equation\n", ueqn->divScheme.c_str());

    // read time discretization scheme
    readDictWord("control.dat", "-dUdtScheme", &(ueqn->ddtScheme));

    if(ueqn->ddtScheme == "CN")
    {
        ueqn->solverType   = "SNES";
        ueqn->relExitTol   = 1e-30;
        ueqn->absExitTol   = 1e-5;
        ueqn->snesMaxIter  = 20;
        ueqn->gmresRestart = 30;

        SNESType snesType  = SNESNEWTONTR;   // SNESNEWTONLS

        // read optional parameters for SNES solver
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolU",          &(ueqn->relExitTol),   PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolU",          &(ueqn->absExitTol),   PETSC_NULL);
        PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-snesMaxItersU",    &(ueqn->snesMaxIter),  PETSC_NULL);
        PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-kspGMRESRestartU", &(ueqn->gmresRestart), PETSC_NULL);

        // KSP type: BiCGStab (default) or GMRES
        PetscBool found;
        char kspTypeBuf[256] = "BiCGStab";
        PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-kspTypeU", kspTypeBuf, sizeof(kspTypeBuf), &found);
        ueqn->kspType        = kspTypeBuf;

        if(ueqn->kspType != "BiCGStab" && ueqn->kspType != "GMRES")
        {
            char error[512];
            sprintf(error,
                "unknown kspTypeU %s for U equation, available types are\n"
                "    BiCGStab (default)\n"
                "    GMRES",
                ueqn->kspType.c_str());
            fatalErrorInFunction("InitializeUEqn", error);
        }

        // create the SNES solver
        SNESCreate(mesh->MESH_COMM, &(ueqn->snesU));
        SNESMonitorSet(ueqn->snesU, SNESMonitorU, (void*)&(mesh->MESH_COMM), PETSC_NULL);

        // set the SNES evaluating function
        SNESSetFunction(ueqn->snesU, ueqn->Rhs, SNESFuncEval, (void *)ueqn);

        // set the newton algorithm used to solve the nonlinear system (SNESNEWTONTR/SNESNEWTONLS)
        SNESSetType(ueqn->snesU, snesType);

        // create matrix-free Jacobian
        MatCreateSNESMF(ueqn->snesU, &(ueqn->JU));
        MatMFFDSetType (ueqn->JU, MATMFFD_WP);
        SNESSetJacobian(ueqn->snesU, ueqn->JU, ueqn->JU, MatMFFDComputeJacobian, (void *)ueqn);

        // set SNES solve and step failures
        SNESSetMaxLinearSolveFailures  (ueqn->snesU,10000);
        SNESSetMaxNonlinearStepFailures(ueqn->snesU,10000);

        // enable Eisenstat-Walker for adaptive KSP tolerances
        SNESKSPSetUseEW       (ueqn->snesU, PETSC_TRUE);
        SNESKSPSetParametersEW(ueqn->snesU, 3, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);

        // SNES tolerances 
        SNESSetTolerances(ueqn->snesU, ueqn->absExitTol, ueqn->relExitTol, 1e-30, ueqn->snesMaxIter, 1000);

        // get KSP and PC
        SNESGetKSP(ueqn->snesU, &(ueqn->ksp));
        KSPGetPC  (ueqn->ksp, &(ueqn->pc));

        // KSP solver
        if(ueqn->kspType == "BiCGStab")
        {
            KSPSetType(ueqn->ksp, KSPBCGS);
        }
        else if(ueqn->kspType == "GMRES")
        {
            KSPSetType(ueqn->ksp, KSPGMRES);
            KSPGMRESSetRestart(ueqn->ksp, ueqn->gmresRestart);
        }

        PCSetType(ueqn->pc, PCNONE);
        PetscReal rtol = ueqn->relExitTol, atol = ueqn->absExitTol, dtol = 1e30;
        KSPSetTolerances(ueqn->ksp, rtol, atol, dtol, 1000);

        // allocate the constant-RHS precomputation buffer
        VecDuplicate(ueqn->Ucont, &(ueqn->bU)); VecSet(ueqn->bU, 0.0);
        
        // print information
        PetscPrintf(mesh->MESH_COMM, "selected %s time scheme for U equation\n", ueqn->ddtScheme.c_str());
        PetscPrintf(mesh->MESH_COMM, " > relTolU = %e\n", ueqn->relExitTol);
        PetscPrintf(mesh->MESH_COMM, " > absTolU = %e\n", ueqn->absExitTol);
        PetscPrintf(mesh->MESH_COMM, " > snesMaxIterU = %d\n", ueqn->snesMaxIter);
        PetscPrintf(mesh->MESH_COMM, " > kspTypeU = %s\n", ueqn->kspType.c_str());
        if(ueqn->kspType == "GMRES")
        {
            PetscPrintf(mesh->MESH_COMM, " > kspGMRESRestartU = %d\n", ueqn->gmresRestart);
        }
    }
    else if
    (
        ueqn->ddtScheme == "ABCN" || 
        ueqn->ddtScheme == "RK3SOCN" || 
        ueqn->ddtScheme == "RK3WCN"
    )
    {
        // KSP type: BiCGStab or GMRES (default)
        //           BiCGStab did not show to work properly
        
        PetscBool found;
        char kspTypeBuf[256] = "GMRES";
        PetscOptionsGetString(PETSC_NULL, PETSC_NULL, "-kspTypeU", kspTypeBuf, sizeof(kspTypeBuf), &found);
        ueqn->kspType      = kspTypeBuf;
        ueqn->solverType   = ueqn->kspType;
        ueqn->snesMaxIter  = 5;
        ueqn->gmresRestart = 30;

        PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-kspGMRESRestartU", &(ueqn->gmresRestart), PETSC_NULL);

        if(ueqn->kspType != "BiCGStab" && ueqn->kspType != "GMRES")
        {
            char error[512];
            sprintf(error,
                "unknown kspTypeU %s for U equation, available types are\n"
                "    BiCGStab\n"
                "    GMRES (default)",
                ueqn->kspType.c_str());
            fatalErrorInFunction("InitializeUEqn", error);
        }

        // read KSP exit tolerances from control.dat
        ueqn->relExitTol = 1e-30;
        ueqn->absExitTol = 1e-5;
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolU", &(ueqn->relExitTol), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolU", &(ueqn->absExitTol), PETSC_NULL);

        // Create a standalone KSP for the linear system A*U^{n+1} = b
        PetscInt n, N;
        VecGetLocalSize(ueqn->Ucont, &n);
        VecGetSize(ueqn->Ucont, &N);
        MatCreateShell(mesh->MESH_COMM, n, n, N, N, (void*)ueqn, &(ueqn->JvIMEX));
        MatShellSetOperation(ueqn->JvIMEX, MATOP_MULT, (void(*)(void))IMEXMatVec);
        //MatShellSetOperation(ueqn->JvIMEX, MATOP_GET_DIAGONAL, (void(*)(void))IMEXMatGetDiagonal);
        KSPCreate(mesh->MESH_COMM, &(ueqn->kspIMEX));
        KSPSetOperators(ueqn->kspIMEX, ueqn->JvIMEX, ueqn->JvIMEX);

        // one day: add a preconditioner (at least Jacobi)
        PC pcIMEX;
        KSPGetPC(ueqn->kspIMEX, &pcIMEX);
        PCSetType(pcIMEX, PCNONE); 

        if(ueqn->kspType == "BiCGStab")
        {
            KSPSetType(ueqn->kspIMEX, KSPBCGSL);
            KSPBCGSLSetEll(ueqn->kspIMEX, 2);
        }
        else if(ueqn->kspType == "GMRES")
        {
            KSPSetType(ueqn->kspIMEX, KSPGMRES);
            KSPGMRESSetRestart(ueqn->kspIMEX, ueqn->gmresRestart);
        }

        // activate initial guess nonzero
        KSPSetInitialGuessNonzero(ueqn->kspIMEX, PETSC_TRUE);
        KSPSetTolerances(ueqn->kspIMEX, ueqn->relExitTol, ueqn->absExitTol, 1e30, 1000);

        // allocate the constant-RHS buffer 
        VecDuplicate(ueqn->Ucont, &(ueqn->bU));        VecSet(ueqn->bU, 0.0);
        VecDuplicate(ueqn->Ucont, &(ueqn->RhsVisc));   VecSet(ueqn->RhsVisc,   0.0);

        // allocate Adams-Bashforth 2 convective RHS history buffers
        if(ueqn->ddtScheme == "ABCN")
        {
            VecDuplicate(ueqn->Ucont, &(ueqn->RhsConv));   VecSet(ueqn->RhsConv,   0.0);
            VecDuplicate(ueqn->Ucont, &(ueqn->RhsConv_o)); VecSet(ueqn->RhsConv_o, 0.0);
        }
        
        // print info 
        PetscPrintf(mesh->MESH_COMM, "selected %s time scheme for U equation\n", ueqn->ddtScheme.c_str());
        PetscPrintf(mesh->MESH_COMM, " > relTolU  = %e\n", ueqn->relExitTol);
        PetscPrintf(mesh->MESH_COMM, " > absTolU  = %e\n", ueqn->absExitTol);
        PetscPrintf(mesh->MESH_COMM, " > kspTypeU = %s\n", ueqn->kspType.c_str());
        if(ueqn->kspType == "GMRES")
        {
            PetscPrintf(mesh->MESH_COMM, " > kspGMRESRestartU = %d\n", ueqn->gmresRestart);
        }
    }
    else if(ueqn->ddtScheme == "FE" || ueqn->ddtScheme == "RK4")
    {
        // print info 
        ueqn->solverType   = "EXP";
        PetscPrintf(mesh->MESH_COMM, "selected %s time scheme for U equation\n", ueqn->ddtScheme.c_str());
    }
    else
    {
        char error[512];
        sprintf(error, "unknown ddtScheme %s for U equation, available schemes are\n"
                       "    - CN (Crank Nicolson)\n"
                       "    - FE (Forward Euler)\n"
                       "    - RK4 (Runge Kutta 4)\n"
                       "    - ABCN (Crank Nicolson - Adam Bashforth IMEX)\n"
                       "    - RK3SOCN (Sor-Osher Runge Kutta 3 (TVD) - Crank Nicolson IMEX)\n"
                       "    - RK3WCN (Wray-Lund Runge Kutta 3 - Crank Nicolson IMEX)\n",
                       ueqn->ddtScheme.c_str());
        fatalErrorInFunction("InitializeUEqn", error);
    }

    // explicit  6th order hyperviscosity.
    ueqn->hyperVisc4i = 0.0; ueqn->hyperVisc4j = 0.0; ueqn->hyperVisc4k = 0.0;
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-hyperViscU_i", &(ueqn->hyperVisc4i), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-hyperViscU_j", &(ueqn->hyperVisc4j), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-hyperViscU_k", &(ueqn->hyperVisc4k), PETSC_NULL);
    if
    (
        ueqn->hyperVisc4i < 0.0 || ueqn->hyperVisc4j < 0.0 || ueqn->hyperVisc4k < 0.0 ||
        ueqn->hyperVisc4i >= 1.0 || ueqn->hyperVisc4j >= 1.0 || ueqn->hyperVisc4k >= 1.0
    )
    {
        char error[512];
        sprintf(error, "-hyperViscU_{i,j,k} must be > 0 and < 1\n");
        fatalErrorInFunction("InitializeUEqn", error);
    }
    if(ueqn->hyperVisc4i > 0.0 || ueqn->hyperVisc4j > 0.0 || ueqn->hyperVisc4k > 0.0)
    {
        PetscPrintf(mesh->MESH_COMM, "explicit hyperviscosity factor for U equation:\n");
        if(ueqn->hyperVisc4i > 0.0) { PetscPrintf(mesh->MESH_COMM, " > hyperViscU_i = %f\n", ueqn->hyperVisc4i); }
        if(ueqn->hyperVisc4j > 0.0) { PetscPrintf(mesh->MESH_COMM, " > hyperViscU_j = %f\n", ueqn->hyperVisc4j); }
        if(ueqn->hyperVisc4k > 0.0) { PetscPrintf(mesh->MESH_COMM, " > hyperViscU_k = %f\n", ueqn->hyperVisc4k); }
    }

    // bulk velocity controller 
    if(ueqn->access->flags->isBulkGradPForcingActive)
    {
        ueqn->uBulk = nSetZero();
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-u_bulk", &(ueqn->meanGradP.x), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-v_bulk", &(ueqn->meanGradP.y), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-w_bulk", &(ueqn->meanGradP.z), PETSC_NULL);

        // at least one component of the bulk velocity controller must be nonzero if the controller is active
        if(ueqn->uBulk.x == 0.0 && ueqn->uBulk.y == 0.0 && ueqn->uBulk.z == 0.0)
        {
            char error[512];
            sprintf(error, "at least one component of the bulk velocity controller for U equation must be nonzero if the controller is active\n");
            fatalErrorInFunction("InitializeUEqn", error);
        }

        PetscPrintf(mesh->MESH_COMM, "bulk velocity controller for U equation:\n");
        if(ueqn->uBulk.x != 0.0) { PetscPrintf(mesh->MESH_COMM, " > u_bulk = %f\n", ueqn->uBulk.x); }
        if(ueqn->uBulk.y != 0.0) { PetscPrintf(mesh->MESH_COMM, " > v_bulk = %f\n", ueqn->uBulk.y); }
        if(ueqn->uBulk.z != 0.0) { PetscPrintf(mesh->MESH_COMM, " > w_bulk = %f\n", ueqn->uBulk.z); }
    }

    // mean constant pressure gradient 
    if(ueqn->access->flags->isMeanGradPForcingActive)
    {
        ueqn->meanGradP = nSetZero();
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-dpdx_mean", &(ueqn->meanGradP.x), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-dpdy_mean", &(ueqn->meanGradP.y), PETSC_NULL);
        PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-dpdz_mean", &(ueqn->meanGradP.z), PETSC_NULL);

        // at least one component of the mean pressure gradient must be nonzero if the controller is active
        if(ueqn->meanGradP.x == 0.0 && ueqn->meanGradP.y == 0.0 && ueqn->meanGradP.z == 0.0)
        {
            char error[512];
            sprintf(error, "at least one component of the mean pressure gradient forcing for U equation must be nonzero if the forcing is active\n");
            fatalErrorInFunction("InitializeUEqn", error);
        }   

        PetscPrintf(mesh->MESH_COMM, "mean pressure gradient forcing for U equation:\n");
        if(ueqn->meanGradP.x != 0.0) { PetscPrintf(mesh->MESH_COMM, " > dpdx_mean = %f\n", ueqn->meanGradP.x); }
        if(ueqn->meanGradP.y != 0.0) { PetscPrintf(mesh->MESH_COMM, " > dpdy_mean = %f\n", ueqn->meanGradP.y); }
        if(ueqn->meanGradP.z != 0.0) { PetscPrintf(mesh->MESH_COMM, " > dpdz_mean = %f\n", ueqn->meanGradP.z); }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveUEqn(ueqn_ *ueqn)
{
    if      (ueqn->ddtScheme == "CN")
        UeqnSNES(ueqn);
    else if (ueqn->ddtScheme == "ABCN")
        UeqnABCN(ueqn);
    else if (ueqn->ddtScheme == "RK3WCN")
        UeqnRK3CN_W(ueqn);
    else if (ueqn->ddtScheme == "RK3SOCN")
        UeqnRK3CN_SO(ueqn);
    else if (ueqn->ddtScheme == "FE")
        UeqnFE(ueqn);
    else if (ueqn->ddtScheme == "RK4")
        UeqnRK4(ueqn);
    
    // reset no penetration fluxes to zero (override numerical errors)
    resetNoPenetrationFluxes(ueqn);

    // reset periodic fluxes to be consistent if the flow is periodic
    resetFacePeriodicFluxesVector(ueqn->access->mesh, ueqn->Ucont, ueqn->lUcont, "globalToLocal");

    // apply IBM boundary condition to reset fluxes at IBM interface
    if(ueqn->access->flags->isIBMActive)
    {
        UpdateImmersedBCs(ueqn->access->ibm);
    }    

    // updates BCs around blanked cells
    if(ueqn->access->flags->isOversetActive)
    {
        setBackgroundBC(ueqn->access->mesh);
    }

    // adjust inflow/outflow fluxes to ensure mass conservation
    adjustFluxesLocal(ueqn);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateFluxLimiter(ueqn_ *ueqn)
{
    mesh_         *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucont, ***limiter;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, mesh->fluxLimiter, &limiter);

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if(j!=0 && k!=0)
                {
                    PetscInt    iL, iR;
                    getFace2Face3StencilCsi(mesh, i, mx, &iL, &iR);

                    limiter[k][j][i].x
                    =
                    vanLeer
                    (
                        ucont[k][j][iL].x,
                        ucont[k][j][i].x,
                        ucont[k][j][iR].x
                    );
                }

                if(i!=0 && k!=0)
                {
                    PetscInt    jL, jR;
                    getFace2Face3StencilEta(mesh, j, my, &jL, &jR);

                    limiter[k][j][i].y
                    =
                    vanLeer
                    (
                        ucont[k][jL][i].y,
                        ucont[k][j][i].y,
                        ucont[k][jR][i].y
                    );
                }

                if(j!=0 && i!=0)
                {
                    PetscInt    kL, kR;
                    getFace2Face3StencilZet(mesh, k, mz, &kL, &kR);

                    limiter[k][j][i].z
                    =
                    vanLeer
                    (
                        ucont[kL][j][i].z,
                        ucont[k][j][i].z,
                        ucont[kR][j][i].z
                    );
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, mesh->fluxLimiter, &limiter);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode contravariantToCartesian(ueqn_ *ueqn)
{
    mesh_           *mesh = ueqn->access->mesh;
    DM               da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    PetscReal        mat[3][3], det, det0, det1, det2;

    PetscReal        ***aj, ***nvert, ***meshTag;
    Cmpnts           ***ucat, ***lucont;

    PetscReal q[3];  // local working array

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    Cmpnts    ***csi,  ***eta,  ***zet;

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);
    DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);

    // interpolate ucont from faces to cell centers
    // transform from contravariant to cartesian

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if ( isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    mat[0][0] = (csi[k][j][i].x);
                    mat[0][1] = (csi[k][j][i].y);
                    mat[0][2] = (csi[k][j][i].z);

                    mat[1][0] = (eta[k][j][i].x);
                    mat[1][1] = (eta[k][j][i].y);
                    mat[1][2] = (eta[k][j][i].z);

                    mat[2][0] = (zet[k][j][i].x);
                    mat[2][1] = (zet[k][j][i].y);
                    mat[2][2] = (zet[k][j][i].z);


                    PetscInt iL = i-1,  iR = i;
                    PetscInt jL = j-1,  jR = j;
                    PetscInt kL = k-1,  kR = k;

                    // interpolate ucont from opposite faces to cell centers
                    q[0] = central( lucont[k][j][iL].x, lucont[k][j][iR].x);
                    q[1] = central( lucont[k][jL][i].y, lucont[k][jR][i].y);
                    q[2] = central( lucont[kL][j][i].z, lucont[kR][j][i].z);

                    // transform to cartesian
                    det =   mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

                    det0 =  q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                            q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                    det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                            q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                    det2 =  q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                            q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                            q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                    ucat[k][j][i].x = det0 / det;
                    ucat[k][j][i].y = det1 / det;
                    ucat[k][j][i].z = det2 / det;
                }
            }
        }
    }

    DMDAVecRestoreArray (fda, ueqn->Ucat,  &ucat);
    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd  (fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMDAVecRestoreArray (fda, ueqn->lUcont, &lucont);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode contravariantToCartesianGeneric(mesh_ *mesh, Vec &lCont, Vec &lCat)
{
    DM               da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    PetscReal        mat[3][3], det, det0, det1, det2;

    PetscReal        ***aj, ***nvert, ***meshTag;
    Cmpnts           ***lcat, ***lcont;
    Cmpnts           ***csi,  ***eta,  ***zet;

    PetscReal q[3];  // local working array

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, lCont, &lcont);
    DMDAVecGetArray(fda, lCat,  &lcat);

    // interpolate ucont from faces to cell centers
    // transform from contravariant to cartesian

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if ( isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    mat[0][0] = (csi[k][j][i].x);
                    mat[0][1] = (csi[k][j][i].y);
                    mat[0][2] = (csi[k][j][i].z);

                    mat[1][0] = (eta[k][j][i].x);
                    mat[1][1] = (eta[k][j][i].y);
                    mat[1][2] = (eta[k][j][i].z);

                    mat[2][0] = (zet[k][j][i].x);
                    mat[2][1] = (zet[k][j][i].y);
                    mat[2][2] = (zet[k][j][i].z);


                    PetscInt iL = i-1,  iR = i;
                    PetscInt jL = j-1,  jR = j;
                    PetscInt kL = k-1,  kR = k;

                    // interpolate ucont from opposite faces to cell centers
                    q[0] = central( lcont[k][j][iL].x, lcont[k][j][iR].x);
                    q[1] = central( lcont[k][jL][i].y, lcont[k][jR][i].y);
                    q[2] = central( lcont[kL][j][i].z, lcont[kR][j][i].z);

                    // transform to cartesian
                    det =   mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

                    det0 =  q[0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
                            q[1] * (mat[0][1] * mat[2][2] - mat[0][2] * mat[2][1]) +
                            q[2] * (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]);

                    det1 = -q[0] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                            q[1] * (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) -
                            q[2] * (mat[0][0] * mat[1][2] - mat[0][2] * mat[1][0]);

                    det2 =  q[0] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) -
                            q[1] * (mat[0][0] * mat[2][1] - mat[0][1] * mat[2][0]) +
                            q[2] * (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]);

                    lcat[k][j][i].x = det0 / det;
                    lcat[k][j][i].y = det1 / det;
                    lcat[k][j][i].z = det2 / det;
                }
            }
        }
    }

    DMDAVecRestoreArray (fda, lCat,  &lcat);
    DMLocalToLocalBegin(fda, lCat, INSERT_VALUES, lCat);
    DMLocalToLocalEnd  (fda, lCat, INSERT_VALUES, lCat);
    DMDAVecRestoreArray (fda, lCont, &lcont);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode adjustFluxesLocal(ueqn_ *ueqn)
{
    mesh_           *mesh = ueqn->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont, ***lucont,
                     ***coor;
    PetscReal        epsilon = 1.e-10, ***nvert;

    PetscScalar      globalFlux;
    PetscScalar      lFluxIn = 0.0, lFluxOut = 0.0,
                     FluxIn  = 0.0, FluxOut  = 0.0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->Ucont,  &ucont);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    // instead of saying that fluxin is on left patches and fluxout is on right patches
    // we can say that fluxin is on the patches with negative velocity and fluxout is on the patches with positive velocity
    // The mass imbalance is then distributed to the cells with positive velocity based on the ratio between their area and
    // the cumulative area of the outflow faces

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // k-left boundary
        if (zs==0 && mesh->boundaryU.kLeftPatchType == 0)
        {
            // k-left boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].z > 0.0 && (!isIBMCell(k+1, j, i, nvert)))
                    {
                        lFluxIn += ucont[k][j][i].z;
                    }
                    else if(ucont[k][j][i].z < 0.0 && (!isIBMCell(k+1, j, i, nvert)))
                    {
                        lFluxOut += fabs(ucont[k][j][i].z);
                    }
                    else 
                    {
                        ucont[k][j][i].z = 0.;
                    }
                }
            }
        }

        // compute outflow flux at k-right boundary
        if (ze==mz && mesh->boundaryU.kRightPatchType == 0)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].z > 0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxOut += ucont[k][j][i].z;
                    }
                    else if(ucont[k][j][i].z < 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxIn += fabs(ucont[k][j][i].z);
                    }
                    else
                    {
                        ucont[k][j][i].z = 0.;
                    }
                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0 && mesh->boundaryU.jLeftPatchType == 0)
        {
            // j-left boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].y > 0.0 && (!isIBMCell(k, j+1, i, nvert)))
                    {
                        lFluxIn += ucont[k][j][i].y;
                    }
                    else if(ucont[k][j][i].y < 0.0 && (!isIBMCell(k, j+1, i, nvert)))
                    {
                        lFluxOut += fabs(ucont[k][j][i].y);
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }
                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my && mesh->boundaryU.jRightPatchType == 0)
        {
            // j-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].y > 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxOut += ucont[k][j][i].y;
                    }
                    else if(ucont[k][j][i].y < 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxIn += fabs(ucont[k][j][i].y);
                    }
                    else
                    {
                        ucont[k][j][i].y = 0.;
                    }
                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0 && mesh->boundaryU.iLeftPatchType == 0)
        {
            // i-left boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].x > 0.0 && (!isIBMCell(k, j, i+1, nvert)))
                    {
                        lFluxIn += ucont[k][j][i].x;
                    }
                    else if(ucont[k][j][i].x < 0.0 && (!isIBMCell(k, j, i+1, nvert)))
                    {
                        lFluxOut += fabs(ucont[k][j][i].x);
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }
                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx && mesh->boundaryU.iRightPatchType == 0)
        {
            // i-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].x > 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxOut += ucont[k][j][i].x;
                    }
                    else if(ucont[k][j][i].x < 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        lFluxIn += fabs(ucont[k][j][i].x);
                    }
                    else
                    {
                        ucont[k][j][i].x = 0.;
                    }
                }
            }
        }
    }

    // cumulate the net influx and net outflux
    MPI_Allreduce(&lFluxIn, &FluxIn, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    //PetscPrintf(mesh->MESH_COMM, "Pre correction: Fluxin = %lf, Fluxout = %lf\n", FluxIn, FluxOut);

    globalFlux = FluxOut - FluxIn;

    lFluxOut = 0.0;

    if
    (
        !mesh->k_periodic && !mesh->kk_periodic
    )
    {
        // k-left boundary
        if (zs==0 && mesh->boundaryU.kLeftPatchType == 0)
        {
            // k-left boundary face
            k = 0;

            // loop on the boundary faces
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].z < 0.0 && (!isIBMCell(k+1, j, i, nvert)))
                    {
                        ucont[k][j][i].z += globalFlux * (fabs(ucont[k][j][i].z) / FluxOut);
                        lFluxOut += fabs(ucont[k][j][i].z);   
                    }
                }
            }
        }

        // correct outflow flux at k-right boundary
        if (ze==mz && mesh->boundaryU.kRightPatchType == 0)
        {
            // k-right boundary face
            k = mz-2;

            // loop on the boundary cells
            for (j=lys; j<lye; j++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].z > 0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        ucont[k][j][i].z -= globalFlux * (ucont[k][j][i].z / FluxOut);
                        lFluxOut += ucont[k][j][i].z;
                    }
                }
            }
        }
    }

    if
    (
        !mesh->j_periodic && !mesh->jj_periodic
    )
    {
        // compute flux at j-left boundary
        if (ys==0 && mesh->boundaryU.jLeftPatchType == 0)
        {
            // j-left boundary face
            j = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].y < 0.0 && (!isIBMCell(k, j+1, i, nvert)))
                    {
                        ucont[k][j][i].y += globalFlux * (fabs(ucont[k][j][i].y) / FluxOut);
                        lFluxOut += fabs(ucont[k][j][i].y);
                    }
                }
            }
        }

        // compute flux at j-right boundary
        if (ye==my && mesh->boundaryU.jRightPatchType == 0)
        {
            // j-right boundary face
            j = my-2;

            // loop on the boundary cells
            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].y > 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        ucont[k][j][i].y -= globalFlux * (ucont[k][j][i].y / FluxOut);
                        lFluxOut += ucont[k][j][i].y;
                    }
                }
            }
        }
    }

    if
    (
        !mesh->i_periodic && !mesh->ii_periodic
    )
    {
        // compute flux at i-left boundary
        if (xs==0 && mesh->boundaryU.iLeftPatchType == 0)
        {
            // i-left boundary face
            i = 0;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].x < 0.0 && (!isIBMCell(k, j, i+1, nvert)))
                    {
                        ucont[k][j][i].x += globalFlux * (fabs(ucont[k][j][i].x) / FluxOut);
                        lFluxOut += fabs(ucont[k][j][i].x);
                    }
                }
            }
        }

        // compute flux at i-right boundary
        if (xe==mx && mesh->boundaryU.iRightPatchType == 0)	
        {
            // i-right boundary face
            i = mx-2;

            // loop on the boundary faces
            for (k=lzs; k<lze; k++)
            {
                for (j=lys; j<lye; j++)
                {
                    // cumulate flux
                    if(ucont[k][j][i].x > 0.0 && (!isIBMCell(k, j, i, nvert)))
                    {
                        ucont[k][j][i].x -= globalFlux * (ucont[k][j][i].x / FluxOut);
                        lFluxOut += ucont[k][j][i].x;
                    }
                }
            }
        }
    }

    // cumulate the net influx and net outflux
    // MPI_Allreduce(&lFluxOut, &FluxOut, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    // PetscPrintf(mesh->MESH_COMM, "Post correction: Fluxin = %lf, Fluxout = %lf\n", FluxIn, FluxOut);

    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &lucont);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

