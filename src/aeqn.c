//! \file  aeqn.c
//! \brief Contains Alpha equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

#include "solvers/aeqn_solvers.h"

//***************************************************************************************************************//

PetscErrorCode SNESMonitorA(SNES snes, PetscInt iter, PetscReal rnorm, void* comm)
{
    MPI_Comm SNES_COMM = *(MPI_Comm*)comm;
    if(iter==1)
    {
        PetscPrintf(SNES_COMM,"%e, ", rnorm);
    }
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeAEqn(aeqn_ *aeqn)
{
    if(aeqn != NULL)
    {
        // set pointer to mesh
        mesh_  *mesh  = aeqn->access->mesh;
        flags_ *flags = aeqn->access->flags;

        // input file
        PetscOptionsInsertFile(mesh->MESH_COMM, PETSC_NULL, "control.dat", PETSC_TRUE);

        VecDuplicate(mesh->Nvert, &(aeqn->AlphaTmp));    VecSet(aeqn->AlphaTmp, 0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Alpha));       VecSet(aeqn->Alpha,    0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Alpha_o));     VecSet(aeqn->Alpha_o,  0.0);
        VecDuplicate(mesh->Nvert, &(aeqn->Rhs));         VecSet(aeqn->Rhs,      0.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lAlpha));      VecSet(aeqn->lAlpha,   0.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lAlpha_o));    VecSet(aeqn->lAlpha_o, 0.0);

        VecDuplicate(mesh->lCent, &(aeqn->lDivA));       VecSet(aeqn->lDivA,    0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lDivAHO));     VecSet(aeqn->lDivAHO,  0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lDivACor));    VecSet(aeqn->lDivACor, 0.0);
        VecDuplicate(mesh->lCent, &(aeqn->lLambdaA));    VecSet(aeqn->lLambdaA, 1.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lAlphaLO));    VecSet(aeqn->lAlphaLO, 0.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lRplusA));     VecSet(aeqn->lRplusA,  1.0);
        VecDuplicate(mesh->lAj,   &(aeqn->lRminusA));    VecSet(aeqn->lRminusA, 1.0);

        // mules: always active, 3 limiter sweeps by default
        aeqn->mulesIter = 3;
        PetscOptionsGetInt(PETSC_NULL, PETSC_NULL, "-mulesIter", &(aeqn->mulesIter), PETSC_NULL);

        // read time discretization scheme
        readDictWord("control.dat", "-dAdtScheme", &(aeqn->ddtScheme));

        // create the SNES solver
        if(aeqn->ddtScheme=="BE" || aeqn->ddtScheme=="BDF2")
        {
            // default parameters
            aeqn->solverType       = "SNES";
            aeqn->absExitTol        = 1e-5;
            aeqn->relExitTol        = 1e-30;

            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-absTolA",  &(aeqn->absExitTol), PETSC_NULL);
            PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-relTolA",  &(aeqn->relExitTol), PETSC_NULL);

            SNESCreate    (mesh->MESH_COMM,&(aeqn->snesA));
            SNESMonitorSet(aeqn->snesA, SNESMonitorA, (void*)&(mesh->MESH_COMM), PETSC_NULL);

            // set the SNES evaluating function
            SNESSetFunction(aeqn->snesA, aeqn->Rhs, SNESFuncEvalA, (void *)aeqn);

            // create jacobian matrix
            MatCreateSNESMF(aeqn->snesA, &(aeqn->JA));
            SNESSetJacobian(aeqn->snesA, aeqn->JA, aeqn->JA, MatMFFDComputeJacobian, (void *)aeqn);

            // set SNES solver type
            SNESSetType(aeqn->snesA, SNESNEWTONTR);           //SNESTR
            //SNESSetType(aeqn->snesA, SNESNEWTONLS);         //SNESLS is better for stiff PDEs such as the one including IB but slower

            // set SNES solve and step failures
            SNESSetMaxLinearSolveFailures  (aeqn->snesA,10000);
            SNESSetMaxNonlinearStepFailures(aeqn->snesA,10000);
            SNESKSPSetUseEW(aeqn->snesA, PETSC_TRUE);

            // set SNES Krylov Sub-Space parameters
            SNESKSPSetParametersEW(aeqn->snesA,3,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

            // SNES tolerances
            // 2nd arg: absolute tolerance
            // 3rd arg: relative tolerance
            // 4th arg: convergene tolerance in terms of the norm of the change in the solution |deltaU| / |U| < tol
            // 5th arg: maximum number of iterations
            // 6th arg: maximum function evaluations
            SNESSetTolerances(aeqn->snesA, aeqn->absExitTol, aeqn->relExitTol, 1e-30, 20, 1000);

            SNESGetKSP(aeqn->snesA, &(aeqn->ksp));
            KSPGetPC  (aeqn->ksp,&(aeqn->pc));

            // set KSP solver type (default to GMRES)
            KSPSetType(aeqn->ksp, KSPGMRES);

            //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);    //2009.09.22 poor performance
            //KSPSetInitialGuessKnoll(ksp, PETSC_TRUE);      //2009.09.22

            //KSPFischerGuess itg;
            //KSPFischerGuessCreate(ksp,1,100,&itg);
            //KSPSetFischerGuess(ksp, itg);                  //2009.09.22

            //KSPGMRESSetPreAllocateVectors(ksp);            --> crazy thing consumes memory

            // maybe one day add a preconditioner (at least Jacobi)
            PCSetType(aeqn->pc, PCNONE);
            PetscReal rtol=aeqn->relExitTol, atol=aeqn->absExitTol, dtol=PETSC_DEFAULT;
            KSPSetTolerances(aeqn->ksp, rtol, atol, dtol, 1000);

            // previous solution for BDF2
            if(aeqn->ddtScheme == "BDF2")
            {
                VecDuplicate(mesh->Nvert, &(aeqn->Alpha_oo));  VecSet(aeqn->Alpha_oo, 0.0);
            }

            // print info
            PetscPrintf(mesh->MESH_COMM, "selected %s time scheme for Alpha-Water equation\n", aeqn->ddtScheme.c_str());
            PetscPrintf(mesh->MESH_COMM, " > relTolA = %e\n", aeqn->relExitTol);
            PetscPrintf(mesh->MESH_COMM, " > absTolA = %e\n", aeqn->absExitTol);
        }
        else if (aeqn->ddtScheme=="RK4")
        {
            aeqn->solverType       = "EXP";
        }
        else
        {
            char error[512];
            sprintf(error,  "unknown ddtScheme %s for Alpha-Water equation, available schemes are\n"
                            "    - RK4 (Runge Kutta 4)\n"
                            "    - BE (Backward Euler - BDF1)\n"
                            "    - BDF2 (2nd order Backward Differentiation Formula)\n",
                            aeqn->ddtScheme.c_str());
            fatalErrorInFunction("InitializeAEqn", error);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolveAEqn(aeqn_ *aeqn)
{
    mesh_          *mesh  = aeqn->access->mesh;
    clock_         *clock = aeqn->access->clock;

    if(aeqn->ddtScheme=="BE" || aeqn->ddtScheme=="BDF2")
        AeqnSNES(aeqn);
    else if (aeqn->ddtScheme=="RK4")
        AeqnRK4(aeqn);

    resetCellPeriodicFluxes(mesh, aeqn->Alpha, aeqn->lAlpha, "scalar", "globalToLocal");

    // compute min max of alpha water for logging
    PetscReal alphaMin, alphaMax;
    VecMin(aeqn->Alpha, NULL, &alphaMin);
    VecMax(aeqn->Alpha, NULL, &alphaMax);
    PetscPrintf(mesh->MESH_COMM, "Alpha-Water min = %e, max = %e\n", alphaMin, alphaMax);

    return(0);
}

//***************************************************************************************************************//
