//! \file  aeqn_solvers.c
//! \brief Alpha-Water-equation solver function implementations.

#include "aeqn_solvers.h"

PetscErrorCode FormA(aeqn_ *aeqn, Vec &Rhs, PetscReal scale, PetscInt formMode)
{
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AeqnSNES(aeqn_ *aeqn)
{
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SNESFuncEvalA(SNES snes, Vec Alpha, Vec Rhs, void *ptr)
{
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode FormExplicitRhsA(aeqn_ *aeqn)
{
    return(0);
}

//***************************************************************************************************************//

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