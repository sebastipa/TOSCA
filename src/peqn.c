//! \file  peqn.c
//! \brief Contains P equation function definitions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"

// connectivity macro
#define matID(i,j,k) (HYPRE_Int)(gid[k][j][i])

// 19 point stencil for building the coefficient matrix
#define CP 0  // center point
#define EP 1  // east point
#define WP 2  // west point
#define NP 3  // north point
#define SP 4  // south point
#define TP 5  // top point
#define BP 6  // bottom point
#define NE 7  // north east point
#define SE 8  // south east point
#define NW 9  // north west point
#define SW 10 // south west point
#define TN 11 // top west point
#define BN 12 // bottom north point
#define TS 13 // top south point
#define BS 14 // bottom south point
#define TE 15 // top east point
#define BE 16 // bottom east point
#define TW 17 // bottom east point
#define BW 18 // bottom west point

//***************************************************************************************************************//

PetscErrorCode InitializePEqn(peqn_ *peqn)
{
    // set pointer to mesh
    mesh_ *mesh = peqn->access->mesh;

    VecDuplicate(mesh->Nvert,  &(peqn->P));         VecSet(peqn->P,   0.0);
    VecDuplicate(mesh->Nvert,  &(peqn->Phi));       VecSet(peqn->Phi, 0.0);

    VecDuplicate(mesh->lNvert, &(peqn->lP));        VecSet(peqn->lP,  0.0);
    VecDuplicate(mesh->lNvert, &(peqn->lPhi));      VecSet(peqn->lPhi,0.0);
    VecDuplicate(mesh->lNvert, &(peqn->lGid));
    VecDuplicate(mesh->lNvert, &(peqn->lLid));

    // default parameters
    peqn->hypreSolverType = 1;
    peqn->amgAgg          = 0;
    peqn->amgThresh       = 0.5;
    peqn->amgCoarsenType  = 10;
    peqn->poissonIt       = 8;
    peqn->poissonTol      = 1.e-8;

    PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-amgAgg",         &(peqn->amgAgg),   PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-amgThresh",      &(peqn->amgThresh), PETSC_NULL);
    PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-amgCoarsenType", &(peqn->amgCoarsenType), PETSC_NULL);
    PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-poissonIt",      &(peqn->poissonIt), PETSC_NULL);
    PetscOptionsGetReal(PETSC_NULL, PETSC_NULL, "-poissonTol",     &(peqn->poissonTol), PETSC_NULL);
    PetscOptionsGetInt (PETSC_NULL, PETSC_NULL, "-hypreSolverType",&(peqn->hypreSolverType),   PETSC_NULL);
    readDictWord       ("control.dat",          "-poissonSolver",  &(peqn->solverType));

    // compute connectivity and create the solver
    SetPoissonConnectivity(peqn);

    if(peqn->solverType == "HYPRE")
    {
        // create parallel coeff. matrix framework
        CreateHypreMatrix(peqn);

        // create parallel sol and rhs vector frameworks
        CreateHypreVector(peqn);

        // create solver framework
        CreateHypreSolver(peqn);
    }
    else if(peqn->solverType == "PETSc")
    {
        // create parallel rhs vector frameworks
        CreatePETScVector(peqn);

        // create parallel coeff. matrix framework
        CreatePETScMatrix(peqn);

    }
    else
    {
        char error[512];
        sprintf(error, "unknown solverType type %s. Known types are HYPRE and PETSc\n", peqn->solverType.c_str());
        fatalErrorInFunction("InitializePEqn",  error);
    }

    SetCoeffMatrix(peqn);

    if(peqn->solverType == "PETSc")
    {
        // create solver framework
        CreatePETScSolver(peqn);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetPoissonConnectivity(peqn_ *peqn)
{
    // Creates connectivity for HYPRE Poisson matrix solver. Writes Gid (Global IDs)
    // and defines the starting label of this processor (p_global_begin), the last
    // label of this processor (p_global_begin + local_Phi2_size - 1) and the number
    // of cells w/o the ghosts cells, namely reduced_p_size.
    // Also creates the parallel vector Phi2 according to the set-up connectivity.

    mesh_         *mesh  = peqn->access->mesh;
    clock_        *clock = peqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k, r;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***nvert, ***meshTag;
    PetscReal     ***gid,         // global cell ID
                  ***lid;         // local cell ID

    PetscInt      nUntilThisRank; // total number of dofs from proc. 0 to this proc. (excluded)

    PetscMPIInt   rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);
    PetscMPIInt   size; MPI_Comm_size(mesh->MESH_COMM, &size);

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // The vector of local IDs will then be transformed to global id by adding
    // the specific processor start label.
    // When building the coefficent matrix of the Poisson iteration the correct
    // cells will be accessed by looping over cells and indexing the gid vector
    // with normal indexing. Gid will return the label of that cell, which is a
    // value going from 0 to ndofs-1.
    // A test will be performed in order to skip IBM cells since gid will be -1
    // there.

    VecSet(peqn->lGid, -1);
    VecSet(peqn->lLid, -1);

    DMDAVecGetArray(da, peqn->lGid,    &gid);
    DMDAVecGetArray(da, peqn->lLid,    &lid);
    DMDAVecGetArray(da, mesh->lNvert,  &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag,  &meshTag);

    // number of pressure dof for each processor
    std::vector<PetscInt> lndof(size), gndof(size);

    for (r=0; r<size; r++)
    {
        lndof[r] = 0;
    }

    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                // accumulate the number of dofs for current rank
                // and set lid to the label of the current node
                lid[k][j][i] = (PetscReal)lndof[rank]++;
            }
        }
    }

    // scatter the number of dofs for each process into ndof_node
    MPI_Allreduce(&lndof[0], &gndof[0], size, MPIU_INT, MPI_MAX, mesh->MESH_COMM);

    // get total accumulated number of dofs belonging to process 0 to this process (excluded)
    // note: the starting label of this processor will be startID, and the last
    //       label will be (startID + nIDs - 1)
    nUntilThisRank = 0;
    for (r=0; r<rank; r++)
    {
        nUntilThisRank = nUntilThisRank + (PetscInt)gndof[r];
    }

    // set local start idx and size into the Peqn context
    peqn->thisRankStart  = (HYPRE_Int)nUntilThisRank;
    peqn->thisRankSize   = (HYPRE_Int)gndof[rank];
    peqn->thisRankEnd    = peqn->thisRankStart + peqn->thisRankSize - 1;

    if (rank == size-1)
    {
        // set the total number of dofs
        peqn->totalSize = (HYPRE_Int)nUntilThisRank + (HYPRE_Int)gndof[rank];
    }

    // broadcast to other processors
    MPI_Bcast(&(peqn->totalSize), 1, MPI_INT, size-1, mesh->MESH_COMM);

    MPI_Barrier(mesh->MESH_COMM);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                // test: only enter if we are not on boundaries
                if((PetscInt)lid[k][j][i] >= 0)
                {
                    // set gid to the current label + the start label on this proc.
                    gid[k][j][i] = (PetscReal)(lid[k][j][i] + nUntilThisRank);
                }
            }
        }
    }

    // destroy Phi2 of previous iteration
    if (clock->it > clock->itStart)
    {
        VecDestroy(&peqn->phi);
    }

    // restore the arrays
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag,  &meshTag);
    DMDAVecRestoreArray(da, peqn->lGid, &gid);
    DMDAVecRestoreArray(da, peqn->lLid, &lid);

    DMLocalToLocalBegin(da, peqn->lGid, INSERT_VALUES, peqn->lGid);
    DMLocalToLocalEnd  (da, peqn->lGid, INSERT_VALUES, peqn->lGid);

    // reset cartesian periodic cells if the flow is periodic
    resetCellPeriodicFluxes(mesh, peqn->lGid, peqn->lGid, "scalar", "localToLocal");

    // create Phi2 of current iteration
    VecCreateMPI(mesh->MESH_COMM, gndof[rank], PETSC_DETERMINE, &peqn->phi);

    // free memeory
    std::vector<PetscInt> ().swap(lndof);
    std::vector<PetscInt> ().swap(gndof);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CreateHypreMatrix(peqn_ *peqn)
{
    // set number of non-zero entries per row
    HYPRE_Int *nonZeroEntriesPerRow;
    PetscMalloc(sizeof(HYPRE_Int)*peqn->thisRankSize, &nonZeroEntriesPerRow);

    for(PetscInt i=0; i<peqn->thisRankSize; i++)
    {
        nonZeroEntriesPerRow[i] = 19;
    }

    HYPRE_IJMatrixCreate
    (
            peqn->access->mesh->MESH_COMM, // communicator
            peqn->thisRankStart,          // i lower
            peqn->thisRankEnd,            // i upper
            peqn->thisRankStart,          // j lower
            peqn->thisRankEnd,            // j upper
            &(peqn->hypreA)
    );

    // set the matrix storage type (only HYPRE_PARCSR available at this time)
    HYPRE_IJMatrixSetObjectType(peqn->hypreA, HYPRE_PARCSR);

    // set the expected number of non-zero in each row (for speedup)
    HYPRE_IJMatrixSetRowSizes(peqn->hypreA, nonZeroEntriesPerRow);

    // set the expected number of addition from current to other procs (for speedup).
    // In HYPRE the matrix is split on the rows and distributed (so that columns are
    // shared by different processors). Since insertions happen on the same row
    // no intra-processor additions/insertions are performed
    HYPRE_IJMatrixSetMaxOffProcElmts (peqn->hypreA, 0);

    // initialize the hypre matrix
    HYPRE_IJMatrixInitialize(peqn->hypreA);

    HYPRE_IJMatrixGetObject(peqn->hypreA, (void**)&(peqn->hypreParA));

    PetscFree(nonZeroEntriesPerRow);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CreatePETScMatrix(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    flags_        *flags = peqn->access->flags;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      localSize;

    MatCreate(mesh->MESH_COMM, &(peqn->petscA));
    VecGetLocalSize(peqn->phi, &localSize);
    MatSetSizes(peqn->petscA, localSize, localSize, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(peqn->petscA,MATMPIAIJ);
    MatMPIAIJSetPreallocation(peqn->petscA, 19, PETSC_NULL, 10, PETSC_NULL);
    MatSetFromOptions(peqn->petscA);

    // this option tells PETSc that if sparsity patterns changes its ok (for dynamic IBM cases, does internal reallocation stuff)
    MatSetOption(peqn->petscA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    // this option tells PETSc to force constant sparsity pattern discarding new entries
    //MatSetOption(peqn->petscA,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CreateHypreVector(peqn_ *peqn)
{
    // p vector
    HYPRE_IJVectorCreate(peqn->access->mesh->MESH_COMM, peqn->thisRankStart, peqn->thisRankEnd, &(peqn->hypreP));
    HYPRE_IJVectorSetObjectType(peqn->hypreP,   HYPRE_PARCSR);

    HYPRE_IJVectorCreate(peqn->access->mesh->MESH_COMM, peqn->thisRankStart, peqn->thisRankEnd, &(peqn->hypreRhs));
    HYPRE_IJVectorSetObjectType(peqn->hypreRhs, HYPRE_PARCSR);

    HYPRE_IJVectorInitialize(peqn->hypreP  );
    HYPRE_IJVectorInitialize(peqn->hypreRhs);

    HYPRE_IJVectorGetObject(peqn->hypreP,   (void **) &(peqn->hypreParP  ));
    HYPRE_IJVectorGetObject(peqn->hypreRhs, (void **) &(peqn->hypreParRhs));

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CreatePETScVector(peqn_ *peqn)
{
    VecDuplicate(peqn->phi, &peqn->petscRhs);
    VecSet(peqn->petscRhs,0.0);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode CreateHypreSolver(peqn_ *peqn)
{
    // COARSENING METHODS FROM HYPRE LIBRARY
    // -------------------------------------
    // 0  CLJP-coarsening (a parallel coarsening algorithm using independent sets.
    // 1  classical Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
    // 3  classical Ruge-Stueben coarsening on each processor, followed by a third pass, which adds coarse
    //    points on the boundaries
    // 6  Falgout coarsening (uses 1 first, followed by CLJP using the interior coarse points
    //    generated by 1 as its first independent set)
    // 7  CLJP-coarsening (using a fixed random vector, for debugging purposes only)
    // 8  PMIS-coarsening (a parallel coarsening algorithm using independent sets, generating
    //    lower complexities than CLJP, might also lead to slower convergence)
    // 9  PMIS-coarsening (using a fixed random vector, for debugging purposes only)
    // 10 HMIS-coarsening (uses one pass Ruge-Stueben on each processor independently, followed
    //    by PMIS using the interior C-points generated as its first independent set)
    // 11 one-pass Ruge-Stueben coarsening on each processor, no boundary treatment (not recommended!)
    // 21 CGC coarsening by M. Griebel, B. Metsch and A. Schweitzer
    // 22 CGC-E coarsening by M. Griebel, B. Metsch and A.Schweitzer

    // PRECONDITIONER
    HYPRE_BoomerAMGCreate(&(peqn->hyprePC));
    HYPRE_BoomerAMGSetInterpType(peqn->hyprePC, 13);                    // 6:ext, 13: FF1

    if(peqn->amgAgg)
    {
        HYPRE_BoomerAMGSetAggNumLevels(peqn->hyprePC, peqn->amgAgg);    // FF1 + aggresive coarsening is good for > 50mil grids
    }

    HYPRE_BoomerAMGSetStrongThreshold(peqn->hyprePC, peqn->amgThresh);  // 0.5 : Cartesian, 0.6 : Distorted

    HYPRE_BoomerAMGSetTol (peqn->hyprePC, peqn->poissonTol);
    HYPRE_BoomerAMGSetPrintLevel(peqn->hyprePC, 0);                     // 1
    HYPRE_BoomerAMGSetCoarsenType(peqn->hyprePC, peqn->amgCoarsenType); // 0:CLJP, 6:Falgout, 8:PMIS, 10:HMIS,
    HYPRE_BoomerAMGSetCycleType(peqn->hyprePC, 1);                      // 1
    HYPRE_BoomerAMGSetMaxIter(peqn->hyprePC, 3);                        // haji 28 oct 2017 set it to 2 (it was 1 before)

    //HYPRE_BoomerAMGSetRelaxType(precon_p, 6);
    //HYPRE_BoomerAMGSetRelaxWt(precon_p,0.8);
    //HYPRE_BoomerAMGSetCycleNumSweeps(precon_p,5,3);         // 5 sweep at coarsest level
    //HYPRE_BoomerAMGSetLevelOuterWt(precon_p,0.5,0);
    HYPRE_BoomerAMGSetSmoothType(peqn->hyprePC, 7);                 // more complex smoother
    HYPRE_BoomerAMGSetSmoothNumSweeps(peqn->hyprePC,3);

    // SOLVER
    if(peqn->hypreSolverType == 1)
    {
        HYPRE_ParCSRGMRESCreate (peqn->access->mesh->MESH_COMM, &(peqn->hypreSlvr));
        HYPRE_ParCSRGMRESSetKDim(peqn->hypreSlvr, 5);
        HYPRE_GMRESSetPrecond   (peqn->hypreSlvr, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, peqn->hyprePC);
        HYPRE_GMRESSetMaxIter   (peqn->hypreSlvr, peqn->poissonIt);
        HYPRE_GMRESSetTol       (peqn->hypreSlvr, peqn->poissonTol);
        HYPRE_GMRESSetPrintLevel(peqn->hypreSlvr, 0);
    }
    else if(peqn->hypreSolverType == 2)
    {
        HYPRE_ParCSRPCGCreate (peqn->access->mesh->MESH_COMM, &(peqn->hypreSlvr));
        HYPRE_PCGSetTwoNorm   (peqn->hypreSlvr, 1);
        HYPRE_PCGSetLogging   (peqn->hypreSlvr, 1);
        HYPRE_PCGSetPrecond   (peqn->hypreSlvr, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, peqn->hyprePC);
        HYPRE_PCGSetMaxIter   (peqn->hypreSlvr, peqn->poissonIt);
        HYPRE_PCGSetTol       (peqn->hypreSlvr, peqn->poissonTol);
        HYPRE_PCGSetPrintLevel(peqn->hypreSlvr, 0 );
    }
    else
    {
       char error[512];
        sprintf(error, "Unknown HYPRE solver type. Available types are:\n    - 1: GMRES\n    - 2: PCG\n");
        fatalErrorInFunction("CreateHypreSolver",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode MyKSPMonitorPoisson(KSP ksp, PetscInt iter, PetscReal rnorm, void* dummy)
{
    peqn_* peqn = (peqn_*)dummy;
    if(iter==1)
    {
        peqn->initialPoissonRes = rnorm;
    }
    return(0);
}

PetscErrorCode CreatePETScSolver(peqn_ *peqn)
{
    mesh_ *mesh = peqn->access->mesh;

    KSPCreate(mesh->MESH_COMM, &(peqn->ksp));

    //PetscViewerAndFormat *vf;
	//PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD,PETSC_VIEWER_DEFAULT,&vf);
	//KSPMonitorSet(peqn->ksp,(PetscErrorCode (*)(KSP,PetscInt,PetscReal,void*))KSPMonitorTrueResidualNorm,vf,(PetscErrorCode (*)(void**))PetscViewerAndFormatDestroy);
    KSPMonitorSet(peqn->ksp, MyKSPMonitorPoisson, (void*)peqn, PETSC_NULL);

    KSPSetOperators(peqn->ksp, peqn->petscA, peqn->petscA);
    // hardcoded option
    //PetscOptionsSetValue(PETSC_NULL, "-ksp_initial_guess_nonzero", "true");

    KSPGetPC(peqn->ksp, &(peqn->petscPC));

    //PCSetType(peqn->petscPC,PCJACOBI);
    //PCSetType(peqn->petscPC, PCMG);
    //PCSetType(peqn->petscPC, PCILU);
    //PCMGSetType(peqn->petscPC,PC_MG_KASKADE);
	PCSetType(peqn->petscPC,PCSOR);
	//PCSetType(peqn->petscPC,PCGAMG);

    // omega = 1 equivalent to Gauss Seidel preconditioner
    PCSORSetOmega(peqn->petscPC,1.0);
    PCSORSetIterations(peqn->petscPC,10,2);

    KSPSetPCSide(peqn->ksp,PC_LEFT);

    // null space
    MatNullSpaceCreate(mesh->MESH_COMM, PETSC_TRUE,0,0, &(peqn->petscNs));
    MatSetNullSpace(peqn->petscA,peqn->petscNs);


    KSPSetType(peqn->ksp,KSPPGMRES);

    KSPSetTolerances(peqn->ksp, 1e-30, peqn->poissonTol, 1e30, peqn->poissonIt);

	KSPGMRESSetRestart(peqn->ksp, 100);
	KSPSetFromOptions(peqn->ksp);

    KSPSetUp(peqn->ksp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DestroyHypreMatrix(peqn_ *peqn)
{
    HYPRE_IJMatrixDestroy(peqn->hypreA);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DestroyHypreVector(peqn_ *peqn)
{
    HYPRE_IJVectorDestroy(peqn->hypreP);
    HYPRE_IJVectorDestroy(peqn->hypreRhs);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DestroyPETScVector(peqn_ *peqn)
{
    VecDestroy(&(peqn->petscRhs));
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DestroyHypreSolver(peqn_ *peqn)
{
    // destroy AMG preconditioner
    HYPRE_BoomerAMGDestroy(peqn->hyprePC);

    if(peqn->hypreSolverType == 1)
    {
        // destroy GMRES solver if created
        HYPRE_ParCSRGMRESDestroy(peqn->hypreSlvr);
    }
    else if(peqn->hypreSolverType == 2)
    {
        // destroy PCG solver if created
        HYPRE_ParCSRPCGDestroy  (peqn->hypreSlvr);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode DestroyPETScSolver(peqn_ *peqn)
{
    // destroy AMG preconditioner
    KSPDestroy(&(peqn->ksp));

    return(0);
}

//***************************************************************************************************************//

cellIds GetIdFromStencil(int stencil, int k, int j, int i)
{
    cellIds SId;
    // center point
    if(stencil == CP)
    {
        SId.k = k;
        SId.j = j;
        SId.i = i;
    }
    // east point
    else if (stencil == EP)
    {
        SId.k = k;
        SId.j = j;
        SId.i = i+1;
    }
    // west point
    else if (stencil == WP)
    {
        SId.k = k;
        SId.j = j;
        SId.i = i-1;
    }
    // north point
    else if (stencil == NP)
    {
        SId.k = k;
        SId.j = j+1;
        SId.i = i;
    }
    // south point
    else if (stencil == SP)
    {
        SId.k = k;
        SId.j = j-1;
        SId.i = i;
    }
    // top point
    else if (stencil == TP)
    {
        SId.k = k+1;
        SId.j = j;
        SId.i = i;
    }
    // south point
    else if (stencil == BP)
    {
        SId.k = k-1;
        SId.j = j;
        SId.i = i;
    }
    // north east
    else if (stencil == NE)
    {
        SId.k = k;
        SId.j = j+1;
        SId.i = i+1;
    }
    // south east
    else if (stencil == SE)
    {
        SId.k = k;
        SId.j = j-1;
        SId.i = i+1;
    }
    // north west
    else if (stencil == NW)
    {
        SId.k = k;
        SId.j = j+1;
        SId.i = i-1;
    }
    // south west
    else if (stencil == SW)
    {
        SId.k = k;
        SId.j = j-1;
        SId.i = i-1;
    }
    // top north
    else if (stencil == TN)
    {
        SId.k = k+1;
        SId.j = j+1;
        SId.i = i;
    }
    // bottom north
    else if (stencil == BN)
    {
        SId.k = k-1;
        SId.j = j+1;
        SId.i = i;
    }
    // top south
    else if (stencil == TS)
    {
        SId.k = k+1;
        SId.j = j-1;
        SId.i = i;
    }
    // bottom south
    else if (stencil == BS)
    {
        SId.k = k-1;
        SId.j = j-1;
        SId.i = i;
    }
    // top east
    else if (stencil == TE)
    {
        SId.k = k+1;
        SId.j = j;
        SId.i = i+1;
    }
    // bottom east
    else if (stencil == BE)
    {
        SId.k = k-1;
        SId.j = j;
        SId.i = i+1;
    }
    // top west
    else if (stencil == TW)
    {
        SId.k = k+1;
        SId.j = j;
        SId.i = i-1;
    }
    // bottom west
    else if (stencil == BW)
    {
        SId.k = k-1;
        SId.j = j;
        SId.i = i-1;
    }
    else
    {
        fatalErrorInFunction("GetIdFromStencil",  "wrong stencil");
    }
    return SId;
}

//***************************************************************************************************************//

PetscErrorCode SetCoeffMatrix(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    flags_        *flags = peqn->access->flags;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt      gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt      gys  = info.gys, gye = info.gys + info.gym;
    PetscInt      gzs  = info.gzs, gze = info.gzs + info.gzm;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;
    Cmpnts        ***cent;

    PetscReal     ***aj, ***iaj, ***jaj, ***kaj;
    PetscReal     ***nvert, ***gid, ***meshTag;

    Vec           G11, G12, G13, G21, G22, G23, G31, G32, G33;
    PetscReal     ***g11, ***g12, ***g13,
                  ***g21, ***g22, ***g23,
                  ***g31, ***g32, ***g33;

    PetscScalar   vol[27];
    HYPRE_Int     row;
    PetscInt      idx[27], m;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    if(peqn->solverType == "PETSc")
    {
        MatZeroEntries(peqn->petscA);
    }

    DMDAVecGetArray(fda, mesh->lCsi,  &csi);
    DMDAVecGetArray(fda, mesh->lEta,  &eta);
    DMDAVecGetArray(fda, mesh->lZet,  &zet);
    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);
    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);
    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);
    DMDAVecGetArray(da,  mesh->lAj,   &aj);
    DMDAVecGetArray(da,  mesh->lIAj,  &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,  &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,  &kaj);

    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, mesh->lCent, &cent);

    DMDAVecGetArray(da, peqn->lGid, &gid);

    // create metric tensor
    VecDuplicate(mesh->lAj, &G11);
    VecDuplicate(mesh->lAj, &G12);
    VecDuplicate(mesh->lAj, &G13);
    VecDuplicate(mesh->lAj, &G21);
    VecDuplicate(mesh->lAj, &G22);
    VecDuplicate(mesh->lAj, &G23);
    VecDuplicate(mesh->lAj, &G31);
    VecDuplicate(mesh->lAj, &G32);
    VecDuplicate(mesh->lAj, &G33);

    DMDAVecGetArray(da, G11, &g11);
    DMDAVecGetArray(da, G12, &g12);
    DMDAVecGetArray(da, G13, &g13);
    DMDAVecGetArray(da, G21, &g21);
    DMDAVecGetArray(da, G22, &g22);
    DMDAVecGetArray(da, G23, &g23);
    DMDAVecGetArray(da, G31, &g31);
    DMDAVecGetArray(da, G32, &g32);
    DMDAVecGetArray(da, G33, &g33);

    for (k=lzs-1; k<lze+1; k++)
    {
        for (j=lys-1; j<lye+1; j++)
        {
            for (i=lxs-1; i<lxe+1; i++)
            {
                PetscInt a=i, b=j, c=k;

                if(i==0)
                {
                    if(mesh->i_periodic)       a=mx-2;
                    else if(mesh->ii_periodic) a=-2;
                    else                      a=1;
                }
                if(i==mx-1)
                {
                    if(mesh->i_periodic)       a=1;
                    else if(mesh->ii_periodic) a=mx+1;
                    else                      a=mx-2;
                }
                if(j==0)
                {
                    if(mesh->j_periodic)       b=my-2;
                    else if(mesh->jj_periodic) b=-2;
                    else                      b=1;
                }
                if(j==my-1)
                {
                    if(mesh->j_periodic)       b=1;
                    else if(mesh->jj_periodic) b=my+1;
                    else                      b=my-2;
                }
                if(k==0)
                {
                    if(mesh->k_periodic)       c=mz-2;
                    else if(mesh->kk_periodic) c=-2;
                    else                      c=1;
                }
                if(k==mz-1)
                {
                    if(mesh->k_periodic)       c=1;
                    else if(mesh->kk_periodic) c=mz+1;
                    else                      c=mz-2;
                }

                // building metrics tensor
                g11[k][j][i] = (icsi[c][b][a].x * icsi[c][b][a].x + icsi[c][b][a].y * icsi[c][b][a].y + icsi[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
                g12[k][j][i] = (ieta[c][b][a].x * icsi[c][b][a].x + ieta[c][b][a].y * icsi[c][b][a].y + ieta[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
                g13[k][j][i] = (izet[c][b][a].x * icsi[c][b][a].x + izet[c][b][a].y * icsi[c][b][a].y + izet[c][b][a].z * icsi[c][b][a].z) * iaj[c][b][a];
                g21[k][j][i] = (jcsi[c][b][a].x * jeta[c][b][a].x + jcsi[c][b][a].y * jeta[c][b][a].y + jcsi[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
                g22[k][j][i] = (jeta[c][b][a].x * jeta[c][b][a].x + jeta[c][b][a].y * jeta[c][b][a].y + jeta[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
                g23[k][j][i] = (jzet[c][b][a].x * jeta[c][b][a].x + jzet[c][b][a].y * jeta[c][b][a].y + jzet[c][b][a].z * jeta[c][b][a].z) * jaj[c][b][a];
                g31[k][j][i] = (kcsi[c][b][a].x * kzet[c][b][a].x + kcsi[c][b][a].y * kzet[c][b][a].y + kcsi[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
                g32[k][j][i] = (keta[c][b][a].x * kzet[c][b][a].x + keta[c][b][a].y * kzet[c][b][a].y + keta[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
                g33[k][j][i] = (kzet[c][b][a].x * kzet[c][b][a].x + kzet[c][b][a].y * kzet[c][b][a].y + kzet[c][b][a].z * kzet[c][b][a].z) * kaj[c][b][a];
            }
        }
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal one  = 1.0;
                PetscReal zero = 0.0;

                row = matID(i, j, k);

                if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
                {
                    if(peqn->solverType == "HYPRE")
                    {
                        HYPRE_Int nrows = 1, ncols = 1;

                        HYPRE_Int col = row;
                        HYPRE_IJMatrixSetValues(peqn->hypreA, nrows, &ncols, &row, &col, &one);
                    }
                    else if(peqn->solverType == "PETSc")
                    {
                        PetscInt rowId = (PetscInt)row;
                        MatSetValues(peqn->petscA, 1, &rowId, 1, &rowId, &one, ADD_VALUES);
                    }
                }
                // i,j,k is a fluid point
                else
                {
                    // initialize contributions to zero
                    for (m=0; m<19; m++)
                    {
                        vol[m] = 0.;
                    }

                    PetscReal r = 1.0;

                    // contribution from east face in i-direction (i+1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (i != mx-2       ||
                        mesh->i_periodic ||
                        mesh->ii_periodic) && (isFluidCell(k,j,i+1,nvert) && isCalculatedCell(k,j,i+1,meshTag))

                    )
                    {
                        // dpdc{i} = (p_{i+1} - p_{i}) * g11_{i}
                        vol[CP] -= g11[k][j][i] / r; // i, j, k
                        vol[EP] += g11[k][j][i] / r; // i+1, j, k

                        // dpde{i} = ({p_{i+1,j+1} + p{i, j+1} - p{i+1, j-1} - p{i, j-1}) * 0.25 * g12[k][j][i]

                        // j-right boundary -> use upwind
                        if( (j == my-2) || isIBMIFace(k, j+1, i, i+1, nvert) || isOversetIFace(k, j+1, i, i+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g12[k][j][i] * 0.5 / r; // i, j, k
                            vol[EP] += g12[k][j][i] * 0.5 / r; // i+1, j, k
                            vol[SP] -= g12[k][j][i] * 0.5 / r; // i, j-1, k
                            vol[SE] -= g12[k][j][i] * 0.5 / r; // i+1, j-1, k
                        }
                        // j-left boundary -> use upwind
                        else if((j == 1) || isIBMIFace(k, j-1, i, i+1, nvert) || isOversetIFace(k, j-1, i, i+1, meshTag))
                        {
                            // upwind differencing
                            vol[NP] += g12[k][j][i] * 0.5 / r; // i, j+1, k
                            vol[NE] += g12[k][j][i] * 0.5 / r; // i+1, j+1, k
                            vol[CP] -= g12[k][j][i] * 0.5 / r; // i, j, k
                            vol[EP] -= g12[k][j][i] * 0.5 / r; // i+1, j, k
                        }
                        else
                        {
                            // central differencing
                            vol[NP] += g12[k][j][i] * 0.25 / r; // i, j+1, k
                            vol[NE] += g12[k][j][i] * 0.25 / r; // i+1, j+1, k
                            vol[SP] -= g12[k][j][i] * 0.25 / r; // i, j-1, k
                            vol[SE] -= g12[k][j][i] * 0.25 / r; // i+1, j-1, k
                        }

                        // dpdz{i}=(p_{i,k+1} + p{i+1,k+1} - p{i, k-1} - p{i+1, k-1}) * 0.25 / r * g13[k][j][i]

                        // k-right boundary -> use upwind
                        if((k==mz-2) || isIBMIFace(k+1, j, i, i+1, nvert) || isOversetIFace(k+1, j, i, i+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g13[k][j][i] * 0.5 / r; // i, j, k
                            vol[EP] += g13[k][j][i] * 0.5 / r; // i+1, j, k
                            vol[BP] -= g13[k][j][i] * 0.5 / r; // i, j, k-1
                            vol[BE] -= g13[k][j][i] * 0.5 / r; // i+1, j, k-1
                        }
                        // k-left boundary  -> use upwind
                        else if((k==1) || isIBMIFace(k-1, j, i, i+1, nvert) || isOversetIFace(k-1, j, i, i+1, meshTag))
                        {
                            // upwind differencing
                            vol[TP] += g13[k][j][i] * 0.5 / r; // i, j, k+1
                            vol[TE] += g13[k][j][i] * 0.5 / r; // i+1, j, k+1
                            vol[CP] -= g13[k][j][i] * 0.5 / r; // i, j, k
                            vol[EP] -= g13[k][j][i] * 0.5 / r; // i+1, j, k
                        }
                        else
                        {
                            // central differencing
                            vol[TP] += g13[k][j][i] * 0.25 / r; //i, j, k+1
                            vol[TE] += g13[k][j][i] * 0.25 / r; //i+1, j, k+1
                            vol[BP] -= g13[k][j][i] * 0.25 / r; //i, j, k-1
                            vol[BE] -= g13[k][j][i] * 0.25 / r; //i+1, j, k-1
                        }
                    }

                    // contribution from west face in i-direction (i-1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (i != 1 ||
                        mesh->i_periodic ||
                        mesh->ii_periodic) && isFluidCell(k, j, i-1, nvert ) && isCalculatedCell(k, j, i-1, meshTag)
                    )
                    {
                        // -dpdc{i-1} = -(p_{i} - p_{i-1}) * g11_{i}
                        vol[CP] -= g11[k][j][i-1] / r;  //i, j, k
                        vol[WP] += g11[k][j][i-1] / r;  //i-1, j, k

                        // -dpde{i-1} = -({p_{i,j+1}+p{i-1, j+1} - p{i, j-1}-p{i-1, j-1}) * 0.25 / r * g12[k][j][i-1]

                        // j-right boundary -> use upwind
                        if((j==my-2) || isIBMIFace(k, j+1, i, i-1, nvert) || isOversetIFace(k, j+1, i, i-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g12[k][j][i-1] * 0.5 / r; // i, j, k
                            vol[WP] -= g12[k][j][i-1] * 0.5 / r; // i-1, j, k
                            vol[SP] += g12[k][j][i-1] * 0.5 / r; // i, j-1, k
                            vol[SW] += g12[k][j][i-1] * 0.5 / r; // i-1, j-1, k
                        }
                        // j-left boundary -> use upwind
                        else if((j==1) || isIBMIFace(k, j-1, i, i-1, nvert) || isOversetIFace(k, j-1, i, i-1, meshTag))
                        {
                            // upwind differencing
                            vol[NP] -= g12[k][j][i-1] * 0.5 / r; // i, j+1, k
                            vol[NW] -= g12[k][j][i-1] * 0.5 / r; // i-1, j+1, k
                            vol[CP] += g12[k][j][i-1] * 0.5 / r; // i, j, k
                            vol[WP] += g12[k][j][i-1] * 0.5 / r; // i-1, j, k
                        }
                        else
                        {
                            // central differencing
                            vol[NP] -= g12[k][j][i-1] * 0.25 / r; // i, j+1, k
                            vol[NW] -= g12[k][j][i-1] * 0.25 / r; // i-1, j+1, k
                            vol[SP] += g12[k][j][i-1] * 0.25 / r; // i, j-1, k
                            vol[SW] += g12[k][j][i-1] * 0.25 / r; // i-1, j-1, k
                        }

                        // -dpdz{i-1}=-(p_{i,k+1}+p{i-1,k+1} - p{i, k-1}-p{i-1, k-1}) * 0.25 / r * g13[k][j][i]

                        // k-right boundary -> use upwind
                        if((k==mz-2) || isIBMIFace(k+1, j, i, i-1, nvert) || isOversetIFace(k+1, j, i, i-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g13[k][j][i-1] * 0.5 / r; // i, j, k
                            vol[WP] -= g13[k][j][i-1] * 0.5 / r; // i-1, j, k
                            vol[BP] += g13[k][j][i-1] * 0.5 / r; // i, j, k-1
                            vol[BW] += g13[k][j][i-1] * 0.5 / r; // i-1, j, k-1
                        }
                        // k-left boundary  -> use upwind
                        else if((k==1) || isIBMIFace(k-1, j, i, i-1, nvert) || isOversetIFace(k-1, j, i, i-1, meshTag))
                        {
                            // upwind differencing
                            vol[TP] -= g13[k][j][i-1] * 0.5 / r; // i, j, k+1
                            vol[TW] -= g13[k][j][i-1] * 0.5 / r; // i-1, j, k+1
                            vol[CP] += g13[k][j][i-1] * 0.5 / r; // i, j, k
                            vol[WP] += g13[k][j][i-1] * 0.5 / r; // i-1, j, k
                        }
                        else
                        {
                            // central differencing
                            vol[TP] -= g13[k][j][i-1] * 0.25 / r; // i, j, k+1
                            vol[TW] -= g13[k][j][i-1] * 0.25 / r; // i-1, j, k+1
                            vol[BP] += g13[k][j][i-1] * 0.25 / r; // i, j, k-1
                            vol[BW] += g13[k][j][i-1] * 0.25 / r; // i-1, j, k-1
                        }
                    }

                    // contribution from north face in j-direction (j+1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (j != my-2  ||
                        mesh->j_periodic ||
                        mesh->jj_periodic) && isFluidCell(k, j+1, i, nvert) && isCalculatedCell(k, j+1, i, meshTag)
                    )
                    {
                        // dpdc{j} = (p_{i+1,j}+p{i+1,j+1} - p{i-1,j}-p{i-1,j+1}) * 0.25

                        // i-right boundary -> use upwind
                        if((i==mx-2) || isIBMJFace(k, j, i+1, j+1, nvert) || isOversetJFace(k, j, i+1, j+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g21[k][j][i] * 0.5 / r; // i, j, k
                            vol[NP] += g21[k][j][i] * 0.5 / r; // i, j+1, k
                            vol[WP] -= g21[k][j][i] * 0.5 / r; // i-1, j, k
                            vol[NW] -= g21[k][j][i] * 0.5 / r; // i-1, j+1, k
                        }
                        // i-left boundary -> use upwind
                        else if((i==1) || isIBMJFace(k, j, i-1, j+1, nvert) || isOversetJFace(k, j, i-1, j+1, meshTag))
                        {
                            // upwind differencing
                            vol[EP] += g21[k][j][i] * 0.5 / r; // i+1, j, k
                            vol[NE] += g21[k][j][i] * 0.5 / r; // i+1, j+1, k
                            vol[CP] -= g21[k][j][i] * 0.5 / r; // i, j, k
                            vol[NP] -= g21[k][j][i] * 0.5 / r; // i, j+1, k
                        }
                        else
                        {
                            // central differencing
                            vol[EP] += g21[k][j][i] * 0.25 / r; // i+1, j, k
                            vol[NE] += g21[k][j][i] * 0.25 / r; // i+1, j+1, k
                            vol[WP] -= g21[k][j][i] * 0.25 / r; // i-1, j, k
                            vol[NW] -= g21[k][j][i] * 0.25 / r; // i-1, j+1, k
                        }

                        // dpde{j} = (p{j+1} - p{j}) * g22[k][j][i]
                        vol[CP] -= g22[k][j][i] / r;
                        vol[NP] += g22[k][j][i] / r;

                        // dpdz{j} = (p{j, k+1}+p{j+1,k+1} - p{j,k-1}-p{j+1,k-1}) *0.25

                        // k-right boundary -> use upwind
                        if((k==mz-2) || isIBMJFace(k+1, j, i, j+1, nvert) || isOversetJFace(k+1, j, i, j+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g23[k][j][i] * 0.5 / r; //i,j,k
                            vol[NP] += g23[k][j][i] * 0.5 / r; //i, j+1, k
                            vol[BP] -= g23[k][j][i] * 0.5 / r;//i, j, k-1
                            vol[BN] -= g23[k][j][i] * 0.5 / r;//i, j+1, k-1
                        }
                        // k-left boundary -> use upwind
                        else if((k==1) || isIBMJFace(k-1, j, i, j+1, nvert) || isOversetJFace(k-1, j, i, j+1, meshTag))
                        {
                            // upwind differencing
                            vol[TP] += g23[k][j][i] * 0.5 / r; //i, j, k+1
                            vol[TN] += g23[k][j][i] * 0.5 / r;//i, j+1, k+1
                            vol[CP] -= g23[k][j][i] * 0.5 / r;//i, j, k
                            vol[NP] -= g23[k][j][i] * 0.5 / r;//i, j+1, k
                        }
                        else
                        {
                            // central differencing
                            vol[TP] += g23[k][j][i] * 0.25 / r; // i, j, k+1
                            vol[TN] += g23[k][j][i] * 0.25 / r; // i, j+1, k+1
                            vol[BP] -= g23[k][j][i] * 0.25 / r; // i, j, k-1
                            vol[BN] -= g23[k][j][i] * 0.25 / r; // i, j+1, k-1
                        }
                    }

                    // contribution from south face in j-direction (j-1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (j!=1       ||
                        mesh->j_periodic ||
                        mesh->jj_periodic) && isFluidCell(k, j-1, i, nvert) && isCalculatedCell(k, j-1, i, meshTag)
                    )
                    {
                        // -dpdc{j-1} = -(p_{i+1,j}+p{i+1,j-1} - p{i-1,j}-p{i-1,j-1}) * 0.25

                        // i-right boundary -> use upwind
                        if((i==mx-2) || isIBMJFace(k, j, i+1, j-1, nvert) || isOversetJFace(k, j, i+1, j-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g21[k][j-1][i] * 0.5 / r; // i, j, k
                            vol[SP] -= g21[k][j-1][i] * 0.5 / r; // i, j-1, k
                            vol[WP] += g21[k][j-1][i] * 0.5 / r; // i-1, j, k
                            vol[SW] += g21[k][j-1][i] * 0.5 / r; // i-1, j-1, k
                        }
                        // i-left boundary -> use upwind
                        else if((i==1) || isIBMJFace(k, j, i-1, j-1, nvert) || isOversetJFace(k, j, i-1, j-1, meshTag))
                        {
                            // upwind differencing
                            vol[EP] -= g21[k][j-1][i] * 0.5 / r; // i+1, j, k
                            vol[SE] -= g21[k][j-1][i] * 0.5 / r; // i+1, j-1, k
                            vol[CP] += g21[k][j-1][i] * 0.5 / r; // i, j, k
                            vol[SP] += g21[k][j-1][i] * 0.5 / r; // i, j-1, k
                        }
                        else
                        {
                            // central differencing
                            vol[EP] -= g21[k][j-1][i] * 0.25 / r; // i+1, j, k
                            vol[SE] -= g21[k][j-1][i] * 0.25 / r; // i+1, j-1, k
                            vol[WP] += g21[k][j-1][i] * 0.25 / r; // i-1, j, k
                            vol[SW] += g21[k][j-1][i] * 0.25 / r; // i-1, j-1, k
                        }

                        // -dpde{j-1} = -(p{j} - p{j-1}) * g22[k][j-1][i]
                        vol[CP] -= g22[k][j-1][i] / r;
                        vol[SP] += g22[k][j-1][i] / r;

                        // -dpdz{j-1} = -(p{j,k+1}+p{j-1,k+1} - p{j,k-1}-p{j-1,k-1}) * 0.25

                        // k-right boundary -> use upwind
                        if((k==mz-2) || isIBMJFace(k+1, j, i, j-1, nvert) || isOversetJFace(k+1, j, i, j-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g23[k][j-1][i] * 0.5 / r; //i, j, k
                            vol[SP] -= g23[k][j-1][i] * 0.5 / r; //i, j-1, k
                            vol[BP] += g23[k][j-1][i] * 0.5 / r; //i, j, k-1
                            vol[BS] += g23[k][j-1][i] * 0.5 / r; //i, j-1, k-1
                        }
                        // k-left boundary -> use upwind
                        else if((k==1) || isIBMJFace(k-1, j, i, j-1, nvert) || isOversetJFace(k-1, j, i, j-1, meshTag))
                        {
                            // upwind differencing
                            vol[TP] -= g23[k][j-1][i] * 0.5 / r; // i, j, k+1
                            vol[TS] -= g23[k][j-1][i] * 0.5 / r; // i, j-1, k+1
                            vol[CP] += g23[k][j-1][i] * 0.5 / r; // i, j, k
                            vol[SP] += g23[k][j-1][i] * 0.5 / r; // i, j-1, k
                        }
                        else
                        {
                            // central differencing
                            vol[TP] -= g23[k][j-1][i] * 0.25 / r; // i, j, k+1
                            vol[TS] -= g23[k][j-1][i] * 0.25 / r; // i, j-1, k+1
                            vol[BP] += g23[k][j-1][i] * 0.25 / r; // i, j, k-1
                            vol[BS] += g23[k][j-1][i] * 0.25 / r; // i, j-1, k-1
                        }
                    }

                    // contribution from top face in k-direction (k+1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (k != mz-2  ||
                        mesh->k_periodic ||
                        mesh->kk_periodic) && isFluidCell(k+1, j, i, nvert) && isCalculatedCell(k+1, j, i, meshTag)
                    )
                    {
                        // dpdc{k} = (p{i+1,k}+p{i+1,k+1} - p{i-1,k}-p{i-1,k+1}) * 0.25

                        // i-right boundary -> use upwind
                        if((i==mx-2) || isIBMKFace(k, j, i+1, k+1, nvert) || isOversetKFace(k, j, i+1, k+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g31[k][j][i] * 0.5 / r; // i, j, k
                            vol[TP] += g31[k][j][i] * 0.5 / r; // i, j, k+1
                            vol[WP] -= g31[k][j][i] * 0.5 / r; // i-1, j, k
                            vol[TW] -= g31[k][j][i] * 0.5 / r; // i-1, j, k+1
                        }
                        // i-left boundary -> use upwind
                        else if((i==1) || isIBMKFace(k, j, i-1, k+1, nvert) || isOversetKFace(k, j, i-1, k+1, meshTag))
                        {
                            // upwind differencing
                            vol[EP] += g31[k][j][i] * 0.5 / r; // i+1, j, k
                            vol[TE] += g31[k][j][i] * 0.5 / r; // i+1, j, k+1
                            vol[CP] -= g31[k][j][i] * 0.5 / r; // i, j, k
                            vol[TP] -= g31[k][j][i] * 0.5 / r; // i, j, k+1
                        }
                        else
                        {
                            // central differencing
                            vol[EP] += g31[k][j][i] * 0.25 / r; // i+1, j, k
                            vol[TE] += g31[k][j][i] * 0.25 / r; // i+1, j, k+1
                            vol[WP] -= g31[k][j][i] * 0.25 / r; // i-1, j, k
                            vol[TW] -= g31[k][j][i] * 0.25 / r; // i-1, j, k+1
                        }

                        // dpde{k} = (p{j+1, k}+p{j+1,k+1} - p{j-1, k}-p{j-1,k+1}) * 0.25

                        // j-right boundary -> use upwind
                        if((j==my-2) || isIBMKFace(k, j+1, i, k+1, nvert) || isOversetKFace(k, j+1, i, k+1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] += g32[k][j][i] * 0.5 / r; // i, j,k
                            vol[TP] += g32[k][j][i] * 0.5 / r; // i, j, k+1
                            vol[SP] -= g32[k][j][i] * 0.5 / r; // i, j-1, k
                            vol[TS] -= g32[k][j][i] * 0.5 / r; // i, j-1, k+1
                        }
                        // j-left boundary -> use upwind
                        else if((j==1) || isIBMKFace(k, j-1, i, k+1, nvert) || isOversetKFace(k, j-1, i, k+1, meshTag))
                        {
                            // upwind differencing
                            vol[NP] += g32[k][j][i] * 0.5 / r; // i, j+1, k
                            vol[TN] += g32[k][j][i] * 0.5 / r; // i, j+1, k+1
                            vol[CP] -= g32[k][j][i] * 0.5 / r; // i, j, k
                            vol[TP] -= g32[k][j][i] * 0.5 / r; // i, j, k+1
                        }
                        else
                        {
                            // central differencing
                            vol[NP] += g32[k][j][i] * 0.25 / r;//i, j+1, k
                            vol[TN] += g32[k][j][i] * 0.25 / r;//i, j+1, k+1
                            vol[SP] -= g32[k][j][i] * 0.25 / r;//i, j-1, k
                            vol[TS] -= g32[k][j][i] * 0.25 / r;//i, j-1, k+1
                        }

                        // dpdz{k} = p{k+1} - p{k}
                        vol[CP] -= g33[k][j][i] / r; // i, j, k
                        vol[TP] += g33[k][j][i] / r; // i, j, k+1
                    }

                    // contribution from bottom face in k-direction (k-1/2)
                    if
                    (
                        // exclude boundary cell for zero gradient BC if non-periodic (crucial also for curvilinear)
                        (k != 1     ||
                        mesh->k_periodic ||
                        mesh->kk_periodic) &&  isFluidCell(k-1, j, i, nvert) && isCalculatedCell(k-1, j, i, meshTag)
                    )
                    {
                        // -dpdc{k-1} = -(p{i+1,k}+p{i+1,k-1} - p{i-1,k}-p{i-1,k-1}) * 0.25

                        // i-right boundary -> use upwind
                        if((i==mx-2) || isIBMKFace(k, j, i+1, k-1, nvert) || isOversetKFace(k, j, i+1, k-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g31[k-1][j][i] * 0.5 / r; // i, j, k
                            vol[BP] -= g31[k-1][j][i] * 0.5 / r; // i, j, k-1
                            vol[WP] += g31[k-1][j][i] * 0.5 / r; // i-1, j, k
                            vol[BW] += g31[k-1][j][i] * 0.5 / r; // i-1, j, k-1
                        }
                        // i-left boundary -> use upwind
                        else if((i==1) || isIBMKFace(k, j, i-1, k-1, nvert) || isOversetKFace(k, j, i-1, k-1, meshTag))
                        {
                            // upwind differencing
                            vol[EP] -= g31[k-1][j][i] * 0.5 / r; // i+1, j, k
                            vol[BE] -= g31[k-1][j][i] * 0.5 / r; // i+1, j, k-1
                            vol[CP] += g31[k-1][j][i] * 0.5 / r; // i, j, k
                            vol[BP] += g31[k-1][j][i] * 0.5 / r; // i, j, k-1
                        }
                        else
                        {
                            // central differencing
                            vol[EP] -= g31[k-1][j][i] * 0.25 / r; // i+1, j, k
                            vol[BE] -= g31[k-1][j][i] * 0.25 / r; // i+1, j, k-1
                            vol[WP] += g31[k-1][j][i] * 0.25 / r; // i-1, j, k
                            vol[BW] += g31[k-1][j][i] * 0.25 / r; // i-1, j, k-1
                        }

                        // -dpde{k-1} = -(p{j+1, k}+p{j+1,k-1} - p{j-1, k}-p{j-1,k-1}) *  0.25 / r * g32[k-1][j][i]

                        // j-right boundary -> use upwind
                        if((j==my-2) || isIBMKFace(k, j+1, i, k-1, nvert) || isOversetKFace(k, j+1, i, k-1, meshTag))
                        {
                            // upwind differencing
                            vol[CP] -= g32[k-1][j][i] * 0.5 / r; // i, j,k
                            vol[BP] -= g32[k-1][j][i] * 0.5 / r; // i, j, k-1
                            vol[SP] += g32[k-1][j][i] * 0.5 / r; // i, j-1, k
                            vol[BS] += g32[k-1][j][i] * 0.5 / r; // i, j-1, k-1
                        }
                        // j-left boundary -> use upwind
                        else if((j==1) || isIBMKFace(k, j-1, i, k-1, nvert) || isOversetKFace(k, j-1, i, k-1, meshTag))
                        {
                            // upwind differencing
                            vol[NP] -= g32[k-1][j][i] * 0.5 / r; // i, j+1, k
                            vol[BN] -= g32[k-1][j][i] * 0.5 / r; // i, j+1, k-1
                            vol[CP] += g32[k-1][j][i] * 0.5 / r; // i, j, k
                            vol[BP] += g32[k-1][j][i] * 0.5 / r; // i, j, k-1
                        }
                        else
                        {
                            // central differencing
                            vol[NP] -= g32[k-1][j][i] * 0.25 / r; // i, j+1, k
                            vol[BN] -= g32[k-1][j][i] * 0.25 / r; // i, j+1, k-1
                            vol[SP] += g32[k-1][j][i] * 0.25 / r; // i, j-1, k
                            vol[BS] += g32[k-1][j][i] * 0.25 / r; // i, j-1, k-1
                        }

                        // -dpdz{k-1} = -(p{k} - p{k-1}) * g33[k-1][j][i]
                        vol[CP] -= g33[k-1][j][i] / r; // i, j, k
                        vol[BP] += g33[k-1][j][i] / r; //i, j, k-1
                    }

                    idx[CP] = matID(i  , j  , k  );
                    idx[EP] = matID(i+1, j  , k  );
                    idx[WP] = matID(i-1, j  , k  );
                    idx[NP] = matID(i  , j+1, k  );
                    idx[SP] = matID(i  , j-1, k  );
                    idx[TP] = matID(i  , j  , k+1);
                    idx[BP] = matID(i  , j  , k-1);
                    idx[NE] = matID(i+1, j+1, k  );
                    idx[SE] = matID(i+1, j-1, k  );
                    idx[NW] = matID(i-1, j+1, k  );
                    idx[SW] = matID(i-1, j-1, k  );
                    idx[TN] = matID(i  , j+1, k+1);
                    idx[BN] = matID(i  , j+1, k-1);
                    idx[TS] = matID(i  , j-1, k+1);
                    idx[BS] = matID(i  , j-1, k-1);
                    idx[TE] = matID(i+1, j  , k+1);
                    idx[BE] = matID(i+1, j  , k-1);
                    idx[TW] = matID(i-1, j  , k+1);
                    idx[BW] = matID(i-1, j  , k-1);

                    for(m=0; m<19; m++)
                    {
                        if
                        (
                            (
                                fabs(vol[m]) > 1.e-10 &&
                                idx[m] >= 0
                            ) ||
                            m == CP
                        )
                        {
                            if(peqn->solverType == "HYPRE")
                            {
                                // the number of values to be set with this call
                                HYPRE_Int nrows = 1, ncols = 1;

                                // the global column index to be set
                                HYPRE_Int col = idx[m];

                                HYPRE_IJMatrixSetValues(peqn->hypreA, nrows, &ncols, &row, &col, &vol[m]);

                            }
                            else if(peqn->solverType == "PETSc")
                            {
                                PetscInt rowId = (PetscInt)row;
                                PetscInt col   = (PetscInt)idx[m];
                                MatSetValues(peqn->petscA, 1, &rowId, 1, &col, &vol[m], ADD_VALUES);
                            }
                        }
                    }
                }
            }
        }
    }

    // assemble the matrix
    if(peqn->solverType == "HYPRE")
    {
        HYPRE_IJMatrixAssemble(peqn->hypreA);
    }
    else if(peqn->solverType == "PETSc")
    {
        MatAssemblyBegin(peqn->petscA, MAT_FINAL_ASSEMBLY);
    	MatAssemblyEnd  (peqn->petscA, MAT_FINAL_ASSEMBLY);
    }

    // restore the arrays
    DMDAVecRestoreArray(da, G11, &g11);
    DMDAVecRestoreArray(da, G12, &g12);
    DMDAVecRestoreArray(da, G13, &g13);
    DMDAVecRestoreArray(da, G21, &g21);
    DMDAVecRestoreArray(da, G22, &g22);
    DMDAVecRestoreArray(da, G23, &g23);
    DMDAVecRestoreArray(da, G31, &g31);
    DMDAVecRestoreArray(da, G32, &g32);
    DMDAVecRestoreArray(da, G33, &g33);

    VecDestroy(&G11);
    VecDestroy(&G12);
    VecDestroy(&G13);
    VecDestroy(&G21);
    VecDestroy(&G22);
    VecDestroy(&G23);
    VecDestroy(&G31);
    VecDestroy(&G32);
    VecDestroy(&G33);

    DMDAVecRestoreArray(fda, mesh->lCsi,  &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,  &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,  &zet);
    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);
    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);
    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);
    DMDAVecRestoreArray(da,  mesh->lAj,   &aj);
    DMDAVecRestoreArray(da,  mesh->lIAj,  &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,  &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,  &kaj);

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    DMDAVecRestoreArray(da, peqn->lGid, &gid);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Petsc2HypreVector(Vec &A, HYPRE_IJVector &B, HYPRE_Int startID)
{
    HYPRE_Int  size;
    PetscReal *a;

    // get Petsc array
    VecGetArray(A, &a);

    // get local Petsc array dimension
    VecGetLocalSize(A, (PetscInt*)&size);

    std::vector<HYPRE_Int> indices(size);

    // get the global indices to be modified by adding the start index for this
    // processor, computed in setup_lidx2
    for(HYPRE_Int i=0; i<size; i++)
    {
        indices[i] = startID + i;
    }

    // insert the values and assemble the vector to be ready to use
    HYPRE_IJVectorSetValues(B, size, &indices[0], &a[0]);
    HYPRE_IJVectorAssemble(B);

    // restore the Petsc array
    VecRestoreArray(A, &a);

    // clear the vector indices
    std::vector<HYPRE_Int> ().swap(indices);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode Hypre2PetscVector(HYPRE_IJVector &B, Vec &A, HYPRE_Int startID)
{
    HYPRE_Int size;

    // get local Petsc array dimension
    VecGetLocalSize(A, (PetscInt*)&size);

    std::vector<PetscReal>    values(size);
    std::vector<HYPRE_Int> indices(size);

    for(HYPRE_Int i=0; i<size; i++)
    {
        indices[i] = startID + i;
    }

    // get HYPRE array
    HYPRE_IJVectorGetValues(B, size, &indices[0], &values[0]);

    // put the HYPRE values into the Petsc array
    PetscReal *a;

    VecGetArray(A, &a);

    for(HYPRE_Int i=0; i<size; i++)
    {
        a[i] = values[i];
    }

    VecRestoreArray(A, &a);

    // free memory
    std::vector<PetscReal>    ().swap(values);
    std::vector<HYPRE_Int> ().swap(indices);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode phiToPhi(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs   = info.xs, xe = info.xs + info.xm;
    PetscInt           ys   = info.ys, ye = info.ys + info.ym;
    PetscInt           zs   = info.zs, ze = info.zs + info.zm;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k, m;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***nvert, ***meshTag, ***phi, ***lphi, *phi1d;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    VecGetArray(peqn->phi, &phi1d);
    DMDAVecGetArray(da, peqn->Phi, &phi);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    PetscInt pos = 0;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(isIBMCell(k,j,i,nvert) || isOversetCell(k,j,i,meshTag))
                {
                    phi[k][j][i] = 0;
                    pos++;
                }
                else
                {
                    phi[k][j][i] = phi1d[pos++];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da, peqn->Phi, &phi);
    VecRestoreArray    (peqn->phi, &phi1d);

    DMGlobalToLocalBegin(da, peqn->Phi, INSERT_VALUES, peqn->lPhi);
    DMGlobalToLocalEnd  (da, peqn->Phi, INSERT_VALUES, peqn->lPhi);

    UpdatePhiBCs(peqn);

    return(0);
}

//***************************************************************************************************************//

PetscReal L2NormHypre(peqn_ *peqn, HYPRE_IJMatrix &A, HYPRE_IJVector &X, HYPRE_IJVector &B)
{
    mesh_         *mesh  = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs   = info.xs, xe = info.xs + info.xm;
    PetscInt           ys   = info.ys, ye = info.ys + info.ym;
    PetscInt           zs   = info.zs, ze = info.zs + info.zm;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k, m;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscInt           idx[19];

    PetscReal        locSumSq = 0;
    PetscReal        glbSumSq;
    PetscReal        locBsumSq = 0;
    PetscReal        glbBsumSq;

    HYPRE_Int     nsets = 1;

    PetscReal     ***gid;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lGid, &gid);

    // compute the relative residual as the L2 norm of r = |Ax - b| / |b|.
    // - for each row evaluate Ax for this process using sparse matrix format
    // - for each row subtract b and take 2nd power
    // - sum up results from all rows belonging to this process
    // - reduce by summing other processors row summations
    // - take the square root, divide by L2 norm of b and return
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if
                (
                        i == 0 || i == mx-1 ||
                        j == 0 || j == my-1 ||
                        k == 0 || k == mz-1
                )
                {
                    continue;
                }

                HYPRE_Int row = matID(i,j,k);

                idx[CP] = matID(i  , j  , k  );
                idx[EP] = matID(i+1, j  , k  );
                idx[WP] = matID(i-1, j  , k  );
                idx[NP] = matID(i  , j+1, k  );
                idx[SP] = matID(i  , j-1, k  );
                idx[TP] = matID(i  , j  , k+1);
                idx[BP] = matID(i  , j  , k-1);
                idx[NE] = matID(i+1, j+1, k  );
                idx[SE] = matID(i+1, j-1, k  );
                idx[NW] = matID(i-1, j+1, k  );
                idx[SW] = matID(i-1, j-1, k  );
                idx[TN] = matID(i  , j+1, k+1);
                idx[BN] = matID(i  , j+1, k-1);
                idx[TS] = matID(i  , j-1, k+1);
                idx[BS] = matID(i  , j-1, k-1);
                idx[TE] = matID(i+1, j  , k+1);
                idx[BE] = matID(i+1, j  , k-1);
                idx[TW] = matID(i-1, j  , k+1);
                idx[BW] = matID(i-1, j  , k-1);

                // dot product corresponding to current row
                PetscReal dot = 0;

                // loop over non-zero elements
                for(m=0; m<19; m++)
                {
                    // get the matrix column
                    HYPRE_Int col = idx[m];

                    // test if can access the elements
                    if(col >= 0 && row >= 0)
                    {
                        // get matrix and vector element corresponding to current row
                        PetscReal matVal, lhsVal;
                        HYPRE_IJMatrixGetValues(A, nsets, &nsets, &row, &col, &matVal);
                        HYPRE_IJVectorGetValues(X, nsets, &col, &lhsVal);

                        // multiply and accumulate on dot product
                        dot += matVal*lhsVal;
                    }
                }

                // get the element corresponding to current row from rhs
                PetscReal rhsVal;
                HYPRE_IJVectorGetValues(B, nsets, &row, &rhsVal);

                // take power of 2 and accumulate on this process summation
                locSumSq += pow(dot - rhsVal,2.0);

                // take the power of 2 of the RHS and accumulate on this process summation
                locBsumSq += pow(rhsVal,2.0);
            }
        }
    }

    DMDAVecRestoreArray(da, peqn->lGid, &gid);

    // sum contributions from other processes
    MPI_Allreduce(&locSumSq,  &glbSumSq,  1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&locBsumSq, &glbBsumSq, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    // take the square root and return
    if(glbSumSq - glbBsumSq < 1e-10)
    {
        return(1.0);
    }
    else
    {
        return(sqrt(glbSumSq)/sqrt(glbBsumSq));
    }
}

//***************************************************************************************************************//

PetscErrorCode SubtractAverageHypre(peqn_ *peqn, HYPRE_IJVector &B)
{
    mesh_         *mesh  = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs   = info.xs, xe = info.xs + info.xm;
    PetscInt           ys   = info.ys, ye = info.ys + info.ym;
    PetscInt           zs   = info.zs, ze = info.zs + info.zm;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***nvert, ***gid, ***meshTag;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lGid, &gid);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    PetscReal lsum = 0, gsum = 0;
    PetscInt    lcount = 0, gcount = 0;

    // local sum the vector values
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal val;

                // ensure gid differs from -1
                if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    HYPRE_Int idx = matID(i, j, k);
                    HYPRE_IJVectorGetValues(B, 1, &idx, &val);

                    lsum += val;
                    lcount++;
                }

            }
        }
    }

    // global number of counts among processors
    MPI_Allreduce(&lsum, &gsum, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    // global number of counts among processors
    MPI_Allreduce(&lcount, &gcount, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    // compute the mean of the values in the vector
    PetscReal val = -gsum/(PetscReal)gcount;

    // subtract the mean to every component of the vector
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    HYPRE_Int idx = matID(i, j, k);
                    HYPRE_IJVectorAddToValues(B, 1, &idx, &val);
                }
            }
        }
    }
    HYPRE_IJVectorAssemble(B);

    DMDAVecRestoreArray(da, peqn->lGid, &gid);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SubtractAveragePETSc(peqn_ *peqn, Vec &B)
{
    mesh_         *mesh  = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt           xs   = info.xs, xe = info.xs + info.xm;
    PetscInt           ys   = info.ys, ye = info.ys + info.ym;
    PetscInt           zs   = info.zs, ze = info.zs + info.zm;
    PetscInt           mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;
    PetscInt           size;
    PetscReal          ***nvert, ***gid, *b;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lGid, &gid);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    VecGetArray(B, &b);

    // get local Petsc array dimension
    VecGetLocalSize(B, &size);

    PetscReal lsum = 0, gsum = 0;
    PetscInt    lcount = 0, gcount = 0;

    for(PetscInt id=0; id<size; id++)
    {
        lsum += b[id];
        lcount++;
    }

    // global number of counts among processors
    MPI_Allreduce(&lsum, &gsum, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    // global number of counts among processors
    MPI_Allreduce(&lcount, &gcount, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

    // compute the mean of the values in the vector
    PetscReal val = -gsum/(PetscReal)gcount;

    for(PetscInt id=0; id<size; id++)
    {
        b[id] += val;
    }

    VecRestoreArray(B, &b);
    DMDAVecRestoreArray(da, peqn->lGid, &gid);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    VecAssemblyBegin(B);
    VecAssemblyEnd  (B);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetRHS(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    clock_        *clock = peqn->access->clock;
    ueqn_         *ueqn  = peqn->access->ueqn;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt           xs     = info.xs, xe = info.xs + info.xm;
    PetscInt           ys     = info.ys, ye = info.ys + info.ym;
    PetscInt           zs     = info.zs, ze = info.zs + info.zm;
    PetscInt           mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt           i, j, k, size;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscReal        dt   = clock->dt;

    PetscReal     ***nvert, ***gid, *rhs, ***meshTag;
    Cmpnts        ***ucont;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(da,  mesh->lNvert,  &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,  &meshTag);
    DMDAVecGetArray(da,  peqn->lGid,    &gid);

    if(peqn->solverType == "PETSc")
    {
        VecGetArray(peqn->petscRhs, &rhs);
        VecGetLocalSize(peqn->petscRhs, &size);
    }

    PetscReal lsum     = 0.0, gsum  = 0.0;
    PetscInt    lcount   = 0, gcount  = 0;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal val;

                // i, j, k is on solid body
                if (isIBMCell(k, j, i, nvert) || isOversetCell(k,j,i,meshTag))
                {
                    val = 0;
                }
                // i, j, k is a fluid point
                else
                {
                    // val is the divergence of the i, j, k cell.
                    //Contravariant fluxes balance on the cell is computed.
                    val = 0;

                    // i+1/2
                    val += ucont[k][j][i].x;

                    // i-1/2
                    val -= ucont[k][j][i-1].x;

                    // j+1/2
                    val += ucont[k][j][i].y;

                    // j-1/2
                    val -= ucont[k][j-1][i].y;

                    // k+1/2
                    val += ucont[k][j][i].z;

                    // k-1/2
                    val -= ucont[k-1][j][i].z;

                    val *=  1.0 / dt;
                }

                if(peqn->solverType == "HYPRE")
                {
                    // get connectivity
                    HYPRE_Int idx = matID(i,j,k);

                    // set the RHS value
                    HYPRE_IJVectorSetValues(peqn->hypreRhs, 1, &idx, &val);

                    // sum up this cell value
                    if ((!isIBMCell(k, j, i, nvert)) && (!isOversetCell(k,j,i,meshTag)))
                    {
                        lsum += val;
                        lcount++;
                    }
                }
                else if(peqn->solverType == "PETSc")
                {
                    // get connectivity
                    PetscInt idx = (PetscInt)matID(i,j,k) - peqn->thisRankStart;

                    // set the RHS value
                    rhs[idx] = val;

                    // sum up this cell value
                    if ((!isIBMCell(k, j, i, nvert)) && (!isOversetCell(k,j,i,meshTag)))
                    {
                        lsum += val;
                        lcount++;
                    }

                }
            }
        }
    }

    MPI_Allreduce(&lsum,   &gsum,   1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lcount, &gcount, 1, MPIU_INT,    MPI_SUM, mesh->MESH_COMM);

    gsum = - gsum / (PetscReal)gcount;

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(matID(i,j,k)>=0)
                {

                    if ((!isIBMCell(k, j, i, nvert)) && (!isOversetCell(k,j,i,meshTag)))
                    {
                        if(peqn->solverType == "HYPRE")
                        {
                            // get connectivity
                            HYPRE_Int idx = matID(i,j,k);

                            // subtract the mean to apply compatibility condition
                            HYPRE_IJVectorAddToValues(peqn->hypreRhs, 1, &idx, &gsum);
                        }
                        else if(peqn->solverType == "PETSc")
                        {
                            // get connectivity
                            PetscInt idx = (PetscInt)matID(i,j,k) - peqn->thisRankStart;

                            // set the RHS value
                            rhs[idx] += gsum;
                        }
                    }

                }
            }
        }
    }


    if(peqn->solverType == "HYPRE")
    {
        HYPRE_IJVectorAssemble(peqn->hypreRhs);
    }
    if(peqn->solverType == "PETSc")
    {
        VecAssemblyBegin(peqn->petscRhs);
    	VecAssemblyEnd  (peqn->petscRhs);

        VecRestoreArray(peqn->petscRhs, &rhs);
    }

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(da,  mesh->lNvert,  &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,  &meshTag);
    DMDAVecRestoreArray(da,  peqn->lGid,    &gid);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode AdjustIBMFlux(peqn_ *peqn)
{
    mesh_         *mesh = peqn->access->mesh;
    ueqn_         *ueqn = peqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     epsilon = 1.e-10; // threshold below which the velocity is zeroed at IBM boundaries

    PetscReal     ***nvert, ***meshTag;
    Cmpnts        ***ucor;

    PetscReal     local_influx = 0.0, local_outflux = 0.0;
    PetscReal     global_influx = 0.0, global_outflux = 0.0, global_flux;

    // indices for internal cells
    lxs = xs; if (lxs==0) lxs++; lxe = xe; if (lxe==mx) lxe--;
    lys = ys; if (lys==0) lys++; lye = ye; if (lye==my) lye--;
    lzs = zs; if (lzs==0) lzs++; lze = ze; if (lze==mz) lze--;

    DMDAVecGetArray(fda, ueqn->Ucont, &ucor);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // First pass: Compute local influx (to IBM) and outflux (from IBM)
    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                // Left boundary (fluid at i,j,k; IBM at i+1,j+1,k+1)
                if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    // i-left boundary
                    if ((isIBMFluidCell(k, j, i+1, nvert) || isInterpolatedCell(k, j, i+1, meshTag)) && i < mx - 2)
                    {
                        if (fabs(ucor[k][j][i].x) > epsilon)
                        {
                            if (ucor[k][j][i].x > 0) // Flow to IBM (u_x > 0)
                            {
                                local_influx += fabs(ucor[k][j][i].x);
                            }
                            else // Flow from IBM (u_x < 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].x);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].x = 0.0;
                        }
                    }

                    // j-left boundary
                    if ((isIBMFluidCell(k, j+1, i, nvert) || isInterpolatedCell(k, j+1, i, meshTag)) && j < my - 2)
                    {
                        if (fabs(ucor[k][j][i].y) > epsilon)
                        {
                            if (ucor[k][j][i].y > 0) // Flow to IBM (u_y > 0)
                            {
                                local_influx += fabs(ucor[k][j][i].y);
                            }
                            else // Flow from IBM (u_y < 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].y);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].y = 0.0;
                        }
                    }

                    // k-left boundary
                    if ((isIBMFluidCell(k+1, j, i, nvert) || isInterpolatedCell(k+1, j, i, meshTag)) && k < mz - 2)
                    {
                        if (fabs(ucor[k][j][i].z) > epsilon)
                        {
                            if (ucor[k][j][i].z > 0) // Flow to IBM (u_z > 0)
                            {
                                local_influx += fabs(ucor[k][j][i].z);
                            }
                            else // Flow from IBM (u_z < 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].z);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].z = 0.0;
                        }
                    }
                }

                // Right boundary (IBM at i,j,k; fluid at i+1,j+1,k+1)
                if (isIBMFluidCell(k, j, i, nvert) || isInterpolatedCell(k, j, i, meshTag))
                {
                    // i-right boundary
                    if (isFluidCell(k, j, i+1, nvert) && isCalculatedCell(k, j, i+1, meshTag) && i < mx - 2)
                    {
                        if (fabs(ucor[k][j][i].x) > epsilon)
                        {
                            if (ucor[k][j][i].x < 0) // Flow to IBM (u_x < 0)
                            {
                                local_influx += fabs(ucor[k][j][i].x);
                            }
                            else // Flow from IBM (u_x > 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].x);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].x = 0.0;
                        }
                    }

                    // j-right boundary
                    if (isFluidCell(k, j+1, i, nvert) && isCalculatedCell(k, j+1, i, meshTag) && j < my - 2)
                    {
                        if (fabs(ucor[k][j][i].y) > epsilon)
                        {
                            if (ucor[k][j][i].y < 0) // Flow to IBM (u_y < 0)
                            {
                                local_influx += fabs(ucor[k][j][i].y);
                            }
                            else // Flow from IBM (u_y > 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].y);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].y = 0.0;
                        }
                    }

                    // k-right boundary
                    if (isFluidCell(k+1, j, i, nvert) && isCalculatedCell(k+1, j, i, meshTag) && k < mz - 2)
                    {
                        if (fabs(ucor[k][j][i].z) > epsilon)
                        {
                            if (ucor[k][j][i].z < 0) // Flow to IBM (u_z < 0)
                            {
                                local_influx += fabs(ucor[k][j][i].z);
                            }
                            else // Flow from IBM (u_z > 0)
                            {
                                local_outflux += fabs(ucor[k][j][i].z);
                            }
                        }
                        else
                        {
                            ucor[k][j][i].z = 0.0;
                        }
                    }
                }
            }
        }
    }

    // Compute global net flux before correction
    MPI_Allreduce(&local_influx, &global_influx, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&local_outflux, &global_outflux, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    global_flux = global_influx - global_outflux;

    // Second pass: Apply global flux-proportional correction
    for (k = lzs; k < lze; k++)
    {
        for (j = lys; j < lye; j++)
        {
            for (i = lxs; i < lxe; i++)
            {
                // Left boundary (fluid at i,j,k; IBM at i+1,j+1,k+1)
                if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                {
                    // i-left boundary
                    if ((isIBMFluidCell(k, j, i+1, nvert) || isInterpolatedCell(k, j, i+1, meshTag)) && i < mx - 2)
                    {
                        if (fabs(ucor[k][j][i].x) > epsilon)
                        {
                            if (ucor[k][j][i].x > 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].x -= 0.5 * global_flux * (ucor[k][j][i].x / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].x += 0.5 * global_flux * (ucor[k][j][i].x / global_outflux);
                            }
                        }
                    }

                    // j-left boundary
                    if ((isIBMFluidCell(k, j+1, i, nvert) || isInterpolatedCell(k, j+1, i, meshTag)) && j < my - 2)
                    {
                        if (fabs(ucor[k][j][i].y) > epsilon)
                        {
                            if (ucor[k][j][i].y > 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].y -= 0.5 * global_flux * (ucor[k][j][i].y / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].y += 0.5 * global_flux * (ucor[k][j][i].y / global_outflux);
                            }
                        }
                    }

                    // k-left boundary
                    if ((isIBMFluidCell(k+1, j, i, nvert) || isInterpolatedCell(k+1, j, i, meshTag)) && k < mz - 2)
                    {
                        if (fabs(ucor[k][j][i].z) > epsilon)
                        {
                            if (ucor[k][j][i].z > 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].z -= 0.5 * global_flux * (ucor[k][j][i].z / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].z += 0.5 * global_flux * (ucor[k][j][i].z / global_outflux);
                            }
                        }
                    }
                }

                // Right boundary (IBM at i,j,k; fluid at i+1,j+1,k+1)
                if (isIBMFluidCell(k, j, i, nvert) || isInterpolatedCell(k, j, i, meshTag))
                {
                    // i-right boundary
                    if (isFluidCell(k, j, i+1, nvert) && isCalculatedCell(k, j, i+1, meshTag) && i < mx - 2)
                    {
                        if (fabs(ucor[k][j][i].x) > epsilon)
                        {
                            if (ucor[k][j][i].x < 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].x -= 0.5 * global_flux * (ucor[k][j][i].x / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].x += 0.5 * global_flux * (ucor[k][j][i].x / global_outflux);
                            }
                        }
                    }

                    // j-right boundary
                    if (isFluidCell(k, j+1, i, nvert) && isCalculatedCell(k, j+1, i, meshTag) && j < my - 2)
                    {
                        if (fabs(ucor[k][j][i].y) > epsilon)
                        {
                            if (ucor[k][j][i].y < 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].y -= 0.5 * global_flux * (ucor[k][j][i].y / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].y += 0.5 * global_flux * (ucor[k][j][i].y / global_outflux);
                            }
                        }
                    }

                    // k-right boundary
                    if (isFluidCell(k+1, j, i, nvert) && isCalculatedCell(k+1, j, i, meshTag) && k < mz - 2)
                    {
                        if (fabs(ucor[k][j][i].z) > epsilon)
                        {
                            if (ucor[k][j][i].z < 0) // Influx (to IBM)
                            {
                                if (global_influx > 1e-15)
                                    ucor[k][j][i].z -= 0.5 * global_flux * (ucor[k][j][i].z / global_influx);
                            }
                            else // Outflux (from IBM)
                            {
                                if (global_outflux > 1e-15)
                                    ucor[k][j][i].z += 0.5 * global_flux * (ucor[k][j][i].z / global_outflux);
                            }
                        }
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucor);

    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    // // Third pass: Recompute flux after correction for debugging
    // local_influx = 0.0; local_outflux = 0.0;

    // DMDAVecGetArray(fda, ueqn->Ucont, &ucor);
    // DMDAVecGetArray(da, mesh->lNvert, &nvert);
    // DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // for (k = lzs; k < lze; k++)
    // {
    //     for (j = lys; j < lye; j++)
    //     {
    //         for (i = lxs; i < lxe; i++)
    //         {
    //             // Left boundary
    //             if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
    //             {
    //                 // i-left boundary
    //                 if ((isIBMFluidCell(k, j, i+1, nvert) || isInterpolatedCell(k, j, i+1, meshTag)) && i < mx - 2)
    //                 {
    //                     if (ucor[k][j][i].x > 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].x);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].x);
    //                     }
    //                 }

    //                 // j-left boundary
    //                 if ((isIBMFluidCell(k, j+1, i, nvert) || isInterpolatedCell(k, j+1, i, meshTag)) && j < my - 2)
    //                 {
    //                     if (ucor[k][j][i].y > 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].y);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].y);
    //                     }
    //                 }

    //                 // k-left boundary
    //                 if ((isIBMFluidCell(k+1, j, i, nvert) || isInterpolatedCell(k+1, j, i, meshTag)) && k < mz - 2)
    //                 {
    //                     if (ucor[k][j][i].z > 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].z);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].z);
    //                     }
    //                 }
    //             }

    //             // Right boundary
    //             if (isIBMFluidCell(k, j, i, nvert) || isInterpolatedCell(k, j, i, meshTag))
    //             {
    //                 // i-right boundary
    //                 if (isFluidCell(k, j, i+1, nvert) && isCalculatedCell(k, j, i+1, meshTag) && i < mx - 2)
    //                 {
    //                     if (ucor[k][j][i].x < 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].x);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].x);
    //                     }
    //                 }

    //                 // j-right boundary
    //                 if (isFluidCell(k, j+1, i, nvert) && isCalculatedCell(k, j+1, i, meshTag) && j < my - 2)
    //                 {
    //                     if (ucor[k][j][i].y < 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].y);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].y);
    //                     }
    //                 }

    //                 // k-right boundary
    //                 if (isFluidCell(k+1, j, i, nvert) && isCalculatedCell(k+1, j, i, meshTag) && k < mz - 2)
    //                 {
    //                     if (ucor[k][j][i].z < 0) // Influx (to IBM)
    //                     {
    //                         local_influx += fabs(ucor[k][j][i].z);
    //                     }
    //                     else // Outflux (from IBM)
    //                     {
    //                         local_outflux += fabs(ucor[k][j][i].z);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    // DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    // DMDAVecRestoreArray(fda, ueqn->Ucont, &ucor);

    // // Compute global net flux after correction
    // MPI_Allreduce(&local_influx, &global_influx, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    // MPI_Allreduce(&local_outflux, &global_outflux, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    // global_flux = global_influx - global_outflux;

    // // Print global flux for debugging
    // PetscPrintf(PETSC_COMM_WORLD, "Post-correction global net flux: %e (influx to IBM: %e, outflux from IBM: %e)\n",
    //             global_flux, global_influx, global_outflux);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ProjectVelocity(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    ueqn_         *ueqn  = peqn->access->ueqn;
    clock_        *clock = peqn->access->clock;
    DM            da     = mesh->da, fda = mesh->fda;
    DMDALocalInfo info   = mesh->info;
    PetscInt      xs     = info.xs, xe = info.xs + info.xm;
    PetscInt      ys     = info.ys, ye = info.ys + info.ym;
    PetscInt      zs     = info.zs, ze = info.zs + info.zm;
    PetscInt      mx     = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    Cmpnts        ***icsi, ***ieta, ***izet;
    Cmpnts        ***jcsi, ***jeta, ***jzet;
    Cmpnts        ***kcsi, ***keta, ***kzet;

    PetscReal     ***iaj, ***jaj, ***kaj;
    PetscReal     ***phi, ***nvert, ***meshTag;

    Cmpnts        ***ucont;

    PetscReal     dpdc, dpde, dpdz;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda,  mesh->lICsi, &icsi);
    DMDAVecGetArray(fda,  mesh->lIEta, &ieta);
    DMDAVecGetArray(fda,  mesh->lIZet, &izet);
    DMDAVecGetArray(fda,  mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda,  mesh->lJEta, &jeta);
    DMDAVecGetArray(fda,  mesh->lJZet, &jzet);
    DMDAVecGetArray(fda,  mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda,  mesh->lKEta, &keta);
    DMDAVecGetArray(fda,  mesh->lKZet, &kzet);
    DMDAVecGetArray(da,   mesh->lIAj, &iaj);
    DMDAVecGetArray(da,   mesh->lJAj, &jaj);
    DMDAVecGetArray(da,   mesh->lKAj, &kaj);
    DMDAVecGetArray(da,   mesh->lNvert, &nvert);
    DMDAVecGetArray(da,   mesh->lmeshTag, &meshTag);

    DMDAVecGetArray(da,   peqn->lPhi, &phi);

    DMDAVecGetArray(fda,  ueqn->Ucont, &ucont);

    PetscInt periodic_i = mesh->i_periodic + mesh->ii_periodic,
             periodic_j = mesh->j_periodic + mesh->jj_periodic,
             periodic_k = mesh->k_periodic + mesh->kk_periodic;

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                if
                (
                    j > 0 && k > 0 &&              // exclude non-defined indices in the i-direction
                    (i>0 && i<mx-2 || periodic_i)  // exclude also first and last face unless periodic (has meaningful ghosts)
                )
                {
                    dpdc = phi[k][j][i + 1] - phi[k][j][i];

                    dpde = 0.;
                    dpdz = 0.;

                    if
                    (
                        // j-right boundary -> use upwind
                        (
                            j==my-2
                        ) || isIBMIFace(k, j+1, i, i+1, nvert) || isOversetIFace(k, j+1, i, i+1, meshTag)
                    )
                    {
                        dpde = (phi[k][j][i] + phi[k][j][i + 1] - phi[k][j - 1][i] - phi[k][j - 1][i + 1]) * 0.5;
                    }
                    else if
                    (
                        // j-left boundary -> use upwind
                        (
                            j == 1
                        ) || isIBMIFace(k, j-1, i, i+1, nvert) || isOversetIFace(k, j-1, i, i+1, meshTag)
                     )
                     {
                         dpde = (phi[k][j+1][i] + phi[k][j+1][i+1] - phi[k][j][i] - phi[k][j][i+1])* 0.5;
                     }
                     else
                     {
                         // central differences
                         dpde = (phi[k][j+1][i] + phi[k][j+1][i+1] - phi[k][j-1][i] - phi[k][j-1][i+1])* 0.25;
                     }

                     if
                     (
                        // k-right boundary -> use upwind
                        (
                            k == mz - 2
                        ) || isIBMIFace(k+1, j, i, i+1, nvert) || isOversetIFace(k+1, j, i, i+1, meshTag)
                    )
                    {
                        dpdz = (phi[k][j][i] + phi[k][j][i+1] - phi[k-1][j][i] - phi[k-1][j][i+1]) * 0.5;
                    }
                    else if
                    (
                        // k-left boundary  -> use upwind
                        (
                            k == 1
                        ) || isIBMIFace(k-1, j, i, i+1, nvert) || isOversetIFace(k-1, j, i, i+1, meshTag)
                    )
                    {
                        dpdz = (phi[k+1][j][i] + phi[k+1][j][i+1]- phi[k][j][i] - phi[k][j][i+1]) * 0.5;
                    }
                    else
                    {
                        dpdz = (phi[k+1][j][i] + phi[k+1][j][i+1] - phi[k-1][j][i] - phi[k-1][j][i+1]) * 0.25;
                    }

                    if
                    (
                        isFluidIFace(k, j, i, i+1, nvert) && isCalculatedIFace(k, j, i, i+1, meshTag)
                    )
                    {
                        ucont[k][j][i].x
                        -=
                        (
                            dpdc *
                            (
                                icsi[k][j][i].x * icsi[k][j][i].x +
                                icsi[k][j][i].y * icsi[k][j][i].y +
                                icsi[k][j][i].z * icsi[k][j][i].z
                            ) * iaj[k][j][i] +
                            dpde *
                            (
                                ieta[k][j][i].x * icsi[k][j][i].x +
                                ieta[k][j][i].y * icsi[k][j][i].y +
                                ieta[k][j][i].z * icsi[k][j][i].z
                            ) * iaj[k][j][i] +
                            dpdz *
                            (
                                izet[k][j][i].x * icsi[k][j][i].x +
                                izet[k][j][i].y * icsi[k][j][i].y +
                                izet[k][j][i].z * icsi[k][j][i].z
                            ) * iaj[k][j][i]
                        ) * clock->dt;
                    }
                }

                if
                (
                    i > 0 && k > 0 &&              // exclude non-defined indices in the j-direction
                    (j>0 && j<my-2 || periodic_j)  // exclude also first and last face unless periodic (has meaningful ghosts)
                )
                {
                    dpdc = 0.;
                    dpdz = 0.;

                    if
                    (
                        // i-right boundary -> use upwind
                        (
                            i == mx-2
                        ) || isIBMJFace(k, j, i+1, j+1, nvert) || isOversetJFace(k, j, i+1, j+1, meshTag)
                    )
                    {
                        dpdc = (phi[k][j][i] + phi[k][j+1][i] - phi[k][j][i-1] - phi[k][j+1][i-1]) * 0.5;
                    }
                    else if
                    (
                        // i-left boundary -> use upwind
                        (
                            i == 1
                        ) || isIBMJFace(k, j, i-1, j+1, nvert) || isOversetJFace(k, j, i-1, j+1, meshTag)
                    )
                    {
                        dpdc = (phi[k][j][i+1] + phi[k][j+1][i+1] - phi[k][j][i] - phi[k][j+1][i]) * 0.5;
                    }
                    else
                    {
                        // central differencing
                        dpdc = (phi[k][j][i+1] + phi[k][j+1][i+1] - phi[k][j][i-1] - phi[k][j+1][i-1]) * 0.25;
                    }

                    dpde = phi[k][j+1][i] - phi[k][j][i];

                    if
                    (
                        // k-right boundary -> use upwind
                        (
                            k == mz-2
                        ) || isIBMJFace(k+1, j, i, j+1, nvert) || isOversetJFace(k+1, j, i, j+1, meshTag)
                    )
                    {
                        dpdz = (phi[k][j][i] + phi[k][j+1][i] - phi[k-1][j][i] - phi[k - 1][j + 1][i]) * 0.5;
                    }
                    else if
                    (
                        // k-left boundary -> use upwind
                        (
                            k == 1
                        ) || isIBMJFace(k-1, j, i, j+1, nvert) || isOversetJFace(k-1, j, i, j+1, meshTag)
                    )
                    {
                        dpdz = (phi[k + 1][j][i] + phi[k + 1][j + 1][i] - phi[k][j][i] - phi[k][j + 1][i]) * 0.5;
                    }
                    else
                    {
                        // central differencing
                        dpdz = (phi[k + 1][j][i] + phi[k + 1][j + 1][i] - phi[k - 1][j][i] - phi[k - 1][j + 1][i]) * 0.25;
                    }

                    if
                    (
                            isFluidJFace(k, j, i, j+1, nvert) && isCalculatedJFace(k, j, i, j+1, meshTag)
                    )
                    {
                        ucont[k][j][i].y
                        -=
                        (
                            dpdc *
                            (
                                jcsi[k][j][i].x * jeta[k][j][i].x +
                                jcsi[k][j][i].y * jeta[k][j][i].y +
                                jcsi[k][j][i].z * jeta[k][j][i].z
                            ) * jaj[k][j][i] +
                            dpde *
                            (
                                jeta[k][j][i].x * jeta[k][j][i].x +
                                jeta[k][j][i].y * jeta[k][j][i].y +
                                jeta[k][j][i].z * jeta[k][j][i].z
                            ) * jaj[k][j][i] +
                            dpdz *
                            (
                                jzet[k][j][i].x * jeta[k][j][i].x +
                                jzet[k][j][i].y * jeta[k][j][i].y +
                                jzet[k][j][i].z * jeta[k][j][i].z
                            ) * jaj[k][j][i]
                        ) * clock->dt;
                    }
                }

                if
                (
                    i > 0 && j > 0 &&              // exclude non-defined indices in the k-direction
                    (k>0 && k<mz-2 || periodic_k)  // exclude also first and last face unless periodic (has meaningful ghosts)
                )
                {
                    dpdc = 0.;
                    dpde = 0.;

                    if
                    (
                        // i-right boundary -> use upwind
                        (
                            i == mx - 2
                        ) || isIBMKFace(k, j, i+1, k+1, nvert) || isOversetKFace(k, j, i+1, k+1, meshTag)
                    )
                    {
                        dpdc = (phi[k][j][i] + phi[k+1][j][i] - phi[k][j][i-1] - phi[k+1][j][i-1]) * 0.5;
                    }
                    else if
                    (
                        // i-left boundary -> use upwind
                        (
                            i == 1
                        ) || isIBMKFace(k, j, i-1, k+1, nvert) || isOversetKFace(k, j, i-1, k+1, meshTag)
                    )
                    {
                        dpdc = (phi[k][j][i+1] + phi[k+1][j][i+1] - phi[k][j][i] - phi[k+1][j][i])* 0.5;
                    }
                    else
                    {
                        dpdc = (phi[k][j][i+1] + phi[k+1][j][i+1] - phi[k][j][i-1] - phi[k+1][j][i-1]) * 0.25;
                    }

                    if
                    (
                        // j-right boundary -> use upwind
                        (
                            j == my - 2
                        ) || isIBMKFace(k, j+1, i, k+1, nvert) || isOversetKFace(k, j+1, i, k+1, meshTag)
                    )
                    {
                        dpde = (phi[k][j][i] + phi[k+1][j][i] - phi[k][j-1][i] - phi[k+1][j-1][i]) * 0.5;
                    }
                    else if
                    (
                        // j-left boundary -> use upwind
                        (
                            j == 1
                        ) || isIBMKFace(k, j-1, i, k+1, nvert) || isOversetKFace(k, j-1, i, k+1, meshTag)
                    )
                    {
                        dpde = (phi[k][j+1][i] + phi[k+1][j+1][i] - phi[k][j][i] - phi[k+1][j][i]) * 0.5;
                    }
                    else
                    {
                        // central differences
                        dpde = (phi[k][j + 1][i] + phi[k + 1][j + 1][i] - phi[k][j - 1][i] - phi[k + 1][j - 1][i]) * 0.25;
                    }

                    dpdz = phi[k + 1][j][i] - phi[k][j][i];

                    if
                    (
                        isFluidKFace(k, j, i, k+1, nvert) && isCalculatedKFace(k, j, i, k+1, meshTag)
                    )
                    {
                        ucont[k][j][i].z
                        -=
                        (
                            dpdc *
                            (
                                kcsi[k][j][i].x * kzet[k][j][i].x +
                                kcsi[k][j][i].y * kzet[k][j][i].y +
                                kcsi[k][j][i].z * kzet[k][j][i].z
                            ) * kaj[k][j][i] +
                            dpde *
                            (
                                keta[k][j][i].x * kzet[k][j][i].x +
                                keta[k][j][i].y * kzet[k][j][i].y +
                                keta[k][j][i].z * kzet[k][j][i].z
                            ) * kaj[k][j][i] +
                            dpdz *
                            (
                                kzet[k][j][i].x * kzet[k][j][i].x +
                                kzet[k][j][i].y * kzet[k][j][i].y +
                                kzet[k][j][i].z * kzet[k][j][i].z
                            ) * kaj[k][j][i]
                        ) * clock->dt;

                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda,  mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda,  mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda,  mesh->lIZet, &izet);
    DMDAVecRestoreArray(fda,  mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda,  mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda,  mesh->lJZet, &jzet);
    DMDAVecRestoreArray(fda,  mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda,  mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda,  mesh->lKZet, &kzet);
    DMDAVecRestoreArray(da,   mesh->lIAj, &iaj);
    DMDAVecRestoreArray(da,   mesh->lJAj, &jaj);
    DMDAVecRestoreArray(da,   mesh->lKAj, &kaj);
    DMDAVecRestoreArray(da,   mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,   mesh->lmeshTag, &meshTag);

    DMDAVecRestoreArray(da,   peqn->lPhi, &phi);

    DMDAVecRestoreArray(fda,  ueqn->Ucont, &ucont);

    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdatePressure(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    clock_        *clock = peqn->access->clock;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    PetscReal     ***p, ***lphi, ***nvert, ***meshTag;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscReal lPsum    = 0.0, gPsum  = 0.0;
    PetscInt  lnPoints = 0, gnPoints = 0;

    DMDAVecGetArray(da, peqn->P, &p);
    DMDAVecGetArray(da, peqn->lPhi, &lphi);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // update pressure and average

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
              if (isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag)) 
              {
                      continue;
              }

              if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
              {
                  // p updated only for fluid cells. For ibm fluid cells it is interpolated from the fluid cells in IBM interpolation function

                  p[k][j][i] += lphi[k][j][i];
                  lPsum += p[k][j][i];
                  lnPoints++;
              }
            }
        }
    }

	// update boundary conditions
	UpdatePressureBCs(peqn);

    MPI_Allreduce(&lPsum,    &gPsum,    1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lnPoints, &gnPoints, 1, MPIU_INT,    MPI_SUM, mesh->MESH_COMM);

    gPsum = gPsum / gnPoints;

    // subtract average
    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
              if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
              {
                  p[k][j][i] -= gPsum;
              }
              else if(isIBMFluidCell(k, j, i, nvert))
              {
                 //do nothing. will be updated in ibm interpolation
              }
              else
              {
                   p[k][j][i] = 0.0;
              }
            }
        }
    }

    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(da, peqn->lPhi, &lphi);
    DMDAVecRestoreArray(da, peqn->P, &p);

    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd  (da, peqn->P, INSERT_VALUES, peqn->lP);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode GradP(peqn_ *peqn)
{
    mesh_         *mesh = peqn->access->mesh;
    ueqn_         *ueqn = peqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***icsi, ***ieta, ***izet,
                  ***jcsi, ***jeta, ***jzet,
                  ***kcsi, ***keta, ***kzet,
                  ***dp;

    PetscReal     ***p, ***nvert, ***meshTag;
    PetscReal     ***iaj, ***jaj, ***kaj;

    PetscInt           i, j, k;
    PetscInt           lxs, lxe, lys, lye, lzs, lze;

    PetscReal        dpdc, dpde, dpdz;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

	PetscInt periodic_i = mesh->i_periodic + mesh->ii_periodic,
             periodic_j = mesh->j_periodic + mesh->jj_periodic,
             periodic_k = mesh->k_periodic + mesh->kk_periodic;

    VecSet(ueqn->dP, 0.0);

    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lIEta, &ieta);
    DMDAVecGetArray(fda, mesh->lIZet, &izet);
    DMDAVecGetArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lJZet, &jzet);
    DMDAVecGetArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecGetArray(fda, mesh->lKEta, &keta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);
    DMDAVecGetArray(da,  mesh->lNvert,&nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag,&meshTag);
    DMDAVecGetArray(da,  mesh->lIAj,  &iaj);
    DMDAVecGetArray(da,  mesh->lJAj,  &jaj);
    DMDAVecGetArray(da,  mesh->lKAj,  &kaj);

    DMDAVecGetArray(da,  peqn->lP,  &p );
    DMDAVecGetArray(fda, ueqn->dP, &dp);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal g11_i = (icsi[k][j][i].x * icsi[k][j][i].x + icsi[k][j][i].y * icsi[k][j][i].y + icsi[k][j][i].z * icsi[k][j][i].z);
                PetscReal g12_i = (ieta[k][j][i].x * icsi[k][j][i].x + ieta[k][j][i].y * icsi[k][j][i].y + ieta[k][j][i].z * icsi[k][j][i].z);
                PetscReal g13_i = (izet[k][j][i].x * icsi[k][j][i].x + izet[k][j][i].y * icsi[k][j][i].y + izet[k][j][i].z * icsi[k][j][i].z);
                PetscReal g21_j = (jcsi[k][j][i].x * jeta[k][j][i].x + jcsi[k][j][i].y * jeta[k][j][i].y + jcsi[k][j][i].z * jeta[k][j][i].z);
                PetscReal g22_j = (jeta[k][j][i].x * jeta[k][j][i].x + jeta[k][j][i].y * jeta[k][j][i].y + jeta[k][j][i].z * jeta[k][j][i].z);
                PetscReal g23_j = (jzet[k][j][i].x * jeta[k][j][i].x + jzet[k][j][i].y * jeta[k][j][i].y + jzet[k][j][i].z * jeta[k][j][i].z);
                PetscReal g31_k = (kcsi[k][j][i].x * kzet[k][j][i].x + kcsi[k][j][i].y * kzet[k][j][i].y + kcsi[k][j][i].z * kzet[k][j][i].z);
                PetscReal g32_k = (keta[k][j][i].x * kzet[k][j][i].x + keta[k][j][i].y * kzet[k][j][i].y + keta[k][j][i].z * kzet[k][j][i].z);
                PetscReal g33_k = (kzet[k][j][i].x * kzet[k][j][i].x + kzet[k][j][i].y * kzet[k][j][i].y + kzet[k][j][i].z * kzet[k][j][i].z);

                // pressure gradient in the i-direction
                if( i==mx-2 && mesh->ii_periodic)
                {
                    dpdc = p[k][j][mx+1] - p[k][j][i];
                }
                else
                {
                    dpdc = p[k][j][i+1] - p[k][j][i];
                }

                if
                (
                    // j-right boundary -> use upwind only at the corner faces
                    (
                        j==my-2 &&
                        (
                            i==mx-2
                        )
                    ) || isOversetIFace(k, j+1, i, i+1, meshTag)
                )
                {
                    dpde = (p[k][j  ][i  ] - p[k][j-1][i  ] + p[k][j  ][i+1] - p[k][j-1][i+1]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (
                        j == 1 &&
                        (
                            i == mx - 2
                        )
                     ) || isOversetIFace(k, j-1, i, i+1, meshTag)
                )
                {
                    dpde = (p[k][j+1][i  ] - p[k][j  ][i  ] + p[k][j+1][i+1] - p[k][j  ][i+1]) * 0.5;
                }
                else
                {
                    dpde = (p[k][j+1][i] - p[k][j-1][i] + p[k][j+1][i+1] - p[k][j-1][i+1]) * 0.25;
                }

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (
                        k == mz - 2 &&
                        (
                            i==mx-2
                        )
                    ) || isOversetIFace(k+1, j, i, i+1, meshTag)
                )
                {
                    dpdz = (p[k][j][i  ] - p[k-1][j][i  ] + p[k][j][i+1] - p[k-1][j][i+1]) * 0.5;
                }
                else if
                (
                    // k-left boundary  -> use upwind  only at the corner faces
                    (
                        k == 1 &&
                        (
                            i==mx-2
                        )
                    ) || isOversetIFace(k-1, j, i, i+1, meshTag)
                )
                {
                    dpdz = (p[k+1][j][i  ] - p[k][j][i  ] + p[k+1][j][i+1] - p[k][j][i+1]) * 0.5;
                }
                else
                {
                    dpdz = (p[k+1][j][i] - p[k-1][j][i] + p[k+1][j][i+1] - p[k-1][j][i+1]) * 0.25;
                }

                dp[k][j][i].x = (dpdc * g11_i + dpde *  g12_i + dpdz * g13_i ) * iaj[k][j][i];

                // pressure gradient in the j-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (
                        i == mx-2 &&
                        (
                            j==my-2
                        )
                    ) || isOversetJFace(k, j, i+1, j+1, meshTag)
                )
                {
                    dpdc = (p[k][j  ][i] - p[k][j  ][i-1] + p[k][j+1][i] - p[k][j+1][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (
                        i == 1 &&
                        (
                            j==my-2
                        )
                    ) || isOversetJFace(k, j, i-1, j+1, meshTag)
                )
                {
                    dpdc = (p[k][j  ][i+1] - p[k][j  ][i] + p[k][j+1][i+1] - p[k][j+1][i]) * 0.5;
                }
                else
                {
                    dpdc = (p[k][j  ][i+1] - p[k][j  ][i-1] + p[k][j+1][i+1] - p[k][j+1][i-1]) * 0.25;
                }

                if( j==my-2 && mesh->jj_periodic)
                {
                    dpde = p[k][my+1][i] - p[k][j][i];
                }
                else
                {
                    dpde = p[k][j+1][i] - p[k][j][i];
                }

                if
                (
                    // k-right boundary -> use upwind  only at the corner faces
                    (
                        k == mz-2 &&
                        (
                            j== my-2
                        )
                    ) || isOversetJFace(k+1, j, i, j+1, meshTag)
                )
                {
                    dpdz = (p[k][j  ][i] - p[k-1][j  ][i] + p[k][j+1][i] - p[k-1][j+1][i]) * 0.5;
                }
                else if
                (
                    // k-left boundary -> use upwind  only at the corner faces
                    (
                        k == 1 &&
                        (
                            j== my-2
                        )
                    ) || isOversetJFace(k-1, j, i, j+1, meshTag)
                )
                {
                    dpdz = (p[k+1][j  ][i] - p[k][j  ][i] + p[k+1][j+1][i] - p[k][j+1][i]) * 0.5;
                }
                else
                {
                    dpdz = (p[k+1][j  ][i] - p[k-1][j  ][i] + p[k+1][j+1][i] - p[k-1][j+1][i]) * 0.25;
                }

                dp[k][j][i].y = (dpdc * g21_j + dpde * g22_j + dpdz * g23_j ) * jaj[k][j][i];

                // pressure gradient in the k-direction
                if
                (
                    // i-right boundary -> use upwind  only at the corner faces
                    (
                        i == mx - 2 &&
                        (
                            k==mz-2
                        )
                    ) || isOversetKFace(k, j, i+1, k+1, meshTag)
                )
                {
                    dpdc = (p[k  ][j][i] - p[k  ][j][i-1] + p[k+1][j][i] - p[k+1][j][i-1]) * 0.5;
                }
                else if
                (
                    // i-left boundary -> use upwind  only at the corner faces
                    (
                        i == 1 &&
                        (
                            k == mz - 2
                        )
                    ) || isOversetKFace(k, j, i-1, k+1, meshTag)
                )
                {
                    dpdc = (p[k  ][j][i+1] - p[k  ][j][i] + p[k+1][j][i+1] - p[k+1][j][i]) * 0.5;
                }
                else
                {
                    dpdc = (p[k  ][j][i+1] - p[k  ][j][i-1] + p[k+1][j][i+1] - p[k+1][j][i-1]) * 0.25;
                }

                if
                (
                    // j-right boundary -> use upwind  only at the corner faces
                    (
                        j == my - 2 &&
                        (
                            k==mz-2
                        )
                    ) || isOversetKFace(k, j+1, i, k+1, meshTag)
                )
                {
                    dpde = (p[k  ][j][i] - p[k  ][j-1][i] + p[k+1][j][i] - p[k+1][j-1][i]) * 0.5;
                }
                else if
                (
                    // j-left boundary -> use upwind  only at the corner faces
                    (
                        j == 1 &&
                        (
                            k==mz-2
                        )
                    )  || isOversetKFace(k, j-1, i, k+1, meshTag)
                )
                {
                    dpde = (p[k  ][j+1][i] - p[k  ][j][i] + p[k+1][j+1][i] - p[k+1][j][i]) * 0.5;
                }
                else
                {
                    dpde = (p[k  ][j+1][i] - p[k  ][j-1][i] + p[k+1][j+1][i] - p[k+1][j-1][i]) * 0.25;
                }

                if( k==mz-2 && mesh->kk_periodic)
                {
                    dpdz = p[mz+1][j][i] - p[k][j][i];
                }
                else
                {
                    dpdz = (p[k+1][j][i] - p[k][j][i]);
                }

                dp[k][j][i].z = (dpdc * g31_k + dpde * g32_k + dpdz * g33_k ) * kaj[k][j][i];

                // periodic: set to zero only on left boundaries since the contrav. velocity is not solved there
                // non-periodic: set to zero also on right boundaries since the contrav. velocity is not solved there
                if
                (
                    i==0 || (!mesh->i_periodic && !mesh->ii_periodic && i==mx-2) || isOversetIFace(k, j, i, i+1, meshTag)
                )
                {
                    dp[k][j][i].x = 0;
                }
                if
                (
                    j==0 || (!mesh->j_periodic && !mesh->jj_periodic && j==my-2) || isOversetJFace(k, j, i, j+1, meshTag)
                )
                {
                    dp[k][j][i].y = 0;
                }
                if
                (
                    k==0 || (!mesh->k_periodic && !mesh->kk_periodic && k==mz-2) || isOversetKFace(k, j, i, k+1, meshTag)
                )
                {
                    dp[k][j][i].z = 0;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lIEta, &ieta);
    DMDAVecRestoreArray(fda, mesh->lIZet, &izet);
    DMDAVecRestoreArray(fda, mesh->lJCsi, &jcsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lJZet, &jzet);
    DMDAVecRestoreArray(fda, mesh->lKCsi, &kcsi);
    DMDAVecRestoreArray(fda, mesh->lKEta, &keta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);
    DMDAVecRestoreArray(da,  mesh->lNvert,&nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag,&meshTag);
    DMDAVecRestoreArray(da,  mesh->lIAj,  &iaj);
    DMDAVecRestoreArray(da,  mesh->lJAj,  &jaj);
    DMDAVecRestoreArray(da,  mesh->lKAj,  &kaj);

    DMDAVecRestoreArray(da,  peqn->lP, &p );
    DMDAVecRestoreArray(fda, ueqn->dP, &dp);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SolvePEqn(peqn_ *peqn)
{
    mesh_         *mesh  = peqn->access->mesh;
    clock_        *clock = peqn->access->clock;
    flags_        *flags = peqn->access->flags;
    PetscReal     ts, te;

    PetscMPIInt   rank, size;

    // get current process
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    PetscTime(&ts);

    // add flux correction
    if(flags->isIBMActive || flags->isOversetActive)
    {
        AdjustIBMFlux(peqn);

        MPI_Barrier(mesh->MESH_COMM);
    }

    if(flags->isIBMActive)
    {
        if( peqn->access->ibm->dynamic )
        {
            if(peqn->solverType == "HYPRE")
            {
                DestroyHypreMatrix(peqn);
                DestroyHypreVector(peqn);


                CreateHypreMatrix(peqn);
                CreateHypreVector(peqn);

                MPI_Barrier(mesh->MESH_COMM);
            }

            // set coefficient matrix
            SetCoeffMatrix(peqn);

            MPI_Barrier(mesh->MESH_COMM);

            if(peqn->solverType == "PETSc")
            {
                // destroy rhs framework
                VecSet(peqn->petscRhs, 0.0);

                // create solver framework
                if(clock->it > clock->itStart)
                {
                    DestroyPETScSolver(peqn);
                }

                CreatePETScSolver(peqn);
            }
        }
    }

    // compute the RHS (divergence of predicted velocity)
    SetRHS(peqn);

    if(peqn->solverType == "HYPRE")
    {
        // transform Phi2 to the unknown in HYPRE (used as initial guess)
        Petsc2HypreVector(peqn->phi, peqn->hypreP, peqn->thisRankStart);

        MPI_Barrier(mesh->MESH_COMM);

        // initialize the solver
        if((clock->it == clock->itStart) || (flags->isIBMActive && peqn->access->ibm->dynamic))
        {
            if(peqn->hypreSolverType == 1)
            {
                PetscPrintf(mesh->MESH_COMM, "BoomerGMRES: Setup = ");
                HYPRE_ParCSRGMRESSetup(peqn->hypreSlvr, peqn->hypreParA, peqn->hypreParRhs, peqn->hypreParP);
                PetscPrintf(mesh->MESH_COMM, "1, Solving for p, ");
            }
            else if (peqn->hypreSolverType == 2)
            {
                PetscPrintf(mesh->MESH_COMM, "BoomerPCG: Setup = ");
                HYPRE_ParCSRPCGSetup  (peqn->hypreSlvr, peqn->hypreParA, peqn->hypreParRhs, peqn->hypreParP);
                PetscPrintf(mesh->MESH_COMM, "1, Solving for p, ");
            }
        }
        else
        {
            if(peqn->hypreSolverType == 1) PetscPrintf(mesh->MESH_COMM, "BoomerGMRES: ");
            else if (peqn->hypreSolverType == 2) PetscPrintf(mesh->MESH_COMM, "BoomerPCG: ");
        }

        MPI_Barrier(mesh->MESH_COMM);

        // compute initial residual
        peqn->initialPoissonRes = L2NormHypre(peqn, peqn->hypreA, peqn->hypreP, peqn->hypreRhs);

        // solve the Poisson equation
        if(peqn->hypreSolverType == 1)
        {
            HYPRE_ParCSRGMRESSolve
            (
                    peqn->hypreSlvr,
                    peqn->hypreParA,
                    peqn->hypreParRhs,
                    peqn->hypreParP
            );

            // get iteration number
            HYPRE_GMRESGetNumIterations(peqn->hypreSlvr, &(peqn->poissonIterations));

            // compute final relative residual norm
            HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm(peqn->hypreSlvr, &(peqn->finalPoissonRes));

            PetscPrintf(mesh->MESH_COMM, "Initial residual = %e, Final residual = %e, Iterations = %ld, Elapsed Time = ", peqn->initialPoissonRes, peqn->finalPoissonRes, peqn->poissonIterations);
        }
        else if (peqn->hypreSolverType == 2)
        {
            HYPRE_ParCSRPCGSolve
            (
                    peqn->hypreSlvr,
                    peqn->hypreParA,
                    peqn->hypreParRhs,
                    peqn->hypreParP
            );

            // get iteration number
            HYPRE_PCGGetNumIterations(peqn->hypreSlvr, &(peqn->poissonIterations));

            // compute final relative residual norm
            HYPRE_PCGGetFinalRelativeResidualNorm(peqn->hypreSlvr, &(peqn->finalPoissonRes));

            PetscPrintf(mesh->MESH_COMM, "Initial residual = %e, Final residual = %e, Iterations = %ld, Elapsed Time = ", peqn->initialPoissonRes, peqn->finalPoissonRes, peqn->poissonIterations);
        }

        // subtract the mean from the solution vector
        SubtractAverageHypre(peqn, peqn->hypreP);

        // transform the HYPRE solution to the Phi2 petsc vector
        Hypre2PetscVector(peqn->hypreP, peqn->phi, peqn->thisRankStart);
    }
    else if(peqn->solverType == "PETSc")
    {
        // compute initial residual
        KSPSolve(peqn->ksp, peqn->petscRhs, peqn->phi);

        // compute info
        KSPGetResidualNorm(peqn->ksp,&(peqn->finalPoissonRes));
        PetscInt numIter;
        KSPGetIterationNumber(peqn->ksp, &numIter);
        peqn->poissonIterations = (HYPRE_Int)numIter;

        PetscPrintf(mesh->MESH_COMM, "MGGMRES: Solving for p, Initial residual = %e, Final residual = %e, Iterations = %ld, Elapsed Time = ",peqn->initialPoissonRes, peqn->finalPoissonRes, peqn->poissonIterations);

        if(flags->isIBMActive)
        {
            if( peqn->access->ibm->dynamic )
            {
                MatMPIAIJSetPreallocation(peqn->petscA, 27, PETSC_NULL, 27, PETSC_NULL);

                MatSetOption(peqn->petscA, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
            }
        }

        SubtractAveragePETSc(peqn, peqn->phi);
    }

    // wait until al processes reach this point
    MPI_Barrier(mesh->MESH_COMM);

    // convert Phi2 to Phi (the actual correction)
    phiToPhi(peqn);

    // wait until al processes reach this point
    MPI_Barrier(mesh->MESH_COMM);

    // update pressure
    UpdatePressure(peqn);

    // set pressure reference
    SetPressureReference(peqn);

    // update velocity
    ProjectVelocity(peqn);

    PetscTime(&te);

    PetscPrintf(mesh->MESH_COMM, "%lf\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetPressureReference(peqn_ *peqn)
{
    mesh_         *mesh = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***p, ***lp, ***nvert, ***meshTag;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lP, &lp);
    DMDAVecGetArray(da, peqn->P, &p);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    // test if this processor contains the 0,0,0 cell
    PetscReal lscale = 0.0, gscale = 0.0;
    if(xs==0 && ys == 0 && zs==0)
    {
        lscale = p[1][1][1];
    }

    MPI_Allreduce(&lscale, &gscale, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if((!isIBMSolidCell(k, j, i, nvert)) && (!isZeroedCell(k, j, i, meshTag)))
                {
                    p[k][j][i] -= gscale;
                }
            }
        }
    }

    DMDAVecRestoreArray(da, peqn->lP, &lp);
    DMDAVecRestoreArray(da, peqn->P, &p);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

    // scatter Phi from global to local
    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd  (da, peqn->P, INSERT_VALUES, peqn->lP);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ContinuityErrors(peqn_ *peqn)
{
    mesh_         *mesh = peqn->access->mesh;
    ueqn_         *ueqn = peqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, r;

    Vec           Div;
    PetscReal     ***div, ***aj, ***nvert, ***meshTag;
    Cmpnts        ***ucont, ***ucat;

    PetscReal     maxdiv;
    PetscReal     lmaxu = 0, gmaxu = 0;

    PetscMPIInt   rank, size;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // get current process
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // get ttal number of processes
    MPI_Comm_size(mesh->MESH_COMM, &size);

    VecDuplicate(peqn->P, &Div);
    VecSet(Div, 1.0);

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecGetArray(da,  mesh->lAj, &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecGetArray(da,  Div, &div);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                maxdiv
                =
                fabs
                (
                    (
                        ucont[k][j][i].x - ucont[k][j][i-1].x +
                        ucont[k][j][i].y - ucont[k][j-1][i].y +
                        ucont[k][j][i].z - ucont[k-1][j][i].z
                    )*aj[k][j][i]
                );

                if(nMag(ucat[k][j][i]) > lmaxu)
                {
                    lmaxu = nMag(ucat[k][j][i]);
                }

                // if
                // (
                //     nvert[k  ][j  ][i  ] +
                //     nvert[k+1][j  ][i  ] +
                //     nvert[k-1][j  ][i  ] +
                //     nvert[k  ][j+1][i  ] +
                //     nvert[k  ][j-1][i  ] +
                //     nvert[k  ][j  ][i+1] +
                //     nvert[k  ][j  ][i-1] > 0.1
                // )
                // {
                //     maxdiv = 0.;
                // }

                if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
                {
                    maxdiv = 0;
                }

                div[k][j][i] = maxdiv;
            }
        }
    }

    std::vector<PetscReal>  lmaxDiv(size);
    std::vector<PetscReal>  gmaxDiv(size);
    cellIds              maxDivIds;

    maxDivIds.i = 1; maxDivIds.j = 1; maxDivIds.k = 1;

    for(r=0; r<size; r++)
    {
        lmaxDiv[r] = 0.0;
        gmaxDiv[r] = 0.0;
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                if(div[k][j][i] > lmaxDiv[rank])
                {
                    lmaxDiv[rank] = div[k][j][i];
                    maxDivIds.i = i;
                    maxDivIds.j = j;
                    maxDivIds.k = k;
                }
            }
        }
    }

    MPI_Allreduce(&lmaxDiv[0], &gmaxDiv[0], size, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);
    MPI_Allreduce(&lmaxu,      &gmaxu,      1,    MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);

    maxdiv      = 0.0;
    PetscInt maxrank = 0;

    for(r=0; r<size; r++)
    {
        if(gmaxDiv[r] > maxdiv)
        {
            maxdiv   = gmaxDiv[r];
            maxrank  = r;
        }
    }

    if(rank==maxrank) printf("Time step continuity error: %e (at cell %ld, %ld, %ld). uMax = %lf\n", maxdiv, maxDivIds.i, maxDivIds.j, maxDivIds.k, gmaxu);

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

    DMDAVecRestoreArray(da,  mesh->lNvert,  &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj,     &aj);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecRestoreArray(da, Div, &div);

    VecDestroy(&Div);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode ContinuityErrorsOptimized(peqn_ *peqn)
{
    mesh_         *mesh = peqn->access->mesh;
    ueqn_         *ueqn = peqn->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k, r;

    PetscReal     ***aj, ***nvert, ***meshTag;
    Cmpnts        ***ucont;

    PetscReal     lmaxdiv = 0.0, gmaxdiv = 0.0;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, ueqn->lUcont, &ucont);
    DMDAVecGetArray(da,  mesh->lAj, &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscReal div
                =
                fabs
                (
                    (
                        ucont[k][j][i].x - ucont[k][j][i-1].x +
                        ucont[k][j][i].y - ucont[k][j-1][i].y +
                        ucont[k][j][i].z - ucont[k-1][j][i].z
                    )*aj[k][j][i]
                );

                // if(isBoxIBMCell(k, j, i, nvert) || isBoxOversetCell(k, j, i, meshTag))
                // {
                //     div = 0;
                // }

                if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
                {
                    div = 0;
                }

                if(div > lmaxdiv) lmaxdiv = div;
            }
        }
    }

    MPI_Allreduce(&lmaxdiv, &gmaxdiv, 1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);

    PetscPrintf(mesh->MESH_COMM, "Time step continuity error: %e\n", gmaxdiv);

    DMDAVecRestoreArray(fda, ueqn->lUcont, &ucont);
    DMDAVecRestoreArray(da,  mesh->lNvert,  &nvert);
    DMDAVecRestoreArray(da,  mesh->lAj,     &aj);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

    MPI_Barrier(mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//
