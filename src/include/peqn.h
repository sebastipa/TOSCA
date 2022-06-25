//! \file  peqn.h
//! \brief P equation solution header file.

#ifndef _PEQN_H_
#define _PEQN_H_

//! \brief struct storing pressure equation
struct peqn_
{
    Vec           phi;                        //!< phi is the actual solution of Poisson equation (then converted in Phi)
    Vec           Phi, lPhi;                  //!< pressure correction for frac. step method
    Vec           P, lP;                      //!< pressure at current and previous time step
    Vec           lLid, lGid;                 //!< matrixFree - matrixBased connectivity (used to build the poisson coeff. matrix)
    Vec           pIBMPt, pIBMCell;           //!< point where the ibm pressure is interpolated and its closest cell id
    HYPRE_Int     thisRankSize,               //!< number of cells owned by this processore (used to build the poisson coeff. matrix)
                  thisRankStart,              //!< first cell ID owned by this processor in global indexing (used to build the poisson coeff. matrix)
                  thisRankEnd;                //!< last cell ID owned by this processor in global indexing (thisRankStart + thisRankSize - 1)
    HYPRE_Int     totalSize;                  //!< total number of cells (excluding physical ghost & IBM points, depending of the Hypre Poisson type, -1, -2, 1)
    PetscReal     initialPoissonRes,
                  finalPoissonRes;            //!< initial and final residual of the Poisson iteration
    HYPRE_Int     poissonIterations;          //!< number of Poisson iterations

    PetscInt      hypreSolverType;            //!< 1: GMRES, 2: PCG
    PetscInt      poissonIt;                  //!< max number of poisson iterations per timestep
    PetscReal     poissonTol;                 //!< relative exit tolerance

    // solver parameters
    PetscInt      amgAgg;                     //!< aggresive coarsening is good for > 50mil grids
    PetscInt      amgCoarsenType;
    PetscReal     amgThresh;                  //!< threshold value - 0.5 : Cartesian, 0.6 : Distorted

    // hypre solver
    HYPRE_Solver       hypreSlvr;             //!< solver
    HYPRE_Solver       PC;                    //!< preconditioner
    HYPRE_IJMatrix     hypreA;                //!< coefficient matrix
    HYPRE_ParCSRMatrix hypreParA;             //!< coefficient matrix
    HYPRE_IJVector     hypreP, hypreRhs;      //!< unknwon and RHS
    HYPRE_ParVector    hypreParP, hypreParRhs;//!< unknown and RHS

    // access
    access_            *access;                 //!< access database

};

#endif

//! \brief Initializes Peqn environment
PetscErrorCode InitializePEqn(peqn_ *peqn);

//! \brief Initializes HYPRE solver
PetscErrorCode CreateHypreSolver(peqn_ *peqn);

//! \brief Initializes HYPRE matrix
PetscErrorCode CreateHypreMatrix(peqn_ *peqn);

//! \brief Initializes HYPRE vector
PetscErrorCode CreateHypreVector(peqn_ *peqn);

//! \brief Destroys HYPRE solver
PetscErrorCode DestroyHypreSolver(peqn_ *peqn);

//! \brief Destroys HYPRE matrix
PetscErrorCode DestroyHypreMatrix(peqn_ *peqn);

//! \brief Destroys HYPRE vector
PetscErrorCode DestroyHypreVector(peqn_ *peqn);

//! \brief Transfer vector values from Petsc 2 Hypre
PetscErrorCode Petsc2HypreVector(Vec &A, HYPRE_IJVector &B, HYPRE_Int startID);

//! \brief Transfer vector values from Hypre 2 Petsc
PetscErrorCode Hypre2PetscVector(HYPRE_IJVector &B, Vec &A, HYPRE_Int startID);

//! \brief Convert phi to Phi (1D to 3D vector)
PetscErrorCode phiToPhi(peqn_ *peqn);

//! \brief Compute L2 norm of system residual
PetscReal         L2NormHypre(peqn_ *peqn, HYPRE_IJMatrix &A, HYPRE_IJVector &X, HYPRE_IJVector &B);

//! \brief Compute pressure gradient term
PetscErrorCode GradP(peqn_ *peqn);

//! \brief Compute cell-matrix connectivity
PetscErrorCode SetPoissonConnectivity(peqn_ *peqn);

//! \brief Compute coefficient matrix
PetscErrorCode SetCoeffMatrix(peqn_ *peqn);

//! \brief Set coefficient matrix to zero
PetscErrorCode ZeroCoeffMatrix(peqn_ *peqn);

//! \brief Set RHS of the pressure equation
PetscErrorCode SetRHS(peqn_ *peqn);

//! \brief Subtract average from the solution
PetscErrorCode SubtractAverage(peqn_ *peqn, HYPRE_IJVector &B);

//! \brief Correct IBM volume flux
PetscErrorCode AdjustIBMFlux(peqn_ *peqn);

//! \brief Update pressure and subtract average
PetscErrorCode UpdatePressure(peqn_ *peqn);

PetscErrorCode updateIBMPhi(ibm_ *ibm);

//! \brief Project Ucont into an incompressible space
PetscErrorCode ProjectVelocity(peqn_ *peqn);

//! \brief Compute pressure gradient term
PetscErrorCode SolvePEqn(peqn_ *peqn);

//! \brief Set pressure = 0 at k = 0, j = 0, i = 0
PetscErrorCode SetPressureReference(peqn_ *peqn);

//! \brief Compute continuity errors (also calculates which cell and processor has the max)
PetscErrorCode ContinuityErrors(peqn_ *peqn);

//! \brief Compute continuity errors (only prints the max)
PetscErrorCode ContinuityErrorsOptimized(peqn_ *peqn);

//! \brief get the cell id from stencil position
cellIds GetIdFromStencil(int stencil, int k, int j, int i);
