//! \file  aeqn.h
//! \brief alpha equation solution header file.

#ifndef _AEQN_H_
#define _AEQN_H_

//! \brief struct storing alpha water equation
struct aeqn_
{
    // solver name 
    word          solverType;                 //!< solver type (kept for logging)
    
    // temperature variables
    SNES          snesA;                      //!< non linear matrix free context
    Mat           JA;                         //!< non linear matrix free preconditioner
    Mat           AA, CA;
    KSP           ksp;                        //!< linear krylov-subspace context (backwardEuler/BDF2 SNES inner KSP)
    PC            pc;
    word          kspType;                    //!< KSP solver type for SNES/IMEX
    PetscInt      gmresRestart;               //!< GMRES restart parameter for SNES/IMEX
    Vec           Rhs;                        //!< rhs of the alpha water equation 

    Vec           AlphaTmp;                   //!< temporary solution
    Vec           Alpha, lAlpha,
                  Alpha_o, lAlpha_o, Alpha_oo; 

    Vec           lDivA;                      //!< low-order (upwind) face fluxes
    Vec           lDivAHO;                    //!< high-order (central4) face fluxes
    Vec           lDivACor;                   //!< correction face fluxes: lDivAHO - lDivA
    Vec           lLambdaA;                   //!< mules face limiter coefficients lambda_f in [0,1]
    Vec           lAlphaLO;                   //!< low-order provisional cell update alpha^{LO}
    Vec           lRplusA;                    //!< per-cell R_plus  (max allowable incoming correction fraction)
    Vec           lRminusA;                   //!< per-cell R_minus (max allowable outgoing correction fraction)

    PetscInt      mulesIter;                  //!< number of mules limiter sweeps (default 3)

    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscReal     relExitTol;                 //!< relative exit tolerance

    word          ddtScheme;                  //!< time derivative scheme

    // initial field
    word          initFieldType;

    // access
    access_       *access;                     //!< access database

};

#endif

//! \brief Initializes alpha equation environment
PetscErrorCode InitializeAEqn(aeqn_ *aeqn);

//! \brief Solve alpha equation
PetscErrorCode SolveAEqn(aeqn_ *aeqn);

//! rief Update alpha water boundary conditions
PetscErrorCode UpdateAlphaWaterBCs(aeqn_ *aeqn);