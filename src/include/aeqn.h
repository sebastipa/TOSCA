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

    Vec           lDivA, lViscA, lViscIBMA;   //!< viscous and divergence alpha water equation fluxes
    Vec           sourceA;                    //!< alpha water sources

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
