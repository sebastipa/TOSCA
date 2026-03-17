//! \file  teqn.h
//! \brief T equation solution header file.

#ifndef _TEQN_H_
#define _TEQN_H_

//! \brief struct storing temperature equation
struct teqn_
{
    // solver name 
    word          solverType;                 //!< solver type (kept for logging)
    
    // temperature variables
    SNES          snesT;                      //!< non linear matrix free context
    Mat           JT;                         //!< non linear matrix free preconditioner
    Mat           AT, CT;
    KSP           ksp;                        //!< linear krylov-subspace context (backwardEuler/BDF2 SNES inner KSP)
    PC            pc;
    KSP           kspIMEX;                    //!< standalone direct KSP for IMEX schemes
    Mat           JvIMEX;                     //!< MatShell for IMEX operator A*v 
    word          kspType;                    //!< KSP solver type for SNES/IMEX
    PetscInt      gmresRestart;               //!< GMRES restart parameter for SNES/IMEX
    Vec           Rhs;                        //!< rhs of the temperature equation 
    Vec           bT;                         //!< rhs of the linear system in IMEX schemes
    Vec           RhsConv;                    //!< convective RHS at current time (for AB2 blending in ABBE)
    Vec           RhsConv_o;                  //!< convective RHS at previous time (for AB2 blending in ABBE)

    Vec           TmprtTmp;                   //!< temporary solution
    Vec           Tmprt, lTmprt,
                  Tmprt_o, lTmprt_o;
    Vec           Tmprt_oo;                   //!< previous solution for BDF2

    Vec           lDivT, lViscT, lViscIBMT;   //!< viscous and divergence temperature equation fluxes
    Vec           sourceT;                    //!< temperature sources

    Vec           lRhoK;                      //!< rhok / rho0 field
    Vec           ghGradRhok;                 //!< buoyancy term for momentum equation (current step)
    Vec           ghGradRhok_o;               //!< buoyancy term for momentum equation (previous step, for AB2)

    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscReal     relExitTol;                 //!< relative exit tolerance

    word          ddtScheme;                  //!< time derivative scheme

    PetscReal     hyperVisc4i;                //!< IMEX biharmonic hyperviscosity along i (ξ) index direction (-imexHyperVisc4U_i, default 0)
    PetscReal     hyperVisc4j;                //!< IMEX biharmonic hyperviscosity along j (η) index direction (-imexHyperVisc4U_j, default 0)
    PetscReal     hyperVisc4k;                //!< IMEX biharmonic hyperviscosity along k (ζ) index direction, typically vertical (-imexHyperVisc4U_k, default 0)

    // wall model patch
    wallModel     *iLWM;                      //!< wall model on the i-left patch
    wallModel     *iRWM;                      //!< wall model on the i-right patch
    wallModel     *jLWM;                      //!< wall model on the j-left patch
    wallModel     *jRWM;                      //!< wall model on the j-right patch

    // t equation flags
    PetscInt      pTildeFormulation;          //!< buoyancy term is expressed as gradient in momentum

    // initial field
    word          initFieldType;

    // access
    access_       *access;                     //!< access database

};

#endif

//! \brief Initializes Teqn environment
PetscErrorCode InitializeTEqn(teqn_ *teqn);

//! \brief Solve T equation
PetscErrorCode SolveTEqn(teqn_ *teqn);

//! \brief Computes g*h times gradient of rho_k / rho_0
PetscErrorCode ghGradRhoK(teqn_ *teqn);

//! \brief Forward-Euler conv-only T predictor: T* = T^n + dt*N(T^n) -> lTmprt (for buoyancy in SolveUEqn)
PetscErrorCode TmprtPredictor(teqn_ *teqn);

//! \brief Restore Tmprt/lTmprt to T^n from Tmprt_o (undo TmprtPredictor before SolveTEqn)
PetscErrorCode TmprtRestoreFromOld(teqn_ *teqn);
