//! \file  teqn.h
//! \brief T equation solution header file.

#ifndef _TEQN_H_
#define _TEQN_H_

//! \brief struct storing temperature equation
struct teqn_
{
    // temperature variables
    SNES          snesT;                      //!< non linear matrix free context
    Mat           JT;                         //!< non linear matrix free preconditioner
    Mat           AT, CT;
    KSP           ksp;                        //!< linear krylov-subspace context (backwardEuler/crankNicholson SNES inner KSP)
    PC            pc;
    KSP           kspIMEX;                    //!< standalone direct KSP for IMEX scheme (no SNES wrapper)
    Mat           JvIMEX;                     //!< MatShell for IMEX operator A*v = v - dt*D(v)
    word          kspType;                    //!< KSP solver type for IMEX: BiCGStab or GMRES (-kspTypeT in control.dat)
    PetscInt      gmresRestart;               //!< GMRES restart parameter for IMEX (-kspGMRESRestartT in control.dat, default 30)
    Vec           Rhs;
    Vec           Rhs_o;
    Vec           bT;                         //!< explicit-half RHS buffer: prebuilt once per step for IMEX (conv) or crankNicholson (FormT+damp)
    Vec           RhsConv;                    //!< convective RHS at T^n (for AB2 blending)
    Vec           RhsConv_o;                  //!< convective RHS at T^{n-1} (for AB2 blending)

    Vec           TmprtTmp;                   //!< temporary solution
    Vec           Tmprt, lTmprt,
                  Tmprt_o, lTmprt_o;
    Vec           lDivT, lViscT, lViscIBMT;              //!< viscous and divergence temperature equation fluxes
    Vec           sourceT;                    //!< temperature sources

    Vec           lRhoK;                      //!< rhok / rho0 field
    Vec           ghGradRhok;                 //!< buoyancy term for momentum equation

    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscReal     relExitTol;                 //!< relative exit tolerance

    word          ddtScheme;                  //!< time derivative scheme

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
