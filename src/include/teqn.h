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
    KSP           ksp;                        //!< linear krylov-subspace context
    PC            pc;
    Vec           Rhs;
    Vec           Rhs_o;
    Vec           TmprtTmp;                   //!< temporary solution
    Vec           Tmprt, lTmprt,
                  Tmprt_o, lTmprt_o;
    Vec           lDivT, lViscT;              //!< viscous and divergence temperature equation fluxes
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

//! \brief SNES evaulation function
PetscErrorCode TeqnSNES(SNES snes, Vec T, Vec Rhs, void *ptr);

//! \brief Computes g*h times gradient of rho_k / rho_0
PetscErrorCode ghGradRhoK(teqn_ *teqn);

//! \brief Compute temperature control source term
PetscErrorCode CorrectSourceTermsT(teqn_ *teqn, PetscInt print);

//! \brief Apply fringe region damping
PetscErrorCode dampingSourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

//! \brief Apply temperature control
PetscErrorCode sourceT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

//! \brief RHS of the potential temperature transport equation
PetscErrorCode FormT(teqn_ *teqn, Vec &Rhs, PetscReal scale);

//! \brief solve Teqn using RungeKutta 4
PetscErrorCode TeqnRK4(teqn_ *teqn);

//! \brief Computed RHS of temperature equation using current lTmprt (updates Rhs), data put in ueqn->Rhs
PetscErrorCode FormExplicitRhsT(teqn_ *teqn);
