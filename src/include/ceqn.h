//! \file  teqn.h
//! \brief T equation solution header file.

#ifndef _CEQN_H_
#define _CEQN_H_

//! \brief struct storing temperature equation
struct ceqn_
{
    // temperature variables
    SNES          snesC;                      //!< non linear matrix free context
    Mat           JC;                         //!< non linear matrix free preconditioner
    Mat           AC, CC;
    KSP           ksp;                        //!< linear krylov-subspace context
    PC            pc;
    Vec           Rhs;
    Vec           Rhs_o;
    Vec           ConcTmp;                   //!< temporary solution
    Vec           Conc, lConc,
                  Conc_o, lConc_o;
    Vec           lDivC, lViscC;              //!< viscous and divergence temperature equation fluxes

    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscReal     relExitTol;                 //!< relative exit tolerance

    word          ddtScheme;                  //!< time derivative scheme

    // wall model patch
    //wallModel     *iLWM;                      //!< wall model on the i-left patch
    //wallModel     *iRWM;                      //!< wall model on the i-right patch
    //wallModel     *jLWM;                      //!< wall model on the j-left patch
    //wallModel     *jRWM;                      //!< wall model on the j-right patch

    // initial field
    word          initFieldType;

    // access
    access_       *access;                     //!< access database

    PetscReal     cRefNoAbl;                    //cRef if no abl

};

#endif

//! \brief Initializes Teqn environment
PetscErrorCode InitializeCEqn(ceqn_ *ceqn);

//! \brief Solve T equation
PetscErrorCode SolveCEqn(ceqn_ *ceqn);

//! \brief SNES evaulation function
PetscErrorCode CeqnSNES(SNES snes, Vec C, Vec Rhs, void *ptr);

//! \brief Apply fringe region damping. Not used in ceqn currently.
//PetscErrorCode dampingSourceC(ceqn_ *ceqn, Vec &Rhs, PetscReal scale);

//! \brief RHS of the potential concentration transport equation
PetscErrorCode FormC(ceqn_ *ceqn, Vec &Rhs, PetscReal scale);

//! \brief solve Teqn using RungeKutta 4
//PetscErrorCode TeqnRK4(teqn_ *teqn);

//! \brief Computed RHS of concentration equation using current lTmprt (updates Rhs), data put in ueqn->Rhs Not used currently for concentration
//PetscErrorCode FormExplicitRhsC(ceqn_ *ceqn);
