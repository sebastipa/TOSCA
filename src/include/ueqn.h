//! \file  ueqn.h
//! \brief U equation solution header file.

#ifndef _UEQN_H_
#define _UEQN_H_

//! \brief structure storing the Shumann wall models information
struct Shumann
{
    word             wfEvalType;
    PetscReal        roughness;
    
};

//! \brief structure storing the Shumann wall models information
struct Cabot
{
    PetscReal        roughness;
};

//! \brief wall models container
struct wallModel
{
    Shumann          *wmShumann;
    Cabot            *wmCabot;
    patchVectorField tauWall;
};

//! \brief structure storing momentum equation
struct ueqn_
{
    // momentum variables
    SNES          snesU;                      //!< non linear matrix free context for momentum equation
    Mat           JU;                         //!< non linear matrix free preconditioner
    Mat           A, C;
    KSP           ksp;                        //!< linear krylov-subspace context
    PC            pc;
    Vec           Utmp;                       //!< temporary solution passed to the SNES evaluation function
    Vec           Rhs;                        //!< rhs of the momentum equation (stores transport and viscous fluxes), low level use in FormU
    Vec           Rhs_o;                      //!< rhs of the momentum equation at previous time step
    Vec           lFp;                        //!< rhs of the momentum equation prior to dotting with curv. coords basis (becomes Rhs)
    Vec           lDiv1, lDiv2, lDiv3;        //!< Components of the convective term in momentum equation
    Vec           lVisc1, lVisc2, lVisc3;     //!< Components of the viscous term in the momentum equation
    Vec           dP;                         //!< pressure term of the momentum equation
    Vec           sourceU;                    //!< source term to drive Uref at Zref with periodic BCs
    Vec           gCont;                      //!< gravity vector in cuvilinear coordinates
    Vec           Ucont, lUcont;              //!< contravariant fluxes (contravariant velocity / J)
    Vec           Ucat, lUcat;                //!< cartesian velocity
    Vec           Ucont_o;                    //!< contravariant fluxes at the previous time step
    Vec           lUstar;

    // momentum settings
    word          divScheme;                  //!< divergence scheme
    PetscReal     relExitTol;                 //!< relative exit tolerance
    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscInt      inviscid;                   //!< inviscid run
    PetscInt      buoyancy;                   //!< buoyancy term
    PetscInt      coriolis;                   //!< coriolis term
    PetscInt      fringe;                     //!< fringe region term

    // divergence schemes
    PetscInt      centralDiv;                 //!< linear divergence scheme
    PetscInt      centralUpwindDiv;           //!< blending between linear and upwind scheme
    PetscInt      centralUpwindWDiv;          //!< blending between linear and upwind scheme for non-uniform mesh
    PetscInt      quickDiv;                   //!< 3rd order QUICK scheme
    PetscInt      weno3Div;                   //!< 3rd order WENO scheme

    // wall model patch
    wallModel     *iLWM;                      //!< wall model on the i-left patch
    wallModel     *iRWM;                      //!< wall model on the i-right patch
    wallModel     *jLWM;                      //!< wall model on the j-left patch
    wallModel     *jRWM;                      //!< wall model on the j-right patch

    // initial field
    word           initFieldType;

    // access
    access_       *access;                    //!< access database

};

#endif

//! \brief Initializes Ueqn environment
PetscErrorCode InitializeUEqn(ueqn_ *ueqn);

//! \brief Updates flux limiter
PetscErrorCode UpdateFluxLimiter(ueqn_ *ueqn);

//! \brief Update driving source terms
PetscErrorCode CorrectSourceTerms(ueqn_ *ueqn, PetscInt print);

//! \brief Correct damping source terms
PetscErrorCode correctDampingSources(ueqn_ *ueqn);

//! \brief Transform contravariant to cartesian
PetscErrorCode contravariantToCartesian(ueqn_ *ueqn);

//! \brief Adjust fluxes to obey mass conservation
PetscErrorCode adjustFluxes(ueqn_ *ueqn);

PetscErrorCode adjustFluxesOverset(ueqn_ *ueqn);

//! \brief Compute driving source term
PetscErrorCode sourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Compute damping source terms (x and z)
PetscErrorCode dampingSourceU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Compute Coriolis source term
PetscErrorCode Coriolis(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Compute Side Force source term
PetscErrorCode SideForce(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Compute buoyancy term
PetscErrorCode Buoyancy(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! Viscous and divergence terms
PetscErrorCode FormU(ueqn_ *ueqn, Vec &Rhs, PetscReal scale);

//! \brief Solve the momentum equation
PetscErrorCode SolveUEqn(ueqn_ *ueqn);

//! \brief SNES evaluation function
PetscErrorCode UeqnSNES(SNES snes, Vec Ucont, Vec Rhs, void *ptr);
