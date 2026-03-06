//! \file  ueqn.h
//! \brief U equation solution header file.

#ifndef _UEQN_H_
#define _UEQN_H_

//! \brief structure storing momentum equation
struct ueqn_
{
    // implicit time stepping
    SNES          snesU;                      //!< non linear matrix free context for momentum equation
    Mat           JU;                         //!< non linear matrix free preconditioner
    Mat           A, C;
    KSP           ksp;                        //!< linear krylov-subspace context (backwardEuler SNES inner KSP)
    PC            pc;
    KSP           kspIMEX;                    //!< standalone direct KSP for IMEX scheme (no SNES wrapper)
    Mat           JvIMEX;                     //!< MatShell for IMEX operator A*v = v - dt*scale*Visc(v)
    word          snesType;                   //!< SNES solver type (kept for logging)
    PetscInt      snesMaxIter;                //!< max SNES outer iterations (-snesMaxItersU in control.dat)
    word          kspType;                    //!< KSP solver type: BiCGStab or GMRES (-kspTypeU in control.dat)
    PetscInt      gmresRestart;               //!< GMRES restart parameter (-kspGMRESRestartU in control.dat, default 30)

    // momentum variables
    Vec           Utmp;                             //!< temporary solution passed to the SNES evaluation function or used for RK4
    Vec           Rhs;                              //!< rhs of the momentum equation (stores transport and viscous fluxes), low level use in FormU
    Vec           Rhs_o;                            //!< rhs of the momentum equation at previous time step
    Vec           lFp;                              //!< rhs of the momentum equation prior to dotting with curv. coords basis (becomes Rhs)
    Vec           lDiv1, lDiv2, lDiv3;              //!< Components of the convective term in momentum equation
    Vec           lVisc1, lVisc2, lVisc3;           //!< Components of the viscous term in the momentum equation
    Vec           lViscIBM1, lViscIBM2, lViscIBM3;  //!< Components of the viscous term in the momentum equation for IBM faces
    Vec           dP;                               //!< pressure term of the momentum equation
	Vec           dPAGW;                            //!< pressure term of the momentum equation as calculated from provided atmopsheric gravity waves pressure
    Vec           bTheta;                           //!< buoyancy field
    Vec           sourceU;                          //!< source term to drive Uref at Zref with periodic BCs
    Vec           gCont;                            //!< gravity vector in cuvilinear coordinates
    Vec           Ucont, lUcont;                    //!< contravariant fluxes (contravariant velocity / J)
    Vec           Ucat, lUcat;                      //!< cartesian velocity
    Vec           Ucont_o;                          //!< contravariant fluxes at the previous time step
    Vec           bU;                               //!< precomputed constant-RHS contributions (backwardEuler/IMEX): built once per time step in SolveUEqn
    Vec           RhsConv;                          //!< convective-only RHS at current step n  (IMEX-CNAB: Adams-Bashforth 2)
    Vec           RhsConv_o;                        //!< convective-only RHS at previous step n-1 (IMEX-CNAB: Adams-Bashforth 2)
    word          teqnPredictorScheme;               //!< T predictor before SolveUEqn for improved buoyancy coupling: "forwardEuler" or "none" (default, -teqnPredictorU)
    Vec           lUstar;
    Cmpnts        meanGradP;
    Cmpnts        uBulk;

    // momentum settings
    word          ddtScheme;                  //!< time derivative scheme
    word          divScheme;                  //!< divergence scheme
    word          viscScheme;
    PetscReal     relExitTol;                 //!< relative exit tolerance
    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscInt      inviscid;                   //!< inviscid run
    PetscInt      buoyancy;                   //!< buoyancy term
    PetscInt      coriolis;                   //!< coriolis term
    PetscInt      fringe;                     //!< fringe region term

    // divergence schemes
    PetscInt      centralDiv;                 //!< linear divergence scheme
    PetscInt      central4Div;                //!< 4th order central scheme
    PetscInt      centralUpwindDiv;           //!< blending between linear and upwind scheme
    PetscInt      centralUpwindWDiv;          //!< blending between linear and upwind scheme for non-uniform mesh
    PetscInt      quickDiv;                   //!< 3rd order QUICK scheme
    PetscInt      weno3Div;                   //!< 3rd order WENO scheme

    PetscReal     hyperVisc;                  //!< hyperviscocity parameter to add artificial diffusion - (b = 0 is blend between 3rd order upwind and 4th order, b = 1 is central 4th order, b = 0.8 is hybrid 3-4 scheme)
    PetscReal     hyperVisc4;                 //!< IMEX explicit biharmonic (4th-order index-space) hyperviscosity coefficient: per-step damping ε₄, contribution to ΔU = -ε₄·δ⁴(lUcont). Targets Nyquist (checkerboard) modes only (-imexHyperVisc4U, default 0)
                                              
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

//! \brief Transform velocity from contravariant to cartesian
PetscErrorCode contravariantToCartesian(ueqn_ *ueqn);

//! \brief Transform generic local vector from contravariant to cartesian
PetscErrorCode contravariantToCartesianGeneric(mesh_ *mesh, Vec &lCont, Vec &lCat);

//! \brief Adjust fluxes to obey mass conservation (on a per-cell basis)
PetscErrorCode adjustFluxesLocal(ueqn_ *ueqn);

//! \brief Solve the momentum equation
PetscErrorCode SolveUEqn(ueqn_ *ueqn);
