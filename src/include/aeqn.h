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
    Vec           Alpha, lAlpha,              //!< volume fraction fields
                  Alpha_o, Alpha_oo; 

    Vec           lDivA;                      //!< combined face fluxes
    Vec           lDivALO;                    //!< low-order (upwind) face fluxes
    Vec           lDivAHO;                    //!< high-order (central4) face fluxes
    Vec           lDivACor;                   //!< correction face fluxes: lDivAHO - lDivA
    Vec           lLambdaA;                   //!< mules face limiter coefficients lambda_f in [0,1]
    Vec           lAlphaLO;                   //!< low-order provisional cell update alpha^{LO}
    Vec           lRplusA;                    //!< per-cell R_plus  (max allowable incoming correction fraction)
    Vec           lRminusA;                   //!< per-cell R_minus (max allowable outgoing correction fraction)

    Vec           dRho;                       //!< density gradient (for multiphase flows)
    Vec           lRho;                       //!< density at cell centers (for multiphase flows)
    Vec           lRhoFace;                   //!< density at faces (for multiphase flows)
    Vec           lRhoFace_o;                 //!< density at faces at previous time step (for multiphase flows)

    PetscInt      nAlphaSubCycles;            //!< number of alpha sub-cycling steps (default 3)
    PetscReal     compCoeff;                  //!< numerical compression term coefficient

    PetscReal     absExitTol;                 //!< absolute exit tolerance
    PetscReal     relExitTol;                 //!< relative exit tolerance

    word          ddtScheme;                  //!< time derivative scheme

    // initial field
    word          initFieldType;

    // damping regions
    PetscReal     kRightAlphaDampingDelta;       //!< damping region thickness (m)
    PetscReal     kRightAlphaDampingCoeff;       //!< damping region coefficient (1/s)
    PetscReal     kRightAlphaDampingWaterLevel;  //!< unperturbed water level (m)

    PetscReal     kLeftAlphaDampingStart;        //!< damping region start (m)
    PetscReal     kLeftAlphaDampingEnd;          //!< damping region end (m)
    PetscReal     kLeftAlphaDampingDelta;        //!< damping region rising and decaying length (m)
    PetscReal     kLeftAlphaDampingCoeff;        //!< damping region coefficient (1/s)

    // wave properties to be used for damping regions
    
    // monochromatic wave (linear or stokes2)
    word         waveType;                    //!< wave theory type: "linear" or "stokes2"
    PetscReal    waveHeight;                  //!< wave height H (m)
    PetscReal    wavePeriod;                  //!< wave period T (s)
    PetscReal    waveLevel;                   //!< MWL position in domain coords (m) - also water depth if bottom is at z=0
    PetscReal    waveDirection;               //!< wave propagation angle (degrees, 0=x-dir)
    PetscReal    wavePhase;                   //!< initial phase offset (radians)
    PetscReal    waveNumber;                  //!< computed wave number k (1/m)
    PetscReal    waveOmega;                   //!< computed angular frequency omega (rad/s)
    PetscReal    waveLambda;                  //!< computed wavelength lambda (m)

    // access
    access_       *access;                     //!< access database
};

#endif

//! \brief Initializes alpha equation environment
PetscErrorCode InitializeAEqn(aeqn_ *aeqn);

//! \brief Solve alpha equation
PetscErrorCode SolveAEqn(aeqn_ *aeqn);

//! \brief clip alpha water to [0,1] and sync ghosts
PetscErrorCode boundAlpha(aeqn_ *aeqn);

//! \brief compute density gradient for momentum right hand side 
PetscErrorCode GradRho(aeqn_ *aeqn);