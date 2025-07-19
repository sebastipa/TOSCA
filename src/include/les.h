//! \file  lesmdl.h
//! \brief LES model header file.

#ifndef _LES_H_
#define _LES_H_

typedef enum {
    SMAGORINSKY,
    STABILITY_BASED,
    DSM,
    DLASI,
    DLASD,
    DPASD,
    AMD,
    VREMAN,
    BV,
    BAMD,
    BDS
} LesModel;

//! \brief struct storing LES model
struct les_
{
    LesModel      model;
    PetscReal     maxCs;
    PetscReal     amdCs;
    PetscInt      turbStat;

    // Dynamic Smagorisnky
    Vec           lSx, lSy, lSz, lS;
    Vec           lLM, lMM;
    Vec           lQN, lNN;
    Vec           lLM_old, lMM_old;
    Vec           lQN_old, lNN_old;
    Vec           lCs;

    // Scale dependent model
    Vec           lN;
    Vec           lCh;
    Vec           L;                       //!< length scale

    // kEquation

    // general
    Vec           lNu_t;                  //!< eddy viscosity
    Vec           lDiff_t;
    Vec           lKsgs;

    Vec           lTau;                   //!< tau sgs tensor for structural models

    // initial nut field
    word          initFieldType;

    // access
    access_       *access;                 //!< access database

};

#endif

//! \brief Initializes LES environment
PetscErrorCode InitializeLES(les_ *les);

//! \brief Update Cs Smagorinky coefficient
PetscErrorCode UpdateCs(les_ *les);

//! \brief Update effective viscosity
PetscErrorCode UpdateNut(les_ *les);

//! \brief Compute turbulent statistics based on the les->turbStat flag
PetscErrorCode ComputeTurbulenceStatistics(les_ *les);

PetscErrorCode updateLESStructuralModelCartesianForm(les_ *les);

PetscErrorCode updateLESStructuralModelContravariantForm(les_ *les);