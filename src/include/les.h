//! \file  lesmdl.h
//! \brief LES model header file.

#ifndef _LES_H_
#define _LES_H_

//! \brief struct storing LES model
struct les_
{
    PetscReal     maxCs;
    // Dynamic Smagorisnky
    Vec           lSx, lSy, lSz, lS;
    Vec           lLM, lMM;
    Vec           lLM_old, lMM_old;
    Vec           lCs;

    // kEquation

    // general
    Vec           lNu_t;                   //!< eddy viscosity

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
