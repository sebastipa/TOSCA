//! \file  pores.h
//! \brief pores model header file. Momentum sink

#ifndef PORES_H
#define PORES_H


typedef struct
{

    PetscReal           xBound1;    //lesser bound value for x axis

    PetscReal           xBound2;    //upper bound value for x axis

    PetscReal           yBound1;    //lower bound value for y axis

    PetscReal           yBound2;    //upper bound value for y axis

    PetscReal           zBound1;    //lower bound value for z axis

    PetscReal           zBound2;    // upper bound for z axis


}porousZone;


struct pores_
{
    PetscInt           numberOfPores;

    porousZone         **poreZone;

    PetscReal          porosity;

    PetscReal          alpha;

    PetscReal          C2;

    PetscReal          filtrationEfficiency;
    // access database
    access_            *access;

};

#endif

// =============================================================================
// FUNCTIONS
// =============================================================================

PetscErrorCode InitializePores(domain_ *domain);

PetscErrorCode readPoresProperties(pores_ *pores);

PetscErrorCode poreSetAndPrint(pores_ *pores);
