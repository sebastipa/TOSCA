//! \file  precursor.h
//! \brief Concurrent precursor header file.

#ifndef _PRECURSOR_H_
#define _PRECURSOR_H_

//! \brief concurrent precursor mapping info
typedef struct
{
    PetscInt kFringeStart, kFringeEnd;             //!< start and end ids of the fringe region in the successr indexing
    PetscInt kStart, kEnd;                         //!< start and end ids of this processor in the k direction in the successor indexing
} mapInfo;

//! \brief concurrent precursor database
struct precursor_
{
    PetscInt      thisProcessorInFringe;      //!< 1: this processors will solve, 0: this processor will idle
    mapInfo       map;                        //!< precursor/successor mapping info

    domain_       *domain;
};

#endif

//! \brief initialize concurrent precursor
PetscErrorCode concurrentPrecursorInitialize(abl_ *abl);

//! \brief Precursor solution flags definition
PetscErrorCode SetSolutionFlagsPrecursor(domain_ *domain);

//! \brief Precursor mesh initialization
PetscErrorCode InitializeMeshPrecursor(abl_ *abl);

//! \brief Precursor boundary conditions initialization function
PetscErrorCode SetBoundaryConditionsPrecursor(mesh_ *mesh);

//! \brief Sets the inflow function on the precursor
PetscErrorCode SetInflowFunctionsPrecursor(mesh_ *mesh);

//! \brief Map fields from successor to precursor
PetscErrorCode MapInitialConditionPrecursor(abl_ *abl);

//! \brief Map vector field from successor to precursor
PetscErrorCode successorPrecursorMapVectorField(abl_ *abl, Vec &Source, Vec &Target, Vec &lTarget);

//! \brief Map scalar tor field from successor to precursor
PetscErrorCode successorPrecursorMapScalarField(abl_ *abl, Vec &Source, Vec &Target, Vec &lTarget);

//! \brief Initialize abl for precursor (x damping layer is not set)
PetscErrorCode ABLInitializePrecursor(domain_ *domain);

//! \brief Solve concurrent precursor
PetscErrorCode concurrentPrecursorSolve(abl_ *abl);
