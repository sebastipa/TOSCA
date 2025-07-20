//! \file  lesScalar.h
//! \brief LES model header file for scalar quantities.

PetscErrorCode InitializeLESScalar(les_ *les);

PetscErrorCode UpdateCsk (les_ *les);

PetscErrorCode UpdatekT(les_ *les);

