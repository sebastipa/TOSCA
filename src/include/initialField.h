//! \file  initialField.h
//! \brief initial field header file.

#ifndef INITIALFIELD_H
#define INITIALFIELD_H

struct ModeData
{
    PetscReal   kappaM; 
    PetscReal   khatx, khaty, khatz;
    PetscReal   amp;
    PetscReal   phi;
    PetscReal   ktx, kty, ktz;
    PetscReal   sx, sy, sz;
};

#endif

// =============================================================================
// FUNCTIONS
// =============================================================================

//!< \brief set the initial internal fields
PetscErrorCode SetInitialField(domain_ *domain);

//!< \brief set the initial internal fields for the concurrent precursor simulation
PetscErrorCode SetInitialFieldPrecursor(abl_ *abl);

//!< \brief set the initial internal contravariant and cartesian velocity field
PetscErrorCode SetInitialFieldU(ueqn_ *ueqn);

//!< \brief set the initial temperature field
PetscErrorCode SetInitialFieldT(teqn_ *teqn);

//!< \brief set the initial pressure field
PetscErrorCode SetInitialFieldP(peqn_ *peqn);

//!< \brief set the initial variables for LES
PetscErrorCode SetInitialFieldLES(les_ *les);

//!< \brief set uniform value for the cartesian velocity (add perturbations if applicable)
PetscErrorCode SetUniformFieldU(ueqn_ *ueqn, Cmpnts &uRef, PetscInt &addPerturbations);

//!< \brief set initial ABL flow U
PetscErrorCode SetABLInitialFlowU(ueqn_ *ueqn);

PetscErrorCode SetABLInitialFlowUZilitinkevich(ueqn_ *ueqn);

//!< \brief set initial field Taylor green vortex problem
PetscErrorCode SetTaylorGreenFieldU(ueqn_ *ueqn, PetscReal &u0, PetscReal &freq);

//!< \brief set initial field by generating synthetic turbulence (homogenous isotropic turbulence)
PetscErrorCode SetHITFieldU(ueqn_ *ueqn,PetscInt  nmodes, PetscReal kMin);

//!< \brief set initial field ABC flow (also called the Arnold–Beltrami–Childress flow) is a well-known three-dimensional, steady, incompressible velocity field
PetscErrorCode SetABCFlow(ueqn_ *ueqn, PetscReal &u0, PetscInt &k0);

//!< \brief set the internal field as the spreaded inlet flow condition
PetscErrorCode SpreadInletFlowU(ueqn_ *ueqn);

//!< \brief set uniform value for the temperature
PetscErrorCode SetUniformFieldT(teqn_ *teqn, PetscReal &tRef);

//!< \brief set linear profile for the temperature
PetscErrorCode SetLinearFieldT(teqn_ *teqn, PetscReal &tRef, PetscReal &tLapse);

//!< \brief set initial ABL flow T
PetscErrorCode SetABLInitialFlowT(teqn_ *teqn);

//!< \brief set the internal field as the spreaded inlet flow condition
PetscErrorCode SpreadInletFlowT(teqn_ *teqn);

//!< \brief set the internal pressure field as the TGV flow
PetscErrorCode SetTaylorGreenFieldP(peqn_ *peqn);
