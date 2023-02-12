//! \file  initialField.h
//! \brief initial field header file.

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
