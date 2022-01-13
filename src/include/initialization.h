//! \file  initialization.h
//! \brief Contains simulation initialization function headers

// =============================================================================
// FUNCTIONS
// =============================================================================

//! \brief Print OkWind logo
PetscErrorCode PrintOkWindLogo();

//! \brief Print number of processors
PetscErrorCode PrintNumberOfProcs();

//! \brief Set simulation flags
PetscErrorCode SetSimulationFlags(flags_ *flags);

//! \brief Set simulation global info
PetscErrorCode SetSimulationInfo(simInfo_ *info);

//! \brief Set number of domains (reads from OversetInput.dat file if overset is active)
PetscErrorCode SetDomainsAndAllocate(domain_ **domain, flags_ *flags, simInfo_ *info);

//! \brief Read time controls
PetscErrorCode ReadTimeControls(clock_ *clock);

//! \brief Read physical constants
PetscErrorCode ReadPhysicalConstants(domain_ *domain);

//! \brief Allocate memory for domain pointers
PetscErrorCode SetDomainMemory(domain_ *domain);

//! \brief Set access database pointers
PetscErrorCode SetAccessPointers(domain_ *domain);

//! \brief Initialize the simulation parameters
PetscErrorCode simulationInitialize(domain_ **domainAddr, clock_ *clock, simInfo_ *info, flags_ *flags);

//!< \brief set the initial internal fields
PetscErrorCode SetInitialField(domain_ *domain);

//!< \brief set the initial internal contravariant and cartesian velocity field
PetscErrorCode SetInitialField(ueqn_ *ueqn);
