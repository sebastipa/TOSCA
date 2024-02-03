//! \file  acquisition.h
//! \brief Acquisition object declaration

#ifndef _ACQUISITION_H_
#define _ACQUISITION_H_

#include "acquisition/averageAcquisition.h"
#include "acquisition/probesAcquisition.h"
#include "acquisition/sectionsAcquisition.h"
#include "acquisition/ablAcquisition.h"
#include "acquisition/keAcquisition.h"
#include "acquisition/perturbAcquisition.h"

//! \brief Struct containing all acquisition data structures
struct acquisition_
{
    rakes         *probes;                    //!< probe rakes
    avgFields     *fields;                    //!< average and turbulence fields
    keFields      *keBudFields;               //!< kinetic energy budjets
    sections      *kSections;                 //!< information about the k-sections
    sections      *jSections;                 //!< information about the j-sections
    sections      *iSections;                 //!< information about the i-sections
    udSections    *userSections;              //!< information about the user defined sections
    dataABL       *statisticsABL;             //!< ABL statistics
    data3LM       *LM3;                       //!< 3LM statistics
    perturbFields *perturbABL;                //!< ABL perturbation fields for gravity waves

    // acquisition flags
    PetscInt      isProbesActive;
    PetscInt      isSectionsActive;
    PetscInt      isAverageABLActive;
    PetscInt      isAverage3LMActive;
    PetscInt      isPerturbABLActive;

    access_       *access;                    //!< access database
};

#endif

// GLOBAL ACQUISITION
// ============================================================================================================= //

//! \brief Reads acquisition settings and calls individual iniitialization functions
PetscErrorCode InitializeAcquisition(domain_ *domain);

//! \brief Reads acquisition settings and calls individual iniitialization functions for concurrent precursor
PetscErrorCode InitializeAcquisitionPrecursor(domain_ *domain);

//! \brief Write acquisition data
PetscErrorCode WriteAcquisition(domain_ *domain);

// PROBES ACQUISITION
// ============================================================================================================= //

//! \brief Initializes the array of probe rakes
PetscErrorCode ProbesInitialize(domain_ *domain, PetscInt postProcessing=0);

//! \brief Initialize probe rake file
PetscErrorCode InitRakeFile(probeRake *rake, const char *fieldName);

//! \brief Writes probes to file
PetscErrorCode writeProbes(domain_ *domain);

// SECTIONS ACQUISITION
// ============================================================================================================= //

//! \brief Initializes sections
PetscErrorCode sectionsInitialize(acquisition_ *acquisition);

//! \brief Writes down section data
PetscErrorCode writeSections(acquisition_ *acquisition);

//! \brief Saves i-section vector data
PetscErrorCode iSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt iplane, Vec &V, const char* fieldName);

//! \brief Saves j-section vector data
PetscErrorCode jSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt jplane, Vec &V, const char* fieldName);

//! \brief Saves k-section vector data
PetscErrorCode kSectionSaveVector(mesh_ *mesh, sections *sec, PetscInt kplane, Vec &V, const char* fieldName);

//! \brief Saves i-section scalar data
PetscErrorCode iSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt iplane, Vec &V, const char* fieldName);

//! \brief Saves j-section scalar data
PetscErrorCode jSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt jplane, Vec &V, const char* fieldName);

//! \brief Saves k-section scalar data
PetscErrorCode kSectionSaveScalar(mesh_ *mesh, sections *sec, PetscInt kplane, Vec &V, const char* fieldName);

// FIELDS AVERAGES ACQUISITION
// ============================================================================================================= //

//! \brief Initialize average fields acquisition
PetscErrorCode averageFieldsInitialize(acquisition_ *acquisition);

//! \brief Average fields
PetscErrorCode averageFields(acquisition_ *acquisition);

//! \brief Compute Q criteria for I/O in cartesian form
PetscErrorCode computeQCritIO(acquisition_ *acquisition);

//! \brief Compute coriolis force for I/O in cartesian form
PetscErrorCode computeCoriolisIO(acquisition_ *acquisition);

//! \brief Compute driving pressure force for I/O in cartesian form
PetscErrorCode computeDrivingSourceIO(acquisition_ *acquisition);

//! \brief Compute x damping force for I/O in cartesian form
PetscErrorCode computeXDampingIO(acquisition_ *acquisition);

//! \brief Compute side force for I/O in cartesian form
PetscErrorCode computeCanopyForceIO(acquisition_ *acquisition);

//! \brief Compute velocity divergence field for I/O
PetscErrorCode computeVelocityDivergence(acquisition_ *acquisition);

// MKE BUDGETS ACQUISITION
// ============================================================================================================= //

//! \brief Initializes ke budgets fields acquisition
PetscErrorCode averageKEBudgetsInitialize(acquisition_ *acquisition);

//! \brief Read box array
PetscErrorCode readKeBoxArray(keFields *ke);

//! \brief Search box bounds and define communicators
PetscErrorCode setKeBoundsAndComms(mesh_ *mesh, keFields *ke);

//! \brief Do box cumulation for mesh. energy budgets
PetscErrorCode boxCumulateKEBudgets(acquisition_ *acquisition);

//! \brief Compute MKE budgets
PetscErrorCode averageKEBudgets(acquisition_ *acquisition);

//! \brief Compute MKE budgets in cartesian form
PetscErrorCode averageKEBudgetsCat(acquisition_ *acquisition);

//! \brief Compute MKE budgets in generalized curvilinaer form w/o mesh stretching
PetscErrorCode averageKEBudgetsCont(acquisition_ *acquisition);

// 3LM ACQUISITION
// ============================================================================================================= //

//! \brief Initialize 3LM data structure
PetscErrorCode averaging3LMInitialize(domain_ *domain);

//! \brief Performs 3LM average and write 3LM fields
PetscErrorCode writeAveraging3LM(domain_ *domain);

//! \brief Write 3LM points to file inside postProcessing/3LM/points
PetscErrorCode write3LMPoints(acquisition_ *acquisition);

//! \brief Write 3LM fields to file inside postProcessing/3LM/
PetscErrorCode write3LMFields(acquisition_ *acquisition);

//! \brief Inizializes the closest cells and the owner of the 3LM mesh points
PetscErrorCode findAvgLineIds(acquisition_ *acquisition);

//! \brief Reads velocity and pressure averages (they don't use auxiliary fields like TBL, LBL, IBL)
PetscErrorCode read3LMFields(acquisition_ *acquisition);

// ABL PERTURBATION ACQUISITION
// ============================================================================================================= //

//! \brief Initializes ABL perturbation fields w.r.t. given reference
PetscErrorCode perturbationABLInitialize(acquisition_ *acquisition);

//! \brief Performs perturbation averaging and writes to memory
PetscErrorCode averagePerturbationABL(acquisition_ *acquisition);

// ABL PLANAR AVERAGING ACQUISITION
// ============================================================================================================= //

//! \brief Initialize ABL data structure
PetscErrorCode averagingABLInitialize(domain_ *domain);

//! \brief Write spatial averaged ABL statistics
PetscErrorCode writeAveragingABL(domain_ *domain);
