#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_openfast.h

#ifndef TRB_OF_H
#define TRB_OF_H

#if USE_OPENFAST

#include "OpenFAST.H"

//! \brief Initialize OpenFAST coupling
PetscErrorCode initOpenFAST(farm_ *farm);

//! \brief Initial condition for OpenFAST solution
PetscErrorCode stepZeroOpenFAST(farm_ *farm);

//! \brief Time step OpenFAST solution
PetscErrorCode stepOpenFAST(farm_ *farm);

//! \brief Compute and end blade velocity to OpenFAST 
PetscErrorCode computeWindVectorsRotorOpenFAST(farm_ *farm);

//! \brief Compute and end tower velocity to OpenFAST
PetscErrorCode computeWindVectorsTowerOpenFAST(farm_ *farm);

//! \brief Get blade force from OpenFAST
PetscErrorCode computeBladeForceOpenFAST(farm_ *farm);

//! \brief Get tower from OpenFAST
PetscErrorCode computeTowerForceOpenFAST(farm_ *farm);

//! \brief Find out which OpenFAST points are controlled by which processor
PetscErrorCode findControlledPointsRotorOpenFAST(farm_ *farm);

//! \brief Find out which OpenFAST points are controlled by which processor
PetscErrorCode findControlledPointsTowerOpenFAST(farm_ *farm);

//! \brief Update turbine points and parameters from OpenFAST after a step
PetscErrorCode getDataFromOpenFAST(farm_ *farm);

    //! \brief Get global turbine parameters from OpenFAST
    PetscErrorCode getGlobParamsOpenFAST(windTurbine *wt, const word actuatorModel);

    //! \brief Get positions of velocity turbine points from OpenFAST
    PetscErrorCode getVelPtsBladeOpenFAST(windTurbine *wt);

    //! \brief Get positions of force turbine points from OpenFAST
    PetscErrorCode getForcePtsBladeOpenFAST(windTurbine *wt);

    //! \brief Get positions of velocity tower points from OpenFAST
    PetscErrorCode getVelPtsTwrOpenFAST(windTurbine *wt);

    //! \brief Get positions of force tower points from OpenFAST
    PetscErrorCode getForcePtsTwrOpenFAST(windTurbine *wt);

//! \brief Translate OpenFAST force points displacements into TOSCA actuator point displacements 
PetscErrorCode mapOFDisplToActPts(farm_ *farm);

//! \brief Translate OpenFAST forces into TOSCA forces depending on the actutor model
PetscErrorCode mapOFForcesToActPts(farm_ *farm);

#endif

#endif