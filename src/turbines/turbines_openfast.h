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

//! \brief Get positions of force tower points from OpenFAST
PetscErrorCode getForcePtsTwrOpenFAST(windTurbine *wt);

//! \brief Send velocity to OpenFAST 
PetscErrorCode computeWindVectorsRotorOpenFAST(farm_ *farm);

//! \brief Get force from OpenFAST
PetscErrorCode computeBladeForceOpenFAST(farm_ *farm);

//! \brief Find out which OpenFAST velocity points are controlled by this processor
PetscErrorCode findControlledVelPointsOpenFAST(farm_ *farm);

//! \brief Initialize processor-controlled vel points 
PetscErrorCode initControlledPointsOpenFAST(windTurbine *wt);

//! \brief Get global turbine parameters from OpenFAST
PetscErrorCode getGlobParamsOpenFAST(windTurbine *wt);

//! \brief Get positions of velocity turbine points from OpenFAST
PetscErrorCode getVelPtsBladeOpenFAST(windTurbine *wt);

//! \brief Get positions of force turbine points from OpenFAST
PetscErrorCode getForcePtsBladeOpenFAST(windTurbine *wt);

//! \brief Get positions of velocity tower points from OpenFAST
PetscErrorCode getVelPtsTwrOpenFAST(windTurbine *wt);

#endif

#endif