#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_phys.h
//! \brief Turbine physics module and control (substituted by OpenFAST if the latter is used)

#ifndef TRB_PHYS_H
#define TRB_PHYS_H

//! \brief Solve rotor dynamics and compute filtered rot speed
PetscErrorCode computeRotSpeed(farm_ *farm);

//! \brief Rotate blades (only for ALM)
PetscErrorCode rotateBlades(windTurbine *wt, PetscReal angle, PetscInt updateAzimuth);

//! \brief Compute generator torque with 5-regions control system model
PetscErrorCode controlGenSpeed(farm_ *farm);

//! \brief Compute blade pitch using the PID control
PetscErrorCode controlBldPitch(farm_ *farm);

//! \brief Compute nacelle yaw
PetscErrorCode controlNacYaw(farm_ *farm);

//! \brief Compute wind turbine control based on wind farm controller
PetscErrorCode windFarmControl(farm_ *farm);

#endif