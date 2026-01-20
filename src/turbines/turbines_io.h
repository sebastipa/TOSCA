#include "../include/base.h"
#include "../include/domain.h"
#include "../include/inline.h"
#include "../include/turbines.h"  

//! \file  turbines_io.h
//! \brief Turbine IO routines

#ifndef TRB_IO_H
#define TRB_IO_H

//! \brief Write output if applicable
PetscErrorCode windTurbinesWrite(farm_ *farm);

//! \brief Write checkpoint file
PetscErrorCode windTurbinesWriteCheckpoint(farm_ *farm);

//! \brief Read checkpoint file and prepare wind turbines
PetscErrorCode windTurbinesReadCheckpoint(farm_ *farm);

//! \brief Print wind farm information
PetscErrorCode printFarmProperties(farm_ *farm);

//! \brief Write the wind farm AD mesh to ucd file
PetscErrorCode writeFarmADMesh(farm_ *wf);

//! \brief Write the wind farm AL mesh to ucd file
PetscErrorCode writeFarmALMesh(farm_ *farm);

//! \brief Write the wind farm tower mesh to ucd file
PetscErrorCode writeFarmTwrMesh(farm_ *wf);

//! \brief Read actuator model parameters and allocate memory 
PetscErrorCode readActuatorModelParameters(const word actuatorModel, windTurbine *wt, const word meshName); 

//! \brief Initializes the ADM reading from files and allocating memory
PetscErrorCode readADM(windTurbine *wt, const word meshName);

//! \brief Creates the ADM mesh 
PetscErrorCode meshADM(windTurbine *wt);

//! \brief Initialize the UADM Model
PetscErrorCode readUADM(windTurbine *wt, const word meshName);

//! \brief Creates the UADM mesh
PetscErrorCode meshUADM(windTurbine *wt);

//! \brief Initializes the ALM reading from files and allocating memory
PetscErrorCode readALM(windTurbine *wt, const word meshName);

//! \brief Creates the ALM mesh
PetscErrorCode meshALM(windTurbine *wt);

//! \brief Initializes the AFM reading from files and allocating memory
PetscErrorCode readAFM(windTurbine *wt, const word meshName);

//! \brief Creates the AFM mesh
PetscErrorCode meshAFM(windTurbine *wt);

//! \brief Initialize the upstream sample points data structure
PetscErrorCode initSamplePoints(windTurbine *wt, const word meshName);

//! \brief Initializes the tower model
PetscErrorCode initTwrModel(windTurbine *wt, Cmpnts &base);

//! \brief Initializes the nacelle model
PetscErrorCode initNacModel(windTurbine *wt, Cmpnts &base);

//! \brief Read the farm_Properties file and fill the structs
PetscErrorCode readFarmProperties(farm_ *farm);

//! \brief Read the wind turbines inside the farm_Properties file
PetscErrorCode readTurbineArray(farm_ *wf);

//! \brief Fill the windTurbine struct reading in the file named as the wind turbine type
PetscErrorCode readTurbineProperties(windTurbine *wt, const char *dictName, const word meshName, const word modelName, Cmpnts &base);

//! \brief Reads the names of the airfoil used in the turbine (given in the airfoils subdict inside the file named as the wind turbine type)
PetscErrorCode readAirfoilProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the blades aero properties used in the turbine (given in the bladeData subdict inside the file named as the wind turbine type)
PetscErrorCode readBladeProperties(windTurbine *wt, const char *dictName, const PetscInt readThickness);

//! \brief Reads the turbine Ct curve used in the turbine (given in the CtTable subdict inside the file named as the wind turbine type)
PetscErrorCode readCtTable(windTurbine *wt, const char *dictName);

//! \brief Reads the turbine blade pitch curve used in the turbine (given in the pitchTable subdict inside the file /turbines/control/bladePitchCurve)
PetscErrorCode readPitchTable(windTurbine *wt, const char *dictName);

//! \brief Reads the turbine rpm curve used in the turbine (given in the rpmTable subdict inside the file /turbines/control/rotorRpmCurve)
PetscErrorCode readRpmTable(windTurbine *wt, const char *dictName);

//! \brief Reads the tower properties used in the turbine (given in the towerData subdict inside the file named as the wind turbine)
PetscErrorCode readTowerProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the nacelle properties used in the turbine (given in the nacelleData subdict inside the file named as the wind turbine)
PetscErrorCode readNacelleProperties(windTurbine *wt, const char *dictName);

//! \brief Reads the 2D airfoil tables given in the file named as the airfoil
PetscErrorCode readAirfoilTable(foilInfo *af, const char *tableName);

//! \brief Reads generator torque controller parameters
PetscErrorCode readGenControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Reads blade pitch controller parameters
PetscErrorCode readPitchControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Reads nacelle yaw controller parameters
PetscErrorCode readYawControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Read wind farm controller table (1 header and time - value list)
PetscErrorCode readWindFarmControlTable(windTurbine *wt);

//! \brief Read discrete induction controller parameters
PetscErrorCode readDipcControllerParameters(windTurbine *wt, const char *dictName, const char *meshName);

//! \brief Check that the mesh is resolved around the turbine
PetscErrorCode checkTurbineMesh(farm_ *farm);

#endif