//! \file  vents.h
//! \brief vents model header file.

#ifndef VENTS_H
#define VENTS_H


typedef struct
{
    word                face;   //face on which vent lies

    word                dir;   //direction of vent flow. either outlet or inlet. Needed for adjustFluxesVents. will affect all backflow codes as well.

    PetscReal           xBound1;    //lesser bound value for x axis

    PetscReal           xBound2;    //upper bound value for x axis

    PetscReal           yBound1;    //lower bound value for y axis

    PetscReal           yBound2;    //upper bound value for y axis

    PetscReal           zBound1;    //lower bound value for z axis

    PetscReal           zBound2;    // upper bound for z axis

    word               ventBC;	   //boundary condition for vent

    Cmpnts             ventBCVec;  //fixed value velocity components, if needed.

    PetscReal          ventDesiredFlux; //flux determined in hospital drawing search. Use this to set velocity cmpnts of ventBCVEc. M3/s

    //word               wallBC;	   //boundary condition for wall portion around vent, may need to add vector read.

    word               ventBCNut; //Nut BC for vent

    //word               wallBCNut; //Nut BC for wall

    PetscReal          ventBCNutReal; //value for Nut if fixed value at wall

    PetscInt           nCellsVent;

    PetscScalar        ventArea;

    PetscReal          lFluxVent; //flux calculated in simulation M3/s

    //PetscReal          fluxToCellRatio;

    word              ventTBC; //Temp BC type

    PetscScalar       ventTBCVal; //value of fixed temp or temp gradient.

    word              ventCBC; //Temp BC type

    PetscScalar       ventCBCVal; //value of fixed temp or temp gradient.

    PetscScalar       ventCBCStart; //start and end time of intermittent C boundary condition

    PetscScalar       ventCBCEnd;

}ventObject;


struct vents_
{
    PetscInt           nCellsVentsOut;

    PetscInt           nCellsVentsIn;

    //PetscInt           nCellsLeak;

    PetscInt           numberOfVents;

    //word               roomPressure;

    ventObject         **vent;

    //PetscReal          desiredLeakFlux;

    // access database
    access_            *access;

};

#endif

// =============================================================================
// FUNCTIONS
// =============================================================================

PetscErrorCode initializeVents(domain_ *domain);

PetscErrorCode readVentsProperties(vents_ *vents);

PetscErrorCode ventSetAndPrint(vents_ *vents);

PetscErrorCode printVentTemp(vents_ *vents);
