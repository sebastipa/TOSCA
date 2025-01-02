//! \file  vents.h
//! \brief vents model header file.

#ifndef VENTS_H
#define VENTS_H


typedef struct
{
    PetscScalar       smBCVal; //value of fixed sm gradient..

}ventSMObject;

typedef struct
{
    // type of inlet function
    PetscInt      typeU;

    //type 6: Synthetic turbulence inlet using random fourier series model.
    PetscInt      FSumNum;                    //!< number of fourier series to include in summations
    PetscInt      iterTKE;                    //!< number of iteration for Energy spectrum adjustment. should be <= 10.
    PetscReal     Urms;                       //!< rms found from TI and desired flux
    PetscReal     kolLScale;                  //!< kolomagrov length scale in m
    PetscReal     intLScale;                //!< length scale in m of largest eddy
    PetscReal     dkn;                        //!< calculated wave vector interval, based on log scale
    PetscReal     kMax;                       //!< largerst wave number to include in driving energy spectrum
    PetscReal     kMin;                       //!< smallest wave number to include in driving energy spectrum
    PetscReal     kKol;                       //!< wave number assoicated with smallest eddies
    PetscReal     kEng;                       //!< wave number assoicated with most energetic eddies (peak of driving energy spectrum)
    Cmpnts        meanU;                      //!< mean velocity vector
    Cmpnts           *kn;                         //!< randomized wave vector
    Cmpnts           *Gn;                         //!< randomized direction unit vector
    PetscReal           *phaseN;                     //!< randomized phase angle
    PetscReal           *uMagN;                      //!< randomized velocity magnitudes
    PetscReal           *knMag;                      //!< wave number magnitudes
    PetscReal           *Ek;                         //!< Driving energy spectra values
    word                genType;                //!< shape function for turbulence. isotropic, transverse, or streamwise.

    PetscReal     TI;                       //!< user defined turbulence intensity
    PetscReal     desiredFlux;                       //!< user defined desired vent flux in m3/s

}ventInletFunction;

typedef struct
{
    word                face;   //face on which vent lies. kLeft, kRight, etc.
    word                dir;   //direction of vent flow. either outlet or inlet.
    word                shape;   //circle or rectangle vent

    PetscReal           xBound1;    //lesser bound value for x axis; rectangle only.
    PetscReal           xBound2;    //upper bound value for x axis
    PetscReal           yBound1;    //lower bound value for y axis
    PetscReal           yBound2;    //upper bound value for y axis
    PetscReal           zBound1;    //lower bound value for z axis
    PetscReal           zBound2;    // upper bound for z axis

    PetscReal           dia;    //diameter of circle vent
    Cmpnts              center;    // center of circle vent in x,y,z

    word               ventBC;	   //boundary condition for vent velocity
    ventInletFunction      *inletF;    //inflow data for vent object if turbulent inflow
    Cmpnts             ventBCVec;  //fixed value velocity components, calculated based on ventDesiredFlux
    PetscReal          ventDesiredFlux; //m3/s of desired vent flux
    word               ventBCNut; //Nut BC for vent
    PetscReal          ventBCNutReal; //value for Nut if fixedValue at wall
    word              ventTBC; //Temp BC type
    PetscScalar       ventTBCVal; //value of fixed temp or temp gradient.
    ventSMObject      **ventSMBC;    //inflow data for vent object
    word              smBC; //sm BC types

    PetscInt           nCellsVent; //number of cells for each vent.
    PetscScalar        ventArea; //area of each vent (m2).
    PetscScalar        ventCentX; //x-coor of vent center
    PetscScalar        ventCentY; //y-coor of vent center
    PetscScalar        ventCentZ; //z-coor of vent center

    PetscReal          lFluxVent; //flux calculated in simulation m3/s. used for adjustFluxesVents


}ventObject;


struct vents_
{
    PetscInt           nCellsVentsOut; //total number of cells for all outlet vents
    PetscInt           nCellsVentsIn; //total number of cells for all inlet vents
    //PetscInt           nCellsLeak;

    PetscInt           numberOfVents;
    ventObject         **vent;

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
