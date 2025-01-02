//! \file  sm.h
//! \brief T equation solution header file.

#ifndef _SMOBJ_H_
#define _SMOBJ_H_

#include <petscts.h>
#include <petscdm.h>
#include <petscdmda.h>
#include <petscdmcomposite.h>

typedef struct
{

    Vec     weight;                //!< weight of particles of size absc
    Vec     absc;                 //!< particle sizes for method of moments
    Vec     tauP;                 //!< responce time of particles with size absc

}WandA;

//! \brief struct storing scalar moment equation
typedef struct
{
    //boundary conditions
    word      iLeft, iRight,  //!< type of boundary condition
              jLeft, jRight,
              kLeft, kRight;
    PetscReal iLval, iRval ,  //!< value (defined if BC prefix is 'fixed' only )
              jLval, jRval ,
              kLval, kRval ;

    // scalar moment variables
    PC            pc;
    Vec           Rhs;
    Vec           Rhs_o;
    Vec           smTmp;                   //!< temporary solution
    Vec           smVal, lsmVal,           //!< global and local sm values of current time step
                  sm_o, lsm_o;             //!< global and local sm values of old time step
    Vec           lDivSM, lViscSM;         //!< viscou and divergence scalar moment equation fluxes
    Vec           lSed, Sed;               //!< sedimatation values at cell faces and centers
    Vec           sedCent, lSedCent;
    Vec           lDev, Dev, DevCell;      //!< deviation values at cell faces and centers
    Vec           devCent, lDevCent;
    Vec           coagSource;              //!< coagulation scalar moment source, only need at cell center
    Vec           lDep, Dep;               //!< deposition scalar moment source, occurs at boundary faces, but external boundary is always 0.


    word          ddtScheme;                  //!< time derivative scheme

    //subaccess to back track to SMOBJ when Needed
    access_       *access;                     //!< access database

}sm_;


struct SMObj_
{
   sm_         **sm;             //!< substracture for scalar moment 0 to 6

   WandA       **weightAbsc;     //!sub structure for weights and abscissi

   Vec          Eps;             //!< turbulent dissipation vector for coagulation source term

   PetscInt     dissWeight;      //!< weight for time average in dissipation calc.

   PetscReal    OGConc;         //!< original particle concentration in #/m3.

   PetscReal    monoDFrac;      //!< monomer size for coagulation in meters.

   PetscReal    rhoPart;        //!< partcile density.

   Vec          quant;
   Vec          Dq;
   Vec          probI;          //!< infection probability estimated based on qunata (quant) and infectious dose (Dq)
   Vec          ExCount, DepCount;

   // initial field type readField or uniform
   word          initFieldType;



};


#endif

//! \brief Initializes scalar moment weight and abscissi vectors
PetscErrorCode InitializeSMObject(SMObj_ *smObject);

//! \brief Initializes scalar moment environments
PetscErrorCode InitializeSM(sm_ *sm);

//! \brief Solve SM equation
PetscErrorCode SolveSM(SMObj_ *smObject);

//! \brief RHS of the potential scalar moment transport equation visc and div terms
PetscErrorCode FormSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii);

//! \brief Computed RHS of scalar moment equation using current lsm (updates Rhs)
PetscErrorCode FormExplicitRhsSM(sm_ *sm, PetscInt ii);

//! \find weights and abscissi for moment equations to add source terms later on
PetscErrorCode quickUpdateWeightsAndAbscissi(SMObj_ *smObject);

//! \find dissipation rate to include in coag source terms
PetscErrorCode DissipationCalc(SMObj_ *smObject);

//! \find flux term dur to sedimentation of particles
PetscErrorCode sedFluxSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii);

//! \find flux term dur to deviation of particles
PetscErrorCode devFluxSM(sm_ *sm, Vec &Rhs, PetscReal scale, PetscInt ii);

//! \esimtate infection probability based on particle concentration (SM0)
PetscErrorCode infectProb(SMObj_ *smObject);

//! \find source term due to coagulation of particles
PetscErrorCode formCoagSourceExp(sm_ *sm, PetscInt ii);

//! \find source term due to deposition of particles
PetscErrorCode formDepSourceExp(sm_ *sm, PetscInt ii);

//! \apply source term due to coagulation of particles
PetscErrorCode sourceSMDep(sm_ *sm, Vec &Rhs, PetscReal scale);

//! \apply source term due to deposition of particles
PetscErrorCode sourceSMCoag(sm_ *sm, Vec &Rhs, PetscReal scale);

//! \track particles deposited or exhausted
PetscErrorCode partCount(SMObj_ *smObject);
