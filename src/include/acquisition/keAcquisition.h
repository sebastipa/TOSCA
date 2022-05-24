//! \file  keAcquisition.h
//! \brief mean kinetic energy budget object declaration

#ifndef _KEACQUISITION_H
#define _KEACQUISITION_H

//! \brief Struct defining the ke fields
struct keBox
{
    word          name;
    Cmpnts        center;
    PetscReal     sizeX;
    PetscReal     sizeY;
    PetscReal     sizeZ;

    PetscInt      minKFace, maxKFace,
                  minJFace, maxJFace,
                  minIFace, maxIFace;

    PetscReal     avgDum, avgDup, avgDpm, avgDpp,
                  avgFI, avgFJ, avgFK, avgEps, avgPf,
                  avgPtheta, avgPinf;

    // communication color
    MPI_Comm          KEBOX_COMM;   //!< communicator for this box
    PetscMPIInt       writerRank;   //!< label of master rank of the KEBOX_COMM communicator in the MPI_COMM_WORLD rank list
    PetscMPIInt        nProcsBox;   //!< size of the KEBOX_COMM communicator

};

//! \brief Struct defining the ke fields
struct keFields
{
    // flags
    PetscInt      cartesian;       //!< perform the budget in cartesian components (GCC otherwise)
    PetscInt      debug;

    PetscInt      nBox;            //!< number of boxes in the array
    keBox         *box;          //!< array of boxes

    // fields
    Vec           Error;           //!< error on the equation
    Vec           D;
    Vec           F;
    Vec           Pinf;            //!< mean background pressure power
    Vec           Pf;              //!< wind turbine power
    Vec           Ptheta;          //!< kinetic to potential energy conversion
    Vec           Eps;             //!< turbulent dissipation
    Vec           lEm;             //!< mechanical energy
    Vec           lDum;            //!< mean kinetic energy advection
    Vec           lDup;            //!< turbulent kinetic energy advection
    Vec           lF;              //!< turbulent fluxes

    // working fields
    Vec           lVmTauSGS;        //!< advection of residual stresses
    Vec           lDpm;             //!< advecton of mean pressure work
    Vec           lDpp;             //!< advection of mean pressure fluct
    Vec           lUm;              //!< average of cartesian velocity

    // working fields at cell faces
    Vec           lVmVpVp;           //!<
    Vec           lVpVpVp;           //!<
    Vec           lVDm;              //!<
    Vec           lVmPmG;            //!<
    Vec           lVpPpG;            //!<
    Vec           lDdVm;             //!<
    Vec           lVmCsi;            //!<
    Vec           lVmEta;            //!<
    Vec           lVmZet;            //!<
    Vec           lVm;               //!<
    Vec           lPm;               //!<
    Vec           lVpVpEta;          //!<
    Vec           lVpVpCsi;          //!<
    Vec           lVpVpZet;          //!<
    Vec           lVpVp;             //!<

};

#endif
