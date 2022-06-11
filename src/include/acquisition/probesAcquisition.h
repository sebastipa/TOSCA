//! \file  probesAcquisition.h
//! \brief probes and rakes object declaration

#ifndef _PROBESACQUISITION_H_
#define _PROBESACQUISITION_H_

//! \brief Struct defining a probe rake
struct probeRake
{
    word          rakeName;                   //!< name of the probes rake
    word          timeName;                   //!< name of the time directory where rake data are written
    PetscInt      probesNumber;               //!< number of probes in this rake
    PetscReal     timeStart;                  //!< start time of acquisition system
    PetscReal     timeInterval;               //!< acquisition time interval (overrides simulation time step if smaller)
    PetscInt      Uflag;                      //!< flag telling if velocity field must be acquired
    PetscInt      Tflag;                      //!< flag telling if fempearature field must be acquired
    PetscInt      Pflag;                      //!< flag telling if pressure field must be acquired

    PetscReal     **locations;                //!< array of [n, 3] where n is the number of probes in the rake and 3 are the probe coords
    PetscInt      **cells;                    //!< array of [n, 3] where n is the number of probes in the rake and 3 are cell indices of the computational pt closes to the probe
    PetscInt      *owner;                     //!< array of [n] where n is the number of probes in the rake and every element is the rank owning the probe

    PetscInt      *domainID;                  //!< number of the domain which controls this probe
    PetscInt      thisRakeControlled;         //!< this processor controls this rake (totally or partially)
    MPI_Comm      RAKE_COMM;                  //!< communicator for this probe rake
};

//! \brief Struct defining all rakes in the mesh
struct rakes
{
    PetscInt      nRakes;                     //!< number of probe rakes
    probeRake     *rakes;                     //!< array storing all probe rakes in the domain
};

#endif
