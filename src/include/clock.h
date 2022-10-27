//! \file  clock.h
//! \brief Simulation time struct header file

#ifndef _CLOCK_H_
#define _CLOCK_H_

struct clock_
{
    PetscReal        time;                      //!< current time value
    PetscReal        startTime;                 //!< simulation start time
    PetscReal        endTime;                   //!< simulation end time

    word             startFrom;                 //!< choose if start from startTime value or latestTime (startTime is ignored)

    PetscReal        dt;                        //!< time step
    PetscReal        startDt;                   //!< time step
    PetscReal        dtOld;                     //!< old time step
    PetscReal        cfl;                       //!< cfl number
    PetscReal        dxMin;                     //!< min cell side size
    PetscReal        acquisitionDt;             //!< uniform dt due to acquistion if applied from the start

    PetscInt         it;                        //!< current iteration value
    PetscInt         itStart;                   //!< start iteration value (zero)

    PetscInt         timePrecision;             //!< time write precision
};

#endif

PetscErrorCode adjustTimeStep(domain_ *domain);
