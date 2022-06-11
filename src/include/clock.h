//! \file  clock.h
//! \brief Simulation time struct header file

#ifndef _CLOCK_H_
#define _CLOCK_H_

struct clock_
{
    PetscReal        time;                      //!< current time value
    PetscReal        startTime;                 //!< simulation start time
    PetscReal        endTime;                   //!< simulation end time

    PetscReal        dt;                        //!< time step
    PetscReal        startDt;                   //!< time step
    PetscReal        dtOld;                     //!< old time step
    PetscReal        cfl;                       //!< cfl number
    PetscReal        dxMin;                     //!< min cell side size

    PetscInt         it;                        //!< current iteration value
    PetscInt         itStart;                   //!< start iteration value (zero)

    PetscInt         timePrecision;             //!< time write precision
};

#endif

PetscErrorCode adjustTimeStep(domain_ *domain);
