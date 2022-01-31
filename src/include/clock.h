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
<<<<<<< HEAD
    PetscReal        startDt;                   //!< time step
=======
    PetscReal        dtOld;                     //!< old time step
>>>>>>> a1ad32b (added alphaOptimized fringe controller and timeStep increment control)
    PetscReal        cfl;                       //!< cfl number

    PetscInt         it;                        //!< current iteration value
    PetscInt         itStart;                   //!< start iteration value (zero)

    PetscInt         timePrecision;             //!< time write precision
};

#endif

PetscErrorCode adjustTimeStep(domain_ *domain);
