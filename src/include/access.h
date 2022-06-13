//! \file  db.h
//! \brief database object declaration

#ifndef _DATABASE_H_
#define _DATABASE_H_

#include "ceqn.h"

//! \brief access database
struct access_
{
    clock_       *clock;
    simInfo_     *info;
    constants_   *constants;
    flags_       *flags;
    mesh_        *mesh;
    io_          *io;
    ueqn_        *ueqn;
    peqn_        *peqn;
    teqn_        *teqn;
    ceqn_        *ceqn;
    les_         *les;
    farm_        *farm;
    overset_     *os;
    abl_         *abl;
    acquisition_ *acquisition;
    ibm_         *ibm;
    PetscInt     *domainID;
    vents_       *vents;
};

#endif
