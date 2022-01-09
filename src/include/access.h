//! \file  db.h
//! \brief database object declaration

#ifndef _DATABASE_H_
#define _DATABASE_H_

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
    les_         *les;
    farm_        *farm;
    overset_     *os;
    abl_         *abl;
    acquisition_ *acquisition;
    ibm_         *ibm;
};

#endif
