//! \file domain.h
//! \brief Domain data structure definition

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "objects.h"
#include "clock.h"
#include "vents.h"
#include "access.h"
#include "flags.h"
#include "mesh.h"
#include "overset.h"
#include "acquisition.h"
#include "ibm.h"
#include "ueqn.h"
#include "peqn.h"
#include "teqn.h"
#include "les.h"
#include "precursor.h"
#include "turbines.h"
#include "abl.h"


//! \brief Domain data structure definition
struct domain_
{
    simInfo_      info;
    constants_    constants;

    PetscInt      domainID;
    clock_        *clock;            //!< time info

    mesh_         *mesh;             //!< mesh data structure

    io_           *io;               //!< io settings data structure

    ueqn_         *ueqn;             //!< momentum equation data structure
    peqn_         *peqn;             //!< pressure equation data structure
    teqn_         *teqn;             //!< potential temperature transport equation
    les_          *les;              //!< LES model data structure

    overset_      *os;               //!< overset data structure

    farm_         *farm;             //!< wind farm data structure

    ibm_          *ibm;              //!< IBM data structure

    vents_        *vents;           //!< vents data structure

    abl_          *abl;              //!< atmospheric boundary layer data structure

    acquisition_  *acquisition;      //!< acquisition data structure

    flags_        flags;             //!< solution flags data structure
    access_       access;            //!< access database data structure
};

#endif
