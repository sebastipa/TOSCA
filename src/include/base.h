//! \file  base.h
//! \brief Base header file including c/c++ libraries and PETSc and HYPRE libs

#ifndef _BASE_H_
#define _BASE_H_

// include base headers
#include <vector>
#include <algorithm>
#include <assert.h>
#include <complex>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iomanip>

// include PETSc headers
#include "petsctime.h"
#include "petscvec.h"
#include "petscdmda.h"
#include "petscksp.h"
#include "petscsnes.h"

// include HYPRE headers
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE_sstruct_ls.h"
#include "HYPRE_IJ_mv.h"

// typedefs
typedef std::string                        word;
typedef std::vector<PetscInt>              labelList;
typedef std::vector<std::vector<PetscInt>> labelListList;
typedef std::complex<PetscReal>            complex;

// constants
#ifndef M_PI
    #define M_PI 3.141592653589793238462643
#endif

struct simInfo_
{
    PetscInt nDomains;
    PetscInt periodic;
};

struct constants_
{
    // dimensionless numbers
    PetscReal        Pr;                         //!< Prantl number
    PetscReal        Sc;                         //!< Schmidt number
    PetscReal        ScT;                         //!< Turbulent Schmidt number

    // physical constants
    PetscReal        nu;                         //!< kinematic viscosity
    PetscReal        rho;                        //!< flow density [Kg/m3]

    // reference temperature
    PetscReal        tRef;                       //!< reference T, required when ABL is not active

};

//! \brief Cell indices
typedef struct
{
   PetscInt i,j,k;
} cellIds;

//! \brief Defines the x, y, z components of a vector
typedef struct
{
    PetscScalar x, y, z;
} Cmpnts;

//! \brief Structure for defining the x, y components of a 2D vector
typedef struct
{
    PetscReal x, y;
} Cpt2D;

//! \brief Defines the components of a symmetric tensor
typedef struct
{
    PetscScalar xx, xy, xz, yy, yz, zz;
} symmTensor;

//! \brief Vector field on the patch (used to store wall models fields)
typedef struct
{
    PetscReal **x;
    PetscReal **y;
    PetscReal **z;
} patchVectorField;

#endif
