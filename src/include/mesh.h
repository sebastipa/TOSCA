//! \file  mesh.h
//! \brief Mesh header file.

#ifndef _MESH_H_
#define _MESH_H_

#include "wallmodel.h"
#include "boundary.h"

//! \brief Structure defining the domain bounding box (easy access to xmin, xmax, ymin, ymax, zmin, zmax)
typedef struct
{
    PetscReal      xmin,                       // domain bounding coordinates in meters
                   xmax,
                   ymin,
                   ymax,
                   zmin,
                   zmax;
    PetscReal      Lx,                         // domain length in meters in x,y,z directions
                   Ly,
                   Lz;
} boundingBox;

//! \brief mesh structure storing mesh and curvilinear coordinates
struct mesh_
{
    word           meshName;                   //!< name of the domain
    word           meshFileType;               //!< type of mesh file

    boundingBox    bounds;                     //!< domain extensions and lengths information
    PetscReal      grndLevel;                  //!< the ground level height - would be = bounds->zmin for normal simulation, but could change if using IBM
      
    PetscInt       IM, JM, KM;                 //!< ncells in the GCC directions

    // distributed arrays
    DM             da;                         //!< data structure for scalars (include the grid geometry informaion, to obtain the mesh information, use DMDAGetCoordinates)
    DM             fda;                        //!< data Structure for vectors
    DM             sda;                        //!< data Structure for symmetric tensors
    DMDALocalInfo  info;                       //!< data struct that contains information about a the Distributed Array (essentially mesh information)

    // metrics
    Vec            Cent, lCent;                //!< cell centers coordinates in the cartesian frame
    Vec            lCsi, lEta, lZet;
    Vec            lAj;                        //!< jacobian evaluated at cell centers (1 / cell volume)
    Vec            lICsi, lIEta, lIZet;
    Vec            lIAj;                       //!< jacobian evaluated at the i-faces
    Vec            lJCsi, lJEta, lJZet;
    Vec            lJAj;                       //!< jacobian evaluated at the j-faces
    Vec            lKCsi, lKEta, lKZet;
    Vec            lKAj;                       //!< jacobian evaluated at the j-faces

    // IBM markup
    Vec            Nvert, Nvert_o;             //!< solid body field for IBM
    Vec            lNvert, lNvert_o;

    // vent and porous zone marker
    Vec            ventMarkers;

    // periodic connectivity
    PetscInt       i_periodic, ii_periodic;
    PetscInt       j_periodic, jj_periodic;
    PetscInt       k_periodic, kk_periodic;

    // boundary conditions
    scalarBC       boundaryNut;
    scalarBC       boundaryT;
    vectorBC       boundaryU;

    // special inflow boundary conditions
    inletFunctions inletF;

    // limiters and numerical fields
    Vec            fluxLimiter;

    // this mesh communicator
    MPI_Comm       MESH_COMM;

    // access database
    access_        *access;

};

#endif

// =============================================================================
// FUNCTIONS
// =============================================================================

//! \brief Initialize mesh
PetscErrorCode InitializeMesh(mesh_ *mesh);

//! \brief Set distributed arrays
PetscErrorCode SetDistributedArrays(mesh_ *mesh);

//! \brief Set curvilinear coordinates metrics
PetscErrorCode SetMeshMetrics(mesh_ *mesh);

//! \brief Set bounding box
PetscErrorCode SetBoundingBox(mesh_ *mesh);

//! \brief Find the cell center for the ghost nodes
PetscErrorCode ghostnodesCellcenter(mesh_ *mesh);

//! \brief Deform the mesh according to prescribed BL Disp
PetscErrorCode DeformMeshBasedOnBLDisp(mesh_ *mesh);
