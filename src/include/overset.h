//! \file  overset.h
//! \brief Definitions for overset objects and functions for two–way coupling

#ifndef _OVERSET_H_
#define _OVERSET_H_

#include <unordered_map>

struct BinIndex {
    int ix, iy, iz;
    bool operator==(const BinIndex& other) const {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

namespace std {
    template <>
    struct hash<BinIndex> {
        std::size_t operator()(const BinIndex& b) const {
            return ((std::hash<int>()(b.ix) ^ (std::hash<int>()(b.iy) << 1)) >> 1) ^ (std::hash<int>()(b.iz) << 1);
        }
    };
}

//! \brief Overset acceptor cell (used for both overset and base mesh acceptor cells)
typedef struct {
    PetscInt      indi, indj, indk;  //!< cell indices
    PetscReal     coorx, coory, coorz; //!< cell center coordinates
    PetscMPIInt   rank;              //!< owning processor
    PetscReal     cell_size;         //!< representative cell size
    PetscInt      face;              //!< face indicator (used for BCs etc.)
} Acell;

//! \brief Overset donor cell (contains donor cell indices, owner, and distance)
typedef struct 
{
    PetscInt      indi, indj, indk;  //!< cell indices
    PetscMPIInt   rank;              //!< donor cell owner
    PetscReal     dist2p;            //!< squared distance from acceptor cell center
} Dcell;

//! \brief Overset motion structure (for dynamic overset motion)
typedef struct 
{
    Cmpnts        prescribedVel; //!< prescribed translational velocity
    word          motionType;    //!< type of motion ("Translation", "Rotation", etc.)
    PetscInt      setMotion;     //!< if nonzero, use the prescribed motion
    PetscInt      ibmAttached;   //!< if nonzero, use IBM–based motion instead
} oversetMotion;

//! \brief Main overset structure.  
//! This structure holds the parameters, connectivity, and interpolation data for
//! both the overset (child) mesh and the base (background) mesh update.
struct overset_ 
{
    // Mesh connectivity information
    labelList    parentMeshId; //!< List of parent (background) domains
    labelList    childMeshId;  //!< List of attached overset (child) domains

    word         interpolationType; //!< Interpolation method ("trilinear", "inverseDist", "LS1", etc.)

    // Dynamic overset parameters
    PetscInt                    dynamicOverset; //!< Flag for dynamic overset (nonzero if active)
    PetscInt                    procChange;     //!< Flag if background processor intersections have changed
    oversetMotion               *oMotion;       //!< Pointer to overset motion data

    // MLS (Least–Squares) interpolation parameters
    PetscReal                   cellAvg;        //!< Average cell size (optional)
    PetscReal                   cellFactor;     //!< Factor to scale cell size for search radius

    // Communication-related arrays for processor connectivity (for donor/acceptor mapping)
    std::vector<MPI_Comm>       oset_comm;  

    std::vector<std::vector<PetscInt>>      AcellProcMat; //!< Processor connectivity matrix for acceptor cells
    std::vector<PetscInt>                   NumAcellPerProc; //!< Number of acceptor cells per processor
    
    std::vector<Acell>                      aCell;         //!< List of acceptor cells
    std::vector<std::vector<Dcell>>         dCell;    //!< List of donor cell  
    std::vector<Dcell>                      closestDonor;  //!< For trilinear interpolation: closest donor cell per acceptor

    std::vector<std::vector<PetscReal>>     DWeights; //!< MLS weights for overset interpolation 

    ibm_    *oibm;
    access_ *access;
};

#endif


// Initialization and update
PetscErrorCode InitializeOverset(domain_ *dm);

PetscErrorCode UpdateOversetInterpolation(domain_ *domain);

PetscErrorCode readOversetProperties(overset_ *os);

PetscErrorCode readBlankingIBMObject(overset_ *os, domain_ *domain);

PetscErrorCode updateAcceptorCoordinates(overset_ *os);

PetscErrorCode createAcceptorCellOverset(overset_ *os);

PetscErrorCode createAcceptorCellBackground(overset_ *os); 

PetscErrorCode findClosestDonor(mesh_ *meshDonor, mesh_ *meshAcceptor);

// Interpolation routines for overset (child) mesh
PetscErrorCode interpolateACellTrilinear(mesh_ *donorMesh, mesh_ *acceptorMesh);

PetscErrorCode oversetContravariantBC(mesh_ *mesh, PetscInt i, PetscInt j, PetscInt k, Cmpnts ucart, PetscInt face);

PetscErrorCode setOversetBC(mesh_ *meshA);

PetscErrorCode setBackgroundBC(mesh_ *meshA);

PetscErrorCode acell1Dcell0Connectivity(mesh_ *meshP, mesh_ *mesh);  

// Least–Squares weights (for overset interpolation from base donor cells)
PetscErrorCode getLSWeights(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode getLSWeights_2nd(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode getLSWeights_3rd(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode interpolateACellInvD(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode interpolateACellLS(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode updateAcellCoordinates(domain_ *domain);

PetscErrorCode updateIntersectingProcessors(domain_ *domain);

PetscErrorCode updateDonorCells(domain_ *domain);

// Dynamic overset motion
PetscErrorCode oversetMeshTranslation(overset_ *os);  

void sum_struct_Acell(void *in, void *inout, int *len, MPI_Datatype *type);
void defineStruct_Acell(MPI_Datatype *tstype);

PetscErrorCode computeOversetIBMElementNormal(ibm_ *ibm);

PetscErrorCode oversetIbmSearch(ibm_ *ibm);


