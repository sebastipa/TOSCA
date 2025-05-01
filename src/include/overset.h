//! \file  overset.h
//! \brief Definitions for overset objects and functions for two–way coupling

#ifndef _OVERSET_H_
#define _OVERSET_H_

#include <unordered_map>

// Structure to hold hole object data
typedef struct {
    word bodyName;
    PetscInt ownerMesh;
    PetscInt donorMesh;
    word fileType;
    Cmpnts baseLocation;
    PetscInt searchCellRatio;
} HoleObject;

//! \brief Overset acceptor cell (used for both overset and base mesh acceptor cells)
typedef struct {
    PetscInt      indi, indj, indk;  //!< cell indices
    PetscReal     coorx, coory, coorz; //!< cell center coordinates
    PetscMPIInt   rank;              //!< owning processor
    PetscReal     cell_size;         //!< representative cell size
    PetscInt      face;              //!< face indicator (used for BCs etc.)
    PetscInt      donorId;
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

    // In general there are 2 sets of acceptor cells, acceptor cell at domain boundary(aCellDb) and acceptor cell at holeCut boundary(aCellHc)
    // All variables need to be created for both sets of acceptors 

    std::vector<std::vector<PetscInt>>      AcellProcMatDb;    //!< Processor connectivity matrix for acceptor cells
    std::vector<std::vector<PetscInt>>      AcellProcMatHc;    //!< Processor connectivity matrix for acceptor cells

    std::vector<PetscInt>                   NumAcellPerProcDb; //!< Number of acceptor cells per processor
    std::vector<PetscInt>                   NumAcellPerProcHc; //!< Number of acceptor cells per processor

    std::vector<Acell>                      aCellDb;           //!< List of acceptor cells
    std::vector<Acell>                      aCellHc;           //!< List of acceptor cells
 
    std::vector<Dcell>                      closestDonorDb;    //!< For trilinear interpolation: closest donor cell per acceptor at domain boundary
    std::vector<Dcell>                      closestDonorHc;    //!< For trilinear interpolation: closest donor cell per acceptor at hole cut boundary

    //for other interpolation methods
    std::vector<std::vector<PetscReal>>     DWeights;        //!< MLS weights for overset interpolation 
    std::vector<std::vector<Dcell>>         dCell;           //!< List of donor cell 

    // Communication-related arrays for processor connectivity (for donor/acceptor mapping)
    std::vector<MPI_Comm>       oset_comm;  

    // MLS (Least–Squares) interpolation parameters
    PetscReal                   cellAvg;        //!< Average cell size (optional)
    PetscReal                   cellFactor;     //!< Factor to scale cell size for search radius
    
    ibm_    *oibm;
    access_ *access;
};

//structs for 3D binning 
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

#endif


// Initialization and update
PetscErrorCode InitializeOverset(domain_ *dm);

PetscErrorCode UpdateOversetInterpolation(domain_ *domain);

PetscErrorCode UpdateDomainInterpolation(PetscInt d, domain_ *domain, PetscInt level);

PetscErrorCode findClosestDomainDonors(PetscInt d, domain_ *domain, PetscInt level, const std::vector<HoleObject> &holeObjects); 

PetscErrorCode findAcceptorCells(PetscInt d, domain_ *domain, PetscInt level, 
    const std::vector<HoleObject> &holeObjects);

PetscErrorCode readOversetProperties(overset_ *os);

PetscErrorCode readBlankingIBMObject(overset_ *os, domain_ *domain, char *holeObjectName, const std::vector<HoleObject> &holeObjects);

PetscErrorCode readHoleObjects(std::vector<HoleObject> &holeObjects, PetscInt numHoleObjects);

PetscErrorCode FindHoleObject(const std::vector<HoleObject> &holeObjects, 
    PetscInt parentId, PetscInt childId, 
    char **holeObjectName);

PetscErrorCode updateAcceptorCoordinates(overset_ *os);

PetscErrorCode createAcceptorCellOverset(overset_ *os);

PetscErrorCode createAcceptorCellBackground(overset_ *os, PetscInt donorMeshId);

PetscErrorCode findClosestDonorC2P(mesh_ *meshDonor, mesh_ *meshAcceptor, PetscInt donorId);

PetscErrorCode findClosestDonorP2C(mesh_ *meshDonor, mesh_ *meshAcceptor);

// Interpolation routines for overset (child) mesh
PetscErrorCode interpolateACellTrilinearParent2Current(mesh_ *donorMesh, mesh_ *acceptorMesh);

PetscErrorCode interpolateACellTrilinearChild2Current(mesh_ *meshD, mesh_ *meshA, PetscInt donorId);

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


