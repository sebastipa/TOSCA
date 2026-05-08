//! \file  overset.h
//! \brief Definitions for overset objects and functions for two–way coupling

#ifndef _OVERSET_H_
#define _OVERSET_H_

#include <unordered_map>
#include <map>

//! \brief Structure to hold hole object data
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
    PetscInt      indi, indj, indk;    //!< cell indices
    PetscReal     coorx, coory, coorz; //!< cell center coordinates
    PetscMPIInt   rank;                //!< owning processor
    PetscReal     cell_size;           //!< representative cell size
    PetscInt      face;                //!< face indicator (used for BCs etc.)
    PetscInt      donorId;
    PetscInt      parentCellId;
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

//! \brief Main overset structure. Holds the parameters, connectivity, and interpolation data between parent and child domains.
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

    // For other interpolation methods (currently deprecated)
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

//! \brief Struct for 3D binning of hole cut boundary elements for efficient donor search
struct BinIndex {
    int ix, iy, iz;
    bool operator==(const BinIndex& other) const {
        return ix == other.ix && iy == other.iy && iz == other.iz;
    }
};

//! \brief Hash function for BinIndex to be used in unordered_map
namespace std {
    template <>
    struct hash<BinIndex> {
        std::size_t operator()(const BinIndex& b) const {
            return ((std::hash<int>()(b.ix) ^ (std::hash<int>()(b.iy) << 1)) >> 1) ^ (std::hash<int>()(b.iz) << 1);
        }
    };
}

#endif

//! \brief Initializes overset environment and data structures
PetscErrorCode InitializeOverset(domain_ *dm);

//! \brief Set initial fields for all overset meshes iteratively
PetscErrorCode SetInitialFieldsOverset(PetscInt d, domain_ *domain);

//! \brief Update overset interpolation for a single time step
PetscErrorCode UpdateOversetInterpolation(domain_ *domain);

//! \brief Synchronize pressure gauge across all domains after overset interpolation
PetscErrorCode SyncPressureAcrossDomains(domain_ *domain);

//! \brief Update overset interpolation for given domains and their dependencies recursively
PetscErrorCode UpdateDomainInterpolation(PetscInt d, domain_ *domain, PetscInt level);

//! \brief Find donor cells from closest domain level
PetscErrorCode findClosestDomainDonors(PetscInt d, domain_ *domain, PetscInt level, const std::vector<HoleObject> &holeObjects); 

//! \brief Find acceptor cells from closest domain level for hole cut boundary
PetscErrorCode findAcceptorCells(PetscInt d, domain_ *domain, PetscInt level, const std::vector<HoleObject> &holeObjects);

//! \brief Read overset properties from input file
PetscErrorCode readOversetProperties(overset_ *os);

//! \brief Read overset blanking IBM object from input file
PetscErrorCode readBlankingIBMObject(overset_ *os, domain_ *domain, char *holeObjectName, const std::vector<HoleObject> &holeObjects);

//! \brief Read overset hole objects from input file
PetscErrorCode readHoleObjects(std::vector<HoleObject> &holeObjects, PetscInt numHoleObjects);

//! \brief Find the hole object for a parent-child pair
PetscErrorCode FindHoleObject(const std::vector<HoleObject> &holeObjects, PetscInt parentId, PetscInt childId, char **holeObjectName);

//! \brief Update coordinates of acceptor cells based on overset motion
PetscErrorCode updateAcceptorCoordinates(overset_ *os);

//! \brief Create acceptor cells at boundaries for a given domain
PetscErrorCode createAcceptorCellOverset(overset_ *os);

//! \brief Create acceptor cells at hole boundary for a given domain
PetscErrorCode createAcceptorCellBackground(overset_ *os, PetscInt donorMeshId);

//! \brief Find donor cells from child to parent 
PetscErrorCode findClosestDonorC2P(mesh_ *meshDonor, mesh_ *meshAcceptor, PetscInt donorId);

//! \brief Find donor cells from parent to child 
PetscErrorCode findClosestDonorP2C(mesh_ *meshDonor, mesh_ *meshAcceptor);

//! \brief Perform trilinear interpolation from parent to child
PetscErrorCode interpolateACellTrilinearP2C(mesh_ *donorMesh, mesh_ *acceptorMesh);

//! \brief Perform trilinear interpolation from child to parent
PetscErrorCode interpolateACellTrilinearC2P(mesh_ *meshD, mesh_ *meshA, PetscInt donorId);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode oversetContravariantBC(mesh_ *mesh, PetscInt i, PetscInt j, PetscInt k, Cmpnts ucart, PetscInt face);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode setOversetBC(mesh_ *meshA);

//! \brief Updates contracariant velocity at interpolated faces
PetscErrorCode setBackgroundBC(mesh_ *meshA);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode acell1Dcell0Connectivity(mesh_ *meshP, mesh_ *mesh);  

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode getLSWeights(mesh_ *meshP, mesh_ *mesh);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode getLSWeights_2nd(mesh_ *meshP, mesh_ *mesh);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode getLSWeights_3rd(mesh_ *meshP, mesh_ *mesh);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode interpolateACellInvD(mesh_ *meshP, mesh_ *mesh);

//! \brief DEPRECATED (COMMENTED OUT IN overset.c)
PetscErrorCode interpolateACellLS(mesh_ *meshP, mesh_ *mesh);

//! \brief DEPRECATED 
PetscErrorCode updateAcellCoordinates(domain_ *domain);

//! \brief DEPRECATED 
PetscErrorCode updateIntersectingProcessors(domain_ *domain);

//! \brief DEPRECATED 
PetscErrorCode updateDonorCells(domain_ *domain);

//! \brief Basic motion for dynamic overset mesh 
PetscErrorCode oversetMeshTranslation(overset_ *os);  

//! \brief MPI operation function to find the sum the elements of the vector of structs
void sum_struct_Acell(void *in, void *inout, int *len, MPI_Datatype *type);

//! \brief Define MPI datatype for Acell struct
void defineStruct_Acell(MPI_Datatype *tstype);

//! \brief Compute the normal vector for overset IBM hole cut boundary elements
PetscErrorCode computeOversetIBMElementNormal(ibm_ *ibm);

//! \brief Perform overset IBM search to find hole cut and the interpolated cells 
PetscErrorCode oversetIbmSearch(ibm_ *ibm);


