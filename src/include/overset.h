//! \file  overset.h
//! \brief overset Objects and functions,

#ifndef _OVERSET_H_
#define _OVERSET_H_

//! \brief overset acceptor cell
typedef struct
{

    PetscInt      indi, indj, indk;
    PetscReal     coorx, coory, coorz;
    PetscMPIInt   rank;
    PetscReal     cell_size;
    PetscInt      face;
} Acell;

//! \brief overset donor cell
typedef struct
{

    PetscInt      indi, indj, indk;
    PetscMPIInt   rank;
    PetscReal     dist2p;

} Dcell;

//! \brief struct with the overset motion
typedef struct
{

  Cmpnts        prescribedVel;
  word          motionType;
  PetscInt      setMotion;                                                       //!< if set to true it uses the prescribed motion to move the overset mesh
  PetscInt      ibmAttached;                                                     //!< set to true if the setMotion is false so motion based on the movement of the IBM attached to mesh

} oversetMotion;


//! \brief Structure for the Overset mesh method
struct overset_
{
    labelList           parentMeshId;                                          //!< Parent node of each level
    labelList           childMeshId;                                           //!< Child node of each level

    word                interpolationType;                                     //!< type of Interpolation

    // dynamic overset parameters
    PetscInt            dynamicOverset;                                        //!< switch for dynamic overset
    PetscInt            procChange;                                            //!< switch to check if the background processors which intersect with overset mesh have changed
    oversetMotion       *oMotion;

    // MLS interpolation parameters
    PetscReal           cellAvg;
    PetscReal           cellFactor;                                            //!< cellFactor scales the cell size by a factor which will be used as the search radius

    // lists variables
    std::vector<MPI_Comm>                oset_comm;                            //!< communicator for overset - background processor interaction

    std::vector<Acell>                   aCell;                                //!< interpolated cells of the overset mesh

    std::vector <std::vector <Dcell>>    dCell;                                //!< MLS interpolation: donor cells of the background mesh
    std::vector <Dcell>                  closestDonor;                         //!< Trilinear interpolation: donor cells of the background mesh

    std::vector <std::vector <PetscInt>>      AcellProcMat;                         //!< rank matrix which indicates the processor connectivity between the Acell1 and Dcell0
    std::vector <PetscInt>                    NumAcellPerProc;                      //!< number of Acceptor cells in each processor
    std::vector <std::vector <PetscReal>>   DWeights;                             //!< MLS interpolation weights

    // access database
    access_ *access;
};

#endif

//! \brief initialize the overset variables
PetscErrorCode InitializeOverset(domain_ *dm);

//! \brief perform the overset interpolation
PetscErrorCode UpdateOversetInterpolation(domain_ *domain);

//! \brief read the overset properties file to set them
PetscErrorCode readOversetProperties(overset_ *overset);

//! \brief In dynamic overset, update the acceptor cell co-ordinates for the next time step
PetscErrorCode updateAcceptorCoordinates(overset_ *os);

//! \brief overset mesh translation based on the prescribed velocity
PetscErrorCode oversetMeshTranslation(overset_ *os);

//! \brief Create the list of acceptor cells which will be interpolated
PetscErrorCode createAcceptorCell(overset_ *os);

//! \brief Find the closest donor cell to the acceptor cell in the overset mesh for trilinear interpolation
PetscErrorCode findClosestDonor(mesh_ *meshP, mesh_ *mesh);

//! \brief Find the  donor cell to the acceptor cell in the overset mesh for Least square interpolation
PetscErrorCode acellDcellConnectivity(mesh_ *meshP, mesh_ *mesh);

//! \brief Find the weight for each donor cell when using least square interpolation first order
PetscErrorCode getLSWeights(mesh_ *meshP, mesh_ *mesh);

//! \brief Find the weight for each donor cell when using least square interpolation second order
PetscErrorCode getLSWeights_2nd(mesh_ *meshP, mesh_ *mesh);

//! \brief Find the weight for each donor cell when using least square interpolation third order
PetscErrorCode getLSWeights_3rd(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode OversetInterpolation(domain_ *domain);

PetscErrorCode interpolateACellInvD(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode interpolateACellTrilinear(mesh_ *meshP, mesh_ *mesh);

PetscErrorCode interpolateACellLS(mesh_ *meshP, mesh_ *mesh);

//! \brief set the contravariant flux at the face of an overset interpolated cell
PetscErrorCode oversetContravariantBC(mesh_ *mesh, PetscInt i, PetscInt j, PetscInt k, Cmpnts ucart, PetscInt face);

// functions not yet defined
PetscErrorCode updateAcellCoordinates(domain_ *domain);

PetscErrorCode updateIntersectingProcessors(domain_ *domain);

PetscErrorCode updateDonorCells(domain_ *domain);

PetscErrorCode oversetMeshTranslation(domain_ *domain);

//***************************************
//MPI Functions
void sum_struct_Acell(void *in, void *inout, int *len, MPI_Datatype *type);
void defineStruct_Acell(MPI_Datatype *tstype);
