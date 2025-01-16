//! \file  ibm.h
//! \brief IBM model header file.

#ifndef IBM_H
#define IBM_H

#define MAX_ELEMENTS_PER_NODE 20

typedef struct Vertex Vertex;
typedef struct HalfEdge HalfEdge;
typedef struct Face Face;

struct Vertex{
    Cmpnts       nCoor;
    PetscInt     vertexId;
    HalfEdge     *edge;
};

struct HalfEdge{
    Vertex       *origin;
    HalfEdge     *next;
    HalfEdge     *twin;
    Face         *face;
};

struct Face{
    PetscInt     faceId;
    HalfEdge     *edge;
};

typedef struct {
    PetscInt elem[MAX_ELEMENTS_PER_NODE];
    PetscInt numConnected;
} ibmNode;

//! \brief Node struct
typedef struct node
{
    PetscInt              Node;             //! node id
    struct node           *next;            //! pointer to the next node
} node;

typedef struct cellNode
{
    cellIds         Node;                   //! cell id
    struct cellNode *next;                  //! pointer to the next cell
} cellNode;

//! \brief Node list
typedef struct list
{
    node         *head;                     //! pointer to the head of the integer node list
} list;

//! \brief cell Node list
typedef struct cellList
{
    cellNode     *head;                     //! pointer to the head of the cell list
} cellList;

typedef struct
{
    PetscInt     ncx, ncy, ncz;                //!< number of search cells in x,y and z direction
    PetscReal    dcx, dcy, dcz;                //!< search cell size in x,y and z direction
} searchBox;


typedef struct
{
    PetscInt     nodes;                        //!< number of nodes in the IBM Body
    PetscInt     elems;                        //!< number of elements in the IBM body

    Vertex       *vertices;
    HalfEdge     *halfEdges;
    Face         *faces;

    Cmpnts       *nCoor;                       //!< pointer to the co-ordinates of the nodes
    PetscInt	 *nID1 , *nID2 , *nID3 ;       //!< pointer to the 3 node ids of a triangular mesh element
    Cmpnts	     *eN;                          //!< pointers to the component of face normal in x, y and z direction of element
    Cmpnts       *eT1;                         //!< pointers to the component of face tangential1 (eT1 = eN x k, where k is unit normal along z, which is taken as a generic direction)
    Cmpnts       *eT2;                         //!< pointers to the component of face tangential2 (eT2 = eN x eT1)

    ibmNode      *ibmMeshNode;
    PetscInt     *eSurface;                    //!< pointer to the surfaceBody that this element belongs to
    PetscReal    *eA;                          //!< area of the element
    Cmpnts       *eCent;                       //!< coordinate of the element center
    Cmpnts       *nU;                          //!< velocity of the nodes
    Cmpnts       *nUPrev;                      //!< velocity of the node at the previous timestep

    Cmpnts       *eQVec;                       //!< center of the smallest bounding sphere of an ibm element
    PetscReal    *eRVec;                       //!< radius of the smallest bounding sphere of an ibm element


}ibmMesh;

typedef struct
{
    PetscReal     angSpeed;
    PetscReal     angAcc;
    PetscReal     rotAngle;
    PetscReal     maxR;

    Cmpnts        rotAxis;
    Cmpnts        rotCenter;

}ibmRotation;

typedef struct
{
    Cmpnts        transVelocity;
}ibmTranslation;

typedef struct
{
    PetscReal     amplitude;
    PetscReal     frequency;
    Cmpnts        motionDir;
    PetscReal     tPrev;

}ibmSineMotion;

typedef struct
{
    PetscReal     amplitude;
    PetscReal     frequency;
    PetscReal     initAngPosition;
    Cmpnts        pitchAxis;
    Cmpnts        pitchCenter;
    PetscReal     tPrev;

}ibmPitchMotion;

typedef struct
{
    word     surfaceName;

    word     surfaceFileType;

    PetscInt surfaceId;

    Cmpnts   baseLocation;

    word     elementSet;

    // element localToGlobalMapping
    PetscInt      *elementMapping;

    ibmMesh  *ibMsh;

    //ibm temp bc
    word          tempBC;
    word          tBCSetType;
    PetscInt      wallFunctionTypeT;
    PetscReal     fixedTemp;
    PetscReal     fixedFluxT;              //characteristic velocity used to find Nusselt number for fixedFLux BC only. Recomended to select average inlet velocoity or IBM tranalation speed.


}surface;

// struct to define the bounding boxes around each process that define the local ibm mesh elements within each processor
typedef struct
{
    PetscInt      *thisElemControlled;                        //!< flag telling if a ibm element is controlled by this processor
    PetscInt      *thisElemTransfered;                   //!< flag telling that the control of this point is being transferred to another processor

    boundingBox   innerZone;
    boundingBox   outerZone;

}elementBox;

typedef struct
{
    //ibm body parameters
    word          bodyName;
    word          bodyType;

    PetscInt      numSurfaces;
    PetscInt      bodyID;
    PetscInt      thinBody;

    word          elementSet;

    ibmMesh       *ibMsh;
    surface       **ibmSurface;

    word          fileType;

    //ibm wall model properties
    wallModel     *ibmWallModelU;
    wallModel     *ibmWallModelT;

    //ibm velocity bc
    word          velocityBC;
    word          uBCSetType;
    PetscInt      wallFunctionTypeU;

    //ibm temp bc
    word          tempBC;
    word          tBCSetType;
    PetscInt      wallFunctionTypeT;
    PetscReal     fixedTemp;
    PetscReal     fixedFluxT;
    PetscReal     uChar;   //global characteristic velocity anf length of simulation.
    PetscReal     lChar;


    //ibm motion variables
    Cmpnts         baseLocation;
    word           bodyMotion;
    ibmRotation    *ibmRot;
    ibmSineMotion  *ibmSine;
    ibmPitchMotion *ibmPitch;
    ibmTranslation *ibmTrans;
    PetscReal      startMove;                                 // start time for movement
    PetscReal      endMove;

    //ibm search parameters
    PetscReal     searchCellRatio;
    boundingBox   *bound;
    list          *searchCellList;

    //ibm force parameters
    Cmpnts        *ibmPForce;                                 //!< pressure force acting on each ibm element

    Cmpnts        procBoundCenter;                            //!< center of the processor bounding box for the ibm body
    Cmpnts        procBoundSize;                              //!< size of the processor bounding box in the x, y and z direction

    //element processor mapping for element normal projection - used in force calculation
    PetscInt      *thisPtControlled;                          //!< flag telling if a ibm mesh point is controlled by this processor based on normal projection
    PetscInt      *thisPtControlTransfer;                     //!< flag telling that the control of this point is being transferred to another processor
    cellIds       *closestCells;                              //!< closest cell id to the outward normal projection from the ibm mesh element

    //ibm source variables
    PetscReal     IBTemp;                                     //!< surface temperature of fixed temp IB sources
    PetscReal     IBTFlux;                                     //!< surface temperature of fixed temp IB sources

    //flags
    PetscInt       ibmControlled;                             //!< flag which tells if this proc controls this IBM body
    PetscInt       tSourceFlag;                               //!< set to 0 if no source, 1 if if fixed surface temp. and 2 if fixed heat flux

    // communication color
    MPI_Comm            IBM_COMM;                             //!< communicator for this IBM

    // element localToGlobalMapping
    PetscInt      *elementMapping;

    //local ibm elements box
    elementBox    *eBox;

    labelList     tSourceFlagSurf;                               //!< set to 0 if no source, 1 if if fixed surface temp. and 2 if fixed heat flux. Set as array for each surface in body
    labelList     uSourceFlagSurf;                               //!< set to 0 if no source, 1 if source
    labelList     smSourceFlagSurf;                               //!< set to 0 if no source, 1 if source

}ibmObject;

typedef struct
{
  cellIds         cellId;                                     //!< cell id of the ibm fluid cell
  cellList        flNodes;                                    //!< list of the cell ids of the fluid cells as support to the ibm fluid cell
  list            slNodes;                                    //!< list of the ibm mesh nodes as support to the ibm fluid cell
  list            bID;                                        //!< list of the body id corresponding to the ibm mesh node in the slNodes list

  PetscInt        numFl;                                      //!< number of fluid cells in the support
  PetscInt        numSl;                                      //!< number of ibm mesh nodes in the support

  PetscReal       rad;                                        //!< support radius of the ibm fluid cell

  PetscInt        closestElem;                                //!< closest ibm mesh element to the ibm fluid cell
  Cmpnts          pMin;                                       //!< projected point on the ibm mesh element from the ibm fluid cell
  PetscReal       minDist;                                    //!< dist to the ibm mesh element from the ibm fluid cell
  PetscInt        bodyID;
  Cmpnts          normal;

  PetscInt        sID;

  PetscReal       cs1, cs2, cs3;                              //!< ibm interpolation coefficient of the projected point from the ibm element nodes
  PetscReal       cr1, cr2, cr3;                              //!< ibm interpolation coefficient of the background mesh point from the background plane triangle nodes
  PetscInt        i1, j1, k1;                                 //!< first node (cell center id of fluid mesh) of the background interpolation plane
  PetscInt        i2, j2, k2;                                 //!< second node (cell center id of fluid mesh) of the background interpolation plane
  PetscInt        i3, j3, k3;                                 //!< third node (cell center id of fluid mesh) of the background interpolation plane
  PetscInt        imode;                                    //!< indicates which triangle (out 0f 8) the interpolated point belongs in the background plane

  PetscReal       dB;                                         //!< distance to the background plane from ibm fluid node along the IBM element normal
} ibmFluidCell;


//! \brief struct storing IBM model
struct ibm_
{
    Vec                lNvertFixed;
    PetscInt           numBodies;                     //!<  number of bodies
    word               IBInterpolationModel;          //!<  interpolation methodology
    word               curvibType;
    word               curvibOrder;
    ibmObject          **ibmBody;                     //!<  array of pointers to ibm objects

    searchBox          *sBox;                         //!< array of searchBox with number of search cells and their size for each ibm object

    ibmFluidCell       *ibmFCells;                    //!< cell ids and supports cells of the ibm fluid cells within each processor

    PetscInt           numIBMFluid;                   //!< number of ibm fluid cells within each processor

    // access database
    access_              *access;

    // flags
    PetscInt                 dbg;
    PetscInt             dynamic;
    PetscInt        computeForce;
    PetscInt         checkNormal;
    PetscInt       averageNormal;
    PetscInt         wallShearOn;
    PetscInt            writeSTL;
    PetscInt              ibmABL;

    PetscReal         interpDist;
    // output parameters
    PetscReal          timeStart;                    //!< start time of acquisition system
    word            intervalType;                    //!< timeStep: sample at every (timeInterval) iter, adjustableTime sample at every (timeInterval) seconds
    PetscReal       timeInterval;                    //!< acquisition time interval (overrides simulation time step if smaller and adjustableTime active)
};

#endif

// =============================================================================
// FUNCTIONS
// =============================================================================

//! \brief initialize ibm: top level function
PetscErrorCode InitializeIBM(ibm_ *ibm);

//! \brief update imb: top level function
PetscErrorCode UpdateIBM(ibm_ *ibm);

//! \brief writes the angular position when the body is rotating
PetscErrorCode writeIBMData(ibm_ *ibm, PetscInt b);

//! \brief writes the IBM force data for element and net force
PetscErrorCode  writeIBMForceData(ibm_ *ibm, PetscInt b, PetscReal *gElemPressure, PetscReal netMoment, PetscReal ibmPower, Cmpnts PresForce, Cmpnts ViscForce, Cmpnts momentVector);

//! \brief find in which processors the ibm body belongs
PetscErrorCode findIBMControlledProcs(ibm_ *ibm);

//! \brief create list of ibm element processors and their closest ibm fluid cell to the ibm element normal projection
PetscErrorCode initElementProjectionProcs(ibm_ *ibm);

//! \brief transfer ibm mesh element across processors based on their movement - for dynamic case
PetscErrorCode IBMProjectionProcessorTransfer(ibm_ *ibm);

//! \brief create inner and outer buffer zones for ibm mesh element parallelization
PetscErrorCode createProcessorBufferZones(ibm_ *ibm);

//! \brief create list of ibm element processors and their closest ibm fluid cell to the ibm element
PetscErrorCode initElementProcs(ibm_ *ibm);

//! \brief transfer ibm mesh element across processors based on their movement - for dynamic case
PetscErrorCode IBMElementProcessorTransfer(ibm_ *ibm);

//! \brief update the ibm mesh when the ibm body is moving
PetscErrorCode UpdateIBMesh(ibm_ *ibm);

//! \brief rotate the ibm mesh based on the angular speed input
PetscErrorCode rotateIBMesh(ibm_ *ibm, PetscInt b);

//! \brief prescribe sinusoidal motion for IBM body
PetscErrorCode sineMotion(ibm_ *ibm, PetscInt b);

//! \brief prescribe pitching oscillation motion for IBM body
PetscErrorCode pitchingMotion(ibm_ *ibm, PetscInt b);

//! \brief translate the ibm mesh based on the speed input
PetscErrorCode translateIBMesh(ibm_ *ibm, PetscInt b);

//! \brief set IBM wall model type and properties
PetscErrorCode setIBMWallModels(ibm_ *ibm);

//! \brief find the bounding box around an ibm body
PetscErrorCode findBodyBoundingBox(ibm_ *ibm);

//! \brief find the search cell (coarser mesh around the ibm) dimension
PetscErrorCode findSearchCellDim(ibm_ *ibm);

//! \brief ibm fluid node search performed using ray tracing algorithm
PetscErrorCode ibmSearch(ibm_ *ibm);

//! \brief find the closest ibm mesh element to a IBM fluid node
PetscErrorCode findClosestIBMElement(ibm_ *ibm);

PetscErrorCode findClosestIBMElement2Solid(ibm_ *ibm);

//! \brief find the interceptionPoint on the background mesh plane for curvib interpolation
PetscErrorCode findInterceptionPoint(ibm_ *ibm);

//! \brief interception point algorithm based on dividing the plane into 8 triangles
PetscErrorCode interceptionPt(Cmpnts pCoor, Cmpnts pc[9], Cmpnts eNorm, ibmFluidCell *ibF);

//! \brief find list of ibm fluid nodes within each processor
PetscErrorCode findIBMFluidCells(ibm_ *ibm);

//! \brief find the fluid nodes that act as support nodes to the given ibm fluid node
PetscErrorCode findFluidSupportNodes(ibm_ *ibm);

//! \brief find the ibm mesh nodes that act as support node to the given ibm fluid node
PetscErrorCode findIBMMeshSupportNodes(ibm_ *ibm);

//! \brief MLS interpolation algorithm
PetscErrorCode MLSInterpolation(ibm_ *ibm);

//! \brief CURVIB normal projection interpolation algorithm
PetscErrorCode CurvibInterpolation(ibm_ *ibm);

PetscErrorCode CurvibInterpolationTriangular(ibm_ *ibm);

PetscErrorCode CurvibInterpolationQuadratic(ibm_ *ibm);

PetscErrorCode CurvibInterpolationInternalCell(ibm_ *ibm);

PetscErrorCode ComputeForceMoment(ibm_ *ibm);

//! \brief compute shear stress at faces close to the IBM
PetscErrorCode findIBMWallShear(ibm_ *ibm);

//! \brief compute shear stress at faces close to the IBM
PetscErrorCode findIBMWallShearChester(ibm_ *ibm);

// check if a given body exists in the mesh based on the IBMProperties.dat file
PetscErrorCode checkIBMexists(ibm_ *ibm);

//! \brief insert ibm mesh elements as a list into the search cell that they belong to
PetscErrorCode createSearchCellList(ibm_ *ibm);

//! \brief destroy the ibm search cell list
PetscErrorCode destroyLists(ibm_ *ibm);

PetscErrorCode computeIBMElementNormal(ibm_ *ibm);

PetscErrorCode recomputeIBMeshProperties(ibm_ *ibm, PetscInt b);

PetscErrorCode createHalfEdgeDataStructure(ibm_ *ibm);

//! \brief ray casting algorithm
PetscReal rayCastingTest(Cmpnts p, ibmMesh *ibMsh, cellIds sCell, searchBox *sBox, boundingBox *ibBox, list *searchCellList);

//! \brief ray casting algorithm in local positive cell neighbourhood
PetscErrorCode rayCastLocal(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts p4, ibmMesh *ibMsh, cellIds sCell, searchBox *sBox, boundingBox *ibBox, list *searchCellList, PetscInt &intersect);
//! \brief find element bounding sphere
PetscErrorCode elementBoundingSphere(ibmObject *ibmBody);

//! \brief inline function to generate a random direction along the z index of the search cell for the ray casting algorithm
inline Cmpnts randomdirection(Cmpnts p, cellIds sCell, boundingBox *ibBox, searchBox *sBox, PetscInt seed);

//! \brief inline function to check if the ray intersects an ibm element
inline PetscInt intsectElement(Cmpnts p, Cmpnts dir, Cmpnts node1, Cmpnts node2, Cmpnts node3, PetscReal *t, PetscReal *u, PetscReal *v);

//! \brief inline function to check if the line intersects an ibm element
inline PetscBool isLineTriangleInt(Cmpnts p1, Cmpnts p2, ibmMesh *ibMsh, PetscInt e);

//! \brief inline function to check if a search cell is a support to an ibm fluid node
inline bool isSearchCellSupport(PetscReal rad, Cmpnts pt, cellIds cID, boundingBox *ibBox, searchBox *sBox);

//! \brief insert ibm nodes into the list of support nodes around an ibm fluid node
inline void insertIBMSupportNodes(Cmpnts pt, cellIds sCell, ibmFluidCell *ibF, ibmObject *ibBody, searchBox *sBox);

//! \brief Check whether 2D point p is located on the same side of line p1p2
inline PetscInt ISSameSide2D(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3);

//! \brief check whether point p is inside the triangle - 2D
inline PetscInt ISInsideTriangle2D(Cpt2D p, Cpt2D pa, Cpt2D pb, Cpt2D pc);

//! \brief check whether point p is inside the triangle - 3D
inline PetscInt isPointInTriangle(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts p3, Cmpnts norm);

//! \brief find the angle averaged normal about a vertex
inline Cmpnts computeVertexAverageNormal(ibmMesh *ibMsh, PetscInt vertexId);

//! \brief find the angle averaged normal about an edge
inline Cmpnts computeEdgeAverageNormal(ibmMesh *ibMsh, PetscInt vertexId, PetscInt faceId);

//! \brief find projection of point on a line and the distance to it
inline void disP2Line(Cmpnts p, Cmpnts p1, Cmpnts p2, Cmpnts *po, PetscReal *d, PetscReal *t);
//! \brief find the interpolation weights of a point inside a triangle from its nodes
inline void triangleIntp(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, ibmFluidCell *ibF);

//! \brief find the interpolation weights of a point inside a triangle from its nodes used for background grid
inline void triangleIntpBg(Cpt2D p, Cpt2D p1, Cpt2D p2, Cpt2D p3, ibmFluidCell *ibF);
