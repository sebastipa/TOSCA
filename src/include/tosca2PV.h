//! \file  tosca2PV.h
//! \brief post processing header file.

#ifndef _TOSCA2PV_H_
#define _TOSCA2PV_H_

//! \brief Structure defining the variables for postProcessing
typedef struct
{
  PetscInt               postProcessFields;
  PetscInt               writeRaster;
  PetscInt               samplingSections;

  PetscInt               postProcessPrecursor;
  precursor_             *precursor;
} postProcess;

#endif

// Top level functions
// ---------------------------------------------------------

//! \brief Reads binary fields and writes paraview data into XMF folder
PetscErrorCode binary3DToXMF(domain_ *domain, postProcess *pp);

//! \brief Reads binary i-section data and writes paraview data into XMF folder
PetscErrorCode binaryISectionsToXMF(domain_ *domain);

//! \brief Reads i-section data from average fields and writes paraview data into XMF folder
PetscErrorCode fieldISectionsToXMF(domain_ *domain);

//! \brief Reads binary i-section pertirbation data and writes paraview data into XMF folder
PetscErrorCode binaryISectionsPerturbToXMF(domain_ *domain);

//! \brief Reads binary j-section data and writes paraview data into XMF folder
PetscErrorCode binaryJSectionsToXMF(domain_ *domain, postProcess *pp);

//! \brief Reads j-section data from average fields and writes paraview data into XMF folder
PetscErrorCode fieldJSectionsToXMF(domain_ *domain);

//! \brief Reads binary j-section perturbation data and writes paraview data into XMF folder
PetscErrorCode binaryJSectionsPerturbToXMF(domain_ *domain, postProcess *pp);

//! \brief Reads binary k-section data and writes paraview data into XMF folder
PetscErrorCode binaryKSectionsToXMF(domain_ *domain);

//! \brief Reads k-section data from average fields and writes paraview data into XMF folder
PetscErrorCode fieldKSectionsToXMF(domain_ *domain);

//! \brief Reads binary k-section perturbation data and writes paraview data into XMF folder
PetscErrorCode binaryKSectionsPerturbToXMF(domain_ *domain);

//! \brief Initialize post processing parameters
PetscErrorCode postProcessInitialize(domain_ **domainAddr, clock_ *clock, simInfo_ *info, flags_ *flags);

//! \brief Initialize concurrent precursor post processing parameters
PetscErrorCode postProcessInitializePrecursor(postProcess *pp, clock_ *clock);

// Operation functions
// ---------------------------------------------------------
PetscErrorCode writeFieldsToXMF(domain_ *domain, const char* filexmf, PetscReal time);

PetscErrorCode getTimeList(const char* dataLoc, std::vector<PetscReal> &timeSeries, PetscInt &ntimes);

//! \brief Reads i,j,k - sections info and allocates moemry. Some info are not necessary thus not read.
PetscErrorCode sectionsReadAndAllocate(domain_ *domain);

//! Generate the section on-the-fly in the post processing phase (for average fields that only have to be done at the end)
PetscErrorCode kSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time);
PetscErrorCode kSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time);

//! \briefReads from k-slices time series and loads the velocity, temperature and nut planes. Important: assumes T and nut databases have the same times of U.
PetscErrorCode kSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time);
PetscErrorCode kSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt kplane, const word &fieldName, PetscReal time);

//! Generate the section on-the-fly in the post processing phase (for average fields that only have to be done at the end)
PetscErrorCode iSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time);
PetscErrorCode iSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time);

//! \briefReads from i-slices time series and loads the velocity, temperature and nut planes. Important: assumes T and nut databases have the same times of U.
PetscErrorCode iSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time);
PetscErrorCode iSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt iplane, const word &fieldName, PetscReal time);

//! Generate the section on-the-fly in the post processing phase (for average fields that only have to be done at the end)
PetscErrorCode jSectionLoadVectorFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time);
PetscErrorCode jSectionLoadScalarFromField(Vec &V, mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time);

//! \briefReads from j-slices time series and loads the velocity, temperature and nut planes. Important: assumes T and nut databases have the same times of U.
PetscErrorCode jSectionLoadVector(mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time);
PetscErrorCode jSectionLoadScalar(mesh_ *mesh, sections *sec, PetscInt jplane, const word &fieldName, PetscReal time);


// The following functions append text the the XMF file
// ----------------------------------------------------

//! \brief Opens a time section in the XMF file
void xmfWriteFileStartTimeSection
(
    FILE *Xmf, const char *FileXmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *Topology,
    PetscReal Time
);

//! \brief Closes a time section in the XMF file
void xmfWriteFileEndTimeSection
(
    FILE *Xmf, const char *FileXmf
);

//! \brief Writes geometry info in the XMF file
void xmfWriteFileGeometry
(
    FILE *Xmf, const char *FileXmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *PathSave
);

//! \brief Writes a symmetric tensor in the XMF files
void xmfWriteFileSymmTensor
(
    FILE *xmf, const char *filexmf,
    PetscInt size_x, PetscInt size_y, PetscInt size_z,
    const char *PathSave, const char *symmTensorName,
    const char * XX, const char * YY, const char *ZZ,
    const char * XY, const char * XZ, const char *YZ,
    const char *center = "Cell"
);

//! \brief Writes a vector in the XMF file
void xmfWriteFileVector
(
    FILE *xmf, const char *filexmf,
    PetscInt size_x, PetscInt size_y, PetscInt size_z,
    const char *PathSave, const char *Vecname,
    const char * V1, const char * V2, const char *V3,
    const char *center = "Cell"
);

//! \brief Writes a scalar in the XMF file
void xmfWriteFileScalar
(
    FILE *Xmf, const char *Filexmf,
    PetscInt Size_x, PetscInt Size_y, PetscInt Size_z,
    const char *PathSave, const char *ScalName,
    const char *Scal,
    const char *Center = "Cell"
);

// The following functions deal with the 1D-data HDF5 file creation
// ----------------------------------------------------------------

//! \brief Writes a dataset to HDF format file
void hdfWriteDataset
(
    hid_t *file_id,
    hid_t *dataspace_id,
    char const *var,
    float *x
);

// Interface functions for writing a single field (3D or 2D)
// ---------------------------------------------------------

//! \brief Writes a scalar field (appends to XMF and creates HDF)
PetscErrorCode writeScalarToXMF
(
    domain_    *domain,
    const char *filexmf,
    const char *hdfilen,
    hid_t	     *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Vec        V
);

//! \brief Writes a vector field (appends to XMF and creates HDF)
PetscErrorCode writeVectorToXMF
(
    domain_    *domain,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Vec        V
);

//! \brief Writes a symmetric tensor field (appends to XMF and creates HDF)
PetscErrorCode writeSymmTensorToXMF
(
    domain_    *domain,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Vec        V
);

//! \brief Writes a scalar defined on an i-section (appends to XMF and creates HDF)
PetscErrorCode writeISectionScalarToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    PetscReal     **field
);

//! \brief Writes a scalar defined on an j-section (appends to XMF and creates HDF)
PetscErrorCode writeJSectionScalarToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    PetscReal     **field
);

//! \brief Writes a scalar defined on an k-section (appends to XMF and creates HDF)
PetscErrorCode writeKSectionScalarToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    PetscReal     **field
);

//! \brief Writes a vector defined on an i-section (appends to XMF and creates HDF)
PetscErrorCode writeISectionVectorToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Cmpnts     **field
);

//! \brief Writes a vector defined on an j-section (appends to XMF and creates HDF)
PetscErrorCode writeJSectionVectorToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Cmpnts     **field
);

//! \brief Writes a vector defined on an k-section (appends to XMF and creates HDF)
PetscErrorCode writeKSectionVectorToXMF
(
    mesh_    *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	   *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    const char *fieldName,
    Cmpnts     **field
);

// Interface functions for writing the mesh (3D or 2D)
// ---------------------------------------------------------

//! \brief Writes 3D mesh (appends to XMF and creates HDF)
PetscErrorCode writePointsToXMF
(
    mesh_      *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	     *file_id,
    hid_t      *dataspace_id,
    PetscReal     time
);

//! \brief Writes 2D i-section mesh (appends to XMF and creates HDF)
PetscErrorCode writeISectionPointsToXMF
(
    mesh_      *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	     *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    PetscInt   iIndex
);

//! \brief Writes 2D j-section mesh (appends to XMF and creates HDF)
PetscErrorCode writeJSectionPointsToXMF
(
    mesh_      *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	     *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    PetscInt   jIndex
);

//! \brief Writes 2D k-section mesh (appends to XMF and creates HDF)
PetscErrorCode writeKSectionPointsToXMF
(
    mesh_      *mesh,
    const char *filexmf,
    const char *hdfilen,
    hid_t	     *file_id,
    hid_t      *dataspace_id,
    PetscReal     time,
    PetscInt   kIndex
);

void setToZero(float *vec, PetscInt n);

// Other file formats
// ---------------------------------------------------------

PetscErrorCode writeJSectionToRaster
(
    mesh_    *mesh,        // user context
    PetscInt        jIndex        // index of the j-section
);
