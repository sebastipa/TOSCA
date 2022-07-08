//! \file  ibm.h
//! \brief IBM input header file.

// =============================================================================
// FUNCTIONS
// =============================================================================

//! \brief read the IBMProperties.dat file , IBM mesh files and allocate memory for the objects
PetscErrorCode readIBMProperties(ibm_ *ibm);

// \brief read the IBM mesh
PetscErrorCode readIBMObjectMesh(ibm_ *ibm, PetscInt b);

//! \brief read the ibm mesh in ASCII format
PetscErrorCode readIBMBodyFileASCIIRaster(ibmObject *ibmBody);

//! \brief read the ibm mesh in ucd format
PetscErrorCode readIBMBodyFileUCD(ibmObject *ibmBody);

//! \brief read the ibm mesh in ucd2 format (without 0 tri entries)
PetscErrorCode readIBMBodyFileUCD2(ibmObject *ibmBody);

//! \brief read the ibm mesh in abaqus inp format
PetscErrorCode readIBMBodyFileAbaqusInp(ibmObject *ibmBody);

//! \brief read the ibm mesh in stl format
PetscErrorCode readIBMBodyFileSTL(ibmObject *ibmBody);

//! \brief read the ibm mesh in ASCII format
PetscErrorCode readIBMSurfaceFileASCIIRaster(surface *ibmSurface);

//! \brief read the ibm mesh in ucd format
PetscErrorCode readIBMSurfaceFileUCD(surface *ibmSurface);

//! \brief read the ibm mesh in ucd2 format (without 0 tri entries)
PetscErrorCode readIBMSurfaceFileUCD2(surface *ibmSurface);

//! \brief read the ibm mesh in abaqus inp format
PetscErrorCode readIBMSurfaceFileAbaqusInp(surface *ibmSurface);

//! \brief read the ibm mesh in stl format
PetscErrorCode readIBMSurfaceFileSTL(surface *ibmSurface);

//!< \brief combine the mesh nodes and elements of different surface bodies into one ibm body
PetscErrorCode combineMesh(ibmObject *ibmBody);

//! \brief write the STL mesh
PetscErrorCode writeSTLFile(ibmObject *ibmBody);
