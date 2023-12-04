//! \file  ibm.h
//! \brief IBM input header file.

// =============================================================================
// FUNCTIONS
// =============================================================================

//! \brief read the IBMProperties.dat file , IBM mesh files and allocate memory for the objects
PetscErrorCode readIBMProperties(ibm_ *ibm);

// \brief read the IBM mesh
PetscErrorCode readIBMObjectMesh(ibm_ *ibm, PetscInt b);

//! \brief read the ibm mesh in GRD format
PetscErrorCode readIBMBodyFileGRD(ibmObject *ibmBody);

//! \brief read the ibm mesh in ASCII format
PetscErrorCode readIBMBodyFileASCII(ibmObject *ibmBody);

//! \brief read the ibm mesh in ucd format
PetscErrorCode readIBMBodyFileUCD(ibmObject *ibmBody);

//! \brief read the ibm mesh in ucd2 format (without 0 tri entries)
PetscErrorCode readIBMBodyFileUCD2(ibmObject *ibmBody);

//! \brief read the ibm mesh in abaqus inp format
PetscErrorCode readIBMBodyFileAbaqusInp(ibmObject *ibmBody);

//! \brief read the ibm mesh in .grd format
PetscErrorCode readIBMSurfaceFileGRD(surface *ibmSurface);

//! \brief read the ibm mesh in ucd format
PetscErrorCode readIBMSurfaceFileUCD(surface *ibmSurface);

//! \brief read the ibm mesh in ucd2 format (without 0 tri entries)
PetscErrorCode readIBMSurfaceFileUCD2(surface *ibmSurface);

//! \brief read the ibm mesh in abaqus inp format
PetscErrorCode readIBMSurfaceFileAbaqusInp(surface *ibmSurface);

//!< \brief combine the mesh nodes and elements of different surface bodies into one ibm body
PetscErrorCode combineMesh(ibmObject *ibmBody);

//! \brief write the STL mesh
PetscErrorCode writeSTLFile(ibm_ *ibm, PetscInt b);

//! \brief find the reverse connectivity which lists the elements connected to a node.
PetscErrorCode nodeElementConnectivity(ibm_ *ibm);
