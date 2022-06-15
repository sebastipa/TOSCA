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
PetscErrorCode readIBMFileASCIIRaster(ibmObject *ibmBody);

//! \brief read the ibm mesh in ucd format
PetscErrorCode readIBMFileUCD(ibmObject *ibmBody);

//! \brief read the ibm mesh in abaqus inp format
PetscErrorCode readIBMFileAbaqusInp(ibmObject *ibmBody);

//! \brief write the STL mesh
PetscErrorCode writeSTLFile(ibmObject *ibmBody);
