//! \file  wallmodel.h
//! \brief Wall models header file

//! \brief structure storing the Shumann wall models information
struct Shumann
{
    word             wfEvalType;
    PetscReal        roughness;
    PetscReal        thetaRef;
    PetscReal        gammaM;
    PetscReal        kappa;
    PetscReal        qWall;
};

//! \brief structure storing the Shumann wall models information
struct Cabot
{
    PetscReal        roughness;
};

//! \brief wall models container
struct wallModel
{
    Shumann          *wmShumann;
    Cabot            *wmCabot;
    patchVectorField tauWall;
    patchVectorField qWall;
};
