//! \file  wallmodel.h
//! \brief Wall models header file

//! \brief structure storing the Shumann wall models information for U and T
struct Shumann
{
    word             wfEvalType;              //!< type of uStar evaluation localized/averaged
    PetscReal        kappa;                   //!< von karman constant (usually 0.4)
    PetscReal        thetaRef;                //!< reference potential temperature
    PetscReal        roughness;               //!< equivalent roughness height
    PetscReal        gammaM;                  //!< momentum gammaM from Paulson 1970 (stable BL)
    PetscReal        betaM;                   //!< momentum betaM from Pauslon 1970 (unstable BL)
    PetscReal        gammaH;                  //!< pot. temp. gammaH from Paulson 1970 (stable BL)
    PetscReal        betaH;                   //!< pot. temp. betaH from Paulson 1970 (stable BL)
    PetscReal        alphaH;                  //!< pot. temp. alphaH from Paulson 1970 (stable BL)
    PetscReal        tLast;                   //!< time at which the last theta update was done
    PetscReal        heatingRate;             //!< surface heating rate
    PetscReal        **surfaceTheta;          //!< surface temperature
    PetscInt         surfaceThetaSet;         //!< surface temperature has been initialized

};

//! \brief structure storing the Shumann wall models information
struct Cabot
{
    PetscReal        roughness;
    PetscReal        kappa;                   //!< von karman constant (usually 0.4)
};

//! \brief structure storing the Shumann wall models information
struct PowerLawAPG
{
    PetscReal        roughness;
    PetscReal        kappa;                   //!< von karman constant (usually 0.4)
};

struct LogLawAPG
{
    PetscReal        roughness;
    PetscReal        kappa;                   //!< von karman constant (usually 0.4)
};

//! \brief wall models container
struct wallModel
{
    Shumann          *wmShumann;
    Cabot            *wmCabot;
    PowerLawAPG      *wmPowerLawAPG;
    LogLawAPG        *wmLogLawAPG;
    
    patchVectorField tauWall;
    patchVectorField qWall;
};
