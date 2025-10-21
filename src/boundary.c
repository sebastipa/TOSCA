//! \file  boundary.c
//! \brief Contains boundary conditions

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/inflow.h"
#include "include/wallfunctions.h"


//***************************************************************************************************************//

PetscErrorCode readScalarBC(const word &location, const word &field, scalarBC *bc)
{
    word filename = location + field;

    PetscPrintf(PETSC_COMM_WORLD, "Reading %s boundary conditions in %s\n\n", field.c_str(), filename.c_str());

    // allocate memory
    PetscMalloc(sizeof(scalarBC), &(*bc));
    PetscMalloc(sizeof(word), &(bc->iLeft));
    PetscMalloc(sizeof(word), &(bc->iRight));
    PetscMalloc(sizeof(word), &(bc->jLeft));
    PetscMalloc(sizeof(word), &(bc->jRight));
    PetscMalloc(sizeof(word), &(bc->kLeft));
    PetscMalloc(sizeof(word), &(bc->kRight));

    // read the file and stores the string
    readDictWordAndDouble(filename.c_str(), "iLeft", &(bc->iLeft), &(bc->iLval));
    readDictWordAndDouble(filename.c_str(), "iRight", &(bc->iRight), &(bc->iRval));
    readDictWordAndDouble(filename.c_str(), "jLeft", &(bc->jLeft), &(bc->jLval));
    readDictWordAndDouble(filename.c_str(), "jRight", &(bc->jRight), &(bc->jRval));
    readDictWordAndDouble(filename.c_str(), "kLeft", &(bc->kLeft), &(bc->kLval));
    readDictWordAndDouble(filename.c_str(), "kRight", &(bc->kRight), &(bc->kRval));

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode readVectorBC(const word &location, const word &field, vectorBC *bc)
{
    word filename = location + field;

    PetscPrintf(PETSC_COMM_WORLD, "Reading %s boundary conditions in %s\n\n", field.c_str(), filename.c_str());

    // allocate memory
    PetscMalloc(sizeof(vectorBC), &(*bc));
    PetscMalloc(sizeof(word), &(bc->iLeft));
    PetscMalloc(sizeof(word), &(bc->iRight));
    PetscMalloc(sizeof(word), &(bc->jLeft));
    PetscMalloc(sizeof(word), &(bc->jRight));
    PetscMalloc(sizeof(word), &(bc->kLeft));
    PetscMalloc(sizeof(word), &(bc->kRight));

    // reads the file and stores the string
    readDictWordAndVector(filename.c_str(), "iLeft", &(bc->iLeft), &(bc->iLval));
    readDictWordAndVector(filename.c_str(), "iRight", &(bc->iRight), &(bc->iRval));
    readDictWordAndVector(filename.c_str(), "jLeft", &(bc->jLeft), &(bc->jLval));
    readDictWordAndVector(filename.c_str(), "jRight", &(bc->jRight), &(bc->jRval));
    readDictWordAndVector(filename.c_str(), "kLeft", &(bc->kLeft), &(bc->kLval));
    readDictWordAndVector(filename.c_str(), "kRight", &(bc->kRight), &(bc->kRval));

    return (0);
}

//***************************************************************************************************************//

PetscErrorCode SetBoundaryConditions(mesh_ *mesh)
{
    word location = "./boundary/" + mesh->meshName + "/";

    // read U boundary conditions
    readVectorBC(location, "U", &(mesh->boundaryU));

    // read nut boundary conditions
    if (mesh->access->flags->isLesActive)
    {
        readScalarBC(location, "nut", &(mesh->boundaryNut));
    }

    // read T boundary conditions
    if (mesh->access->flags->isTeqnActive)
    {
        readScalarBC(location, "T", &(mesh->boundaryT));
    }

    // check boundary conditions
    checkBCsAndSetPatchTypes(mesh);

    PetscBarrier(NULL);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkBCsAndSetPatchTypes(mesh_ *mesh)
{

    word location = "./boundary/" + mesh->meshName + "/";

    vectorBC boundaryU   = mesh->boundaryU;

    std::vector<PetscInt> flagU(6,0), flagT(6,0), flagNut(6,0), patchType(6,0);

    std::vector<word> UAvailableBC   = {"inletFunction","inletFunctionkLeft", "inletFunctionkRight", "inletFunctionjLeft","inletFunctionjRight",
                                        "inletFunctioniLeft", "inletFunctioniRight", "noSlip", "slip", "velocityWallFunction",
                                        "fixedValue", "zeroGradient", "periodic", "oversetInterpolate"};

    std::vector<word> TAvailableBC   = {"inletFunction","inletFunctionkLeft", "inletFunctionkRight", "inletFunctionjLeft","inletFunctionjRight",
                                        "inletFunctioniLeft", "inletFunctioniRight", "zeroGradient", "fixedValue", "thetaWallFunction",
                                        "fixedGradient", "periodic", "oversetInterpolate"};

    std::vector<word> nutAvailableBC = {"inletFunction","inletFunctionkLeft", "inletFunctionkRight", "inletFunctionjLeft","inletFunctionjRight",
                                        "inletFunctioniLeft", "inletFunctioniRight", "zeroGradient", "fixedValue",
                                        "periodic", "oversetInterpolate"};

    std::vector<word> wallPatchTypes = {"noSlip", "slip", "velocityWallFunction"};

    for(PetscInt i = 0; i < UAvailableBC.size(); i++)
    {
        if( boundaryU.kLeft  == UAvailableBC[i])
            flagU[0] = 1;

        if( boundaryU.kRight == UAvailableBC[i]) 
            flagU[1] = 1;

        if( boundaryU.jLeft  == UAvailableBC[i])
            flagU[2] = 1;

        if( boundaryU.jRight == UAvailableBC[i]) 
            flagU[3] = 1;

        if( boundaryU.iLeft  == UAvailableBC[i])
            flagU[4] = 1;

        if( boundaryU.iRight == UAvailableBC[i])
            flagU[5] = 1;
    }

    if(flagU[0] == 0)
    {
        char error[512];
        sprintf(error, "In %s/U, U boundary condition at kLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.kLeft.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(flagU[1] == 0)
    {
        char error[512];
        sprintf(error, "In %sU, U boundary condition at kRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.kRight.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(flagU[2] == 0)
    {
        char error[512];
        sprintf(error, "In %sU, U boundary condition at jLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.jLeft.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(flagU[3] == 0)
    {
        char error[512];
        sprintf(error, "In %sU, U boundary condition at jRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.jRight.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(flagU[4] == 0)
    {
        char error[512];
        sprintf(error, "In %sU, U boundary condition at iLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.iLeft.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(flagU[5] == 0)
    {
        char error[512];
        sprintf(error, "In %s/U, U boundary condition at iRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryU.iRight.c_str());
        fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
    }

    if(mesh->access->flags->isTeqnActive)
    {
        scalarBC boundaryT   = mesh->boundaryT;

        for(PetscInt i = 0; i < TAvailableBC.size(); i++)
        {
            if( boundaryT.kLeft  == TAvailableBC[i])
                flagT[0] = 1;

            if( boundaryT.kRight == TAvailableBC[i]) 
                flagT[1] = 1;

            if( boundaryT.jLeft  == TAvailableBC[i])
                flagT[2] = 1;

            if( boundaryT.jRight == TAvailableBC[i])
                flagT[3] = 1;

            if( boundaryT.iLeft  == TAvailableBC[i]) 
                flagT[4] = 1;

            if( boundaryT.iRight == TAvailableBC[i])
                flagT[5] = 1;
        }

        if(flagT[0] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at kLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryT.kLeft.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }

        if(flagT[1] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at kRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryT.kRight.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }

        if(flagT[2] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at jLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryT.jLeft.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }

        if(flagT[3] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at jRight = '%s' does not match with available BC\n", location.c_str(), boundaryT.jRight.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }

        if(flagT[4] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at iLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryT.iLeft.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }

        if(flagT[5] == 0)
        {
            char error[512];
            sprintf(error, "In %s/T, T boundary condition at iRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryT.iRight.c_str());
            fatalErrorInFunction("checkBCsAndSetPatchTypes", error);
        }
    }

    if(mesh->access->flags->isLesActive)
    {
        scalarBC boundaryNut = mesh->boundaryNut;

        for(PetscInt i = 0; i < nutAvailableBC.size(); i++)
        {
            if( boundaryNut.kLeft == nutAvailableBC[i]) 
                flagNut[0] = 1;

            if( boundaryNut.kRight == nutAvailableBC[i]) 
                flagNut[1] = 1;

            if( boundaryNut.jLeft == nutAvailableBC[i]) 
                flagNut[2] = 1;

            if( boundaryNut.jRight == nutAvailableBC[i]) 
                flagNut[3] = 1;

            if( boundaryNut.iLeft == nutAvailableBC[i]) 
                flagNut[4] = 1;

            if( boundaryNut.iRight == nutAvailableBC[i]) 
                flagNut[5] = 1;
        }

        if(flagNut[0] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at kLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryNut.kLeft.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }

        if(flagNut[1] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at kRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryNut.kRight.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }

        if(flagNut[2] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at jLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryNut.jLeft.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }

        if(flagNut[3] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at jRight = '%s' does not match with available BC\n", location.c_str(), boundaryNut.jRight.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }

        if(flagNut[4] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at iLeft = '%s' does not match with available BCs.\n", location.c_str(), boundaryNut.iLeft.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }

        if(flagNut[5] == 0)
        {
            char error[512];
            sprintf(error, "In %s/nut, nut boundary condition at iRight = '%s' does not match with available BCs.\n", location.c_str(), boundaryNut.iRight.c_str());
            fatalErrorInFunction("checkBoundaryConditions", error);
        }
    }

    // set patch types to zero
    mesh->boundaryU.iLeftPatchType  = 0;
    mesh->boundaryU.iRightPatchType = 0;
    mesh->boundaryU.jLeftPatchType  = 0;
    mesh->boundaryU.jRightPatchType = 0;      
    mesh->boundaryU.kLeftPatchType  = 0;
    mesh->boundaryU.kRightPatchType = 0;

    // set patch types based on velocity BCs
    for(PetscInt i = 0; i < wallPatchTypes.size(); i++)
    {
        if( boundaryU.kLeft  == wallPatchTypes[i])
            mesh->boundaryU.kLeftPatchType = 1;

        if( boundaryU.kRight == wallPatchTypes[i])
            mesh->boundaryU.kRightPatchType = 1;

        if( boundaryU.jLeft  == wallPatchTypes[i])
            mesh->boundaryU.jLeftPatchType = 1;

        if( boundaryU.jRight == wallPatchTypes[i])
            mesh->boundaryU.jRightPatchType = 1;

        if( boundaryU.iLeft  == wallPatchTypes[i])
            mesh->boundaryU.iLeftPatchType = 1;

        if( boundaryU.iRight == wallPatchTypes[i])
            mesh->boundaryU.iRightPatchType = 1;
    }

    return 0;
}

//***************************************************************************************************************//

PetscErrorCode SetPeriodicConnectivity(mesh_ *mesh, word &meshFileName)
{
    // initialize all types to zero
    mesh->i_periodic = mesh->ii_periodic = 0;
    mesh->j_periodic = mesh->jj_periodic = 0;
    mesh->k_periodic = mesh->kk_periodic = 0;

    PetscInt iType = 0, jType = 0, kType = 0;

    // for each boundary check and set
    if(mesh->boundaryU.iLeft == "periodic")
    {
        if(mesh->boundaryU.iRight!="periodic")
        {
           char error[512];
            sprintf(error, "i-left patch is periodic but opposite patch is not. Try setting periodic in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
        else
        {

            readDictInt(meshFileName.c_str(), "-iPeriodicType", &iType);

            if(iType == 1)
            {
                mesh->i_periodic = 1;
            }
            else if(iType == 2)
            {
                mesh->ii_periodic = 1;
            }
            else
            {
               char error[512];
                sprintf(error, "unknown periodic i connectivity type. Known types are 1, 2");
                fatalErrorInFunction("SetPeriodicConnectivity", error);
            }
        }
    }

    if (mesh->boundaryU.iRight=="periodic")
    {
        if(mesh->boundaryU.iLeft!="periodic")
        {
           char error[512];
            sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
    }

    if(mesh->boundaryU.jLeft == "periodic")
    {
        if(mesh->boundaryU.jRight!="periodic")
        {
           char error[512];
            sprintf(error, "j-left patch is periodic but opposite patch is not. Try setting periodic in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
        else
        {
            readDictInt(meshFileName.c_str(), "-jPeriodicType", &jType);

            if(jType == 1)
            {
                mesh->j_periodic = 1;
            }
            else if(jType == 2)
            {
                mesh->jj_periodic = 1;
            }
            else
            {
               char error[512];
                sprintf(error, "unknown periodic j connectivity type. Known types are 1, 2");
                fatalErrorInFunction("SetPeriodicConnectivity", error);
            }
        }
    }
    if (mesh->boundaryU.jRight=="periodic")
    {
        if(mesh->boundaryU.jLeft!="periodic")
        {
           char error[512];
            sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
    }

    if(mesh->boundaryU.kLeft == "periodic")
    {
        if(mesh->boundaryU.kRight!="periodic")
        {
           char error[512];
            sprintf(error, "k-left patch is periodic but opposite patch is not. Try setting periodic in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
        else
        {
            readDictInt(meshFileName.c_str(), "-kPeriodicType", &kType);

            if(kType == 1)
            {
                mesh->k_periodic = 1;
            }
            else if(kType == 2)
            {
                mesh->kk_periodic = 1;
            }
            else
            {
               char error[512];
                sprintf(error, "unknown periodic k connectivity type. Known types are 1, 2");
                fatalErrorInFunction("SetPeriodicConnectivity", error);
            }
        }
    }
    if (mesh->boundaryU.kRight=="periodic")
    {
        if(mesh->boundaryU.kLeft!="periodic")
        {
           char error[512];
            sprintf(error, "k-right patch is periodic but opposite patch is not. Try setting periodic in in boundary/U");
            fatalErrorInFunction("SetPeriodicConnectivity",  error);
        }
    }

    // check on sgs viscosity
    if(mesh->access->flags->isLesActive)
    {
        if(mesh->boundaryNut.iLeft == "periodic")
        {
            if(mesh->boundaryNut.iRight!="periodic")
            {
               char error[512];
                sprintf(error, "i-left patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryNut.iRight=="periodic")
        {
            if(mesh->boundaryNut.iLeft!="periodic")
            {
               char error[512];
                sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }

        if(mesh->boundaryNut.jLeft == "periodic")
        {
            if(mesh->boundaryNut.jRight!="periodic")
            {
               char error[512];
                sprintf(error, "j-left patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryNut.jRight=="periodic")
        {
            if(mesh->boundaryNut.jLeft!="periodic")
            {
               char error[512];
                sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }

        if(mesh->boundaryNut.kLeft == "periodic")
        {
            if(mesh->boundaryNut.kRight!="periodic")
            {
               char error[512];
                sprintf(error, "k-left patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryNut.kRight=="periodic")
        {
            if(mesh->boundaryNut.kLeft!="periodic")
            {
               char error[512];
                sprintf(error, "k-right patch is periodic but opposite patch is not. Try setting periodic in boundary/Nut");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
    }

    // check on potential temperature
    if(mesh->access->flags->isTeqnActive)
    {
        if(mesh->boundaryT.iLeft == "periodic")
        {
            if(mesh->boundaryT.iRight!="periodic")
            {
               char error[512];
                sprintf(error, "i-left patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryT.iRight=="periodic")
        {
            if(mesh->boundaryT.iLeft!="periodic")
            {
               char error[512];
                sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }

        if(mesh->boundaryT.jLeft == "periodic")
        {
            if(mesh->boundaryT.jRight!="periodic")
            {
               char error[512];
                sprintf(error, "j-left patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryT.jRight=="periodic")
        {
            if(mesh->boundaryT.jLeft!="periodic")
            {
               char error[512];
                sprintf(error, "i-right patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }

        if(mesh->boundaryT.kLeft == "periodic")
        {
            if(mesh->boundaryT.kRight!="periodic")
            {
               char error[512];
                sprintf(error, "k-left patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
        if (mesh->boundaryT.kRight=="periodic")
        {
            if(mesh->boundaryT.kLeft!="periodic")
            {
               char error[512];
                sprintf(error, "k-right patch is periodic but opposite patch is not. Try setting periodic in boundary/T");
                fatalErrorInFunction("SetPeriodicConnectivity",  error);
            }
        }
    }

    mesh->access->info->periodic
    =
    mesh->kk_periodic + mesh->k_periodic +
    mesh->jj_periodic + mesh->j_periodic +
    mesh->ii_periodic + mesh->i_periodic;

    return(0);
}
//***************************************************************************************************************//

PetscErrorCode UpdateContravariantBCs(ueqn_ *ueqn)
{
    mesh_        *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***ucont, ***lucont,
                  ***lucat,
                  ***icsi, ***jeta, ***kzet;
    PetscReal     ***nvert;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);

    DMDAVecGetArray(fda, ueqn->Ucont,  &ucont);
    DMDAVecGetArray(fda, ueqn->lUcont, &lucont);
    DMDAVecGetArray(fda, ueqn->lUcat,  &lucat);

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                // fixedValue type: interpolate the contravariant velocity
                if
                (
                    (mesh->boundaryU.kLeftPatchType == 0 && k == 0) ||
                    (mesh->boundaryU.kRightPatchType == 0 && k == mz-2)
                )
                {
                    ucont[k][j][i].z
                    =
                    0.5 * (lucat[k+1][j][i].x + lucat[k][j][i].x) * kzet[k][j][i].x +
                    0.5 * (lucat[k+1][j][i].y + lucat[k][j][i].y) * kzet[k][j][i].y +
                    0.5 * (lucat[k+1][j][i].z + lucat[k][j][i].z) * kzet[k][j][i].z;
                }

                if 
                (
                    (mesh->boundaryU.jLeftPatchType == 0 && j == 0) ||
                    (mesh->boundaryU.jRightPatchType == 0 && j == my-2)
                )
                {
                    ucont[k][j][i].y
                    =
                    0.5 * (lucat[k][j+1][i].x + lucat[k][j][i].x) * jeta[k][j][i].x +
                    0.5 * (lucat[k][j+1][i].y + lucat[k][j][i].y) * jeta[k][j][i].y +
                    0.5 * (lucat[k][j+1][i].z + lucat[k][j][i].z) * jeta[k][j][i].z;
                }

                if 
                (
                    (mesh->boundaryU.iLeftPatchType == 0 && i == 0) ||
                    (mesh->boundaryU.iRightPatchType == 0 && i == mx-2)
                )
                {
                    ucont[k][j][i].x
                    =
                    0.5 * (lucat[k][j][i+1].x + lucat[k][j][i].x) * icsi[k][j][i].x +
                    0.5 * (lucat[k][j][i+1].y + lucat[k][j][i].y) * icsi[k][j][i].y +
                    0.5 * (lucat[k][j][i+1].z + lucat[k][j][i].z) * icsi[k][j][i].z;
                }

                // no penetration type: set velocity at the boundary faces to zero
                if
                (
                    (mesh->boundaryU.iLeftPatchType == 1 && i==0) || 
                    (mesh->boundaryU.iRightPatchType == 1 && i==mx-2) 
                )
                {
                    ucont[k][j][i].x = 0.0;
                }
                if
                (
                    (mesh->boundaryU.jLeftPatchType == 1 && j==0) || 
                    (mesh->boundaryU.jRightPatchType == 1 && j==my-2) 
                )
                {
                    ucont[k][j][i].y = 0.0;
                }
                if
                (
                    (mesh->boundaryU.kLeftPatchType == 1 && k==0) || 
                    (mesh->boundaryU.kRightPatchType == 1 && k==mz-2) 
                )
                {
                    ucont[k][j][i].z = 0.0;
                }

                // zero gradient type: velocity is solved, but set flux to zero if reverse flow
                if (mesh->boundaryU.iLeft=="zeroGradient" && i==0)
                {
                    if(ucont[k][j][i].x > 0.0) ucont[k][j][i].x = 0.0;
                }
                if (mesh->boundaryU.iRight=="zeroGradient" && i==mx-2)
                {
                    if(ucont[k][j][i].x < 0.0) ucont[k][j][i].x = 0.0;
                }
                if (mesh->boundaryU.jLeft=="zeroGradient" && j==0)
                {
                    if(ucont[k][j][i].y > 0.0) ucont[k][j][i].y = 0.0;
                }
                if (mesh->boundaryU.jRight=="zeroGradient" && j==my-2)
                {
                    if(ucont[k][j][i].y < 0.0) ucont[k][j][i].y = 0.0;
                }
                if (mesh->boundaryU.kLeft=="zeroGradient" && k==0)
                {
                    if(ucont[k][j][i].z > 0.0) ucont[k][j][i].z = 0.0;
                }
                if (mesh->boundaryU.kRight=="zeroGradient" && k==mz-2)
                {
                    if(ucont[k][j][i].z < 0.0) ucont[k][j][i].z = 0.0;
                }

                // i,j,k  or ii, jj, kk periodic type: the right
                // boundary velocity has been solved in this case: put it on the
                // left boundary
                if
                (
                    mesh->boundaryU.iLeft=="periodic" &&
                    i==0
                )
                {
                    if(mesh->i_periodic) ucont[k][j][i].x = lucont[k][j][mx-2].x;
                    else if(mesh->ii_periodic) ucont[k][j][i].x = lucont[k][j][-2].x;
                }

                if
                (
                    mesh->boundaryU.iRight=="periodic" &&
                    i==mx-2
                )
                {
                    // do nothing, already have the fluxes
                }

                if
                (
                    mesh->boundaryU.jLeft=="periodic" &&
                    j==0
                )
                {
                    if(mesh->j_periodic) ucont[k][j][i].y = lucont[k][my-2][i].y;
                    else if(mesh->jj_periodic) ucont[k][j][i].y = lucont[k][-2][i].y;
                }

                if
                (
                    mesh->boundaryU.jRight=="periodic" &&
                    j==my-2
                )
                {
                    // do nothing, already have the fluxes
                }

                if
                (
                    mesh->boundaryU.kLeft=="periodic" &&
                    k==0
                )
                {
                    if(mesh->k_periodic)       ucont[k][j][i].z = lucont[mz-2][j][i].z;
                    else if(mesh->kk_periodic) ucont[k][j][i].z = lucont[-2][j][i].z;
                }

                if
                (
                    mesh->boundaryU.kRight=="periodic" &&
                    k==mz-2
                )
                {
                    // do nothing, already have the fluxes
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, ueqn->lUcont, &lucont);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &lucat);

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);

    // scatter new contravariant velocity values
    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateCartesianBCs(ueqn_ *ueqn)
{
    mesh_        *mesh = ueqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/U";

    Cmpnts        ***icsi, ***jeta, ***kzet;
    Cmpnts        ***csi,  ***eta,  ***zet;
    Cmpnts        ***cent;

    PetscReal     ***aj, ***iaj;
    PetscReal     ***nvert, ***ustar, ***meshTag;
    Cmpnts        ***ucat,  ***lucat,
                  ***lucont;

    PetscMPIInt   rank;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // initialize random number generator
    srand(time(NULL));

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);
    DMDAVecGetArray(fda, mesh->lCsi, &csi);
    DMDAVecGetArray(fda, mesh->lEta, &eta);
    DMDAVecGetArray(fda, mesh->lZet, &zet);
    DMDAVecGetArray(da,  mesh->lAj,  &aj);
    DMDAVecGetArray(da,  mesh->lIAj,  &iaj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecGetArray(fda, ueqn->Ucat,  &ucat);
    DMDAVecGetArray(fda, ueqn->lUcat,  &lucat);
    DMDAVecGetArray(fda, ueqn->lUcont,  &lucont);
    DMDAVecGetArray(da,  ueqn->lUstar, &ustar);

    // read inflow if necessary
    if(mesh->boundaryU.kLeft == "inletFunction")
    {
        inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

        // read the inflow data if necessary
        if (ifPtr->typeU == 3 || ifPtr->typeU == 4)
        {
            readInflowU(ifPtr, mesh->access->clock);

            if(ifPtr->typeU == 4)
            {
                // for each y coordinate, find the intex for the right shifted interpolation point
                if(ifPtr->shift2)
                {
                    PetscReal refRatio = ifPtr->shiftSpeed*ueqn->access->clock->time / mesh->bounds.Ly;
                    PetscReal refShift = (refRatio - floor(refRatio))*mesh->bounds.Ly;

                    // print information
                    /*
                    PetscPrintf(ifPtr->IFFCN_COMM,"Shift info:\n");
                    PetscPrintf(ifPtr->IFFCN_COMM," distance        = %.2f\n", refRatio);
                    PetscPrintf(ifPtr->IFFCN_COMM," shift (Ly=%.1f) = %.2f\n", mesh->bounds.Ly, refShift);
                    PetscPrintf(ifPtr->IFFCN_COMM," point specific\n");
                    */

                    for (i=lxs; i<lxe; i++)
                    {
                        PetscReal locRatio = (refShift + (cent[lzs][lys][i].y - mesh->bounds.ymin)) / mesh->bounds.Ly;
                        PetscReal locShift = (locRatio - floor(locRatio))*mesh->bounds.Ly + mesh->bounds.ymin;

                        PetscReal minDist = 1e20;
                        PetscInt  iClose  = 0;

                        for (PetscInt ii=1; ii<mx-1; ii++)
                        {
                            PetscReal dist = fabs(ifPtr->ycent[ii] - locShift);
                            if(dist < minDist)
                            {
                                minDist = dist;
                                iClose  = ii;
                            }
                        }

                        // make sure iClose is the right index
                        if(ifPtr->ycent[iClose] <= locShift)
                        {
                            iClose++;
                        }

                        // only set the ones belonging to this processor. This is the index that has to
                        // be sourced from the inflow data to apply the shift at index i
                        ifPtr->yIDs[i] = iClose;

                        // set the right weight
                        PetscReal delta    = ifPtr->ycent[ifPtr->yIDs[i]] - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        PetscReal dist     = locShift - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        ifPtr->yWeights[i] = dist / delta;

                        /*
                        PetscPrintf(ifPtr->IFFCN_COMM,"  distance = %.2f\n", locRatio);
                        PetscPrintf(ifPtr->IFFCN_COMM,"  shift    = %.2f\n", locShift);
                        PetscPrintf(ifPtr->IFFCN_COMM,"  iRight   = %ld\n", iClose);
                        PetscPrintf(ifPtr->IFFCN_COMM,"  weight   = %.2f\n", dist / delta);
                        */
                    }
                }
            }
        }
    }

    // set velocity boundary conditions
    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                PetscInt solid_flag = 0;

                if(isIBMSolidCell(k,j,i,nvert))
                {
                    mSetValue(ucat[k][j][i], 0);
                    continue;
                }

                if(isZeroedCell(k, j, i, meshTag))
                {
                    mSetValue(ucat[k][j][i], 0);
                    continue;
                }

                // wall functions: directly look at the wall function type,
                //                 note: zero means not allocated
                if
                (
                    (mesh->boundaryU.iLWF==-1 && i==1) ||
                    (mesh->boundaryU.iRWF==-1 && i==mx-2)
                )
                {
                    PetscReal roughness;

                    if(mesh->boundaryU.iLWF==-1)      roughness = ueqn->iLWM->wmCabot->roughness;
                    else if(mesh->boundaryU.iRWF==-1) roughness = ueqn->iRWM->wmCabot->roughness;

                    PetscReal area = nMag(csi[k][j][i]);

                    // get wall distance from 1st cell center and velocity
                    PetscReal sb = 0.5/aj[k][j][i]/area;
                    Cmpnts Ub    = nSet(lucat[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);
                    Cmpnts n = nSetFromComponents(ni[0], ni[1], ni[2]);

                    // get cartesian components of i-normal
                    PetscReal nx = ni[0], ny = ni[1], nz = ni[2];

                    // compute wall-normal velocity at point b
                    Cmpnts UbNormal = nScale(nDot(Ub, n), n);

                    // compute wall-parallel velocity at points b and c
                    Cmpnts UbParallel = nSub(Ub, UbNormal);

                    // compute magnitude
                    PetscReal UbParallelMag = nMag(UbParallel);

                    // get kinematic viscosity
                    PetscReal nu = mesh->access->constants->nu;

                    // compute uStar
                    if (roughness > 1.e-19) ustar[k][j][i] = utau_wf(nu, roughness, sb, UbParallelMag);
                    else ustar[k][j][i] = uTauCabot(nu, UbParallelMag, sb, 0.01, 0);

                    Cmpnts Ughost;

                    // set uTau and U ghost to 0 if U parallel is small
                    if(UbParallelMag<1.e-20)
                    {
                        ustar[k][j][i] = 0.0;
                        mSetValue(Ughost, 0.0);
                    }
                    // compute U ghost
                    else
                    {
                        if(isIBMCell(k,j,i,nvert))
                        {
                            mSetValue(Ughost, 0.0);
                        }

                        if(isOversetCell(k, j, i, meshTag))
                        {
                            mSetValue(Ughost, 0.0);
                        }

                        else
                        {
                            PetscReal tau_w = ustar[k][j][i]*ustar[k][j][i];

                            // Seba: added correction based on the fact that if we use central scheme the velocity on the face
                            //       can't be lower than zero. For coarse mesh is equivalent to no slip
                            PetscReal UghostParallelMag = PetscMax(UbParallelMag - 2.0*sb*tau_w/nu, -1.0*UbParallelMag);

                            // define UghostParallelDir
                            Cmpnts UghostParallelDir = nUnit(UbParallel);

                            // set Ughost parallel component
                            Ughost = nScale(UghostParallelMag, UghostParallelDir);

                            // add opposite normal component
                            mSub(Ughost, UbNormal);
                        }
                    }

                    // set U on the physical ghost cells
                    if(i==1)
                    {
                        mSet(ucat[k][j][i-1], Ughost);
                    }
                    else
                    {
                        mSet(ucat[k][j][i+1], Ughost);
                    }

                    solid_flag=1;
                }

                if
                (
                    (mesh->boundaryU.jLWF==-1 && j==1) ||
                    (mesh->boundaryU.jRWF==-1 && j==my-2)
                )
                {
                    PetscReal roughness;

                    if(mesh->boundaryU.jLWF==-1)      roughness = ueqn->jLWM->wmCabot->roughness;
                    else if(mesh->boundaryU.jRWF==-1) roughness = ueqn->jRWM->wmCabot->roughness;

                    // get the area at the boundary cell
                    PetscReal area = nMag(eta[k][j][i]);

                    // distance from wall and velocity at point b
                    PetscReal sb = 0.5/aj[k][j][i]/area;
                    Cmpnts    Ub = nSet(lucat[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // compute wall-normal velocity at points b and c
                    Cmpnts UbNormal = nScale(nDot(Ub, n), n);

                    // compute wall-parallel velocity at points b and c
                    Cmpnts UbParallel = nSub(Ub, UbNormal);

                    // compute magnitude
                    PetscReal UbParallelMag = nMag(UbParallel);

                    // get kinematic viscosity
                    PetscReal nu = mesh->access->constants->nu;

                    if (roughness > 1.e-19) ustar[k][j][i] = utau_wf (nu, roughness, sb, UbParallelMag);
                    else ustar[k][j][i] = uTauCabot(nu, UbParallelMag, sb, 0.01, 0);

                    Cmpnts Ughost;

                    if(UbParallelMag < 1.e-20)
                    {
                        ustar[k][j][i] = 0.0;
                        mSetValue(Ughost, 0.0);
                    }
                    else
                    {
                        if(isIBMCell(k,j,i,nvert))
                        {
                            mSetValue(Ughost, 0.0);
                        }

                        if(isOversetCell(k,j,i,meshTag))
                        {
                            mSetValue(Ughost, 0.0);
                        }
                        else
                        {
                            PetscReal tau_w = ustar[k][j][i]*ustar[k][j][i];

                            // seba: added correction based on the fact that if we use central scheme the velocity on the cell
                            //       can't be lower than zero. For coarse mesh is equivalent to no slip
                            PetscReal UghostParallelMag = PetscMax(UbParallelMag - 2.0*sb*tau_w/nu, -1.0*UbParallelMag);

                            // define UghostParallelDir
                            Cmpnts UghostParallelDir = nUnit(UbParallel);

                            // set Ughost parallel component
                            Ughost = nScale(UghostParallelMag, UghostParallelDir);

                            // add opposite normal component
                            mSub(Ughost, UbNormal);
                        }
                    }

                    if(j==1)
                    {
                        mSet(ucat[k][j-1][i], Ughost);
                    }
                    else
                    {
                        mSet(ucat[k][j+1][i], Ughost);
                    }

                    solid_flag=1;
                }

                // slip preserved normal gradient boundary condition
                // to be used in conjunction with Shumann ABL Wall Shear Stress model
                if
                (
                    (mesh->boundaryU.jLWF==-3 && j==1) ||
                    (mesh->boundaryU.jRWF==-3 && j==my-2)
                )
                {
                    // get the area at the boundary cell
                    PetscReal area = nMag(eta[k][j][i]);

                    // b is the boundary cell and c is first internal cell
                    PetscReal sb, sc;
                    Cmpnts    Uc, Ub;

                    // distance from wall and velocity at point b
                    sb = 0.5/aj[k][j][i]/area;
                    Ub = nSet(lucat[k][j][i]);

                    // distance from wall and velocity at point c
                    if (j==1)
                    {
                        Uc = nSet(lucat[k][j+1][i]);
                        sc = 1.0/aj[k][j][i]/area + 0.5/aj[k][j+1][i]/area;
                    }
                    else
                    {
                        Uc = nSet(lucat[k][j-1][i]);
                        sc = 1.0/aj[k][j][i]/area + 0.5/aj[k][j-1][i]/area;
                    }

                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // compute wall-normal velocity at points b and c
                    Cmpnts UbNormal = nScale(nDot(Ub, n), n);
                    Cmpnts UcNormal = nScale(nDot(Uc, n), n);

                    // compute wall-parallel velocity at points b and c
                    Cmpnts UbParallel = nSub(Ub, UbNormal);
                    Cmpnts UcParallel = nSub(Uc, UcNormal);

                    // compute magnitude
                    PetscReal UbParallelMag = nMag(UbParallel);

                    // compute parallel velocity at ghost cell
                    Cmpnts Ughost;

                    // uniform mesh at the wall
                    Ughost.x = 2. * UbParallel.x - UcParallel.x;
                    Ughost.y = 2. * UbParallel.y - UcParallel.y;
                    Ughost.z = 2. * UbParallel.z - UcParallel.z;

                    // compute normal gradient of the parallel component (to be preserved)
                    Cmpnts dudn1 = nSub(UcParallel, UbParallel); mScale(1.0/(sc-sb),  dudn1);
                    Cmpnts dudn0 = nSub(UbParallel, Ughost);     mScale(1.0/(2.0*sb), dudn0);

                    // add normal component to enforce no-penetration
                    mSub(Ughost, UbNormal);

                    if(j==1)
                    {
                        mSet(ucat[k][j-1][i], Ughost);
                    }
                    else
                    {
                        mSet(ucat[k][j+1][i], Ughost);
                    }

                    solid_flag=1;
                }

                if (k==1 && mesh->boundaryU.kLeft == "inletFunction")
                {
                    inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

                    // power law profile
                    if (ifPtr->typeU == 1)
                    {
                        PetscReal h = cent[k][j][i].z - mesh->grndLevel;
                        PetscReal x = cent[k][j][i].x;
                        PetscReal y = cent[k][j][i].y;
                        PetscReal z = cent[k][j][i].z;

                        // shear exponent
                        PetscReal alpha = 0.107027;

                        // power low contribution
                        ucat[k-1][j][i] = nScale(pow(h/ifPtr->Href, alpha), ifPtr->Uref);

                        // fluctuations contribution
                        PetscReal f          = ifPtr->uPrimeRMS;
                        PetscInt  n1         = rand() % 20000 - 10000; // RAND_MAX = 65535
                        PetscReal turbulence = ( 1 + ( (PetscReal)n1 ) / 10000. * f );

                        mScale(turbulence, ucat[k-1][j][i]);
                    }
                    // log-law profile
                    else if (ifPtr->typeU == 2)
                    {
                        // Compute scalar velocity

                        PetscReal h = cent[k][j][i].z - mesh->grndLevel;
                        PetscReal x = cent[k][j][i].x;
                        PetscReal y = cent[k][j][i].y;
                        PetscReal z = cent[k][j][i].z;

                        PetscReal vkConstant = 0.4;

                        PetscReal uMag;

                        if(h <= ifPtr->hInv)
                        {
                            uMag
                            =
                            PetscMax
                            (
                                (ifPtr->uTau/vkConstant)*std::log(h/ifPtr->roughness),
                                1e-5
                            );
                        }
                        else
                        {
                            uMag
                            =
                            (ifPtr->uTau/vkConstant)*std::log(ifPtr->hInv/ifPtr->roughness);
                        }

                        // set velocity according to log law
                        ucat[k-1][j][i] = nScale(uMag, ifPtr->Udir);
                    }
                    // Nieuwstadt model
                    else if (ifPtr->typeU == 5)
                    {
                        PetscReal h = cent[k][j][i].z - mesh->grndLevel;
                        ucat[k-1][j][i] = NieuwstadtInflowEvaluate(ifPtr, h);
                    }
                    // sinusoidal inflow
                    else if (ifPtr->typeU == 6)
                    {
                        PetscReal fraction = (cent[k][j][i].y - mesh->bounds.ymin)/mesh->bounds.Ly;
                        PetscReal uMag     = (1.0 + ifPtr->amplitude*std::cos(ifPtr->periods * 2.0 * M_PI * fraction)) * nMag(ifPtr->Uref);
                        ucat[k-1][j][i]    = nScale(uMag, ifPtr->Udir);
                    }
                    // unsteady mapped inflow
                    else if (ifPtr->typeU == 3)
                    {
                        // periodize inflow according to input

                        // compute period fraction (handle index = n case)
                        PetscInt jif = j % ifPtr->n1 == 0 ? ifPtr->n1 : j % ifPtr->n1;
                        PetscInt iif = i % ifPtr->n2 == 0 ? ifPtr->n2 : i % ifPtr->n2;

                        // index is less than nPrds times inflow points: have data
                        if
                        (
                            j<=ifPtr->n1*ifPtr->prds1 &&
                            i<=ifPtr->n2*ifPtr->prds2
                        )
                        {
                            if(ifPtr->merge1)
                            {
                                PetscReal height = cent[k][j][i].z - mesh->grndLevel;
                                PetscInt  IDs[2];
                                PetscReal Wg [2];

                                findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                                ucat[k-1][j][i].x = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    ifPtr->ucat_plane[jif][iif].x +
                                                    scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    (
                                                        ifPtr->uBarAvgTopX[IDs[0]].x * Wg[0] +
                                                        ifPtr->uBarAvgTopX[IDs[1]].x * Wg[1]
                                                    );
                                ucat[k-1][j][i].y = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    ifPtr->ucat_plane[jif][iif].y +
                                                    scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    (
                                                        ifPtr->uBarAvgTopX[IDs[0]].y * Wg[0] +
                                                        ifPtr->uBarAvgTopX[IDs[1]].y * Wg[1]
                                                    );
                                ucat[k-1][j][i].z = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    ifPtr->ucat_plane[jif][iif].z +
                                                    scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                                    (
                                                        ifPtr->uBarAvgTopX[IDs[0]].z * Wg[0] +
                                                        ifPtr->uBarAvgTopX[IDs[1]].z * Wg[1]
                                                    );
                            }
                            else
                            {
                                ucat[k-1][j][i].x = ifPtr->ucat_plane[jif][iif].x;
                                ucat[k-1][j][i].y = ifPtr->ucat_plane[jif][iif].y;
                                ucat[k-1][j][i].z = ifPtr->ucat_plane[jif][iif].z;
                            }
                        }
                        // index is more than nPrds times inflow points: extrapolate
                        else
                        {
                            if(ifPtr->merge1)
                            {
                                ucat[k-1][j][i].x = ifPtr->uBarAvgTopX[9].x;
                                ucat[k-1][j][i].y = ifPtr->uBarAvgTopX[9].y;
                                ucat[k-1][j][i].z = ifPtr->uBarAvgTopX[9].z;
                            }
                            else
                            {
                                // extrapolate along j
                                if(j>ifPtr->n1*ifPtr->prds1) jif = ifPtr->n1;

                                // extrapolate along i
                                if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                                ucat[k-1][j][i].x = ifPtr->ucat_plane[jif][iif].x;
                                ucat[k-1][j][i].y = ifPtr->ucat_plane[jif][iif].y;
                                ucat[k-1][j][i].z = ifPtr->ucat_plane[jif][iif].z;
                            }
                        }
                    }
                    else if (ifPtr->typeU == 4)
                    {
                        Cmpnts uGhost;

                        if(ifPtr->shift2)
                        {
                            PetscInt  iLeft    = ifPtr->yIDs[i]-1,
                                      iRight   = ifPtr->yIDs[i];

                            Cmpnts uLeft;
                            uLeft.x =
                                ifPtr->inflowWeights[j][iLeft][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][0].j][ifPtr->closestCells[j][iLeft][0].i].x +
                                ifPtr->inflowWeights[j][iLeft][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][1].j][ifPtr->closestCells[j][iLeft][1].i].x +
                                ifPtr->inflowWeights[j][iLeft][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][2].j][ifPtr->closestCells[j][iLeft][2].i].x +
                                ifPtr->inflowWeights[j][iLeft][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][3].j][ifPtr->closestCells[j][iLeft][3].i].x;

                            if(ifPtr->interpMethod == "spline")
                            {
                                uLeft.y =
                                    ifPtr->inflowWeights_2[j][iLeft][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][0].j][ifPtr->closestCells_2[j][iLeft][0].i].y +
                                    ifPtr->inflowWeights_2[j][iLeft][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][1].j][ifPtr->closestCells_2[j][iLeft][1].i].y +
                                    ifPtr->inflowWeights_2[j][iLeft][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][2].j][ifPtr->closestCells_2[j][iLeft][2].i].y +
                                    ifPtr->inflowWeights_2[j][iLeft][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][3].j][ifPtr->closestCells_2[j][iLeft][3].i].y +
                                    ifPtr->inflowWeights_2[j][iLeft][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][4].j][ifPtr->closestCells_2[j][iLeft][4].i].y +
                                    ifPtr->inflowWeights_2[j][iLeft][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iLeft][5].j][ifPtr->closestCells_2[j][iLeft][5].i].y;

                                uLeft.z =
                                    ifPtr->inflowWeights_1[j][iLeft][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][0].j][ifPtr->closestCells_1[j][iLeft][0].i].z +
                                    ifPtr->inflowWeights_1[j][iLeft][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][1].j][ifPtr->closestCells_1[j][iLeft][1].i].z +
                                    ifPtr->inflowWeights_1[j][iLeft][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][2].j][ifPtr->closestCells_1[j][iLeft][2].i].z +
                                    ifPtr->inflowWeights_1[j][iLeft][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][3].j][ifPtr->closestCells_1[j][iLeft][3].i].z +
                                    ifPtr->inflowWeights_1[j][iLeft][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][4].j][ifPtr->closestCells_1[j][iLeft][4].i].z +
                                    ifPtr->inflowWeights_1[j][iLeft][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iLeft][5].j][ifPtr->closestCells_1[j][iLeft][5].i].z;
                            }
                            else
                            {
                                uLeft.y =
                                    ifPtr->inflowWeights[j][iLeft][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][0].j][ifPtr->closestCells[j][iLeft][0].i].y +
                                    ifPtr->inflowWeights[j][iLeft][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][1].j][ifPtr->closestCells[j][iLeft][1].i].y +
                                    ifPtr->inflowWeights[j][iLeft][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][2].j][ifPtr->closestCells[j][iLeft][2].i].y +
                                    ifPtr->inflowWeights[j][iLeft][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][3].j][ifPtr->closestCells[j][iLeft][3].i].y;

                                uLeft.z =
                                    ifPtr->inflowWeights[j][iLeft][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][0].j][ifPtr->closestCells[j][iLeft][0].i].z +
                                    ifPtr->inflowWeights[j][iLeft][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][1].j][ifPtr->closestCells[j][iLeft][1].i].z +
                                    ifPtr->inflowWeights[j][iLeft][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][2].j][ifPtr->closestCells[j][iLeft][2].i].z +
                                    ifPtr->inflowWeights[j][iLeft][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iLeft][3].j][ifPtr->closestCells[j][iLeft][3].i].z;
                            }

                            Cmpnts uRight;
                            uRight.x =
                                ifPtr->inflowWeights[j][iRight][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][0].j][ifPtr->closestCells[j][iRight][0].i].x +
                                ifPtr->inflowWeights[j][iRight][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][1].j][ifPtr->closestCells[j][iRight][1].i].x +
                                ifPtr->inflowWeights[j][iRight][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][2].j][ifPtr->closestCells[j][iRight][2].i].x +
                                ifPtr->inflowWeights[j][iRight][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][3].j][ifPtr->closestCells[j][iRight][3].i].x;

                            if(ifPtr->interpMethod == "spline")
                            {
                                uRight.y =
                                    ifPtr->inflowWeights_2[j][iRight][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][0].j][ifPtr->closestCells_2[j][iRight][0].i].y +
                                    ifPtr->inflowWeights_2[j][iRight][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][1].j][ifPtr->closestCells_2[j][iRight][1].i].y +
                                    ifPtr->inflowWeights_2[j][iRight][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][2].j][ifPtr->closestCells_2[j][iRight][2].i].y +
                                    ifPtr->inflowWeights_2[j][iRight][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][3].j][ifPtr->closestCells_2[j][iRight][3].i].y +
                                    ifPtr->inflowWeights_2[j][iRight][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][4].j][ifPtr->closestCells_2[j][iRight][4].i].y +
                                    ifPtr->inflowWeights_2[j][iRight][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][iRight][5].j][ifPtr->closestCells_2[j][iRight][5].i].y;

                                uRight.z =
                                    ifPtr->inflowWeights_1[j][iRight][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][0].j][ifPtr->closestCells_1[j][iRight][0].i].z +
                                    ifPtr->inflowWeights_1[j][iRight][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][1].j][ifPtr->closestCells_1[j][iRight][1].i].z +
                                    ifPtr->inflowWeights_1[j][iRight][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][2].j][ifPtr->closestCells_1[j][iRight][2].i].z +
                                    ifPtr->inflowWeights_1[j][iRight][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][3].j][ifPtr->closestCells_1[j][iRight][3].i].z +
                                    ifPtr->inflowWeights_1[j][iRight][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][4].j][ifPtr->closestCells_1[j][iRight][4].i].z +
                                    ifPtr->inflowWeights_1[j][iRight][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][iRight][5].j][ifPtr->closestCells_1[j][iRight][5].i].z;
                            }
                            else
                            {
                                uRight.y =
                                    ifPtr->inflowWeights[j][iRight][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][0].j][ifPtr->closestCells[j][iRight][0].i].y +
                                    ifPtr->inflowWeights[j][iRight][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][1].j][ifPtr->closestCells[j][iRight][1].i].y +
                                    ifPtr->inflowWeights[j][iRight][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][2].j][ifPtr->closestCells[j][iRight][2].i].y +
                                    ifPtr->inflowWeights[j][iRight][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][3].j][ifPtr->closestCells[j][iRight][3].i].y;

                                uRight.z =
                                    ifPtr->inflowWeights[j][iRight][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][0].j][ifPtr->closestCells[j][iRight][0].i].z +
                                    ifPtr->inflowWeights[j][iRight][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][1].j][ifPtr->closestCells[j][iRight][1].i].z +
                                    ifPtr->inflowWeights[j][iRight][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][2].j][ifPtr->closestCells[j][iRight][2].i].z +
                                    ifPtr->inflowWeights[j][iRight][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][iRight][3].j][ifPtr->closestCells[j][iRight][3].i].z;
                            }

                            uGhost = nSum(nScale(1.0-ifPtr->yWeights[i], uLeft), nScale(ifPtr->yWeights[i], uRight));
                        }
                        else
                        {
                            uGhost.x =
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].x +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].x +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].x +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].x;

                            if(ifPtr->interpMethod == "spline")
                            {
                                uGhost.y =
                                    ifPtr->inflowWeights_2[j][i][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][0].j][ifPtr->closestCells_2[j][i][0].i].y +
                                    ifPtr->inflowWeights_2[j][i][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][1].j][ifPtr->closestCells_2[j][i][1].i].y +
                                    ifPtr->inflowWeights_2[j][i][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][2].j][ifPtr->closestCells_2[j][i][2].i].y +
                                    ifPtr->inflowWeights_2[j][i][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][3].j][ifPtr->closestCells_2[j][i][3].i].y +
                                    ifPtr->inflowWeights_2[j][i][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][4].j][ifPtr->closestCells_2[j][i][4].i].y +
                                    ifPtr->inflowWeights_2[j][i][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_2[j][i][5].j][ifPtr->closestCells_2[j][i][5].i].y;

                                uGhost.z =
                                    ifPtr->inflowWeights_1[j][i][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][0].j][ifPtr->closestCells_1[j][i][0].i].z +
                                    ifPtr->inflowWeights_1[j][i][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][1].j][ifPtr->closestCells_1[j][i][1].i].z +
                                    ifPtr->inflowWeights_1[j][i][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][2].j][ifPtr->closestCells_1[j][i][2].i].z +
                                    ifPtr->inflowWeights_1[j][i][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][3].j][ifPtr->closestCells_1[j][i][3].i].z +
                                    ifPtr->inflowWeights_1[j][i][4] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][4].j][ifPtr->closestCells_1[j][i][4].i].z +
                                    ifPtr->inflowWeights_1[j][i][5] *
                                    ifPtr->ucat_plane[ifPtr->closestCells_1[j][i][5].j][ifPtr->closestCells_1[j][i][5].i].z;
                            }
                            else
                            {
                                uGhost.y =
                                    ifPtr->inflowWeights[j][i][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].y +
                                    ifPtr->inflowWeights[j][i][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].y +
                                    ifPtr->inflowWeights[j][i][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].y +
                                    ifPtr->inflowWeights[j][i][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].y;

                                uGhost.z =
                                    ifPtr->inflowWeights[j][i][0] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i].z +
                                    ifPtr->inflowWeights[j][i][1] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i].z +
                                    ifPtr->inflowWeights[j][i][2] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i].z +
                                    ifPtr->inflowWeights[j][i][3] *
                                    ifPtr->ucat_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i].z;
                            }
                        }

                        if(ifPtr->merge1)
                        {
                            PetscReal height = cent[k][j][i].z - mesh->grndLevel;
                            PetscInt  IDs[2];
                            PetscReal Wg [2];

                            findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                            ucat[k-1][j][i].x
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                uGhost.x
                            ) +
                            scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->uBarAvgTopX[IDs[0]].x * Wg[0] +
                                ifPtr->uBarAvgTopX[IDs[1]].x * Wg[1]
                            );

                            ucat[k-1][j][i].y
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                uGhost.y
                            ) +
                            scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->uBarAvgTopX[IDs[0]].y * Wg[0] +
                                ifPtr->uBarAvgTopX[IDs[1]].y * Wg[1]
                            );

                            ucat[k-1][j][i].z
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                uGhost.z
                            ) +
                            scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->uBarAvgTopX[IDs[0]].z * Wg[0] +
                                ifPtr->uBarAvgTopX[IDs[1]].z * Wg[1]
                            );
                        }
                        else
                        {
                            ucat[k-1][j][i] = nSet(uGhost);
                        }
                    }
                }

                // fixedValue boundary condition on i-left patch
                if (mesh->boundaryU.iLeft=="fixedValue" && i==1)
                {
                    ucat[k][j][i-1].x = mesh->boundaryU.iLval.x;
                    ucat[k][j][i-1].y = mesh->boundaryU.iLval.y;
                    ucat[k][j][i-1].z = mesh->boundaryU.iLval.z;
                }
                // fixedValue boundary condition on i-right patch
                if (mesh->boundaryU.iRight=="fixedValue" && i==mx-2)
                {
                    ucat[k][j][i+1].x = mesh->boundaryU.iRval.x;
                    ucat[k][j][i+1].y = mesh->boundaryU.iRval.y;
                    ucat[k][j][i+1].z = mesh->boundaryU.iRval.z;
                }
                // fixedValue boundary condition on j-left patch
                if (mesh->boundaryU.jLeft=="fixedValue" && j==1)
                {
                    ucat[k][j-1][i].x = mesh->boundaryU.jLval.x;
                    ucat[k][j-1][i].y = mesh->boundaryU.jLval.y;
                    ucat[k][j-1][i].z = mesh->boundaryU.jLval.z;
                }
                // fixedValue boundary condition on j-right patch
                if (mesh->boundaryU.jRight=="fixedValue" && j==my-2)
                {
                    ucat[k][j+1][i].x = mesh->boundaryU.jRval.x;
                    ucat[k][j+1][i].y = mesh->boundaryU.jRval.y;
                    ucat[k][j+1][i].z = mesh->boundaryU.jRval.z;
                }
                // fixedValue boundary condition on k-left patch
                if (mesh->boundaryU.kLeft=="fixedValue" && k==1)
                {
                    ucat[k-1][j][i].x = mesh->boundaryU.kLval.x;
                    ucat[k-1][j][i].y = mesh->boundaryU.kLval.y;
                    ucat[k-1][j][i].z = mesh->boundaryU.kLval.z;
                }
                // fixedValue boundary condition on k-right patch
                if (mesh->boundaryU.kRight=="fixedValue" && k==mz-2)
                {
                    ucat[k+1][j][i].x = mesh->boundaryU.kRval.x;
                    ucat[k+1][j][i].y = mesh->boundaryU.kRval.y;
                    ucat[k+1][j][i].z = mesh->boundaryU.kRval.z;
                }

                // slip on X left boundary - set on the physical ghost cell
                // works for arbitrary structured meshes
                if (mesh->boundaryU.iLeft=="slip" && i==1)
                {
                    // get cell normals (only direction is important)
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get i normal
                    Cmpnts n = nSetFromComponents(ni[0], ni[1], ni[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k][j][i-1], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j][i-1],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j][i-1],0);
                }

                // slip on X right boundary - set on the physical ghost cell
                // seba: works for arbitrary structured meshes
                if (mesh->boundaryU.iRight=="slip" && i==mx-2)
                {
                    // get cell normals (only direction is important)
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get i normal
                    Cmpnts n = nSetFromComponents(ni[0], ni[1], ni[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k][j][i+1], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j][i+1],0.0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j][i+1],0.0);

                }

                // slip on Y left boundary - set on the physical ghost cell
                // seba: works for arbitrary structured meshes
                if (mesh->boundaryU.jLeft=="slip" && j==1)
                {
                    // get cell normals (only direction is important)
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k][j-1][i], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j-1][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j-1][i],0);

                }

                // slip on Y right boundary - set on the physical ghost cell
                // seba: works for arbitrary structured meshes
                if (mesh->boundaryU.jRight=="slip" && j==my-2)
                {
                    // get cell normals and flip them
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k][j+1][i], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j+1][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j+1][i],0);

                }

                // slip on Z left boundary - set on the physical ghost cell
                // seba: works for arbitrary structured meshes
                if (mesh->boundaryU.kLeft=="slip" && k==1)
                {
                    // get cell normals (only direction is important)
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get i normal
                    Cmpnts n = nSetFromComponents(nk[0], nk[1], nk[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k-1][j][i], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k-1][j][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k-1][j][i],0);

                }

                // slip on Y right boundary - set on the physical ghost cell
                // seba: works for arbitrary structured meshes
                if (mesh->boundaryU.kRight=="slip" && k==mz-2)
                {
                    // get cell normals (only direction is important)
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get i normal
                    Cmpnts n = nSetFromComponents(nk[0], nk[1], nk[2]);

                    // get U at internal cells
                    Cmpnts Ucell = nSet(lucat[k][j][i]);

                    // compute 2x wall-normal velocity
                    Cmpnts UcellNormal2 = nScale(2.0*nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellGhost = nSub(Ucell, UcellNormal2);

                    // set normal wall component opposite on ghosts
                    mSet(ucat[k+1][j][i], UcellGhost);

                    // check IBM and eventually set to zero
                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k+1][j][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k+1][j][i],0);

                }

                // no-slip on X left boundary - set on the physical ghost cell
                if(i==1 && mesh->boundaryU.iLeft=="noSlip")
                {
                    mSetScale(-1.0, ucat[k][j][i-1], lucat[k][j][i]);
                    solid_flag=1;
                }

                // no-slip on X right boundary - set on the physical ghost cell
                if(i==mx-2 && mesh->boundaryU.iRight=="noSlip")
                {
                    mSetScale(-1.0, ucat[k][j][i+1], lucat[k][j][i]);
                    solid_flag=1;
                }

                // no-slip on Y left boundary - set on the physical ghost cell
                if(j==1 && mesh->boundaryU.jLeft=="noSlip")
                {
                    mSetScale(-1.0, ucat[k][j-1][i], lucat[k][j][i]);
                    solid_flag=1;
                }

                // no-slip on Y right boundary - set on the physical ghost cell
                if(j==my-2 && mesh->boundaryU.jRight=="noSlip")
                {
                    mSetScale(-1.0, ucat[k][j+1][i], lucat[k][j][i]);
                    solid_flag=1;
                }

                // no-slip on Z left boundary - set on the physical ghost cell
                if(k==1 && mesh->boundaryU.kLeft=="noSlip")
                {
                    mSetScale(-1.0, ucat[k-1][j][i], lucat[k][j][i]);
                    solid_flag=1;
                }

                // no-slip on Z right boundary - set on the physical ghost cell
                if(k==mz-2 && mesh->boundaryU.kRight=="noSlip")
                {
                    mSetScale(-1.0, ucat[k+1][j][i], lucat[k][j][i]);
                    solid_flag=1;
                }

                // outflow (zero gradient) in X left boundary
                if(mesh->boundaryU.iLeft=="zeroGradient" && i==1)
                {
                    ucat[k][j][i-1].x = 2*lucat[k][j][i].x - lucat[k][j][i+1].x;
                    ucat[k][j][i-1].y = 2*lucat[k][j][i].y - lucat[k][j][i+1].y;
                    ucat[k][j][i-1].z = 2*lucat[k][j][i].z - lucat[k][j][i+1].z;

                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j][i-1],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j][i-1],0);

                }
                // outflow (zero gradient) in X right boundary
                if (mesh->boundaryU.iRight=="zeroGradient" && i==mx-2)
                {
                    ucat[k][j][i+1].x = 2*lucat[k][j][i].x - lucat[k][j][i-1].x;
                    ucat[k][j][i+1].y = 2*lucat[k][j][i].y - lucat[k][j][i-1].y;
                    ucat[k][j][i+1].z = 2*lucat[k][j][i].z - lucat[k][j][i-1].z;

                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j][i+1],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j][i+1],0);

                }
                // outflow (zero gradient) in Y left boundary
                if (mesh->boundaryU.jLeft=="zeroGradient" && j==1)
                {
                    ucat[k][j-1][i].x = 2*lucat[k][j][i].x - lucat[k][j+1][i].x;
                    ucat[k][j-1][i].y = 2*lucat[k][j][i].y - lucat[k][j+1][i].y;
                    ucat[k][j-1][i].z = 2*lucat[k][j][i].z - lucat[k][j+1][i].z;

                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j-1][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j-1][i],0);

                }
                // outflow (zero gradient) in Y right boundary
                if (mesh->boundaryU.jRight=="zeroGradient" && j==my-2)
                {
                    ucat[k][j+1][i].x = 2*lucat[k][j][i].x - lucat[k][j-1][i].x;
                    ucat[k][j+1][i].y = 2*lucat[k][j][i].y - lucat[k][j-1][i].y;
                    ucat[k][j+1][i].z = 2*lucat[k][j][i].z - lucat[k][j-1][i].z;


                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j+1][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k][j+1][i],0);

                }
                // outflow (zero gradient) in Z left boundary
                if (mesh->boundaryU.kLeft=="zeroGradient" && k==1)
                {
                    ucat[k-1][j][i].x = 2*lucat[k][j][i].x - lucat[k+1][j][i].x;
                    ucat[k-1][j][i].y = 2*lucat[k][j][i].y - lucat[k+1][j][i].y;
                    ucat[k-1][j][i].z = 2*lucat[k][j][i].z - lucat[k+1][j][i].z;


                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k-1][j][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k-1][j][i],0);

                }
                // outflow (zero gradient) in Z right boundary
                if (mesh->boundaryU.kRight=="zeroGradient" && k==mz-2)
                {
                    ucat[k+1][j][i].x = 2*lucat[k][j][i].x - lucat[k-1][j][i].x;
                    ucat[k+1][j][i].y = 2*lucat[k][j][i].y - lucat[k-1][j][i].y;
                    ucat[k+1][j][i].z = 2*lucat[k][j][i].z - lucat[k-1][j][i].z;

                    if(isIBMCell(k,j,i,nvert)) mSetValue(ucat[k+1][j][i],0);

                    if(isOversetCell(k,j,i,meshTag)) mSetValue(ucat[k+1][j][i],0);

                }

                // i-periodic boundary condition on i-left patch
                if (mesh->boundaryU.iLeft=="periodic" && i==1)
                {
                    if(mesh->i_periodic) ucat[k][j][i-1] = lucat[k][j][mx-2];
                    else if (mesh->ii_periodic) ucat[k][j][i-1] = lucat[k][j][-2];

                    if ( isIBMCell(k,j,i,nvert)) mSetValue(ucat[k][j][i-1],0);
                }
                // i-periodic boundary condition on i-right patch
                if (mesh->boundaryU.iRight=="periodic" && i==mx-2)
                {

                    if(mesh->i_periodic) ucat[k][j][i+1] = lucat[k][j][1];
                    else if (mesh->ii_periodic) ucat[k][j][i+1] = lucat[k][j][mx+1];

                    if ( isIBMCell(k,j,i,nvert) ) mSetValue(ucat[k][j][i+1],0);
                }
                // j-periodic boundary condition on j-left patch
                if (mesh->boundaryU.jLeft=="periodic" && j==1)
                {

                    if(mesh->j_periodic) ucat[k][j-1][i] = lucat[k][my-2][i];
                    else if(mesh->jj_periodic) ucat[k][j-1][i] = lucat[k][-2][i];

                    if ( isIBMCell(k,j,i,nvert) ) mSetValue(ucat[k][j-1][i],0);
                }
                // j-periodic boundary condition on j-right patch
                if (mesh->boundaryU.jRight=="periodic" && j==my-2)
                {
                    if(mesh->j_periodic) ucat[k][j+1][i] = lucat[k][1][i];
                    else if(mesh->jj_periodic) ucat[k][j+1][i] = lucat[k][my+1][i];

                    if ( isIBMCell(k,j,i,nvert) ) mSetValue(ucat[k][j+1][i],0);
                }
                // k-periodic boundary condition on k-left patch
                if (mesh->boundaryU.kLeft=="periodic" && k==1)
                {

                    if(mesh->k_periodic) ucat[k-1][j][i] = lucat[mz-2][j][i];
                    else if(mesh->kk_periodic) ucat[k-1][j][i] = lucat[-2][j][i];

                    if ( isIBMCell(k,j,i,nvert)) mSetValue(ucat[k-1][j][i],0);
                }
                // k-periodic boundary condition on k-right patch
                if (mesh->boundaryU.kRight=="periodic" && k==mz-2)
                {
                    if(mesh->k_periodic) ucat[k+1][j][i] = lucat[1][j][i];
                    else if(mesh->kk_periodic) ucat[k+1][j][i] = lucat[mz+1][j][i];

                    if ( isIBMCell(k,j,i,nvert) ) mSetValue(ucat[k+1][j][i],0);
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);
    DMDAVecRestoreArray(fda, mesh->lCsi, &csi);
    DMDAVecRestoreArray(fda, mesh->lEta, &eta);
    DMDAVecRestoreArray(fda, mesh->lZet, &zet);
    DMDAVecRestoreArray(da,  mesh->lAj,  &aj);
    DMDAVecRestoreArray(da,  mesh->lIAj,  &iaj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecRestoreArray(fda, ueqn->Ucat,  &ucat);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &lucat);
    DMDAVecRestoreArray(fda, ueqn->lUcont,  &lucont);

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd  (fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    DMDAVecRestoreArray(da, ueqn->lUstar, &ustar);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateTemperatureBCs(teqn_ *teqn)
{
    mesh_          *mesh = teqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/T";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***t, ***lt, ***nvert, ***meshTag;
    PetscReal     ***aj, ***iaj;
    Cmpnts        ***csi, ***eta, ***zet, ***icsi, ***cent;

    // variables to recover inletFunction 3 lapse rate
    PetscReal        ldataHeight = 0, gdataHeight = 0;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    // read inflow if necessary
    if(mesh->boundaryT.kLeft == "inletFunction")
    {
        inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

        if (ifPtr->typeT == 3 || ifPtr->typeT == 4)
        {
            readInflowT(ifPtr, mesh->access->clock);

            // compute hight at which inflow data ends

            // type 3: inflow and actual meshes have same cell dimensions
            if (ifPtr->typeT == 3)
            {
                PetscInt lcount = 0, gcount = 0;

                // make sure this processor can access data
                if(lys <= ifPtr->n1*ifPtr->prds1 && ifPtr->n1*ifPtr->prds1 <= lye)
                {
                    // compute end of data height (j is vertical direction)
                    i = std::floor(0.5*(lxe-lxs) + lxs);
                    k = std::floor(0.5*(lze-lzs) + lzs);
                    ldataHeight = cent[k][ifPtr->n1*ifPtr->prds1][i].z;
                    lcount      = 1;
                }

                MPI_Allreduce(&ldataHeight, &gdataHeight, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
                MPI_Allreduce(&lcount, &gcount, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                gdataHeight = gdataHeight / gcount;
            }

            // type 4: inflow and actual meshes are different, use inflow width info
            else if (ifPtr->typeT == 4)
            {
                gdataHeight = ifPtr->inflowHeigth;

                // for each y coordinate, find the intex for the right shifted interpolation point
                if(ifPtr->shift2)
                {
                    PetscReal refRatio = ifPtr->shiftSpeed*teqn->access->clock->time / mesh->bounds.Ly;
                    PetscReal refShift = (refRatio - floor(refRatio))*mesh->bounds.Ly;

                    for (i=lxs; i<lxe; i++)
                    {
                        PetscReal locRatio = (refShift + (cent[lzs][lys][i].y - mesh->bounds.ymin)) / mesh->bounds.Ly;
                        PetscReal locShift = (locRatio - floor(locRatio))*mesh->bounds.Ly + mesh->bounds.ymin;

                        PetscReal minDist = 1e20;
                        PetscInt  iClose  = 0;

                        for (PetscInt ii=1; ii<mx-1; ii++)
                        {
                            PetscReal dist = fabs(ifPtr->ycent[ii] - locShift);
                            if(dist < minDist)
                            {
                                minDist = dist;
                                iClose  = ii;
                            }
                        }

                        // make sure iClose is the right index
                        if(ifPtr->ycent[iClose] <= locShift)
                        {
                            iClose++;
                        }

                        // only set the ones belonging to this processor. This is the index that has to
                        // be sourced from the inflow data to apply the shift at index i
                        ifPtr->yIDs[i] = iClose;

                        // set the right weight
                        PetscReal delta    = ifPtr->ycent[ifPtr->yIDs[i]] - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        PetscReal dist     = locShift - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        ifPtr->yWeights[i] = dist / delta;
                    }
                }
            }
        }
    }

    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);
    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  mesh->lIAj,   &iaj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecGetArray(da, teqn->lTmprt, &lt);
    DMDAVecGetArray(da, teqn->Tmprt,  &t);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // set to zero if solid
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    t[k][j][i] = teqn->access->abl->tRef;
                    continue;
                }

                // special boundary condition where inflow is mapped from precursor
                if(mesh->boundaryT.kLeft=="inletFunction" && k==1)
                {
                    inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

                    // rampanelli and zardi model
                    if(ifPtr->typeT == 2)
                    {
                        PetscReal b      = ifPtr->smear * ifPtr->gTop * ifPtr->dInv;
                        PetscReal a      = ifPtr->gInv - b;
                        PetscReal c      = ifPtr->smear * ifPtr->gABL * ifPtr->dInv;
                        PetscReal h      = cent[k][j][i].z - mesh->grndLevel;
                        PetscReal etaLim = ifPtr->hInv / ifPtr->smear / ifPtr->dInv;

                        // non dimensional height eta
                        PetscReal eta = (h - ifPtr->hInv) / ifPtr->smear / ifPtr->dInv;

                        // below BL and capping
                        if(eta < etaLim)
                        {
                            // non dimensional functions
                            PetscReal f_eta = (std::tanh(eta) + 1.0) / 2.0;
                            PetscReal g_eta = (std::log(2.0 * std::cosh(eta)) + eta) / 2.0;
                            PetscReal h_eta = (eta - std::log(2.0 * std::cosh(eta))) / 2.0;

                            // potential temperature
                            t[k-1][j][i] = ifPtr->tRef + a * f_eta + b * g_eta + c * h_eta + ifPtr->gABL*ifPtr->hInv;
                        }
                        // asymptotic behavior
                        else
                        {
                            // non dimensional functions
                            PetscReal f_eta = (std::tanh(eta) + 1.0) / 2.0;
                            PetscReal g_eta = (std::log(2.0 * std::cosh(eta)) + eta) / 2.0;
                            PetscReal h_eta = (eta - std::log(2.0 * std::cosh(eta))) / 2.0;

                            // potential temperature
                            t[k-1][j][i] = ifPtr->tRef + a * f_eta + b * g_eta + c * h_eta + ifPtr->gABL*ifPtr->hInv;

                            // Dries implementation (to add limit from below)
                            // gLim = (abs(eta) + eta)/2;

                            // potential temperature
                            // t[k-1][j][i] = ifPtr->tRef + a + b * eta;
                        }
                    }
                    // periodized mapped inflow
                    else if (ifPtr->typeT == 3)
                    {
                        // periodize inflow according to input

                        // compute period fraction (handle index = n case)
                        PetscInt jif = j % ifPtr->n1 == 0 ? ifPtr->n1 : j % ifPtr->n1;
                        PetscInt iif = i % ifPtr->n2 == 0 ? ifPtr->n2 : i % ifPtr->n2;

                        // index is less than nPrds times inflow points: have data
                        if
                        (
                            j<=ifPtr->n1*ifPtr->prds1 &&
                            i<=ifPtr->n2*ifPtr->prds2
                        )
                        {
                            if(ifPtr->merge1)
                            {
                                PetscReal height = cent[k][j][i].z - mesh->grndLevel;
                                PetscInt  IDs[2];
                                PetscReal Wg [2];

                                findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                                t[k-1][j][i] = scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                               ifPtr->t_plane[jif][iif] +
                                               scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                                               (
                                                   ifPtr->tBarAvgTopX[IDs[0]] * Wg[0] +
                                                   ifPtr->tBarAvgTopX[IDs[1]] * Wg[1]
                                               );
                            }
                            else
                            {
                                t[k-1][j][i] = ifPtr->t_plane[jif][iif];
                            }
                        }
                        // index is more than nPrds times inflow points: apply lapse rate
                        else
                        {
                            PetscReal delta = 0;

                            // extrapolate along j
                            if(j>ifPtr->n1*ifPtr->prds1)
                            {
                                jif   = ifPtr->n1;

                                delta = cent[k][j][i].z - ifPtr->avgTopLength;
                            }

                            // extrapolate along i
                            if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                            if(ifPtr->merge1)
                            {
                                t[k-1][j][i] = ifPtr->tBarAvgTopX[9] + delta * teqn->access->abl->gTop;
                            }
                            else
                            {
                                t[k-1][j][i] = ifPtr->t_plane[jif][iif] + delta * teqn->access->abl->gTop;;
                            }
                        }
                    }

                    // interpolated periodized mapped inflow
                    else if (ifPtr->typeT == 4)
                    {
                        PetscReal delta  = PetscMax(0.0, cent[k][j][i].z - ifPtr->avgTopLength);
                        PetscReal tGhost;

                        if(ifPtr->shift2)
                        {
                            PetscInt  iLeft    = ifPtr->yIDs[i]-1,
                                      iRight   = ifPtr->yIDs[i];

                            PetscReal tLeft =
                                ifPtr->inflowWeights[j][iLeft][0] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iLeft][0].j][ifPtr->closestCells[j][iLeft][0].i] +
                                ifPtr->inflowWeights[j][iLeft][1] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iLeft][1].j][ifPtr->closestCells[j][iLeft][1].i] +
                                ifPtr->inflowWeights[j][iLeft][2] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iLeft][2].j][ifPtr->closestCells[j][iLeft][2].i] +
                                ifPtr->inflowWeights[j][iLeft][3] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iLeft][3].j][ifPtr->closestCells[j][iLeft][3].i];

                            PetscReal tRight =
                                ifPtr->inflowWeights[j][iRight][0] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iRight][0].j][ifPtr->closestCells[j][iRight][0].i] +
                                ifPtr->inflowWeights[j][iRight][1] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iRight][1].j][ifPtr->closestCells[j][iRight][1].i] +
                                ifPtr->inflowWeights[j][iRight][2] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iRight][2].j][ifPtr->closestCells[j][iRight][2].i] +
                                ifPtr->inflowWeights[j][iRight][3] *
                                ifPtr->t_plane[ifPtr->closestCells[j][iRight][3].j][ifPtr->closestCells[j][iRight][3].i];

                            tGhost   = (1.0 - ifPtr->yWeights[i]) * tLeft + ifPtr->yWeights[i] * tRight;
                        }
                        else
                        {
                            tGhost =
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i] +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i] +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i] +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->t_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i];
                        }

                        if(ifPtr->merge1)
                        {
                            PetscReal height = cent[k][j][i].z - mesh->grndLevel;

                            PetscInt  IDs[2];
                            PetscReal Wg [2];
                            findInterpolationWeightsWithExtrap(Wg, IDs, ifPtr->avgTopPointCoords, 10, height);

                            t[k-1][j][i]
                            =
                            scaleHyperTangBot(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                tGhost
                            ) +
                            scaleHyperTangTop(height, ifPtr->avgTopLength, ifPtr->avgTopDelta) *
                            (
                                ifPtr->tBarAvgTopX[IDs[0]] * Wg[0] +
                                ifPtr->tBarAvgTopX[IDs[1]] * Wg[1]
                            ) +
                            delta * teqn->access->abl->gTop;
                        }
                        else
                        {
                            t[k-1][j][i] = tGhost + delta * teqn->access->abl->gTop;
                        }
                    }
                }

                // zeroGradient boundary condition on i-left patch
                if ((mesh->boundaryT.iLeft=="zeroGradient" || mesh->boundaryT.iLeft=="thetaWallFunction") && i==1)
                {
                    t[k][j][i-1] = lt[k][j][i];
                }
                // zeroGradient boundary condition on i-right patch
                if ((mesh->boundaryT.iRight=="zeroGradient" || mesh->boundaryT.iRight=="thetaWallFunction") && i==mx-2)
                {
                    t[k][j][i+1] = lt[k][j][i];
                }
                // zeroGradient boundary condition on j-left patch
                if ((mesh->boundaryT.jLeft=="zeroGradient" || mesh->boundaryT.jLeft=="thetaWallFunction") && j==1)
                {
                    if(mesh->boundaryT.jLeft=="thetaWallFunction" && mesh->boundaryT.jLWF==-4)
                    {
                        if(teqn->access->clock->it==teqn->access->clock->itStart)
                        {
                            t[k][j-1][i] = lt[k][j][i];
                        }
                        else 
                        {
                            Shumann *wm  = teqn->jLWM->wmShumann;

                            // find interpolation weights for surface temp
                            PetscReal w[2];
                            PetscInt  l[2];

                            findInterpolationWeights(w, l, wm->timeVec, wm->numT, teqn->access->clock->time);

                            PetscReal surfaceTemp = w[0] * wm->surfTemp[l[0]] + w[1] * wm->surfTemp[l[1]];

                            t[k][j-1][i] = 2.0*surfaceTemp - lt[k][j][i];
                        }
                    }
                    else 
                    {
                        t[k][j-1][i] = lt[k][j][i];
                    }
                }
                // zeroGradient boundary condition on j-right patch
                if ((mesh->boundaryT.jRight=="zeroGradient" || mesh->boundaryT.jRight=="thetaWallFunction") && j==my-2)
                {
                    t[k][j+1][i] = lt[k][j][i];
                }
                // zeroGradient boundary condition on k-left patch
                if (mesh->boundaryT.kLeft=="zeroGradient" && k==1)
                {
                    t[k-1][j][i] = lt[k][j][i];
                }
                // zeroGradient boundary condition on k-right patch
                if (mesh->boundaryT.kRight=="zeroGradient" && k==mz-2)
                {
                    t[k+1][j][i] = lt[k][j][i];
                }

                // fixedValue boundary condition on i-left patch
                if (mesh->boundaryT.iLeft=="fixedValue" && i==1)
                {
                    t[k][j][i-1] = mesh->boundaryT.iLval;
                }
                // fixedValue boundary condition on i-right patch
                if (mesh->boundaryT.iRight=="fixedValue" && i==mx-2)
                {
                    t[k][j][i+1] = mesh->boundaryT.iRval;
                }
                // fixedValue boundary condition on j-left patch
                if (mesh->boundaryT.jLeft=="fixedValue" && j==1)
                {
                    t[k][j-1][i] = mesh->boundaryT.jLval;
                }
                // fixedValue boundary condition on j-right patch
                if (mesh->boundaryT.jRight=="fixedValue" && j==my-2)
                {
                    t[k][j+1][i] = mesh->boundaryT.jRval;
                }
                // fixedValue boundary condition on k-left patch
                if (mesh->boundaryT.kLeft=="fixedValue" && k==1)
                {
                    t[k-1][j][i] = mesh->boundaryT.kLval;
                }
                // fixedValue boundary condition on k-right patch
                if (mesh->boundaryT.kRight=="fixedValue" && k==mz-2)
                {
                    t[k+1][j][i] = mesh->boundaryT.kRval;
                }

                // fixedGradient boundary condition on i-left patch
                if (mesh->boundaryT.iLeft=="fixedGradient" && i==1)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        csi[k][j][i].x*csi[k][j][i].x +
                        csi[k][j][i].y*csi[k][j][i].y +
                        csi[k][j][i].z*csi[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k][j][i-1] = d * mesh->boundaryT.iLval + lt[k][j][i];
                }
                // fixedGradient boundary condition on i-right patch
                if (mesh->boundaryT.iRight=="fixedGradient" && i==mx-2)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        csi[k][j][i].x*csi[k][j][i].x +
                        csi[k][j][i].y*csi[k][j][i].y +
                        csi[k][j][i].z*csi[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k][j][i+1] = d * mesh->boundaryT.iRval + lt[k][j][i];
                }
                // fixedGradient boundary condition on j-left patch
                if (mesh->boundaryT.jLeft=="fixedGradient" && j==1)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        eta[k][j][i].x*eta[k][j][i].x +
                        eta[k][j][i].y*eta[k][j][i].y +
                        eta[k][j][i].z*eta[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k][j-1][i] = d * mesh->boundaryT.jLval + lt[k][j][i];
                }
                // fixedGradient boundary condition on j-right patch
                if (mesh->boundaryT.jRight=="fixedGradient" && j==my-2)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        eta[k][j][i].x*eta[k][j][i].x +
                        eta[k][j][i].y*eta[k][j][i].y +
                        eta[k][j][i].z*eta[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k][j+1][i] = d * mesh->boundaryT.jRval + lt[k][j][i];
                }
                // fixedGradient boundary condition on k-left patch
                if (mesh->boundaryT.kLeft=="fixedGradient" && k==1)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        zet[k][j][i].x*zet[k][j][i].x +
                        zet[k][j][i].y*zet[k][j][i].y +
                        zet[k][j][i].z*zet[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k-1][j][i] = d * mesh->boundaryT.kLval + lt[k][j][i];
                }
                // fixedGradient boundary condition on k-right patch
                if (mesh->boundaryT.kRight=="fixedGradient" && k==mz-2)
                {
                    // get the area at the boundary cell
                    PetscReal area
                    =
                    sqrt
                    (
                        zet[k][j][i].x*zet[k][j][i].x +
                        zet[k][j][i].y*zet[k][j][i].y +
                        zet[k][j][i].z*zet[k][j][i].z
                    );

                    PetscReal d = (1.0/aj[k][j][i])/area;

                    t[k+1][j][i] = d * mesh->boundaryT.kRval + lt[k][j][i];
                }

                // periodic boundary condition on i-left patch
                if (mesh->boundaryT.iLeft=="periodic" && i==1)
                {
                    if(mesh->i_periodic)       t[k][j][i-1] = lt[k][j][mx-2];
                    else if(mesh->ii_periodic) t[k][j][i-1] = lt[k][j][-2];
                }
                // periodic boundary condition on i-right patch
                if (mesh->boundaryT.iRight=="periodic" && i==mx-2)
                {
                    if(mesh->i_periodic)        t[k][j][i+1] = lt[k][j][1];
                    else if (mesh->ii_periodic) t[k][j][i+1] = lt[k][j][mx+1];
                }
                // periodic boundary condition on j-left patch
                if (mesh->boundaryT.jLeft=="periodic" && j==1)
                {
                    if(mesh->j_periodic)       t[k][j-1][i] = lt[k][my-2][i];
                    else if(mesh->jj_periodic) t[k][j-1][i] = lt[k][-2][i];
                }
                // periodic boundary condition on j-right patch
                if (mesh->boundaryT.jRight=="periodic" && j==my-2)
                {
                    if(mesh->j_periodic)       t[k][j+1][i] = lt[k][1][i];
                    else if(mesh->jj_periodic) t[k][j+1][i] = lt[k][my+1][i];
                }
                // periodic boundary condition on k-left patch
                if (mesh->boundaryT.kLeft=="periodic" && k==1)
                {
                    if(mesh->k_periodic)       t[k-1][j][i] = lt[mz-2][j][i];
                    else if(mesh->kk_periodic) t[k-1][j][i] = lt[-2][j][i];
                }
                // periodic boundary condition on k-right patch
                if (mesh->boundaryT.kRight=="periodic" && k==mz-2)
                {
                    if(mesh->k_periodic)       t[k+1][j][i] = lt[1][j][i];
                    else if(mesh->kk_periodic) t[k+1][j][i] = lt[mz+1][j][i];
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);
    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(da,  mesh->lIAj,   &iaj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);

    DMDAVecRestoreArray(da, teqn->lTmprt, &lt);
    DMDAVecRestoreArray(da, teqn->Tmprt,  &t);

    // scatter global to local
    DMGlobalToLocalBegin(da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);
    DMGlobalToLocalEnd  (da, teqn->Tmprt, INSERT_VALUES, teqn->lTmprt);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateNutBCs(les_ *les)
{
    mesh_          *mesh = les->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/nut";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***nut, ***nvert, ***meshTag, ***iaj;
    Cmpnts        ***icsi, ***cent;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, les->lNu_t, &nut);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(da,  mesh->lIAj,   &iaj);
    DMDAVecGetArray(fda, mesh->lCent,  &cent);

    // read inflow if necessary
    if(mesh->boundaryNut.kLeft == "inletFunction")
    {
        inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

        // read the inflow data if necessary
        if (ifPtr->typeNut == 3 || ifPtr->typeNut == 4)
        {
            readInflowNut(ifPtr, mesh->access->clock);

            if(ifPtr->typeU == 4)
            {
                // for each y coordinate, find the intex for the right shifted interpolation point
                if(ifPtr->shift2)
                {
                    PetscReal refRatio = ifPtr->shiftSpeed*les->access->clock->time / mesh->bounds.Ly;
                    PetscReal refShift = (refRatio - floor(refRatio))*mesh->bounds.Ly;

                    for (i=lxs; i<lxe; i++)
                    {
                        PetscReal locRatio = (refShift + (cent[lzs][lys][i].y - mesh->bounds.ymin)) / mesh->bounds.Ly;
                        PetscReal locShift = (locRatio - floor(locRatio))*mesh->bounds.Ly + mesh->bounds.ymin;

                        PetscReal minDist = 1e20;
                        PetscInt  iClose  = 0;

                        for (PetscInt ii=1; ii<mx-1; ii++)
                        {
                            PetscReal dist = fabs(ifPtr->ycent[ii] - locShift);
                            if(dist < minDist)
                            {
                                minDist = dist;
                                iClose  = ii;
                            }
                        }

                        // make sure iClose is the right index
                        if(ifPtr->ycent[iClose] <= locShift)
                        {
                            iClose++;
                        }

                        // only set the ones belonging to this processor. This is the index that has to
                        // be sourced from the inflow data to apply the shift at index i
                        ifPtr->yIDs[i] = iClose;

                        // set the right weight
                        PetscReal delta    = ifPtr->ycent[ifPtr->yIDs[i]] - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        PetscReal dist     = locShift - ifPtr->ycent[ifPtr->yIDs[i]-1];
                        ifPtr->yWeights[i] = dist / delta;
                    }
                }
            }
        }
    }

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // set to zero at solid internal cells and skip
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    nut[k][j][i] = 0.0;
                    continue;
                }

                // special boundary condition where inflow is mapped from precursor,
                // contains a check to ensure that also velocity is mapped. Only
                // avilable for k-left patches. Values are applied at the ghost
                if(mesh->boundaryNut.kLeft=="inletFunction" && k==1)
                {
                    inletFunctionTypes *ifPtr = mesh->inletF.kLeft;

                    // power law profile
                    if (ifPtr->typeNut == 3)
                    {
                        // periodize inflow according to input

                        // compute period fraction (handle index = n case)
                        PetscInt jif = j % ifPtr->n1 == 0 ? ifPtr->n1 : j % ifPtr->n1;
                        PetscInt iif = i % ifPtr->n2 == 0 ? ifPtr->n2 : i % ifPtr->n2;

                        // index is less than nPrds times inflow points: have data
                        if
                        (
                            j<=ifPtr->n1*ifPtr->prds1 &&
                            i<=ifPtr->n2*ifPtr->prds2
                        )
                        {
                            nut[k-1][j][i] = ifPtr->nut_plane[jif][iif];
                        }
                        // index is more than nPrds times inflow points: extrapolate
                        else
                        {
                            // extrapolate along j
                            if(j>ifPtr->n1*ifPtr->prds1) jif = ifPtr->n1;

                            // extrapolate along i
                            if(i>ifPtr->n2*ifPtr->prds2) iif = ifPtr->n2;

                            nut[k-1][j][i] = ifPtr->nut_plane[jif][iif];
                        }
                    }
                    else if (ifPtr->typeNut == 4)
                    {
                        PetscReal nutGhost;

                        if(ifPtr->shift2)
                        {
                            PetscInt  iLeft    = ifPtr->yIDs[i]-1,
                                      iRight   = ifPtr->yIDs[i];

                            PetscReal nutLeft =
                                ifPtr->inflowWeights[j][iLeft][0] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iLeft][0].j][ifPtr->closestCells[j][iLeft][0].i] +
                                ifPtr->inflowWeights[j][iLeft][1] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iLeft][1].j][ifPtr->closestCells[j][iLeft][1].i] +
                                ifPtr->inflowWeights[j][iLeft][2] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iLeft][2].j][ifPtr->closestCells[j][iLeft][2].i] +
                                ifPtr->inflowWeights[j][iLeft][3] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iLeft][3].j][ifPtr->closestCells[j][iLeft][3].i];

                            PetscReal nutRight =
                                ifPtr->inflowWeights[j][iRight][0] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iRight][0].j][ifPtr->closestCells[j][iRight][0].i] +
                                ifPtr->inflowWeights[j][iRight][1] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iRight][1].j][ifPtr->closestCells[j][iRight][1].i] +
                                ifPtr->inflowWeights[j][iRight][2] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iRight][2].j][ifPtr->closestCells[j][iRight][2].i] +
                                ifPtr->inflowWeights[j][iRight][3] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][iRight][3].j][ifPtr->closestCells[j][iRight][3].i];

                            nutGhost   = (1.0-ifPtr->yWeights[i]) * nutLeft + ifPtr->yWeights[i]*nutRight;
                        }
                        else
                        {
                            nutGhost =
                                ifPtr->inflowWeights[j][i][0] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][i][0].j][ifPtr->closestCells[j][i][0].i] +
                                ifPtr->inflowWeights[j][i][1] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][i][1].j][ifPtr->closestCells[j][i][1].i] +
                                ifPtr->inflowWeights[j][i][2] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][i][2].j][ifPtr->closestCells[j][i][2].i] +
                                ifPtr->inflowWeights[j][i][3] *
                                ifPtr->nut_plane[ifPtr->closestCells[j][i][3].j][ifPtr->closestCells[j][i][3].i];
                        }

                        nut[k-1][j][i] = nutGhost;
                    }
                }

                // zeroGradient boundary condition on i-left patch
                if (mesh->boundaryNut.iLeft=="zeroGradient" && i==1)
                {
                    nut[k][j][i-1] = nut[k][j][i];
                }
                // zeroGradient boundary condition on i-right patch
                if (mesh->boundaryNut.iRight=="zeroGradient" && i==mx-2)
                {
                    nut[k][j][i+1] = nut[k][j][i];
                }
                // zeroGradient boundary condition on j-left patch
                if (mesh->boundaryNut.jLeft=="zeroGradient" && j==1)
                {
                    nut[k][j-1][i] = nut[k][j][i];
                }
                // zeroGradient boundary condition on j-right patch
                if (mesh->boundaryNut.jRight=="zeroGradient" && j==my-2)
                {
                    nut[k][j+1][i] = nut[k][j][i];
                }
                // zeroGradient boundary condition on k-left patch
                if (mesh->boundaryNut.kLeft=="zeroGradient" && k==1)
                {
                    nut[k-1][j][i] = nut[k][j][i];
                }
                // zeroGradient boundary condition on k-right patch
                if (mesh->boundaryNut.kRight=="zeroGradient" && k==mz-2)
                {
                    nut[k+1][j][i] = nut[k][j][i];
                }

                // fixedValue boundary condition on i-left patch
                if (mesh->boundaryNut.iLeft=="fixedValue" && i==1)
                {
                    nut[k][j][i-1] = mesh->boundaryNut.iLval;
                }
                // fixedValue boundary condition on i-right patch
                if (mesh->boundaryNut.iRight=="fixedValue" && i==mx-2)
                {
                    nut[k][j][i+1] = mesh->boundaryNut.iRval;
                }
                // fixedValue boundary condition on j-left patch
                if (mesh->boundaryNut.jLeft=="fixedValue" && j==1)
                {
                    nut[k][j-1][i] = mesh->boundaryNut.jLval;
                }
                // fixedValue boundary condition on j-right patch
                if (mesh->boundaryNut.jRight=="fixedValue" && j==my-2)
                {
                    nut[k][j+1][i] = mesh->boundaryNut.jRval;
                }
                // fixedValue boundary condition on k-left patch
                if (mesh->boundaryNut.kLeft=="fixedValue" && k==1)
                {
                    nut[k-1][j][i] = mesh->boundaryNut.kLval;
                }
                // fixedValue boundary condition on k-right patch
                if (mesh->boundaryNut.kRight=="fixedValue" && k==mz-2)
                {
                    nut[k+1][j][i] = mesh->boundaryNut.kRval;
                }

                // periodic boundary condition on i-left patch
                if (mesh->boundaryNut.iLeft=="periodic" && i==1)
                {
                    if(mesh->i_periodic)       nut[k][j][i-1] = nut[k][j][mx-2];
                    else if(mesh->ii_periodic) nut[k][j][i-1] = nut[k][j][-2];
                }
                // periodic boundary condition on i-right patch
                if (mesh->boundaryNut.iRight=="periodic" && i==mx-2)
                {
                    if(mesh->i_periodic)        nut[k][j][i+1] = nut[k][j][1];
                    else if (mesh->ii_periodic) nut[k][j][i+1] = nut[k][j][mx+1];
                }
                // periodic boundary condition on j-left patch
                if (mesh->boundaryNut.jLeft=="periodic" && j==1)
                {
                    if(mesh->j_periodic)       nut[k][j-1][i] = nut[k][my-2][i];
                    else if(mesh->jj_periodic) nut[k][j-1][i] = nut[k][-2][i];
                }
                // periodic boundary condition on j-right patch
                if (mesh->boundaryNut.jRight=="periodic" && j==my-2)
                {
                    if(mesh->j_periodic)       nut[k][j+1][i] = nut[k][1][i];
                    else if(mesh->jj_periodic) nut[k][j+1][i] = nut[k][my+1][i];
                }
                // periodic boundary condition on k-left patch
                if (mesh->boundaryNut.kLeft=="periodic" && k==1)
                {
                    if(mesh->k_periodic)       nut[k-1][j][i] = nut[mz-2][j][i];
                    else if(mesh->kk_periodic) nut[k-1][j][i] = nut[-2][j][i];
                }
                // periodic boundary condition on k-right patch
                if (mesh->boundaryNut.kRight=="periodic" && k==mz-2)
                {
                    if(mesh->k_periodic)       nut[k+1][j][i] = nut[1][j][i];
                    else if(mesh->kk_periodic) nut[k+1][j][i] = nut[mz+1][j][i];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, les->lNu_t, &nut);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(da,  mesh->lIAj,   &iaj);
    DMDAVecRestoreArray(fda, mesh->lCent,  &cent);

    // scatter nut from global to local
    DMLocalToLocalBegin(da, les->lNu_t, INSERT_VALUES, les->lNu_t);
    DMLocalToLocalEnd  (da, les->lNu_t, INSERT_VALUES, les->lNu_t);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdatektBCs(les_ *les)
{
    mesh_          *mesh = les->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/kt";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***kt, ***nvert, ***meshTag;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, les->lk_t, &kt);
    DMDAVecGetArray(da, mesh->lNvert, &nvert);
    DMDAVecGetArray(da, mesh->lmeshTag, &meshTag);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                // set to zero at solid internal cells and skip
                if(isIBMSolidCell(k, j, i, nvert) || isZeroedCell(k, j, i, meshTag))
                {
                    kt[k][j][i] = 0.0;
                    continue;
                }

                // periodic boundary condition on i-left patch
                if (mesh->boundaryNut.iLeft=="periodic" && i==1)
                {
                    if(mesh->i_periodic)       kt[k][j][i-1] = kt[k][j][mx-2];
                    else if(mesh->ii_periodic) kt[k][j][i-1] = kt[k][j][-2];
                }
                else if (i==1)
                {
                    kt[k][j][i-1] = kt[k][j][i];
                }

                // periodic boundary condition on i-right patch
                if (mesh->boundaryNut.iRight=="periodic" && i==mx-2)
                {
                    if(mesh->i_periodic)        kt[k][j][i+1] = kt[k][j][1];
                    else if (mesh->ii_periodic) kt[k][j][i+1] = kt[k][j][mx+1];
                }
                else if (i==mx-2)
                {
                    kt[k][j][i+1] = kt[k][j][i];
                }

                // periodic boundary condition on j-left patch
                if (mesh->boundaryNut.jLeft=="periodic" && j==1)
                {
                    if(mesh->j_periodic)       kt[k][j-1][i] = kt[k][my-2][i];
                    else if(mesh->jj_periodic) kt[k][j-1][i] = kt[k][-2][i];
                }
                else if (j==1)
                {
                    kt[k][j-1][i] = kt[k][j][i];
                }

                // periodic boundary condition on j-right patch
                if (mesh->boundaryNut.jRight=="periodic" && j==my-2)
                {
                    if(mesh->j_periodic)       kt[k][j+1][i] = kt[k][1][i];
                    else if(mesh->jj_periodic) kt[k][j+1][i] = kt[k][my+1][i];
                }
                else if (j==my-2)
                {
                    kt[k][j+1][i] = kt[k][j][i];
                }

                // periodic boundary condition on k-left patch
                if (mesh->boundaryNut.kLeft=="periodic" && k==1)
                {
                    if(mesh->k_periodic)       kt[k-1][j][i] = kt[mz-2][j][i];
                    else if(mesh->kk_periodic) kt[k-1][j][i] = kt[-2][j][i];
                }
                else if (k==1)
                {
                    kt[k-1][j][i] = kt[k][j][i];
                }

                // periodic boundary condition on k-right patch
                if (mesh->boundaryNut.kRight=="periodic" && k==mz-2)
                {
                    if(mesh->k_periodic)       kt[k+1][j][i] = kt[1][j][i];
                    else if(mesh->kk_periodic) kt[k+1][j][i] = kt[mz+1][j][i];
                }
                else if (k==mz-2)
                {
                    kt[k+1][j][i] = kt[k][j][i];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, les->lk_t, &kt);
    DMDAVecRestoreArray(da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da, mesh->lmeshTag, &meshTag);

    // scatter nut from local to local
    DMLocalToLocalBegin(da, les->lk_t, INSERT_VALUES, les->lk_t);
    DMLocalToLocalEnd  (da, les->lk_t, INSERT_VALUES, les->lk_t);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdatePressureBCs(peqn_ *peqn)
{
    mesh_          *mesh = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/P";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***p, ***lp;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lP, &lp);
    DMDAVecGetArray(da, peqn->P, &p);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                PetscInt a=i, b=j, c=k, flag=0;

                if(i==0)
                {
                    if(mesh->i_periodic)       a=mx-2, flag=1;
                    else if(mesh->ii_periodic) a=-2, flag=1;
                    else                      a=1, flag=1;
                }
                if(i==mx-1)
                {
                    if(mesh->i_periodic)       a=1, flag=1;
                    else if(mesh->ii_periodic) a=mx+1, flag=1;
                    else                      a=mx-2, flag=1;
                }
                if(j==0)
                {
                    if(mesh->j_periodic)       b=my-2, flag=1;
                    else if(mesh->jj_periodic) b=-2, flag=1;
                    else                      b=1, flag=1;
                }
                if(j==my-1)
                {
                    if(mesh->j_periodic)       b=1, flag=1;
                    else if(mesh->jj_periodic) b=my+1, flag=1;
                    else                      b=my-2, flag=1;
                }
                if(k==0)
                {
                    if(mesh->k_periodic)       c=mz-2, flag=1;
                    else if(mesh->kk_periodic) c=-2, flag=1;
                    else                      c=1, flag=1;
                }
                if(k==mz-1)
                {
                    if(mesh->k_periodic)       c=1, flag=1;
                    else if(mesh->kk_periodic) c=mz+1, flag=1;
                    else                      c=mz-2, flag=1;
                }

                if(flag)
                {
                    p[k][j][i] = lp[c][b][a];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, peqn->lP, &lp);
    DMDAVecRestoreArray(da, peqn->P, &p);

    // scatter Phi from global to local
    DMGlobalToLocalBegin(da, peqn->P, INSERT_VALUES, peqn->lP);
    DMGlobalToLocalEnd  (da, peqn->P, INSERT_VALUES, peqn->lP);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdatePhiBCs(peqn_ *peqn)
{
    mesh_          *mesh = peqn->access->mesh;
    DM            da   = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs   = info.xs, xe = info.xs + info.xm;
    PetscInt      ys   = info.ys, ye = info.ys + info.ym;
    PetscInt      zs   = info.zs, ze = info.zs + info.zm;
    PetscInt      mx   = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/Phi";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscReal     ***phi, ***lphi;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da, peqn->lPhi, &lphi);
    DMDAVecGetArray(da, peqn->Phi, &phi);

    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                PetscInt a=i, b=j, c=k, flag=0;

                if(i==0)
                {
                    if(mesh->i_periodic)       a=mx-2, flag=1;
                    else if(mesh->ii_periodic) a=-2, flag=1;
                    else                      a=1, flag=1;
                }
                if(i==mx-1)
                {
                    if(mesh->i_periodic)       a=1, flag=1;
                    else if(mesh->ii_periodic) a=mx+1, flag=1;
                    else                      a=mx-2, flag=1;
                }
                if(j==0)
                {
                    if(mesh->j_periodic)       b=my-2, flag=1;
                    else if(mesh->jj_periodic) b=-2, flag=1;
                    else                      b=1, flag=1;
                }
                if(j==my-1)
                {
                    if(mesh->j_periodic)       b=1, flag=1;
                    else if(mesh->jj_periodic) b=my+1, flag=1;
                    else                      b=my-2, flag=1;
                }
                if(k==0)
                {
                    if(mesh->k_periodic)       c=mz-2, flag=1;
                    else if(mesh->kk_periodic) c=-2, flag=1;
                    else                      c=1, flag=1;
                }
                if(k==mz-1)
                {
                    if(mesh->k_periodic)       c=1, flag=1;
                    else if(mesh->kk_periodic) c=mz+1, flag=1;
                    else                      c=mz-2, flag=1;
                }

                if(flag)
                {
                    phi[k][j][i] = lphi[c][b][a];
                }
            }
        }
    }

    DMDAVecRestoreArray(da, peqn->lPhi, &lphi);
    DMDAVecRestoreArray(da, peqn->Phi, &phi);

    // scatter Phi from global to local
    DMGlobalToLocalBegin(da, peqn->Phi, INSERT_VALUES, peqn->lPhi);
    DMGlobalToLocalEnd  (da, peqn->Phi, INSERT_VALUES, peqn->lPhi);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateWallModelsU(ueqn_ *ueqn)
{
    mesh_         *mesh  = ueqn->access->mesh;
    flags_        *flags = ueqn->access->flags;
    clock_        *clock = ueqn->access->clock;
    teqn_         *teqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/WallModels";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscScalar   ***aj, ***nvert, ***meshTag, ***t;
    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***jeta;
    Cmpnts        ***ucat;

    Cmpnts        e1, e2, e3;                          // local wall normal co-ordinate system

    //rotation tensor direction cosines
    PetscReal        a11, a21, a31,
                     a12, a22, a32,
                     a13, a23, a33;

    PetscReal       tau11, tau21, tau31,
                    tau12, tau22, tau32,
                    tau13, tau23, tau33;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // only applies if nut is defined
    if(!flags->isLesActive) return(0);

    // get teqn pointer if teqn is active
    if(flags->isTeqnActive)
    {
        teqn = ueqn->access->teqn;
    }

    // compute averaging weights
    PetscReal aN = (PetscReal)(ueqn->access->clock->it);
    PetscReal m1 = aN  / (aN + 1.0);
    PetscReal m2 = 1.0 / (aN + 1.0);

    // this function applies boundary conditions for those cases in which
    // the wall shear stress tensor is specified.

    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);

    if(flags->isTeqnActive)
    {
        DMDAVecGetArray(da, teqn->lTmprt,  &t);
    }

    // wall function -3: opposite mechanism, whatever value of nut is applied at the wall
    // it doesen't matter since the whole viscous flux at the wall is modeled. The velocity
    // boundary condition it is just necessary to set a reasonable value at the ghost nodes,
    // which is calculated by preserving the wall normal gradient.
    if(mesh->boundaryU.jLWF==-3)
    {
        Shumann *wm = ueqn->jLWM->wmShumann;

        // compute some averaged global parameters first depending on the wall model
        PetscReal UParallelMeanMag  = 0.0;
        PetscReal dist              = 0.0;
        PetscInt  nCells            = 0;
        PetscReal frictionVel       = 0.0;
        PetscReal phiM, L, surfaceL;

        PetscReal lUParallelMeanMag = 0.0;
        PetscReal ldist             = 0.0;
        PetscInt  lnCells           = 0;

        PetscInt  computeqWall      = 0;
        PetscReal qWall             = 0.0;
        PetscReal lqWall            = 0.0;

        // check if must compute deltaTheta (if teqn is not active or another bc is used on T, qWall is 0.0 - neutral)
        if(flags->isTeqnActive)
        {
            // check if same wall function is used
            if(mesh->boundaryT.jLWF == -3 || mesh->boundaryT.jLWF == -2 || mesh->boundaryT.jLWF == -4)
            {
                computeqWall = 1;
            }
        }

        // compute average friction velocity
        if(wm->wfEvalType=="averaged")
        {
            // compute qWall (otherwise neutral wallfunction is used)
            if(computeqWall)
            {
                // compute average deltaTheta
                if(ys==0)
                {
                    for (k=lzs; k<lze; k++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            lqWall += teqn->jLWM->qWall.z[k-zs][i-xs];
                        }
                    }
                }

                MPI_Allreduce(&lqWall, &qWall, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            }

            // select processors next to the j-left patch
            if(ys==0)
            {
                j=1;

                // loop over boundary cells and compute UParallelMeanMag
                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // get cell velocity
                        Cmpnts    Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                        {
                            // increment velocity
                            lUParallelMeanMag += nMag(UcellParallel);

                            // increment distance
                            ldist += s;

                            // increment cell count
                            lnCells++;
                        }
                    }
                }
            }

            MPI_Allreduce(&lUParallelMeanMag, &UParallelMeanMag, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldist, &dist, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnCells, &nCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            UParallelMeanMag = UParallelMeanMag / nCells;
            qWall            = qWall / nCells;
            dist             = dist / nCells;

            uStarShumann
            (
                UParallelMeanMag, dist, wm->roughness,
                wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                frictionVel, phiM, L, k, j, i
            );

            // print information (debugging)
            // PetscPrintf(mesh->MESH_COMM, "ShumannGrotzbach U: uStar = %lf, <U_||> = %lf, L = %lf, phiM = %lf, qWall = %lf\n", frictionVel, UParallelMeanMag, L, phiM,qWall);
        }

        // initialize filtered velocity
        if(!ueqn->jLWM->uFiltSet)
        {
            if(ys==0)
            {
                j=1;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // get cell normals
                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // velocity at point b
                        Cmpnts Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        ueqn->jLWM->uFilt.x[k-zs][i-xs] = UcellParallel.x;
                        ueqn->jLWM->uFilt.y[k-zs][i-xs] = UcellParallel.y;
                        ueqn->jLWM->uFilt.z[k-zs][i-xs] = UcellParallel.z;
                    }
                }
            }

            ueqn->jLWM->uFiltSet = 1;
        }

        // compute the wall shear stress

        // select processors next to the j-left patch
        if(ys==0)
        {
            j=1;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // get the cell area
                    PetscReal area = nMag(eta[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // velocity at point b
                    Cmpnts Ucell = nSet(ucat[k][j][i]);

                    // compute wall-normal velocity
                    Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                    //local co-ordinate system
                    e1 = nUnit(UcellParallel);

                    e3 = nSet(n);

                    e2 = nCross(e3, e1);

                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;
                    // average velocity
                    ueqn->jLWM->uMeanMag[k_wm][i_wm] = m1 * ueqn->jLWM->uMeanMag[k_wm][i_wm]  + m2 * nMag(UcellParallel);

                    // compute filter scale
                    Cmpnts        l = nSetFromComponents(1.0 / aj[k][j][i] / nMag(zet[k][j][i]), 1.0 / aj[k][j][i] / nMag(csi[k][j][i]), 0.0);
                    PetscReal     T = std::max(nMag(l) / ueqn->jLWM->uMeanMag[k_wm][i_wm], clock->dt);

                    // filter velocity
                    ueqn->jLWM->uFilt.x[k_wm][i_wm] = (1.0 - clock->dt/T) * ueqn->jLWM->uFilt.x[k_wm][i_wm] + (clock->dt/T) * UcellParallel.x;
                    ueqn->jLWM->uFilt.y[k_wm][i_wm] = (1.0 - clock->dt/T) * ueqn->jLWM->uFilt.y[k_wm][i_wm] + (clock->dt/T) * UcellParallel.y;

                    // compute local friction velocity
                    if(wm->wfEvalType=="localized")
                    {
                        // compute qWall (otherwise neutral wallfunction is used)
                        if(computeqWall)
                        {
                            qWall = teqn->jLWM->qWall.z[k-zs][i-xs];
                        }

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        UParallelMeanMag = nMag(nSetFromComponents(ueqn->jLWM->uFilt.x[k_wm][i_wm], ueqn->jLWM->uFilt.y[k_wm][i_wm], 0.0));

                        uStarShumann
                        (
                            UParallelMeanMag, s, wm->roughness,
                            wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                            frictionVel, phiM, L, k, j, i
                        );
                    }

                    //create the transformation vector for rotation to global axis.
                    a11 = e1.x; a12 = e2.x, a13 = e3.x;
                    a21 = e1.y; a22 = e2.y, a23 = e3.y;
                    a31 = e1.z; a32 = e2.z, a33 = e3.z;

                    //transform it to original co-ordinate system
                    tau11 = -2.0 * a11 * a13 * frictionVel*frictionVel;
                    tau12 = -(a13 * a21 + a11 * a23) * frictionVel*frictionVel;
                    tau13 = -(a13 * a31 + a11 * a33) * frictionVel*frictionVel;

                    tau21 = -(a23 * a11 + a21 * a13) * frictionVel*frictionVel;
                    tau22 = -2.0 * a23 * a21 * frictionVel*frictionVel;
                    tau23 = -(a23 * a31 + a21 * a33) * frictionVel*frictionVel;

                    tau31 = -(a33 * a11 + a31 * a13) * frictionVel*frictionVel;
                    tau32 = -(a33 * a21 + a31 * a23) * frictionVel*frictionVel;
                    tau33 = -2.0 * a33 * a31 * frictionVel*frictionVel;

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        //printf("uFilt.x = %f, u.x = %f\n", ueqn->jLWM->uFilt.x[k_wm][i_wm], UcellParallel.x);

                        if(mesh->meshFileType == "curvilinear")
                        {
                            ueqn->jLWM->tauWall.x[k_wm][i_wm] = (tau11* jeta[k][j-1][i].x + tau12 * jeta[k][j-1][i].y + tau13 * jeta[k][j-1][i].z);
                            ueqn->jLWM->tauWall.y[k_wm][i_wm] = (tau21* jeta[k][j-1][i].x + tau22 * jeta[k][j-1][i].y + tau23 * jeta[k][j-1][i].z);
                            ueqn->jLWM->tauWall.z[k_wm][i_wm] = (tau31* jeta[k][j-1][i].x + tau32 * jeta[k][j-1][i].y + tau33 * jeta[k][j-1][i].z);
                        }
                        else if(mesh->meshFileType == "cartesian")
                        {
                            PetscReal TauXZ = - frictionVel*frictionVel * (ueqn->jLWM->uFilt.x[k_wm][i_wm] / PetscMax(UParallelMeanMag, 1e-5));
                            PetscReal TauYZ = - frictionVel*frictionVel * (ueqn->jLWM->uFilt.y[k_wm][i_wm] / PetscMax(UParallelMeanMag, 1e-5));

                            // transform to contravariant coords
                            ueqn->jLWM->tauWall.x[k_wm][i_wm] = jeta[k][j-1][i].z * TauXZ;
                            ueqn->jLWM->tauWall.y[k_wm][i_wm] = jeta[k][j-1][i].z * TauYZ;
                            ueqn->jLWM->tauWall.z[k_wm][i_wm] = jeta[k][j-1][i].x * TauXZ + jeta[k][j-1][i].y * TauYZ;
                        }
                    }
                }
            }
        }
    }

    if(mesh->boundaryU.jRWF==-3)
    {
        Shumann *wm = ueqn->jRWM->wmShumann;

        // compute some averaged global parameters first depending on the wall model
        PetscReal UParallelMeanMag  = 0.0;
        PetscReal dist              = 0.0;
        PetscInt  nCells            = 0;
        PetscReal frictionVel       = 0.0;
        PetscReal phiM, L;

        PetscReal lUParallelMeanMag = 0.0;
        PetscReal ldist             = 0.0;
        PetscInt  lnCells           = 0;

        PetscInt  computeqWall      = 0;
        PetscReal qWall             = 0.0;
        PetscReal lqWall            = 0.0;

        // check if must compute deltaTheta (if teqn is not active or another bc is used on T, qWall is 0.0 - neutral)
        if(flags->isTeqnActive)
        {
            // check if same wall function is used
            if(mesh->boundaryT.jRWF == -3 || mesh->boundaryT.jRWF == -2)
            {
                computeqWall = 1;
            }
        }

        // compute average friction velocity
        if(wm->wfEvalType=="averaged")
        {
            // compute qWall (otherwise neutral wallfunction is used)
            if(computeqWall)
            {
                // compute average deltaTheta
                if(ye==my)
                {
                    for (k=lzs; k<lze; k++)
                    {
                        for (i=lxs; i<lxe; i++)
                        {
                            lqWall += teqn->jRWM->qWall.z[k-zs][i-xs];
                        }
                    }
                }

                MPI_Allreduce(&lqWall, &qWall, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            }

            // select processors next to the j-right patch
            if(ye==my)
            {
                j=my-2;

                // loop over boundary cells and compute UParallelMeanMag
                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // get cell normals and flip them
                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // velocity at point b
                        Cmpnts Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                        {
                            // increment velocity
                            lUParallelMeanMag += nMag(UcellParallel);

                            // increment distance
                            ldist += s;

                            // increment cell count
                            lnCells++;
                        }
                    }
                }
            }

            MPI_Allreduce(&lUParallelMeanMag, &UParallelMeanMag, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldist, &dist, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnCells, &nCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            UParallelMeanMag = UParallelMeanMag / nCells;
            qWall            = qWall / nCells;
            dist             = dist / nCells;

            uStarShumann
            (
                UParallelMeanMag, dist, wm->roughness,
                wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                frictionVel, phiM, L, k, j, i
            );
            // print information (debugging)
            // PetscPrintf(PETSC_COMM_WORLD, "ShumannGrotzbach: uStar = %lf, <U_||> = %lf\n", frictionVel, UParallelMeanMag);
        }

        // initialize filtered velocity
        if(!ueqn->jRWM->uFiltSet)
        {
            if(ys==mx)
            {
                j=my-2;

                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // get cell normals
                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // velocity at point b
                        Cmpnts Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        ueqn->jRWM->uFilt.x[k-zs][i-xs] = UcellParallel.x;
                        ueqn->jRWM->uFilt.y[k-zs][i-xs] = UcellParallel.y;
                        ueqn->jRWM->uFilt.z[k-zs][i-xs] = UcellParallel.z;
                    }
                }
            }

            ueqn->jRWM->uFiltSet = 1;
        }

        // compute the wall shear stress

        // select processors next to the j-left patch
        if(ye==my)
        {
            j=my-2;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    // get the cell area
                    PetscReal area = nMag(eta[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // velocity at point b
                    Cmpnts Ucell = nSet(ucat[k][j][i]);

                    // compute wall-normal velocity
                    Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;

                    // average velocity
                    ueqn->jRWM->uMeanMag[k_wm][i_wm] = m1 * ueqn->jRWM->uMeanMag[k_wm][i_wm]  + m2 * nMag(UcellParallel);

                    // compute filter scale
                    Cmpnts        l = nSetFromComponents(1.0 / aj[k][j][i] / nMag(zet[k][j][i]), 1.0 / aj[k][j][i] / nMag(csi[k][j][i]), 0.0);
                    PetscReal     T = nMag(l) / ueqn->jRWM->uMeanMag[k_wm][i_wm];

                    // filter velocity
                    ueqn->jRWM->uFilt.x[k_wm][i_wm] = (1.0 - clock->dt/T) * ueqn->jRWM->uFilt.x[k_wm][i_wm] + (clock->dt/T) * UcellParallel.x;
                    ueqn->jRWM->uFilt.y[k_wm][i_wm] = (1.0 - clock->dt/T) * ueqn->jRWM->uFilt.y[k_wm][i_wm] + (clock->dt/T) * UcellParallel.y;

                    // compute local friction velocity
                    if(wm->wfEvalType=="localized")
                    {
                        // compute qWall (otherwise neutral wallfunction is used)
                        if(computeqWall)
                        {
                            qWall = teqn->jRWM->qWall.z[k-zs][i-xs];
                        }

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        UParallelMeanMag = nMag(nSetFromComponents(ueqn->jRWM->uFilt.x[k_wm][i_wm], ueqn->jRWM->uFilt.y[k_wm][i_wm], 0.0));

                        uStarShumann
                        (
                            UParallelMeanMag, s, wm->roughness,
                            wm->gammaM, wm->kappa, qWall, wm->thetaRef,
                            frictionVel, phiM, L, k, j, i
                        );
                    }

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        // wall shear stress in cartesian coords
                        PetscReal TauXZ = - frictionVel*frictionVel * (ueqn->jRWM->uFilt.x[k_wm][i_wm] / PetscMax(UParallelMeanMag, 1e-5));
                        PetscReal TauYZ = - frictionVel*frictionVel * (ueqn->jRWM->uFilt.y[k_wm][i_wm] / PetscMax(UParallelMeanMag, 1e-5));

                        // wall shear stress in curvilinear coords
                        ueqn->jRWM->tauWall.x[k_wm][i_wm] = jeta[k][j][i].z * TauXZ;
                        ueqn->jRWM->tauWall.y[k_wm][i_wm] = jeta[k][j][i].z * TauYZ;
                        ueqn->jRWM->tauWall.z[k_wm][i_wm] = jeta[k][j][i].x * TauXZ + jeta[k][j][i].y * TauYZ;
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);

    if(flags->isTeqnActive)
    {
        DMDAVecRestoreArray(da, teqn->lTmprt,  &t);
    }

    return(0);
}

//***************************************************************************************************************//


PetscErrorCode UpdateWallModelsT(teqn_ *teqn)
{
    mesh_        *mesh  = teqn->access->mesh;
    ueqn_        *ueqn  = teqn->access->ueqn;
    clock_       *clock = teqn->access->clock;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;
    PetscInt      xs    = info.xs, xe = info.xs + info.xm;
    PetscInt      ys    = info.ys, ye = info.ys + info.ym;
    PetscInt      zs    = info.zs, ze = info.zs + info.zm;
    PetscInt      mx    = info.mx, my = info.my, mz = info.mz;

    word          typeName = "boundary/WallModels";

    PetscInt      lxs, lxe, lys, lye, lzs, lze;
    PetscInt      i, j, k;

    PetscInt      mx_wm = xe-xs;
    PetscInt      my_wm = ye-ys;
    PetscInt      mz_wm = ze-zs;

    PetscScalar   ***aj, ***nvert, ***meshTag, ***t;
    Cmpnts        ***csi, ***eta, ***zet;
    Cmpnts        ***jeta;
    Cmpnts        ***ucat;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // only applies if nut is defined
    if(!teqn->access->flags->isLesActive) return(0);

    // this function wall heat flux to be added in the teqn

    DMDAVecGetArray(fda, mesh->lCsi,   &csi);
    DMDAVecGetArray(fda, mesh->lEta,   &eta);
    DMDAVecGetArray(fda, mesh->lZet,   &zet);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);
    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecGetArray(da,  teqn->lTmprt, &t);

    if(mesh->boundaryT.jLWF==-3)
    {
        Shumann *wm = teqn->jLWM->wmShumann;

        PetscInt updateTemp = 0;

        if(clock->time != wm->tLast)
        {
            updateTemp = 1;
            wm->tLast = clock->time;
        }

        // set surface temperature
        if(!wm->surfaceThetaSet)
        {
            if(ys==0)
            {
                j=1;

                PetscMalloc(sizeof(PetscReal*)*mz_wm, &(wm->surfaceTheta));

                for (k=lzs; k<lze; k++)
                {
                    PetscMalloc(sizeof(PetscReal)*mx_wm, &(wm->surfaceTheta[k-zs]));

                    for (i=lxs; i<lxe; i++)
                    {
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j-1][i]);
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j-1][i]);
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j-1][i]);
                    }

                }
            }

            wm->surfaceThetaSet = 1;
        }

        // compute some averaged global parameters first depending on the wall model
        PetscReal UParallelMeanMag  = 0.0;
        PetscReal deltaTheta        = 0.0;
        PetscReal dist              = 0.0;
        PetscInt  nCells            = 0;
        PetscReal frictionVel       = 0.0;
        PetscReal phiM, phiH, L;
        PetscReal qWall;

        PetscReal lUParallelMeanMag = 0.0;
        PetscReal ldeltaTheta       = 0.0;
        PetscReal ldist             = 0.0;
        PetscInt  lnCells           = 0;

        // compute average friction velocity
        if(wm->wfEvalType=="averaged")
        {
            // select processors next to the j-left patch
            if(ys==0)
            {
                j=1;

                // loop over boundary cells and compute UParallelMeanMag
                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        PetscInt k_wm = k - zs;
                        PetscInt i_wm = i - xs;

                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // get cell temperature
                        PetscReal cellTheta = t[k][j][i];

                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // get cell velocity
                        Cmpnts    Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                        {
                            // update surface temperature
                            if(updateTemp) wm->surfaceTheta[k_wm][i_wm] += wm->heatingRate * clock->dt;

                            ldeltaTheta += (cellTheta - wm->surfaceTheta[k_wm][i_wm]);

                            // increment velocity
                            lUParallelMeanMag += nMag(UcellParallel);

                            // increment distance
                            ldist += s;

                            // increment cell count
                            lnCells++;
                        }
                    }
                }
            }

            MPI_Allreduce(&lUParallelMeanMag, &UParallelMeanMag, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldeltaTheta, &deltaTheta, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldist, &dist, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnCells, &nCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            UParallelMeanMag = UParallelMeanMag / nCells;
            deltaTheta       = deltaTheta / nCells;
            dist             = dist / nCells;

            qWallShumann
            (
                UParallelMeanMag, dist, wm->roughness,
                wm->gammaM, wm->gammaH, wm->alphaH,
                wm->thetaRef, deltaTheta, wm->kappa,
                qWall, frictionVel, phiM, phiH, L, k, j, i
            );

            // print information (debugging)
            // PetscPrintf(mesh->MESH_COMM, "ShumannGrotzbach T: uStar = %lf, <U_||> = %lf, L = %lf, phiM = %lf, phiH = %lf, deltaT = %lf, qWall = %lf\n", frictionVel, UParallelMeanMag, L, phiM, phiH, deltaTheta, qWall);
        }

        // compute the wall shear stress

        // select processors next to the j-left patch
        if(ys==0)
        {
            j=1;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;

                    // get the cell area
                    PetscReal area = nMag(eta[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // velocity at point b
                    Cmpnts Ucell = nSet(ucat[k][j][i]);

                    // get cell temperature
                    PetscReal cellTheta = t[k][j][i];

                    // compute wall-normal velocity
                    Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                    // compute local friction velocity
                    if(wm->wfEvalType=="localized")
                    {
                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // update surface temperature
                        if(updateTemp) wm->surfaceTheta[k_wm][i_wm] += wm->heatingRate * clock->dt;

                        deltaTheta = (cellTheta - wm->surfaceTheta[k_wm][i_wm]);

                        UParallelMeanMag = nMag(UcellParallel);

                        qWallShumann
                        (
                            UParallelMeanMag, s, wm->roughness,
                            wm->gammaM, wm->gammaH, wm->alphaH,
                            wm->thetaRef, deltaTheta, wm->kappa,
                            qWall, frictionVel, phiM, phiH, L, k, j, i
                        );
                    }

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        teqn->jLWM->qWall.x[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.y[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.z[k_wm][i_wm] = qWall;
                    }
                }
            }
        }
    }
    else if(mesh->boundaryT.jLWF==-2)
    {
        Shumann *wm = teqn->jLWM->wmShumann;

        // select processors next to the j-left patch
        if(ys==0)
        {
            j=1;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        teqn->jLWM->qWall.x[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.y[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.z[k_wm][i_wm] = wm->qWall;
                    }
                }
            }
        }
        
    }
    else if(mesh->boundaryT.jLWF==-4)
    {
        Shumann *wm  = teqn->jLWM->wmShumann;

        // compute some averaged global parameters first depending on the wall model
        PetscReal UParallelMeanMag  = 0.0;
        PetscReal deltaTheta        = 0.0;
        PetscReal dist              = 0.0;
        PetscInt  nCells            = 0;
        PetscReal frictionVel       = 0.0;
        PetscReal surfaceTemp;
        PetscReal phiM, phiH, L;
        PetscReal qWall;

        PetscReal lUParallelMeanMag = 0.0;
        PetscReal ldeltaTheta       = 0.0;
        PetscReal ldist             = 0.0;
        PetscInt  lnCells           = 0;

        // find interpolation weights
        PetscReal w[2];
        PetscInt  l[2];

        findInterpolationWeights(w, l, wm->timeVec, wm->numT, clock->time);

        //interpolate the surface temperature from closest available mesoscale times
        surfaceTemp = w[0] * wm->surfTemp[l[0]] + w[1] * wm->surfTemp[l[1]];

        // compute average friction velocity
        if(wm->wfEvalType=="averaged")
        {
            // select processors next to the j-left patch
            if(ys==0)
            {
                j=1;

                // loop over boundary cells and compute UParallelMeanMag
                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        PetscInt k_wm = k - zs;
                        PetscInt i_wm = i - xs;

                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // get cell temperature
                        PetscReal cellTheta = t[k][j][i];

                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // get cell velocity
                        Cmpnts    Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                        {
                            ldeltaTheta += (cellTheta - surfaceTemp);

                            // increment velocity
                            lUParallelMeanMag += nMag(UcellParallel);

                            // increment distance
                            ldist += s;

                            // increment cell count
                            lnCells++;
                        }
                    }
                }
            }

            MPI_Allreduce(&lUParallelMeanMag, &UParallelMeanMag, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldeltaTheta, &deltaTheta, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldist, &dist, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnCells, &nCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            UParallelMeanMag = UParallelMeanMag / nCells;
            deltaTheta       = deltaTheta / nCells;
            dist             = dist / nCells;

            qWallShumann
            (
                UParallelMeanMag, dist, wm->roughness,
                wm->gammaM, wm->gammaH, wm->alphaH,
                wm->thetaRef, deltaTheta, wm->kappa,
                qWall, frictionVel, phiM, phiH, L, k, j, i
            );

            // print information (debugging)
            // PetscPrintf(mesh->MESH_COMM, "ShumannGrotzbach T: uStar = %lf, <U_||> = %lf, L = %lf, phiM = %lf, phiH = %lf, deltaT = %lf, qWall = %lf\n", frictionVel, UParallelMeanMag, surfaceL, phiM, phiH, deltaTheta, qWall);
        }

        if(ys==0)
        {
            j=1;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;

                    // get the cell area
                    PetscReal area = nMag(eta[k][j][i]);

                    // get cell normals
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // velocity at point b
                    Cmpnts Ucell = nSet(ucat[k][j][i]);

                    // get cell temperature
                    PetscReal cellTheta = t[k][j][i];

                    // compute wall-normal velocity
                    Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                    // compute local friction velocity
                    if(wm->wfEvalType=="localized")
                    {
                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        deltaTheta = (cellTheta - surfaceTemp);

                        UParallelMeanMag = nMag(UcellParallel);

                        qWallShumann
                        (
                            UParallelMeanMag, s, wm->roughness,
                            wm->gammaM, wm->gammaH, wm->alphaH,
                            wm->thetaRef, deltaTheta, wm->kappa,
                            qWall, frictionVel, phiM, phiH, L, k, j, i
                        );
                    }

                    if(i  == 5 && k == 5 && !teqn->access->flags->isIBMActive)
                        PetscPrintf(PETSC_COMM_SELF, "surfTemp = %lf, L = %lf, bPtTemp = %lf, deltaTheta = %lf, qWall = %lf, ustar = %lf, utmag = %lf, walldist = %lf\n", surfaceTemp, L, cellTheta, deltaTheta, qWall, frictionVel, UParallelMeanMag, 0.5/aj[k][j][i]/area);

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        teqn->jLWM->qWall.x[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.y[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.z[k_wm][i_wm] = qWall;
                    }
                }
            }
        }
    }

    if(mesh->boundaryT.jRWF==-3)
    {
        Shumann *wm = teqn->jRWM->wmShumann;

        PetscInt updateTemp = 0;

        if(clock->time != wm->tLast)
        {
            updateTemp = 1;
            wm->tLast = clock->time;
        }

        // set surface temperature
        if(!wm->surfaceThetaSet)
        {
            if(ye==my)
            {
                j=my-2;

                PetscMalloc(sizeof(PetscReal*)*mz_wm, &(wm->surfaceTheta));

                for (k=lzs; k<lze; k++)
                {
                    PetscMalloc(sizeof(PetscReal)*mx_wm, &(wm->surfaceTheta[k-zs]));

                    for (i=lxs; i<lxe; i++)
                    {
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j+1][i]);
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j+1][i]);
                        wm->surfaceTheta[k-zs][i-xs] = 0.5*(t[k][j][i] + t[k][j+1][i]);
                    }

                }
            }

            wm->surfaceThetaSet = 1;
        }

        // compute some averaged global parameters first depending on the wall model
        PetscReal UParallelMeanMag  = 0.0;
        PetscReal deltaTheta        = 0.0;
        PetscReal dist              = 0.0;
        PetscInt  nCells            = 0;
        PetscReal frictionVel       = 0.0;
        PetscReal phiM, phiH, L;
        PetscReal qWall;

        PetscReal lUParallelMeanMag = 0.0;
        PetscReal ldeltaTheta       = 0.0;
        PetscReal ldist             = 0.0;
        PetscInt  lnCells           = 0;

        // compute average friction velocity
        if(wm->wfEvalType=="averaged")
        {
            // select processors next to the j-right patch
            if(ye==my)
            {
                j=my-2;

                // loop over boundary cells and compute UParallelMeanMag
                for (k=lzs; k<lze; k++)
                {
                    for (i=lxs; i<lxe; i++)
                    {
                        PetscInt k_wm = k - zs;
                        PetscInt i_wm = i - xs;

                        // get the cell area
                        PetscReal area = nMag(eta[k][j][i]);

                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // get cell temperature
                        PetscReal cellTheta = t[k][j][i];

                        // get cell normals and flip them
                        PetscReal ni[3], nj[3], nk[3];
                        calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                        // get j normal
                        Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                        // velocity at point b
                        Cmpnts Ucell = nSet(ucat[k][j][i]);

                        // compute wall-normal velocity
                        Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                        // compute wall-parallel velocity
                        Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                        if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                        {
                            // update surface temperature
                            if(updateTemp) wm->surfaceTheta[k_wm][i_wm] += wm->heatingRate * clock->dt;

                            ldeltaTheta += (cellTheta - wm->surfaceTheta[k_wm][i_wm]);

                            // increment velocity
                            lUParallelMeanMag += nMag(UcellParallel);

                            // increment distance
                            ldist += s;

                            // increment cell count
                            lnCells++;
                        }
                    }
                }
            }

            MPI_Allreduce(&lUParallelMeanMag, &UParallelMeanMag, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldeltaTheta, &deltaTheta, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&ldist, &dist, 1, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
            MPI_Allreduce(&lnCells, &nCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            UParallelMeanMag = UParallelMeanMag / nCells;
            deltaTheta       = deltaTheta / nCells;
            dist             = dist / nCells;

            qWallShumann
            (
                UParallelMeanMag, dist, wm->roughness,
                wm->gammaM, wm->gammaH, wm->alphaH,
                wm->thetaRef, deltaTheta, wm->kappa,
                qWall, frictionVel, phiM, phiH, L, k, j, i
            );

            // print information (debugging)
            // PetscPrintf(PETSC_COMM_WORLD, "ShumannGrotzbach: uStar = %lf, <U_||> = %lf\n", frictionVel, UParallelMeanMag);
        }

        // compute the wall shear stress

        // select processors next to the j-left patch
        if(ye==my)
        {
            j=my-2;

            for (k=lzs; k<lze; k++)
            {
                for (i=lxs; i<lxe; i++)
                {
                    PetscInt k_wm = k - zs;
                    PetscInt i_wm = i - xs;

                    // get the cell area
                    PetscReal area = nMag(eta[k][j][i]);

                    // get cell normals and flip them
                    PetscReal ni[3], nj[3], nk[3];
                    calculateNormal(csi[k][j][i], eta[k][j][i], zet[k][j][i], ni, nj, nk);

                    // get j normal
                    Cmpnts n = nSetFromComponents(nj[0], nj[1], nj[2]);

                    // velocity at point b
                    Cmpnts Ucell = nSet(ucat[k][j][i]);

                    // get cell temperature
                    PetscReal cellTheta = t[k][j][i];

                    // compute wall-normal velocity
                    Cmpnts UcellNormal = nScale(nDot(Ucell, n), n);

                    // compute wall-parallel velocity
                    Cmpnts UcellParallel = nSub(Ucell, UcellNormal);

                    // compute local friction velocity
                    if(wm->wfEvalType=="localized")
                    {
                        // distance from wall and velocity at point b
                        PetscReal s = 0.5/aj[k][j][i]/area;

                        // update surface temperature
                        if(updateTemp) wm->surfaceTheta[k_wm][i_wm] += wm->heatingRate * clock->dt;

                        deltaTheta = (cellTheta - wm->surfaceTheta[k_wm][i_wm]);

                        UParallelMeanMag = nMag(UcellParallel);

                        qWallShumann
                        (
                            UParallelMeanMag, s, wm->roughness,
                            wm->gammaM, wm->gammaH, wm->alphaH,
                            wm->thetaRef, deltaTheta, wm->kappa,
                            qWall, frictionVel, phiM, phiH, L, k, j, i
                        );
                    }

                    if (isFluidCell(k, j, i, nvert) && isCalculatedCell(k, j, i, meshTag))
                    {
                        teqn->jLWM->qWall.x[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.y[k_wm][i_wm] = 0.0;
                        teqn->jLWM->qWall.z[k_wm][i_wm] = qWall;
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCsi,   &csi);
    DMDAVecRestoreArray(fda, mesh->lEta,   &eta);
    DMDAVecRestoreArray(fda, mesh->lZet,   &zet);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(da,  mesh->lmeshTag, &meshTag);
    DMDAVecRestoreArray(fda, ueqn->lUcat,  &ucat);
    DMDAVecRestoreArray(da,  teqn->lTmprt, &t);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode UpdateImmersedBCs(ibm_ *ibm)
{
    mesh_         *mesh = ibm->access->mesh;
    ueqn_         *ueqn = ibm->access->ueqn;
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo info  = mesh->info;

    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    Cmpnts        ***lucat, ***icsi, ***jeta, ***kzet;
    Cmpnts        ***ucont, ***ucat;
    PetscReal     ***nvert, ***t;
    PetscReal     ucx, ucy, ucz;

    PetscInt      i, j, k;
    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(fda, ueqn->lUcat, &lucat);
    DMDAVecGetArray(fda, ueqn->Ucat, &ucat);
    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);
    DMDAVecGetArray(fda, mesh->lICsi, &icsi);
    DMDAVecGetArray(fda, mesh->lJEta, &jeta);
    DMDAVecGetArray(fda, mesh->lKZet, &kzet);

    for (PetscInt k=zs; k<lze; k++)
    {
        for (PetscInt j=ys; j<lye; j++)
        {
            for (PetscInt i=xs; i<lxe; i++)
            {

                if (isIBMFluidIFace(k, j, i, i+1, nvert))
                {
                    ucx = (lucat[k][j][i].x + lucat[k][j][i+1].x) * 0.5;
                    ucy = (lucat[k][j][i].y + lucat[k][j][i+1].y) * 0.5;
                    ucz = (lucat[k][j][i].z + lucat[k][j][i+1].z) * 0.5;

                    if(ibm->wallShearOn)
                    {
                        ucont[k][j][i].x = 0.0;
                    }
                    else 
                    {
                        ucont[k][j][i].x = (ucx * icsi[k][j][i].x + ucy * icsi[k][j][i].y + ucz * icsi[k][j][i].z);
                    }
                }

                if (isIBMFluidJFace(k, j, i, j+1, nvert))
                {
                    ucx = (lucat[k][j+1][i].x + lucat[k][j][i].x) * 0.5;
                    ucy = (lucat[k][j+1][i].y + lucat[k][j][i].y) * 0.5;
                    ucz = (lucat[k][j+1][i].z + lucat[k][j][i].z) * 0.5;
                    
                    if(ibm->wallShearOn)
                    {
                        ucont[k][j][i].y = 0.0;
                    }
                    else
                    {
                        ucont[k][j][i].y = (ucx * jeta[k][j][i].x + ucy * jeta[k][j][i].y + ucz * jeta[k][j][i].z);
                    }
                }

                if (isIBMFluidKFace(k, j, i, k+1, nvert))
                {
                    ucx = (lucat[k+1][j][i].x + lucat[k][j][i].x) * 0.5;
                    ucy = (lucat[k+1][j][i].y + lucat[k][j][i].y) * 0.5;
                    ucz = (lucat[k+1][j][i].z + lucat[k][j][i].z) * 0.5;

                    if(ibm->wallShearOn)
                    {
                        ucont[k][j][i].z = 0.0;
                    }
                    else
                    {
                        ucont[k][j][i].z = (ucx * kzet[k][j][i].x + ucy * kzet[k][j][i].y + ucz * kzet[k][j][i].z);
                    }
                }

                if(isIBMSolidIFace(k, j, i, i+1, nvert))
                {
                  ucont[k][j][i].x = 0;
                }

                if(isIBMSolidJFace(k, j, i, j+1, nvert))
                {
                  ucont[k][j][i].y = 0;
                }

                if(isIBMSolidKFace(k, j, i, k+1, nvert))
                {
                  ucont[k][j][i].z = 0;
                }

                if(isIBMSolidCell(k, j, i, nvert))
                {
                  mSetValue(ucat[k][j][i], 0);
                }

            }
        }
    }

    DMDAVecRestoreArray(da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    DMDAVecRestoreArray(fda, ueqn->Ucat, &ucat);
    DMDAVecRestoreArray(fda, ueqn->Ucont, &ucont);
    DMDAVecRestoreArray(fda, mesh->lICsi, &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta, &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet, &kzet);

    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    DMGlobalToLocalBegin(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);
    DMGlobalToLocalEnd(fda, ueqn->Ucat, INSERT_VALUES, ueqn->lUcat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode SetWallModels(ueqn_ *ueqn)
{
    // set useful pointers
    mesh_         *mesh  = ueqn->access->mesh;

    flags_        *flags = ueqn->access->flags;
    teqn_         *teqn;

    DMDALocalInfo info = mesh->info;
    PetscInt      xs  = info.xs, xe = info.xs + info.xm;
    PetscInt      ys  = info.ys, ye = info.ys + info.ym;
    PetscInt      zs  = info.zs, ze = info.zs + info.zm;
    PetscInt      mx  = info.mx, my = info.my, mz = info.mz;

    PetscInt      i, j, k;

    PetscInt      mx_wm = xe-xs;
    PetscInt      my_wm = ye-ys;
    PetscInt      mz_wm = ze-zs;

    // set dictionaries
    word          fileNameU = "./boundary/" + mesh->meshName + "/U";
    word          fileNameT = "./boundary/" + mesh->meshName + "/T";

    // read wall functions subdictionaries for wall function bc and set type
    // set type to zero otherwise so that its allocated for the cartesian BC update
    PetscPrintf(mesh->MESH_COMM, "Reading wall model...");

    // velocity wall models
    {
        // i-left boundary wall function
        if (mesh->boundaryU.iLeft == "velocityWallFunction")
        {
            readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(mesh->boundaryU.iLWF));

            // allocate memory and connect pointers
            PetscMalloc(sizeof(wallModel), &(ueqn->iLWM));

            // set to zero the filtered velocity initialization flag
            ueqn->iLWM->uFiltSet = 0;

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->tauWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->tauWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->tauWall.z));

            // initialize velocity average for wall model filtering LLM correction
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->uMeanMag));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->uFilt.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->uFilt.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iLWM->uFilt.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->tauWall.x[k]));
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->tauWall.y[k]));
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->tauWall.z[k]));
                PetscMalloc(sizeof(PetscReal) *mz_wm, &(ueqn->iLWM->uMeanMag[k]));
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->uFilt.x[k]));
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->uFilt.y[k]));
                PetscMalloc(sizeof(PetscReal) *my_wm, &(ueqn->iLWM->uFilt.z[k]));

                for (j = 0; j < my_wm; j++)
                {
                    ueqn->iLWM->tauWall.x[k][j] = 0.0;
                    ueqn->iLWM->tauWall.y[k][j] = 0.0;
                    ueqn->iLWM->tauWall.z[k][j] = 0.0;
                    ueqn->iLWM->uFilt.x  [k][j] = 0.0;
                    ueqn->iLWM->uFilt.y  [k][j] = 0.0;
                    ueqn->iLWM->uFilt.z  [k][j] = 0.0;
                    ueqn->iLWM->uMeanMag [k][j] = 0.0;
                }
            }

            if(mesh->boundaryU.iLWF == -1)
            {
                PetscMalloc(sizeof(Cabot), &(ueqn->iLWM->wmCabot));

                Cabot *wm = ueqn->iLWM->wmCabot;
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));

            }
            else if (mesh->boundaryU.iLWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(ueqn->iLWM->wmShumann));

                Shumann *wm = ueqn->iLWM->wmShumann;
                readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryU.iLWF = 0;
        }

        // i-right boundary wall function
        if (mesh->boundaryU.iRight == "velocityWallFunction")
        {
            readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(mesh->boundaryU.iRWF));

            PetscMalloc(sizeof(wallModel), &(ueqn->iRWM));

            // set to zero the filtered velocity initialization flag
            ueqn->iRWM->uFiltSet = 0;

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->tauWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->tauWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->tauWall.z));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->uMeanMag));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->uFilt.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->uFilt.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->iRWM->uFilt.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->tauWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->tauWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->tauWall.z[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->uMeanMag [k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->uFilt.x[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->uFilt.y[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(ueqn->iRWM->uFilt.z[k]));

                for (j = 0; j < my_wm; j++)
                {
                    ueqn->iRWM->tauWall.x[k][j] =
                    ueqn->iRWM->tauWall.y[k][j] =
                    ueqn->iRWM->tauWall.z[k][j] =
                    ueqn->iRWM->uMeanMag [k][j] =
                    ueqn->iRWM->uFilt.x  [k][j] =
                    ueqn->iRWM->uFilt.y  [k][j] =
                    ueqn->iRWM->uFilt.z  [k][j] =
                    0.0;
                }
            }

            if(mesh->boundaryU.iRWF == -1)
            {
                PetscMalloc(sizeof(Cabot), &(ueqn->iRWM->wmCabot));

                Cabot *wm = ueqn->iRWM->wmCabot;
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
            }
            else if (mesh->boundaryU.iRWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(ueqn->iRWM->wmShumann));

                Shumann *wm = ueqn->iLWM->wmShumann;
                readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryU.iRWF = 0;
        }

        // j-left boundary wall function
        if (mesh->boundaryU.jLeft == "velocityWallFunction")
        {
            readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(mesh->boundaryU.jLWF));

            PetscMalloc(sizeof(wallModel), &(ueqn->jLWM));

            // set to zero the filtered velocity initialization flag
            ueqn->jLWM->uFiltSet = 0;

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->tauWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->tauWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->tauWall.z));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->uMeanMag));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->uFilt.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->uFilt.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jLWM->uFilt.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->tauWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->tauWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->tauWall.z[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->uMeanMag [k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->uFilt.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->uFilt.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jLWM->uFilt.z[k]));

                for (i = 0; i < mx_wm; i++)
                {
                    ueqn->jLWM->tauWall.x[k][i] = 0.0;
                    ueqn->jLWM->tauWall.y[k][i] = 0.0;
                    ueqn->jLWM->tauWall.z[k][i] = 0.0;
                    ueqn->jLWM->uMeanMag [k][i] = 0.0;
                    ueqn->jLWM->uFilt.x  [k][i] = 0.0;
                    ueqn->jLWM->uFilt.y  [k][i] = 0.0;
                    ueqn->jLWM->uFilt.z  [k][i] = 0.0;
                }
            }

            if(mesh->boundaryU.jLWF == -1)
            {
                PetscMalloc(sizeof(Cabot), &(ueqn->jLWM->wmCabot));

                Cabot *wm = ueqn->jLWM->wmCabot;
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));

            }
            else if (mesh->boundaryU.jLWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(ueqn->jLWM->wmShumann));

                Shumann *wm = ueqn->jLWM->wmShumann;
                readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
                //readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "betaM",    &(wm->betaM));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryU.jLWF = 0;
        }

        // j-right boundary wall function
        if (mesh->boundaryU.jRight == "velocityWallFunction")
        {
            readSubDictInt(fileNameU.c_str(), "velocityWallFunction", "type", &(mesh->boundaryU.jRWF));

            PetscMalloc(sizeof(wallModel), &(ueqn->jRWM));

            // set to zero the filtered velocity initialization flag
            ueqn->jRWM->uFiltSet = 0;

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->tauWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->tauWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->tauWall.z));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->uMeanMag));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->uFilt.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->uFilt.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(ueqn->jRWM->uFilt.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->tauWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->tauWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->tauWall.z[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->uMeanMag [k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->uFilt.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->uFilt.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(ueqn->jRWM->uFilt.z[k]));

                for (i = 0; i < mx_wm; i++)
                {
                    ueqn->jRWM->tauWall.x[k][i] = 0.0;
                    ueqn->jRWM->tauWall.y[k][i] = 0.0;
                    ueqn->jRWM->tauWall.z[k][i] = 0.0;
                    ueqn->jRWM->uMeanMag [k][i] = 0.0;
                    ueqn->jRWM->uFilt.x  [k][i] = 0.0;
                    ueqn->jRWM->uFilt.y  [k][i] = 0.0;
                    ueqn->jRWM->uFilt.z  [k][i] = 0.0;
                }
            }

            if(mesh->boundaryU.jRWF == -1)
            {
                PetscMalloc(sizeof(Cabot), &(ueqn->jRWM->wmCabot));

                Cabot *wm = ueqn->jRWM->wmCabot;
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",  &(wm->roughness));
            }
            else if (mesh->boundaryU.jRWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(ueqn->jRWM->wmShumann));

                Shumann *wm = ueqn->jLWM->wmShumann;
                readSubDictWord  (fileNameU.c_str(), "velocityWallFunction", "uStarEval", &(wm->wfEvalType));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kappa",     &(wm->kappa));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "thetaRef",  &(wm->thetaRef));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "kRough",    &(wm->roughness));
                readSubDictDouble(fileNameU.c_str(), "velocityWallFunction", "gammaM",    &(wm->gammaM));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryU.jRWF = 0;
        }
    }

    // temperature wall models
    if(flags->isTeqnActive)
    {
        // get teqn pointer
        teqn = ueqn->access->teqn;

        // i-left boundary wall function
        if (mesh->boundaryT.iLeft == "thetaWallFunction")
        {
            readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(mesh->boundaryT.iLWF));

            // allocate memory
            PetscMalloc(sizeof(wallModel), &(teqn->iLWM));

            // initialize the patch field for Visc term at the wall in the temperature eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iLWM->qWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iLWM->qWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iLWM->qWall.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iLWM->qWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iLWM->qWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iLWM->qWall.z[k]));

                for (j = 0; j < my_wm; j++)
                {
                    teqn->iLWM->qWall.x[k][j] = 0.0;
                    teqn->iLWM->qWall.y[k][j] = 0.0;
                    teqn->iLWM->qWall.z[k][j] = 0.0;
                }
            }

            if (mesh->boundaryT.iLWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->iLWM->wmShumann));

                Shumann *wm = teqn->iLWM->wmShumann;
                readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));
            }
            else if(mesh->boundaryT.iLWF == -2)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->iLWM->wmShumann));

                Shumann *wm = teqn->iLWM->wmShumann;

                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryT.iLWF = 0;
        }

        // i-right boundary wall function
        if (mesh->boundaryT.iRight == "thetaWallFunction")
        {
            readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(mesh->boundaryT.iRWF));

            PetscMalloc(sizeof(wallModel), &(teqn->iRWM));

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iRWM->qWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iRWM->qWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->iRWM->qWall.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iRWM->qWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iRWM->qWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*my_wm, &(teqn->iRWM->qWall.z[k]));

                for (j = 0; j < my_wm; j++)
                {
                    teqn->iRWM->qWall.x[k][j] = 0.0;
                    teqn->iRWM->qWall.y[k][j] = 0.0;
                    teqn->iRWM->qWall.z[k][j] = 0.0;
                }
            }

            if (mesh->boundaryT.iRWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->iRWM->wmShumann));

                Shumann *wm = teqn->iRWM->wmShumann;
                readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));
            }
            else if(mesh->boundaryT.iRWF == -2)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->iRWM->wmShumann));

                Shumann *wm = teqn->iRWM->wmShumann;

                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryT.iRWF = 0;
        }

        // j-left boundary wall function
        if (mesh->boundaryT.jLeft == "thetaWallFunction")
        {
            readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(mesh->boundaryT.jLWF));

            PetscMalloc(sizeof(wallModel), &(teqn->jLWM));

            // initialize the patch field for heat flux term at the wall in the temperature eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jLWM->qWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jLWM->qWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jLWM->qWall.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jLWM->qWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jLWM->qWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jLWM->qWall.z[k]));

                for (i = 0; i < mx_wm; i++)
                {
                    teqn->jLWM->qWall.x[k][i] = 0.0;
                    teqn->jLWM->qWall.y[k][i] = 0.0;
                    teqn->jLWM->qWall.z[k][i] = 0.0;
                }
            }

            if (mesh->boundaryT.jLWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->jLWM->wmShumann));

                Shumann *wm = teqn->jLWM->wmShumann;
                readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));

                // force surface temperature update at first time step
                wm->tLast = -1.0;

                // set surface theta flag to zero
                wm->surfaceThetaSet = 0;
            }
            else if(mesh->boundaryT.jLWF == -2)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->jLWM->wmShumann));
                
                Shumann *wm = teqn->jLWM->wmShumann;

                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
            }
            else if(mesh->boundaryT.jLWF == -4)
            {

                PetscMalloc(sizeof(Shumann), &(teqn->jLWM->wmShumann));

                Shumann *wm = teqn->jLWM->wmShumann;
                readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));

                //read the surface temp and Obhukhov length 
                readSurfaceTempData(wm);
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -1 or -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryT.jLWF = 0;
        }

        // j-right boundary wall function
        if (mesh->boundaryT.jRight == "thetaWallFunction")
        {
            readSubDictInt(fileNameT.c_str(), "thetaWallFunction", "type", &(mesh->boundaryT.jRWF));

            PetscMalloc(sizeof(wallModel), &(teqn->jRWM));

            // initialize the patch field for Visc term at the wall in the momentum eqn.
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jRWM->qWall.x));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jRWM->qWall.y));
            PetscMalloc(sizeof(PetscReal*)*mz_wm, &(teqn->jRWM->qWall.z));

            for (k = 0; k < mz_wm; k++)
            {
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jRWM->qWall.x[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jRWM->qWall.y[k]));
                PetscMalloc(sizeof(PetscReal)*mx_wm, &(teqn->jRWM->qWall.z[k]));

                for (i = 0; i < mx_wm; i++)
                {
                    teqn->jRWM->qWall.x[k][i] = 0.0;
                    teqn->jRWM->qWall.y[k][i] = 0.0;
                    teqn->jRWM->qWall.z[k][i] = 0.0;
                }
            }

            if (mesh->boundaryT.jRWF == -3)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->jRWM->wmShumann));

                Shumann *wm = teqn->jRWM->wmShumann;
                readSubDictWord  (fileNameT.c_str(), "thetaWallFunction", "uStarEval",   &(wm->wfEvalType));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kappa",       &(wm->kappa));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "thetaRef",    &(wm->thetaRef));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "kRough",      &(wm->roughness));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaM",      &(wm->gammaM));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaM",       &(wm->betaM));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "gammaH",      &(wm->gammaH));
                //readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "betaH",       &(wm->betaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "alphaH",      &(wm->alphaH));
                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "heatingRate", &(wm->heatingRate));
            }
            else if(mesh->boundaryT.jRWF == -2)
            {
                PetscMalloc(sizeof(Shumann), &(teqn->jRWM->wmShumann));
                
                Shumann *wm = teqn->jRWM->wmShumann;

                readSubDictDouble(fileNameT.c_str(), "thetaWallFunction", "qWall",       &(wm->qWall));
            }
            else
            {
                char error[512];
                sprintf(error, "invalid wall model chosen. Please use option -3 \n");
                fatalErrorInFunction("SetWallModels", error);
            }
        }
        else
        {
            // zero value on the wall function means no wall function
            mesh->boundaryT.jRWF = 0;
        }
    }

    PetscPrintf(mesh->MESH_COMM, "done\n\n");

    MPI_Barrier(mesh->MESH_COMM);

    return(0);
}

PetscErrorCode readSurfaceTempData(Shumann *wm)
{
    word variableName;
    word fileName;
    word filePath = "inflowDatabase/mesoscaleData";

    PetscInt  dim1, dim2;

    PetscInt  numtV;

    std::vector<word> fileList;
    PetscInt nFiles = 0;

    getFileList(filePath.c_str(), fileList, nFiles);

    if(nFiles == 0)
    {
        char error[512];
        sprintf(error, "no file found in folder %s\n", filePath.c_str());
        fatalErrorInFunction("readSurfaceTempData",  error);
    }

    PetscPrintf(PETSC_COMM_WORLD, "   reading surface temperature and Obhukhov Length data\n");

    // Check if "surfTemp" file exists in fileList
    bool surfTempExists = false;
    for(PetscInt i = 0; i < nFiles; i++)
    {
        if (fileList[i] == "surfTemp")
        {
            surfTempExists = true;
            break;
        }
    }

    if (!surfTempExists)
    {
        char error[512];
        sprintf(error, "Error: 'surfTemp' file not found in folder %s\n", filePath.c_str());
        fatalErrorInFunction("readSurfaceTempData", error);
    }
    
    for(PetscInt i=0; i<nFiles; i++)
    {
        fileName = filePath + "/" + fileList[i];

        std::ifstream indata;
        indata.open(fileName.c_str());

        if(!indata)
        {
            char error[512];
            sprintf(error, "cannot open file %s\n", fileName.c_str());
            fatalErrorInFunction("readSurfaceTempData",  error);
        }
        else
        {
            indata >> dim1 >> dim2;

            if(fileList[i] == "velTime" || fileList[i] == "time")
            {
                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2; 
                }
                else 
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(wm->timeVec));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> wm->timeVec[j];
                }

                wm->numT = arrSize;
            }

            if(fileList[i] == "surfTemp")
            {
                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2; 
                }
                else 
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(wm->surfTemp));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> wm->surfTemp[j];
                }
            }

            if(fileList[i] == "L")
            {
                PetscInt arrSize;

                if(dim1 == 1)
                {
                    arrSize = dim2; 
                }
                else 
                {
                    arrSize = dim1;
                }

                PetscMalloc(sizeof(PetscReal) * arrSize, &(wm->surfL));

                for(PetscInt j=0; j<arrSize; j++)
                {
                    indata >> wm->surfL[j];
                }
            }
                        
            indata.close();
        }
    }

    return(0);
}
