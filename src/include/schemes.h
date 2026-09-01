//! \file  schemes.h
//! \brief Numerical schemes for advection/divergence operators

#ifndef _SCHEMES_
#define _SCHEMES_

// NOTE: This file assumes that base types and inline helper functions (nScale, nSub, nSum, 
// isIBMCell, etc.) have already been defined. It should be included from inline.h after 
// those definitions.

// NUMERICAL SCHEMES
// ============================================================================================================= //

inline PetscReal central(PetscReal f0, PetscReal f1)
{
    return 0.5 * (f0 + f1);
}

//***************************************************************************************************************//

inline Cmpnts centralVec(Cmpnts f0, Cmpnts f1)
{
    Cmpnts result;
    result.x = 0.5 * (f0.x + f1.x);
    result.y = 0.5 * (f0.y + f1.y);
    result.z = 0.5 * (f0.z + f1.z);

    return result;
}

//***************************************************************************************************************//

inline symmTensor centralSymmT(symmTensor f0, symmTensor f1)
{
    symmTensor result;
    result.xx = 0.5 * (f0.xx + f1.xx);
    result.yy = 0.5 * (f0.yy + f1.yy);
    result.zz = 0.5 * (f0.zz + f1.zz);
    result.xy = 0.5 * (f0.xy + f1.xy);
    result.xz = 0.5 * (f0.xz + f1.xz);
    result.yz = 0.5 * (f0.yz + f1.yz);

    return result;
}

//***************************************************************************************************************//

inline PetscReal wCentral(PetscReal f0, PetscReal f1, PetscReal d0, PetscReal d1)
{
    // distance btw upwind and downwind cell centers
    PetscReal d  = 0.5 * (d0 + d1);

    // distance btw upwind cell center and face center
    PetscReal dF = 0.5 * d0;

    return((f1 - f0)/d * dF + f0);
}

//***************************************************************************************************************//

inline Cmpnts wCentralVec(Cmpnts f0, Cmpnts f1, PetscReal d0, PetscReal d1)
{
    // distance btw upwind and downwind cell centers
    PetscReal d  = 0.5 * (d0 + d1);

    // distance btw upwind cell center and face center
    PetscReal dF = 0.5 * d0;

    return
    (
        nSum
        (
            f0,
            nScale(dF/d, nSub(f1, f0))
        )
    );
}

//***************************************************************************************************************//

inline Cmpnts centralVec4thCsi(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt mx, 
        PetscReal ***nvert, PetscReal ***meshTag, Cmpnts ***ucat, PetscReal wavespeed, PetscReal hyperVisc)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    PetscInt im1, i0, ip1, ip2;
    Cmpnts result;
    if(i==0 || i==mx-2)
    {
        if(mesh->i_periodic)
        {
            im1 = mx-3, i0  = i, ip1  = i+1, ip2 = 2;

            result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
            result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
            result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].x - ucat[k][j][im1].x) - 3.0 * (ucat[k][j][ip1].x - ucat[k][j][i0].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].y - ucat[k][j][im1].y) - 3.0 * (ucat[k][j][ip1].y - ucat[k][j][i0].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].z - ucat[k][j][im1].z) - 3.0 * (ucat[k][j][ip1].z - ucat[k][j][i0].z));
        }
        else if(mesh->ii_periodic && i==mx-2)
        {
            im1 = i-1, i0  = i, ip1  = i+1, ip2 = mx+2;
            
            result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
            result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
            result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].x - ucat[k][j][im1].x) - 3.0 * (ucat[k][j][ip1].x - ucat[k][j][i0].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].y - ucat[k][j][im1].y) - 3.0 * (ucat[k][j][ip1].y - ucat[k][j][i0].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].z - ucat[k][j][im1].z) - 3.0 * (ucat[k][j][ip1].z - ucat[k][j][i0].z));
        }
        else if(mesh->ii_periodic && i==0)
        {
            im1 = -3, i0  = i, ip1  = i+1, ip2 = i+2;

            result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
            result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
            result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].x - ucat[k][j][im1].x) - 3.0 * (ucat[k][j][ip1].x - ucat[k][j][i0].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].y - ucat[k][j][im1].y) - 3.0 * (ucat[k][j][ip1].y - ucat[k][j][i0].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].z - ucat[k][j][im1].z) - 3.0 * (ucat[k][j][ip1].z - ucat[k][j][i0].z));
        }
        else
        {
            i0 = i, ip1 = i+1;

            result.x = 0.5 * (ucat[k][j][i0].x + ucat[k][j][ip1].x);
            result.y = 0.5 * (ucat[k][j][i0].y + ucat[k][j][ip1].y);
            result.z = 0.5 * (ucat[k][j][i0].z + ucat[k][j][ip1].z);
        }
    }
    else if ( (i == 1 || i == mx-3) && (!isIBMIFace(k, j, i, i+1, nvert)) && (!isOversetIFace(k, j, i, i+1, nvert)) && (!(mesh->i_periodic || mesh->ii_periodic))) // Quick/Second-order central scheme for j=1 and j=my-3 depending on direction of flow
    {
        im1 = i-1, i0 = i, ip1 = i+1, ip2 = i+2;
        
        if (i == 1)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[k][j][i0].x + ucat[k][j][ip1].x);
                result.y = 0.5 * (ucat[k][j][i0].y + ucat[k][j][ip1].y);
                result.z = 0.5 * (ucat[k][j][i0].z + ucat[k][j][ip1].z);
            }
            else
            {
                result.x = 0.5 * (ucat[k][j][ip1].x + ucat[k][j][i0].x) + (2.0*ucat[k][j][ip1].x - ucat[k][j][i0].x - ucat[k][j][ip2].x) / 8.0;
                result.y = 0.5 * (ucat[k][j][ip1].y + ucat[k][j][i0].y) + (2.0*ucat[k][j][ip1].y - ucat[k][j][i0].y - ucat[k][j][ip2].y) / 8.0;
                result.z = 0.5 * (ucat[k][j][ip1].z + ucat[k][j][i0].z) + (2.0*ucat[k][j][ip1].z - ucat[k][j][i0].z - ucat[k][j][ip2].z) / 8.0;
            }
        }

        if (i == mx-3)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[k][j][ip1].x + ucat[k][j][i0].x) + (2.0*ucat[k][j][i0].x - ucat[k][j][ip1].x - ucat[k][j][im1].x) / 8.0;
                result.y = 0.5 * (ucat[k][j][ip1].y + ucat[k][j][i0].y) + (2.0*ucat[k][j][i0].y - ucat[k][j][ip1].y - ucat[k][j][im1].y) / 8.0;
                result.z = 0.5 * (ucat[k][j][ip1].z + ucat[k][j][i0].z) + (2.0*ucat[k][j][i0].z - ucat[k][j][ip1].z - ucat[k][j][im1].z) / 8.0;
            }
            else
            {
                result.x = 0.5 * (ucat[k][j][i0].x + ucat[k][j][ip1].x);
                result.y = 0.5 * (ucat[k][j][i0].y + ucat[k][j][ip1].y);
                result.z = 0.5 * (ucat[k][j][i0].z + ucat[k][j][ip1].z);
            }
        }
    }
    else if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        im1 = i-1, i0 = i, ip1 = i+1, ip2 = i+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[k][j][i0].x + ucat[k][j][ip1].x);
            result.y = 0.5 * (ucat[k][j][i0].y + ucat[k][j][ip1].y);
            result.z = 0.5 * (ucat[k][j][i0].z + ucat[k][j][ip1].z);
        }
        else
        {
            result.x = 0.5 * (ucat[k][j][ip1].x + ucat[k][j][i0].x) + (2.0*ucat[k][j][ip1].x - ucat[k][j][i0].x - ucat[k][j][ip2].x) / 8.0;
            result.y = 0.5 * (ucat[k][j][ip1].y + ucat[k][j][i0].y) + (2.0*ucat[k][j][ip1].y - ucat[k][j][i0].y - ucat[k][j][ip2].y) / 8.0;
            result.z = 0.5 * (ucat[k][j][ip1].z + ucat[k][j][i0].z) + (2.0*ucat[k][j][ip1].z - ucat[k][j][i0].z - ucat[k][j][ip2].z) / 8.0;
        }
    }
    else if(isIBMCell(k, j, i+1, nvert) || isOversetCell(k, j, i+1, meshTag))
    {
        im1 = i-1, i0 = i, ip1 = i+1, ip2 = i+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[k][j][ip1].x + ucat[k][j][i0].x) + (2.0*ucat[k][j][i0].x - ucat[k][j][ip1].x - ucat[k][j][im1].x) / 8.0;
            result.y = 0.5 * (ucat[k][j][ip1].y + ucat[k][j][i0].y) + (2.0*ucat[k][j][i0].y - ucat[k][j][ip1].y - ucat[k][j][im1].y) / 8.0;
            result.z = 0.5 * (ucat[k][j][ip1].z + ucat[k][j][i0].z) + (2.0*ucat[k][j][i0].z - ucat[k][j][ip1].z - ucat[k][j][im1].z) / 8.0;
        }
        else
        {
            result.x = 0.5 * (ucat[k][j][i0].x + ucat[k][j][ip1].x);
            result.y = 0.5 * (ucat[k][j][i0].y + ucat[k][j][ip1].y);
            result.z = 0.5 * (ucat[k][j][i0].z + ucat[k][j][ip1].z);
        }
    }
    else
    {
        im1 = i-1, i0  = i, ip1  = i+1, ip2 = i+2;
            
        result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
        result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
        result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;

        result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].x - ucat[k][j][im1].x) - 3.0 * (ucat[k][j][ip1].x - ucat[k][j][i0].x));
        result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].y - ucat[k][j][im1].y) - 3.0 * (ucat[k][j][ip1].y - ucat[k][j][i0].y));
        result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][j][ip2].z - ucat[k][j][im1].z) - 3.0 * (ucat[k][j][ip1].z - ucat[k][j][i0].z));
    }

    return result;
}

//***************************************************************************************************************//

inline Cmpnts centralVec4thEta(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt my, 
                               PetscReal ***nvert, PetscReal ***meshTag, Cmpnts ***ucat, PetscReal wavespeed, PetscReal hyperVisc)
{
    PetscInt jm1, j0, jp1, jp2;
    Cmpnts result;

    if (j == 0 || j == my-2)
    {
        if (mesh->j_periodic) // Assuming j_periodic for y-direction periodicity
        {
            jm1 = my-3, j0 = j, jp1 = j+1, jp2 = 2;

            result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
            result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
            result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].x - ucat[k][jm1][i].x) - 3.0 * (ucat[k][jp1][i].x - ucat[k][j0][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].y - ucat[k][jm1][i].y) - 3.0 * (ucat[k][jp1][i].y - ucat[k][j0][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].z - ucat[k][jm1][i].z) - 3.0 * (ucat[k][jp1][i].z - ucat[k][j0][i].z));
        }
        else if (mesh->jj_periodic && j == my-2) // Assuming jj_periodic for y-direction
        {
            jm1 = j-1, j0 = j, jp1 = j+1, jp2 = my+2;

            result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
            result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
            result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].x - ucat[k][jm1][i].x) - 3.0 * (ucat[k][jp1][i].x - ucat[k][j0][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].y - ucat[k][jm1][i].y) - 3.0 * (ucat[k][jp1][i].y - ucat[k][j0][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].z - ucat[k][jm1][i].z) - 3.0 * (ucat[k][jp1][i].z - ucat[k][j0][i].z));
        }
        else if (mesh->jj_periodic && j == 0)
        {
            jm1 = -3, j0 = j, jp1 = j+1, jp2 = j+2;

            result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
            result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
            result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].x - ucat[k][jm1][i].x) - 3.0 * (ucat[k][jp1][i].x - ucat[k][j0][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].y - ucat[k][jm1][i].y) - 3.0 * (ucat[k][jp1][i].y - ucat[k][j0][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].z - ucat[k][jm1][i].z) - 3.0 * (ucat[k][jp1][i].z - ucat[k][j0][i].z));
        }
        else
        {
            // use second order central at boundary
            j0 = j, jp1 = j+1;

            result.x = 0.5 * (ucat[k][jp1][i].x + ucat[k][j0][i].x);
            result.y = 0.5 * (ucat[k][jp1][i].y + ucat[k][j0][i].y);
            result.z = 0.5 * (ucat[k][jp1][i].z + ucat[k][j0][i].z);
        }
    }
    else if ( (j == 1 || j == my-3) && (!isIBMJFace(k, j, i, j+1, nvert)) && (!isOversetJFace(k, j, i, j+1, nvert)) && (!(mesh->j_periodic || mesh->jj_periodic))) // Quick/Second-order central scheme for j=1 and j=my-3 depending on direction of flow
    {
        jm1 = j-1, j0 = j, jp1 = j+1, jp2 = j+2;
        
        if (j == 1)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[k][j0][i].x + ucat[k][jp1][i].x);
                result.y = 0.5 * (ucat[k][j0][i].y + ucat[k][jp1][i].y);
                result.z = 0.5 * (ucat[k][j0][i].z + ucat[k][jp1][i].z);
            }
            else
            {
                result.x = 0.5 * (ucat[k][jp1][i].x + ucat[k][j0][i].x) + (2.0*ucat[k][jp1][i].x - ucat[k][j0][i].x - ucat[k][jp2][i].x) / 8.0;
                result.y = 0.5 * (ucat[k][jp1][i].y + ucat[k][j0][i].y) + (2.0*ucat[k][jp1][i].y - ucat[k][j0][i].y - ucat[k][jp2][i].y) / 8.0;
                result.z = 0.5 * (ucat[k][jp1][i].z + ucat[k][j0][i].z) + (2.0*ucat[k][jp1][i].z - ucat[k][j0][i].z - ucat[k][jp2][i].z) / 8.0;
            }
        }

        if (j == my-3)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[k][jp1][i].x + ucat[k][j0][i].x) + (2.0*ucat[k][j0][i].x - ucat[k][jp1][i].x - ucat[k][jm1][i].x) / 8.0;
                result.y = 0.5 * (ucat[k][jp1][i].y + ucat[k][j0][i].y) + (2.0*ucat[k][j0][i].y - ucat[k][jp1][i].y - ucat[k][jm1][i].y) / 8.0;
                result.z = 0.5 * (ucat[k][jp1][i].z + ucat[k][j0][i].z) + (2.0*ucat[k][j0][i].z - ucat[k][jp1][i].z - ucat[k][jm1][i].z) / 8.0;
            }
            else
            {
                result.x = 0.5 * (ucat[k][j0][i].x + ucat[k][jp1][i].x);
                result.y = 0.5 * (ucat[k][j0][i].y + ucat[k][jp1][i].y);
                result.z = 0.5 * (ucat[k][j0][i].z + ucat[k][jp1][i].z);
            }
        }
    }
    else if (isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        jm1 = j-1, j0 = j, jp1 = j+1, jp2 = j+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[k][j0][i].x + ucat[k][jp1][i].x);
            result.y = 0.5 * (ucat[k][j0][i].y + ucat[k][jp1][i].y);
            result.z = 0.5 * (ucat[k][j0][i].z + ucat[k][jp1][i].z);
        }
        else
        {
            result.x = 0.5 * (ucat[k][jp1][i].x + ucat[k][j0][i].x) + (2.0*ucat[k][jp1][i].x - ucat[k][j0][i].x - ucat[k][jp2][i].x) / 8.0;
            result.y = 0.5 * (ucat[k][jp1][i].y + ucat[k][j0][i].y) + (2.0*ucat[k][jp1][i].y - ucat[k][j0][i].y - ucat[k][jp2][i].y) / 8.0;
            result.z = 0.5 * (ucat[k][jp1][i].z + ucat[k][j0][i].z) + (2.0*ucat[k][jp1][i].z - ucat[k][j0][i].z - ucat[k][jp2][i].z) / 8.0;
        }
    }
    else if (isIBMCell(k, j+1, i, nvert) || isOversetCell(k, j+1, i, meshTag))
    {
        jm1 = j-1, j0 = j, jp1 = j+1, jp2 = j+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[k][jp1][i].x + ucat[k][j0][i].x) + (2.0*ucat[k][j0][i].x - ucat[k][jp1][i].x - ucat[k][jm1][i].x) / 8.0;
            result.y = 0.5 * (ucat[k][jp1][i].y + ucat[k][j0][i].y) + (2.0*ucat[k][j0][i].y - ucat[k][jp1][i].y - ucat[k][jm1][i].y) / 8.0;
            result.z = 0.5 * (ucat[k][jp1][i].z + ucat[k][j0][i].z) + (2.0*ucat[k][j0][i].z - ucat[k][jp1][i].z - ucat[k][jm1][i].z) / 8.0;
        }
        else
        {
            result.x = 0.5 * (ucat[k][j0][i].x + ucat[k][jp1][i].x);
            result.y = 0.5 * (ucat[k][j0][i].y + ucat[k][jp1][i].y);
            result.z = 0.5 * (ucat[k][j0][i].z + ucat[k][jp1][i].z);
        }
    }
    else
    {
        jm1 = j-1, j0 = j, jp1 = j+1, jp2 = j+2;

        result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
        result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
        result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;

        result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].x - ucat[k][jm1][i].x) - 3.0 * (ucat[k][jp1][i].x - ucat[k][j0][i].x));
        result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].y - ucat[k][jm1][i].y) - 3.0 * (ucat[k][jp1][i].y - ucat[k][j0][i].y));
        result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[k][jp2][i].z - ucat[k][jm1][i].z) - 3.0 * (ucat[k][jp1][i].z - ucat[k][j0][i].z));
    }

    return result;
}

//***************************************************************************************************************//

inline Cmpnts centralVec4thZet(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt mz, 
                               PetscReal ***nvert, PetscReal ***meshTag, Cmpnts ***ucat, PetscReal wavespeed, PetscReal hyperVisc)
{
    PetscInt km1, k0, kp1, kp2;
    Cmpnts result;

    if (k == 0 || k == mz-2)
    {
        if (mesh->k_periodic) // Assuming k_periodic for z-direction periodicity
        {
            km1 = mz-3, k0 = k, kp1 = k+1, kp2 = 2;

            result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
            result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
            result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].x - ucat[km1][j][i].x) - 3.0 * (ucat[kp1][j][i].x - ucat[k0][j][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].y - ucat[km1][j][i].y) - 3.0 * (ucat[kp1][j][i].y - ucat[k0][j][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].z - ucat[km1][j][i].z) - 3.0 * (ucat[kp1][j][i].z - ucat[k0][j][i].z));
        }
        else if (mesh->kk_periodic && k == mz-2) // Assuming kk_periodic for z-direction
        {
            km1 = k-1, k0 = k, kp1 = k+1, kp2 = mz+2;

            result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
            result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
            result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].x - ucat[km1][j][i].x) - 3.0 * (ucat[kp1][j][i].x - ucat[k0][j][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].y - ucat[km1][j][i].y) - 3.0 * (ucat[kp1][j][i].y - ucat[k0][j][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].z - ucat[km1][j][i].z) - 3.0 * (ucat[kp1][j][i].z - ucat[k0][j][i].z));
        }
        else if (mesh->kk_periodic && k == 0)
        {
            km1 = -3, k0 = k, kp1 = k+1, kp2 = k+2;

            result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
            result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
            result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;

            result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].x - ucat[km1][j][i].x) - 3.0 * (ucat[kp1][j][i].x - ucat[k0][j][i].x));
            result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].y - ucat[km1][j][i].y) - 3.0 * (ucat[kp1][j][i].y - ucat[k0][j][i].y));
            result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].z - ucat[km1][j][i].z) - 3.0 * (ucat[kp1][j][i].z - ucat[k0][j][i].z));
        }
        else
        {
            k0 = k, kp1 = k+1;

            result.x = 0.5 * (ucat[k0][j][i].x + ucat[kp1][j][i].x);
            result.y = 0.5 * (ucat[k0][j][i].y + ucat[kp1][j][i].y);
            result.z = 0.5 * (ucat[k0][j][i].z + ucat[kp1][j][i].z);
        }
    }
    else if ( (k == 1 || k == mz-3) && (!isIBMKFace(k, j, i, k+1, nvert)) && (!isOversetKFace(k, j, i, k+1, nvert)) && (!(mesh->k_periodic || mesh->kk_periodic))) // Quick/Second-order central scheme for k=1 and k=mz-3 depending on direction of flow
    {
        km1 = k-1, k0 = k, kp1 = k+1, kp2 = k+2;
        
        if (k == 1)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[k0][j][i].x + ucat[kp1][j][i].x);
                result.y = 0.5 * (ucat[k0][j][i].y + ucat[kp1][j][i].y);
                result.z = 0.5 * (ucat[k0][j][i].z + ucat[kp1][j][i].z);
            }
            else
            {
                result.x = 0.5 * (ucat[kp1][j][i].x + ucat[k0][j][i].x) + (2.0*ucat[kp1][j][i].x - ucat[k0][j][i].x - ucat[kp2][j][i].x) / 8.0;
                result.y = 0.5 * (ucat[kp1][j][i].y + ucat[k0][j][i].y) + (2.0*ucat[kp1][j][i].y - ucat[k0][j][i].y - ucat[kp2][j][i].y) / 8.0;
                result.z = 0.5 * (ucat[kp1][j][i].z + ucat[k0][j][i].z) + (2.0*ucat[kp1][j][i].z - ucat[k0][j][i].z - ucat[kp2][j][i].z) / 8.0;
            }
        }

        if (k == mz-3)
        {
            if(wavespeed>0)
            {
                result.x = 0.5 * (ucat[kp1][j][i].x + ucat[k0][j][i].x) + (2.0*ucat[k0][j][i].x - ucat[kp1][j][i].x - ucat[km1][j][i].x) / 8.0;
                result.y = 0.5 * (ucat[kp1][j][i].y + ucat[k0][j][i].y) + (2.0*ucat[k0][j][i].y - ucat[kp1][j][i].y - ucat[km1][j][i].y) / 8.0;
                result.z = 0.5 * (ucat[kp1][j][i].z + ucat[k0][j][i].z) + (2.0*ucat[k0][j][i].z - ucat[kp1][j][i].z - ucat[km1][j][i].z) / 8.0;
            }
            else
            {
                result.x = 0.5 * (ucat[k0][j][i].x + ucat[kp1][j][i].x);
                result.y = 0.5 * (ucat[k0][j][i].y + ucat[kp1][j][i].y);
                result.z = 0.5 * (ucat[k0][j][i].z + ucat[kp1][j][i].z);
            }
        }
    }
    else if (isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        km1 = k-1, k0 = k, kp1 = k+1, kp2 = k+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[k0][j][i].x + ucat[kp1][j][i].x);
            result.y = 0.5 * (ucat[k0][j][i].y + ucat[kp1][j][i].y);
            result.z = 0.5 * (ucat[k0][j][i].z + ucat[kp1][j][i].z);
        }
        else
        {
            result.x = 0.5 * (ucat[kp1][j][i].x + ucat[k0][j][i].x) + (2.0*ucat[kp1][j][i].x - ucat[k0][j][i].x - ucat[kp2][j][i].x) / 8.0;
            result.y = 0.5 * (ucat[kp1][j][i].y + ucat[k0][j][i].y) + (2.0*ucat[kp1][j][i].y - ucat[k0][j][i].y - ucat[kp2][j][i].y) / 8.0;
            result.z = 0.5 * (ucat[kp1][j][i].z + ucat[k0][j][i].z) + (2.0*ucat[kp1][j][i].z - ucat[k0][j][i].z - ucat[kp2][j][i].z) / 8.0;
        }
    }
    else if (isIBMCell(k+1, j, i, nvert) || isOversetCell(k+1, j, i, meshTag))
    {
        km1 = k-1, k0 = k, kp1 = k+1, kp2 = k+2;

        if(wavespeed>0)
        {
            result.x = 0.5 * (ucat[kp1][j][i].x + ucat[k0][j][i].x) + (2.0*ucat[k0][j][i].x - ucat[kp1][j][i].x - ucat[km1][j][i].x) / 8.0;
            result.y = 0.5 * (ucat[kp1][j][i].y + ucat[k0][j][i].y) + (2.0*ucat[k0][j][i].y - ucat[kp1][j][i].y - ucat[km1][j][i].y) / 8.0;
            result.z = 0.5 * (ucat[kp1][j][i].z + ucat[k0][j][i].z) + (2.0*ucat[k0][j][i].z - ucat[kp1][j][i].z - ucat[km1][j][i].z) / 8.0;
        }
        else
        {
            result.x = 0.5 * (ucat[k0][j][i].x + ucat[kp1][j][i].x);
            result.y = 0.5 * (ucat[k0][j][i].y + ucat[kp1][j][i].y);
            result.z = 0.5 * (ucat[k0][j][i].z + ucat[kp1][j][i].z);
        }
    }
    else
    {
        km1 = k-1, k0 = k, kp1 = k+1, kp2 = k+2;

        result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
        result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
        result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;

        result.x += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].x - ucat[km1][j][i].x) - 3.0 * (ucat[kp1][j][i].x - ucat[k0][j][i].x));
        result.y += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].y - ucat[km1][j][i].y) - 3.0 * (ucat[kp1][j][i].y - ucat[k0][j][i].y));
        result.z += PetscSign(wavespeed) * ((1.0-hyperVisc)/12.0)*((ucat[kp2][j][i].z - ucat[km1][j][i].z) - 3.0 * (ucat[kp1][j][i].z - ucat[k0][j][i].z));
    }

    return result;
}

//***************************************************************************************************************//
inline PetscReal central4(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3)
{
    return 0.0625 * (- f0 + 9. * f1 + 9.* f2 - f3);
}

//***************************************************************************************************************//

inline PetscReal tvd12(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3, PetscReal wavespeed)
{
    // 2nd-order TVD face value with Van Leer limiter
    PetscReal n, d, r, lim;
    PetscReal fU, fD;
    const PetscReal eps = 1e-5;

    if(wavespeed > 0)
    {
        n = f1 - f0;
        d = f2 - f1;
        fU = f1;
        fD = f2;
    }
    else
    {
        n = f2 - f3;
        d = f1 - f2;
        fU = f2;
        fD = f1;
    }

    if(fabs(d) < eps)
    {
        return fU;
    }

    r = n / d;

    // Van Leer limiter: phi(r) = (r + |r|)/(1 + |r|)
    lim = (r + fabs(r)) / (1.0 + fabs(r));

    return fU + 0.5 * lim * (fD - fU);
}

//***************************************************************************************************************//

inline PetscReal tvd34(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3, PetscReal wavespeed)
{
    // 3rd-to-4th order TVD-like face value with Van Leer limiter
    // low-order branch: QUICK (3rd-order upwind-biased)
    // high-order branch: 4th-order central interpolation
    PetscReal n, d, r, lim;
    PetscReal fU;
    const PetscReal eps = 1e-5;

    if(wavespeed > 0)
    {
        n = f1 - f0;
        d = f2 - f1;
        fU = f1;
    }
    else
    {
        n = f2 - f3;
        d = f1 - f2;
        fU = f2;
    }

    if(fabs(d) < eps)
    {
        return fU;
    }

    r = n / d;

    // Van Leer limiter: phi(r) = (r + |r|)/(1 + |r|)
    lim = (r + fabs(r)) / (1.0 + fabs(r));

    // QUICK (3rd-order upwind-biased) face value
    PetscReal q3;
    if(wavespeed > 0)
    {
        q3 = (-f0 + 6.0 * f1 + 3.0 * f2) / 8.0;
    }
    else
    {
        q3 = (-f3 + 6.0 * f2 + 3.0 * f1) / 8.0;
    }

    // 4th-order central face value
    PetscReal c4 = central4(f0, f1, f2, f3);

    // lim=0 uses QUICK (3rd), lim=1 uses central4 (4th)
    return q3 + lim * (c4 - q3);
}

//***************************************************************************************************************//

inline PetscReal centralUpwind(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3, PetscReal wavespeed)
{
    PetscReal r, n, d;
    PetscReal fUU, fU, fD;

    if(wavespeed>0)
    {
        n = f1 - f0; d = PetscMax(f2 - f1, 1e-5); r = n/d;
        fUU = f0; fU  = f1; fD  = f2;
	}
	else
    {
        n = f2 - f3; d = PetscMax(f1 - f2, 1e-5); r = n/d;
        fUU = f3; fU  = f2; fD  = f1;
	}

    // Van Leer limiter function
    PetscReal lim = 0.5 * (r + fabs(r)) / (1.0 + fabs(r));

    // central scheme
    PetscReal C = 0.5 * (fU + fD);

    // quickDiv correction to the central scheme
    PetscReal corr = (2.0*fU - fD - fUU) / 8.0;

	return
    (
       C + (1.0-lim)*corr
    );
}

//***************************************************************************************************************//

inline Cmpnts centralUpwindVec(Cmpnts f0, Cmpnts f1, Cmpnts f2, Cmpnts f3, PetscReal wavespeed)
{
    Cmpnts r, n, d;
    Cmpnts fUU, fU, fD;

    if(wavespeed>0)
    {
        n = nSub(f1,f0); d = nMax(nSub(f2,f1), 1e-5); r = nDivElWise(n,d);
        fUU = nSet(f0); fU  = nSet(f1); fD  = nSet(f2);
	}
	else
    {
        n = nSub(f2,f3); d = nMax(nSub(f1,f2), 1e-5); r = nDivElWise(n,d);
        fUU = nSet(f3); fU  = nSet(f2); fD  = nSet(f1);
	}

    // Van Leer limiter function
    Cmpnts lim = nScale(0.5, nDivElWise(nSum(r, nAbs(r)), nSum(nSetFromComponents(1.0, 1.0, 1.0), nAbs(r))));

    // central scheme
    Cmpnts C = nScale(0.5, nSum(fU, fD));

    // quickDiv correction to the central scheme
    Cmpnts corr = nScale(1./8., nSub(nScale(2.0, fU), nSum(fD, fUU)));

	return
    (
       nSum(C, nMultElWise(nSub(nSetFromComponents(1.0, 1.0, 1.0),lim),corr))
    );
}


//***************************************************************************************************************//

inline PetscReal centralUpwind
(
    PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3,
    PetscReal wavespeed, PetscReal limiter
)
{
    PetscReal fUU, fU, fD;

    if(wavespeed>0)
    {
        fUU = f0; fU  = f1; fD  = f2;
	}
	else
    {
        fUU = f3; fU  = f2; fD  = f1;
	}

    // central scheme
    PetscReal C = 0.5 * (fU + fD);

    // quickDiv correction to the central scheme
    PetscReal corr = (2.0*fU - fD - fUU) / 8.0;

	return
    (
       C + (1.0-limiter)*corr
    );
}

//***************************************************************************************************************//

inline PetscReal upwind
(
    PetscReal f1, PetscReal f2,
    PetscReal wavespeed
)
{
    PetscReal fU, fD;

    if(wavespeed>0)
    {
        fU  = f1;
	}
	else
    {
        fU  = f2;
	}

	return
    (
       fU
    );
}

//***************************************************************************************************************//

inline Cmpnts upwindVec
(
    Cmpnts f1, Cmpnts f2,
    PetscReal wavespeed
)
{
    Cmpnts fU, fD;

    if(wavespeed>0)
    {
        fU  = f1;
	}
	else
    {
        fU  = f2;
	}

	return
    (
       fU
    );
}

//***************************************************************************************************************//

inline Cmpnts centralUpwindVec
(
    Cmpnts f0, Cmpnts f1, Cmpnts f2, Cmpnts f3,
    PetscReal wavespeed, PetscReal limiter
)
{
    Cmpnts fUU, fU, fD;

    if(wavespeed>0)
    {
        fUU = nSet(f0); fU  = nSet(f1); fD  = nSet(f2);
	}
	else
    {
        fUU = nSet(f3); fU  = nSet(f2); fD  = nSet(f1);
	}

    // central scheme
    Cmpnts C = nScale(0.5, nSum(fU,fD));

    // quickDiv correction to the central scheme
    Cmpnts corr = nScale(1./8., nSub(nScale(2.0,fU), nSum(fD,fUU)));

	return
    (
       nSum(C, nScale(1.0-limiter,corr))
    );
}

//***************************************************************************************************************//

inline PetscReal wCentralUpwind
(
    PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3,
    PetscReal d0, PetscReal d1, PetscReal d2, PetscReal d3,
    PetscReal wavespeed, PetscReal limiter
)
{
    PetscReal fUU, fU, fD;
    PetscReal dUU, dU, dD;

    if(d0 != d0 || d1 != d1 || d2 != d2 || d3 != d3)
    {
        char error[512];
        sprintf(error, "nan-distances provided\n");
        fatalErrorInFunction("wCentralUpwind",  error);
    }

    if(wavespeed > 0)
    {
        fUU = f0; fU  = f1; fD  = f2;
        dUU = d0; dU  = d1; dD  = d2;
	}
	else
    {
        fUU = f3; fU  = f2; fD  = f1;
        dUU = d3; dU  = d2; dD  = d1;
	}

    // compute cell to cell and cell to face distances
    PetscReal da  = 0.5 * (dUU + dU);
    PetscReal db  = 0.5 * (dU  + dD);
    PetscReal dc  = 0.5 * (dUU + 2.0*dU + dD);
    PetscReal dfQ = 0.5 * (2.0*dU + dUU);
    PetscReal dfC = 0.5 * dU;

    // central scheme
    PetscReal C = (fD - fU) / db * dfC + fU;

    // quickDiv correction a, b, c terms
    PetscReal a = (da*fD - dc*fU + db*fUU);
    PetscReal b = ((dc*dc+dc*da)*fU - (da*da+dc*da)*fD - (dc*dc-da*da)*fUU);
    PetscReal c = ((fD - fU) / db * da + fUU - fU ) * (da*db*dc);

    // quadratic quickDiv correction to the central scheme
    PetscReal corr = (a*dfQ*dfQ + b*dfQ + c ) / (da*db*dc);

	return
    (
       C + (1.0 - limiter) * corr
    );
}

//***************************************************************************************************************//

inline Cmpnts wCentralUpwindVec
(
    Cmpnts f0, Cmpnts f1, Cmpnts f2, Cmpnts f3,
    PetscReal d0, PetscReal d1, PetscReal d2, PetscReal d3,
    PetscReal wavespeed, PetscReal limiter
)
{
    Cmpnts    fUU, fU, fD;
    PetscReal dUU, dU, dD;

    if(d0 != d0 || d1 != d1 || d2 != d2 || d3 != d3)
    {
        char error[512];
        sprintf(error, "nan-distances provided\n");
        fatalErrorInFunction("wCentralUpwind",  error);
    }

    if(wavespeed > 0)
    {
        fUU = nSet(f0); fU  = nSet(f1); fD  = nSet(f2);
        dUU = d0; dU  = d1; dD  = d2;
	}
	else
    {
        fUU = nSet(f3); fU  = nSet(f2); fD  = nSet(f1);
        dUU = d3; dU  = d2; dD  = d1;
	}

    // compute cell to cell and cell to face distances
    PetscReal da  = 0.5 * (dUU + dU);
    PetscReal db  = 0.5 * (dU  + dD);
    PetscReal dc  = 0.5 * (dUU + 2.0*dU + dD);
    PetscReal dfQ = 0.5 * (2.0*dU + dUU);
    PetscReal dfC = 0.5 * dU;

    // central scheme
    Cmpnts C = nSum(nScale(dfC/ db, nSub(fD, fU)), fU);

    // quickDiv correction a, b, c terms
    Cmpnts a = nSum(nSub(nScale(da,fD), nScale(dc,fU)), nScale(db,fUU));
    Cmpnts b = nSub(nScale(dc*dc+dc*da,fU), nSum(nScale(da*da+dc*da,fD), nScale(dc*dc-da*da,fUU)));
    Cmpnts c = nScale(da*db*dc, nSub(nSum(nScale(da/db, nSub(fD,fU)), fUU), fU));

    // quadratic quickDiv correction to the central scheme
    Cmpnts corr = nScale(1.0/(da*db*dc), nSum(nSum(nScale(dfQ*dfQ, a), nScale(dfQ, b)), c));

	return
    (
       nSum(C, nScale(1.0 - limiter, corr))
    );
}

//***************************************************************************************************************//

inline PetscReal quadraticUpwind
(
    PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3,
    PetscReal wavespeed
)
{
	PetscReal fUU, fU, fD;

	if(wavespeed>0)
    {
		fUU = f0; fU  = f1; fD  = f2;
	}
	else
    {
		fUU = f3; fU  = f2; fD  = f1;
	}

    return 0.5 * (fU + fD) + (2.0*fU - fD - fUU) / 8.0;
}

//***************************************************************************************************************//

inline Cmpnts quadraticUpwindVec
(
    Cmpnts f0, Cmpnts f1, Cmpnts f2, Cmpnts f3,
    PetscReal wavespeed
)
{
	Cmpnts fUU, fU, fD;

	if(wavespeed>0)
    {
		fUU = nSet(f0); fU  = nSet(f1); fD  = nSet(f2);
	}
	else
    {
		fUU = nSet(f3); fU  = nSet(f2); fD  = nSet(f1);
	}

    return nSum(nScale(0.5, nSum(fU, fD)), nScale(1./8., nSub(nScale(2.0,fU), nSum(fD, fUU))));
}

//***************************************************************************************************************//

inline PetscReal wQuadraticUpwind
(
    PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3,
    PetscReal d0, PetscReal d1, PetscReal d2, PetscReal d3,
    PetscReal wavespeed
)
{
    PetscReal fUU, fU, fD;
    PetscReal dUU, dU, dD;

    if(wavespeed > 0)
    {
        fUU = f0; fU  = f1; fD  = f2;
        dUU = d0; dU  = d1; dD  = d2;
	}
	else
    {
        fUU = f3; fU  = f2; fD  = f1;
        dUU = d3; dU  = d2; dD  = d1;
	}

    // compute cell to cell and cell to face distances
    PetscReal da  = 0.5 * (dUU + dU);
    PetscReal db  = 0.5 * (dU  + dD);
    PetscReal dc  = 0.5 * (dUU + 2.0*dU + dD);
    PetscReal dfQ = 0.5 * (2.0*dU + dUU);

    // quickDiv scheme a, b terms (c = fUU)
    PetscReal a = (da*fD - dc*fU + db*fUU) / (da*db*dc);
    PetscReal b = (dc*dc*fU - da*da*fD - (dc*dc-da*da)*fUU) / (da*db*dc);

    return(a*dfQ*dfQ + b*dfQ + fUU);
}

//***************************************************************************************************************//

inline PetscReal weno3(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3, PetscReal wavespeed)
{
	PetscReal fL, fC, fR;

	if(wavespeed>0)
    {
		fL = f0; fC = f1; fR = f2;
	}
	else
    {
		fL = f3; fC = f2; fR = f1;
	}

	PetscReal d0=2./3., d1=1./3.;

	const PetscReal eps=1.e-6;

	PetscReal beta0 = pow( fC - fR,  2. );
	PetscReal beta1 = pow( fL - fC,  2. );

	PetscReal alpha0 = d0 / pow( eps + beta0, 2. );
	PetscReal alpha1 = d1 / pow( eps + beta1, 2. );

	PetscReal sumalpha = alpha0 + alpha1;

	PetscReal w0 = alpha0 / sumalpha;
	PetscReal w1 = alpha1 / sumalpha;

	PetscReal u0 = ( fC*0.5  + fR*0.5 );
	PetscReal u1 = ( -fL*0.5 + fC*1.5 );

	return w0*u0 + w1*u1;
}

//***************************************************************************************************************//

inline Cmpnts weno3Vec(Cmpnts f0, Cmpnts f1, Cmpnts f2, Cmpnts f3, PetscReal wavespeed)
{
	Cmpnts fL, fC, fR;

	if(wavespeed>0)
    {
		fL = nSet(f0); fC = nSet(f1); fR = nSet(f2);
	}
	else
    {
		fL = nSet(f3); fC = nSet(f2); fR = nSet(f1);
	}

	PetscReal d0=2./3., d1=1./3.;

	const Cmpnts eps = nSetFromComponents(1e-6, 1e-6, 1e-6);

	Cmpnts beta0 = nPow(nSub(fC,fR),  2.);
	Cmpnts beta1 = nPow(nSub(fL,fC),  2.);

	Cmpnts alpha0 = nScale(d0, nInv(nPow( nSum(eps,beta0), 2. )));
	Cmpnts alpha1 = nScale(d1, nInv(nPow( nSum(eps,beta1), 2. )));

	Cmpnts sumalpha = nSum(alpha0, alpha1);

	Cmpnts w0 = nMultElWise(alpha0, nInv(sumalpha));
	Cmpnts w1 = nMultElWise(alpha1, nInv(sumalpha));

	Cmpnts u0 = nSum(nScale(0.5,fC),  nScale(0.5,fR));
	Cmpnts u1 = nSum(nScale(-0.5,fL), nScale(1.5,fC));

	return nSum(nMultElWise(w0,u0), nMultElWise(w1,u1));
}

//***************************************************************************************************************//

inline PetscReal weno5(PetscReal f0, PetscReal f1, PetscReal f2, PetscReal f3, PetscReal f4, PetscReal f5, PetscReal wavespeed)
{
	PetscReal A, B, C, D, E;

	if(wavespeed>0)
    {
		A = f0; B = f1; C = f2; D = f3; E = f4;
	}
	else
    {
		A = f5; B = f4; C = f3; D = f2; E = f1;
	}

	PetscReal eps = 1.e-6;
	PetscReal d0, d1, d2;

	d0=0.3, d1=0.6, d2=0.1;

	PetscReal beta0 = 13./12. * pow( A - 2. * B + C, 2. ) + 1./4. * pow ( A - 4. * B  + 3. * C, 2. );
	PetscReal beta1 = 13./12. * pow( B - 2. * C + D, 2. ) + 1./4. * pow ( B - D, 2. );
	PetscReal beta2 = 13./12. * pow( C - 2. * D + E, 2. ) + 1./4. * pow ( 3. * C - 4. * D  + E, 2. );

	PetscReal alpha0 = d0 / pow( eps + beta0, 2.);
	PetscReal alpha1 = d1 / pow( eps + beta1, 2.);
	PetscReal alpha2 = d2 / pow( eps + beta2, 2.);

	PetscReal sumalpha = alpha0 + alpha1 + alpha2;

	PetscReal w0 = alpha0 / sumalpha;
	PetscReal w1 = alpha1 / sumalpha;
	PetscReal w2 = alpha2 / sumalpha;

	PetscReal u0 = 2./6. * A - 7./6. * B + 11./6. * C;
	PetscReal u1 = -1./6. * B + 5./6. * C + 2./6. * D;
	PetscReal u2 = 2./6. * C + 5./6. * E - 1./6. * E;

	return w0*u0 + w1*u1 + w2*u2;

}

// ============================================================================================================= //
// DIVERGENCE SCHEME WRAPPER FUNCTIONS
// These inline functions wrap the numerical schemes to work with the function pointer interface
// ============================================================================================================= //

//! \brief Central scheme wrapper
inline Cmpnts divScheme_central(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts ucat_l, ucat_r;
    
    if(dir == 0) { // csi
        ucat_l = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        ucat_r = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
    } else if(dir == 1) { // eta
        ucat_l = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        ucat_r = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
    } else { // zet
        ucat_l = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        ucat_r = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
    }
    
    return centralVec(ucat_l, ucat_r);
}

//! \brief Central 4th order scheme wrapper
inline Cmpnts divScheme_central4(DivSchemeAccess *divaccess, PetscInt dir)
{
    if(dir == 0) { // csi
        return centralVec4thCsi(divaccess->mesh, divaccess->k, divaccess->j, divaccess->i, divaccess->mx, 
                                divaccess->nvert, divaccess->meshTag, divaccess->ucat, divaccess->wavespeed, divaccess->hyperVisc);
    } else if(dir == 1) { // eta
        return centralVec4thEta(divaccess->mesh, divaccess->k, divaccess->j, divaccess->i, divaccess->my, 
                                divaccess->nvert, divaccess->meshTag, divaccess->ucat, divaccess->wavespeed, divaccess->hyperVisc);
    } else { // zet
        return centralVec4thZet(divaccess->mesh, divaccess->k, divaccess->j, divaccess->i, divaccess->mz, 
                                divaccess->nvert, divaccess->meshTag, divaccess->ucat, divaccess->wavespeed, divaccess->hyperVisc);
    }
}

//! \brief Central-Upwind scheme wrapper
inline Cmpnts divScheme_centralUpwind(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts f0, f1, f2, f3;
    PetscReal limiterVal;
    
    if(dir == 0) { // csi
        f0 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iL];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
        f3 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iR];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].x;
    } else if(dir == 1) { // eta
        f0 = divaccess->ucat[divaccess->k][divaccess->jL][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
        f3 = divaccess->ucat[divaccess->k][divaccess->jR][divaccess->i];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].y;
    } else { // zet
        f0 = divaccess->ucat[divaccess->kL][divaccess->j][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
        f3 = divaccess->ucat[divaccess->kR][divaccess->j][divaccess->i];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].z;
    }
    
    return centralUpwindVec(f0, f1, f2, f3, divaccess->wavespeed, limiterVal);
}

//! \brief Central-Upwind weighted scheme wrapper
inline Cmpnts divScheme_centralUpwindW(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts f0, f1, f2, f3;
    PetscReal d0, d1, d2, d3, limiterVal;
    
    if(dir == 0) { // csi
        d0 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j][divaccess->iL ] * nMag(divaccess->csi[divaccess->k][divaccess->j][divaccess->iL ]));
        d1 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j][divaccess->i  ] * nMag(divaccess->csi[divaccess->k][divaccess->j][divaccess->i  ]));
        d2 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j][divaccess->i+1] * nMag(divaccess->csi[divaccess->k][divaccess->j][divaccess->i+1]));
        d3 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j][divaccess->iR ] * nMag(divaccess->csi[divaccess->k][divaccess->j][divaccess->iR ]));
        f0 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iL];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
        f3 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iR];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].x;
    } else if(dir == 1) { // eta
        d0 = 1.0 / (divaccess->aj[divaccess->k][divaccess->jL ][divaccess->i] * nMag(divaccess->eta[divaccess->k][divaccess->jL ][divaccess->i]));
        d1 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j  ][divaccess->i] * nMag(divaccess->eta[divaccess->k][divaccess->j  ][divaccess->i]));
        d2 = 1.0 / (divaccess->aj[divaccess->k][divaccess->j+1][divaccess->i] * nMag(divaccess->eta[divaccess->k][divaccess->j+1][divaccess->i]));
        d3 = 1.0 / (divaccess->aj[divaccess->k][divaccess->jR ][divaccess->i] * nMag(divaccess->eta[divaccess->k][divaccess->jR ][divaccess->i]));
        f0 = divaccess->ucat[divaccess->k][divaccess->jL][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
        f3 = divaccess->ucat[divaccess->k][divaccess->jR][divaccess->i];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].y;
    } else { // zet
        d0 = 1.0 / (divaccess->aj[divaccess->kL ][divaccess->j][divaccess->i] * nMag(divaccess->zet[divaccess->kL ][divaccess->j][divaccess->i]));
        d1 = 1.0 / (divaccess->aj[divaccess->k  ][divaccess->j][divaccess->i] * nMag(divaccess->zet[divaccess->k  ][divaccess->j][divaccess->i]));
        d2 = 1.0 / (divaccess->aj[divaccess->k+1][divaccess->j][divaccess->i] * nMag(divaccess->zet[divaccess->k+1][divaccess->j][divaccess->i]));
        d3 = 1.0 / (divaccess->aj[divaccess->kR ][divaccess->j][divaccess->i] * nMag(divaccess->zet[divaccess->kR ][divaccess->j][divaccess->i]));
        f0 = divaccess->ucat[divaccess->kL][divaccess->j][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
        f3 = divaccess->ucat[divaccess->kR][divaccess->j][divaccess->i];
        limiterVal = divaccess->limiter[divaccess->k][divaccess->j][divaccess->i].z;
    }
    
    return wCentralUpwindVec(f0, f1, f2, f3, d0, d1, d2, d3, divaccess->wavespeed, limiterVal);
}

//! \brief QUICK scheme wrapper
inline Cmpnts divScheme_quick(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts f0, f1, f2, f3;
    
    if(dir == 0) { // csi
        f0 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iL];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
        f3 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iR];
    } else if(dir == 1) { // eta
        f0 = divaccess->ucat[divaccess->k][divaccess->jL][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
        f3 = divaccess->ucat[divaccess->k][divaccess->jR][divaccess->i];
    } else { // zet
        f0 = divaccess->ucat[divaccess->kL][divaccess->j][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
        f3 = divaccess->ucat[divaccess->kR][divaccess->j][divaccess->i];
    }
    
    return quadraticUpwindVec(f0, f1, f2, f3, divaccess->wavespeed);
}

//! \brief WENO3 scheme wrapper
inline Cmpnts divScheme_weno3(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts f0, f1, f2, f3;
    
    if(dir == 0) { // csi
        f0 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iL];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
        f3 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->iR];
    } else if(dir == 1) { // eta
        f0 = divaccess->ucat[divaccess->k][divaccess->jL][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
        f3 = divaccess->ucat[divaccess->k][divaccess->jR][divaccess->i];
    } else { // zet
        f0 = divaccess->ucat[divaccess->kL][divaccess->j][divaccess->i];
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
        f3 = divaccess->ucat[divaccess->kR][divaccess->j][divaccess->i];
    }
    
    return weno3Vec(f0, f1, f2, f3, divaccess->wavespeed);
}

//! \brief Upwind scheme wrapper (1st order)
inline Cmpnts divScheme_upwind(DivSchemeAccess *divaccess, PetscInt dir)
{
    Cmpnts f1, f2;
    
    if(dir == 0) { // csi
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i+1];
    } else if(dir == 1) { // eta
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k][divaccess->j+1][divaccess->i];
    } else { // zet
        f1 = divaccess->ucat[divaccess->k][divaccess->j][divaccess->i];
        f2 = divaccess->ucat[divaccess->k+1][divaccess->j][divaccess->i];
    }
    
    return upwindVec(f1, f2, divaccess->wavespeed);
}

#endif
