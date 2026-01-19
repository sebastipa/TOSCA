#include "turbines_sampling.h"

//***************************************************************************************************************//

PetscErrorCode computeWindVectorsRotor(farm_ *farm)
{
    clock_           *clock = farm->access->clock;
    mesh_            *mesh  = farm->access->mesh;
    ueqn_            *ueqn  = farm->access->ueqn;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)
    Cmpnts           ***ucat;   // local vector (no ambiguity in this context)
    PetscReal        ***aj;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);
    DMDAVecGetArray(da,  mesh->lAj,   &aj);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if this processor controls this turbine
        if(wt->turbineControlled && !wt->useOpenFAST)
        {
            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // this turbine angular velocity
                Cmpnts omega_t = nScale(wt->rtrOmega, wt->omega_hat);

                // number of points in the AD mesh
                PetscInt npts_t = wt->adm.nPoints;

                // local relative velocity for this processor
                std::vector<Cmpnts> lU(npts_t);

                // global inflow velocity for this processor
                std::vector<Cmpnts> lWind(npts_t);

                // global inflow velocity for this processor
                std::vector<Cmpnts> gWind(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize temporary variables
                    lU[p].x = 0.0;
                    lU[p].y = 0.0;
                    lU[p].z = 0.0;

                    lWind[p].x = 0.0;
                    lWind[p].y = 0.0;
                    lWind[p].z = 0.0;

                    // save this point locally for speed
                    Cmpnts point_p = wt->adm.points[p];

                    // this point position from COR
                    Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                    // compute this blade point blade velocity
                    Cmpnts u_p  = nCross(omega_t, r_p);

                    // compute this blade point flow velocity
                    Cmpnts uf_p = nScale(-1.0, u_p);

                    // get the closest cell center
                    PetscInt i = wt->adm.closestCells[p].i,
                             j = wt->adm.closestCells[p].j,
                             k = wt->adm.closestCells[p].k;

                    if(wt->adm.thisPtControlled[p])
                    {
                        // get velocity at that point
                        Cmpnts uc_p = nSet(ucat[k][j][i]);

                        // now we have to sample the velocity from the background mesh,
                        // uc_p could also be behind the rotor so the estimation
                        // would be unstable. We use the velocity info to go back along
                        // the local streamline at a distance equal to uc_p*dt from the rotor,
                        // and use that as an estimate. This is supported by the fact
                        // that upon exiting this time iteration, the particle being at a distance
                        // from the rotor of uc_p*dt will likely be at this AD mesh point.

                        // reverse sign to the velocity (we go backward along streamline)
                        mScale(-1.0, uc_p);

                        // find the point at which velocity must be sampled
                        Cmpnts sample = nScale(clock->dt, uc_p);
                                        mSum(sample, point_p);

                        // find the closest cell indices to the sample point,
                        // allow for max 2 delta cells to stay in this processor
                        PetscReal  r_c_minMag = 1e20;
                        cellIds    closestCell;
                        PetscInt   k1, j1, i1;

                        for (k1=k-2; k1<k+3; k1++)
                        for (j1=j-2; j1<j+3; j1++)
                        for (i1=i-2; i1<i+3; i1++)
                        {
                            // compute distance from mesh cell to AD point
                            Cmpnts r_c = nSub(sample, cent[k1][j1][i1]);

                            // compute magnitude
                            PetscReal r_c_mag = nMag(r_c);

                            if(r_c_mag < r_c_minMag)
                            {
                                r_c_minMag = r_c_mag;
                                closestCell.i = i1;
                                closestCell.j = j1;
                                closestCell.k = k1;
                            }
                        }

                        // trilinear interpolate
                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            sample.x, sample.y, sample.z,
                            closestCell.i, closestCell.j, closestCell.k,
                            cent, ucat, uc_p
                        );

                        // compute the relative velocity at the AD point
                        Cmpnts ur_p  = nSet(uc_p);
                                       mSum(ur_p, uf_p);
                        lU[p]        = nSet(ur_p);

                        // compute the inflow wind at the AD point
                        lWind[p]     = nSet(uc_p);
                    }
                }

                MPI_Allreduce(&(lU[0]),  &(wt->adm.U[0]),  wt->adm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);
                MPI_Allreduce(&(lWind[0]), &(gWind[0]), wt->adm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                PetscReal rtrAvgMagU = 0.0;
                PetscReal areaSum    = 0.0;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // save this point locally for speed
                    Cmpnts point_p = wt->adm.points[p];

                    // this point position from COR
                    Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                    // build the blade reference frame based on rotation type:
                    // 1. if counter-clockwise: z from tip to root,
                    //    x as rotor axis from nacelle cone to back,
                    //    y blade tangent, directed as wind due to rotation that blade sees.
                    // 2. if clockwise: z from root to tip,
                    //    x as rotor axis from nacelle cone to back,
                    //    y blade tangent, directed as wind due to rotation that blade sees.
                    // This is done in order to have the wind vector lying in the positive
                    // quadrant in each case (it is a standard practice).

                    Cmpnts xb_hat, yb_hat, zb_hat;

                    if(wt->rotDir == "cw")
                    {
                        zb_hat = nUnit(r_p);
                        xb_hat = nScale(-1.0, wt->rtrAxis);
                        yb_hat = nCross(zb_hat, xb_hat);
                    }
                    else if(wt->rotDir == "ccw")
                    {
                        zb_hat = nUnit(r_p);
                                 mScale(-1.0, zb_hat);
                        xb_hat = nScale(-1.0, wt->rtrAxis);
                        yb_hat = nCross(zb_hat, xb_hat);
                    }

                    // transform the velocity in the bladed reference frame and
                    // remove radial (z) component
                    PetscReal ub_x = nDot(gWind[p], xb_hat);
                    PetscReal ub_y = nDot(gWind[p], yb_hat);

                    PetscReal dA;
                    if(!wt->useOpenFAST) dA = wt->adm.dr[p] * wt->adm.chord[p];
                    else dA = M_PI * wt->rTip * wt->rTip / (PetscReal)(wt->alm.nPoints);
                    
                    rtrAvgMagU += sqrt(ub_x*ub_x + ub_y*ub_y) * dA;

                    areaSum += dA;
                }

                wt->adm.rtrAvgMagU = rtrAvgMagU / areaSum;

                // clean memory
                std::vector<Cmpnts> ().swap(lU);
                std::vector<Cmpnts> ().swap(lWind);
                std::vector<Cmpnts> ().swap(gWind);
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // number of points in the AD mesh
                PetscInt npts_t = wt->uadm.nPoints;

                // local relative velocity for this processor
                std::vector<Cmpnts> lU(npts_t);

                // global inflow velocity for this processor
                std::vector<Cmpnts> lWind(npts_t);

                // global inflow velocity for this processor
                std::vector<Cmpnts> gWind(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize temporary variables
                    lU[p].x = 0.0;
                    lU[p].y = 0.0;
                    lU[p].z = 0.0;

                    lWind[p].x = 0.0;
                    lWind[p].y = 0.0;
                    lWind[p].z = 0.0;

                    // save this point locally for speed
                    Cmpnts point_p = wt->uadm.points[p];

                    // get the closest cell center
                    PetscInt i = wt->uadm.closestCells[p].i,
                             j = wt->uadm.closestCells[p].j,
                             k = wt->uadm.closestCells[p].k;

                    if(wt->uadm.thisPtControlled[p])
                    {
                        // get velocity at that point
                        Cmpnts uc_p = nSet(ucat[k][j][i]);

                        // reverse sign to the velocity (we go backward along streamline)
                        mScale(-1.0, uc_p);

                        // find the point at which velocity must be sampled
                        Cmpnts sample = nScale(clock->dt, uc_p);
                                        mSum(sample, point_p);

                        // find the closest cell indices to the sample point,
                        // allow for max 2 delta cells to stay in this processor
                        PetscReal  r_c_minMag = 1e20;
                        cellIds closestCell;
                        PetscInt     k1, j1, i1;

                        for (k1=k-2; k1<k+3; k1++)
                        for (j1=j-2; j1<j+3; j1++)
                        for (i1=i-2; i1<i+3; i1++)
                        {
                            // compute distance from mesh cell to AD point
                            Cmpnts r_c = nSub(sample, cent[k1][j1][i1]);

                            // compute magnitude
                            PetscReal r_c_mag = nMag(r_c);

                            if(r_c_mag < r_c_minMag)
                            {
                                r_c_minMag = r_c_mag;
                                closestCell.i = i1;
                                closestCell.j = j1;
                                closestCell.k = k1;
                            }
                        }

                        // trilinear interpolate
                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            sample.x, sample.y, sample.z,
                            closestCell.i, closestCell.j, closestCell.k,
                            cent, ucat, uc_p
                        );

                        // compute relative velocity at the AD point (UADM does not have rotation)
                        lU[p]    = nSet(uc_p);

                        // compute the inflow wind at the AD point
                        lWind[p] = nSet(uc_p);
                    }
                }

                MPI_Allreduce(&(lU[0]), &(wt->uadm.U[0]), wt->uadm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);
                MPI_Allreduce(&(lWind[0]), &(gWind[0]), wt->uadm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                PetscReal rtrAvgMagU = 0.0;
                PetscReal areaSum    = 0.0;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // save this point locally for speed
                    Cmpnts point_p = wt->uadm.points[p];

                    // this point position from COR
                    Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                    // form blade reference frame (as 'cw' in ADM since no rotation is present)
                    Cmpnts xb_hat, yb_hat, zb_hat;
                    zb_hat = nUnit(r_p);
                    xb_hat = nScale(-1.0, wt->rtrAxis);
                    yb_hat = nCross(zb_hat, xb_hat);

                    // transform the velocity in the bladed reference frame and
                    // remove radial (z) component
                    PetscReal ub_x = nDot(gWind[p], xb_hat);
                    PetscReal ub_y = nDot(gWind[p], yb_hat);

                    rtrAvgMagU += sqrt(ub_x*ub_x + ub_y*ub_y) * wt->uadm.dA[p];

                    areaSum += wt->uadm.dA[p];
                }

                wt->uadm.rtrAvgMagU = rtrAvgMagU / areaSum;

                // clean memory
                std::vector<Cmpnts> ().swap(lU);
                std::vector<Cmpnts> ().swap(lWind);
                std::vector<Cmpnts> ().swap(gWind);
            }
            // actuator line model
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // this turbine angular velocity
                Cmpnts omega_t = nScale(wt->rtrOmega, wt->omega_hat);

                // number of points in the AD mesh
                PetscInt npts_t = wt->alm.nPoints;

                // local relative velocity for this processor
                std::vector<Cmpnts> lU(npts_t);

                // local inflow velocity for this processor
                std::vector<Cmpnts> lWind(npts_t);

                if(wt->alm.sampleType == "rotorDisk")
                {
                    // loop over the AL mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        // initialize temporary variables
                        lU[p].x = 0.0;
                        lU[p].y = 0.0;
                        lU[p].z = 0.0;

                        lWind[p].x = 0.0;
                        lWind[p].y = 0.0;
                        lWind[p].z = 0.0;

                        // save this point locally for speed
                        Cmpnts point_p = wt->alm.points[p];

                        // this point position from COR
                        Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                        // compute this blade point blade velocity
                        Cmpnts u_p  = nCross(omega_t, r_p);

                        // compute this blade point flow velocity
                        Cmpnts uf_p = nScale(-1.0, u_p);

                        if(wt->alm.thisPtControlled[p])
                        {
                            // now we have to sample the velocity at the actuator line point
                            // from the background mesh, we use trilinear interpolation knowing
                            // that one of the closest cells to actuator line point is i,j,k


                            // get the closest cell center
                            PetscInt i = wt->alm.closestCells[p].i,
                                     j = wt->alm.closestCells[p].j,
                                     k = wt->alm.closestCells[p].k;

                            Cmpnts uc_p;

                            // trilinear interpolate
                            vectorPointLocalVolumeInterpolation
                            (
                                mesh,
                                point_p.x, point_p.y, point_p.z,
                                i, j, k,
                                cent, ucat, uc_p
                            );

                            // compute the relative velocity at the AL point
                            Cmpnts ur_p  = nSet(uc_p);
                                        mSum(ur_p, uf_p);
                            lU[p]        = nSet(ur_p);

                            // compute the inflow wind at the AL point
                            lWind[p]     = nSet(uc_p);
                        }
                    }

                    MPI_Allreduce(&(lU[0]),  &(wt->alm.U[0]),  wt->alm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);
                    MPI_Allreduce(&(lWind[0]), &(wt->alm.gWind[0]), wt->alm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);
                }

                else if(wt->alm.sampleType == "integral")
                {
                    // projection type
                    PetscInt projectionType = 0;

                    if(wt->alm.projectionType == "anisotropic")
                    {
                        projectionType = 1;
                    }

                    // loop over the AL mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        // initialize temporary variables
                        lWind[p].x = 0.0;
                        lWind[p].y = 0.0;
                        lWind[p].z = 0.0;

                        // save this point locally for speed
                        Cmpnts point_p = wt->alm.points[p];

                        // this point position from COR
                        Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                        // loop over the sphere cells
                        for(c=0; c<wt->nControlled; c++)
                        {
                            // cell indices
                            PetscInt i = wt->controlledCells[c].i,
                                    j = wt->controlledCells[c].j,
                                    k = wt->controlledCells[c].k;

                            // compute distance from mesh cell to AL point
                            Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                            // compute magnitude
                            PetscReal r_c_mag = nMag(r_c);

                            Cmpnts xb_hat, yb_hat, zb_hat;

                            // define the bladed coordinate system
                            if(wt->rotDir == "cw")
                            {
                                zb_hat = nUnit(r_p);
                                xb_hat = nScale(-1.0, wt->rtrAxis);
                                yb_hat = nCross(zb_hat, xb_hat);
                            }
                            else if(wt->rotDir == "ccw")
                            {
                                zb_hat = nUnit(r_p);
                                         mScale(-1.0, zb_hat);
                                xb_hat = nScale(-1.0, wt->rtrAxis);
                                yb_hat = nCross(zb_hat, xb_hat);
                            }

                            PetscReal pf, rPrj;

                            if(projectionType==0)
                            {
                                // get projection epsilon
                                PetscReal eps  = wt->eps;

                                // projection parameters
                                rPrj = eps * wt->prjNSigma;

                                // compute projection factor
                                pf
                                =
                                std::exp
                                (
                                    -(r_c_mag / eps)*
                                     (r_c_mag / eps)
                                ) /
                                (
                                    pow(eps,    3) *
                                    pow(M_PI, 1.5)
                                );

                            }

                            if(projectionType==1)
                            {

                                // get projection epsilon
                                PetscReal eps_x =     wt->alm.chord[p] * wt->eps_x,
                                          eps_y =     wt->alm.thick[p] * wt->eps_y,
                                          eps_z =     0.5 * wt->alm.chord[p];

                                // apply wind farm control
                                PetscReal wfControlPitch = 0.0;

                                // apply wind farm controller
                                if(farm->farmControlActive[t])
                                {
                                   wfControlPitch = wt->wfControlCollPitch;
                                }

                                // form a reference blade starting from the
                                // bladed frame, but rotate x and y around z by
                                // the twist angle

                                PetscReal sectionAngle = wt->alm.twist[p]*wt->deg2rad + wt->collPitch + wfControlPitch;

                                sectionAngle = M_PI/2.0 - sectionAngle;

                                Cmpnts xb_af_hat = nRot(zb_hat, xb_hat,  sectionAngle);
                                Cmpnts yb_af_hat = nRot(zb_hat, yb_hat,  sectionAngle);

                                PetscReal rcx = nDot(r_c,xb_af_hat),
                                          rcy = nDot(r_c,yb_af_hat),
                                          rcz = nDot(r_c,zb_hat);

                                // projection parameters
                                rPrj = std::max(std::max(eps_x, eps_y),eps_z) * wt->prjNSigma;

                                // compute projection factor
                                pf
                                =
                                std::exp
                                (
                                    -rcx * rcx / (eps_x * eps_x)
                                    -rcy * rcy / (eps_y * eps_y)
                                    -rcz * rcz / (eps_z * eps_z)
                                ) /
                                (
                                    eps_x * eps_y * eps_z *
                                    pow(M_PI, 1.5)
                                );
                            }

                            if(r_c_mag<rPrj)
                            {
                                mSum(lWind[p], nScale(pf/aj[k][j][i], ucat[k][j][i]));
                            }
                        }
                    }

                    MPI_Allreduce(&(lWind[0]), &(wt->alm.gWind[0]), wt->alm.nPoints*3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                    // loop over the AL mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        // save this point locally for speed
                        Cmpnts point_p = wt->alm.points[p];

                        // this point position from COR
                        Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                        // compute this blade point blade velocity
                        Cmpnts u_p  = nCross(omega_t, r_p);

                        // compute this blade point flow velocity
                        Cmpnts uf_p = nScale(-1.0, u_p);

                        wt->alm.U[p] = nSum(wt->alm.gWind[p], uf_p);
                    }
                }

                PetscReal rtrAvgMagU = 0.0;
                PetscReal areaSum    = 0.0;

                // loop over the AL mesh points
                for(p=0; p<npts_t; p++)
                {
                    // save this point locally for speed
                    Cmpnts point_p = wt->alm.points[p];

                    // this point position from COR
                    Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                    // build the blade reference frame based on rotation type:
                    // 1. if counter-clockwise: z from tip to root,
                    //    x as rotor axis from nacelle cone to back,
                    //    y blade tangent, directed as wind due to rotation that blade sees.
                    // 2. if clockwise: z from root to tip,
                    //    x as rotor axis from nacelle cone to back,
                    //    y blade tangent, directed as wind due to rotation that blade sees.
                    // This is done in order to have the wind vector lying in the positive
                    // quadrant in each case (it is a standard practice).

                    Cmpnts xb_hat, yb_hat, zb_hat;

                    if(wt->rotDir == "cw")
                    {
                        zb_hat = nUnit(r_p);
                        xb_hat = nScale(-1.0, wt->rtrAxis);
                        yb_hat = nCross(zb_hat, xb_hat);
                    }
                    else if(wt->rotDir == "ccw")
                    {
                        zb_hat = nUnit(r_p);
                                mScale(-1.0, zb_hat);
                        xb_hat = nScale(-1.0, wt->rtrAxis);
                        yb_hat = nCross(zb_hat, xb_hat);
                    }

                    // transform the velocity in the bladed reference frame and
                    // remove radial (z) component
                    PetscReal ub_x = nDot(wt->alm.gWind[p], xb_hat);
                    PetscReal ub_y = nDot(wt->alm.gWind[p], yb_hat);

                    PetscReal dA;
                    if(!wt->useOpenFAST) dA = wt->alm.dr[p] * wt->alm.chord[p];
                    else dA = M_PI * wt->rTip * wt->rTip / (PetscReal)(wt->alm.nPoints);

                    rtrAvgMagU += sqrt(ub_x*ub_x + ub_y*ub_y) * dA;

                    areaSum += dA;
                }

                wt->alm.rtrAvgMagU = rtrAvgMagU / areaSum;

                // clean memory
                std::vector<Cmpnts> ().swap(lU);
                std::vector<Cmpnts> ().swap(lWind);

            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // this processor velocity
                Cmpnts lU = nSetZero(), gU = nSetZero();

                if
                (
                    wt->afm.sampleType == "momentumTheory" ||
                    wt->afm.sampleType == "rotorDisk"
                )
                {
                    if(wt->afm.thisPtControlled)
                    {
                        // save this point locally for speed
                        Cmpnts point_p = wt->afm.point;

                        // get the closest cell center
                        PetscInt i = wt->afm.closestCell.i,
                                 j = wt->afm.closestCell.j,
                                 k = wt->afm.closestCell.k;

                        // get velocity at that point
                        Cmpnts uc_p = nSet(ucat[k][j][i]);

                        // now we have to sample the velocity from the background mesh,
                        // uc_p could also be behind the rotor so the estimation
                        // would be unstable. We use the velocity info to go back along
                        // the local streamline at a distance equal to uc_p*dt from the rotor,
                        // and use that as an estimate. This is supported by the fact
                        // that upon exiting this time iteration, the particle being at a distance
                        // from the rotor of uc_p*dt will likely be at this AF mesh point.

                        // reverse sign to the velocity (we go backward along streamline)
                        mScale(-1.0, uc_p);

                        // find the point at which velocity must be sampled
                        Cmpnts sample = nScale(clock->dt, uc_p);
                                        mSum(sample, point_p);

                        // find the closest cell indices to the sample point,
                        // allow for max 2 delta cells to stay in this processor
                        PetscReal  r_c_minMag = 1e20;
                        cellIds    closestCell;
                        PetscInt   k1, j1, i1;

                        for (k1=k-2; k1<k+3; k1++)
                        for (j1=j-2; j1<j+3; j1++)
                        for (i1=i-2; i1<i+3; i1++)
                        {
                            // make sure not to interpolate from ghost: AFM control cells could be near ground
                            if
                            (
                                (i1>0 && i1<mx-1) &&
                                (j1>0 && j1<my-1) &&
                                (k1>0 && k1<mz-1)
                            )
                            {
                                // compute distance from mesh cell to AF point
                                Cmpnts r_c = nSub(sample, cent[k1][j1][i1]);

                                // compute magnitude
                                PetscReal r_c_mag = nMag(r_c);

                                if(r_c_mag < r_c_minMag)
                                {
                                    r_c_minMag = r_c_mag;
                                    closestCell.i = i1;
                                    closestCell.j = j1;
                                    closestCell.k = k1;
                                }
                            }
                        }

                        // trilinear interpolate
                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            sample.x, sample.y, sample.z,
                            closestCell.i, closestCell.j, closestCell.k,
                            cent, ucat, uc_p
                        );

                        // compute the inflow at the AF point
                        lU = nSet(uc_p);
                    }
                    else
                    {
                        lU.x = 0.0;
                        lU.y = 0.0;
                        lU.z = 0.0;
                    }

                    MPI_Allreduce(&lU, &gU, 3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                    // this is an actual sampling so simulate turbine acquisition with filter
                    // signal filtering: 1-pole low pass recursive filter
                    PetscReal a = std::exp(-clock->dt * wt->afm.rtrUFilterFreq);

                    if(clock->it == clock->itStart)
                    {
                        wt->afm.U = nSet(gU);
                    }
                    else
                    {
                        wt->afm.U = nSum(nScale(a, wt->afm.U), nScale(1.0-a, gU));
                    }

                }
                else if (wt->afm.sampleType == "integral")
                {
                    // projection type
                    PetscInt projectionType;

                    if(wt->afm.projectionType == "anisotropic")
                    {
                        projectionType = 0;
                    }
                    else if(wt->afm.projectionType == "gaussexp")
                    {
                        projectionType = 1;
                    }

                    // save this point locally for speed
                    Cmpnts point_p = wt->afm.point;

                    // define local reference frame where
                    // xa_hat is aligned with the wind
                    // za_hat is vertical
                    // ya_hat defined to form right handed frame
                    Cmpnts xa_hat, ya_hat, za_hat;
                    xa_hat = nScale(-1.0, wt->rtrAxis);
                    za_hat = wt->twrDir;
                    ya_hat = nCross(za_hat, xa_hat);

                    // loop in sphere points
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute distance from mesh cell to AF point
                        Cmpnts r_c;
                        r_c.x = nDot(nSub(point_p, cent[k][j][i]), xa_hat);
                        r_c.y = nDot(nSub(point_p, cent[k][j][i]), ya_hat);
                        r_c.z = nDot(nSub(point_p, cent[k][j][i]), za_hat);

                        // compute projection factor
                        PetscReal pf;

                        if(projectionType==0)
                        {
                            pf
                            =
                            std::exp
                            (
                                -r_c.x * r_c.x / (wt->eps_x * wt->eps_x)
                                -r_c.y * r_c.y / (wt->eps_y * wt->eps_y)
                                -r_c.z * r_c.z / (wt->eps_z * wt->eps_z)
                            ) /
                            (
                                wt->eps_x * wt->eps_y * wt->eps_z *
                                pow(M_PI, 1.5)
                            );
                        }
                        else if(projectionType==1)
                        {
                            pf
                            =
                            (
                                std::exp(-r_c.x * r_c.x / (wt->eps_x * wt->eps_x)) /
                                (
                                    std::exp
                                    (
                                        (
                                            sqrt(r_c.y*r_c.y + r_c.z*r_c.z) - wt->r12
                                        ) / wt->flat
                                    ) + 1.0
                                )
                            )/wt->I;
                        }

                        mSum(lU, nScale(pf/aj[k][j][i], ucat[k][j][i]));
                    }

                    MPI_Allreduce(&lU, &gU, 3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                    if(clock->it != clock->itStart)
                    {
                        PetscReal a = std::exp(-clock->dt * wt->afm.rtrUFilterFreq);
                        wt->afm.U   = nSum(nScale(a, wt->afm.U), nScale(1.0-a, gU));
                    }
                    else
                    {
                        wt->afm.U = nSet(gU);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  mesh->lAj,   &aj);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeWindVectorsTower(farm_ *farm)
{
    clock_           *clock = farm->access->clock;
    mesh_            *mesh  = farm->access->mesh;
    ueqn_            *ueqn  = farm->access->ueqn;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)
    Cmpnts           ***ucat;   // local vector (no ambiguity in this context)

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if tower is modeled in this turbine
        if(wt->includeTwr && !wt->useOpenFAST)
        {
            // test if this processor controls this tower
            if(wt->twr.nControlled)
            {
                // number of points in the tower mesh
                PetscInt npts_t = wt->twr.nPoints;

                // local velocity for this processor
                std::vector<Cmpnts> lU(npts_t);

                // loop over the tower mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize to zero
                    lU[p].x = 0.0;
                    lU[p].y = 0.0;
                    lU[p].z = 0.0;

                    // save this point locally for speed
                    Cmpnts point_p = wt->twr.points[p];

                    if(wt->twr.thisPtControlled[p])
                    {
                        // get the closest cell center
                        PetscInt i = wt->twr.closestCells[p].i,
                                 j = wt->twr.closestCells[p].j,
                                 k = wt->twr.closestCells[p].k;

                        // find the point at which velocity must be sampled
                        Cmpnts sample = nSet(point_p);

                        // trilinear interpolate
                        vectorPointLocalVolumeInterpolation
                        (
                            mesh,
                            sample.x, sample.y, sample.z,
                            i, j, k,
                            cent, ucat, lU[p]
                        );
                    }
                }

                MPI_Allreduce(&(lU[0]), &(wt->twr.U[0]), wt->twr.nPoints*3, MPIU_REAL, MPIU_SUM, wt->twr.TWR_COMM);

                // subtract axial velocity
                for(p=0; p<npts_t; p++)
                {
                    PetscReal uAxMag = nDot(wt->twr.U[p], wt->twrDir);

                    Cmpnts uAx    = nScale(uAxMag, wt->twrDir);

                    mSub(wt->twr.U[p], uAx);
                }

                // clean memory
                std::vector<Cmpnts> ().swap(lU);
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeWindVectorsNacelle(farm_ *farm)
{
    clock_           *clock = farm->access->clock;
    mesh_            *mesh  = farm->access->mesh;
    ueqn_            *ueqn  = farm->access->ueqn;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)
    Cmpnts           ***ucat;   // local vector (no ambiguity in this context)

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if nacelle is modeled in this turbine
        if(wt->includeNacelle && !wt->useOpenFAST)
        {
            // test if this processor controls this nacelle
            if(wt->nac.nControlled)
            {
                // local velocity for this processor
                Cmpnts lU = nSetZero();

                // save this point locally for speed
                Cmpnts point_p = wt->nac.point;

                if(wt->nac.thisPtControlled)
                {
                    PetscInt i = wt->nac.closestCell.i,
                             j = wt->nac.closestCell.j,
                             k = wt->nac.closestCell.k;

                    // get velocity at that point
                    Cmpnts uc_p = nSet(ucat[k][j][i]);

                    // reverse sign to the velocity (we go backward along streamline)
                    mScale(-1.0, uc_p);

                    // find the point at which velocity must be sampled
                    Cmpnts sample = nScale(clock->dt, uc_p);
                                    mSum(sample, point_p);

                    // find the closest cell indices to the sample point,
                    // allow for max 2 delta cells to stay in this processor
                    PetscReal  r_c_minMag = 1e20;
                    cellIds    closestCell;
                    PetscInt   k1, j1, i1;

                    for (k1=k-2; k1<k+3; k1++)
                    for (j1=j-2; j1<j+3; j1++)
                    for (i1=i-2; i1<i+3; i1++)
                    {
                        // make sure not to interpolate from ghost: tower is usually next to boundaries
                        if
                        (
                            (i1>0 && i1<mx-1) &&
                            (j1>0 && j1<my-1) &&
                            (k1>0 && k1<mz-1)
                        )
                        {
                            // compute distance from mesh cell to AD point
                            Cmpnts r_c = nSub(sample, cent[k1][j1][i1]);

                            // compute magnitude
                            PetscReal r_c_mag = nMag(r_c);

                            if(r_c_mag < r_c_minMag)
                            {
                                r_c_minMag = r_c_mag;
                                closestCell.i = i1;
                                closestCell.j = j1;
                                closestCell.k = k1;
                            }
                        }
                    }

                    // trilinear interpolate
                    vectorPointLocalVolumeInterpolation
                    (
                        mesh,
                        sample.x, sample.y, sample.z,
                        closestCell.i, closestCell.j, closestCell.k,
                        cent, ucat, uc_p
                    );

                    lU = nSet(uc_p);
                }

                MPI_Allreduce(&lU, &(wt->nac.U), 3, MPIU_REAL, MPIU_SUM, wt->nac.NAC_COMM);
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeWindVectorsSample(farm_ *farm)
{
    clock_           *clock = farm->access->clock;
    mesh_            *mesh  = farm->access->mesh;
    ueqn_            *ueqn  = farm->access->ueqn;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)
    Cmpnts           ***ucat;   // local vector (no ambiguity in this context)

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, ueqn->lUcat, &ucat);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        upSampling *upPoints = farm->wt[t]->upPoints;

        // test if this processor controls this rig
        if(upPoints->thisRigControlled && !farm->wt[t]->useOpenFAST)
        {
            // number of points in the sample mesh
            PetscInt npts_t = upPoints->nPoints;

            // global inflow velocity for this processor
            std::vector<Cmpnts> lWind(npts_t);

            // global inflow velocity for this processor
            std::vector<Cmpnts> gWind(npts_t);

            // loop over the sample mesh points
            for(p=0; p<npts_t; p++)
            {
                // initialize to zero
                lWind[p].x = 0.0;
                lWind[p].y = 0.0;
                lWind[p].z = 0.0;

                gWind[p].x = 0.0;
                gWind[p].y = 0.0;
                gWind[p].z = 0.0;

                // save this point locally for speed
                Cmpnts point_p = upPoints->points[p];

                // get the closest cell center
                PetscInt i = upPoints->closestCells[p].i,
                         j = upPoints->closestCells[p].j,
                         k = upPoints->closestCells[p].k;

                if(upPoints->thisPtControlled[p])
                {
                    // find the point at which velocity must be sampled
                    Cmpnts sample = nSet(point_p);

                    // trilinear interpolate (do not go upstream for the sample)
                    vectorPointLocalVolumeInterpolation
                    (
                        mesh,
                        sample.x, sample.y, sample.z,
                        i, j, k,
                        cent, ucat, lWind[p]
                    );
                }
            }

            MPI_Allreduce(&(lWind[0]), &(gWind[0]), upPoints->nPoints*3, MPIU_REAL, MPIU_SUM, upPoints->UPW_COMM);

            PetscReal rtrAvgMagU = 0.0;
            PetscReal areaSum    = 0.0;

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                // save this point locally for speed
                Cmpnts point_p = upPoints->points[p];

                // this point position from COR
                Cmpnts r_p  = nSub(point_p, upPoints->center);

                // form blade reference frame (as 'cw' in ADM since no rotation is present)
                Cmpnts xb_hat, yb_hat, zb_hat;
                zb_hat = nUnit(r_p);
                xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                yb_hat = nCross(zb_hat, xb_hat);

                // transform the velocity in the bladed reference frame and
                // remove radial (z) component
                PetscReal ub_x = nDot(gWind[p], xb_hat);
                PetscReal ub_y = nDot(gWind[p], yb_hat);

                PetscReal dA   = upPoints->dA[p];

                rtrAvgMagU += sqrt(ub_x*ub_x + ub_y*ub_y) * dA;

                areaSum += dA;
            }

            // this will be on all procs belonging to UPW_COMM and TRB_COMM
            upPoints->Uref = rtrAvgMagU / areaSum;

            // clean memory
            std::vector<Cmpnts> ().swap(lWind);
            std::vector<Cmpnts> ().swap(gWind);
        }
    }

    DMDAVecRestoreArray(fda, ueqn->lUcat, &ucat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}