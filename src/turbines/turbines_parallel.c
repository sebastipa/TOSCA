#include "turbines_parallel.h"

//***************************************************************************************************************//

PetscErrorCode initControlledCells(farm_ *farm)
{
    // initialize the controlled cells labels.
    // Loop over wind turbines and check if this processor
    // controls that wind turbine by defining a sphere of cells around the
    // rotor where this processor have access. If the size of the sphere cell
    // is zero this processor does not control the wind turbine, so set
    // turbineControlled flag to zero. A wind turbine can be shared between
    // two or more processors, in that case the sphere cells of the two or more processors
    // will sum up to an actual sphere of cell labels. In that case turbineControlled is
    // set to one since must compute the body force on the sub-set of the sphere.

    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, p, c;

    PetscInt         hasCells;

    Cmpnts           ***cent;

    PetscMPIInt      rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // define communicator color
        PetscInt commColor;

        // define this turbine sphere cells
        std::vector<cellIds> sphereCells;

        // set flag determining if this proc has cells for this turbine to zero
        hasCells = 0;

        // this turbine sphere radius
        PetscReal radius_t;

        if((*farm->turbineModels[t]) != "AFM")
        {
            if((*farm->turbineModels[t]) == "ALM")
            {
                if(wt->alm.projectionType == "anisotropic")
                {
                    // use 5m as reference length for chord and thickness
                    PetscReal   eps_x =     5.0 * wt->eps_x,
                                eps_y =     5.0 * wt->eps_y;


                    PetscReal eps = PetscMax( eps_x, eps_y );

                    radius_t = wt->rTip +  eps * wt->prjNSigma;
                }
                else
                {
                    radius_t = wt->rTip +  wt->eps * wt->prjNSigma;
                }
            }
            else
            {
                radius_t = wt->rTip +  wt->eps * wt->prjNSigma;
            }
        }
        else
        {
            if(wt->afm.projectionType == "anisotropic")
            {
                radius_t = PetscMax(wt->eps_x, PetscMax(wt->eps_y, wt->eps_z)) * wt->prjNSigma;
            }
            else if(wt->afm.projectionType == "gaussexp")
            {
                radius_t = PetscMax(wt->eps_x*wt->prjNSigma, 2.0*wt->r12);
            }
        }

        // this turbine rotor center
        Cmpnts rotor_c = wt->rotCenter;

        // loop over this processor mesh points
        for(k=lzs; k<lze; k++)
        for(j=lys; j<lye; j++)
        for(i=lxs; i<lxe; i++)
        {
            // compute cell center to point distance
            Cmpnts dist    = nSub(cent[k][j][i], rotor_c);

            // compute distance magnitude
            PetscReal distMag = nMag(dist);

            // test if inside sphere
            if(distMag <= radius_t)
            {
                cellIds thisCell;

                thisCell.i = i;
                thisCell.j = j;
                thisCell.k = k;

                sphereCells.push_back(thisCell);

                // this proc has cells: activate flag
                if(hasCells == 0) hasCells = 1;
            }
        }

        // allocate memory and store in the AD model
        if(hasCells)
        {
            // set number of controlled cells
            PetscInt nControlled = sphereCells.size();
            wt->nControlled = nControlled;

            // allocate memory
            PetscMalloc(nControlled*sizeof(cellIds), &(wt->controlledCells));

            // save cell ids
            for(c=0; c<nControlled; c++)
            {
                wt->controlledCells[c].i = sphereCells[c].i;
                wt->controlledCells[c].j = sphereCells[c].j;
                wt->controlledCells[c].k = sphereCells[c].k;
            }

            // this processor controls this turbine
            wt->turbineControlled = 1;

            // set the communicator color
            commColor = 1;
        }
        else
        {
            // set number of controlled cells to zero
            wt->nControlled = 0;

            // this processor does not control this turbine
            wt->turbineControlled = 0;

            // set the communicator color to undefined
            commColor = MPI_UNDEFINED;
        }

        // clean memory
        std::vector<cellIds> ().swap(sphereCells);

        // print processor control if debug
        if(wt->dbg)
        {
            printf("    > rank %d: turbine %s controlled flag = %ld\n", rank, (*farm->turbineIds[t]).c_str(), wt->turbineControlled);
            MPI_Barrier(mesh->MESH_COMM);
        }

        // create communicator
        MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(wt->TRB_COMM));

        // the master rank of this TRB_COMM communicator will write this turbine I/O file
        PetscMPIInt thisTrbRank = 10, lwriterRank;
        MPI_Comm_rank(wt->TRB_COMM, &thisTrbRank);

        if(!thisTrbRank && wt->turbineControlled)
        {
            lwriterRank = rank;
        }
        else
        {
            lwriterRank = 0;
        }

        // scatter this info among all processors in the mesh->MESH_COMM
        MPI_Allreduce(&lwriterRank, &(wt->writerRank), 1, MPI_INT, MPI_SUM, mesh->MESH_COMM);

        // get the number of processor that control this turbine to average turbine data
        PetscMPIInt lnProcsTrb, gnProcsTrb;
        if(!thisTrbRank && wt->turbineControlled)
        {
            MPI_Comm_size(wt->TRB_COMM, &lnProcsTrb);
        }
        else
        {
            lnProcsTrb = 0;
        }

        // scatter this info among all processors in the mesh->MESH_COMM
        MPI_Allreduce(&lnProcsTrb, &gnProcsTrb, 1, MPI_INT, MPI_SUM, mesh->MESH_COMM);

        wt->nProcsTrb = (PetscReal)gnProcsTrb;

        if(wt->dbg && rank == wt->writerRank)
        {
            printf("    > turbine %s: writer rank = %d, controlling ranks = %d\n", (*farm->turbineIds[t]).c_str(), rank, gnProcsTrb);
        }
        MPI_Barrier(mesh->MESH_COMM);
    }

    // initialize the tower cells from which velocity is sampled
    // loop over each tower point and search the cells wich are
    // inside the projection radius.

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if(wt->includeTwr)
        {
            // define communicator color
            PetscInt commColor;

            // cells where the velocity will be affected by the tower
            std::vector<cellIds> towerCells;

            // set this tower projection radius
            // projection parameters
            PetscReal rPrj = wt->twr.eps * wt->twr.prjNSigma;

            // set flag determining if this proc has cells for this turbine to zero
            PetscInt hasCells = 0;

            // set tower base point
            Cmpnts base = wt->twr.points[0];

            // compute a dummy reference frame as follows:
            // x: in the non tilted rotor direction
            // y: as the tower from bottom to top
            // z: right hand rule

            Cmpnts xd_hat = wt->rtrAxis;
            Cmpnts yd_hat = wt->twrDir;
            Cmpnts zd_hat = nCross(xd_hat, yd_hat);

            // loop over this processor mesh points
            for(k=lzs; k<lze; k++)
            for(j=lys; j<lye; j++)
            for(i=lxs; i<lxe; i++)
            {
                // compute distance from point to wind turbine base
                Cmpnts dist = nSub(cent[k][j][i], base);

                // compute distance components in the dummy directions
                PetscReal xdist = fabs(nDot(dist, xd_hat));
                PetscReal ydist = fabs(nDot(dist, yd_hat));
                PetscReal zdist = fabs(nDot(dist, zd_hat));

                // compute radial surface parallel distance from tower
                PetscReal radDist = sqrt(xdist*xdist+zdist*zdist);

                // test if inside a cylinder which has the tower at the center
                if
                (
                    radDist < wt->hTwr / 2 + rPrj &&
                    ydist < wt->hTwr + rPrj
                )
                {
                    cellIds thisCell;

                    thisCell.i = i;
                    thisCell.j = j;
                    thisCell.k = k;

                    towerCells.push_back(thisCell);

                    // this proc has cells: activate flag
                    if(hasCells == 0) hasCells = 1;
                }
            }

            // if this processor has tower cells allocate and save
            if(hasCells)
            {
                // set number of controlled cells
                PetscInt nControlled = towerCells.size();
                wt->twr.nControlled = nControlled;

                // allocate memory
                PetscMalloc(nControlled*sizeof(cellIds), &(wt->twr.controlledCells));

                // save cell ids
                for(c=0; c<nControlled; c++)
                {
                    wt->twr.controlledCells[c].i = towerCells[c].i;
                    wt->twr.controlledCells[c].j = towerCells[c].j;
                    wt->twr.controlledCells[c].k = towerCells[c].k;
                }

                // set the communicator color
                commColor = 1;
            }
            else
            {
                // set number of controlled cells to zero
                wt->twr.nControlled = 0;

                // set the communicator color
                commColor = MPI_UNDEFINED;
            }

            // clean memory
            std::vector<cellIds> ().swap(towerCells);

            // create communicator
            MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(wt->twr.TWR_COMM));

            PetscMPIInt thisTwrRank = 10;
            MPI_Comm_rank(wt->twr.TWR_COMM, &thisTwrRank);

            // get the number of processor that control this turbine to average turbine data
            PetscMPIInt lnProcsTwr, gnProcsTwr;
            if(!thisTwrRank && wt->twr.nControlled)
            {
                MPI_Comm_size(wt->twr.TWR_COMM, &lnProcsTwr);
            }
            else
            {
                lnProcsTwr = 0;
            }

            // scatter this info among all processors in the mesh->MESH_COMM
            MPI_Allreduce(&lnProcsTwr, &gnProcsTwr, 1, MPI_INT, MPI_SUM, mesh->MESH_COMM);

            wt->twr.nProcsTwr = (PetscReal)gnProcsTwr;

            if(wt->dbg && rank == wt->writerRank)
            {
                printf("    > turbine %s: controlling tower ranks = %d\n", (*farm->turbineIds[t]).c_str(), gnProcsTwr);
            }

            MPI_Barrier(mesh->MESH_COMM);
        }

        if(wt->includeNacelle)
        {
            // define communicator color
            PetscInt commColor;

            // cells where the velocity will be affected by the nacelle
            std::vector<cellIds> nacelleCells;

            // set this nacelle projection radius
            PetscReal rPrj = std::min(wt->nac.eps * wt->nac.prjNSigma, wt->hTwr);

            // set flag determining if this proc has cells for this nacelle to zero
            PetscInt hasCells = 0;

            // this turbine rotor center
            Cmpnts rotor_c = wt->rotCenter;

            // loop over this processor mesh points
            for(k=lzs; k<lze; k++)
            for(j=lys; j<lye; j++)
            for(i=lxs; i<lxe; i++)
            {
                // compute distance from point to wind turbine base
                Cmpnts dist = nSub(cent[k][j][i], rotor_c);

                // compute distance magnitude
                PetscReal distMag = nMag(dist);

                // test if inside sphere
                if(distMag < rPrj)
                {
                    cellIds thisCell;

                    thisCell.i = i;
                    thisCell.j = j;
                    thisCell.k = k;

                    nacelleCells.push_back(thisCell);

                    // this proc has cells: activate flag
                    if(hasCells == 0) hasCells = 1;
                }
            }

            // if this processor has tower cells allocate and save
            if(hasCells)
            {
                // set number of controlled cells
                PetscInt nControlled = nacelleCells.size();
                wt->nac.nControlled  = nControlled;

                // allocate memory
                PetscMalloc(nControlled*sizeof(cellIds), &(wt->nac.controlledCells));

                // save cell ids
                for(c=0; c<nControlled; c++)
                {
                    wt->nac.controlledCells[c].i = nacelleCells[c].i;
                    wt->nac.controlledCells[c].j = nacelleCells[c].j;
                    wt->nac.controlledCells[c].k = nacelleCells[c].k;
                }

                // set the communicator color
                commColor = 1;
            }
            else
            {
                // set number of controlled cells to zero
                wt->nac.nControlled = 0;

                // set the communicator color
                commColor = MPI_UNDEFINED;
            }

            // clean memory
            std::vector<cellIds> ().swap(nacelleCells);

            // create communicator
            MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(wt->nac.NAC_COMM));

            PetscMPIInt thisNacRank = 10;
            MPI_Comm_rank(wt->nac.NAC_COMM, &thisNacRank);

            // get the number of processor that control this turbine to average turbine data
            PetscMPIInt lnProcsNac, gnProcsNac;
            if(!thisNacRank && wt->nac.nControlled)
            {
                MPI_Comm_size(wt->nac.NAC_COMM, &lnProcsNac);
            }
            else
            {
                lnProcsNac = 0;
            }

            // scatter this info among all processors in the mesh->MESH_COMM
            MPI_Allreduce(&lnProcsNac, &gnProcsNac, 1, MPI_INT, MPI_SUM, mesh->MESH_COMM);

            wt->nac.nProcsNac = (PetscReal)gnProcsNac;

            if(wt->dbg && rank == wt->writerRank)
            {
                printf("    > turbine %s: controlling nacelle ranks = %d\n", (*farm->turbineIds[t]).c_str(), gnProcsNac);
            }

            MPI_Barrier(mesh->MESH_COMM);
        }
    }

    // check that there are tower and nacelle cells if tower and nacelle are active, otherwise throws
    // errors and suggests to increase epsilon
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if(wt->includeTwr)
        {
            PetscInt lnCells, gnCells;

            lnCells = wt->twr.nControlled;

            // scatter on MESH_COMM because there might not be a TWR_COMM if no cells were selected
            MPI_Allreduce(&lnCells, &gnCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            if(gnCells == 0)
            {
                char error[512];
                sprintf(error, "could not find any tower cells when searching, try increase tower epsilon\n");
                fatalErrorInFunction("initControlledCells",  error);
            }
        }

        if(wt->includeNacelle)
        {
            PetscInt lnCells, gnCells;

            lnCells = wt->nac.nControlled;

            // scatter on MESH_COMM because there might not be a TWR_COMM if no cells were selected
            MPI_Allreduce(&lnCells, &gnCells, 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            if(gnCells == 0)
            {
                char error[512];
                sprintf(error, "could not find any nacelle cells when searching, try increase nacelle epsilon\n");
                fatalErrorInFunction("initControlledCells",  error);
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initSampleControlledCells(farm_ *farm)
{
    // initialize the controlled cells labels.
    // Loop over wind turbine sample points and check if this processor
    // controls that wind turbine by defining a sphere of cells around the
    // rotor where this processor have access. If the size of the sphere cell
    // is zero, this processor does not control any sample point, so set
    // turbineControlled flag to zero. A wind turbine can be shared between
    // two processors, in that case the sphere cells of the two or more processors
    // will sum up to an actual sphere of cell labels. In that case turbineControlled is
    // set to one since must compute the body force on the sub-set of the sphere.

    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, p, c;

    PetscInt         hasCells;

    Cmpnts           ***cent;

    PetscMPIInt      rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        // define this turbine communicator color
        PetscInt commColor;

        // define this turbine sphere cells
        std::vector<cellIds> sphereCells;

        // set a pointer to this turbine sample points
        upSampling* upPoints = farm->wt[t]->upPoints;

        // set flag determining if this proc has cells for this turbine to zero
        hasCells = 0;

        // this sphere radius
        PetscReal radius_t = 3.0 * (2.0*farm->wt[t]->rTip);

        // this turbine rotor center
        Cmpnts rotor_c = farm->wt[t]->rotCenter;

        // loop over this processor mesh points
        for(k=lzs; k<lze; k++)
        for(j=lys; j<lye; j++)
        for(i=lxs; i<lxe; i++)
        {
            // compute cell center to point distance
            Cmpnts dist    = nSub(cent[k][j][i], rotor_c);

            // compute distance magnitude
            PetscReal distMag = nMag(dist);

            // test if inside sphere
            if(distMag < radius_t)
            {
                cellIds thisCell;

                thisCell.i = i;
                thisCell.j = j;
                thisCell.k = k;

                sphereCells.push_back(thisCell);

                // this proc has cells: activate flag
                if(hasCells == 0) hasCells = 1;
            }
        }

        // allocate memory and store in the AD model
        if(hasCells)
        {
            // set number of controlled cells
            PetscInt nControlled = sphereCells.size();
            upPoints->nControlled = nControlled;

            // allocate memory
            PetscMalloc(nControlled*sizeof(cellIds), &(upPoints->controlledCells));

            // save cell ids
            for(c=0; c<nControlled; c++)
            {
                upPoints->controlledCells[c].i = sphereCells[c].i;
                upPoints->controlledCells[c].j = sphereCells[c].j;
                upPoints->controlledCells[c].k = sphereCells[c].k;
            }

            // this processor controls this sampling rig
            upPoints->thisRigControlled = 1;

            // set the communicator color
            commColor = 1;
        }
        else
        {
            // set number of controlled cells to zero
            upPoints->nControlled = 0;

            // this processor does not control this sampling rig
            upPoints->thisRigControlled = 0;

            // set the communicator color
            commColor = MPI_UNDEFINED;
        }

        // clean memory
        std::vector<cellIds> ().swap(sphereCells);

        // create communicator
        MPI_Comm_split(mesh->MESH_COMM, commColor, rank, &(upPoints->UPW_COMM));
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledPointsRotor(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // create point perturbation range for breaking ties in determining the
    // owner processor. Note: the actual position of the points is not changed,
    // this is only to make sure that each AD points is only controlled by a
    // single processor.

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM" && wt->yawChanged)
            {
                // number of points in the AD mesh
                PetscInt npts_t = wt->adm.nPoints;

                // create temporary vectors
                std::vector<PetscReal> lminDist(npts_t);
                std::vector<PetscReal> gminDist(npts_t);
                std::vector<Cmpnts>    perturb(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize min dists to a big value
                    lminDist[p] = 1e20;
                    gminDist[p] = 1e20;

                    // set point perturbation
                    perturb[p].x =  procContrib;
                    perturb[p].y =  procContrib;
                    perturb[p].z =  procContrib;

                    // save this point locally for speed
                    Cmpnts point_p = wt->adm.points[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute distance from mesh cell to AD point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute magnitude
                        PetscReal r_c_mag = nMag(r_c);

                        if(r_c_mag < r_c_minMag)
                        {
                            r_c_minMag = r_c_mag;
                            closestCell.i = i;
                            closestCell.j = j;
                            closestCell.k = k;
                        }
                    }

                    // save closest cell indices
                    wt->adm.closestCells[p].i = closestCell.i;
                    wt->adm.closestCells[p].j = closestCell.j;
                    wt->adm.closestCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->adm.nPoints, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

                // now compare the lists, if they have the same min distance then the point
                // is controlled by this processor, otherwise not. Ties are break by making
                // the perturbation different among the processors:
                // distance was not equal: it can't become equal
                // distance was equal    : it is made different
                for(p=0; p<npts_t; p++)
                {
                    // point is controlled
                    if(lminDist[p] == gminDist[p])
                    {
                        wt->adm.thisPtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->adm.thisPtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts> ().swap(perturb);
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM" && wt->yawChanged)
            {
                // number of points in the AD mesh
                PetscInt npts_t = wt->uadm.nPoints;

                // create temporary vectors
                std::vector<PetscReal> lminDist(npts_t);
                std::vector<PetscReal> gminDist(npts_t);
                std::vector<Cmpnts> perturb(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize min dists to a big value
                    lminDist[p] = 1e20;
                    gminDist[p] = 1e20;

                    // set point perturbation
                    perturb[p].x =  procContrib;
                    perturb[p].y =  procContrib;
                    perturb[p].z =  procContrib;

                    // save this point locally for speed
                    Cmpnts point_p = wt->uadm.points[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute distance from mesh cell to AD point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute magnitude
                        PetscReal r_c_mag = nMag(r_c);

                        if(r_c_mag < r_c_minMag)
                        {
                            r_c_minMag = r_c_mag;
                            closestCell.i = i;
                            closestCell.j = j;
                            closestCell.k = k;
                        }
                    }

                    // save closest cell indices
                    wt->uadm.closestCells[p].i = closestCell.i;
                    wt->uadm.closestCells[p].j = closestCell.j;
                    wt->uadm.closestCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->uadm.nPoints, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

                // now compare the lists, if they have the same min distance then the point
                // is controlled by this processor, otherwise not. Ties are break by making
                // the perturbation different among the processors:
                // distance was not equal: it can't become equal
                // distance was equal    : it is made different
                for(p=0; p<npts_t; p++)
                {
                    // point is controlled
                    if(lminDist[p] == gminDist[p])
                    {
                        wt->uadm.thisPtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->uadm.thisPtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts> ().swap(perturb);

            }
            // actuator line model (always do this as blades rotate)
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // number of points in the AL mesh
                PetscInt npts_t = wt->alm.nPoints;

                // create temporary vectors
                std::vector<PetscReal> lminDist(npts_t);
                std::vector<PetscReal> gminDist(npts_t);
                std::vector<Cmpnts> perturb(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // initialize min dists to a big value
                    lminDist[p] = 1e20;
                    gminDist[p] = 1e20;

                    // set point perturbation
                    perturb[p].x =  procContrib;
                    perturb[p].y =  procContrib;
                    perturb[p].z =  procContrib;

                    // save this point locally for speed
                    Cmpnts point_p = wt->alm.points[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute distance from mesh cell to AD point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute magnitude
                        PetscReal r_c_mag = nMag(r_c);

                        if(r_c_mag < r_c_minMag)
                        {
                            r_c_minMag = r_c_mag;
                            closestCell.i = i;
                            closestCell.j = j;
                            closestCell.k = k;
                        }
                    }

                    // save closest cell indices
                    wt->alm.closestCells[p].i = closestCell.i;
                    wt->alm.closestCells[p].j = closestCell.j;
                    wt->alm.closestCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                // this call can be in the turbineControlled test as long as the communicator is TRB_COMM (will hang otherwise)
                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->alm.nPoints, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

                // now compare the lists, if they have the same min distance then the point
                // is controlled by this processor, otherwise not. Ties are break by making
                // the perturbation different among the processors:
                // distance was not equal: it can't become equal
                // distance was equal    : it is made different
                for(p=0; p<npts_t; p++)
                {
                    // point is controlled
                    if(lminDist[p] == gminDist[p])
                    {
                        wt->alm.thisPtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->alm.thisPtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts> ().swap(perturb);
            }
            // actuator farm model
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // do the search only once even if yaw changes for this
                if(!wt->afm.searchDone)
                {
                    Cmpnts perturbVec;
                           perturbVec.x = procContrib;
                           perturbVec.y = procContrib;
                           perturbVec.z = procContrib;

                    // we must find the closest cell to the AF point
                    Cmpnts rotor_c = wt->afm.point;

                    PetscReal gmindist, lmindist = 1e10;
                    cellIds   closestCell;

                    // loop over the sphere cells
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute cell center to point distance
                        Cmpnts dist = nSub(cent[k][j][i], rotor_c);
                                      mSub(dist, perturbVec);

                        // compute distance magnitude
                        PetscReal distMag = nMag(dist);

                        // test if inside sphere
                        if(distMag < lmindist)
                        {
                            closestCell.i = i;
                            closestCell.j = j;
                            closestCell.k = k;

                            lmindist = distMag;
                        }
                    }

                    // find the overall minimum distance
                    MPI_Allreduce(&lmindist, &gmindist, 1, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

                    // compare the lists, where they agree the point is in this processor
                    if(lmindist == gmindist)
                    {
                        // set closest cell indices
                        wt->afm.closestCell.i = closestCell.i;
                        wt->afm.closestCell.j = closestCell.j;
                        wt->afm.closestCell.k = closestCell.k;

                        // set controlled flag
                        wt->afm.thisPtControlled = 1;
                    }
                    else
                    {
                        // set controlled flag
                        wt->afm.thisPtControlled = 0;
                    }

                    wt->afm.searchDone = 1;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    // check that discrimination algorithm worked properly
    checkPointDiscriminationRotor(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledPointsTower(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // create point perturbation range for breaking ties in determining the
    // owner processor. Note: the actual position of the points is not changed,
    // this is only to make sure that each twr point is only controlled by a
    // single processor.

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if tower is modeled in this turbine
        if(wt->includeTwr)
        {
            // test if this processor controls this turbine
            if(wt->twr.nControlled)
            {
                // number of points in the tower mesh
                PetscInt npts_t = wt->twr.nPoints;

                // create temporary vectors
                std::vector<PetscReal> lminDist(npts_t);
                std::vector<PetscReal> gminDist(npts_t);
                std::vector<Cmpnts>    perturb(npts_t);

                // loop over the tower mesh points
                for(p=0; p<npts_t; p++)
                {
                    // set min dists to a big value
                    lminDist[p] = 1e20;
                    gminDist[p] = 1e20;

                    // set point perturbation
                    perturb[p].x =  procContrib;
                    perturb[p].y =  procContrib;
                    perturb[p].z =  procContrib;

                    // save this point locally for speed
                    Cmpnts point_p = wt->twr.points[p];

                    // perturb the point position
                    mSum(point_p, perturb[p]);

                    // find the closest cell center
                    PetscReal  r_c_minMag = 1e20;
                    cellIds closestCell;

                    // loop over the tower cells
                    for(c=0; c<wt->twr.nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->twr.controlledCells[c].i,
                                 j = wt->twr.controlledCells[c].j,
                                 k = wt->twr.controlledCells[c].k;

                        // compute distance from mesh cell to AD point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute magnitude
                        PetscReal r_c_mag = nMag(r_c);

                        if(r_c_mag < r_c_minMag)
                        {
                            r_c_minMag = r_c_mag;
                            closestCell.i = i;
                            closestCell.j = j;
                            closestCell.k = k;
                        }
                    }

                    // save closest cell indices
                    wt->twr.closestCells[p].i = closestCell.i;
                    wt->twr.closestCells[p].j = closestCell.j;
                    wt->twr.closestCells[p].k = closestCell.k;

                    // save min dist
                    lminDist[p] = r_c_minMag;
                }

                MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->twr.nPoints, MPIU_REAL, MPIU_MIN, wt->twr.TWR_COMM);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // point is controlled
                    if(lminDist[p] == gminDist[p])
                    {
                        wt->twr.thisPtControlled[p] = 1;
                    }
                    // point is not controlled
                    else
                    {
                        wt->twr.thisPtControlled[p] = 0;
                    }
                }

                // clean memory
                std::vector<PetscReal> ().swap(lminDist);
                std::vector<PetscReal> ().swap(gminDist);
                std::vector<Cmpnts>    ().swap(perturb);
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    // check that discrimination algorithm worked
    checkPointDiscriminationTower(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledPointsNacelle(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;   // local vector (no ambiguity in this context)

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // create point perturbation range for breaking ties in determining the
    // owner processor. Note: the actual position of the points is not changed,
    // this is only to make sure that the nac point is only controlled by a
    // single processor.

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if nacelle is modeled in this turbine
        if(wt->includeNacelle)
        {
            // test if this processor controls this turbine
            if(wt->nac.nControlled)
            {
                // create temporary vectors
                PetscReal lminDist = 1e20;
                PetscReal gminDist = 1e20;
                Cmpnts    perturb  = nSetFromComponents(procContrib,procContrib,procContrib);

                // save this point locally
                Cmpnts point_p = wt->nac.point;

                // perturb the point position
                mSum(point_p, perturb);

                // find the closest cell center
                PetscReal  r_c_minMag = 1e20;
                cellIds    closestCell;

                // loop over the tower cells
                for(c=0; c<wt->nac.nControlled; c++)
                {
                    // cell indices
                    PetscInt i = wt->nac.controlledCells[c].i,
                             j = wt->nac.controlledCells[c].j,
                             k = wt->nac.controlledCells[c].k;

                    // compute distance from mesh cell to AD point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    // compute magnitude
                    PetscReal r_c_mag = nMag(r_c);

                    if(r_c_mag < r_c_minMag)
                    {
                        r_c_minMag = r_c_mag;
                        closestCell.i = i;
                        closestCell.j = j;
                        closestCell.k = k;
                    }
                }

                // save closest cell indices
                wt->nac.closestCell.i = closestCell.i;
                wt->nac.closestCell.j = closestCell.j;
                wt->nac.closestCell.k = closestCell.k;

                // save min dist
                lminDist = r_c_minMag;

                MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, wt->nac.NAC_COMM);

                // point is controlled
                if(lminDist == gminDist)
                {
                    wt->nac.thisPtControlled = 1;
                }
                // point is not controlled
                else
                {
                    wt->nac.thisPtControlled = 0;
                }
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    // check that discrimination algorithm worked
    checkPointDiscriminationNacelle(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode findControlledPointsSample(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         t, c, p;

    Cmpnts           ***cent;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // create point perturbation range for breaking ties in determining the
    // owner processor. Note: the actual position of the points is not changed,
    // this is only to make sure that each sample point is only controlled by a
    // single processor.

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    DMDAVecGetArray(fda, mesh->lCent, &cent);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        upSampling *upPoints = farm->wt[t]->upPoints;

        // test if this processor controls this rig
        if(upPoints->thisRigControlled && farm->wt[t]->yawChanged)
        {
            // number of points in the sampling mesh
            PetscInt npts_t = upPoints->nPoints;

            // create temporary vectors
            std::vector<PetscReal> lminDist(npts_t);
            std::vector<PetscReal> gminDist(npts_t);
            std::vector<Cmpnts>    perturb(npts_t);

            // loop over the sample mesh points
            for(p=0; p<npts_t; p++)
            {
                // set min dists to a big value
                lminDist[p] = 1e20;
                gminDist[p] = 1e20;

                // set point perturbation
                perturb[p].x =  procContrib;
                perturb[p].y =  procContrib;
                perturb[p].z =  procContrib;

                // save this point locally for speed
                Cmpnts point_p = upPoints->points[p];

                // perturb the point position
                mSum(point_p, perturb[p]);

                // find the closest cell center
                PetscReal  r_c_minMag = 1e20;
                cellIds    closestCell;

                // loop over the sphere cells
                for(c=0; c<upPoints->nControlled; c++)
                {
                    // cell indices
                    PetscInt i = upPoints->controlledCells[c].i,
                             j = upPoints->controlledCells[c].j,
                             k = upPoints->controlledCells[c].k;

                    // compute distance from mesh cell to sample point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    // compute magnitude
                    PetscReal r_c_mag = nMag(r_c);

                    if(r_c_mag < r_c_minMag)
                    {
                        r_c_minMag = r_c_mag;
                        closestCell.i = i;
                        closestCell.j = j;
                        closestCell.k = k;
                    }
                }

                // save closest cell indices
                upPoints->closestCells[p].i = closestCell.i;
                upPoints->closestCells[p].j = closestCell.j;
                upPoints->closestCells[p].k = closestCell.k;

                // save min dist
                lminDist[p] = r_c_minMag;
            }

            MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), upPoints->nPoints, MPIU_REAL, MPIU_MIN, upPoints->UPW_COMM);

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                // point is controlled
                if(lminDist[p] == gminDist[p])
                {
                    upPoints->thisPtControlled[p] = 1;
                }
                // point is not controlled
                else
                {
                    upPoints->thisPtControlled[p] = 0;
                }
            }

            // clean memory
            std::vector<PetscReal> ().swap(lminDist);
            std::vector<PetscReal> ().swap(gminDist);
            std::vector<Cmpnts>    ().swap(perturb);
        }

        // set yaw changed to zero (will be overwritten by controlNacYaw if applicable)
        farm->wt[t]->yawChanged = 0;
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    checkPointDiscriminationSample(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkPointDiscriminationRotor(farm_ *farm)
{
    if(farm->dbg)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscInt t, c, p;

        // test for points accounted for twice or none
        std::vector<std::vector<PetscInt>> laccounted(farm->size);
        std::vector<std::vector<PetscInt>> gaccounted(farm->size);

        // set to zero
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];
            PetscInt npts_t;

            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // number of points in the AD mesh
                npts_t = wt->adm.nPoints;
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // number of points in the AD mesh
                npts_t = wt->uadm.nPoints;
            }
            // actuator line model
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // number of points in the AD mesh
                npts_t = wt->alm.nPoints;
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // number of points in the AD mesh
                npts_t = 1;
            }

            laccounted[t].resize(npts_t);
            gaccounted[t].resize(npts_t);

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                laccounted[t][p] = 0;
                gaccounted[t][p] = 0;
            }
        }

        // loop over each wind turbine
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // test if this processor controls this turbine
            if(wt->turbineControlled)
            {
                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    // number of points in the AD mesh
                    PetscInt npts_t = wt->adm.nPoints;

                    // loop over the AD mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->adm.thisPtControlled[p])
                        {
                            // set this point as accounted for. This is to check that
                            // the discrimination algorithm works properly.
                            laccounted[t][p]++;
                        }
                    }
                }
                // uniform actuator disk model
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    // number of points in the AD mesh
                    PetscInt npts_t = wt->uadm.nPoints;

                    // loop over the AD mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->uadm.thisPtControlled[p])
                        {
                            // set this point as accounted for. This is to check that
                            // the discrimination algorithm works properly.
                            laccounted[t][p]++;
                        }
                    }
                }
                // actuator line model
                else if((*farm->turbineModels[t]) == "ALM")
                {
                    // number of points in the AD mesh
                    PetscInt npts_t = wt->alm.nPoints;

                    // loop over the AD mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->alm.thisPtControlled[p])
                        {
                            // set this point as accounted for. This is to check that
                            // the discrimination algorithm works properly.
                            laccounted[t][p]++;
                        }
                    }
                }
                else if((*farm->turbineModels[t]) == "AFM")
                {
                    // number of points in the AD mesh
                    if(wt->afm.thisPtControlled)
                    {
                        laccounted[t][0]++;
                    }
                }
            }
        }

        // scatter the accounted variable
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];
            PetscInt npts_t;

            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // number of points in the AD mesh
                npts_t = wt->adm.nPoints;
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // number of points in the AD mesh
                npts_t = wt->uadm.nPoints;
            }
            // actuator line model
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // number of points in the AL mesh
                npts_t = wt->alm.nPoints;
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // number of points in the AL mesh
                npts_t = 1;
            }

            // scatter 'accounted' for processor point owning test
            MPI_Allreduce(&(laccounted[t][0]), &(gaccounted[t][0]), npts_t, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            // check discrimination algorithm
            for(p=0; p<npts_t; p++)
            {
                Cmpnts point_p;

                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    point_p = wt->adm.points[p];
                }
                // uniform actuator disk model
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    point_p = wt->uadm.points[p];
                }
                // actuator line model
                else if((*farm->turbineModels[t]) == "ALM")
                {
                    point_p = wt->alm.points[p];
                }
                else if((*farm->turbineModels[t]) == "AFM")
                {
                    point_p = wt->afm.point;
                }

                // not accounted
                if(gaccounted[t][p] == 0)
                {
                    char warning[256];
                    sprintf(warning, "turbine point %ld at location (%.2f %.2f %.2f) was not accounted for\n", p, point_p.x, point_p.y, point_p.z);
                    warningInFunction("checkPointDiscriminationRotor",  warning);
                }
                // accounted twice
                else if(gaccounted[t][p] == 2)
                {
                    char warning[256];
                    sprintf(warning, "turbine point %ld at location (%.2f %.2f %.2f) was accounted twice\n", p, point_p.x, point_p.y, point_p.z);
                    warningInFunction("checkPointDiscriminationRotor",  warning);
                }
            }

            // clean memory
            std::vector<PetscInt> ().swap(laccounted[t]);
            std::vector<PetscInt> ().swap(gaccounted[t]);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkPointDiscriminationTower(farm_ *farm)
{
    if(farm->dbg)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscInt t, c, p;

        // test for points accounted for twice or none
        std::vector<std::vector<PetscInt>> laccounted(farm->size);
        std::vector<std::vector<PetscInt>> gaccounted(farm->size);

        // set to zero
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeTwr)
            {
                PetscInt npts_t = wt->twr.nPoints;

                laccounted[t].resize(npts_t);
                gaccounted[t].resize(npts_t);

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    laccounted[t][p] = 0;
                    gaccounted[t][p] = 0;
                }
            }
        }

        // loop over each wind turbine
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeTwr)
            {
                // test if this processor controls this turbine
                if(wt->twr.nControlled)
                {
                    // number of points in the AD mesh
                    PetscInt npts_t = wt->twr.nPoints;

                    // loop over the AD mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->twr.thisPtControlled[p])
                        {
                            // set this point as accounted for. This is to check that
                            // the discrimination algorithm works properly.
                            laccounted[t][p]++;
                        }
                    }
                }
            }
        }

        // scatter the accounted variable
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeTwr)
            {
                PetscInt npts_t = wt->twr.nPoints;

                // scatter 'accounted' for processor point owning test
                MPI_Allreduce(&(laccounted[t][0]), &(gaccounted[t][0]), npts_t, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                // check discrimination algorithm
                for(p=0; p<npts_t; p++)
                {
                    Cmpnts point_p = wt->twr.points[p];

                    // not accounted
                    if(gaccounted[t][p] == 0)
                    {
                        char warning[256];
                        sprintf(warning, "tower point %ld at location (%.2f %.2f %.2f) was not accounted for\n", p, point_p.x, point_p.y, point_p.z);
                        warningInFunction("checkPointDiscriminationTower",  warning);
                    }
                    // accounted twice
                    else if(gaccounted[t][p] == 2)
                    {
                        char warning[256];
                        sprintf(warning, "tower point %ld at location (%.2f %.2f %.2f) was accounted twice\n", p, point_p.x, point_p.y, point_p.z);
                        warningInFunction("checkPointDiscriminationTower",  warning);
                    }
                }

                // clean memory
                std::vector<PetscInt> ().swap(laccounted[t]);
                std::vector<PetscInt> ().swap(gaccounted[t]);
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkPointDiscriminationNacelle(farm_ *farm)
{
    if(farm->dbg)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscInt t, c;

        // test for points accounted for twice or none
        std::vector<PetscInt> laccounted(farm->size);
        std::vector<PetscInt> gaccounted(farm->size);

        // set to zero
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeNacelle)
            {
                laccounted[t] = 0;
                gaccounted[t] = 0;
            }
        }

        // loop over each wind turbine
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeNacelle)
            {
                // test if this processor controls this turbine
                if(wt->nac.nControlled)
                {
                    if(wt->nac.thisPtControlled) laccounted[t]++;
                }
            }
        }

        // scatter the accounted variable
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->includeNacelle)
            {
                // scatter 'accounted' for processor point owning test
                MPI_Allreduce(&laccounted[t], &gaccounted[t], 1, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

                // not accounted
                if(gaccounted[t] == 0)
                {
                    char warning[256];
                    sprintf(warning, "nacelle point from turbine %ld at location (%.2f %.2f %.2f) was not accounted for\n", t, wt->nac.point.x, wt->nac.point.y, wt->nac.point.z);
                    warningInFunction("checkPointDiscriminationNacelle",  warning);
                }
                // accounted twice
                else if(gaccounted[t] == 2)
                {
                    char warning[256];
                    sprintf(warning, "nacelle point from turbine %ld at location (%.2f %.2f %.2f) was accounted twice\n", t, wt->nac.point.x, wt->nac.point.y, wt->nac.point.z);
                    warningInFunction("checkPointDiscriminationNacelle",  warning);
                }
            }
        }

        // clean memory
        std::vector<PetscInt> ().swap(laccounted);
        std::vector<PetscInt> ().swap(gaccounted);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkPointDiscriminationSample(farm_ *farm)
{
    if(farm->dbg)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscInt t, c, p;

        // test for points accounted for twice or none
        std::vector<std::vector<PetscInt>> laccounted(farm->size);
        std::vector<std::vector<PetscInt>> gaccounted(farm->size);

        // set to zero
        for(t=0; t<farm->size; t++)
        {
            upSampling *upPoints = farm->wt[t]->upPoints;

            PetscInt npts_t = upPoints->nPoints;

            laccounted[t].resize(npts_t);
            gaccounted[t].resize(npts_t);

            // loop over the AD mesh points
            for(p=0; p<npts_t; p++)
            {
                laccounted[t][p] = 0;
                gaccounted[t][p] = 0;
            }
        }

        // loop over each wind turbine
        for(t=0; t<farm->size; t++)
        {
            upSampling *upPoints = farm->wt[t]->upPoints;

            // test if this processor controls this rig
            if(upPoints->thisRigControlled)
            {
                // number of points in the sample mesh
                PetscInt npts_t = upPoints->nPoints;

                // loop over the sample mesh points
                for(p=0; p<npts_t; p++)
                {
                    if(upPoints->thisPtControlled[p])
                    {
                        // set this point as accounted for. This is to check that
                        // the discrimination algorithm works properly.
                        laccounted[t][p]++;
                    }
                }
            }
        }

        // scatter the accounted variable
        for(t=0; t<farm->size; t++)
        {
            upSampling *upPoints = farm->wt[t]->upPoints;

            PetscInt npts_t = upPoints->nPoints;

            // scatter 'accounted' for processor point owning test
            MPI_Allreduce(&(laccounted[t][0]), &(gaccounted[t][0]), npts_t, MPIU_INT, MPI_SUM, mesh->MESH_COMM);

            // check discrimination algorithm
            for(p=0; p<npts_t; p++)
            {
                Cmpnts point_p = upPoints->points[p];

                // not accounted
                if(gaccounted[t][p] == 0)
                {
                    char warning[256];
                    sprintf(warning, "sample point %ld at location (%.2f %.2f %.2f) was not accounted for\n", p, point_p.x, point_p.y, point_p.z);
                    warningInFunction("checkPointDiscriminationSample",  warning);
                }
                // accounted twice
                else if(gaccounted[t][p] == 2)
                {
                    char warning[256];
                    sprintf(warning, "sample point %ld at location (%.2f %.2f %.2f) was accounted twice\n", p, point_p.x, point_p.y, point_p.z);
                    warningInFunction("checkPointDiscriminationSample",  warning);
                }
            }

            // clean memory
            std::vector<PetscInt> ().swap(laccounted[t]);
            std::vector<PetscInt> ().swap(gaccounted[t]);
        }
    }

    return(0);
}

