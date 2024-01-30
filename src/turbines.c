//! \file  turbines.c
//! \brief Contains top to bottom level routines for the wind farm modeling.

#include "include/base.h"
#include "include/domain.h"
#include "include/io.h"
#include "include/inline.h"
#include "include/turbines.h"

// Our ADM/ALM models are detailed actuator disk/line
// models with blade pitch and generator speed controls, blade and airfoil properties
// linearly interpolated where necessary. UniformADM is a simple AD version, where the
// thrust is weighted based on the rotor mesh element area. Yaw control
// can be activated for all models and is based upon 1D upstream velocity sampling. Turbines
// rotate with a prescribed nacelle rotation speed.
// We have fully optimized parallel communication: each turbine, tower and up-sampling
// rig (it s just a rotor defined 2.5D upstream of each rotor for velocty sampling)
// has its own communicator, and turbine solution is done in parallel by each communicator.
// Each turbine has an influence sphere of mesh cells, where the rotor can possibly
// be located in any yaw condition, and those cells define the turbine communicator based
// on the processor which they belong to. In this way communications are performed between
// max 8 processors (that is the maximum number of processors that can control a wind
// turbine, most of the time is 2 or 1). Turbine I/O is performed also in parallel:
// the master rank of each communicator is responsible for writing down turbine log
// and checkpoint file. In normal operation (no debug), a global communication between
// all turbines is never performed. There is only 1 communication of type MPI_Reduce
// for wind turbine yaw syncronization on the master node, which will then write
// the mesh to file. That is done only at global write time (not at turbine log write time).
// The only function call for turbine update is "PetscErrorCode UpdateWindTurbines(farm_ *farm)"".


PetscErrorCode UpdateWindTurbines(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    PetscReal ts, te;

    PetscTime(&ts);

    if(!farm->dbg) PetscPrintf(mesh->MESH_COMM,"Updating wind turbines, ");

    // solve rotor dynamics and compute rot speeds
    computeRotSpeed(farm);

    // apply the generator speed controller
    controlGenSpeed(farm);

    // apply the pitch PID controller
    controlBldPitch(farm);

    // apply the wind farm controller
    windFarmControl(farm);

    // apply the yaw controller
    controlNacYaw(farm);

    // find which AD points this processor controls
    findControlledPointsRotor(farm);

    // find which sample points this processor controls
    findControlledPointsSample(farm);

    // compute wind velocity at the rotor mesh points
    computeWindVectorsRotor(farm);

    // compute wind velocity at the tower mesh points
    computeWindVectorsTower(farm);

    // compute wind velocity at the nacelle mesh point
    computeWindVectorsNacelle(farm);

    // compute wind velocity at the sample mesh points
    computeWindVectorsSample(farm);

    // compute aerodynamic forces at the turbine mesh points
    computeBladeForce(farm);

    // project the wind turbine forces on the background mesh
    projectBladeForce(farm);

    // project the tower forces on the background mesh
    projectTowerForce(farm);

    // project the nacelle forces on the background mesh
    projectNacelleForce(farm);

    // transform forces from cartesian to contravariant
    bodyForceCartesian2Contravariant(farm);

    // write output if runTimeWrite == 1
    windTurbinesWrite(farm);

    PetscTime(&te);

    if(!farm->dbg) PetscPrintf(mesh->MESH_COMM,"Elapsed Time = %lf\n", te-ts);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode InitializeWindFarm(farm_ *farm)
{
    if(farm != NULL)
    {
        mesh_ *mesh = farm->access->mesh;

        PetscPrintf(mesh->MESH_COMM, "Initializing wind turbines...\n");

        // read farm propeties and initialize all parameters
        PetscPrintf(mesh->MESH_COMM, "   reading wind farm properties...\n");
        readFarmProperties(farm);

        PetscPrintf(mesh->MESH_COMM, "   initializing wind farm models...\n");
        // initialize the wind turbine models
        PetscInt nT = farm->size;

        for(PetscInt t=0; t<nT; t++)
        {
            if((*farm->turbineModels[t]) == "ADM")
            {
                // initialize actuator disk model parameters
                initADM(farm->wt[t], farm->base[t], mesh->meshName);
            }
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // initialize uniform actuator disk model parameters
                initUADM(farm->wt[t], farm->base[t], mesh->meshName);
            }
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // initialize actuator line model parameters
                initALM(farm->wt[t], farm->base[t], mesh->meshName);
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // initialize actuator farm model parameters
                initAFM(farm->wt[t], farm->base[t], mesh->meshName);
            }
            else
            {
               char error[512];
               sprintf(error, "unknown wind turbine model for turbine %s\n", (*farm->turbineIds[t]).c_str());
               fatalErrorInFunction("InitializeWindFarm",  error);
            }

            // initialize the tower model
            if(farm->wt[t]->includeTwr)
            {
                initTwrModel(farm->wt[t], farm->base[t]);
            }

            // initialize the tower model
            if(farm->wt[t]->includeNacelle)
            {
                initNacModel(farm->wt[t], farm->base[t]);
            }

            // initialize upstream velocity sampling
            initSamplePoints(farm->wt[t], farm->base[t], mesh->meshName);
        }

        PetscPrintf(mesh->MESH_COMM, "   pre-calculating influence sphere cells...\n");

        // initialize turbine and tower sphere cells and see what processor controls which wind turbine
        initControlledCells(farm);

        // ensure best practices with respect to mesh size and turbine dimension are followed
        checkTurbineMesh(farm);

        // initialize sample sphere cells and see what processor controls which sampling rig
        initSampleControlledCells(farm);

        // determine which points in each tower are controlled by which processor (also finds closest cell)
        findControlledPointsTower(farm);

        // determine which points in each nacelle are controlled by which processor (also finds closest cell)
        findControlledPointsNacelle(farm);

        // read the checkpoint file if present and prepare wind turbines for restart
        windTurbinesReadCheckpoint(farm);

        PetscPrintf(mesh->MESH_COMM, "done\n\n");
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode checkTurbineMesh(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, p, c;

    PetscMPIInt      nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt      rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    Cmpnts           ***cent;
    PetscReal        ***aj;
    PetscReal        lMaxCell, gMaxCell, lMinCell, gMinCell;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    // max perturbation amplitude
    PetscReal maxPerturb  = 1e-10;

    // processor perturbation (changes between processors)
    PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        Cmpnts rtrAxis  = farm->wt[t]->rtrAxis;

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            lMaxCell = 0.0, gMaxCell = 0.0;
            lMinCell = 1e20, gMinCell = 1e20;

            // loop in sphere points
            for(c=0; c<wt->nControlled; c++)
            {
                // cell indices
                PetscInt i = wt->controlledCells[c].i,
                    j = wt->controlledCells[c].j,
                    k = wt->controlledCells[c].k;

                if(lMaxCell < 1/aj[k][j][i])
                {
                    lMaxCell = 1/aj[k][j][i];
                }

                if(lMinCell > 1/aj[k][j][i])
                {
                    lMinCell = 1/aj[k][j][i];
                }

            }

            MPI_Allreduce(&lMaxCell, &gMaxCell, 1, MPIU_REAL, MPIU_MAX, wt->TRB_COMM);
            MPI_Allreduce(&lMinCell, &gMinCell, 1, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

            gMaxCell = std::pow(gMaxCell, 1.0/3.0);
            gMinCell = std::pow(gMinCell, 1.0/3.0);


            if(2.0 * wt->rTip < 8.0 * gMaxCell)
            {
                char warning[512];
                sprintf(warning, "turbine diameter (%lf) < 10 mesh cells (%lf), not resolved properly. Revise mesh to improve resolution.\n", 2.0 * wt->rTip, 10.0 * gMaxCell);
                warningInFunction("checkTurbineMesh",  warning);
            }

            //additional checks for Advanced actuator line method
            if(((*farm->turbineModels[t]) == "ALM"))
            {
                if(wt->alm.projectionType == "anisotropic")
                {
                    //check that there are enough radial elements
                    PetscReal cellAvg;

                    cellAvg = 0.5 * (gMaxCell + gMinCell);

                    PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                    PetscReal drvalOpt = 2.0 *  cellAvg;

                    PetscInt nRadial = PetscInt((wt->rTip - wt->rHub)/drvalOpt) + 1;

                    if(drval > drvalOpt)
                    {
                        char error[512];
                        sprintf(error, "Not enough radial elements in the Actuator line mesh. Increase the radial resolution to have %ld points\n", nRadial);
                        fatalErrorInFunction("checkTurbineMesh",  error);
                    }

                    if(drval < gMinCell)
                    {
                        char error[512];
                        sprintf(error, "Too many radial elements. Radial element size smaller than mesh size. Decrease the radial resolution to have %ld points\n", nRadial);
                        fatalErrorInFunction("checkTurbineMesh",  error);
                    }

                    //check the mesh size with respect to the projection radius
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

                    MPI_Allreduce(&(lminDist[0]), &(gminDist[0]), wt->alm.nPoints, MPIU_REAL, MPIU_MIN, wt->TRB_COMM);

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

                    // loop over the AL mesh points
                    for(p=0; p<npts_t; p++)
                    {
                        if(wt->alm.thisPtControlled[p])
                        {
                            // get the closest cell center
                            PetscInt i = wt->alm.closestCells[p].i,
                                     j = wt->alm.closestCells[p].j,
                                     k = wt->alm.closestCells[p].k;

                            PetscReal cellsize = std::pow(1/aj[k][j][i], 1.0/3.0);

                            PetscInt nRadial = PetscInt(0.8 * wt->alm.nRadial);

                            PetscInt radPt = PetscInt (p /wt->alm.nAzimuth);

                            PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                            PetscReal eps_x =     wt->alm.chord[p]*wt->eps_x,
                                      eps_y =     wt->alm.thick[p]*wt->eps_y,
                                      eps_z =     drval * wt->eps_z;

                            PetscPrintf(PETSC_COMM_SELF, "radial point = %ld, nRadial = %ld, eps_x = %lf \n", radPt, nRadial, eps_x);
                            // excluding the tip points for this check
                            if(radPt <= nRadial)
                            {
                                if(eps_x < 1.2 * cellsize)
                                {
                                    char error[512];
                                    sprintf(error, "Fluid Mesh size not optimal for AALM simulation. At radial distance %lf m (%0.2lf %%) from blade root, %lf * chordLength (%lf) < 1.5 * cell size (%lf). Refine mesh or increase epsilonFactor_x (Note: epsilonFactor_x optimal <= 1)\n", drval * radPt, drval * radPt * 100/(wt->rTip - wt->rHub), wt->eps_x, wt->alm.chord[p], cellsize);
                                    fatalErrorInFunction("checkTurbineMesh",  error);
                                }
                            }
                        }
                    }
                }
            }

        }
    }

    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    return(0);
}
//***************************************************************************************************************//
PetscErrorCode computeRotSpeed(farm_ *farm)
{
    // turbine and model mesh point indices
    PetscInt t;

    // set physical time clock pointer
    clock_ *clock = farm->access->clock;

    // set io pointer
    io_    *io    = farm->access->io;

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
        {
            // test if this processor controls this turbine
            if(wt->turbineControlled)
            {
                // do not distinguish between ADM or ALM, this is a
                // global wind turbine property
                if(wt->genControllerType == "none")
                {
                    // set to zero, controller deactivated
                    wt->genOmega = 0.0;
                    wt->genPwr   = 0.0;
                }
                else
                {
                    PetscReal rtrTorque;

                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        rtrTorque = wt->adm.rtrTorque;
                    }
                    else if((*farm->turbineModels[t]) == "ALM")
                    {
                        rtrTorque = wt->alm.rtrTorque;
                    }

                    // solve rotor dynamics and update rotor speed
                    wt->rtrOmega
                    +=
                    (clock->dt / wt->driveTrainInertia) *
                    (
                        wt->gbxEff * rtrTorque -
                        wt->gbxRatioG2R * wt->genTorque
                    );

                    // limit rotational speed if applicable
                    if(wt->rtrSpdLimiter)
                    {
                        wt->rtrOmega
                        =
                        PetscMin
                        (
                            PetscMax(0.0, wt->rtrOmega),
                            wt->ratedRotorSpd
                        );
                    }

                    // simulate acquisition system of turbine velocity and
                    // signal filtering: 1-pole low pass recursive filter
                    PetscReal a = std::exp(-clock->dt * wt->rtrSpdFilterFreq);

                    wt->rtrOmegaFilt = (1-a)*wt->rtrOmega + a*wt->rtrOmegaFilt;

                    // compute the generator speed (rpm) with the filtered rotor speed
                    wt->genOmega = wt->rtrOmegaFilt * wt->gbxRatioG2R;
                }

                // rotate points if turbine model is ALM
                if((*farm->turbineModels[t]) == "ALM")
                {
                    PetscReal angle = clock->dt * wt->rtrOmega;
                    rotateBlades(wt, angle);
                }
            }
            else
            {
                // zero the rotor and generator speeds in those processors which do not
                // control this wind turbine (in order obtain correct values when summing
                // among processors and dividing by the number of controlling processors)
                wt->genPwr       = 0.0;
                wt->rtrOmega     = 0.0;
                wt->genOmega     = 0.0;
                wt->rtrOmegaFilt = 0.0;

                if((*farm->turbineModels[t]) == "ALM")
                {
                    wt->alm.azimuth = 0.0;
                }
            }

            // sync rotor azimuts on the master rank for writing (only for ALM)
            if(io->runTimeWrite && (*farm->turbineModels[t]) == "ALM")
            {
                // set azimut to zero before scattering
                PetscReal gazimuth;

                // set mesh pointer
                mesh_ *mesh = farm->access->mesh;

                PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

                // scatter on the master rank
                MPI_Reduce(&(wt->alm.azimuth), &gazimuth, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

                // divide by nprocs in turbine communicator
                gazimuth = gazimuth / wt->nProcsTrb;

                // rotate turbines
                if(!rank)
                {
                    PetscReal angle = (gazimuth - wt->alm.azimuth) * wt->deg2rad;
                    rotateBlades(wt, angle);
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode rotateBlades(windTurbine *wt, PetscReal angle)
{
    // number of points in the AL mesh
    PetscInt p, npts_t = wt->alm.nPoints;

    // loop over the AL mesh points
    for(p=0; p<npts_t; p++)
    {
        // save this point locally for speed
        Cmpnts point_p = wt->alm.points[p];

        // this point position from COR
        mSub(wt->alm.points[p], wt->rotCenter);

        // rotate points in COR ref frame
        mRot(wt->omega_hat, wt->alm.points[p], angle);

        // add COR ref frame
        mSum(wt->alm.points[p], wt->rotCenter);
    }

    wt->alm.azimuth += angle * wt->rad2deg;

    // bound azimuth between 360 and 0
    if(wt->alm.azimuth >= 360.0)
    {
        wt->alm.azimuth -= 360.0;
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode controlGenSpeed(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    // turbine and AD mesh point indices
    PetscInt t;

    // set physical time clock pointer
    clock_ *clock = farm->access->clock;

    PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
        {
            // test if this processor controls this turbine
            if(wt->turbineControlled)
            {
                // do not distinguish between ADM or ALM, this is a
                // global wind turbine property
                if(wt->genControllerType == "none")
                {
                    // do nothing
                }
                else
                {
                    PetscReal commGenTq = wt->genTorque;

                    // region 1: disconnected
                    if(wt->genOmega < wt->cutInGenSpd)
                    {
                        commGenTq = 0.0;
                    }
                    // region 1 and 1/2: linearly interpolate
                    else if
                    (
                        wt->genOmega >= wt->cutInGenSpd &&
                        wt->genOmega < wt->regTwoStartGenSpd
                    )
                    {
                        PetscReal deltaGenOmega    = wt->genOmega - wt->cutInGenSpd;
                        PetscReal regTwoStartGenTq = wt->omegaKP * wt->regTwoStartGenSpd * wt->regTwoStartGenSpd;
                        PetscReal torqueSlope      = (regTwoStartGenTq - wt->cutInGenTq) / (wt->regTwoStartGenSpd - wt->cutInGenSpd);
                        commGenTq = wt->cutInGenTq + torqueSlope * deltaGenOmega;
                    }
                    // region 2: proportional controller
                    else if
                    (
                        wt->genOmega >= wt->regTwoStartGenSpd &&
                        wt->genOmega < wt->regTwoEndGenSpd
                    )
                    {
                        commGenTq = wt->omegaKP * wt->genOmega * wt->genOmega;
                    }
                    // region 2 and 1/2: linearly extrapolate
                    else if
                    (
                        wt->genOmega >= wt->regTwoEndGenSpd &&
                        wt->genOmega < (wt->ratedRotorSpd*wt->gbxRatioG2R)
                    )
                    {
                        PetscReal ratedGenSpd      = wt->ratedRotorSpd*wt->gbxRatioG2R;
                        PetscReal deltaGenOmega    = wt->genOmega - wt->regTwoEndGenSpd;
                        PetscReal regTwoEndGenTq   = wt->omegaKP * wt->regTwoEndGenSpd * wt->regTwoEndGenSpd;
                        PetscReal torqueSlope      = (wt->ratedGenTq - regTwoEndGenTq) / (ratedGenSpd - wt->regTwoEndGenSpd);
                        commGenTq = regTwoEndGenTq + torqueSlope * deltaGenOmega;
                    }
                    // region 3: fix the rotor torque
                    else if(wt->genOmega >= (wt->ratedRotorSpd*wt->gbxRatioG2R))
                    {
                        commGenTq = wt->ratedGenTq;
                    }

                    // apply torque rate limiter
                    if(wt->tqRateLimiter)
                    {
                        PetscReal tqRate = (commGenTq - wt->genTorque) / PetscMax(clock->dt, 1e-5);
                        tqRate    = PetscMin(PetscMax(tqRate, -wt->tqMaxRate), wt->tqMaxRate);
                        commGenTq = wt->genTorque + (tqRate * clock->dt);
                    }

                    if(wt->dbg && rank == wt->writerRank)
                    {
                        printf("Turbine %s: commGenTq = %.3f KNm, commGenTqRate = %.3f KNm/s, rotorSpeed = %.3f rpm\n", (*farm->turbineIds[t]).c_str(), commGenTq / 1000, (commGenTq - wt->genTorque) / 1000 / clock->dt, wt->rtrOmega/wt->rpm2RadSec);
                    }

                    // update generator torque
                    wt->genTorque = commGenTq;

                }
            }
            else
            {
                wt->genTorque = 0.0;
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode controlBldPitch(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    // turbine and AD mesh point indices
    PetscInt t;

    // set physical time clock pointer
    clock_ *clock = farm->access->clock;

    PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
        {
            // test if this processor controls this turbine
            if(wt->turbineControlled)
            {
                // do not distinguish between ADM or ALM, this is a
                // global wind turbine property
                if(wt->pitchControllerType == "none")
                {
                    // set collective pitch to zero
                    wt->collPitch = 0.0;
                }
                else
                {
                    // compute gain scheduling
                    PetscReal G = 1.0 / (1.0 + wt->collPitch / wt->pitchS2R);

                    // save old PID error
                    PetscReal errPID_old = wt->errPID;

                    // activation above rated speed is ensured by pitch saturation
                    wt->errPID = wt->rtrOmega - wt->ratedRotorSpd;

                    // cumulate integral PID error
                    wt->intErrPID    += wt->errPID * clock->dt;

                    // compute derivative PID error
                    PetscReal derErrPID  = (wt->errPID - errPID_old) / PetscMax(clock->dt, 1e-5);

                    // saturate the integrated PID error based on pitch min/max
                    PetscReal intMinErr  = wt->pitchMin / (G * wt->pitchKI);
                    PetscReal intMaxErr  = wt->pitchMax / (G * wt->pitchKI);
                    wt->intErrPID     = PetscMin(PetscMax(wt->intErrPID, intMinErr), intMaxErr);

                    // compute PID contributions
                    PetscReal pitchP     = G * wt->pitchKP * wt->errPID;
                    PetscReal pitchI     = G * wt->pitchKI * wt->intErrPID;
                    PetscReal pitchD     = G * wt->pitchKD * derErrPID;

                    // compute commanded pitch
                    PetscReal pitchComm  = pitchP + pitchI + pitchD;

                    // saturate pitch based on pitch min/max
                    pitchComm = PetscMin(PetscMax(pitchComm, wt->pitchMin), wt->pitchMax);

                    // add single blade components for specific control systems
                    // (not implemented)

                    if(wt->dbg && rank == wt->writerRank)
                    {
                        printf("Turbine %s: commPitch = %.3f deg, commPitchRate = %.3f deg/s, rotorSpeed = %.3f rpm\n", (*farm->turbineIds[t]).c_str(), wt->collPitch*wt->rad2deg, (pitchComm-wt->collPitch)*wt->rad2deg/clock->dt, wt->rtrOmega/wt->rpm2RadSec);
                    }

                    wt->collPitch = pitchComm;
                }
            }
            else
            {
                wt->collPitch = 0.0;
                wt->intErrPID = 0.0;
                wt->errPID    = 0.0;
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode controlNacYaw(farm_ *farm)
{
    // Seba: each turbine is yawed only if this processor actually controls it. This means that the same turbine
    //       has different yaw angles depending on the processor, and only the controlling processor will have
    //       the true yaw angle. They are put back in sync only on the master node, which is the one who writes
    //       the wind turbine .inp file, but only when runTimeWrite is one (when the actual file must be written).
    //       This strongly reduces the parallel communications only on those processors which control the turbines.

    mesh_ *mesh = farm->access->mesh;

    PetscInt t, c, p, nYawed = 0;

    // get global parallel info
    PetscMPIInt nprocs; MPI_Comm_size(mesh->MESH_COMM, &nprocs);
    PetscMPIInt rank;   MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // count how many turbines have the yaw controller and avoid the search a priori
    // in case no turbine has the yaw controller
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        if(wt->yawControllerType != "none")
        {
            nYawed++;
        }
    }

    // yaw angle
    PetscReal angle;

    if(nYawed)
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

        Cmpnts           ***cent, ***lucat;

        DMDAVecGetArray(fda, mesh->lCent, &cent);
        DMDAVecGetArray(fda, ueqn->lUcat, &lucat);

        // max perturbation amplitude
        PetscReal maxPerturb  = 1e-10;

        // processor perturbation for search (changes between processors)
        PetscReal procContrib = maxPerturb * ((PetscReal)rank + 1) / (PetscReal)nprocs;

        Cmpnts perturbVec;
               perturbVec.x = procContrib;
               perturbVec.y = procContrib;
               perturbVec.z = procContrib;

        // update sampling points
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            upSampling* upPoints = farm->wt[t]->upPoints;

            // test if this turbine has yaw control
            if(wt->yawControllerType != "none")
            {
                // test if this processor controls this turbine
                if(upPoints->thisRigControlled)
                {
                    // hubUpDistance sampling type
                    if(wt->yawSamplingType == "hubUpDistance")
                    {
                        // local and global min distances
                        PetscReal  lminDist = 1e20;
                        PetscReal  gminDist = 1e20;
                        cellIds lclosestIds; lclosestIds.i = lclosestIds.j = lclosestIds.k = 0;

                        // upstream distance is of 1 diameter
                        Cmpnts upDist  = nScale(2.0*wt->rTip, wt->rtrDir);

                        // this turbine rotor center
                        Cmpnts rotor_c = farm->wt[t]->rotCenter;

                        // this turbine upstream point
                        Cmpnts upPoint = nSum(rotor_c, upDist);

                        // loop in sphere points
                        for(c=0; c<upPoints->nControlled; c++)
                        {
                            // cell indices
                            PetscInt i = upPoints->controlledCells[c].i,
                                     j = upPoints->controlledCells[c].j,
                                     k = upPoints->controlledCells[c].k;

                            Cmpnts distVec = nSub(upPoint, cent[k][j][i]);
                                             mSub(distVec, perturbVec);
                            PetscReal distMag = nMag(distVec);

                            if(distMag < lminDist)
                            {
                                lminDist      = distMag;
                                lclosestIds.i = i;
                                lclosestIds.j = j;
                                lclosestIds.k = k;
                            }
                        }

                        // scatter the minimum dist to all processors
                        MPI_Allreduce(&lminDist, &gminDist, 1, MPIU_REAL, MPIU_MIN, upPoints->UPW_COMM);

                        PetscReal  lyawAngle  = 0.0;
                        PetscReal  lflowAngle = 0.0;

                        if(gminDist == lminDist)
                        {
                            wt->yawSampleIds.i = lclosestIds.i;
                            wt->yawSampleIds.j = lclosestIds.j;
                            wt->yawSampleIds.k = lclosestIds.k;

                            // sample the velocity upstream
                            Cmpnts upU    = lucat[wt->yawSampleIds.k][wt->yawSampleIds.j][wt->yawSampleIds.i];

                            // compute flow angle
                            lflowAngle = wt->flowAngle * std::exp(-clock->dt/wt->yawAverageWindow) + clock->dt / wt->yawAverageWindow * std::atan2(upU.y, upU.x) * wt->rad2deg;

                            // compute turbine yaw angle
                            lyawAngle  = std::atan2(-1.0 * wt->rtrDir.y, -1.0 * wt->rtrDir.x) * wt->rad2deg;
                        }

                        // scatter info to all processors from the rig comm: if the point is oustude of TRB_COMM but inside UPW_COMM
                        // scattering on TRB_COMM would result in zero flow and yaw angles. We have to use UPW_COMM
                        MPI_Allreduce(&lyawAngle,  &(wt->yawAngle),  1, MPIU_REAL, MPIU_SUM, upPoints->UPW_COMM);
                        MPI_Allreduce(&lflowAngle, &(wt->flowAngle), 1, MPIU_REAL, MPIU_SUM, upPoints->UPW_COMM);

                        wt->yawError  = wt->flowAngle - wt->yawAngle;

                        PetscReal rotationSign;

                        if(wt->yawError > 1.0*wt->yawAllowedError)       rotationSign =  1.0;
                        else if(wt->yawError < -1.0*wt->yawAllowedError) rotationSign = -1.0;
                        else                                             rotationSign =  0.0;

                        // saturate if yaw angle is out of bounds
                        if(wt->yawAngle >= wt->yawMax || wt->yawAngle <= wt->yawMin)
                        {
                            rotationSign = 0.0;
                        }

                        // compute rotation angle in radiants
                        angle = rotationSign * clock->dt * wt->yawSpeed * wt->deg2rad;

                        // rotate turbine
                        if(angle != 0.0)
                        {
                            // rotate common parameters
                            mRot(wt->twrDir,  wt->rtrDir, angle);
                            mRot(wt->twrDir,  wt->rtrAxis, angle);

                            // actuator disk model
                            if((*farm->turbineModels[t]) == "ADM")
                            {
                                // number of points in the AD mesh
                                PetscInt npts_t = wt->adm.nPoints;

                                // loop over the AD mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(wt->adm.points[p], wt->rotCenter);
                                    mRot(wt->twrDir, wt->adm.points[p], angle);
                                    mSum(wt->adm.points[p], wt->rotCenter);
                                }

                                // rotate other parameters
                                mRot(wt->twrDir, wt->omega_hat, angle);
                            }
                            // uniform actuator disk model
                            else if((*farm->turbineModels[t]) == "uniformADM")
                            {
                                // number of points in the AD mesh
                                PetscInt npts_t = wt->uadm.nPoints;

                                // loop over the AD mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(wt->uadm.points[p], wt->rotCenter);
                                    mRot(wt->twrDir, wt->uadm.points[p], angle);
                                    mSum(wt->uadm.points[p], wt->rotCenter);
                                }
                            }
                            // actuator line model
                            else if((*farm->turbineModels[t]) == "ALM")
                            {
                                // number of points in the AL mesh
                                PetscInt npts_t = wt->alm.nPoints;

                                // loop over the AL mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(wt->alm.points[p], wt->rotCenter);
                                    mRot(wt->twrDir, wt->alm.points[p], angle);
                                    mSum(wt->alm.points[p], wt->rotCenter);
                                }

                                // rotate other parameters
                                mRot(wt->twrDir, wt->omega_hat, angle);
                            }
                            else if((*farm->turbineModels[t]) == "AFM")
                            {
                                // nothing to do
                            }

                            // rotate up-sampling points
                            {
                                upSampling *upPoints = farm->wt[t]->upPoints;

                                // rotate center
                                mSub(upPoints->center, wt->rotCenter);
                                mRot(wt->twrDir, upPoints->center, angle);
                                mSum(upPoints->center, wt->rotCenter);

                                // number of points in the sample mesh
                                PetscInt npts_t = upPoints->nPoints;

                                // loop over the sample mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(upPoints->points[p], wt->rotCenter);
                                    mRot(wt->twrDir, upPoints->points[p], angle);
                                    mSum(upPoints->points[p], wt->rotCenter);
                                }
                            }
                        }

                        // turbine point search will be performed
                        wt->yawChanged = 1;
                    }
                    else if(wt->yawSamplingType == "anemometer")
                    {
                        // not implemented
                    }
                }
            }
        }

        DMDAVecRestoreArray(fda, mesh->lCent, &cent);
        DMDAVecRestoreArray(fda, ueqn->lUcat, &lucat);
    }

    // for each processor, the wind turbines are yawed only if they are controlled, so that the
    // solution can proceed smoothly. This causes problems when the .inp file is written. Since it
    // is written by the master rank of the mesh->MESH_COMM communicator, we have to sync the
    // turbines there prior to write, so that the mesh will reflect the actual yaw of each turbine

    io_ *io = farm->access->io;

    if(io->runTimeWrite)
    {
        // update sampling points
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // test if this turbine has yaw control
            if(wt->yawControllerType != "none")
            {
                // true yaw value that the turbine should have
                PetscReal ltrueYaw = 0.0, gtrueYaw = 0.0;

                if (wt->turbineControlled) ltrueYaw = std::atan2(-1.0 * wt->rtrDir.y, -1.0 * wt->rtrDir.x) * wt->rad2deg;

                MPI_Reduce(&ltrueYaw, &gtrueYaw, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

                gtrueYaw = gtrueYaw / wt->nProcsTrb;

                // current yaw value that the turbine has
                PetscReal currYaw = std::atan2(-1.0 * wt->rtrDir.y, -1.0 * wt->rtrDir.x) * wt->rad2deg;

                if(wt->dbg)
                {
                    PetscPrintf(mesh->MESH_COMM, " > turbine %s: true yaw = %lf, current yaw = %lf, nProcs = %ld\n",(*farm->turbineIds[t]).c_str(), gtrueYaw, currYaw, (PetscInt)wt->nProcsTrb);
                }

                // compute rotation angle
                PetscReal yawError = gtrueYaw - currYaw;

                // recast rotation angle in [0 - 360 interval]
                if(yawError < 0.0) yawError = 360 + yawError;

                // compute rotation angle in radiants (rotation sign is always positive)
                angle = yawError * wt->deg2rad;

                // perform the rotation
                if(angle != 0.0 && !rank)
                {
                    // rotate common parameters
                    mRot(wt->twrDir,  wt->rtrDir, angle);
                    mRot(wt->twrDir,  wt->rtrAxis, angle);

                    // actuator disk model
                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->adm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->adm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->adm.points[p], angle);
                            mSum(wt->adm.points[p], wt->rotCenter);
                        }

                        // rotate other parameters
                        mRot(wt->twrDir, wt->omega_hat, angle);
                    }
                    // uniform actuator disk model
                    else if((*farm->turbineModels[t]) == "uniformADM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->uadm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->uadm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->uadm.points[p], angle);
                            mSum(wt->uadm.points[p], wt->rotCenter);
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
                            mSub(wt->alm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->alm.points[p], angle);
                            mSum(wt->alm.points[p], wt->rotCenter);
                        }

                        // rotate other parameters
                        mRot(wt->twrDir, wt->omega_hat, angle);
                    }
                    else if((*farm->turbineModels[t]) == "AFM")
                    {
                        // nothing to do
                    }

                    // rotate up-sampling points
                    {
                        upSampling *upPoints = farm->wt[t]->upPoints;

                        // rotate center
                        mSub(upPoints->center, wt->rotCenter);
                        mRot(wt->twrDir, upPoints->center, angle);
                        mSum(upPoints->center, wt->rotCenter);

                        // number of points in the sample mesh
                        PetscInt npts_t = upPoints->nPoints;

                        // loop over the sample mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(upPoints->points[p], wt->rotCenter);
                            mRot(wt->twrDir, upPoints->points[p], angle);
                            mSum(upPoints->points[p], wt->rotCenter);
                        }
                    }
                }
            }
        }
    }

    MPI_Barrier(mesh->MESH_COMM);
    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windFarmControl(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    // turbine and AD mesh point indices
    PetscInt t;

    // set physical time clock pointer
    clock_ *clock = farm->access->clock;

    PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            if(farm->farmControlActive[t])
            {
                // get the 2 time values closest to the current time
                PetscInt idx_1 = wt->currentCloseIdx;
                PetscInt idx_2 = wt->currentCloseIdx + 1;

                PetscInt lwrBound = 0;
                PetscInt uprBound = wt->wfControlNData;

                // if past first iteration do the search on a subset to speed up the process
                if(clock->it > clock->itStart)
                {
                    lwrBound = PetscMax(0, (wt->currentCloseIdx - 50));
                    uprBound = PetscMin(wt->wfControlNData, (wt->currentCloseIdx + 50));
                }

                // build error vector for the time search
                PetscReal  diff[wt->wfControlNData];

                for(PetscInt i=lwrBound; i<uprBound; i++)
                {
                    diff[i] = fabs(wt->wfControlTimes[i] - clock->time);
                }

                // find the two closest times
                for(PetscInt i=lwrBound; i<uprBound; i++)
                {
                    if(diff[i] < diff[idx_1])
                    {
                        idx_2 = idx_1;
                        idx_1 = i;
                    }
                    if(diff[i] < diff[idx_2] && i != idx_1)
                    {
                        idx_2 = i;
                    }
                }

                // always put the lower time at idx_1 and higher at idx_2
                if(wt->wfControlTimes[idx_2] < wt->wfControlTimes[idx_1])
                {
                    PetscInt idx_tmp = idx_2;
                    idx_2 = idx_1;
                    idx_1 = idx_tmp;
                }

                // find interpolation weights
                PetscReal idx = (idx_2 - idx_1) / (wt->wfControlTimes[idx_2] - wt->wfControlTimes[idx_1]) * (clock->time - wt->wfControlTimes[idx_1]) + idx_1;
                PetscReal w1 = (idx_2 - idx) / (idx_2 - idx_1);
                PetscReal w2 = (idx - idx_1) / (idx_2 - idx_1);

                // reset the closest index for nex iteration
                wt->currentCloseIdx = idx_1;

                // these models are controlled by changing the Ct coefficient
                if((*farm->turbineModels[t]) == "uniformADM" || (*farm->turbineModels[t]) == "AFM")
                {
                    // save commanded Ct delta
                    wt->wfControlCt = w1* wt->wfControlValues[idx_1] + w2 * wt->wfControlValues[idx_2];
                }
                // these models are controlled by changing the blade pitch
                else if ((*farm->turbineModels[t]) == "ADM" || (*farm->turbineModels[t]) == "ALM")
                {
                    // save commanded pitch delta
                    wt->wfControlCollPitch = w1* wt->wfControlValues[idx_1] + w2 * wt->wfControlValues[idx_2];
                }
            }
        }
    }

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
            // actuator line model
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
        if(wt->turbineControlled)
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

                    PetscReal dA   = wt->adm.dr[p] * wt->adm.chord[p];

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

                    PetscReal dA   = wt->alm.dr[p] * wt->alm.chord[p];

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
                    // projection distance
                    PetscReal eps_x  = wt->eps_x,
                              eps_y  = wt->eps_y,
                              eps_z  = wt->eps_z;

                    // save this point locally for speed
                    Cmpnts point_p    = wt->afm.point;

                    // loop in sphere points
                    for(c=0; c<wt->nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->controlledCells[c].i,
                                 j = wt->controlledCells[c].j,
                                 k = wt->controlledCells[c].k;

                        // compute distance from mesh cell to AF point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute projection factor
                        PetscReal pf
                        =
                        std::exp
                        (
                            -r_c.x * r_c.x / (eps_x * eps_x)
                            -r_c.y * r_c.y / (eps_y * eps_y)
                            -r_c.z * r_c.z / (eps_z * eps_z)
                        ) /
                        (
                            eps_x * eps_y * eps_z *
                            pow(M_PI, 1.5)
                        );

                        mSum(lU, nScale(pf/aj[k][j][i], ucat[k][j][i]));
                    }

                    MPI_Allreduce(&lU, &gU, 3, MPIU_REAL, MPIU_SUM, wt->TRB_COMM);

                    wt->afm.U = nSet(gU);

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
        if(wt->includeTwr)
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
                        PetscInt i = wt->twr.closestCells[p].i,
                                 j = wt->twr.closestCells[p].j,
                                 k = wt->twr.closestCells[p].k;

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

                        lU[p]     = nSet(uc_p);
                    }
                }

                MPI_Allreduce(&(lU[0]), &(wt->twr.U[0]), wt->twr.nPoints*3, MPIU_REAL, MPIU_SUM, wt->twr.TWR_COMM);



                // loop over the tower mesh points
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
        if(wt->includeNacelle)
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
        if(upPoints->thisRigControlled)
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

//***************************************************************************************************************//

PetscErrorCode computeBladeForce(farm_ *farm)
{
    // turbine and AD mesh point indices
    PetscInt t, p;

    // simulation constants
    constants_ *constants = farm->access->constants;

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        Cmpnts rtrAxis  = farm->wt[t]->rtrAxis;

        // test if this processor controls this turbine
        if(wt->turbineControlled)
        {
            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // apply wind farm control
                PetscReal wfControlPitch = 0.0;

                // apply wind farm controller
                if(farm->farmControlActive[t])
                {
                    wfControlPitch = wt->wfControlCollPitch;
                }

                // zero the rotor torque at this time step
                wt->adm.rtrTorque = 0.0;

                // cumulate rotor thrust at this time step
                wt->adm.rtrThrust = 0.0;

                // number of points in the AD mesh
                PetscInt npts_t = wt->adm.nPoints;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // do not test if this point is inside this processor.
                    // If at least 1 AD point is inside this processor we have
                    // to compute the parameters at all AD points.
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
                            xb_hat = nScale(-1.0, rtrAxis);
                            yb_hat = nCross(zb_hat, xb_hat);
                        }
                        else if(wt->rotDir == "ccw")
                        {
                            zb_hat = nUnit(r_p);
                                     mScale(-1.0, zb_hat);
                            xb_hat = nScale(-1.0, rtrAxis);
                            yb_hat = nCross(zb_hat, xb_hat);
                        }

                        // transform the velocity in the bladed reference frame in
                        // order to compute the angle of attack. Radial (z) component
                        // is not considered
                        PetscReal ub_x = nDot(wt->adm.U[p], xb_hat);
                        PetscReal ub_y = nDot(wt->adm.U[p], yb_hat);

                        // compute angle of angle of attack in degrees: alpha = phi - twist - pitch
                        wt->adm.alpha[p] = wt->rad2deg * std::atan2(ub_x, ub_y) - wt->adm.twist[p] - wt->collPitch * wt->rad2deg - wfControlPitch;

                        // The idea is to interpolate the aero coeffs for the computed
                        // angle of attack, but the airfoils are discrete in radius.
                        // So interpolate the cl and cd from the two closest airfoils,
                        // previously stored.

                        // get the neighboring airfoils interpolation weights and labels
                        PetscReal af_w1 = wt->adm.iw[p][0],
                                  af_w2 = wt->adm.iw[p][1];
                        PetscInt  af_l1 = wt->adm.foilIds[p][0],
                                  af_l2 = wt->adm.foilIds[p][1];

                        // get the pointers to the aero tables
                        PetscReal *an1 = wt->foils[af_l1]->aoa;
                        PetscReal *an2 = wt->foils[af_l2]->aoa;
                        PetscReal *cl1 = wt->foils[af_l1]->cl;
                        PetscReal *cl2 = wt->foils[af_l2]->cl;
                        PetscReal *cd1 = wt->foils[af_l1]->cd;
                        PetscReal *cd2 = wt->foils[af_l2]->cd;

                        // get the size of the aero tables
                        PetscInt    s1   = wt->foils[af_l1]->size;
                        PetscInt    s2   = wt->foils[af_l2]->size;

                        // find interpolation weights for the aero coeffs
                        PetscReal w[2];
                        PetscInt  l[2];

                        // interpolate coeffs for the first airfoil
                        findInterpolationWeigths(w, l, an1, s1, wt->adm.alpha[p]);
                        PetscReal cl_1 = w[0]*cl1[l[0]] + w[1]*cl1[l[1]];
                        PetscReal cd_1 = w[0]*cd1[l[0]] + w[1]*cd1[l[1]];

                        // interpolate coeffs for the second airfoil
                        findInterpolationWeigths(w, l, an2, s2, wt->adm.alpha[p]);
                        PetscReal cl_2 = w[0]*cl2[l[0]] + w[1]*cl2[l[1]];
                        PetscReal cd_2 = w[0]*cd2[l[0]] + w[1]*cd2[l[1]];

                        // interpolate between the airfoils and save the values
                        wt->adm.Cl[p] = af_w1 * cl_1 + af_w2 * cl_2;
                        wt->adm.Cd[p] = af_w1 * cd_1 + af_w2 * cd_2;

                        // apply tip/root loss correction factors (from AeroDyn theory manual)
                        PetscReal F;
                        {
                            PetscReal g       = 1.0;
                            PetscReal r_p_mag = nMag(r_p);
                            PetscReal inflow  = fabs(wt->adm.alpha[p] + wt->adm.twist[p])*wt->deg2rad;

                            PetscReal ftip    = (wt->rTip - r_p_mag) / (r_p_mag * std::sin(inflow));
                            PetscReal Ftip    = (2.0 / M_PI) * std::acos(PetscMin(1.0, std::exp(-g * ((PetscReal)wt->nBlades/2.0) * ftip)));

                            PetscReal froot   = (r_p_mag - wt->rHub) / (r_p_mag * std::sin(inflow));
                            PetscReal Froot   = (2.0 / M_PI) * std::acos(PetscMin(1.0, std::exp(-g * ((PetscReal)wt->nBlades/2.0) * froot)));

                            F = Ftip * Froot;

                            wt->adm.Cl[p] *= F;
                            wt->adm.Cd[p] *= F;
                        }

                        // get velocity magnitude
                        PetscReal uMag = sqrt(ub_x*ub_x + ub_y*ub_y);

                        // compute scalar lift and drag per density
                        PetscReal lift_s = 0.5 * uMag * uMag * wt->adm.chord[p] * wt->adm.dr[p] * wt->adm.Cl[p];
                        PetscReal drag_s = 0.5 * uMag * uMag * wt->adm.chord[p] * wt->adm.dr[p] * wt->adm.Cd[p];

                        // define the aerodynamic reference frame (x as the flow, y along the blade, z will be lift)
                        Cmpnts xa_hat = nUnit(wt->adm.U[p]);
                        Cmpnts ya_hat = nUnit(r_p);
                        Cmpnts za_hat = nUnit(nCross(xa_hat, ya_hat));

                        // compute lift and drag vectors in the cartesian frame
                        Cmpnts lift_v = nScale(lift_s, za_hat);
                        Cmpnts drag_v = nScale(drag_s, xa_hat);

                        // compute body force (flow on blade)
                        wt->adm.B[p]  = nSum(lift_v, drag_v);

                        // reverse the body force (blade on flow)
                        mScale(-1.0, wt->adm.B[p]);

                        // compute axial force
                        wt->adm.axialF[p] = nDot(wt->adm.B[p], wt->rtrAxis);

                        // compute tangential force
                        wt->adm.tangtF[p] = nDot(wt->adm.B[p], yb_hat);

                        // cumulate rotor thrust
                        wt->adm.rtrThrust += wt->adm.axialF[p] * wt->adm.solidity[p] * constants->rho;

                        // cumulate rotor torque
                        wt->adm.rtrTorque += wt->adm.tangtF[p] * wt->adm.solidity[p] * constants->rho *
                                             nMag(r_p) * std::cos(wt->deg2rad*wt->precone);
                    }
                }

                // aerodynamic power
                wt->adm.aeroPwr = wt->adm.rtrTorque * wt->rtrOmega;
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // get thrust coefficient
                PetscReal Ct_t = wt->uadm.Ct;

                // apply wind farm controller
                if(farm->farmControlActive[t])
                {
                    Ct_t += wt->wfControlCt;
                }

                // get induction factor
                PetscReal a_t = wt->uadm.axiInd;

                // get reference velocity
                PetscReal Uref;

                if(wt->uadm.sampleType == "rotorUpstream")
                {
                    Uref = wt->upPoints->Uref;
                }
                else if(wt->uadm.sampleType == "givenVelocity")
                {
                    Uref = wt->uadm.Uref;
                }
                else if(wt->uadm.sampleType == "rotorDisk")
                {
                    Uref = wt->uadm.rtrAvgMagU;
                }

                // cumulate rotor thrust at this time step
                wt->uadm.rtrThrust = 0.0;

                // number of points in the AD mesh
                PetscInt npts_t = wt->uadm.nPoints;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // compute blade force / density
                    PetscReal bMag = 0.5 * Uref * Uref * wt->uadm.dA[p] * Ct_t;

                    // compute body force (rotor on flow)
                    wt->uadm.B[p]  = nScale(bMag, wt->rtrAxis);

                    // compute axial force
                    wt->uadm.axialF[p] = nDot(wt->uadm.B[p], wt->rtrAxis);

                    // cumulate rotor thrust
                    wt->uadm.rtrThrust += wt->uadm.axialF[p] * constants->rho;
                }

                // aerodynamic power
                wt->uadm.aeroPwr = wt->uadm.rtrThrust * (1-a_t) * Uref;
            }
            // actuator line model
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // apply wind farm control
                PetscReal wfControlPitch = 0.0;

                // apply wind farm controller
                if(farm->farmControlActive[t])
                {
                    wfControlPitch = wt->wfControlCollPitch;
                }

                // zero the rotor torque at this time step
                wt->alm.rtrTorque = 0.0;

                // cumulate rotor thrust at this time step
                wt->alm.rtrThrust = 0.0;

                // number of points in the AL mesh
                PetscInt npts_t = wt->alm.nPoints;

                // loop over the AL mesh points
                for(p=0; p<npts_t; p++)
                {
                    // do not test if this point is inside this processor.
                    // If at least 1 AL point is inside this processor we have
                    // to compute the parameters at all AD points.
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
                            xb_hat = nScale(-1.0, rtrAxis);
                            yb_hat = nCross(zb_hat, xb_hat);
                        }
                        else if(wt->rotDir == "ccw")
                        {
                            zb_hat = nUnit(r_p);
                                     mScale(-1.0, zb_hat);
                            xb_hat = nScale(-1.0, rtrAxis);
                            yb_hat = nCross(zb_hat, xb_hat);
                        }

                        // transform the velocity in the bladed reference frame in
                        // order to compute the angle of attack. Radial (z) component
                        // is not considered
                        PetscReal ub_x = nDot(wt->alm.U[p], xb_hat);
                        PetscReal ub_y = nDot(wt->alm.U[p], yb_hat);

                        // compute angle of angle of attack in degrees: alpha = phi - twist - pitch
                        wt->alm.alpha[p] = wt->rad2deg * std::atan2(ub_x, ub_y) - wt->alm.twist[p] - wt->collPitch * wt->rad2deg - wfControlPitch;

                        // The idea is to interpolate the aero coeffs for the computed
                        // angle of attack, but the airfoils are discrete in radius.
                        // So interpolate the cl and cd from the two closest airfoils,
                        // previously stored.

                        // get the neighboring airfoils interpolation weights and labels
                        PetscReal af_w1 = wt->alm.iw[p][0],
                                  af_w2 = wt->alm.iw[p][1];
                        PetscInt  af_l1 = wt->alm.foilIds[p][0],
                                  af_l2 = wt->alm.foilIds[p][1];

                        // get the pointers to the aero tables
                        PetscReal *an1 = wt->foils[af_l1]->aoa;
                        PetscReal *an2 = wt->foils[af_l2]->aoa;
                        PetscReal *cl1 = wt->foils[af_l1]->cl;
                        PetscReal *cl2 = wt->foils[af_l2]->cl;
                        PetscReal *cd1 = wt->foils[af_l1]->cd;
                        PetscReal *cd2 = wt->foils[af_l2]->cd;

                        // get the size of the aero tables
                        PetscInt    s1   = wt->foils[af_l1]->size;
                        PetscInt    s2   = wt->foils[af_l2]->size;

                        // find interpolation weights for the aero coeffs
                        PetscReal w[2];
                        PetscInt  l[2];

                        // interpolate coeffs for the first airfoil
                        findInterpolationWeigths(w, l, an1, s1, wt->alm.alpha[p]);
                        PetscReal cl_1 = w[0]*cl1[l[0]] + w[1]*cl1[l[1]];
                        PetscReal cd_1 = w[0]*cd1[l[0]] + w[1]*cd1[l[1]];

                        // interpolate coeffs for the second airfoil
                        findInterpolationWeigths(w, l, an2, s2, wt->alm.alpha[p]);
                        PetscReal cl_2 = w[0]*cl2[l[0]] + w[1]*cl2[l[1]];
                        PetscReal cd_2 = w[0]*cd2[l[0]] + w[1]*cd2[l[1]];

                        // interpolate between the airfoils and save the values
                        wt->alm.Cl[p] = af_w1 * cl_1 + af_w2 * cl_2;
                        wt->alm.Cd[p] = af_w1 * cd_1 + af_w2 * cd_2;

                        // apply tip/root loss correction factors (from AeroDyn theory manual)
                        PetscReal F;
                        {
                            PetscReal g       = 1.0;
                            PetscReal r_p_mag = nMag(r_p);
                            PetscReal inflow  = fabs(wt->alm.alpha[p] + wt->alm.twist[p])*wt->deg2rad;

                            PetscReal ftip    = (wt->rTip - r_p_mag) / (r_p_mag * std::sin(inflow));
                            PetscReal Ftip    = (2.0 / M_PI) * std::acos(PetscMin(1.0, std::exp(-g * ((PetscReal)wt->nBlades/2.0) * ftip)));

                            PetscReal froot   = (r_p_mag - wt->rHub) / (r_p_mag * std::sin(inflow));
                            PetscReal Froot   = (2.0 / M_PI) * std::acos(PetscMin(1.0, std::exp(-g * ((PetscReal)wt->nBlades/2.0) * froot)));

                            F = Ftip * Froot;

                            wt->alm.Cl[p] *= F;
                            wt->alm.Cd[p] *= F;
                        }

                        // get velocity magnitude
                        PetscReal uMag = sqrt(ub_x*ub_x + ub_y*ub_y);

                        // compute scalar lift and drag per density
                        PetscReal lift_s = 0.5 * uMag * uMag * wt->alm.chord[p] * wt->alm.dr[p] * wt->alm.Cl[p];
                        PetscReal drag_s = 0.5 * uMag * uMag * wt->alm.chord[p] * wt->alm.dr[p] * wt->alm.Cd[p];

                        // define the aerodynamic reference frame (x as the flow, y along the blade, z will be lift)
                        Cmpnts xa_hat = nUnit(wt->alm.U[p]);
                        Cmpnts ya_hat = nUnit(r_p);
                        Cmpnts za_hat = nUnit(nCross(xa_hat, ya_hat));

                        // compute lift and drag vectors in the cartesian frame
                        Cmpnts lift_v = nScale(lift_s, za_hat);
                        Cmpnts drag_v = nScale(drag_s, xa_hat);

                        // compute body force (flow on blade)
                        wt->alm.B[p]  = nSum(lift_v, drag_v);

                        // reverse the body force (blade on flow)
                        mScale(-1.0, wt->alm.B[p]);

                        // compute axial force
                        wt->alm.axialF[p] = nDot(wt->alm.B[p], wt->rtrAxis);

                        // compute tangential force
                        wt->alm.tangtF[p] = nDot(wt->alm.B[p], yb_hat);

                        // cumulate rotor thrust
                        wt->alm.rtrThrust += wt->alm.axialF[p] * wt->alm.solidity[p] * constants->rho;

                        // cumulate rotor torque
                        wt->alm.rtrTorque += wt->alm.tangtF[p] * wt->alm.solidity[p] * constants->rho *
                                             nMag(r_p) * std::cos(wt->deg2rad*wt->precone);
                    }
                }

                // aerodynamic power
                wt->alm.aeroPwr = wt->alm.rtrTorque * wt->rtrOmega;

            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                PetscReal a_t, Uref, bMag, Ct_t = wt->afm.Ct;

                if(farm->farmControlActive[t])
                {
                    Ct_t += wt->wfControlCt;
                }

                // reference velocity is extrapolated using induction coefficient - Ct relation
                if(wt->afm.sampleType == "momentumTheory")
                {
                    // compute induction factor
                    a_t     = (1.0 - std::sqrt(1.0 - Ct_t)) / 2.0;

                    // compute freestream velocity
                    Uref  = nMag(wt->afm.U) / (1.0 - a_t);

                    // compute body force
                    bMag = 0.5 * Uref * Uref * M_PI * wt->rTip * wt->rTip * Ct_t;
                }
                // disk based Ct is used
                else if(wt->afm.sampleType == "rotorDisk" || wt->afm.sampleType == "integral")
                {
                    // compute freestream velocity
                    Uref  = nMag(wt->afm.U);

                    // compute body force
                    bMag = 0.5 * Uref * Uref * M_PI * wt->rTip * wt->rTip * Ct_t;

                    // compute induction factor (zero so that power is correct bcs Uref is already at disk)
                    a_t = 0.0;
                }

                wt->afm.B = nScale(bMag, wt->rtrAxis);

                // aerodynamic thrust
                wt->afm.rtrThrust = bMag * constants->rho;

                // aerodynamic power
                wt->afm.aeroPwr   = wt->afm.rtrThrust * (1-a_t) * Uref;
            }

            // electric power (only for AD/AL models)
            if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
            {
                // second test inside because this variable is not defined for UADM
                if(wt->genControllerType != "none")
                {
                    wt->genPwr = wt->genTorque * (wt->rtrOmega * wt->gbxRatioG2R) * wt->genEff;
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode projectBladeForce(farm_ *farm)
{
    constants_       *constants = farm->access->constants;
    mesh_            *mesh  = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, c, p;

    Cmpnts           ***sCat;

    Cmpnts           ***cent;

    PetscReal        ***aj;

    PetscMPIInt           rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // global sums of thrust and torque over the entire wind farm
    // used for the computation of the body force projection error
    PetscReal lThrustSum   = 0.0, gThrustSum   = 0.0;
    PetscReal lTorqueSum   = 0.0, gTorqueSum   = 0.0;
    PetscReal lThrustBFSum = 0.0, gThrustBFSum = 0.0;
    PetscReal lTorqueBFSum = 0.0, gTorqueBFSum = 0.0;

    // zero out the body force
    VecSet(farm->lsourceFarmCat, 0.0);

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

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
                // projection parameters
                PetscReal rPrj = wt->eps * wt->prjNSigma;

                // projection distance
                PetscReal eps  = wt->eps;

                // number of points in the AD mesh
                PetscInt npts_t = wt->adm.nPoints;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // do not test if this point is controlled since body force can
                    // come from other processor, the test is done on the projection radius
                    {
                        // save this point locally for speed
                        Cmpnts point_p    = wt->adm.points[p];

                        // this point position from COR
                        Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                        // save this point solidity for speed
                        PetscReal solidity_p = wt->adm.solidity[p];

                        // loop in sphere points
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

                            if(r_c_mag<rPrj)
                            {
                                // compute projection factor
                                PetscReal pf
                                =
                                /*
                                std::exp
                                (
                                    -(r_c_mag *r_c_mag) /
                                     (2.0 * eps * eps)
                                ) /
                                (
                                    pow(eps,    3) *
                                    pow(2.0*M_PI, 1.5)
                                );
                                */
                                std::exp
                                (
                                    -(r_c_mag / eps)*
                                     (r_c_mag / eps)
                                ) /
                                (
                                    pow(eps,    3) *
                                    pow(M_PI, 1.5)
                                );

                                Cmpnts bfCell = nScale(solidity_p * pf, wt->adm.B[p]);

                                sCat[k][j][i].x += bfCell.x;
                                sCat[k][j][i].y += bfCell.y;
                                sCat[k][j][i].z += bfCell.z;

                                // cumulate wind farm BF for projection error
                                PetscReal vCell = 1.0 / aj[k][j][i];
                                Cmpnts thrustBF = nScale(vCell*constants->rho, bfCell);
                                Cmpnts torqueBF = nScale(vCell*constants->rho*nMag(r_p)*std::cos(wt->deg2rad*wt->precone), bfCell);

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

                                // cumulate contribution from this AD point at this cell
                                lThrustBFSum += nDot(thrustBF, wt->rtrAxis);
                                lTorqueBFSum += nDot(torqueBF, yb_hat);
                            }
                        }
                    }
                }
            }
            // uniform actuator disk model
            else if((*farm->turbineModels[t]) == "uniformADM")
            {
                // projection parameters
                PetscReal rPrj = wt->eps * wt->prjNSigma;

                // projection distance
                PetscReal eps  = wt->eps;

                // number of points in the AD mesh
                PetscInt npts_t = wt->uadm.nPoints;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // do not test if this point is controlled since body force can
                    // come from other processor, the test is done on the projection radius
                    {
                        // save this point locally for speed
                        Cmpnts point_p    = wt->uadm.points[p];

                        // loop in sphere points
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

                            if(r_c_mag<rPrj)
                            {
                                // compute projection factor
                                PetscReal pf
                                =
                                /*
                                std::exp
                                (
                                    -(r_c_mag *r_c_mag) /
                                     (2.0 * eps * eps)
                                ) /
                                (
                                    pow(eps,    3) *
                                    pow(2.0*M_PI, 1.5)
                                );
                                */
                                std::exp
                                (
                                    -(r_c_mag / eps)*
                                     (r_c_mag / eps)
                                ) /
                                (
                                    pow(eps,    3) *
                                    pow(M_PI, 1.5)
                                );

                                Cmpnts bfCell = nScale(pf, wt->uadm.B[p]);

                                sCat[k][j][i].x += bfCell.x;
                                sCat[k][j][i].y += bfCell.y;
                                sCat[k][j][i].z += bfCell.z;

                                // cumulate wind farm BF for projection error
                                PetscReal vCell    = 1.0 / aj[k][j][i];
                                Cmpnts thrustBF = nScale(vCell*constants->rho, bfCell);

                                // cumulate contribution from this AD point at this cell
                                lThrustBFSum += nDot(thrustBF, wt->rtrAxis);
                                lTorqueBFSum += 0.0;
                            }
                        }
                    }
                }
            }
            // actuator line model
            else if((*farm->turbineModels[t]) == "ALM")
            {
                // projection type
                PetscInt projectionType = 0;

                if(wt->alm.projectionType == "anisotropic")
                {
                    projectionType = 1;
                }

                // number of points in the AL mesh
                PetscInt npts_t = wt->alm.nPoints;

                // loop over the AD mesh points
                for(p=0; p<npts_t; p++)
                {
                    // do not test if this point is controlled since body force can
                    // come from other processor, the test is done on the projection radius
                    {
                        // save this point locally for speed
                        Cmpnts point_p    = wt->alm.points[p];

                        // this point position from COR
                        Cmpnts r_p  = nSub(point_p, wt->rotCenter);

                        // save this point solidity for speed
                        PetscReal solidity_p = wt->alm.solidity[p];

                        // loop in sphere points
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

                                // radial mesh cell size size
                                PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                                // get projection epsilon
                                PetscReal eps_x =     wt->alm.chord[p] * wt->eps_x,
                                          eps_y =     wt->alm.thick[p] * wt->eps_y,
                                          eps_z =     drval * wt->eps_z;

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

                                Cmpnts bfCell = nScale(solidity_p * pf, wt->alm.B[p]);

                                sCat[k][j][i].x += bfCell.x;
                                sCat[k][j][i].y += bfCell.y;
                                sCat[k][j][i].z += bfCell.z;

                                // cumulate wind farm BF for projection error
                                PetscReal vCell    = 1.0 / aj[k][j][i];
                                Cmpnts thrustBF = nScale(vCell*constants->rho, bfCell);
                                Cmpnts torqueBF = nScale(vCell*constants->rho*nMag(r_p)*std::cos(wt->deg2rad*wt->precone), bfCell);

                                // cumulate contribution from this AD point at this cell
                                lThrustBFSum += nDot(thrustBF, wt->rtrAxis);
                                lTorqueBFSum += nDot(torqueBF, yb_hat);
                            }
                        }
                    }
                }
            }
            else if((*farm->turbineModels[t]) == "AFM")
            {
                // projection distance
                PetscReal eps_x  = wt->eps_x,
                          eps_y  = wt->eps_y,
                          eps_z  = wt->eps_z;

                // save this point locally for speed
                Cmpnts point_p    = wt->afm.point;

                // loop in sphere points
                for(c=0; c<wt->nControlled; c++)
                {
                    // cell indices
                    PetscInt i = wt->controlledCells[c].i,
                             j = wt->controlledCells[c].j,
                             k = wt->controlledCells[c].k;

                    // compute distance from mesh cell to AF point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    //if(r_c_mag < 2.0*wt->rTip)
                    {
                        // compute projection factor
                        PetscReal pf
                        =
                        std::exp
                        (
                            -r_c.x * r_c.x / (eps_x * eps_x)
                            -r_c.y * r_c.y / (eps_y * eps_y)
                            -r_c.z * r_c.z / (eps_z * eps_z)
                        ) /
                        (
                            eps_x * eps_y * eps_z *
                            pow(M_PI, 1.5)
                        );

                        Cmpnts bfCell = nScale(pf, wt->afm.B);

                        sCat[k][j][i].x += bfCell.x;
                        sCat[k][j][i].y += bfCell.y;
                        sCat[k][j][i].z += bfCell.z;

                        // cumulate wind farm BF for projection error
                        PetscReal vCell = 1.0 / aj[k][j][i];
                        Cmpnts thrustBF = nScale(vCell*constants->rho, bfCell);

                        // cumulate contribution from this AF point at this cell
                        lThrustBFSum += nDot(thrustBF, wt->rtrAxis);
                        lTorqueBFSum += 0.0;
                    }
                }
            }
        }
    }

    // restore the array
    DMDAVecRestoreArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);
    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);

    // scatter local to local
    DMLocalToLocalBegin(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);
    DMLocalToLocalEnd(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);

    // output screen information: global parameters and projection errors
    if(farm->dbg)
    {
        for(t=0; t<farm->size; t++)
        {
            windTurbine  *wt = farm->wt[t];
            PetscReal    aeroPwr = 0.0, genPwr = 0.0,
                         Thr     = 0.0, Trq    = 0.0;

            if(wt->turbineControlled)
            {
                // reduce turbine parameters
                if((*farm->turbineModels[t]) == "ADM")
                {
                    MPI_Reduce(&(wt->adm.aeroPwr),   &aeroPwr, 1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.rtrThrust), &Thr,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.rtrTorque), &Trq,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    MPI_Reduce(&(wt->uadm.aeroPwr),   &aeroPwr, 1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->uadm.rtrThrust), &Thr,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                }
                else if((*farm->turbineModels[t]) == "ALM")
                {
                    MPI_Reduce(&(wt->alm.aeroPwr),   &aeroPwr, 1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.rtrThrust), &Thr,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.rtrTorque), &Trq,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                }
                else if((*farm->turbineModels[t]) == "AFM")
                {
                    MPI_Reduce(&(wt->afm.aeroPwr),   &aeroPwr, 1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->afm.rtrThrust), &Thr,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                }

                // electric power
                if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                {
                    // second test inside because this variable is not defined for UADM
                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->genPwr), &genPwr, 1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }
                }


                if(rank == wt->writerRank)
                {
                    // divide the summed up turbine parameters by the n of processors
                    // that control this turbine since all had the data in the reduce operation

                    aeroPwr /= wt->nProcsTrb;
                    Thr     /= wt->nProcsTrb;
                    Trq     /= wt->nProcsTrb;

                    if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                    {
                        // second test inside because this variable is not defined for UADM
                        if(wt->genControllerType != "none")
                        {
                            genPwr  /= wt->nProcsTrb;
                        }
                    }

                    // print turbine level information
                    if(wt->dbg)
                    {
                        if((*farm->turbineModels[t]) == "uniformADM" || (*farm->turbineModels[t]) == "AFM")
                        {
                            printf("Turbine %s: aero pwr = %.5f MW, Thrust on rtr axis = %.5f KN\n",(*farm->turbineIds[t]).c_str(), aeroPwr/1.0e6, Thr/1000.0);
                        }
                        else
                        {
                            printf("Turbine %s: electrical pwr = %.5f MW, aero pwr = %.5f MW, Thrust on rtr axis = %.5f KN, Torque on rot dir = %.5f KNm\n",(*farm->turbineIds[t]).c_str(), genPwr/1.0e6, aeroPwr/1.0e6, Thr/1000.0, Trq/1000.0);
                        }
                    }

                    // cumulate forces from this turbine
                    lThrustSum += Thr;
                    lTorqueSum += Trq;
                }
            }
        }

        // reduce global force
        MPI_Reduce(&lThrustSum, &gThrustSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        MPI_Reduce(&lTorqueSum, &gTorqueSum , 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // reduce the global body force
        MPI_Reduce(&lThrustBFSum, &gThrustBFSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);
        MPI_Reduce(&lTorqueBFSum, &gTorqueBFSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // compute projection errors
        if(!rank)
        {
            // relative percentage projection errors
            PetscReal errThr = fabs(gThrustBFSum - gThrustSum) / gThrustSum * 100.0;

            word percent = "%";

            // all models are uniformADM/AFM (no torque)
            if(gTorqueSum <= 1e-10)
            {
                // print projection error info
                printf("Wind turbine rotors: global projection error on thrust = %.5f %s\n", errThr, percent.c_str());
            }
            // mixed turbine models
            else
            {
                // compute error on torque
                PetscReal errTrq = fabs(gTorqueBFSum - gTorqueSum) / gTorqueSum * 100.0;

                // print projection error info
                printf("Wind turbine rotors: global projection error on thrust = %.5f %s, global projection error on torque = %.5f %s\n", errThr, percent.c_str(), errTrq, percent.c_str());
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode projectTowerForce(farm_ *farm)
{
    constants_       *constants = farm->access->constants;
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, c, p;

    Cmpnts           ***sCat;

    Cmpnts           ***cent;

    PetscReal        ***aj;

    PetscMPIInt           rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // global sums of tangential tower force over the entire wind farm
    // used for the computation of the body force projection error
    PetscReal lTangSum   = 0, gTangSum   = 0;
    PetscReal lTangBFSum = 0, gTangBFSum = 0;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if tower is modeled in this turbine
        if(wt->includeTwr)
        {
            // test if this processor controls this tower
            if(wt->twr.nControlled)
            {
                // zero this tower thrust
                wt->twr.twrThrust = 0.0;

                // projection parameters
                PetscReal rPrj = wt->twr.eps * wt->twr.prjNSigma,
                          eps  = wt->twr.eps;

                // number of points in the tower mesh
                PetscInt npts_t = wt->twr.nPoints;

                // loop over the tower mesh points
                for(p=0; p<npts_t; p++)
                {
                    // get this point velocity
                    PetscReal uMag = nMag(wt->twr.U[p]);

                    // compute scalar body force
                    PetscReal magB = 0.5 * uMag * uMag * wt->twr.dA[p] * wt->twr.Cd;

                    // get force direction unit vector (x in the local aero frame)
                    Cmpnts xa_hat = nUnit(wt->twr.U[p]);

                    // compute vector body force
                    wt->twr.B[p] = nScale(-magB, xa_hat);

                    // compute axial force
                    wt->twr.tangF[p] = nDot(wt->twr.B[p], wt->rtrAxis);

                    // cumulate this tower thrust
                    wt->twr.twrThrust += wt->twr.tangF[p] * constants->rho;

                    // save this point locally for speed
                    Cmpnts point_p    = wt->twr.points[p];

                    // loop in tower cells
                    for(c=0; c<wt->twr.nControlled; c++)
                    {
                        // cell indices
                        PetscInt i = wt->twr.controlledCells[c].i,
                            j = wt->twr.controlledCells[c].j,
                            k = wt->twr.controlledCells[c].k;

                        // compute distance from mesh cell to tower point
                        Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                        // compute magnitude
                        PetscReal r_c_mag = nMag(r_c);

                        if(r_c_mag<rPrj)
                        {
                            // compute projection factor
                            PetscReal pf
                            =
                            /*
                            std::exp
                            (
                                -(r_c_mag *r_c_mag) /
                                 (2.0 * eps * eps)
                            ) /
                            (
                                pow(eps,    3) *
                                pow(2.0*M_PI, 1.5)
                            );
                            */
                            std::exp
                            (
                                -(r_c_mag / eps)*
                                 (r_c_mag / eps)
                            ) /
                            (
                                pow(eps,    3) *
                                pow(M_PI, 1.5)
                            );

                            Cmpnts bfCell = nScale(pf, wt->twr.B[p]);

                            sCat[k][j][i].x += bfCell.x;
                            sCat[k][j][i].y += bfCell.y;
                            sCat[k][j][i].z += bfCell.z;

                            // cumulate tower BF for projection error
                            PetscReal vCell    = 1.0 / aj[k][j][i];
                            Cmpnts tangBF = nScale(vCell*constants->rho, bfCell);

                            // cumulate contribution from this AD point at this cell
                            lTangBFSum += nDot(tangBF, wt->rtrAxis);
                        }
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    // scatter local to local
    DMLocalToLocalBegin(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);
    DMLocalToLocalEnd(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);

    if(farm->dbg)
    {
        // output screen information: global parameters and projection errors
        for(t=0; t<farm->size; t++)
        {
            windTurbine  *wt = farm->wt[t];

            // all procs enter this test
            if(wt->includeTwr)
            {
                PetscReal       Tang        = 0.0;

                // reduce tower thrust
                MPI_Reduce(&(wt->twr.twrThrust), &Tang, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

                if(!rank)
                {
                    // divide the summed up turbine parameters by the n of processors
                    // that control this turbine since all had the data in the reduce operation
                    Tang /= wt->twr.nProcsTwr;

                    // cumulate forces from this turbine
                    lTangSum += Tang;
                }
            }
        }

        // reduce the global tangential force
        MPI_Reduce(&lTangSum, &gTangSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // reduce the global body force
        MPI_Reduce(&lTangBFSum, &gTangBFSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // compute projection errors if overall tower force is present
        if(!rank && gTangSum > 1e-5)
        {
            // relative percentage projection errors
            PetscReal errTang = fabs(gTangBFSum - gTangSum) / gTangSum * 100.0;

            word percent = "%";

            // print projection error info
            printf("Wind turbine towers: global projection error on tang. force = %.5f %s, tangBFSum = %.5f KN, tangFSum = %.5f KN\n", errTang, percent.c_str(), gTangBFSum/1000.0, gTangSum/1000.0);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode projectNacelleForce(farm_ *farm)
{
    constants_       *constants = farm->access->constants;
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, c, p;

    Cmpnts           ***sCat;

    Cmpnts           ***cent;

    PetscReal        ***aj;

    PetscMPIInt           rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // global sums of tangential tower force over the entire wind farm
    // used for the computation of the body force projection error
    PetscReal lTangSum   = 0, gTangSum   = 0;
    PetscReal lTangBFSum = 0, gTangBFSum = 0;

    DMDAVecGetArray(fda, mesh->lCent, &cent);
    DMDAVecGetArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecGetArray(da,  mesh->lAj,    &aj);

    // loop over each wind turbine
    for(t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        // test if nacelle is modeled in this turbine
        if(wt->includeNacelle)
        {
            // test if this processor controls this nacelle
            if(wt->nac.nControlled)
            {
                // zero this nacelle thrust
                wt->nac.nacThrust = 0.0;

                // projection parameters
                PetscReal rPrj = wt->nac.eps * wt->nac.prjNSigma,
                          eps  = wt->nac.eps;

                // get this point velocity
                PetscReal uMag = nMag(wt->nac.U);

                // compute scalar body force
                PetscReal magB = 0.5 * uMag * uMag * wt->nac.A * wt->nac.Cd;

                // get force direction unit vector (x in the local aero frame)
                Cmpnts xa_hat = nUnit(wt->nac.U);

                // compute vector body force
                wt->nac.B = nScale(-magB, xa_hat);

                // compute axial force
                wt->nac.tangF = nDot(wt->nac.B, wt->rtrAxis);

                // cumulate this nacelle thrust
                wt->nac.nacThrust += wt->nac.tangF * constants->rho;

                // save this point locally for speed
                Cmpnts point_p = wt->nac.point;

                // loop in nacelle cells
                for(c=0; c<wt->nac.nControlled; c++)
                {
                    // cell indices
                    PetscInt i = wt->nac.controlledCells[c].i,
                             j = wt->nac.controlledCells[c].j,
                             k = wt->nac.controlledCells[c].k;

                    // compute distance from mesh cell to tower point
                    Cmpnts r_c = nSub(point_p, cent[k][j][i]);

                    // compute magnitude
                    PetscReal r_c_mag = nMag(r_c);

                    if(r_c_mag<rPrj)
                    {
                        // compute projection factor
                        PetscReal pf
                        =
                        /*
                        std::exp
                        (
                            -(r_c_mag *r_c_mag) /
                             (2.0 * eps * eps)
                        ) /
                        (
                            pow(eps,    3) *
                            pow(2.0*M_PI, 1.5)
                        );
                        */
                        std::exp
                        (
                            -(r_c_mag / eps)*
                             (r_c_mag / eps)
                        ) /
                        (
                            pow(eps,    3) *
                            pow(M_PI, 1.5)
                        );

                        Cmpnts bfCell = nScale(pf, wt->nac.B);

                        sCat[k][j][i].x += bfCell.x;
                        sCat[k][j][i].y += bfCell.y;
                        sCat[k][j][i].z += bfCell.z;

                        // cumulate tower BF for projection error
                        PetscReal vCell    = 1.0 / aj[k][j][i];
                        Cmpnts tangBF = nScale(vCell*constants->rho, bfCell);

                        // cumulate contribution from this AD point at this cell
                        lTangBFSum += nDot(tangBF, wt->rtrAxis);
                    }
                }
            }
        }
    }

    DMDAVecRestoreArray(da,  mesh->lAj,    &aj);
    DMDAVecRestoreArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecRestoreArray(fda, mesh->lCent, &cent);

    // scatter local to local
    DMLocalToLocalBegin(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);
    DMLocalToLocalEnd(fda, farm->lsourceFarmCat, INSERT_VALUES, farm->lsourceFarmCat);

    if(farm->dbg)
    {
        // output screen information: global parameters and projection errors
        for(t=0; t<farm->size; t++)
        {
            windTurbine  *wt = farm->wt[t];

            // all procs enter this test
            if(wt->includeNacelle)
            {
                PetscReal       Tang        = 0.0;

                // reduce tower thrust
                MPI_Reduce(&(wt->nac.nacThrust), &Tang, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

                if(!rank)
                {
                    // divide the summed up turbine parameters by the n of processors
                    // that control this turbine since all had the data in the reduce operation
                    Tang /= wt->nac.nProcsNac;

                    // cumulate forces from this turbine
                    lTangSum += Tang;
                }
            }
        }

        // reduce the global tangential force
        MPI_Reduce(&lTangSum, &gTangSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // reduce the global body force
        MPI_Reduce(&lTangBFSum, &gTangBFSum, 1, MPIU_REAL, MPIU_SUM, 0, mesh->MESH_COMM);

        // compute projection errors if overall tower force is present
        if(!rank && gTangSum > 1e-5)
        {
            // relative percentage projection errors
            PetscReal errTang = fabs(gTangBFSum - gTangSum) / gTangSum * 100.0;

            word percent = "%";

            // print projection error info
            printf("Wind turbine nacelles: global projection error on tang. force = %.5f %s, tangBFSum = %.5f  KN, tangFSum = %.5f KN\n", errTang, percent.c_str(), gTangBFSum/1000.0, gTangSum/1000.0);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode bodyForceCartesian2Contravariant(farm_ *farm)
{
    mesh_            *mesh = farm->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k, t, c, p;

    Cmpnts           ***sCat, ***sCont;

    Cmpnts           ***icsi, ***jeta, ***kzet;

    lxs = xs; if (xs==0) lxs = xs+1; lxe = xe; if (xe==mx) lxe = xe-1;
    lys = ys; if (ys==0) lys = ys+1; lye = ye; if (ye==my) lye = ye-1;
    lzs = zs; if (zs==0) lzs = zs+1; lze = ze; if (ze==mz) lze = ze-1;

    // zero out the body force
    VecSet(farm->sourceFarmCont, 0.0);

    DMDAVecGetArray(fda, farm->lsourceFarmCat, &sCat);
    DMDAVecGetArray(fda, farm->sourceFarmCont, &sCont);

    DMDAVecGetArray(fda, mesh->lICsi,  &icsi);
    DMDAVecGetArray(fda, mesh->lJEta,  &jeta);
    DMDAVecGetArray(fda, mesh->lKZet,  &kzet);

    // interpolate BF to cell faces
    for(k=lzs; k<lze; k++)
    {
        for(j=lys; j<lye; j++)
        {
            for(i=lxs; i<lxe; i++)
            {
                // interpolate to the i-faces
                Cmpnts sFaceI;
                sFaceI.x = 0.5 * (sCat[k][j][i].x + sCat[k][j][i+1].x);
                sFaceI.y = 0.5 * (sCat[k][j][i].y + sCat[k][j][i+1].y);
                sFaceI.z = 0.5 * (sCat[k][j][i].z + sCat[k][j][i+1].z);

                // body force i-flux
                sCont[k][j][i].x
                =
                (
                    sFaceI.x * icsi[k][j][i].x +
                    sFaceI.y * icsi[k][j][i].y +
                    sFaceI.z * icsi[k][j][i].z
                );

                // interpolate to the j-faces
                Cmpnts sFaceJ;
                sFaceJ.x = 0.5 * (sCat[k][j][i].x + sCat[k][j+1][i].x);
                sFaceJ.y = 0.5 * (sCat[k][j][i].y + sCat[k][j+1][i].y);
                sFaceJ.z = 0.5 * (sCat[k][j][i].z + sCat[k][j+1][i].z);

                // body force j-flux
                sCont[k][j][i].y
                =
                (
                    sFaceJ.x * jeta[k][j][i].x +
                    sFaceJ.y * jeta[k][j][i].y +
                    sFaceJ.z * jeta[k][j][i].z
                );

                // interpolate to the k-faces
                Cmpnts sFaceK;
                sFaceK.x = 0.5 * (sCat[k][j][i].x + sCat[k+1][j][i].x);
                sFaceK.y = 0.5 * (sCat[k][j][i].y + sCat[k+1][j][i].y);
                sFaceK.z = 0.5 * (sCat[k][j][i].z + sCat[k+1][j][i].z);

                // body force k-flux
                sCont[k][j][i].z
                =
                (
                    sFaceK.x * kzet[k][j][i].x +
                    sFaceK.y * kzet[k][j][i].y +
                    sFaceK.z * kzet[k][j][i].z
                );
            }
        }
    }

    DMDAVecRestoreArray(fda, mesh->lICsi,  &icsi);
    DMDAVecRestoreArray(fda, mesh->lJEta,  &jeta);
    DMDAVecRestoreArray(fda, mesh->lKZet,  &kzet);

    DMDAVecRestoreArray(fda, farm->sourceFarmCont,  &sCont);
    DMDAVecRestoreArray(fda, farm->lsourceFarmCat, &sCat);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode computeMaxTipSpeed(farm_ *farm)
{
    mesh_    *mesh = farm->access->mesh;

    PetscReal lmaxTipSpeed = 0.0;
    PetscReal gmaxTipSpeed = 0.0;

    // loop over each wind turbine
    for(PetscInt t=0; t<farm->size; t++)
    {
        windTurbine *wt = farm->wt[t];

        PetscReal   tipSpeed = wt->rTip * wt->rtrOmega;

        if(tipSpeed > lmaxTipSpeed)
        {
            lmaxTipSpeed = tipSpeed;
        }
    }

    // compute the maximum among all processors
    MPI_Allreduce(&lmaxTipSpeed, &gmaxTipSpeed, 1, MPIU_REAL, MPIU_MAX, mesh->MESH_COMM);

    farm->maxTipSpeed = gmaxTipSpeed;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windTurbinesWrite(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    PetscMPIInt rank, t;

    clock_ *clock = farm->access->clock;
    io_    *io    = farm->access->io;

    word   turbineFolderName     = "./postProcessing/" + mesh->meshName + "/turbines";
    word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // create/initialize turbines directory (at simulation start only)
    if
    (
        clock->it == clock->itStart
    )
    {
        if(!rank)
        {
            errno = 0;
            PetscInt dirRes;

            dirRes = mkdir(turbineFolderName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create %s directory\n", turbineFolderName.c_str());
                fatalErrorInFunction("windTurbinesWrite",  error);
            }

            dirRes = mkdir(turbineFolderTimeName.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
                char error[512];
                sprintf(error, "could not create %s directory\n", turbineFolderTimeName.c_str());
                fatalErrorInFunction("windTurbinesWrite",  error);
            }

            // if directory already exist remove everything inside
            if(errno == EEXIST)
            {
                remove_subdirs(farm->access->mesh->MESH_COMM, turbineFolderTimeName.c_str());
            }
        }

        // ensure folder is there for every processor
        MPI_Barrier(mesh->MESH_COMM);

        // create turbines files
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(rank == wt->writerRank)
            {
                FILE *f;
                char fileName[80];
                sprintf(fileName, "%s/%s", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                f = fopen(fileName, "w");

                if(!f)
                {
                   char error[512];
                    sprintf(error, "cannot open file %s\n", fileName);
                    fatalErrorInFunction("windTurbinesWrite",  error);
                }
                else
                {

                    int width = -20;

                    word w0 = "time [s]";

                    word w1  = "rtrAvgMagU [m/s]";
                    word w2  = "rtrAvgUpMagU [m/s]";
                    word w3  = "rtrThrust [kN]";
                    word w4  = "aeroPwr [MW]";

                    word w5  = "CtInf [-]";
                    word w6  = "CtLoc [-]";
                    word w7  = "CtUp [-]";

                    word w8  = "rtrTorque [kNm]";
                    word w9  = "rtrOmega [rpm]";

                    word w10 = "genTorque [kNm]";
                    word w11 = "genPwr [MW]";
                    word w12 = "genOmega [rpm]";

                    word w13 = "collPitch [deg]";

                    word w14 = "flowAngle [deg]";
                    word w15 = "yawAngle [deg]";

                    word w16 = "azimuth [deg]";

                    // actuator disk model
                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str(), width, w8.c_str(), width, w9.c_str());

                        if(wt->genControllerType != "none")
                        {
                            fprintf(f, "%*s %*s %*s ", width, w10.c_str(), width, w11.c_str(), width, w12.c_str());
                        }

                        if(wt->pitchControllerType != "none")
                        {
                            fprintf(f, "%*s ", width, w13.c_str());
                        }
                    }
                    // uniform actuator disk model
                    else if((*farm->turbineModels[t]) == "uniformADM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str());
                    }
                    // actuator farm model
                    else if((*farm->turbineModels[t]) == "AFM")
                    {
                        fprintf(f, "%*s %*s %*s ", width, w0.c_str(), width, w3.c_str(), width, w4.c_str());
                    }
                    else if((*farm->turbineModels[t]) == "ALM")
                    {
                        fprintf(f, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s ", width, w0.c_str(), width, w1.c_str(), width, w2.c_str(), width, w3.c_str(), width, w4.c_str(), width, w5.c_str(), width, w6.c_str(), width, w7.c_str(), width, w8.c_str(), width, w9.c_str(), width, w16.c_str());

                        if(wt->genControllerType != "none")
                        {
                            fprintf(f, "%*s %*s %*s ", width, w10.c_str(), width, w11.c_str(), width, w12.c_str());
                        }

                        if(wt->pitchControllerType != "none")
                        {
                            fprintf(f, "%*s ", width, w13.c_str());
                        }
                    }

                    // yaw controller is not model specific
                    if(wt->yawControllerType != "none")
                    {
                        fprintf(f, "%*s %*s ", width, w14.c_str(), width, w15.c_str());
                    }

                    fprintf(f, "\n");

                    fclose(f);
                }

                if((*farm->turbineModels[t]) == "ALM")
                {
                    //create file to write the turbine radial point angle of attack values
                    sprintf(fileName, "%s/%s_aoa", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point relative velocity magnitude
                    sprintf(fileName, "%s/%s_relVelMag", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point rel velocity angle phi
                    sprintf(fileName, "%s/%s_phi", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in axial direction
                    sprintf(fileName, "%s/%s_uAxl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in radial direction
                    sprintf(fileName, "%s/%s_uRad", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point inflow velocity in tangential direction
                    sprintf(fileName, "%s/%s_uTan", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point Cl
                    sprintf(fileName, "%s/%s_Cl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point Cd
                    sprintf(fileName, "%s/%s_Cd", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point axial force
                    sprintf(fileName, "%s/%s_axialF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }

                    //create file to write the turbine radial point tangential force
                    sprintf(fileName, "%s/%s_tangF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = 20;

                        word w0 = "time [s]";

                        fprintf(f, "%s ", w0.c_str());

                        // radial mesh cell size size
                        PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

                        PetscReal rMag;
                        word w1 = "R";

                        for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
                        {
                            // compute radius
                            rMag = drval*ri;

                            // add hub radius (as if all was aligned)
                            rMag += wt->rHub;

                            fprintf(f, "%*s%ld (%lf)", width, w1.c_str(), ri, rMag);
                        }

                        fprintf(f, "\n");
                        fclose(f);
                    }
                }

            }
        }
    }

    // see if must write to file
    PetscInt  turbinesWrite = 0;
    word      intervalType  = farm->intervalType;
    PetscReal timeInterval  = farm->timeInterval;
    PetscReal timeStart     = farm->timeStart;

    // write every "timeInterval" seconds
    if
    (
        (intervalType == "adjustableTime") &&
        (
            (clock->time - timeStart) / timeInterval -
            std::floor
            (
                (clock->time - timeStart) / timeInterval
            ) < 1e-10
        )
    )
    {
        turbinesWrite = 1;
    }
    // write every "timeInterval" iterations
    else if
    (
        (clock->it > 0) &&
        (intervalType == "timeStep") &&
        (
            clock->it / timeInterval -
            std::floor
            (
                clock->it / timeInterval
            ) < 1e-10
        )
    )
    {
        turbinesWrite = 1;
    }

    if(turbinesWrite)
    {
        constants_ *constants = farm->access->constants;

        // loop over wind turbines
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->turbineControlled)
            {
                PetscReal Ar          = M_PI * pow(wt->rTip,2.0);

                // UADM/ADM
                PetscReal rtrAvgMagU  = 0.0;
                PetscReal rtrAvgUpMagU= 0.0;
                PetscReal rtrThrust   = 0.0;
                PetscReal aeroPwr     = 0.0;
                PetscReal ctInf       = 0.0;
                PetscReal ctLoc       = 0.0;
                PetscReal ctUp        = 0.0;

                // ADM
                PetscReal rtrTorque   = 0.0;
                PetscReal rtrOmega    = 0.0;

                // ALM
                PetscReal azimuth     = 0.0;

                // controls
                PetscReal genTorque   = 0.0;
                PetscReal genPwr      = 0.0;
                PetscReal genOmega    = 0.0;
                PetscReal collPitch   = 0.0;

                // do turbine level scatter first
                if((*farm->turbineModels[t]) == "ADM")
                {
                    MPI_Reduce(&(wt->adm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->adm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->adm.rtrTorque),  &rtrTorque,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->rtrOmega),       &rtrOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genOmega),    &genOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        genTorque    = genTorque  / wt->nProcsTrb / 1.0e3;
                        genPwr       = genPwr     / wt->nProcsTrb / 1.0e6;
                        genOmega     = genOmega   / wt->nProcsTrb / wt->rpm2RadSec;
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        collPitch    = collPitch  / wt->nProcsTrb;
                    }

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    rtrTorque    = rtrTorque  / wt->nProcsTrb / 1.0e3;
                    rtrOmega     = rtrOmega   / wt->nProcsTrb / wt->rpm2RadSec;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->adm.Uref,   2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    MPI_Reduce(&(wt->uadm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->uadm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->uadm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->uadm.Uref,  2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }
                else if((*farm->turbineModels[t]) == "AFM")
                {
                    MPI_Reduce(&(wt->afm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->afm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;
                }
                else if((*farm->turbineModels[t]) == "ALM")
                {
                    MPI_Reduce(&(wt->alm.rtrAvgMagU), &rtrAvgMagU,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.rtrThrust),  &rtrThrust,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->alm.aeroPwr),    &aeroPwr,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->alm.rtrTorque),  &rtrTorque,     1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    MPI_Reduce(&(wt->rtrOmega),       &rtrOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    MPI_Reduce(&(wt->alm.azimuth),    &azimuth,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genOmega),    &genOmega,      1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        genTorque    = genTorque  / wt->nProcsTrb / 1.0e3;
                        genPwr       = genPwr     / wt->nProcsTrb / 1.0e6;
                        genOmega     = genOmega   / wt->nProcsTrb / wt->rpm2RadSec;
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                        collPitch    = collPitch  / wt->nProcsTrb;
                    }

                    rtrAvgMagU   = rtrAvgMagU / wt->nProcsTrb;
                    rtrAvgUpMagU = wt->upPoints->Uref;
                    rtrThrust    = rtrThrust  / wt->nProcsTrb / 1.0e3;
                    aeroPwr      = aeroPwr    / wt->nProcsTrb / 1.0e6;

                    rtrTorque    = rtrTorque  / wt->nProcsTrb / 1.0e3;
                    rtrOmega     = rtrOmega   / wt->nProcsTrb / wt->rpm2RadSec;

                    azimuth      = azimuth    / wt->nProcsTrb;

                    ctInf        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(wt->alm.Uref,   2.0) * Ar );
                    ctLoc        = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgMagU,     2.0) * Ar );
                    ctUp         = (2.0 * rtrThrust * 1.0e3) / (constants->rho * pow(rtrAvgUpMagU,   2.0) * Ar );
                }

                // write to file
                if(rank == wt->writerRank)
                {
                    FILE *f;
                    char fileName[80];
                    sprintf(fileName, "%s/%s", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "a");

                    if(!f)
                    {
                        char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWrite",  error);
                    }
                    else
                    {
                        int width = -20;

                        // actuator disk model
                        if((*farm->turbineModels[t]) == "ADM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp, width, rtrTorque, width, rtrOmega);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "%*.4f %*.4f %*.4f ", width, genTorque, width, genPwr, width, genOmega);
                            }

                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "%*.4f ", width, collPitch*wt->rad2deg);
                            }
                        }
                        // uniform actuator disk model
                        else if((*farm->turbineModels[t]) == "uniformADM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp);
                        }
                        // actuator farm model
                        else if((*farm->turbineModels[t]) == "AFM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f ", width, clock->time, width, rtrThrust, width, aeroPwr);
                        }
                        else if((*farm->turbineModels[t]) == "ALM")
                        {
                            fprintf(f, "%*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f %*.4f ", width, clock->time, width, rtrAvgMagU, width, rtrAvgUpMagU, width, rtrThrust, width, aeroPwr, width, ctInf, width, ctLoc, width, ctUp, width, rtrTorque, width, rtrOmega, width, azimuth);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "%*.4f %*.4f %*.4f ", width, genTorque, width, genPwr, width, genOmega);
                            }

                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "%*.4f ", width, collPitch*wt->rad2deg);
                            }
                        }

                        // yaw controller is not model specific
                        if(wt->yawControllerType != "none")
                        {
                            fprintf(f, "%*.4f %*.4f ", width, wt->flowAngle, width, wt->yawAngle);
                        }

                        fprintf(f, "\n");

                        fclose(f);
                    }

                    if((*farm->turbineModels[t]) == "ALM")
                    {
                        sprintf(fileName, "%s/%s_aoa", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sumAlpha = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sumAlpha += wt->alm.alpha[nAz*ri + ai];
                                }

                                sumAlpha/=nAz;

                                fprintf(f, "%*.4f", width, sumAlpha);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_relVelMag", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += nMag(wt->alm.U[nAz*ri + ai]);
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_phi", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai, p;
                            Cmpnts point_p, r_p;
                            Cmpnts xb_hat, yb_hat, zb_hat;
                            PetscReal ub_x, ub_y, phi;

                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    p =  nAz*ri + ai;

                                    // save this point locally for speed
                                    point_p = wt->alm.points[p];

                                    // this point position from COR
                                    r_p  = nSub(point_p, wt->rotCenter);

                                    if(wt->rotDir == "cw")
                                    {
                                        zb_hat = nUnit(r_p);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }
                                    else if(wt->rotDir == "ccw")
                                    {
                                        zb_hat = nUnit(r_p);
                                                 mScale(-1.0, zb_hat);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }

                                    ub_x = nDot(wt->alm.U[p], xb_hat);
                                    ub_y = nDot(wt->alm.U[p], yb_hat);

                                    phi = wt->rad2deg * std::atan2(ub_x, ub_y);
                                    sum += phi;
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        FILE *f1, *f2, *f3;
                        char fileName1[80], fileName2[80], fileName3[80];
                        sprintf(fileName1, "%s/%s_uAxl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f1 = fopen(fileName1, "a");
                        sprintf(fileName2, "%s/%s_uRad", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f2 = fopen(fileName2, "a");
                        sprintf(fileName3, "%s/%s_uTan", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f3 = fopen(fileName3, "a");
                        if(!f1 || !f2 || !f3)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName1);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f1, "%.4f", clock->time);
                            fprintf(f2, "%.4f", clock->time);
                            fprintf(f3, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai, p;
                            Cmpnts point_p, r_p;
                            Cmpnts xb_hat, yb_hat, zb_hat;
                            PetscReal ub_x, ub_y, ub_z, phi;

                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sumAxl = 0.0, sumTan = 0.0, sumRad = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    p =  nAz*ri + ai;

                                    // save this point locally for speed
                                    point_p = wt->alm.points[p];

                                    // this point position from COR
                                    r_p  = nSub(point_p, wt->rotCenter);

                                    if(wt->rotDir == "cw")
                                    {
                                        zb_hat = nUnit(r_p);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }
                                    else if(wt->rotDir == "ccw")
                                    {
                                        zb_hat = nUnit(r_p);
                                                 mScale(-1.0, zb_hat);
                                        xb_hat = nScale(-1.0, farm->wt[t]->rtrAxis);
                                        yb_hat = nCross(zb_hat, xb_hat);
                                    }

                                    ub_x = nDot(wt->alm.gWind[p], xb_hat);
                                    ub_y = nDot(wt->alm.gWind[p], yb_hat);
                                    ub_z = nDot(wt->alm.gWind[p], zb_hat);

                                    sumAxl += ub_x;
                                    sumTan += ub_y;
                                    sumRad += ub_z;
                                }

                                sumAxl/=nAz;
                                sumTan/=nAz;
                                sumRad/=nAz;

                                fprintf(f1, "%*.4f", width, sumAxl);
                                fprintf(f2, "%*.4f", width, sumRad);
                                fprintf(f3, "%*.4f", width, sumTan);
                            }

                            fprintf(f1, "\n");
                            fprintf(f2, "\n");
                            fprintf(f3, "\n");

                            fclose(f1);
                            fclose(f2);
                            fclose(f3);
                        }

                        sprintf(fileName, "%s/%s_Cl", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.Cl[nAz*ri + ai];
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_Cd", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.Cd[nAz*ri + ai];
                                }

                                sum/=nAz;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_axialF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;
                                PetscReal sumDr = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.axialF[nAz*ri + ai];
                                    sumDr += wt->alm.dr[nAz*ri + ai];
                                }

                                sum/=nAz;
                                sumDr/=nAz;

                                //to compute the force per unit length divide by element radial length
                                sum/=sumDr;

                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }

                        sprintf(fileName, "%s/%s_tangF", turbineFolderTimeName.c_str(), (*farm->turbineIds[t]).c_str());
                        f = fopen(fileName, "a");

                        if(!f)
                        {
                            char error[512];
                            sprintf(error, "cannot open file %s\n", fileName);
                            fatalErrorInFunction("windTurbinesWrite",  error);
                        }
                        else
                        {
                            int width = 20;

                            fprintf(f, "%.4f", clock->time);

                            PetscInt nRad = wt->alm.nRadial;
                            PetscInt nAz = wt->alm.nAzimuth;

                            PetscInt ri, ai;
                            //average the aoa over the azimuthal points
                            for(ri=0; ri<nRad; ri++)
                            {
                                PetscReal sum = 0.0;
                                PetscReal sumDr = 0.0;

                                for(ai=0; ai<nAz; ai++)
                                {
                                    sum += wt->alm.tangtF[nAz*ri + ai];
                                    sumDr += wt->alm.dr[nAz*ri + ai];
                                }

                                sum/=nAz;
                                sumDr/=nAz;

                                //to compute the force per unit length divide by element radial length
                                sum/=sumDr;
                                                                
                                fprintf(f, "%*.4f", width, sum);
                            }

                            fprintf(f, "\n");

                            fclose(f);
                        }
                    }

                }
            }
        }
    }

    // write wind farm mesh according to runTimeWrite option
    if(io->runTimeWrite)
    {
        if(farm->writeNumber == 0)
        {
            // write the wind farm mesh
            writeFarmTwrMesh(farm);
        }

        // write AD mesh
        writeFarmADMesh(farm);

        // write AL mesh
        writeFarmALMesh(farm);

        // increase write number counter
        farm->writeNumber++;
    }

    // write checkpoint file
    windTurbinesWriteCheckpoint(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windTurbinesWriteCheckpoint(farm_ *farm)
{
    // writes in fields directory which has been already created at sim start

    clock_      *clock = farm->access->clock;
    io_         *io    = farm->access->io;
    mesh_       *mesh  = farm->access->mesh;

    word        timeName;
    word        path, turbineFolderName;

    PetscInt    t;

    PetscMPIInt rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    // set time folder name
    timeName          = getTimeName(clock);
    path              = "./fields/" + mesh->meshName + "/turbines/" + timeName;
    turbineFolderName = "./fields/" + mesh->meshName + "/turbines";

    // create/initialize fields/turbines directory (at simulation start only)
    if
    (
        clock->it == clock->itStart && !rank
    )
    {
        errno = 0;
        PetscInt dirRes = mkdir(turbineFolderName.c_str(), 0777);
        if(dirRes != 0 && errno != EEXIST)
        {
           char error[512];
            sprintf(error, "could not create fields/turbines directory");
            fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
        }

        // if directory already exist remove everything inside except the start time (safe)
        if(errno == EEXIST)
        {
            word startTimeName = getStartTimeName(clock);
            remove_subdirs_except(mesh->MESH_COMM, turbineFolderName.c_str(), startTimeName.c_str());
        }
    }

    // test if must write at this time step
    if(io->runTimeWrite)
    {
        // creates time folder
        if(!rank)
        {
            PetscInt dirRes = mkdir(path.c_str(), 0777);
            if(dirRes != 0 && errno != EEXIST)
            {
               char error[512];
                sprintf(error, "could not create %s directory\n", path.c_str());
                fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
            }
        }

        // ensure folder is there for every processor
        MPI_Barrier(mesh->MESH_COMM);

        // write the checkpoint for each wind turbine
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            if(wt->turbineControlled)
            {
                // first scatter the necessary information on the master node of this TRB_COMM
                PetscReal rtrOmega     = 0.0;
                PetscReal rtrOmegaFilt = 0.0;
                PetscReal genTorque    = 0.0;
                PetscReal genPwr       = 0.0;
                PetscReal collPitch    = 0.0;
                PetscReal errPID       = 0.0;
                PetscReal intErrPID    = 0.0;
                PetscReal azimuth      = 0.0;

                if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                {
                    MPI_Reduce(&(wt->rtrOmega), &rtrOmega,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);

                    if(wt->genControllerType != "none")
                    {
                        MPI_Reduce(&(wt->rtrOmegaFilt),&rtrOmegaFilt,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genTorque),   &genTorque,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->genPwr),      &genPwr,          1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }
                    if(wt->pitchControllerType != "none")
                    {
                        MPI_Reduce(&(wt->collPitch),   &collPitch,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->errPID),      &errPID,       1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                        MPI_Reduce(&(wt->intErrPID),   &intErrPID,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }

                    if((*farm->turbineModels[t]) == "ALM")
                    {
                        MPI_Reduce(&(wt->alm.azimuth), &azimuth,    1, MPIU_REAL, MPIU_SUM, 0, wt->TRB_COMM);
                    }
                }

                // now only the master of TRB_COMM writes
                if(rank == wt->writerRank)
                {
                    FILE *f;
                    char fileName[80];
                    sprintf(fileName, "%s/%s", path.c_str(), (*farm->turbineIds[t]).c_str());
                    f = fopen(fileName, "w");

                    if(!f)
                    {
                       char error[512];
                        sprintf(error, "cannot open file %s\n", fileName);
                        fatalErrorInFunction("windTurbinesWriteCheckpoint",  error);
                    }
                    else
                    {
                        // turbine level properties
                        fprintf(f, "turbineLevelProperties\n{\n");
                        fprintf(f, "    rtrDir             (%lf %lf %lf)\n", wt->rtrDir.x, wt->rtrDir.y, wt->rtrDir.z);
                        fprintf(f, "    rtrAxis            (%lf %lf %lf)\n", wt->rtrAxis.x, wt->rtrAxis.y, wt->rtrAxis.z);

                        if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
                        {
                            fprintf(f, "    omega_hat          (%lf %lf %lf)\n", wt->omega_hat.x, wt->omega_hat.y, wt->omega_hat.z);
                            fprintf(f, "    rtrOmega           %lf\n",         rtrOmega / wt->nProcsTrb);

                            if(wt->genControllerType != "none")
                            {
                                fprintf(f, "    rtrOmegaFilt       %lf\n",         rtrOmegaFilt / wt->nProcsTrb);
                                fprintf(f, "    genTorque          %lf\n",         genTorque    / wt->nProcsTrb);
                                fprintf(f, "    genPwr             %lf\n",         genPwr       / wt->nProcsTrb);
                            }
                            if(wt->pitchControllerType != "none")
                            {
                                fprintf(f, "    collPitch          %lf\n",         collPitch    / wt->nProcsTrb);
                                fprintf(f, "    errPID             %lf\n",         errPID       / wt->nProcsTrb);
                                fprintf(f, "    intErrPID          %lf\n",         intErrPID    / wt->nProcsTrb);
                            }

                            if((*farm->turbineModels[t]) == "ALM")
                            {
                                fprintf(f, "    azimuth          %lf\n",           azimuth    / wt->nProcsTrb);
                            }

                        }
                        if(wt->yawControllerType != "none")
                        {
                            fprintf(f, "    yawAngle           %lf\n",         wt->yawAngle);
                            fprintf(f, "    flowAngle          %lf\n",         wt->flowAngle);
                            fprintf(f, "    yawError           %lf\n",         wt->yawError);
                        }

                        fprintf(f, "}\n\n");
                    }

                    fclose(f);
                }
            }
        }

        // ensure all files are closed
        MPI_Barrier(mesh->MESH_COMM);

        // remove old checkpoint files except last one after all files are written (safe)
        if(!rank && mesh->access->io->purgeWrite) remove_subdirs_except(mesh->MESH_COMM, turbineFolderName.c_str(), timeName);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode windTurbinesReadCheckpoint(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    clock_      *clock = farm->access->clock;

    word        timeName;
    PetscInt         t, p;

    PetscReal      initialYaw, initialYawError;

    timeName = "./fields/" + mesh->meshName + "/turbines/" + getTimeName(clock);

    // check if this folder exists
    if(dir_exist(timeName.c_str()))
    {
        PetscPrintf(mesh->MESH_COMM, "   found checkpoint data %s, restoring configuration...\n",timeName.c_str());

        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // each processor reads data and stores
            char dictName[80];
            sprintf(dictName, "%s/%s", timeName.c_str(), (*farm->turbineIds[t]).c_str());

            // save initial yaw before overwriting it
            initialYaw = std::atan2(-1.0 * wt->rtrDir.y, -1.0 * wt->rtrDir.x) * wt->rad2deg;

            readSubDictVector(dictName, "turbineLevelProperties", "rtrDir", &(wt->rtrDir));
            readSubDictVector(dictName, "turbineLevelProperties", "rtrAxis", &(wt->rtrAxis));

            if((*farm->turbineModels[t]) != "uniformADM" && (*farm->turbineModels[t]) != "AFM")
            {
                readSubDictVector(dictName, "turbineLevelProperties", "omega_hat", &(wt->omega_hat));
                readSubDictDouble(dictName, "turbineLevelProperties", "rtrOmega", &(wt->rtrOmega));

                if(wt->genControllerType != "none")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "rtrOmegaFilt", &(wt->rtrOmegaFilt));
                    readSubDictDouble(dictName, "turbineLevelProperties", "genTorque", &(wt->genTorque));
                    readSubDictDouble(dictName, "turbineLevelProperties", "genPwr", &(wt->genPwr));
                }
                if(wt->pitchControllerType != "none")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "collPitch", &(wt->collPitch));
                    readSubDictDouble(dictName, "turbineLevelProperties", "errPID", &(wt->errPID));
                    readSubDictDouble(dictName, "turbineLevelProperties", "intErrPID", &(wt->intErrPID));
                }
                if((*farm->turbineModels[t]) == "ALM")
                {
                    readSubDictDouble(dictName, "turbineLevelProperties", "azimuth", &(wt->alm.azimuth));
                }
            }
            if(wt->yawControllerType != "none")
            {
                readSubDictDouble(dictName, "turbineLevelProperties", "yawAngle", &(wt->yawAngle));
                readSubDictDouble(dictName, "turbineLevelProperties", "flowAngle", &(wt->flowAngle));
                readSubDictDouble(dictName, "turbineLevelProperties", "yawError", &(wt->yawError));
            }
        }

        // yaw turbines
        for(t=0; t<farm->size; t++)
        {
            windTurbine *wt = farm->wt[t];

            // test if this turbine has yaw control
            if(wt->yawControllerType != "none")
            {
                PetscReal yawError = wt->yawAngle - initialYaw;
                PetscReal rotationSign;

                if(yawError > 0.0)       rotationSign =  1.0;
                else if(yawError < -0.0) rotationSign = -1.0;
                else                     rotationSign =  0.0;

                // compute rotation angle in radiants
                PetscReal angle = rotationSign * yawError * wt->deg2rad;

                if(angle != 0.0)
                {
                    // actuator disk model
                    if((*farm->turbineModels[t]) == "ADM")
                    {
                        // number of points in the AD mesh
                        PetscInt npts_t = wt->adm.nPoints;

                        // loop over the AD mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(wt->adm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->adm.points[p], angle);
                            mSum(wt->adm.points[p], wt->rotCenter);
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
                            mSub(wt->uadm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->uadm.points[p], angle);
                            mSum(wt->uadm.points[p], wt->rotCenter);
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
                            mSub(wt->alm.points[p], wt->rotCenter);
                            mRot(wt->twrDir, wt->alm.points[p], angle);
                            mSum(wt->alm.points[p], wt->rotCenter);
                        }
                    }
                    else if((*farm->turbineModels[t]) == "AFM")
                    {
                        // nothing to rotate
                    }

                    // rotate up-sampling points
                    {
                        upSampling *upPoints = farm->wt[t]->upPoints;

                        // rotate center
                        mSub(upPoints->center, wt->rotCenter);
                        mRot(wt->twrDir, upPoints->center, angle);
                        mSum(upPoints->center, wt->rotCenter);

                        // number of points in the sample mesh
                        PetscInt npts_t = upPoints->nPoints;

                        // loop over the sample mesh points
                        for(p=0; p<npts_t; p++)
                        {
                            mSub(upPoints->points[p], wt->rotCenter);
                            mRot(wt->twrDir, upPoints->points[p], angle);
                            mSum(upPoints->points[p], wt->rotCenter);
                        }
                    }
                }
            }

            // rotate blades
            if((*farm->turbineModels[t]) == "ALM")
            {
                rotateBlades(wt, wt->alm.azimuth * wt->deg2rad);
            }
        }
    }
    else
    {
        PetscPrintf(mesh->MESH_COMM, "   no checkpoint data found, setting reference configuration...\n");
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initADM(windTurbine *wt, Cmpnts &base, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(ADM), &(wt->adm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AD parameters
    readDictInt(descrFile.c_str(), "nRadPts", &(wt->adm.nRadial));
    readDictInt(descrFile.c_str(), "nAziPts", &(wt->adm.nAzimuth));
    readDictDouble(descrFile.c_str(), "Uref",    &(wt->adm.Uref));

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->adm.dbg));

    // set total numer of points in the mesh
    wt->adm.nPoints = wt->adm.nRadial * wt->adm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->adm.rtrThrust = 0.0;
    wt->adm.rtrTorque = 0.0;
    wt->adm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->adm.rtrAvgMagU = 0.0;

    // build the AD mesh

    // allocate memory for the ADM parameters
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.points));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.dr));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.chord));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.twist));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.solidity));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscInt*), &(wt->adm.foilIds));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal*), &(wt->adm.iw));

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->adm.nPoints*sizeof(PetscInt), &(wt->adm.thisPtControlled));
    PetscMalloc(wt->adm.nPoints*sizeof(cellIds),&(wt->adm.closestCells));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.Cd));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.Cl));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.alpha));
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.U));
    PetscMalloc(wt->adm.nPoints*sizeof(Cmpnts), &(wt->adm.B));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.axialF));
    PetscMalloc(wt->adm.nPoints*sizeof(PetscReal), &(wt->adm.tangtF));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);

    // set rotor center point
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);
    wt->rotCenter = center;

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);

    // set rtrAxis and omega_hat turbine param. here
    {
        // set rotor axis (up-tilted, from nacell back to cone)
        wt->rtrAxis = xr_hat;

        // set the rotor rotation unit vector
        if(wt->rotDir == "cw")
        {
            wt->omega_hat = nScale(-1.0, xr_hat);
        }
        else if(wt->rotDir == "ccw")
        {
            wt->omega_hat = nSet(xr_hat);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown rotationDir, avilable options are 'cw' or 'ccw'\n");
            fatalErrorInFunction("initADM",  error);
        }
    }

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->adm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->adm.nAzimuth;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr;

    for(PetscInt ri=0; ri<wt->adm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->adm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

		// compute radius
		PetscReal  rMag = nMag(rvec);

        // add the initial hub radius to rvec (not aligned if precone != 0)
        mSum(rvec, rHub);

		// add hub radius (as if all was aligned)
		rMag += nMag(rHub);

        // interpolate blade propertes (only depend on r)
        PetscReal  weights[2];
        PetscInt   labels[2];

        findInterpolationWeigths(weights, labels, wt->blade.radius, wt->blade.size, rMag);

        for(PetscInt ai=0; ai<wt->adm.nAzimuth; ai++)
        {
            // set interpolation variables
            PetscReal w1 = weights[0]; PetscInt l1 = labels[0];
            PetscReal w2 = weights[1]; PetscInt l2 = labels[1];

            // allocate memory for the 2 closest foil ids
            PetscMalloc(2*sizeof(PetscInt), &(wt->adm.foilIds[pi]));

            // allocate memory for the 2 interpolation weights
            PetscMalloc(2*sizeof(PetscReal), &(wt->adm.iw[pi]));

            // set airfoil chord
            wt->adm.chord[pi] = w1*wt->blade.chord[l1] + w2*wt->blade.chord[l2];

            // set airfoil twist
            wt->adm.twist[pi] = w1*wt->blade.twist[l1] + w2*wt->blade.twist[l2];

            // set airfoil ids
            wt->adm.foilIds[pi][0] = wt->blade.foilIds[l1];
            wt->adm.foilIds[pi][1] = wt->blade.foilIds[l2];

            // set interpolation weights
            wt->adm.iw[pi][0] = w1;
            wt->adm.iw[pi][1] = w2;

            // set the rotor solidity (only depends on radius)
            wt->adm.solidity[pi] = (PetscReal)wt->nBlades / (PetscReal)wt->adm.nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(xr_hat, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, center);

            // set the point value
            mSet(wt->adm.points[pi], point);

            // set the dr value (uniform for now)
            wt->adm.dr[pi] = dr;

            pi++;
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initUADM(windTurbine *wt, Cmpnts &base, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(UADM), &(wt->uadm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AD parameters
    readDictInt(descrFile.c_str(),    "nRadPts",    &(wt->uadm.nRadial));
    readDictInt(descrFile.c_str(),    "nAziPts",    &(wt->uadm.nAzimuth));
    readDictDouble(descrFile.c_str(), "Ct",         &(wt->uadm.Ct));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->uadm.sampleType));
    readDictDouble(descrFile.c_str(), "Uref",       &(wt->uadm.Uref));

    // check sample type
    if
    (
        wt->uadm.sampleType != "rotorUpstream" &&  // samples velocity on 2.5D upstream disk (needs Ct)
        wt->uadm.sampleType != "givenVelocity" &&  // uses input velocity (needs Ct, good for isolated turbine)
        wt->uadm.sampleType != "rotorDisk"         // samples velocity on rotor disk (needs CtPrime)
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are givenVelocity, rotorUpstream or rotorDisk");
        fatalErrorInFunction("initUADM",  error);
    }

    // compute axial induction factor
    if(wt->uadm.sampleType != "rotorDisk")
    {
        // check that Ct does not make induction complex or negative
        if(wt->uadm.Ct <= 0.0 || wt->uadm.Ct >= 1.0)
        {
            char error[512];
            sprintf(error, "provided thrust coefficient ouside of bounds ([0 1] excluded). Change or switch to sampleType = rotorDisk");
            fatalErrorInFunction("initUADM",  error);
        }

        // Ct
        wt->uadm.axiInd = (1.0 - sqrt(1 - wt->uadm.Ct)) / 2.0;
    }
    else
    {
        // CtPrime
        wt->uadm.axiInd = 0.0;
    }

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->uadm.dbg));

    // set total numer of points in the mesh
    wt->uadm.nPoints = wt->uadm.nRadial * wt->uadm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->uadm.rtrThrust = 0.0;
    wt->uadm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->uadm.rtrAvgMagU = 0.0;

    // build the AD mesh

    // allocate memory for the ADM parameters
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.points));
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscReal), &(wt->uadm.dA));

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscInt), &(wt->uadm.thisPtControlled));
    PetscMalloc(wt->uadm.nPoints*sizeof(cellIds),&(wt->uadm.closestCells));
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.U));
    PetscMalloc(wt->uadm.nPoints*sizeof(Cmpnts), &(wt->uadm.B));
    PetscMalloc(wt->uadm.nPoints*sizeof(PetscReal), &(wt->uadm.axialF));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);

    // set rotor center point
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);
    wt->rotCenter = center;

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);

    // set rotor axis (up-tilted, from nacell back to cone)
    wt->rtrAxis = xr_hat;

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->uadm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->uadm.nAzimuth;

    // area of the circular crown from hub to tip radius
    PetscReal bladeSweptAreaNoHub = M_PI * (wt->rTip*wt->rTip - wt->rHub*wt->rHub);

    // hub area
    PetscReal hubFrontArea = M_PI * wt->rHub * wt->rHub;

    // hub area distribution coeff
    PetscReal hCoeff = hubFrontArea / bladeSweptAreaNoHub;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr, crownArea, r = wt->rHub;

    for(PetscInt ri=0; ri<wt->uadm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->uadm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // compute r+dr crown area
        crownArea = M_PI * ((r+dr)*(r+dr) - r*r);

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

        // add the initial hub radius
        mSum(rvec, rHub);

        for(PetscInt ai=0; ai<wt->uadm.nAzimuth; ai++)
        {
            // set cell area
            wt->uadm.dA[pi] = (1.0 + hCoeff) * crownArea / wt->uadm.nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(xr_hat, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, center);

            // set the point value
            mSet(wt->uadm.points[pi], point);

            pi++;
        }

        r += dr;
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode initSamplePoints(windTurbine *wt, Cmpnts &base, const word meshName)
{
    // allocate memory for this turbine sample points struct
    wt->upPoints = (upSampling*)malloc(sizeof(upSampling));
    upSampling* upPoints = wt->upPoints;

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // discretization is not very fine since only area-weighted
    // average must be performed
    PetscInt nAzimuth = 12;
    PetscInt nRadial  = 12;

    upPoints->nPoints = nAzimuth * nRadial;

    // allocate objects memory
    PetscMalloc(upPoints->nPoints*sizeof(Cmpnts),    &(upPoints->points));
    PetscMalloc(upPoints->nPoints*sizeof(cellIds),   &(upPoints->closestCells));
    PetscMalloc(upPoints->nPoints*sizeof(PetscInt),  &(upPoints->thisPtControlled));
    PetscMalloc(upPoints->nPoints*sizeof(PetscReal), &(upPoints->dA));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);

    // set rotor center point
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    // (do not rotate with uptilt)
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);

    // set the upstream sample distance
    Cmpnts upDist = nScale(2.5*2.0*wt->rTip, xr_hat);

    // set the sampling rotor center
    mSum(center, upDist);
    upPoints->center = center;

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / nAzimuth;

    // area of the circular crown from hub to tip radius
    PetscReal bladeSweptAreaNoHub = M_PI * (wt->rTip*wt->rTip - wt->rHub*wt->rHub);

    // hub area
    PetscReal hubFrontArea = M_PI * wt->rHub * wt->rHub;

    // hub area distribution coeff
    PetscReal hCoeff = hubFrontArea / bladeSweptAreaNoHub;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr, crownArea, r = wt->rHub;

    for(PetscInt ri=0; ri<nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==nRadial-1) dr = drval / 2;
        else dr = drval;

        // compute r+dr crown area
        crownArea = M_PI * ((r+dr)*(r+dr) - r*r);

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, yr_hat);

        // add the initial hub radius
        mSum(rvec, rHub);

        for(PetscInt ai=0; ai<nAzimuth; ai++)
        {
            // set cell area
            upPoints->dA[pi] = (1.0 + hCoeff) * crownArea / nAzimuth;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(xr_hat, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, center);

            // set the point value
            mSet(upPoints->points[pi], point);

            pi++;
        }

        r += dr;
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode initALM(windTurbine *wt, Cmpnts &base, const word meshName)
{
    // allocate memory for the ADM
    PetscMalloc(sizeof(ALM), &(wt->alm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    // read from file AL parameters
    readDictWord(descrFile.c_str(),   "projection", &(wt->alm.projectionType));
    readDictInt(descrFile.c_str(),    "nRadPts", &(wt->alm.nRadial));
    readDictDouble(descrFile.c_str(), "Uref",    &(wt->alm.Uref));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->alm.sampleType));

        // check sample type
    if
    (
        wt->alm.sampleType != "rotorDisk" &&
        wt->alm.sampleType != "integral"
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are rotorDisk or integral");
        fatalErrorInFunction("initALM",  error);
    }

    // check projection type
    if(wt->alm.projectionType!="isotropic" && wt->alm.projectionType!="anisotropic")
    {
        char error[512];
        sprintf(error, "unknown ALM projection, available possibilities are:\n    1. isotropic\n    2. anisotropic\n");
        fatalErrorInFunction("initALM",  error);
    }

    //read projection radius
    if(wt->alm.projectionType=="anisotropic")
    {
        readDictDouble(descrFile.c_str(), "epsilonFactor_x",     &(wt->eps_x));
        readDictDouble(descrFile.c_str(), "epsilonFactor_y",     &(wt->eps_y));
        readDictDouble(descrFile.c_str(), "epsilonFactor_z",     &(wt->eps_z));
    }
    else
    {
        readDictDouble(descrFile.c_str(), "epsilon",     &(wt->eps));
    }

    wt->alm.nAzimuth = wt->nBlades;

    // debug switch
    readDictInt(descrFile.c_str(), "debug", &(wt->alm.dbg));

    // set total numer of points in the mesh
    wt->alm.nPoints = wt->alm.nRadial * wt->alm.nAzimuth;

    // set rotor torque and power to zero (will remain zero in the processors
    // that do not control the turbine for parallel scatter/gather)
    wt->alm.rtrThrust = 0.0;
    wt->alm.rtrTorque = 0.0;
    wt->alm.aeroPwr   = 0.0;

    // set average rotor velocity mag to zero
    wt->alm.rtrAvgMagU = 0.0;

    // set initial azimuth
    wt->alm.azimuth = 0.0;

    // build the AL mesh

    // allocate memory for the ALM parameters
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.points));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.dr));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.chord));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.twist));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.solidity));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscInt*), &(wt->alm.foilIds));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal*), &(wt->alm.iw));

    if(wt->alm.projectionType=="anisotropic")
    {
        PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.thick));
    }

    // allocate memory for the variables used during the simulation
    PetscMalloc(wt->alm.nPoints*sizeof(PetscInt), &(wt->alm.thisPtControlled));
    PetscMalloc(wt->alm.nPoints*sizeof(cellIds),&(wt->alm.closestCells));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.Cd));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.Cl));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.alpha));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.U));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.gWind));
    PetscMalloc(wt->alm.nPoints*sizeof(Cmpnts), &(wt->alm.B));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.axialF));
    PetscMalloc(wt->alm.nPoints*sizeof(PetscReal), &(wt->alm.tangtF));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);

    // set rotor center point
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);
    wt->rotCenter = center;

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

	// preconed r_hat
	Cmpnts r_hat = nSet(yr_hat);
	mRot(zr_hat, r_hat, wt->precone*wt->deg2rad);


    // set rtrAxis and omega_hat turbine param. here
    {
        // set rotor axis (up-tilted, from nacell back to cone)
        wt->rtrAxis = xr_hat;

        // set the rotor rotation unit vector
        if(wt->rotDir == "cw")
        {
            wt->omega_hat = nScale(-1.0, xr_hat);
        }
        else if(wt->rotDir == "ccw")
        {
            wt->omega_hat = nSet(xr_hat);
        }
        else
        {
            char error[512];
            sprintf(error, "unknown rotationDir, avilable options are 'cw' or 'ccw'\n");
            fatalErrorInFunction("initALM",  error);
        }
    }

    // set the hub radius vector
    Cmpnts rHub = nScale(wt->rHub, yr_hat);

    // radial mesh cell size size
    PetscReal drval = (wt->rTip - wt->rHub) / (wt->alm.nRadial - 1);

    // delta angle in radiants
    PetscReal daval = 2 * M_PI / wt->alm.nAzimuth;

    // points counter
    PetscInt    pi = 0;

    // varying variables
    PetscReal dr;

    for(PetscInt ri=0; ri<wt->alm.nRadial; ri++)
    {
        // delta radius (beware start and end points)
        if(ri==0 || ri==wt->alm.nRadial-1) dr = drval / 2;
        else dr = drval;

        // this station vector radius
        Cmpnts rvec = nScale(drval*ri, r_hat);

		// compute radius
		PetscReal  rMag = nMag(rvec);

        // add the initial hub radius to rvec (not aligned if precone != 0)
        mSum(rvec, rHub);

		// add hub radius (as if all was aligned)
		rMag += nMag(rHub);

        // interpolate blade propertes (only depend on r)
        PetscReal  weights[2];
        PetscInt   labels[2];


        findInterpolationWeigths(weights, labels, wt->blade.radius, wt->blade.size, rMag);

        for(PetscInt ai=0; ai<wt->alm.nAzimuth; ai++)
        {
            // set interpolation variables
            PetscReal w1 = weights[0]; PetscInt l1 = labels[0];
            PetscReal w2 = weights[1]; PetscInt l2 = labels[1];

            // allocate memory for the 2 closest foil ids
            PetscMalloc(2*sizeof(PetscInt), &(wt->alm.foilIds[pi]));

            // allocate memory for the 2 interpolation weights
            PetscMalloc(2*sizeof(PetscReal), &(wt->alm.iw[pi]));

            // set airfoil chord
            wt->alm.chord[pi] = w1*wt->blade.chord[l1] + w2*wt->blade.chord[l2];

            // set airfoil twist
            wt->alm.twist[pi] = w1*wt->blade.twist[l1] + w2*wt->blade.twist[l2];

            // set airfoil thickness
            if(wt->alm.projectionType=="anisotropic")
            {
                wt->alm.thick[pi] = w1*wt->blade.thick[l1] + w2*wt->blade.thick[l2];
            }

            // set airfoil ids
            wt->alm.foilIds[pi][0] = wt->blade.foilIds[l1];
            wt->alm.foilIds[pi][1] = wt->blade.foilIds[l2];

            // set interpolation weights
            wt->alm.iw[pi][0] = w1;
            wt->alm.iw[pi][1] = w2;

            // set the rotor solidity (it is one for the ALM)
            wt->alm.solidity[pi] = 1.0;

            // new mesh point
            Cmpnts point = nSet(rvec);

            // rotate the point
            mRot(xr_hat, point, daval*ai);

            // add the rotor center vector from origin
            mSum(point, center);

            // set the point value
            mSet(wt->alm.points[pi], point);

            // set the dr value (uniform for now)
            wt->alm.dr[pi] = dr;

            pi++;
        }
    }

    return(0);
};

//***************************************************************************************************************//

PetscErrorCode initAFM(windTurbine *wt, Cmpnts &base, const word meshName)
{
    // allocate memory for the AFM
    PetscMalloc(sizeof(AFM), &(wt->afm));

    // read necessary properties from file
    word descrFile = "./turbines/" + meshName + "/" + wt->type;

    readDictDouble(descrFile.c_str(), "Ct",         &(wt->afm.Ct));
    readDictDouble(descrFile.c_str(), "Uref",       &(wt->afm.Uref));
    readDictWord(descrFile.c_str(),   "sampleType", &(wt->afm.sampleType));

    // check sample type
    if
    (
        wt->afm.sampleType != "rotorDisk" &&
        wt->afm.sampleType != "momentumTheory" &&
        wt->afm.sampleType != "integral"
    )
    {
        char error[512];
        sprintf(error, "unknown velocity sampling type. Available types are momentumTheory, rotorDisk or integral");
        fatalErrorInFunction("initAFM",  error);
    }

    // check input
    if(wt->afm.sampleType == "momentumTheory" && wt->afm.Ct > 0.999)
    {
        char error[512];
        sprintf(error, "Ct coefficient must be less than one if sampling is of momentumTheory type");
        fatalErrorInFunction("initAFM",  error);
    }

    // read sampling frequency
    PetscReal windSpeedFilterPrd;
    readDictDouble(descrFile.c_str(), "windSpeedFilterPrd", &windSpeedFilterPrd);
    wt->afm.rtrUFilterFreq = 1.0 / windSpeedFilterPrd;

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);
                    mSum(tower, base);

    // set rotor center point
    Cmpnts overH  = nScale(wt->ovrHang, wt->rtrDir);
    Cmpnts center = nSum(tower, overH);

    wt->afm.point  = center;
    wt->rotCenter  = center;

    // set the rotor reference frame
    // x from nacelle back to cone,
    // y on the rotor at zero azimuth,
    // z as the right hand rule
    Cmpnts xr_hat = nUnit(wt->rtrDir);
    Cmpnts yr_hat = nUnit(wt->twrDir);
    Cmpnts zr_hat = nCross(xr_hat, yr_hat);
    mRot(zr_hat, xr_hat, wt->upTilt*wt->deg2rad);  // rotate xr_hat with uptilt
    mRot(zr_hat, yr_hat, wt->upTilt*wt->deg2rad);  // rotate yr_hat with uptilt

    // set rotor axis (up-tilted, from nacell back to cone)
    wt->rtrAxis = xr_hat;

    wt->afm.searchDone = 0;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initTwrModel(windTurbine *wt, Cmpnts &base)
{
    // set tower thrust to zero
    wt->twr.twrThrust = 0.0;

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->twr.prjNSigma = sqrt(std::log(1.0/0.001));

    // allocate memory
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.points));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscReal), &(wt->twr.dA));
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.U));
    PetscMalloc(wt->twr.nPoints*sizeof(Cmpnts), &(wt->twr.B));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscReal), &(wt->twr.tangF));
    PetscMalloc(wt->twr.nPoints*sizeof(PetscInt), &(wt->twr.thisPtControlled));
    PetscMalloc(wt->twr.nPoints*sizeof(cellIds), &(wt->twr.closestCells));

    // set tower top point
    Cmpnts tower  = nScale(wt->hTwr, wt->twrDir);

    // linear element distance
    PetscReal dhval  = nMag(tower) / (wt->twr.nPoints - 1);

    // rastremation coefficient
    PetscReal rastCoeff = (wt->twr.rTop - wt->twr.rBase) / wt->hTwr;

    // varying variables
    PetscReal dh, h = 0;

    for(PetscInt pi=0; pi<wt->twr.nPoints; pi++)
    {
        if(pi==0 || pi==(wt->twr.nPoints-1)) dh = dhval / 2;
        else dh = dhval;

        // bottom and top radiuses of this tower segment
        PetscReal r_lo = rastCoeff*h + wt->twr.rBase;
        PetscReal r_hi = rastCoeff*(h+dh) + wt->twr.rBase;

        // compute area using trapezoidal formula for rastremation
        wt->twr.dA[pi] = 0.5 * (r_lo + r_hi) * dh;

        // set current tower point
        Cmpnts twr_point = nScale(h, wt->twrDir);
                           mSum(twr_point, base);

        // save tower point
        wt->twr.points[pi] = nSet(twr_point);

        h += dh;
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initNacModel(windTurbine *wt, Cmpnts &base)
{
    // set tower thrust to zero
    wt->nac.nacThrust = 0.0;

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->nac.prjNSigma = sqrt(std::log(1.0/0.001));

    // set nacelle point
    wt->nac.point     = nSum(nScale(wt->hTwr, wt->twrDir), base);

    // set nacelle area
    wt->nac.A         = M_PI*wt->rHub*wt->rHub;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode initControlledCells(farm_ *farm)
{
    // initialize the controlled cells labels.
    // Loop over wind turbines and check if this processor
    // controls that wind turbine by defining a sphere of cells around the
    // rotor where this processor have access. If the size of the sphere cell
    // is zero this processor does not control the wind turbine, so set
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

                    PetscReal   eps_x =     wt->alm.chord[0] * wt->eps_x,
                                eps_y =     wt->alm.thick[0 ]* wt->eps_y;

                    
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
            radius_t = PetscMax(wt->eps_x, PetscMax(wt->eps_y, wt->eps_z)) * wt->prjNSigma;
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

PetscErrorCode writeFarmADMesh(farm_ *farm)
{
    mesh_     *mesh = farm->access->mesh;
    clock_    *clock = farm->access->clock;

    PetscMPIInt rank; MPI_Comm_rank(mesh->MESH_COMM, &rank);

    if(!rank)
    {
        // compute number of points and cells
        PetscInt npts = 0;
        PetscInt ncll = 0;
        for(PetscInt t=0; t<farm->size; t++)
        {
            PetscInt npts_t = 0, nrc_t = 0, nac_t = 0;

            // actuator disk model
            if((*farm->turbineModels[t]) == "ADM")
            {
                // number of points in this turbine's AD mesh
                npts_t = farm->wt[t]->adm.nPoints;

                // number of radial cells in this turbine's AD mesh
                nrc_t  = farm->wt[t]->adm.nRadial - 1;

                // number of azimuthal cells in this turbine's AD mesh
                nac_t  = farm->wt[t]->adm.nAzimuth;
            }
            // uniform actuator disk model
            if((*farm->turbineModels[t]) == "uniformADM")
            {
                // number of points in this turbine's AD mesh
                npts_t = farm->wt[t]->uadm.nPoints;

                // number of radial cells in this turbine's AD mesh
                nrc_t  = farm->wt[t]->uadm.nRadial - 1;

                // number of azimuthal cells in this turbine's AD mesh
                nac_t  = farm->wt[t]->uadm.nAzimuth;
            }

            npts += npts_t;
            ncll += nrc_t*nac_t;
        }

        if(npts>0)
        {
            word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

            char fileName[256];
            sprintf(fileName, "%s/ADMesh_%ld.inp", turbineFolderTimeName.c_str(), farm->writeNumber);

            PetscInt width = -20;

            FILE *f = fopen(fileName, "w");

            // header
            PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");

            // number of points, number of cells
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

            // write coordinates
            npts = 1;

            for(PetscInt t=0; t<farm->size; t++)
            {
                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    // number of points in this turbine's AD mesh
                    PetscInt npts_t = farm->wt[t]->adm.nPoints;

                    for(PetscInt pi=0; pi<npts_t; pi++)
                    {
                        Cmpnts point_i;

                        point_i.x = farm->wt[t]->adm.points[pi].x;
                        point_i.y = farm->wt[t]->adm.points[pi].y;
                        point_i.z = farm->wt[t]->adm.points[pi].z;

                        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, npts, width, point_i.x, width, point_i.y, width, point_i.z);

                        npts++;
                    }
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    // number of points in this turbine's AD mesh
                    PetscInt npts_t = farm->wt[t]->uadm.nPoints;

                    for(PetscInt pi=0; pi<npts_t; pi++)
                    {
                        Cmpnts point_i;

                        point_i.x = farm->wt[t]->uadm.points[pi].x;
                        point_i.y = farm->wt[t]->uadm.points[pi].y;
                        point_i.z = farm->wt[t]->uadm.points[pi].z;

                        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, npts, width, point_i.x, width, point_i.y, width, point_i.z);

                        npts++;
                    }
                }
            }

            // write connectivity

            // storage for the labels
            PetscInt cellPtLabels[ncll][4];

            npts = 0;
            ncll = 0;

            // build connectivity
            for(PetscInt t=0; t<farm->size; t++)
            {
                PetscInt nRadCells_t, nAziCells_t;
                PetscInt nRadPts_t, nAziPts_t;

                // actuator disk model
                if((*farm->turbineModels[t]) == "ADM")
                {
                    // radial and azimuthal cells
                    nRadCells_t = farm->wt[t]->adm.nRadial - 1;
                    nAziCells_t = farm->wt[t]->adm.nAzimuth;

                    // radial and azimuthal points
                    nRadPts_t = farm->wt[t]->adm.nRadial;
                    nAziPts_t = farm->wt[t]->adm.nAzimuth;
                }
                else if((*farm->turbineModels[t]) == "uniformADM")
                {
                    // radial and azimuthal cells
                    nRadCells_t = farm->wt[t]->uadm.nRadial - 1;
                    nAziCells_t = farm->wt[t]->uadm.nAzimuth;

                    // radial and azimuthal points
                    nRadPts_t = farm->wt[t]->uadm.nRadial;
                    nAziPts_t = farm->wt[t]->uadm.nAzimuth;
                }

                // number of points per turbine
                PetscInt npt = nRadPts_t*nAziPts_t;

                for(PetscInt ri=0; ri<nRadPts_t; ri++)
                {
                    for(PetscInt ai=0; ai<nAziPts_t; ai++)
                    {
                        if(ri < nRadCells_t)
                        {
                            // all azi cells except last
                            if(ai < (nAziCells_t-1))
                            {
                                cellPtLabels[ncll][0] = (t*npt) + (ri * nAziPts_t + ai + 1);
                                cellPtLabels[ncll][1] = (t*npt) + (ri * nAziPts_t + ai + 2);
                                cellPtLabels[ncll][2] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 2);
                                cellPtLabels[ncll][3] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 1);

                                ncll++;
                            }
                            // last cell: pts 0,3 must be connected to pts 0,3 of first azi cell
                            //            so pt 1,2 become 0,3 of first azi cell
                            else if (ai == (nAziCells_t-1))
                            {
                                cellPtLabels[ncll][0] = (t*npt) + (ri * nAziPts_t + ai + 1);
                                cellPtLabels[ncll][1] = (t*npt) + (ri * nAziPts_t + 1);
                                cellPtLabels[ncll][2] = (t*npt) + ((ri + 1) * nAziPts_t + 1);
                                cellPtLabels[ncll][3] = (t*npt) + ((ri + 1) * nAziPts_t + ai + 1);

                                ncll++;
                            }
                        }
                        npts++;
                    }
                }
            }

            width = -5;

            // write
            for(PetscInt ci=0; ci<ncll; ci++)
            {
                PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d %*d %*d\n", width, ci+1, width, 0, width, "quad", width, cellPtLabels[ci][0], width, cellPtLabels[ci][1], width, cellPtLabels[ci][2], width, cellPtLabels[ci][3]);
            }

            fclose(f);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeFarmTwrMesh(farm_ *farm)
{
    mesh_  *mesh  = farm->access->mesh;
    clock_ *clock = farm->access->clock;

    PetscMPIInt           rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    if(!rank)
    {
        // compute number of points and edges
        // each tower is only one line
        PetscInt npts = 2*farm->size;
        PetscInt ncll = farm->size;

        word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

        char fileName[256];
        sprintf(fileName, "%s/twrMesh.inp", turbineFolderTimeName.c_str());

        PetscInt width = -20;

        FILE *f = fopen(fileName, "w");

        // header
        PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
        PetscFPrintf(mesh->MESH_COMM, f, "#\n");
        PetscFPrintf(mesh->MESH_COMM, f, "#\n");

        // number of points, number of cells
        PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

        for(PetscInt t=0; t<farm->size; t++)
        {
            // set tower top point
            Cmpnts towerPt  = nScale(farm->wt[t]->hTwr, farm->wt[t]->twrDir);
                              mSum(towerPt, farm->base[t]);

            // set tower base point
            Cmpnts basePt   = nSet(farm->base[t]);

            // set point labels
            PetscInt id_base = 2*(t + 1) - 1;
            PetscInt id_twr  = 2*(t + 1);

            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, id_base, width, basePt.x, width, basePt.y, width, basePt.z);
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, id_twr, width, towerPt.x, width, towerPt.y, width, towerPt.z);
        }

        // write connectivity
        width = -5;

        // build connectivity
        for(PetscInt t=0; t<farm->size; t++)
        {
            // set point labels
            PetscInt id_base = 2*(t + 1) - 1;
            PetscInt id_twr  = 2*(t + 1);

            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d\n", width, t+1, width, 0, width, "line", width, id_base, width, id_twr);
        }

        fclose(f);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode writeFarmALMesh(farm_ *farm)
{
    mesh_      *mesh = farm->access->mesh;
    clock_    *clock = farm->access->clock;

    PetscMPIInt rank;
    MPI_Comm_rank(mesh->MESH_COMM, &rank);

    PetscInt    point_id, cell_id;

    if(!rank)
    {
        // compute number of points and cells
        PetscInt npts = 0;
        PetscInt ncll = 0;
        for(PetscInt t=0; t<farm->size; t++)
        {
            if((*farm->turbineModels[t]) == "ALM")
            {
                npts += 6;
                ncll += 3;
            }
        }

        if(npts>0)
        {
            word   turbineFolderTimeName = "./postProcessing/" + mesh->meshName + "/turbines/" + getStartTimeName(clock);

            char fileName[256];
            sprintf(fileName, "%s/ALMesh_%ld.inp", turbineFolderTimeName.c_str(), farm->writeNumber);

            PetscInt width = -20;

            FILE *f = fopen(fileName, "w");

            // header
            PetscFPrintf(mesh->MESH_COMM, f, "#UCD geometry file from TOSCA: Toolbox fOr Stratified Convective Atmospheres\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");
            PetscFPrintf(mesh->MESH_COMM, f, "#\n");

            // number of points, number of cells
            PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*d %*d %*d\n", width, npts, width, ncll, width, 0, width, 0, width, 0);

            // write coordinates
            point_id = 0;

            for(PetscInt t=0; t<farm->size; t++)
            {
                // define tip/root blade points in reference (vertical) configuration
                Cmpnts bldTip  = nScale(farm->wt[t]->rTip, farm->wt[t]->twrDir);
                Cmpnts bldRoot = nScale(farm->wt[t]->rHub, farm->wt[t]->twrDir);

                Cmpnts *points;
                PetscMalloc(6*sizeof(Cmpnts), &points);

                PetscReal angle = farm->wt[t]->alm.azimuth;

                for(PetscInt b=0; b<3; b++)
                {
                    Cmpnts bldTip_b  = nRot(farm->wt[t]->omega_hat, bldTip,  angle*farm->wt[t]->deg2rad);
                    Cmpnts bldRoot_b = nRot(farm->wt[t]->omega_hat, bldRoot, angle*farm->wt[t]->deg2rad);

                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, point_id, width, bldTip_b.x, width, bldTip_b.y, width, bldTip_b.z);
                    point_id++;
                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*.7f %*.7f %*.7f\n", width, point_id, width, bldRoot_b.x, width, bldRoot_b.y, width, bldRoot_b.z);
                    point_id++;

                    angle = angle + 120;
                }
            }

            // write connectivity
            width = -5;

            point_id = 0;
            cell_id  = 1;

            // build connectivity
            for(PetscInt t=0; t<farm->size; t++)
            {
                for(PetscInt b=0; b<3; b++)
                {
                    PetscFPrintf(mesh->MESH_COMM, f, "%*d %*d %*s %*d %*d\n", width, cell_id, width, 0, width, "line", width, point_id, width, point_id+1);
                    cell_id++;
                    point_id = point_id + 2;
                }
            }

            fclose(f);
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readFarmProperties(farm_ *farm)
{
    // The dictionary to read is windFarmProperties,
    // located inside ./turbines/ folder. It has a
    // first subdictionary for the write settings:
    //
    // writeSettings
    // {
    //     timeStart           10.0     // time at which the output starts to be written
    //     intervalType        timeStep // adjustableTime
    //     timeInterval        1        // output every second if adjustableTime or iteration if timeStep
    // }
    //
    // Then it contans the wind farm name and the wind farm
    // specification (onebyone or grid), if onebyone each turbine is
    // specified as follows
    //
    // WT1
    // {
    //     turbineType         NREL5MW
    //     turbineModel        ADM
    //     base                (5.0191 0.0 -90.0)
    //     windFarmController  1
    // }
    //
    // If grid we still have to implement it.

    // get mesh pointer
    mesh_ *mesh = farm->access->mesh;

    word   windFarmPropertiesFile = "./turbines/" + mesh->meshName + "/windFarmProperties";

    // allocate memory for wind farm body force
    VecDuplicate(mesh->lCent, &(farm->lsourceFarmCat));
    VecDuplicate(mesh->Cent,  &(farm->sourceFarmCont));

    // read the write settings
    readSubDictDouble(windFarmPropertiesFile.c_str(),"writeSettings","timeStart", &(farm->timeStart));
    readSubDictWord  (windFarmPropertiesFile.c_str(),"writeSettings","intervalType", &(farm->intervalType));
    readSubDictDouble(windFarmPropertiesFile.c_str(),"writeSettings","timeInterval", &(farm->timeInterval));

    // initialize write number
    farm->writeNumber = 0;

    // read debug flag
    readDictInt(windFarmPropertiesFile.c_str(), "debug", &(farm->dbg));

    // read wind farm name
    readDictWord(windFarmPropertiesFile.c_str(), "windFarmName", &(farm->name));

    // read array specification
    word arraySpec;
    readDictWord(windFarmPropertiesFile.c_str(), "arraySpecification", &arraySpec);

    // read the wind turbines one by one until end of file
    if(arraySpec=="onebyone")
    {
        readTurbineArray(farm);
    }
    else if(arraySpec=="grid")
    {
        char error[512];
        sprintf(error, "array specification type %s not yet implemented\n", arraySpec.c_str());
        fatalErrorInFunction("readFarmProperties",  error);
    }
    else
    {
       char error[512];
        sprintf(error, "unknown array specification type %s\n", arraySpec.c_str());
        fatalErrorInFunction("readFarmProperties",  error);
    }

    // Allocate memory for the wind turbines, this is an array of pointers.
    // Each pointer will point to a turbine. Advantages: the size is the
    // one of a pointer for the allocation and the objects can be not
    // linear in memory.

    farm->wt = new windTurbine* [farm->size];

    // go through each wind turbine and read its description file
    for(PetscInt t=0; t<farm->size; t++)
    {
        // dynamic memory allocation (make the t-th pointer point to the struct)
        farm->wt[t] = new windTurbine;

        // set this turbine type and ID
        farm->wt[t]->type = *(farm->turbineTypes[t]);
        farm->wt[t]->id   = *(farm->turbineIds[t]);

        // allocate memory for the rotor axis (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->rtrAxis));

        // allocate memory for the rotor angular velocity unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->omega_hat));

        // allocate memory for the tower direction unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->twrDir));

        // allocate memory for the rotor direction unit vector (set when building the AD/AL turbine mesh)
        PetscMalloc(sizeof(Cmpnts), &(farm->wt[t]->rtrDir));

        // this turbine description file (can be shared if type is the same)
        word descrFile = "./turbines/" + mesh->meshName + "/" + (*farm->turbineTypes[t]);

        // read turbine properties
        readTurbineProperties(farm->wt[t], descrFile.c_str(), mesh->meshName, (*farm->turbineModels[t]));

        // read wind farm control table for this wind turbine
        if(farm->farmControlActive[t])
        {
            readWindFarmControlTable(farm->wt[t]);
        }
    }

    // set CFL checking flag to zero
    farm->checkCFL = 0;

    printFarmProperties(farm);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTurbineArray(farm_ *farm)
{
    // get mesh pointer
    mesh_ *mesh = farm->access->mesh;

    // define the local variables
    std::vector<std::string>   turbineTypes;
    std::vector<std::string>   turbineIds;
    std::vector<std::string>   turbineModels;
    std::vector<Cmpnts>        base;
    std::vector<PetscInt>      windFarmController;

    std::string turbineTypes_i,
                turbineIds_i,
                turbineModels_i;
    Cmpnts      base_i;
    PetscInt    controller_i;

    PetscInt    nturbines = 0;

    // pointer for strtod and strtol
    char        *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];
    std::string token;

    // dictionary name
    std::string dictName = "./turbines/" + mesh->meshName + "/windFarmProperties" ;

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName.c_str());
        fatalErrorInFunction("readTurbineArray",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "turbineArray",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read all turbines until end of file (should find "}" before eof)
                    while(!indata.eof())
                    {
                        // from here start to read the wind turbines

                        // read WT ID
                        indata >> word;
                        turbineIds_i = word;

                        // check if have hit the end of the list
                        if(trim(turbineIds_i)=="}")
                        {
                            if(nturbines>0)
                            {
                                // close the file
                                indata.close();

                                // allocate memory and store the variables
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineTypes));
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineIds));
                                PetscMalloc(nturbines*sizeof(std::string*), &(farm->turbineModels));
                                PetscMalloc(nturbines*sizeof(Cmpnts),  &(farm->base));
                                PetscMalloc(nturbines*sizeof(PetscInt),&(farm->farmControlActive));

                                // base locations
                                for(PetscInt p=0; p<nturbines; p++)
                                {
                                    PetscMalloc(sizeof(std::string), &(farm->turbineTypes[p]));
                                    PetscMalloc(sizeof(std::string), &(farm->turbineIds[p]));
                                    PetscMalloc(sizeof(std::string), &(farm->turbineModels[p]));

                                    // assign the pointers to the singly created variables in memory
                                    farm->turbineTypes[p]  = new std::string(turbineTypes[p]);
                                    farm->turbineIds[p]    = new std::string(turbineIds[p]);
                                    farm->turbineModels[p] = new std::string(turbineModels[p]);

                                    //PetscPrintf(mesh->MESH_COMM,"%s\n",(*farm->turbineTypes[p]).c_str());

                                    farm->base[p].x = base[p].x;
                                    farm->base[p].y = base[p].y;
                                    farm->base[p].z = base[p].z;

                                    farm->farmControlActive[p] = windFarmController[p];
                                }

                                // wind farm size
                                farm->size = nturbines;

                                // clear the local variables
                                std::vector<std::string>   ().swap(turbineTypes);
                                std::vector<std::string>   ().swap(turbineIds);
                                std::vector<std::string>   ().swap(turbineModels);
                                std::vector<Cmpnts>        ().swap(base);

                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected at least one turbine in subdictionary turbineArray of %s dictionary\n", dictName.c_str());
                                fatalErrorInFunction("readTurbineArray",  error);
                            }
                        }

                        // check if have hit another list (this was not closed with "}")
                        if(trim(turbineIds_i)=="}")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token in subdictionary turbineArray of %s dictionary\n", dictName.c_str());
                            fatalErrorInFunction("readTurbineArray",  error);
                        }

                        // store the turbine ID
                        turbineIds.push_back(turbineIds_i);

                        // parameters are enclosed by '()', so read the first
                        indata >> word;
                        token = word;
                        if(trim(token)=="(")
                        {
                            // read turbineType keyword
                            indata >> word;
                            token = word;
                            if(trim(token)=="turbineType")
                            {
                                // read the turbine type
                                indata >> word;
                                turbineTypes_i = word;
                                turbineTypes.push_back(turbineTypes_i);

                                // read turbineModel keyword
                                indata >> word;
                                token = word;
                                if(trim(token)=="turbineModel")
                                {
                                    // read the turbine model
                                    indata >> word;
                                    turbineModels_i = word;
                                    turbineModels.push_back(turbineModels_i);

                                    // read baseLocation keyword
                                    indata >> word;
                                    token = word;
                                    if(trim(token)=="baseLocation")
                                    {
                                        // read the turbine base vector

                                        // start reading vector ------------------------------------------------------------------------------------------------------------------------------------------------

                                        // the vector is in (x y z) format, so
                                        // 1. read the first parenthesis
                                        // 2. read the 3 doubles
                                        // 3. look for the closing parethesis

                                        // read the first component (contains "(" character)
                                        indata >> word;

                                        std::string first(word);
                                        if (first.find ("(") != std::string::npos)
                                        {
                                           // remove "("" character from the first component
                                           PetscInt l1 = first.size();
                                           for(PetscInt i=0;i<l1;i++)
                                           {
                                               // save the first component
                                               word[i] = word[i+1];
                                           }

                                           base_i.x = std::strtod(word, &eptr);

                                           // check if the first component is a PetscReal, throw error otherwise
                                           std::string cmp1(word);

                                           if(isNumber(cmp1))
                                           {
                                               if (cmp1.find ('.') == std::string::npos)
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                              char error[512];
                                               sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }

                                           // read the second component
                                           indata >> word;
                                           base_i.y = std::strtod(word, &eptr);

                                           // check if the second component is a PetscReal, throw error otherwise
                                           std::string cmp2(word);

                                           if(isNumber(cmp2))
                                           {
                                               if (cmp2.find ('.') == std::string::npos)
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                              char error[512];
                                               sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }

                                           // read the third component (contains ")" character)
                                           indata >> word;

                                           std::string last(word);
                                           if (last.find (")") != std::string::npos)
                                           {
                                               // remove ") character from the last component and store
                                               base_i.z = std::strtod(word, &eptr);

                                               // check if the first component is a PetscReal, throw error otherwise
                                               std::string cmp3(word);

                                               if(isNumber(cmp3))
                                               {
                                                   if (cmp3.find ('.') == std::string::npos)
                                                   {
                                                      char error[512];
                                                       sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                       fatalErrorInFunction("readTurbineArray",  error);
                                                   }
                                               }
                                               else
                                               {
                                                  char error[512];
                                                   sprintf(error, "expected <PetscReal> in vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }

                                               // save the base location
                                               base.push_back(base_i);

                                               // read windFarmController keyword
                                               indata >> word;
                                               token = word;
                                               if(trim(token)=="windFarmController")
                                               {
                                                   // read wether controller is active or not
                                                   indata >> word;
                                                   controller_i = (PetscInt)std::strtol(word, &eptr, 10);
                                                   windFarmController.push_back(controller_i);

                                                   // increase turbine counter
                                                   nturbines++;

                                                   // read the closing ')'
                                                   indata >> word;
                                                   token = word;
                                                   if(trim(token)!=")")
                                                   {
                                                       char error[512];
                                                       sprintf(error, "expected <)>  at end of subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                                       fatalErrorInFunction("readTurbineArray",  error);
                                                   }
                                               }
                                               else
                                               {
                                                   char error[512];
                                                   sprintf(error, "expected windFarmController keyword after %s in subdictionary %s of %s dictionary, found '%s'\n", last.c_str(), turbineIds_i.c_str(), dictName.c_str(), word);
                                                   fatalErrorInFunction("readTurbineArray",  error);
                                               }
                                           }
                                           else
                                           {
                                               char error[512];
                                               sprintf(error, "expected <(>  after vector defined by keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                               fatalErrorInFunction("readTurbineArray",  error);
                                           }
                                       }
                                       else
                                       {
                                           char error[512];
                                           sprintf(error, "expected <(>  after keyword baseLocation in subdictionary %s of %s dictionary, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                                           fatalErrorInFunction("readTurbineArray",  error);
                                       }

                                       // End of reading vector ------------------------------------------------------------------------------------------------------------------------------------------------
                                    }
                                    else
                                    {
                                        char error[512];
                                        sprintf(error, "expected turbineModel keyword after %s in subdictionary %s of %s dictionary, found '%s'\n", turbineModels_i.c_str(), turbineIds_i.c_str(), dictName.c_str(), word);
                                        fatalErrorInFunction("readTurbineArray",  error);
                                    }
                                }
                                else
                                {
                                   char error[512];
                                    sprintf(error, "expected turbineModel keyword after %s in sub-dictionary %s, found '%s'\n", turbineTypes_i.c_str(), turbineIds_i.c_str(), word);
                                    fatalErrorInFunction("readTurbineArray",  error);
                                }
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected turbineType keyword after ( in sub-dictionary %s, found '%s'\n", turbineIds_i.c_str(), word);
                                fatalErrorInFunction("readTurbineArray",  error);
                            }
                        }
                        else
                        {
                           char error[512];
                            sprintf(error, "expected <(> token after keyword %s in dictionary %s, found '%s'\n", turbineIds_i.c_str(), dictName.c_str(), word);
                            fatalErrorInFunction("readTurbineArray",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of turbineArray subdictionary in %s dictionary\n", dictName.c_str());
                    fatalErrorInFunction("readTurbineArray",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword turbineArray in dictionary %s, found '%s'\n", dictName.c_str(), word);
                    fatalErrorInFunction("readTurbineArray",  error);
                }

            }
        }

       char error[512];
        sprintf(error, "could not find subdictionary turbineArray in dictionary %s\n", dictName.c_str());
        fatalErrorInFunction("readTurbineArray",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTurbineProperties(windTurbine *wt, const char *dictName, const word meshName, const word modelName)
{
    // local vectors (to be normalised)
    Cmpnts twrDirVec, rtrDirVec;

    // read main parameters
    readDictDouble(dictName, "rTip",        &(wt->rTip));
    readDictDouble(dictName, "rHub",        &(wt->rHub));
    readDictDouble(dictName, "hTower",      &(wt->hTwr));
    readDictDouble(dictName, "overHang",    &(wt->ovrHang));
    readDictDouble(dictName, "precone",     &(wt->precone));
    readDictVector(dictName, "towerDir",    &(twrDirVec));
    readDictVector(dictName, "rotorDir",    &(rtrDirVec));
    readDictDouble(dictName, "upTilt",      &(wt->upTilt));

    if(modelName == "ADM" || modelName == "uniformADM")
    {
        readDictDouble(dictName, "epsilon",     &(wt->eps));
    }
    else if(modelName == "AFM")
    {
        readDictDouble(dictName, "epsilon_x",     &(wt->eps_x));
        readDictDouble(dictName, "epsilon_y",     &(wt->eps_y));
        readDictDouble(dictName, "epsilon_z",     &(wt->eps_z));
    }
    else if(modelName == "ALM")
    {
        //do nothing - read later based on projection type
    }

    // set projection confidence interval as number of standard deviations (harcoded to 2.7)
    wt->prjNSigma = sqrt(std::log(1.0/0.001));

    // tower model flag
    readDictInt(dictName,    "includeTower",&(wt->includeTwr));

    // nacelle model flag
    readDictInt(dictName,    "includeNacelle",&(wt->includeNacelle));

    // debug switch
    readDictInt(dictName,    "debug",     &(wt->dbg));

    // store normalized vectors
    wt->twrDir = nUnit(twrDirVec);
    wt->rtrDir = nUnit(rtrDirVec);

    // conversion factors
    wt->deg2rad = M_PI / 180.0;
    wt->rad2deg = 180.0 / M_PI;
    wt->rpm2RadSec = 2 * M_PI / 60.0;

    // initialize rotor omega to zero (for CFL check access)
    wt->rtrOmega = 0.0;

    // read parameters for AD/AL models
    if(modelName == "ADM" || modelName == "ALM")
    {
        PetscReal initOmega;

        readDictInt(dictName,    "nBlades",     &(wt->nBlades));
        readDictWord(dictName,   "rotationDir", &(wt->rotDir));
        readDictDouble(dictName, "initialOmega",&(initOmega));
        wt->rtrOmega = initOmega * wt->rpm2RadSec;

        // read torque controller
        readDictWord(dictName,   "genControllerType", &(wt->genControllerType));

        // read pitch controller type
        readDictWord(dictName,   "pitchControllerType", &(wt->pitchControllerType));

        // read controllers input parameters
        if(wt->genControllerType   != "none") readGenControllerParameters(wt,   wt->genControllerType.c_str(), meshName.c_str());
        if(wt->pitchControllerType != "none") readPitchControllerParameters(wt, wt->pitchControllerType.c_str(), meshName.c_str());

        // read airfoil types used in this turbine
        readAirfoilProperties(wt, dictName);

        // Allocate memory for the foil info, this is an array of pointers.
        // Each pointer will point to an airfoil. Advantages: the size is the
        // one of a pointer for the allocation and the objects can be not
        // linear in memory.
        wt->foils = (foilInfo**)malloc(wt->nFoils*sizeof(foilInfo*));

        // read the property file for each airfoil
        for(PetscInt f=0; f<wt->nFoils; f++)
        {
            // dynamic memory allocation (make the f-th pointer point to the struct)
            wt->foils[f] = (foilInfo*)malloc(sizeof(foilInfo));

            // set the name of the airfoil
            PetscMalloc(sizeof(word), &(wt->foils[f]->name));
            wt->foils[f]->name = *(wt->foilNames[f]);

            // set the path to the airfoil data file
            word name2af = "./turbines/" + meshName + "/airfoils/" + *(wt->foilNames[f]);

            // set the reset of the variables by reading the table
            readAirfoilTable(wt->foils[f], name2af.c_str());
        }

        PetscInt readThickness = 0;

        if(modelName == "ALM")
        {
            word projectionType;
            readDictWord(dictName, "projection", &projectionType);

            if(projectionType=="anisotropic")
            {
                readThickness = 1;
            }
        }

        // read blade properties
        readBladeProperties(wt, dictName, readThickness);

        // check that the max airfoil label in the blade properties actually
        // matches the number if provided airfoils
        PetscInt max_id = 0;

        // find max id
        for(PetscInt i=0; i<wt->blade.size; i++)
        {
            if(wt->blade.foilIds[i] > max_id)
            {
                max_id = wt->blade.foilIds[i];
            }
        }

        if(max_id+1>wt->nFoils)
        {
           char error[512];
            sprintf(error, "requested more airfoils than the number provided in turbine type %s (airfoils and bladeData mismatch)\n",wt->type.c_str());
            fatalErrorInFunction("readTurbineProperties",  error);
        }
    }

    // read yaw controller parameters (all models: ADM/ALM/UADM/AFM)
    readDictWord(dictName,   "yawControllerType", &(wt->yawControllerType));

    if(wt->yawControllerType   != "none") readYawControllerParameters(wt,   wt->yawControllerType.c_str(), meshName.c_str());

    if(wt->includeTwr)
    {
        readTowerProperties(wt, dictName);
    }

    if(wt->includeNacelle)
    {
        readNacelleProperties(wt, dictName);
    }

    // set yaw changed parameter to 1 at initialization
    // this makes sure all models do the first cell to point search
    // ALM always does the search as blades are rotating
    wt->yawChanged = 1;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readAirfoilProperties(windTurbine *wt, const char *dictName)
{
    // First a list of names is read from the dictionary passed as argument.
    // This must be defined as follows:
    //
    // airfoils
    // {
    //     name1
    //     ...
    //     nameN
    // }
    //
    // Then the each airfoil look up tables are read using the readAirfoilTable
    // function.

    // define the local variables
    std::vector<std::string>   foilNames;

    std::string foilNames_i;
    PetscInt         nfoils = 0;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readAirfoilProperties",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "airfoils",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the airfoil name
                        indata >> word;

                        // check that did not hit end of subdictionary "}"
                        std::string token(word);
                        if(trim(token)=="}")
                        {
                            if(nfoils>0)
                            {
                                // close file
                                indata.close();

                                // allocate memory
                                PetscMalloc(nfoils*sizeof(std::string*), &(wt->foilNames));

                                // store the data
                                for(PetscInt f=0; f<nfoils; f++)
                                {
                                    PetscMalloc(sizeof(std::string), &(wt->foilNames[f]));

                                    wt->foilNames[f] = new std::string(foilNames[f]);
                                }

                                wt->nFoils = nfoils;

                                // clear memory
                                std::vector<std::string> ().swap(foilNames);

                                // exit
                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 1 entry as reading airfoils subdictionary in %s dictionary\n", dictName);
                                fatalErrorInFunction("readAirfoilProperties",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of airfoils subdictionary in %s dictionary\n", dictName);
                            fatalErrorInFunction("readAirfoilProperties",  error);
                        }

                        // store the airfoil name
                        foilNames_i = word;
                        foilNames.push_back(foilNames_i);

                        nfoils++;

                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of airfoils subdictionary in %s dictionary\n", dictName);
                    fatalErrorInFunction("readAirfoilProperties",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword airfoils in dictionary %s, found '%s'\n", dictName, word);
                    fatalErrorInFunction("readAirfoilProperties",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword airfoils in dictionary %s\n", dictName);
        fatalErrorInFunction("readAirfoilProperties",  error);

    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readAirfoilTable(foilInfo *af, const char *tableName)
{
    // Airfoil tables must be provided as follows:
    //
    // table
    // {
    //     (alpha_deg1  cl1  cd1)
    //               ...
    //     (alpha_degN  clN  cdN)
    // }
    //
    // Suggestion: in order for the AD an AL models to work
    // properly, provide the curves between -180 to +180 degs.
    // The number of points is arbitrary, minimum is 3.

    // define the local variables
    std::vector<PetscReal>   aoa;
    std::vector<PetscReal>   cl;
    std::vector<PetscReal>   cd;

    PetscReal aoa_i, cl_i, cd_i;
    PetscInt    nPoints = 0;

    // pointer for the strtod function
    char   *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(tableName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s file\n", tableName);
        fatalErrorInFunction("readAirfoilTable",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "table",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the first component, contains "(" character
                        indata >> word;

                        // put the word in a string for performing checks
                        std::string token(word);

                        std::string first(word);
                        if (first.find ("(") != std::string::npos)
                        {
                            // remove "( character from the first component
                            PetscInt l1 = first.size();
                            for(PetscInt i=0;i<l1;i++)
                            {
                                word[i] = word[i+1];
                            }

                            // store the first variable (aoa)
                            aoa_i = std::strtod(word, &eptr);
                            aoa.push_back(aoa_i);

                            // read the second component (cl)
                            indata >> word;

                            // store second component
                            cl_i = std::strtod(word, &eptr);
                            cl.push_back(cl_i);

                            // read the third component (cd), contains ")" character
                            indata >> word;

                            std::string last(word);
                            if (last.find (")") != std::string::npos)
                            {
                                // remove ") character from the last component and store
                                cd_i = std::strtod(word, &eptr);
                                cd.push_back(cd_i);

                                // increament point counter
                                nPoints++;
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected <)>  at end of line as reading table in %s dictionary\n", tableName);
                                fatalErrorInFunction("readAirfoilTable",  error);
                            }
                        }
                        else if(trim(token)=="}")
                        {
                            if(nPoints >= 3)
                            {
                                // close file
                                indata.close();

                                // allocate memory
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->aoa));
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->cl));
                                PetscMalloc(nPoints*sizeof(PetscReal), &(af->cd));

                                // store the table
                                for(PetscInt i=0; i<nPoints; i++)
                                {
                                    af->aoa[i] = aoa[i];
                                    af->cl[i]  = cl[i];
                                    af->cd[i]  = cd[i];
                                }

                                // store number of points in the table
                                af->size = nPoints;

                                // clear memory
                                std::vector<PetscReal>   ().swap(aoa);
                                std::vector<PetscReal>   ().swap(cl);
                                std::vector<PetscReal>   ().swap(cd);

                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 3 data points as reading table in %s dictionary\n", tableName);
                                fatalErrorInFunction("readAirfoilTable",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of table in %s dictionary\n", tableName);
                            fatalErrorInFunction("readAirfoilTable",  error);
                        }
                        // we are at a new line and neither '(' nor '}' were found: throws error
                        else
                        {
                             char error[512];
                              sprintf(error, "expected either <(>  or <}> at new line as reading table in %s dictionary\n", tableName);
                              fatalErrorInFunction("readAirfoilTable",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of table in %s dictionary\n", tableName);
                    fatalErrorInFunction("readAirfoilTable",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword table in dictionary %s, found '%s'\n", tableName, word);
                    fatalErrorInFunction("readAirfoilTable",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword 'table' in file %s\n", tableName);
        fatalErrorInFunction("readAirfoilTable",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readBladeProperties(windTurbine *wt, const char *dictName, const PetscInt readThickness)
{
    // The dictionary must be located inside a
    // file named as the turbine type and located
    // inside ./turbines/, e.g.: ./turbines/NREL5MW.dat.
    // It must be defined as follows:
    //
    // bladeData
    // {
    //     (radius1 chord1 twist1 foilIds1)
    //     (radius2 chord2 twist2 foilIds2)
    //                   ...
    //     (radiusN chordN twistN foilIdsN)
    // }
    //
    // two is the minimum number of data points.
    // no checks for doubles or integers are performed,
    // only the format is checked.
    // For anisotropic gaussian ALM projection also blade thickness should be provided in the last column

    // define the local variables
    std::vector<PetscReal> radius;
    std::vector<PetscReal> chord;
    std::vector<PetscReal> twist;
    std::vector<PetscReal> thick;
    std::vector<PetscInt> foilIds;

    PetscReal   radius_i, chord_i, twist_i, thick_i;
    PetscInt    foilIds_i;
    PetscInt    nlines = 0;

    // pointer for strtod and strtol
    char   *eptr;

    // file stream
    std::ifstream indata;

    // word by word read
    char word[256];

    // open dictionary
    indata.open(dictName);

    if(!indata)
    {
       char error[512];
        sprintf(error, "could not open %s dictionary\n", dictName);
        fatalErrorInFunction("readBladeProperties",  error);
    }
    else
    {
        // get word by word till end of dictionary
        while(!indata.eof())
        {
            indata >> word;

            // test if found subdictionary
            if
            (
                strcmp
                (
                    "bladeData",
                    word
                ) == 0
            )
            {
                // read the first "{"
                indata >> word;

                std::string token1(word);
                std::string token2;

                // test if braket is the first word after the subdictionary entry
                if(trim(token1)=="{")
                {
                    // read until end of file (should find "}" before)
                    while(!indata.eof())
                    {
                        // read the first component, contains "(" character
                        indata >> word;

                        // put the word in a string for performing checks
                        std::string token3(word);

                        std::string first(word);
                        if (first.find ("(") != std::string::npos)
                        {
                            // remove "( character from the first component
                            PetscInt l1 = first.size();
                            for(PetscInt i=0;i<l1;i++)
                            {
                                word[i] = word[i+1];
                            }

                            // store the first variable (radius)
                            radius_i = std::strtod(word, &eptr);
                            radius.push_back(radius_i);

                            // read the second component (chord)
                            indata >> word;

                            // store second component
                            chord_i = std::strtod(word, &eptr);
                            chord.push_back(chord_i);

                            // read the third component (twist)
                            indata >> word;

                            // store third component
                            twist_i = std::strtod(word, &eptr);
                            twist.push_back(twist_i);

                            if(readThickness)
                            {
                                // read the thickness
                                indata >> word;

                                // store thickness component
                                thick_i = std::strtod(word, &eptr);
                                thick.push_back(thick_i);
                            }

                            // read the fourth component (foilIds), contains ")" character
                            indata >> word;

                            std::string last(word);
                            if (last.find (")") != std::string::npos)
                            {
                                // remove ") character from the last component and store
                                foilIds_i = std::strtol(word, &eptr, 10);
                                foilIds.push_back(foilIds_i);

                                // increament line counter
                                nlines++;
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "expected <)>  at end of line as reading readData table in %s dictionary\n", dictName);
                                fatalErrorInFunction("readBladeProperties",  error);
                            }

                        }
                        // look for the terminating "}" if found: close the file, store the data and exit
                        // if not found: may be another line
                        else if(trim(token3)=="}")
                        {
                            if(nlines >= 2)
                            {
                                // close file
                                indata.close();

                                // allocate memory for the blade properties in the current wind turbine
                                PetscMalloc(sizeof(bladeAeroInfo), &(wt->blade));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.radius));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.chord));
                                PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.twist));
                                PetscMalloc(nlines*sizeof(PetscInt), &(wt->blade.foilIds));

                                if(readThickness) PetscMalloc(nlines*sizeof(PetscReal), &(wt->blade.thick));

                                // store the blade properties
                                wt->blade.size = nlines;

                                for(PetscInt p=0; p<wt->blade.size; p++)
                                {
                                    wt->blade.radius[p]  = radius[p];
                                    wt->blade.chord[p]   = chord[p];
                                    wt->blade.twist[p]   = twist[p];
                                    wt->blade.foilIds[p] = foilIds[p];

                                    if(readThickness) wt->blade.thick[p]   = thick[p];
                                }

                                // clean the local variables
                                std::vector<PetscReal> ().swap(radius);
                                std::vector<PetscReal> ().swap(chord);
                                std::vector<PetscReal> ().swap(twist);
                                std::vector<PetscReal> ().swap(thick);
                                std::vector<PetscInt>    ().swap(foilIds);

                                // exit
                                return(0);
                            }
                            else
                            {
                               char error[512];
                                sprintf(error, "Required at least 2 data points as reading readData table in %s dictionary\n", dictName);
                                fatalErrorInFunction("readBladeProperties",  error);
                            }
                        }
                        // if find another "{" means another subdict is entered: throws error
                        else if(trim(token3)=="{")
                        {
                           char error[512];
                            sprintf(error, "missing '}' token at end of readData table in %s dictionary\n", dictName);
                            fatalErrorInFunction("readBladeProperties",  error);
                        }
                        // we are at a new line and neither '(' nor '}' were found: throws error
                        else
                        {
                             char error[512];
                              sprintf(error, "expected either <(>  or <}> at new line as reading readData table in %s dictionary\n", dictName);
                              fatalErrorInFunction("readBladeProperties",  error);
                        }
                    }

                    // have reached this point without finding }: throws error
                   char error[512];
                    sprintf(error, "missing '}' token at end of readData table in %s dictionary\n", dictName);
                    fatalErrorInFunction("readBladeProperties",  error);
                }
                else
                {
                   char error[512];
                    sprintf(error, "expected '{' token after keyword bladeData in dictionary %s, found '%s'\n", dictName, word);
                    fatalErrorInFunction("readBladeProperties",  error);
                }
            }
        }

       char error[512];
        sprintf(error, "could not find keyword bladeData in dictionary %s\n", dictName);
        fatalErrorInFunction("readBladeProperties",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readTowerProperties(windTurbine *wt, const char *dictName)
{
    // allocate memory for the tower model
    PetscMalloc(sizeof(towerModel), &(wt->twr));

    readSubDictDouble(dictName, "towerData", "Cd",      &(wt->twr.Cd));
    readSubDictDouble(dictName, "towerData", "epsilon", &(wt->twr.eps));
    readSubDictInt(dictName,    "towerData", "nLinPts", &(wt->twr.nPoints));
    readSubDictDouble(dictName, "towerData", "rBase",   &(wt->twr.rBase));
    readSubDictDouble(dictName, "towerData", "rTop",    &(wt->twr.rTop));

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readNacelleProperties(windTurbine *wt, const char *dictName)
{
    // allocate memory for the tower model
    PetscMalloc(sizeof(nacelleModel), &(wt->nac));

    readSubDictDouble(dictName, "nacelleData", "Cd",      &(wt->nac.Cd));
    readSubDictDouble(dictName, "nacelleData", "epsilon", &(wt->nac.eps));

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readGenControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // rotor dynamics parameters
    readDictDouble(path2dict, "genInertia",         &(wt->genInertia));
    readDictDouble(path2dict, "hubInertia",         &(wt->hubInertia));
    readDictDouble(path2dict, "bldIntertia",        &(wt->bldInertia));
    readDictDouble(path2dict, "gbxRatioG2R",        &(wt->gbxRatioG2R));
    readDictDouble(path2dict, "gbxEfficiency",      &(wt->gbxEff));

    // initial conditions
    readDictDouble(path2dict, "gbxRatioG2R",        &(wt->gbxRatioG2R));
    readDictDouble(path2dict, "initialGenTorque",   &(wt->genTorque));

    // generator electrical efficiency
    readDictDouble(path2dict,  "genEff",            &(wt->genEff));

    // controller parameters
    readSubDictDouble(path2dict, "genTqControllerParameters","genSpeedFilterFreq", &(wt->rtrSpdFilterFreq));
    readSubDictDouble(path2dict, "genTqControllerParameters","cutInGenSpeed",      &(wt->cutInGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","cutInGenTorque",     &(wt->cutInGenTq));
    readSubDictDouble(path2dict, "genTqControllerParameters","regTwoStartGenSpeed",&(wt->regTwoStartGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","regTwoEndGenSpeed",  &(wt->regTwoEndGenSpd));
    readSubDictDouble(path2dict, "genTqControllerParameters","ratedGenTorque",     &(wt->ratedGenTq));
    readSubDictDouble(path2dict, "genTqControllerParameters","controllerPGain",    &(wt->omegaKP));

    // limits
    readSubDictInt(path2dict,    "genTqControllerParameters","torqueRateLimiter",  &(wt->tqRateLimiter));
    readSubDictInt(path2dict,    "genTqControllerParameters","rtrSpeedLimiter",    &(wt->rtrSpdLimiter));
    readSubDictDouble(path2dict, "genTqControllerParameters","torqueMaxRate",      &(wt->tqMaxRate));
    readSubDictDouble(path2dict, "genTqControllerParameters","ratedRotorSpeed",    &(wt->ratedRotorSpd));

    // convert rpm to rad/s
    wt->ratedRotorSpd     *= wt->rpm2RadSec;
    wt->regTwoEndGenSpd   *= wt->rpm2RadSec;
    wt->regTwoStartGenSpd *= wt->rpm2RadSec;
    wt->cutInGenSpd       *= wt->rpm2RadSec;
    wt->omegaKP           /= (wt->rpm2RadSec*wt->rpm2RadSec);
    wt->rtrSpdFilterFreq  /= 60;

    // total inertia on HS LS shafts
    wt->driveTrainInertia
    =
    wt->nBlades*wt->bldInertia +
    wt->hubInertia +
    wt->gbxRatioG2R*wt->gbxRatioG2R*wt->genInertia;

    // generator speed in rpm
    wt->genOmega  = wt->rtrOmega * wt->gbxRatioG2R;

    // fitered rotor omega
    wt->rtrOmegaFilt = wt->rtrOmega;

    // generator power
    wt->genPwr    = 0.0;

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readPitchControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // allocate memory for collective pitch
    PetscMalloc(wt->nBlades*sizeof(PetscReal), &(wt->pitch));

    // PID controller paramters
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerPGain",   &(wt->pitchKP));
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerIGain",   &(wt->pitchKI));
    readSubDictDouble(path2dict, "pitchControllerParameters","controllerDGain",   &(wt->pitchKD));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchS2R",          &(wt->pitchS2R));

    // limits
    readSubDictInt(path2dict,    "pitchControllerParameters","pitchRateLimiter",  &(wt->pitchRateLimiter));
    readSubDictInt(path2dict,    "pitchControllerParameters","pitchAngleLimiter", &(wt->pitchAngleLimiter));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMaxRate",      &(wt->pitchMaxRate));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMin",          &(wt->pitchMin));
    readSubDictDouble(path2dict, "pitchControllerParameters","pitchMax",          &(wt->pitchMax));

    // transform all in radiants
    wt->pitchMaxRate *= wt->deg2rad;
    wt->pitchMax     *= wt->deg2rad;
    wt->pitchMin     *= wt->deg2rad;
    wt->pitchKP      *= wt->deg2rad;
    wt->pitchKI      *= wt->deg2rad;
    wt->pitchKD      *= wt->deg2rad;
    wt->pitchS2R     *= wt->deg2rad;

    // check and validate
    if(fabs(wt->pitchKI) < 1e-5)
    {
       char error[512];
        sprintf(error, "integral gain 'controlledIGain' cannot be zero\n");
        fatalErrorInFunction("readPitchControllerParameters",  error);
    }
    if(fabs(wt->pitchS2R) < 1e-5)
    {
       char error[512];
        sprintf(error, "point at which the power sensitivity to pitch is twice the one at rated conditions 'pitchS2R' cannot be zero\n");
        fatalErrorInFunction("readPitchControllerParameters",  error);
    }

    // set collective and individual blade pitch
    wt->collPitch = 0.0;
    for(PetscInt bld_i=0; bld_i<wt->nBlades; bld_i++)
    {
        wt->pitch[bld_i] = 0.0;
    }

    // set initial speed error to zero
    wt->errPID = 0.0;

    // set integrated error such that if the pitch is correct at the startTime
    // it will not change at the next time step
    PetscReal G   = 1.0 / (1.0 + wt->collPitch / wt->pitchS2R);
    wt->intErrPID = wt->collPitch / (G * wt->pitchKI);

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readYawControllerParameters(windTurbine *wt, const char *dictName, const char *meshName)
{
    char path2dict[256];
    sprintf(path2dict, "./turbines/%s/control/%s", meshName, dictName);

    // controller paramters
    readSubDictWord(path2dict,   "yawControllerParameters","sampleType",  &(wt->yawSamplingType));
    readSubDictDouble(path2dict, "yawControllerParameters","avgWindow",   &(wt->yawAverageWindow));
    readSubDictDouble(path2dict, "yawControllerParameters","yawMin",      &(wt->yawMin));
    readSubDictDouble(path2dict, "yawControllerParameters","yawMax",      &(wt->yawMax));
    readSubDictDouble(path2dict, "yawControllerParameters","yawSpeed",    &(wt->yawSpeed));
    readSubDictDouble(path2dict, "yawControllerParameters","allowedError",&(wt->yawAllowedError));

    // set initial flow angle
    readSubDictDouble(path2dict, "yawControllerParameters","initialFlowAngle", &(wt->flowAngle));

    // make sure yawSpeed is positive
    wt->yawSpeed = fabs(wt->yawSpeed);

    // check that sampling type is known
    if
    (
        wt->yawSamplingType != "hubUpDistance" &&
        wt->yawSamplingType != "anemometer"
    )
    {
       char error[512];
        sprintf(error, "unknown yaw sampleType %s. Known types are hubUpDistance and anemometer", wt->yawSamplingType.c_str());
        fatalErrorInFunction("readYawControllerParameters",  error);
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode readWindFarmControlTable(windTurbine *wt)
{
    word tableName = "turbines/control/" + wt->id;

    std::vector<std::vector<PetscReal>> table;

    // open file stream
    std::ifstream indata;
    indata.open(tableName);

    // word by word read
    char word[256];

    // buffer for read data
    PetscReal buffer;

    // time counter
    PetscInt ntimes;

    if(!indata)
    {
        char error[512];
        sprintf(error, "cannot open file %s\n", tableName.c_str());
        fatalErrorInFunction("readWindFarmControlTable",  error);
    }
    else
    {
        std::string tmpStr;

        // read lines and get number of saved times
        ntimes = 0;
        for (int t = 0; std::getline(indata, tmpStr); t++)
        {
            if (!tmpStr.empty())
            {
                ntimes++;
            }
        }

        // first line is header
        ntimes--;

        // save the number of times
        wt->wfControlNData  = ntimes;
        wt->currentCloseIdx = 0;

        // go back on top of file
        indata.close();
        indata.open(tableName);

        // skip header line
        std::getline(indata, tmpStr);

        // resize the source table
        table.resize(ntimes);

        for(PetscInt t=0; t<ntimes; t++)
        {
            // read along the line: time | value
            for(PetscInt i=0; i<2; i++)
            {
                table[t].resize(2);

                indata >> word;
                std::sscanf(word, "%lf", &buffer);

                table[t][i] = buffer;
            }

        }

        indata.close();
    }

    // now store the source  and free the temporary variable
    PetscMalloc(sizeof(PetscReal) * ntimes, &(wt->wfControlTimes));
    PetscMalloc(sizeof(PetscReal) * ntimes, &(wt->wfControlValues));

	// hard-coded difference between the sim. start and control action start
	PetscReal initialShift = 100000;

    for(PetscInt t=0; t<ntimes; t++)
    {
        wt->wfControlTimes[t]  = table[t][0] + initialShift;
        wt->wfControlValues[t] = table[t][1];

		//printf("time %.1f, ct: %.5f\n", wt->wfControlTimes[t], wt->wfControlValues[t]);
    }

    // clean the temporary variables
    for(PetscInt t=0; t<ntimes; t++)
    {
        std::vector<PetscReal> ().swap(table[t]);
    }

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

//***************************************************************************************************************//

PetscErrorCode printFarmProperties(farm_ *farm)
{
    mesh_ *mesh = farm->access->mesh;

    // check if at least one turbine has its debug flag activated
    PetscInt atLeastOneDebug = 0;

    for(PetscInt i=0;i<farm->size;i++)
    {
        if(farm->wt[i]->dbg)
        {
            atLeastOneDebug++;
        }
    }

    if(atLeastOneDebug)
    {
        PetscPrintf(mesh->MESH_COMM,"Wind Farm Properties for debugging\n\n");
    }

    for(PetscInt i=0;i<farm->size;i++)
    {
        if(farm->wt[i]->dbg)
        {
            // main turbine properties
            PetscPrintf(mesh->MESH_COMM,"Turbine %ld\n\n", i);
            PetscPrintf(mesh->MESH_COMM," ID      : %s\n", (*farm->turbineIds[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Type    : %s\n", (*farm->turbineTypes[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Model   : %s\n", (*farm->turbineModels[i]).c_str());
            PetscPrintf(mesh->MESH_COMM," Base    : (%lf, %lf, %lf)\n", farm->base[i].x, farm->base[i].y, farm->base[i].z);

            if((*farm->turbineModels[i]) != "uniformADM" && (*farm->turbineModels[i]) != "AFM")
            {
                PetscPrintf(mesh->MESH_COMM," Airfoils: ");

                for(PetscInt j=0;j<farm->wt[i]->nFoils; j++)
                {
                    PetscPrintf(mesh->MESH_COMM,"%s, ", (*farm->wt[i]->foilNames[j]).c_str());
                }

                PetscPrintf(mesh->MESH_COMM,"\n\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                PetscPrintf(mesh->MESH_COMM," Blade properties for this turbine:\n");

                // blade properties
                {
                    PetscInt n_bldpts = farm->wt[i]->blade.size;

                    PetscInt nch = -15;

                    PetscPrintf(mesh->MESH_COMM,"    | %*s | %*s | %*s| %*s|\n", nch, "radius [m]", nch, "chord [m]", nch, "twist [deg]", nch, "foil ID [-]");
                    for(PetscInt pti=0; pti<n_bldpts; pti++)
                    {
                        PetscReal radius   = farm->wt[i]->blade.radius[pti];
                        PetscReal chord    = farm->wt[i]->blade.chord[pti];
                        PetscReal twist    = farm->wt[i]->blade.twist[pti];
                        PetscInt    foilIds  = farm->wt[i]->blade.foilIds[pti];
                        word   foilName = (*farm->wt[i]->foilNames[foilIds]);

                        PetscPrintf(mesh->MESH_COMM,"    | %*.4f | %*.4f | %*.4f| %*s|\n", nch, radius, nch, chord, nch, twist, nch, foilName.c_str());
                    }
                }

                PetscPrintf(mesh->MESH_COMM,"\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                PetscPrintf(mesh->MESH_COMM," Breakdown of this turbine airfoil tables:\n");

                // airfoil properties
                for(PetscInt j=0;j<farm->wt[i]->nFoils; j++)
                {
                    // this airfoil table size
                    PetscInt n_afpts = farm->wt[i]->foils[j]->size;

                    PetscInt nch = -15;

                    PetscPrintf(mesh->MESH_COMM,"  \n%s:\n", (*farm->wt[i]->foilNames[j]).c_str());
                    PetscPrintf(mesh->MESH_COMM,"    | %*s | %*s | %*s|\n", nch, "alpha [deg]", nch, "Cl [-]", nch, "Cd [-]");
                    for(PetscInt pti=0; pti<n_afpts; pti++)
                    {
                        PetscReal alpha = farm->wt[i]->foils[j]->aoa[pti];
                        PetscReal cl    = farm->wt[i]->foils[j]->cl[pti];
                        PetscReal cd    = farm->wt[i]->foils[j]->cd[pti];
                        PetscPrintf(mesh->MESH_COMM,"    | %*.4f | %*.4f | %*.4f|\n", nch, alpha, nch, cl, nch, cd);
                    }
                }

                PetscPrintf(mesh->MESH_COMM,"\n");
                PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");

                // controller properties
                if(farm->wt[i]->genControllerType != "none")
                {
                    PetscPrintf(mesh->MESH_COMM," =============================================================================\n\n");
                    PetscPrintf(mesh->MESH_COMM," Generator torque controller:\n");
                    PetscPrintf(mesh->MESH_COMM, "  genControllerType : %s\n", farm->wt[i]->genControllerType.c_str());
                    PetscPrintf(mesh->MESH_COMM, "  rtrOmega          : %lf\n", farm->wt[i]->rtrOmega);
                    PetscPrintf(mesh->MESH_COMM, "  genOmega          : %lf\n", farm->wt[i]->genOmega);
                    PetscPrintf(mesh->MESH_COMM, "  rtrOmegaFilt      : %lf\n", farm->wt[i]->rtrOmegaFilt);
                    PetscPrintf(mesh->MESH_COMM, "  rtrSpdFilterFreq  : %lf\n", farm->wt[i]->rtrSpdFilterFreq);
                    PetscPrintf(mesh->MESH_COMM, "  cutInGenSpd       : %lf\n", farm->wt[i]->cutInGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  cutInGenTq        : %lf\n", farm->wt[i]->cutInGenTq);
                    PetscPrintf(mesh->MESH_COMM, "  regTwoStartGenSpd : %lf\n", farm->wt[i]->regTwoStartGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  regTwoEndGenSpd   : %lf\n", farm->wt[i]->regTwoEndGenSpd);
                    PetscPrintf(mesh->MESH_COMM, "  ratedGenTq        : %lf\n", farm->wt[i]->ratedGenTq);
                    PetscPrintf(mesh->MESH_COMM, "  omegaKP           : %lf\n", farm->wt[i]->omegaKP);
                    PetscPrintf(mesh->MESH_COMM, "  genTorque         : %lf\n", farm->wt[i]->genTorque);
                    PetscPrintf(mesh->MESH_COMM, "  genPwr            : %lf\n", farm->wt[i]->genPwr);
                    PetscPrintf(mesh->MESH_COMM, "  tqRateLimiter     : %ld\n", farm->wt[i]->tqRateLimiter);
                    PetscPrintf(mesh->MESH_COMM, "  rtrSpdLimiter     : %ld\n", farm->wt[i]->rtrSpdLimiter);
                    PetscPrintf(mesh->MESH_COMM, "  tqMaxRate         : %lf\n", farm->wt[i]->tqMaxRate);
                    PetscPrintf(mesh->MESH_COMM, "  ratedRotorSpd     : %lf\n", farm->wt[i]->ratedRotorSpd);
                    PetscPrintf(mesh->MESH_COMM, "  driveTrainInertia : %lf\n", farm->wt[i]->driveTrainInertia);
                    PetscPrintf(mesh->MESH_COMM, "  genInertia        : %lf\n", farm->wt[i]->genInertia);
                    PetscPrintf(mesh->MESH_COMM, "  hubInertia        : %lf\n", farm->wt[i]->hubInertia);
                    PetscPrintf(mesh->MESH_COMM, "  bldInertia        : %lf\n", farm->wt[i]->bldInertia);
                    PetscPrintf(mesh->MESH_COMM, "  gbxRatioG2R       : %lf\n", farm->wt[i]->gbxRatioG2R);
                    PetscPrintf(mesh->MESH_COMM, "  gbxEff            : %lf\n", farm->wt[i]->gbxEff);

                }
            }

            PetscPrintf(mesh->MESH_COMM,"\n");
        }
    }

    MPI_Barrier(mesh->MESH_COMM);

    return(0);
}

//***************************************************************************************************************//
