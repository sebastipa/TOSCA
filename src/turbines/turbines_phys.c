#include "turbines_phys.h"

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

        if((*farm->turbineModels[t]) == "ALM" || (*farm->turbineModels[t]) == "ADM")
        {
            // test if this processor controls this turbine
            if(wt->turbineControlled && !wt->useOpenFAST)
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
                    PetscReal angle        = clock->dt * wt->rtrOmega;
                    PetscInt updateAzimuth = 1;
                    rotateBlades(wt, angle, updateAzimuth);
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

                // sync turbines on overall master node
                if(!rank)
                {
                    PetscReal angle        = (gazimuth - wt->alm.azimuth) * wt->deg2rad;
                    PetscInt updateAzimuth = 1;
                    rotateBlades(wt, angle, updateAzimuth);
                }
            }
        }
    }

    return(0);
}

//***************************************************************************************************************//

PetscErrorCode rotateBlades(windTurbine *wt, PetscReal angle, PetscInt updateAzimuth)
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

    if(updateAzimuth)
    {
        wt->alm.azimuth += angle * wt->rad2deg;

        // bound azimuth between 360 and 0
        if(wt->alm.azimuth >= 360.0)
        {
            wt->alm.azimuth -= 360.0;
        }
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

        if((*farm->turbineModels[t]) == "ALM" || (*farm->turbineModels[t]) == "ADM")
        {
            // test if this processor controls this turbine
            if(wt->turbineControlled && !wt->useOpenFAST)
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
            if(wt->turbineControlled && !wt->useOpenFAST) 
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
    //       This strongly reduces the parallel communications, performed only on the turbine controlling processors.

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

        // test if this turbine has yaw control (either via TOSCA yaw controller or OpenFAST)
        if(wt->yawControllerType != "none" || wt->useOpenFAST)
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
            if(wt->yawControllerType != "none" && !wt->useOpenFAST)
            {
                // re-initialize moved flag if yaw control is active
                wt->trbMoved = 0; 

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
                            mRot(wt->twrDir, wt->rtrDir,    angle);
                            mRot(wt->twrDir, wt->rtrAxis,   angle);
                            mRot(wt->twrDir, wt->rotCenter, angle);

                            // actuator disk model
                            if((*farm->turbineModels[t]) == "ADM")
                            {
                                // number of points in the AD mesh
                                PetscInt npts_t = wt->adm.nPoints;

                                // loop over the AD mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(wt->adm.points[p], wt->twrTop);
                                    mRot(wt->twrDir, wt->adm.points[p], angle);
                                    mSum(wt->adm.points[p], wt->twrTop);
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
                                    mSub(wt->uadm.points[p], wt->twrTop);
                                    mRot(wt->twrDir, wt->uadm.points[p], angle);
                                    mSum(wt->uadm.points[p], wt->twrTop);
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
                                    mSub(wt->alm.points[p], wt->twrTop);
                                    mRot(wt->twrDir, wt->alm.points[p], angle);
                                    mSum(wt->alm.points[p], wt->twrTop);
                                }

                                // rotate other parameters
                                mRot(wt->twrDir, wt->omega_hat, angle);
                            }
                            else if((*farm->turbineModels[t]) == "AFM")
                            {
                                // nothing to do (nothing to rotate)
                            }

                            // rotate up-sampling points
                            {
                                upSampling *upPoints = farm->wt[t]->upPoints;

                                // rotate center
                                mSub(upPoints->center, wt->twrTop);
                                mRot(wt->twrDir, upPoints->center, angle);
                                mSum(upPoints->center, wt->twrTop);

                                // number of points in the sample mesh
                                PetscInt npts_t = upPoints->nPoints;

                                // loop over the sample mesh points
                                for(p=0; p<npts_t; p++)
                                {
                                    mSub(upPoints->points[p], wt->twrTop);
                                    mRot(wt->twrDir, upPoints->points[p], angle);
                                    mSum(upPoints->points[p], wt->twrTop);
                                }
                            }
                        }

                        // turbine point search will be performed
                        wt->trbMoved = 1;
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
            if(wt->yawControllerType != "none" && !wt->useOpenFAST)
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
        if(wt->turbineControlled && !wt->useOpenFAST)
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