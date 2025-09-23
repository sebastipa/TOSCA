#include "turbines_force.h"

//***************************************************************************************************************//

PetscErrorCode computeBladeForce(farm_ *farm)
{
    // turbine and mesh point indices
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
                        findInterpolationWeights(w, l, an1, s1, wt->adm.alpha[p]);
                        PetscReal cl_1 = w[0]*cl1[l[0]] + w[1]*cl1[l[1]];
                        PetscReal cd_1 = w[0]*cd1[l[0]] + w[1]*cd1[l[1]];

                        // interpolate coeffs for the second airfoil
                        findInterpolationWeights(w, l, an2, s2, wt->adm.alpha[p]);
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
                        findInterpolationWeights(w, l, an1, s1, wt->alm.alpha[p]);
                        PetscReal cl_1 = w[0]*cl1[l[0]] + w[1]*cl1[l[1]];
                        PetscReal cd_1 = w[0]*cd1[l[0]] + w[1]*cd1[l[1]];

                        // interpolate coeffs for the second airfoil
                        findInterpolationWeights(w, l, an2, s2, wt->alm.alpha[p]);
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
                // second test inside because this variable is not defined for UADM and AFM
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
                Cmpnts point_p    = wt->afm.point;

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