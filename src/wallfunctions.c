#include "include/base.h"
#include "include/domain.h"
#include "include/inline.h"
#include "include/wallfunctions.h"

PetscInt     pre_integrate_flag = 0;
PetscReal    *integration_buffer;
PetscReal    *integration_buffer_rough;

//***************************************************************************************************************//

void uStarShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &z0,
    PetscReal &gammaM, PetscReal &kappa, PetscReal &qwall, PetscReal &thetaRef,
    PetscReal &uStar, PetscReal &phiM, PetscReal &L, PetscReal nu
)
{
    PetscReal uStar0 = (kappa * UParallelMeanMag) / std::log(wallDist / z0);
    PetscReal uStar1 = uStar0 + 1e-5;
    PetscReal f0, f1;

    PetscReal g      = 9.81;

    PetscInt  iter = 0, iterMax = 100;
    PetscReal tol  = 1e-5;

    if(qwall == 0.0)
    {
        uStar = uStar0;
        L     = 1e30;
        phiM  = 1.0;
        return;
    }
    else if(qwall < 0.0)
    {
        do
        {
            // set initial guesses at L
            PetscReal L0 = -(std::pow(uStar0,3.0)) / (kappa * (g/thetaRef) * qwall);
            PetscReal L1 = -(std::pow(uStar1,3.0)) / (kappa * (g/thetaRef) * qwall);

            // limit L to always be positive and finite
            L0 = std::max(1e-10,L0);
            L1 = std::max(1e-10,L1);

            // form the "zeta" variable
            PetscReal zeta0 = wallDist/L0;
            PetscReal zeta1 = wallDist/L1;

            // form psiM
            PetscReal psiM0 = -gammaM * zeta0;
            PetscReal psiM1 = -gammaM * zeta1;

            // form the function that we are driving to zero
            PetscReal denom0 = std::max(std::log(wallDist/z0) - psiM0, 1e-10);
            PetscReal denom1 = std::max(std::log(wallDist/z0) - psiM1, 1e-10);

            f0 = uStar0 - ((UParallelMeanMag * kappa) / denom0);
            f1 = uStar1 - ((UParallelMeanMag * kappa) / denom1);

            // update uStar
            PetscReal slope = std::max((f1 - f0) / (uStar1 - uStar0), 1e-10);
            PetscReal uStarTmp = uStar1;
            uStar1 = std::max(0.0, uStar0 - (f0 / slope));
            uStar0 = uStarTmp;

            uStar  = uStar1;
            L      = L1;
            phiM   = 1.0 + (gammaM * zeta1);

            iter++;

        } while(fabs(f1)>tol && uStar0 != uStar1 && iter < iterMax);

        if(iter == iterMax)
        {
            printf("max uStar iteration reached\n");
        }

        // if uStar is close zero it means that the zero does not exist, but a maximum exists
        if(uStar < 1e-2)
        {
            PetscReal duStar = 0.5 / 100.0;
            uStar0 = 1e-10;
            uStar1 = uStar0 + duStar;

            do
            {
                uStar0 = uStar1;
                uStar1 = uStar0 + duStar;

                // set initial guesses at L
                PetscReal L0 = -(std::pow(uStar0,3.0)) / (kappa * (g/thetaRef) * qwall);
                PetscReal L1 = -(std::pow(uStar1,3.0)) / (kappa * (g/thetaRef) * qwall);

                // limit L to always be positive and finite
                L0 = std::max(1e-10,L0);
                L1 = std::max(1e-10,L1);

                // form the "zeta" variable
                PetscReal zeta0 = wallDist/L0;
                PetscReal zeta1 = wallDist/L1;

                // form psiM
                PetscReal psiM0 = -gammaM * zeta0;
                PetscReal psiM1 = -gammaM * zeta1;

                // form the function that we are driving to zero
                f0 = UParallelMeanMag - (uStar0/kappa)*(std::log(wallDist/z0) - psiM0);
                f1 = UParallelMeanMag - (uStar1/kappa)*(std::log(wallDist/z0) - psiM1);

                uStar = uStar0;
                L     = L0;
                phiM  = 1.0 + (gammaM * zeta0);

            } while (f1 > f0);

            printf("uStar is zero, looking for maximum: uStar = %lf\n", uStar);
        }

        return;
    }
    else
    {
        char error[512];
        sprintf(error, "wall model for unstable flow not yet implemented");
        fatalErrorInFunction("uStarShumann", error);
    }
}

//***************************************************************************************************************//

void qWallShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &z0,
    PetscReal &gammaM, PetscReal &gammaH, PetscReal &alphaH,
    PetscReal &thetaRef, PetscReal &deltaTheta, PetscReal &kappa,
    PetscReal &qWall, PetscReal &uStar, PetscReal &phiM, PetscReal &phiH, PetscReal &L
)
{
    PetscReal uStar0 = std::max(0.0, (kappa * UParallelMeanMag) / std::log(wallDist / z0));
    PetscReal qWall0 = (-deltaTheta * uStar0 * kappa) / (alphaH * std::log(wallDist / z0));

    if(fabs(qWall0) < 1e-6)
    {
        qWall0 = signNonZero(qWall) * 1e-6;
    }

    PetscReal g      = 9.81;

    PetscInt  iter = 0, iterMax = 100;
    PetscReal tol  = 1e-5;

    if(deltaTheta == 0.0)
    {
        uStar = uStar0;
        qWall = qWall0;
        L     = 1e30;
        phiM  = 1.0;
        phiH  = alphaH;
        return;
    }
    else if(deltaTheta > 0.0)
    {
        PetscReal uStarDiff = 1e30;
        PetscReal qWallDiff = 1e30;

        do
        {
            // set initial guesses at L
            PetscReal L0 = -(std::pow(uStar0,3.0)) / (kappa * (g/thetaRef) * qWall0);

            // limit L to always be positive and finite
            L0 = std::max(1e-10,L0);

            // form the "zeta" variable
            PetscReal zeta0 = wallDist/L0;

            // form psiM
            PetscReal psiM0 = -gammaM * zeta0;

            // form psiH
            PetscReal psiH0 = -gammaH * zeta0;

            // update uStar
            PetscReal uStarOld = uStar0;
            uStar0 = std::max(0.0, (kappa * UParallelMeanMag) / (std::log(wallDist / z0) - psiM0));

            // update qWall
            PetscReal qWallOld = qWall0;
            qWall0 = (-deltaTheta * uStar0 * kappa) / (alphaH * std::log(wallDist / z0) - psiH0);

            // compute change in uStar and qWall
            uStarDiff = uStar0 - uStarOld;
            qWallDiff = qWall0 - qWallOld;

            // evaluate in case it is converged
            uStar  = uStar0;
            qWall  = qWall0;
            L      = L0;
            phiM   = 1.0 + (gammaM * zeta0);
            phiH   = alphaH + (gammaH * zeta0);

            iter++;

        } while((fabs(uStarDiff)>tol || abs(qWallDiff)>tol) && iter < iterMax);

        if(iter == iterMax)
        {
            printf("max uStar iteration reached\n");
        }

        return;
    }
    else
    {
        char error[512];
        sprintf(error, "wall model for unstable flow not yet implemented");
        fatalErrorInFunction("qWallShumann", error);
    }
}

//***************************************************************************************************************//

inline PetscReal near_wall_eddy_viscosity(PetscReal yplus)
{
    PetscReal   kappa = 0.41;
    return kappa * yplus * pow(1. - exp(-yplus / 19.), 2.0);
}

//***************************************************************************************************************//

inline PetscReal near_wall_eddy_viscosity(PetscReal yplus, PetscReal yp_shift)
{
    PetscReal   kappa = 0.41;
	return kappa * (yplus+yp_shift) * pow ( 1. - exp( - (yplus+yp_shift) / 19.0 ) , 2.0 );
};

//***************************************************************************************************************//

inline void pre_integrate
(
    PetscInt     &n_yp,
    PetscInt     &interval_yp,
    PetscInt     &max_yp
)
{
    if (pre_integrate_flag)
        return;
    else
        pre_integrate_flag = 1;

    n_yp = (max_yp / interval_yp);

    integration_buffer       = new PetscReal [ n_yp + 1 ];
    integration_buffer_rough = new PetscReal [ n_yp + 1 ];

    integration_buffer[0] = 0.;
    integration_buffer_rough[0] = 0.;

    for (PetscInt i = 1; i <= n_yp; i++)
    {
        PetscInt N = 24;

        PetscReal ya_p = (PetscReal) (i - 1) * interval_yp;
        PetscReal yb_p = (PetscReal) (i) * interval_yp;
        PetscReal ydiff = yb_p - ya_p, dy_p = ydiff / (PetscReal) N;

        std::vector<PetscReal> E(N+1);

        PetscReal val = 0, ybegin_p = ya_p;

        for (PetscInt k = 0; k <= N; k++)
        {
            E[k] = 1. / (1. + near_wall_eddy_viscosity(ybegin_p + dy_p * k));
        }

        for (PetscInt k = 0; k < N; k++)
        {
            PetscReal F[4];
            F[0] = E[k];
            F[1] = 1. /
            (
                1. + near_wall_eddy_viscosity
                (ybegin_p + dy_p * 1. / 3.)
            );

            F[2] = 1. /
            (
                1. + near_wall_eddy_viscosity
                (ybegin_p + dy_p * 2. / 3.)
            );

            F[3] = E[k + 1];

            val += dy_p / 3. * (3 * F[0] + 9 * F[1] + 9 * F[2] + 3 * F[3]) / 8.;

            ybegin_p += dy_p;
        }

        integration_buffer[i] = integration_buffer[i - 1] + val;
        integration_buffer_rough[i] = integration_buffer_rough[i - 1] + val;
    }

    return;
}

//***************************************************************************************************************//

inline PetscReal integrate_F(PetscReal nu, PetscReal utau, PetscReal yb)
{
    PetscReal val = 0;

    PetscInt     n_yp = 0;
    PetscInt     interval_yp = 2;
    PetscInt     max_yp = 1e7;

    pre_integrate
    (
        n_yp,
        interval_yp,
        max_yp
    );

    PetscReal ya_plus = 0 * utau / nu;
    PetscReal yb_plus = yb * utau / nu;

    PetscInt ib = (PetscInt) (yb_plus / (PetscReal) interval_yp);
    PetscInt N = 4;

    if (yb_plus <= (PetscReal) max_yp)
    {
        PetscReal int_a = 0;
        PetscReal int_b = (integration_buffer[ib + 1] - integration_buffer[ib])
                        / (PetscReal) (interval_yp)
                        * (yb_plus - (PetscReal) (ib) * interval_yp)
                        + integration_buffer[ib];

        return (int_b - int_a) / utau;
    }
    else
    {
        val = integration_buffer[n_yp];
        ya_plus = max_yp;
    }

    PetscReal ydiff = yb_plus - ya_plus, dy = ydiff / (PetscReal) N;

    std::vector<PetscReal> E(N+1);
    PetscReal ybegin = ya_plus;

    for (PetscInt i = 0; i <= N; i++)
    {
        E[i] = 1. / (1. + near_wall_eddy_viscosity(ybegin + dy * i));
    }

    for (PetscInt i = 0; i < N; i++)
    {
        PetscReal F[4];
        F[0] = E[i];
        F[1] = 1. / (1. + near_wall_eddy_viscosity(ybegin + dy * 1. / 3.));
        F[2] = 1. / (1. + near_wall_eddy_viscosity(ybegin + dy * 2. / 3.));
        F[3] = E[i + 1];
        val += dy / 3. * (3 * F[0] + 9 * F[1] + 9 * F[2] + 3 * F[3]) / 8.;
        ybegin += dy;
    }

    val /= utau;
    return val;
}

//***************************************************************************************************************//

inline PetscReal integrate_F(PetscReal nu, PetscReal utau, PetscReal yb, PetscReal ks)
{
	PetscReal ks_plus = utau * ks / nu;
  	PetscReal yp_shift = 0.9*(sqrt(ks_plus)-(ks_plus)*exp(-ks_plus/6.0));

	if(yp_shift<0) yp_shift = 0;

    PetscReal val = 0;

	PetscReal ya_plus = 0 * utau / nu;
	PetscReal yb_plus = yb * utau / nu;

	int N=25;
	val=0;

	PetscReal ydiff = yb_plus - ya_plus, dy = ydiff / (PetscReal)N;
	std::vector<PetscReal> E(N+1);
	PetscReal ybegin=ya_plus;

	for(PetscInt i=0; i<=N; i++)
    {
		E[i] =  1./ ( 1. + near_wall_eddy_viscosity(ybegin + dy*i, yp_shift ) );
	}

	for(PetscInt i=0; i<N; i++)
    {
	    PetscReal F[4];
	    F[0] = E[i];
		F[1] = 1./ ( 1. + near_wall_eddy_viscosity(ybegin+dy*1./3., yp_shift) );
		F[2] = 1./ ( 1. + near_wall_eddy_viscosity(ybegin+dy*2./3., yp_shift) );
		F[3] = E[i+1];
		val += dy / 3.* ( 3 * F[0] + 9 * F[1] + 9 * F[2] + 3 *F[3] ) / 8.;
		ybegin += dy;
	}

	val /= utau;
	return val;
};

//***************************************************************************************************************//

inline PetscReal f_Cabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    return (utau * utau * integrate_F(nu, utau, y) - u);
}

//***************************************************************************************************************//

inline PetscReal f_Cabot_roughness(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn, PetscReal ks)
{
	return utau * utau * integrate_F( nu, utau, y, ks ) - u;
}

//***************************************************************************************************************//

inline PetscReal df_Cabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    PetscReal eps = 1.e-7;

    return
    (
        f_Cabot(nu, u, y, utau + eps, dpdn) -
        f_Cabot(nu, u, y, utau - eps, dpdn)
    ) / (2.0 * eps);
}

//***************************************************************************************************************//

inline PetscReal df_Cabot_roughness(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn, PetscReal ks)
{
	PetscReal eps=1.e-7;
	return ( f_Cabot_roughness(nu, u, y, utau+eps, dpdn, ks) - f_Cabot_roughness (nu, u, y, utau-eps, dpdn, ks) ) / ( 2*eps ) ;
}

//***************************************************************************************************************//

inline PetscReal u_Cabot(PetscReal nu, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    return utau * utau * integrate_F(nu, utau, y); // + dpdn * integrate_Fy( nu, utau, 0, y );
}

//***************************************************************************************************************//

inline PetscReal u_Cabot_roughness(PetscReal nu, PetscReal y, PetscReal utau, PetscReal dpdn, PetscReal ks)
{
  return utau * utau * integrate_F( nu, utau, y, ks );// + dpdn * integrate_Fy( nu, utau, 0, y );
}

//***************************************************************************************************************//

PetscReal uTauCabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn)
{
    PetscReal x, x0 = guess;

    PetscInt i;

    for (i = 0; i < 30; i++)
    {
        x = x0 - f_Cabot(nu, u, y, x0, dpdn) / df_Cabot(nu, u, y, x0, dpdn);
        if (fabs(x0 - x) < 1.e-10)
            break;
        x0 = x;
    }

    if (fabs(x0 - x) > 1.e-5 && i >= 29)
    {
        char wrngMsg[256];
        sprintf(wrngMsg, "Cabot newton iteration failed\n");
        warningInFunction("uTauCabot", wrngMsg);
    }

    return x;
}

//***************************************************************************************************************//

PetscReal uTauCabotRoughness(PetscReal nu, PetscReal u, PetscReal y, PetscReal guess, PetscReal dpdn, PetscReal ks)
{
	PetscReal x, x0=guess;

	int i;

	for(i=0; i<30; i++) {
		x = x0 - f_Cabot_roughness(nu, u, y, x0, dpdn, ks)/df_Cabot_roughness(nu, u, y, x0, dpdn, ks);
		if( fabs(x0 - x) < 1.e-10) break;
		x0 = x;
	}

	if( fabs(x0 - x) > 1.e-5 && i>=29 )
    {
        char wrngMsg[256];
        sprintf(wrngMsg, "Cabot newton iteration failed\n");
        warningInFunction("uTauCabotRoughness", wrngMsg);
    }
	return x;
};

//***************************************************************************************************************//

PetscReal utau_wf(PetscReal nu, PetscReal ks, PetscReal sb, PetscReal Ut_mag)
{
    PetscReal _ks = ks * 0.033;
    PetscReal   kappa = 0.41;

    if (_ks > 1.e-13)
    {
        return Ut_mag * kappa / log(sb / _ks);
    }
    else
    {
        PetscReal A = 8.3;
        PetscReal B = 1.0 / 7.0;
        return pow(Ut_mag * pow(nu, B) / (A * pow(sb, B)), 1.0 / (1.0 + B)); ///

    }
}

//***************************************************************************************************************//

void wallFunctionCabot(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc,
        Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{

    PetscReal u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
    PetscReal un = u_c * nf.x + v_c * nf.y + w_c * nf.z;
    PetscReal ut = u_c - un * nf.x, vt = v_c - un * nf.y, wt = w_c - un * nf.z;
    PetscReal ut_mag = sqrt(ut * ut + vt * vt + wt * wt);
    PetscReal ut_mag_modeled;

    *ustar = uTauCabot(nu, ut_mag, sc, 0.01, 0);
    ut_mag_modeled = u_Cabot(nu, sb, *ustar, 0);

    PetscReal Ua_abs = sqrt(Ua.x * Ua.x + Ua.y * Ua.y + Ua.z * Ua.z);

    if (ut_mag > 1.e-10) {
        ut *= ut_mag_modeled / ut_mag;
        vt *= ut_mag_modeled / ut_mag;
        wt *= ut_mag_modeled / ut_mag;
    }
    else
    {
        ut = vt = wt = 0;
    }

    // u = ut + (u.n)n
    (*Ub).x = ut + sb / sc * un * nf.x;
    (*Ub).y = vt + sb / sc * un * nf.y;
    (*Ub).z = wt + sb / sc * un * nf.z;

    (*Ub).x += Ua.x;
    (*Ub).y += Ua.y;
    (*Ub).z += Ua.z;

    // (*Ub).x = (sb/sc) * Uc.x + (1.0 - (sb/sc)) * Ua.x;
    // (*Ub).y= (sb/sc) * Uc.y + (1.0 - (sb/sc)) * Ua.y;
    // (*Ub).z = (sb/sc) * Uc.z + (1.0 - (sb/sc)) * Ua.z;

    return;
}

//***************************************************************************************************************//

void wallFunctionCabotRoughness(PetscReal nu, PetscReal ks, PetscReal sc, PetscReal sb,
    Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    PetscReal u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
    PetscReal un = u_c * nf.x + v_c * nf.y + w_c * nf.z;
    PetscReal ut = u_c - un * nf.x, vt = v_c - un * nf.y, wt = w_c - un * nf.z;
    PetscReal ut_mag = sqrt(ut * ut + vt * vt + wt * wt);

	*ustar = uTauCabotRoughness(nu, ut_mag, sc, 0.01, 0, ks);
	PetscReal ut_mag_modeled = u_Cabot_roughness(nu, sb, *ustar, 0, ks);

	if(ut_mag>1.e-10) {
		ut *= ut_mag_modeled/ut_mag;
		vt *= ut_mag_modeled/ut_mag;
		wt *= ut_mag_modeled/ut_mag;
	}
	else ut=vt=wt=0;

	// u = ut + (u.n)n
    (*Ub).x = ut + sb / sc * un * nf.x;
    (*Ub).y = vt + sb / sc * un * nf.y;
    (*Ub).z = wt + sb / sc * un * nf.z;

	(*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;

    // (*Ub).x = (sb/sc) * Uc.x + (1.0 - (sb/sc)) * Ua.x;
    // (*Ub).y= (sb/sc) * Uc.y + (1.0 - (sb/sc)) * Ua.y;
    // (*Ub).z = (sb/sc) * Uc.z + (1.0 - (sb/sc)) * Ua.z;


    return;
}

//***************************************************************************************************************//
// Wall function based on the log law model
void wallFunctionSchumann(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness,
    PetscReal kappa, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    unc  = nScale(nDot(u_c, nf), nf);
    Cmpnts    utc  = nSub(u_c, unc);
    Cmpnts    ut_b;

    Cmpnts    et  = nUnit(utc);

    PetscReal ut_mag = nMag(utc), ut_magb;

    *ustar = ut_mag * kappa / log(sc/roughness);

    // ut_magb = (*ustar/kappa) * log(sb/roughness);
    
    // if (ut_magb < 1.e-10)
    // {
    //     ut_magb = 0.0;
    // }
    
    // ut_b  = nScale(ut_magb, et);
    
    // (*Ub) = nSum(ut_b, nScale( (sb/sc), unc));
    
    // mSum(*Ub, Ua);

    (*Ub).x = (sb/sc) * Uc.x + (1.0 - (sb/sc)) * Ua.x;
    (*Ub).y= (sb/sc) * Uc.y + (1.0 - (sb/sc)) * Ua.y;
    (*Ub).z = (sb/sc) * Uc.z + (1.0 - (sb/sc)) * Ua.z;

    return;
}

//***************************************************************************************************************//

// Wall function based on the power law model
void wallFunctionPowerlaw(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua,
        Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    PetscReal u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
    PetscReal un = u_c * nf.x + v_c * nf.y + w_c * nf.z;
    PetscReal ut = u_c - un * nf.x, vt = v_c - un * nf.y, wt = w_c - un * nf.z;
    PetscReal ut_mag = sqrt(ut * ut + vt * vt + wt * wt);

    PetscReal A = 8.3;
    PetscReal B = 1.0 / 7.0;
    *ustar = pow(ut_mag * pow(nu, B) / (A * pow(sc, B)), 1.0 / (1.0 + B));  ///

    PetscReal ybp = sb * (*ustar) / nu;
    PetscReal ycp = sc * (*ustar) / nu;
    PetscReal ut_mag_modeled;
    if (ybp > 12.0)
    {
        ut_mag_modeled = A * (*ustar) * pow(ybp, B);
    }
    else
    {
        ut_mag_modeled = ut_mag * sb / sc;
    }

    if (ut_mag > 1.e-10)
    {
        ut *= ut_mag_modeled / ut_mag;
        vt *= ut_mag_modeled / ut_mag;
        wt *= ut_mag_modeled / ut_mag;
    }
    else
    {
        ut = vt = wt = 0;
    }

    // u = ut + (u.n)n
    (*Ub).x = ut + sb / sc * un * nf.x;
    (*Ub).y = vt + sb / sc * un * nf.y;
    (*Ub).z = wt + sb / sc * un * nf.z;

    (*Ub).x += Ua.x;
    (*Ub).y += Ua.y;
    (*Ub).z += Ua.z;

    return;
}

//***************************************************************************************************************//

// Wall function - Adverse pressure gradient Power Law (APGPL) - A New Explicit Algebraic Wall Model for LES of Turbulent Flows
// Under Adverse Pressure Gradient, Sylvia Wilhelm, Pierre Sagaut - 2020
void wallFunctionPowerlawAPG(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness, PetscReal kappa, Cmpnts Ua,
        Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf, PetscReal dpdx, PetscReal dpdy, PetscReal dpdz)
{
    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    un  = nScale(nDot(u_c, nf), nf);
    Cmpnts    ut  = nSub(u_c, un), ut_b;
    Cmpnts    et  = nUnit(ut);

    PetscReal ut_mag = nMag(ut), ut_magb;

    PetscReal A = 8.3;
    PetscReal B = 1.0 / 7.0;
    PetscReal dpds = dpdx*et.x + dpdy*et.y + dpdz*et.z;               //the pressure gradient in the tangential direction to the local body normal

    PetscReal ustarNoslip = sqrt(nu*ut_mag/sc);
    PetscReal ybp, ycp;

    if(dpds <= 0)       //favourable pressure gradient
    {
        if(roughness > 1.e-9)
        {
            *ustar = ut_mag * kappa / log(sc/roughness);
        }
        else
        {
            *ustar = pow(ut_mag * pow(nu, B) / (A * pow(sc, B)), 1.0 / (1.0 + B));
        }

        ybp = sb * std::max((*ustar), ustarNoslip)/nu;
    	ycp = sc * std::max((*ustar), ustarNoslip)/nu;

        if(ycp>12.0)
        {
            if(roughness > 1.e-9)
            {
                ut_magb = (*ustar) * log(sb/roughness)/kappa;
            }
            else
            {
                ut_magb = A * (*ustar) * pow(ybp, B);
            }
        }
        else
        {
            ut_magb = ut_mag * sb/sc;
            *ustar = ustarNoslip;
        }
    }
    else //adverse pressure gradient
    {
        PetscReal D;
        PetscReal alpha = 7.5789;
        PetscReal beta  = -1.4489;
        PetscReal gamma = 191.1799;

        // D determines regions of separation and reattachment
        D = fabs(ut_mag - alpha * sqrt (sc * dpds) - beta * pow(nu * dpds, 1.0/3.0) * log (gamma * pow(sc, 3.0) * dpds / pow(nu, 2.0) ) );

        ycp = sc * ustarNoslip / nu;

        if (D >= 0) // use Adverse pressure gradient power law
        {
            if(ycp>12.0)
            {
                *ustar  = pow(A, -1.0 / (1.0 + B)) * pow( (nu/sc), B / (1.0 + B)) * pow(D, 1.0 / (1.0 + B));

                ut_magb = A * pow( (*ustar), B + 1.0) * pow( (sb/nu), B) + (ut_mag-D);
            }
            else
            {
                ut_magb = ut_mag * sb/sc;;
                *ustar = ustarNoslip;
            }
        }
        else //separation or reattachment zone - no wall model used
        {
            ut_magb = ut_mag * sb/sc;;
            *ustar = ustarNoslip;
        }
    }

    if (ut_magb < 1.e-10)
    {
        ut_magb = 0.0;
    }

    ut_b  = nScale(ut_magb, et);

    (*Ub) = nSum(ut_b, nScale( (sb/sc), un));

    mSum(*Ub, Ua);

    return;
}


//***************************************************************************************************************//

void wallFunctionLogLawAPG(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness, PetscReal kappa, Cmpnts Ua,
        Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf, PetscReal dpdx, PetscReal dpdy, PetscReal dpdz)
{
    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    un  = nScale(nDot(u_c, nf), nf);
    Cmpnts    ut  = nSub(u_c, un), ut_b;
    Cmpnts    et  = nUnit(ut);

    PetscReal ut_mag = nMag(ut), ut_magb;

    PetscReal A = 8.3;
    PetscReal B = 1.0 / 7.0;
    PetscReal dpds = dpdx*et.x + dpdy*et.y + dpdz*et.z;               //the pressure gradient in the tangential direction to the local body normal

    PetscReal ustarNoslip = sqrt(nu*ut_mag/sc);
    PetscReal ybp, ycp;

    if(dpds <= 0)       //favourable pressure gradient
    {
        if(roughness > 1.e-9)
        {
            *ustar = ut_mag * kappa / log(sc/roughness);
        }
        else
        {
            *ustar = pow(ut_mag * pow(nu, B) / (A * pow(sc, B)), 1.0 / (1.0 + B));
        }

        ybp = sb * std::max((*ustar), ustarNoslip)/nu;
    	ycp = sc * std::max((*ustar), ustarNoslip)/nu;

        if(ycp>12.0)
        {
            if(roughness > 1.e-9)
            {
                ut_magb = (*ustar) * log(sb/roughness)/kappa;
            }
            else
            {
                ut_magb = A * (*ustar) * pow(ybp, B);
            }
        }
        else
        {
            ut_magb = ut_mag * sb/sc;
            *ustar = ustarNoslip;
        }
    }
    else //adverse pressure gradient
    {

        ycp = sc * ustarNoslip / nu;
        ybp = sb * ustarNoslip / nu;

        if(ycp>12.0)
        {
            PetscReal up = pow(nu * dpds, 1.0/3.0);
            PetscReal lambda = 0.0;

            if(roughness > 1.e-9)
            {
                *ustar = uTauCabotRoughness(nu, ut_mag, sc, 0.01, 0, roughness);

                if(*ustar < 1e-5)
                {
                    ut_magb = (*ustar) * log(sb/roughness)/kappa;
                }
                else
                {
                    lambda = pow(up/ (*ustar), 3.0);

                    PetscReal ct = pow(1.0 + lambda * (sb/roughness), 0.5);

                    ut_magb = ((*ustar)/kappa) * (log(sb/roughness) - 2.0 * log((ct+1.0)/2.0) + 2.0 * (ct - 1.0));
                }
            }
            else
            {
                *ustar = uTauCabot(nu, ut_mag, sc, 0.01, 0);

                if(*ustar < 1e-5)
                {
                    ut_magb = (*ustar) * log(ybp)/kappa;
                }
                else
                {
                    lambda = pow( up/(*ustar), 3.0);

                    PetscReal ct = pow(1.0 + lambda * ybp, 0.5);

                    ut_magb = ((*ustar)/kappa) * (log(ybp) - 2.0 * log((ct+1.0)/2.0) + 2.0 * (ct - 1.0));
                }
            }

        }
        else
        {
            ut_magb = ut_mag * sb/sc;;
            *ustar = ustarNoslip;
        }

    }

    if (ut_magb < 1.e-10)
    {
        ut_magb = 0.0;
    }

    ut_b  = nScale(ut_magb, et);

    (*Ub) = nSum(ut_b, nScale( (sb/sc), un));

    mSum(*Ub, Ua);

    return;
}

void slipBC(PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, Cmpnts nf)
{

    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    un  = nScale(nDot(u_c, nf), nf);
    Cmpnts    ut  = nSub(u_c, un);

    *Ub = nSum(ut, nScale(sb/sc, un));

    mSum(*Ub, Ua);

	return;
}
