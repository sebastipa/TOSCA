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
    PetscReal &uStar, PetscReal &phiM, PetscReal &L,
    PetscInt k, PetscInt j, PetscInt i
)
{
    if (wallDist <= z0) 
    { 
        char error[512];
        sprintf(error, "wall distance is smaller than the roughness length for cell (%ld %ld %ld)\n", k, j, i);
        fatalErrorInFunction("uStarShumann",  error); 
    }

    PetscReal logLaw = std::log(wallDist / z0);
    PetscReal uStar_neutral = std::max(0.0, (kappa * UParallelMeanMag) / logLaw);

    if (fabs(qwall) < 1e-6) {
        uStar = uStar_neutral;
        L = 1e30;
        phiM = 1.0;
        return;
    }

    PetscReal g = 9.81;
    PetscInt iter = 0, iterMax = 500;
    PetscReal tol = 1e-5;
    PetscReal omega = 0.5;  // Under-relaxation factor

    // Start with neutral for iteration
    PetscReal uStar0 = uStar_neutral;

    PetscReal uStarDiff = 1e30;
    do {
        // Initial L guess
        PetscReal L0 = -std::pow(uStar0, 3.0) / (kappa * (g / thetaRef) * qwall);
        // Clamp magnitude to avoid zero/NaN (sign follows qwall)
        L0 = (L0 > 0) ? std::max(1e-10, L0) : std::min(-1e-10, L0);

        PetscReal zeta0 = wallDist / L0;
        PetscReal psiM0;

        if (L0 > 0) {  // Stable (qwall < 0)
            psiM0 = -gammaM * zeta0;
            phiM = 1.0 + gammaM * zeta0;
        } else {  // Unstable (qwall > 0)
            // Cap for numerical stability
            if (zeta0 < -1.0 / 15.0) zeta0 = -1.0 / 15.0 + 1e-8;
            PetscReal xi = std::pow(1.0 - 15.0 * zeta0, 0.25);
            psiM0 = 2.0 * std::log((1.0 + xi) / 2.0) + std::log((1.0 + xi * xi) / 2.0) - 2.0 * std::atan(xi) + M_PI / 2.0;
            phiM = std::pow(1.0 - 15.0 * zeta0, -0.25);
        }

        // Update with relaxation (fixed qwall, iterate u*)
        PetscReal uStarNew = std::max(0.0, (kappa * UParallelMeanMag) / (logLaw - psiM0));
        uStar0 = omega * uStarNew + (1.0 - omega) * uStar0;

        uStarDiff = fabs(uStarNew - uStar0);

        iter++;
    } while (uStarDiff > tol && iter < iterMax);

    if (iter == iterMax) {
        char error[512];
        sprintf(error, "Max iterations reached for cell (%ld %ld %ld)\n", k, j, i);
        printf("%s", error);
        uStar = uStar_neutral; L = 1e30; phiM = 1.0;
        return;
    }

    uStar = uStar0;
    L = -std::pow(uStar, 3.0) / (kappa * (g / thetaRef) * qwall);  // Final L
    L = (L > 0) ? std::max(1e-10, L) : std::min(-1e-10, L);

    // Stable low-u* check (only if L>0 and no solution)
    if (L > 0 && uStar < 1e-2) {
        // Bisection for min U_model (df/du*=0, f = U - U_model)
        PetscReal uLow = 1e-6, uHigh = uStar_neutral;
        PetscReal uMid, dfMid, dfLow, dfHigh;
        PetscInt bisIter = 0;
        while (bisIter < 50) {
            uMid = 0.5 * (uLow + uHigh);

            // Use fixed qwall for L_mid
            PetscReal L_mid = std::max(1e-10, -std::pow(uMid, 3.0) / (kappa * (g / thetaRef) * qwall));
            PetscReal zeta_mid = wallDist / L_mid;
            PetscReal psiM_mid = -gammaM * zeta_mid;
            PetscReal U_model_mid = (uMid / kappa) * (logLaw - psiM_mid);
            dfMid = (U_model_mid - UParallelMeanMag) / uMid;  // Approx slope proxy; seek ~0 for extremum

            if (fabs(uHigh - uLow) < 1e-6 || fabs(dfMid) < tol) break;

            // Compute dfLow, dfHigh with fixed qwall
            PetscReal L_low = std::max(1e-10, -std::pow(uLow, 3.0) / (kappa * (g / thetaRef) * qwall));
            PetscReal zeta_low = wallDist / L_low;
            PetscReal psiM_low = -gammaM * zeta_low;
            PetscReal U_model_low = (uLow / kappa) * (logLaw - psiM_low);
            dfLow = (U_model_low - UParallelMeanMag) / uLow;

            PetscReal L_high = std::max(1e-10, -std::pow(uHigh, 3.0) / (kappa * (g / thetaRef) * qwall));
            PetscReal zeta_high = wallDist / L_high;
            PetscReal psiM_high = -gammaM * zeta_high;
            PetscReal U_model_high = (uHigh / kappa) * (logLaw - psiM_high);
            dfHigh = (U_model_high - UParallelMeanMag) / uHigh;

            if ((dfLow > 0 && dfMid > 0) || (dfLow < 0 && dfMid < 0)) uLow = uMid;
            else uHigh = uMid;
            bisIter++;
        }
        
        uStar = uMid;
        PetscReal L_final = std::max(1e-10, -std::pow(uStar, 3.0) / (kappa * (g / thetaRef) * qwall));
        PetscReal zeta_final = wallDist / L_final;
        L = L_final;
        phiM = 1.0 + gammaM * zeta_final;
    }
}

//***************************************************************************************************************//

void qWallShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &z0,
    PetscReal &gammaM, PetscReal &gammaH, PetscReal &alphaH,
    PetscReal &thetaRef, PetscReal &deltaTheta, PetscReal &kappa,
    PetscReal &qWall, PetscReal &uStar, PetscReal &phiM, PetscReal &phiH, PetscReal &L,
    PetscInt k, PetscInt j, PetscInt i
)
{
    if (wallDist <= z0) 
    { 
        char error[512];
        sprintf(error, "wall distance is smaller than the roughness length for cell (%ld %ld %ld)\n", k, j, i);
        fatalErrorInFunction("qWallShumann",  error); 
    }

    PetscReal uStar_neutral = std::max(0.0, (kappa * UParallelMeanMag) / std::log(wallDist / z0));
    PetscReal qWall_neutral = (-deltaTheta * uStar_neutral * kappa) / (alphaH * std::log(wallDist / z0));

    if (fabs(qWall_neutral) < 1e-6) 
    {
        qWall_neutral = signNonZero(qWall_neutral) * 1e-6;  
    }

    PetscReal g = 9.81;
    PetscInt  iter = 0, iterMax = 500;
    PetscReal tol = 1e-5;
    PetscReal omega = 0.5;  

    //neutral case
    if (deltaTheta == 0.0) 
    {
        uStar = uStar_neutral;
        qWall = qWall_neutral;
        L = 1e30;
        phiM = 1.0;
        phiH = alphaH;
        return;
    }

    // Start with neutral for iteration
    PetscReal uStar0 = uStar_neutral;
    PetscReal qWall0 = qWall_neutral;

    PetscReal uStarDiff = 1e30, qWallDiff = 1e30;
    do {
        // Initial L guess
        PetscReal L0 = -std::pow(uStar0, 3.0) / (kappa * (g / thetaRef) * qWall0);
        // Clamp magnitude to avoid zero/NaN 
        L0 = (L0 > 0) ? std::max(1e-10, L0) : std::min(-1e-10, L0);

        PetscReal zeta0 = wallDist / L0;
        PetscReal psiM0, psiH0;

        if (L0 > 0) 
        {   // Stable: 
            psiM0 = -gammaM * zeta0;
            psiH0 = -gammaH * zeta0;
            phiM = 1.0 + gammaM * zeta0;
            phiH = alphaH + gammaH * zeta0;  
        } 
        else 
        {   // Unstable: 
            if (zeta0 < -1.0 / 15.0) zeta0 = -1.0 / 15.0 + 1e-8;
            PetscReal xi = std::pow(1.0 - 15.0 * zeta0, 0.25);
            psiM0 = 2.0 * std::log((1.0 + xi) / 2.0) + std::log((1.0 + xi * xi) / 2.0) - 2.0 * std::atan(xi) + M_PI / 2.0;
            psiH0 = 2.0 * std::log((1.0 + xi * xi) / 2.0);
            phiM = std::pow(1.0 - 15.0 * zeta0, -0.25);
            phiH = alphaH * std::pow(1.0 - 15.0 * zeta0, -0.5);
        }

        // Update with relaxation
        PetscReal uStarNew = std::max(0.0, (kappa * UParallelMeanMag) / (std::log(wallDist / z0) - psiM0));
        uStar0 = omega * uStarNew + (1.0 - omega) * uStar0;

        PetscReal qWallNew = (-deltaTheta * uStar0 * kappa) / (alphaH * std::log(wallDist / z0) - psiH0);
        qWall0 = omega * qWallNew + (1.0 - omega) * qWall0;

        uStarDiff = fabs(uStarNew - uStar0);
        qWallDiff = fabs(qWallNew - qWall0);

        iter++;
    } while ((uStarDiff > tol || qWallDiff > tol) && iter < iterMax);

    if (iter == iterMax) {
        char error[512];
        sprintf(error, "Max iterations reached for cell (%ld %ld %ld)\n", k, j, i);
        printf("%s", error);
        uStar = uStar_neutral; qWall = qWall_neutral; L = 1e30; phiM = 1.0; phiH = alphaH;
        return;
    }

    uStar = uStar0; qWall = qWall0; L = -std::pow(uStar, 3.0) / (kappa * (g / thetaRef) * qWall);  // Final L
    L = (L > 0) ? std::max(1e-10, L) : std::min(-1e-10, L);

    // Stable low-u*, u* is close to zero it means that the zero does not exist, but a maximum exists
    if (L > 0 && uStar < 1e-2) 
    {
        // Neutral u* as upper bound for bisection
        PetscReal uLow = 1e-6, uHigh = uStar_neutral;
        PetscReal uMid, dfMid, dfLow, dfHigh;
        PetscInt  bisIter = 0;
        
        while (bisIter < 50) 
        {
            uMid = 0.5 * (uLow + uHigh);

            // Use fixed neutral qWall for L_mid to decouple
            PetscReal L_mid = std::max(1e-10, -std::pow(uMid, 3.0) / (kappa * (g / thetaRef) * qWall_neutral));
            PetscReal zeta_mid = wallDist / L_mid;
            PetscReal psiM_mid = -gammaM * zeta_mid;
            PetscReal U_model_mid = (uMid / kappa) * (std::log(wallDist / z0) - psiM_mid);
            dfMid = (U_model_mid - UParallelMeanMag) / uMid;  // Approx slope 

            if (fabs(uHigh - uLow) < 1e-6 || fabs(dfMid) < tol) break;

            //sign-based bisection on df 
            PetscReal L_low = std::max(1e-10, -std::pow(uLow, 3.0) / (kappa * (g / thetaRef) * qWall_neutral));
            PetscReal zeta_low = wallDist / L_low;
            PetscReal psiM_low = -gammaM * zeta_low;
            PetscReal U_model_low = (uLow / kappa) * (std::log(wallDist / z0) - psiM_low);
            dfLow = (U_model_low - UParallelMeanMag) / uLow;

            PetscReal L_high = std::max(1e-10, -std::pow(uHigh, 3.0) / (kappa * (g / thetaRef) * qWall_neutral));
            PetscReal zeta_high = wallDist / L_high;
            PetscReal psiM_high = -gammaM * zeta_high;
            PetscReal U_model_high = (uHigh / kappa) * (std::log(wallDist / z0) - psiM_high);
            dfHigh = (U_model_high - UParallelMeanMag) / uHigh;

            if ((dfLow > 0 && dfMid > 0) || (dfLow < 0 && dfMid < 0)) uLow = uMid;
            else uHigh = uMid;
            bisIter++;
        }
        
        uStar = uMid;
        PetscReal L_final = std::max(1e-10, -std::pow(uStar, 3.0) / (kappa * (g / thetaRef) * qWall_neutral));
        PetscReal zeta_final = wallDist / L_final;
        PetscReal psiH_final = -gammaH * zeta_final;
        qWall = (-deltaTheta * uStar * kappa) / (alphaH * std::log(wallDist / z0) - psiH_final);
        L = std::max(1e-10, -std::pow(uStar, 3.0) / (kappa * (g / thetaRef) * qWall));
        phiM = 1.0 + gammaM * zeta_final;
        phiH = alphaH + gammaH * zeta_final;
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

    ut_magb = (*ustar/kappa) * log(sb/roughness);

    if (ut_magb < 1.e-10)
    {
        ut_magb = 0.0;
    }
    
    ut_b  = nScale(ut_magb, et);
    
    (*Ub) = nSum(ut_b, nScale( (sb/sc), unc));
    
    mSum(*Ub, Ua);

    return;
}

//***************************************************************************************************************//
// velocity set to match the background fluid velocity.
void wallShearVelocityBC(PetscReal nu, PetscReal sc, PetscReal sb, PetscReal roughness,
    PetscReal kappa, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    unc  = nScale(nDot(u_c, nf), nf);
    Cmpnts    utc  = nSub(u_c, unc);
    Cmpnts    ut_b;

    Cmpnts    et  = nUnit(utc);

    PetscReal ut_mag = nMag(utc);

    *ustar = ut_mag * kappa / log(sc/roughness);

    (*Ub) = nSum(utc, nScale( (sb/sc), unc));

	(*Ub).x = (*Ub).x + Ua.x;
	(*Ub).y = (*Ub).y + Ua.y;
	(*Ub).z = (*Ub).z + Ua.z;

    return;
}

//***************************************************************************************************************//

void wallShearVelocityBCQuadratic(PetscReal nu,  PetscReal sd, PetscReal sc, PetscReal sb, PetscReal roughness,
    PetscReal kappa, Cmpnts Ua, Cmpnts Ud, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    unc  = nScale(nDot(u_c, nf), nf);
    Cmpnts    utc  = nSub(u_c, unc);

    Cmpnts    u_d = nSub(Ud, Ua);
    Cmpnts    und  = nScale(nDot(u_d, nf), nf);
    Cmpnts    utd  = nSub(u_d, und);    
    Cmpnts    ut_b;

    Cmpnts    et  = nUnit(utd);

    PetscReal ut_mag = nMag(utc);
    PetscReal coeff  = (sc - sb)/(sd - sc);

    //allign utc along utd direction
    utc = nScale(ut_mag, et);

    ut_b.x = utc.x - coeff * (utd.x - utc.x);
    ut_b.y = utc.y - coeff * (utd.y - utc.y);
    ut_b.z = utc.z - coeff * (utd.z - utc.z);

    Cmpnts un_b = nScale((-sb / sc), unc);

    // Enforce no-penetration
    mSub(ut_b, un_b);

    (*Ub).x = ut_b.x + Ua.x;
    (*Ub).y = ut_b.y + Ua.y;
    (*Ub).z = ut_b.z + Ua.z;

    return;
}

//***************************************************************************************************************//

void wallShearGhostVelocityBC
(
    PetscReal sd, PetscReal sc, PetscReal sb, PetscReal roughness,
    PetscReal kappa, Cmpnts Ua, Cmpnts Ud, Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf
)
{
    Cmpnts u_d = nSub(Ud, Ua);
    Cmpnts und = nScale(nDot(u_d, nf), nf);
    Cmpnts utd = nSub(u_d, und);
    PetscReal ut_mag = nMag(utd);
    Cmpnts et = nUnit(utd);

    *ustar = ut_mag * kappa / std::log(sd / roughness);

    PetscReal ut_magc = (*ustar / kappa) * std::log(sc / roughness);
    if (ut_magc < 1.e-10) ut_magc = 0.0;

    Cmpnts ut_c = nScale(ut_magc, et);
    Cmpnts un_c = nScale((sc / sd), und);

    PetscReal coeff = (sc - sb) / (sd - sc);

    Cmpnts ut_b;
    ut_b.x = ut_c.x - coeff * (utd.x - ut_c.x);
    ut_b.y = ut_c.y - coeff * (utd.y - ut_c.y);
    ut_b.z = ut_c.z - coeff * (utd.z - ut_c.z);

    Cmpnts un_b = nScale((-sb / sd), und);

    // Enforce no-penetration
    mSub(ut_b, un_b);

    (*Ub).x = ut_b.x + Ua.x;
    (*Ub).y = ut_b.y + Ua.y;
    (*Ub).z = ut_b.z + Ua.z;

    return;
}

//***************************************************************************************************************//

void ghostTempVelocityBCShumann
(
    Shumann *wm, PetscReal sd, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Ud, Cmpnts *Ub, 
    PetscReal Td, PetscReal surfaceTemp, Cmpnts nf, PetscInt k, PetscInt j, PetscInt i
)
{
    Cmpnts u_d = nSub(Ud, Ua);
    Cmpnts und = nScale(nDot(u_d, nf), nf);
    Cmpnts utd = nSub(u_d, und);
    PetscReal ut_mag = nMag(utd);
    Cmpnts et = nUnit(utd);

    PetscReal surfaceL;

    PetscReal deltaTheta = Td - surfaceTemp;
    PetscReal phiM, phiH;
    PetscReal qWall, ustar;

    qWallShumann
    (
        ut_mag, sd, wm->roughness,
        wm->gammaM, wm->gammaH, wm->alphaH,
        wm->thetaRef, deltaTheta, wm->kappa,
        qWall, ustar, phiM, phiH, surfaceL, k, j, i
    );

    PetscReal psiM;
    PetscReal zeta = sc/surfaceL;

    if(surfaceL > 0)
    {
        psiM = -wm->gammaM * zeta;
    }
    else
    {
        if (zeta < -1.0 / 15.0) zeta = -1.0 / 15.0 + 1e-8;
        PetscReal xi = std::pow(1.0 - 15.0 * zeta, 0.25);
        psiM = 2.0 * std::log((1.0 + xi) / 2.0) + std::log((1.0 + xi * xi) / 2.0) - 2.0 * std::atan(xi) + M_PI / 2.0;
    }

    //set velocity at pt C and interpolate at B
    PetscReal ut_magc = (ustar / wm->kappa) * (std::log(sc / wm->roughness) - psiM);
    if (ut_magc < 1.e-10) ut_magc = 0.0;

    Cmpnts ut_c = nScale(ut_magc, et);

    PetscReal coeff = (sc - sb) / (sd - sc);

    Cmpnts ut_b;
    ut_b.x = ut_c.x - coeff * (utd.x - ut_c.x);
    ut_b.y = ut_c.y - coeff * (utd.y - ut_c.y);
    ut_b.z = ut_c.z - coeff * (utd.z - ut_c.z);

    Cmpnts un_b = nScale((-sb / sd), und);

    // Enforce no-penetration
    mSub(ut_b, un_b);

    (*Ub).x = ut_b.x + Ua.x;
    (*Ub).y = ut_b.y + Ua.y;
    (*Ub).z = ut_b.z + Ua.z;

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

    (*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;
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

    (*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;

    return;
}

void slipBC(PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc, Cmpnts *Ub, Cmpnts nf)
{

    Cmpnts    u_c = nSub(Uc, Ua);
    Cmpnts    un  = nScale(nDot(u_c, nf), nf);
    Cmpnts    ut  = nSub(u_c, un);

    *Ub = nSum(ut, nScale(sb/sc, un));

    (*Ub).x += Ua.x;
	(*Ub).y += Ua.y;
	(*Ub).z += Ua.z;

	return;
}
