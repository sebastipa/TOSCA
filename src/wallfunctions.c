#include "include/base.h"
#include "include/domain.h"
#include "include/wallfunctions.h"

PetscInt     pre_integrate_flag = 0;
PetscReal    *integration_buffer;
PetscReal    *integration_buffer_rough;

void uStarShumann
(
    PetscReal &UParallelMeanMag, PetscReal &wallDist, PetscReal &z0,
    PetscReal &gammaM, PetscReal &kappa, PetscReal &qwall, PetscReal &thetaRef,
    PetscReal &uStar, PetscReal &phiM, PetscReal &L
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
            uStar1 = uStar + duStar;

            printf("uStar is zero, looking for maximum...");

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
                f0 = UParallelMeanMag - (uStar0/kappa)*(std::log(wallDist/z0) - psiM0);
                f1 = UParallelMeanMag - (uStar1/kappa)*(std::log(wallDist/z0) - psiM1);

                uStar = uStar0;
                L     = L0;
                phiM  = 1.0 + (gammaM * zeta0);

            } while (f1 > f0);

            printf("uStar = %f\n", uStar);
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

inline PetscReal near_wall_eddy_viscosity(PetscReal yplus)
{
    PetscReal   kappa = 0.41;
    return kappa * yplus * pow(1. - exp(-yplus / 19.), 2.0);
}

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
    std
    ::vector<PetscReal> E(N+1);
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

inline PetscReal f_Cabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    return (utau * utau * integrate_F(nu, utau, y) - u);
}

inline PetscReal df_Cabot(PetscReal nu, PetscReal u, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    PetscReal eps = 1.e-7;

    return
    (
        f_Cabot(nu, u, y, utau + eps, dpdn) -
        f_Cabot(nu, u, y, utau - eps, dpdn)
    ) / (2.0 * eps);
}

PetscReal u_Cabot(PetscReal nu, PetscReal y, PetscReal utau, PetscReal dpdn)
{
    return utau * utau * integrate_F(nu, utau, y); // + dpdn * integrate_Fy( nu, utau, 0, y );
}

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

void wallFunctionCabot(PetscReal nu, PetscReal sc, PetscReal sb, Cmpnts Ua, Cmpnts Uc,
        Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{

    double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
    double un = u_c * nf.x + v_c * nf.y + w_c * nf.z;
    double ut = u_c - un * nf.x, vt = v_c - un * nf.y, wt = w_c - un * nf.z;
    double ut_mag = sqrt(ut * ut + vt * vt + wt * wt);
    double ut_mag_modeled;

    *ustar = uTauCabot(nu, ut_mag, sc, 0.01, 0);
    ut_mag_modeled = u_Cabot(nu, sb, *ustar, 0);

    double Ua_abs = sqrt(Ua.x * Ua.x + Ua.y * Ua.y + Ua.z * Ua.z);

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

// Wall function based on the power law model
void wallFunctionPowerlaw(double nu, double sc, double sb, Cmpnts Ua,
        Cmpnts Uc, Cmpnts *Ub, PetscReal *ustar, Cmpnts nf)
{
    double u_c = Uc.x - Ua.x, v_c = Uc.y - Ua.y, w_c = Uc.z - Ua.z;
    double un = u_c * nf.x + v_c * nf.y + w_c * nf.z;
    double ut = u_c - un * nf.x, vt = v_c - un * nf.y, wt = w_c - un * nf.z;
    double ut_mag = sqrt(ut * ut + vt * vt + wt * wt);

    double A = 8.3;
    double B = 1.0 / 7.0;
    *ustar = pow(ut_mag * pow(nu, B) / (A * pow(sc, B)), 1.0 / (1.0 + B));  ///

    double ybp = sb * (*ustar) / nu;
    double ycp = sc * (*ustar) / nu;
    double ut_mag_modeled;
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
