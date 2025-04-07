//! \file  inline.h
//! \brief Inline functions

#ifndef _INLINE_
#define _INLINE_

#include "io.h"
#include "ibm.h"

// DYNAMIC TIME STEPPING
// ============================================================================================================= //

inline PetscReal currentDistanceToWriteTime(clock_ *clock, PetscReal timeStart, PetscReal timeInterval)
{
    PetscReal currentAcquistionPeriodFraction = ((clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval + 1e-10));
    return((1.0-currentAcquistionPeriodFraction)*timeInterval);
}

inline void timeStepSet(clock_ *clock, PetscReal timeStart, PetscReal timeInterval, PetscReal dByU, PetscInt &flag, PetscReal &cfl)
{
    PetscReal initialTimeStep  = clock->dt;

    // compute current fraction of the acquisition period
    PetscReal currentAcquistionPeriodFraction = (clock->time - timeStart) / timeInterval - std::floor((clock->time - timeStart) / timeInterval + 1e-10);

    // compute current distance to write time
    PetscReal currentDistanceToWriteTime = (1.0-currentAcquistionPeriodFraction)*timeInterval;

    // compute what will be the two CFL if adjusting with current or next distance
    PetscReal currentCFL = currentDistanceToWriteTime / dByU;

    if(currentCFL < clock->cfl)
    {
        if(clock->time + initialTimeStep >= timeStart)
        {
            if(currentDistanceToWriteTime < clock->dt)
            {
                clock->dt = currentDistanceToWriteTime;
                flag = 1;
                cfl  = currentCFL;
            }

        }
    }

    return;
}

inline PetscInt mustWrite(PetscReal time, PetscReal startTime, PetscReal timeInterval)
{
    return
    (
        (
            (time - startTime) / timeInterval -
            std::floor
            (
                (time - startTime) / timeInterval + 1e-10
            ) < 1e-10
        ) &&
        (time >= startTime)
    );
}

// MATH FUNCTIONS
// ============================================================================================================= //

//find the matrix product of two matrices 

inline void matMatProduct(PetscReal **A, PetscReal **B, PetscReal **C, PetscInt numRowA, PetscInt numColA, PetscInt numRowB, PetscInt numColB)
{
    if( numColA != numRowB )
    {
        char error[512];
        sprintf(error, "number of columns of matrix 1 not equal to number of rows of matrix 2\n");
        fatalErrorInFunction("matMatProduct",  error);
    }

        for (PetscInt i = 0; i < numRowA; i++) 
        {
            for (PetscInt j = 0; j < numColB; j++) 
            {
                PetscReal dotProduct = 0.0;

                for (PetscInt k = 0; k < numColA; k++) 
                {
                    dotProduct += A[i][k] * B[k][j];
                }
                C[i][j] = dotProduct;
            }
        }

    return;
}

inline void matVecProduct(PetscReal **A, PetscReal *b, PetscReal *c, PetscInt numRowA, PetscInt numColA, PetscInt numRowB)
{
    if( numColA != numRowB )
    {
        char error[512];
        sprintf(error, "number of columns of matrix 1 not equal to number of rows of vector 2\n");
        fatalErrorInFunction("matVecProduct",  error);
    }

    for (PetscInt i = 0; i < numRowA; i++) 
    {
        PetscReal dotProduct = 0.0;

        for (PetscInt k = 0; k < numColA; k++) 
        {
            dotProduct += A[i][k] * b[k];
        }
        c[i] = dotProduct;
    }

    return;
}

//=============================================================================================================

inline void gaussianSmooth1D(PetscReal *input, PetscReal *output, PetscInt n, PetscInt window)
{
    PetscReal weights[window];
    PetscReal sigma = window/3.0;
    PetscReal sum = 0;
    
    // Create Gaussian kernel
    for(PetscInt i=0; i<window; i++) 
    {
        weights[i] = exp(-0.5*pow((i-window/2)/sigma, 2));
        sum += weights[i];
    }

    //normalize
    for(PetscInt i=0; i<window; i++) weights[i] /= sum;
    
    // Apply convolution
    for(PetscInt i=0; i<n; i++) 
    {
        output[i] = 0;
        for(PetscInt j=-window/2; j<=window/2; j++) 
        {
            PetscInt idx = i + j;
            if(idx < 0) idx = 0;
            if(idx >= n) idx = n-1;
            output[i] += input[idx] * weights[j + window/2];
        }
    }   
    return;
}

//=============================================================================================================
//check for duplicates in an integer array
inline PetscInt isPresent(PetscInt arr[], PetscInt n, PetscInt elem)
{
    for (PetscInt i = 0; i < n; i++)
    {
        if (arr[i] == elem)
        {
            return 1;
        }
    }

    return 0;
}

// ============================================================================================================= //

inline PetscReal sign(PetscReal a)
{
    if (a > 0)      return  1.0;
    else if (a < 0) return -1.0;
    else            return  0.0;

}

// ============================================================================================================= //

inline PetscReal signNonZero(PetscReal a)
{
    if (a > 0)      return  1.0;
    else if (a < 0) return -1.0;
    else            return  1.0;

}

// ============================================================================================================= //

inline PetscReal gcd(PetscReal a, PetscReal b)
{
    if (a < b) return gcd(b, a);

    // base case
    if (fabs(b) < 1e-8) return a;

    else return (gcd(b, a - floor(a / b) * b));

}

// ============================================================================================================= //

inline PetscReal gcdN(PetscReal *v, PetscInt n)
{
    PetscReal result = v[0];

    for (PetscInt i = 1; i < n; i++)
    {
        result = gcd(v[i], result);
    }

    return result;
}

// ============================================================================================================= //

inline complex gamma(complex x)
{
    complex f, z = x;

    PetscReal c[11] =
    {
        1.000000000000000174663,
        5716.400188274341379136,
        -14815.30426768413909044,
        14291.49277657478554025,
        -6348.160217641458813289,
        1301.608286058321874105,
        -108.1767053514369634679,
        2.605696505611755827729,
       -0.7423452510201416151527e-2,
        0.5384136432509564062961e-7,
       -0.4023533141268236372067e-8
    };

    if(x.real() < 0.0) z = -1.0*z;

    complex   ss, g = complex(9.0,0.0), t = z + g, s = complex(0.0, 0.0);

    for(PetscInt k=10; k>0; k--)
    {
        s = s + (complex)c[k] / t;
        t = t - 1.0;
    }

    s  = s + (complex)c[0];
    ss = (z + g - 0.5);
    s  = std::log(s * std::sqrt(2.0*M_PI)) + (z - 0.5) * std::log(ss) - ss;
    f  = std::exp(s);

    // recursive formula
    if(x.real() < 0.0) f = -1.0 * M_PI / (x * f * std::sin(M_PI*x));

    // negative integers
    if(std::round(x.real()) == x.real() && x.imag() == 0.0 && x.real() <= 0.0) f = 1e20;

    // exact result for integers arguments
    if(std::round(x.real()) == x.real() && x.imag() == 0.0 && x.real() > 0.0)
    {
        complex pp = complex(1.0, 0.0);

        for(PetscInt m = 1; m<(PetscInt)x.real(); m++)
        {
            pp = pp*complex((PetscReal)m,0.0);
        }
        f = pp;
    }

    return(f);
}

// ============================================================================================================= //

inline complex digamma(complex x)
{
    complex f, z = x;

    PetscReal c[15] =
    {
        0.99999999999999709182,
       57.156235665862923517,
       -59.597960355475491248,
       14.136097974741747174,
       -0.49191381609762019978,
       .33994649984811888699e-4,
       .46523628927048575665e-4,
       -.98374475304879564677e-4,
       .15808870322491248884e-3,
       -.21026444172410488319e-3,
       .21743961811521264320e-3,
       -.16431810653676389022e-3,
       .84418223983852743293e-4,
       -.26190838401581408670e-4,
       .36899182659531622704e-5
    };

    if(x.real() < 0.5) z = 1.0 - z;

    complex dz, dd, gg, g = complex(607.0/128.0, 0.0), n = complex(0.0,0.0), d = n;

    for(PetscInt k=14; k>0; k--)
    {
        dz = 1.0 / (z + (PetscReal)(k+1) - 2.0);
        dd = c[k] * dz;
        d  = d + dd;
        n  = n - dd * dz;
    }

    d  = d + c[0];
    gg = z + g - 0.5;

    f = std::log(gg) + (n/d - g/gg);

    if(x.real() < 0.5) f = f - M_PI * (1.0 / tan(M_PI * x));

    if(std::round(x.real()) == x.real() && x.imag() == 0.0 && x.real() <= 0.0) f = 1e20;

    return(f);
}

// ============================================================================================================= //

inline complex hypergeom(complex a, complex b, complex c, PetscReal z)
{
	if(std::round(c.real()) == c.real() && c.imag() == 0 && c.real() <= 0)
	{
		return(complex(1e20,0.0));
	}

    PetscReal epsilon = 1e-15;
    complex   term    = complex(1.0,0.0),
              result  = complex(1.0,0.0);
    PetscInt  termination = 0, k;

    for(k=0; k<30000; k++)
    {
        term = term * ((a + (PetscReal)k) * (b + (PetscReal)k) * z / ( (c + (PetscReal)k) * ((PetscReal)k + 1.0)));

        if(std::abs(term) < abs(result*epsilon))
        {
            termination ++;
            if(termination >= 3)
            {
                return(result);
            }
        }
        else
        {
            termination = 0;
        }

        if(term == 0.0) break;

        result = result + term;

        if(std::isnan(result.real()) || std::isnan(result.imag()))
        {
            char warning[256];
            sprintf(warning, "NaN value detected in hypergeometric function at z = %.5f\n", z);
            warningInFunction("hypergeom",  warning);

            return(complex(std::nan(""),std::nan("")));
        }
    }

    return(result);
}

// ============================================================================================================= //

inline PetscReal polyLog2(PetscReal x)
{
    const PetscReal PI = 3.1415926535897932;
    const PetscReal P[] =
    {
      0.9999999999999999502e+0,
     -2.6883926818565423430e+0,
      2.6477222699473109692e+0,
     -1.1538559607887416355e+0,
      2.0886077795020607837e-1,
     -1.0859777134152463084e-2
    };
    const PetscReal Q[] =
    {
      1.0000000000000000000e+0,
     -2.9383926818565635485e+0,
      3.2712093293018635389e+0,
     -1.7076702173954289421e+0,
      4.1596017228400603836e-1,
     -3.9801343754084482956e-2,
      8.2743668974466659035e-4
    };

    PetscReal y = 0, r = 0, s = 1;

    // transform to [0, 1/2]
    if (x < -1) {
      const PetscReal l = log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - log(-x));
      s = 1;
    } else if (x == -1) {
      return -PI*PI/12;
    } else if (x < 0) {
      const PetscReal l = log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
    } else if (x == 0) {
      return 0;
    } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
    } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - log(x)*log1p(-x);
      s = -1;
    } else if (x == 1) {
      return PI*PI/6;
    } else if (x < 2) {
      const PetscReal l = log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(log(y) + 0.5*l);
      s = 1;
    } else {
      const PetscReal l = log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
    }

    const PetscReal y2 = y*y;
    const PetscReal y4 = y2*y2;
    const PetscReal p = P[0] + y * P[1] + y2 * (P[2] + y * P[3]) +
                    y4 * (P[4] + y * P[5]);
    const PetscReal q = Q[0] + y * Q[1] + y2 * (Q[2] + y * Q[3]) +
                    y4 * (Q[4] + y * Q[5] + y2 * Q[6]);

    return r + s*y*p/q;
}

// ============================================================================================================= //

inline PetscReal splineB1(PetscReal a)
{
    if (fabs(a) <= 1)      return 1.0 - fabs(a);
    else                   return 0.0;
}

// ============================================================================================================= //

inline PetscReal splineB2(PetscReal a)
{
    if      (fabs(a) <= 0.5)      return 3.0/4.0 - a*a;
    else if (fabs(a) <= 1.5 &&
             fabs(a) > 0.5)       return 1.0/8.0 * (2.0*fabs(a) - 3.0)*(2.0*fabs(a) - 3.0);
    else                          return 0.0;

}

// VECTOR ALGEBRA
// ============================================================================================================= //

inline Cmpnts nMax(Cmpnts v, PetscReal c)
{
    Cmpnts vMax;

    vMax.x = PetscMax(v.x, c);
    vMax.y = PetscMax(v.y, c);
    vMax.z = PetscMax(v.z, c);

    return(vMax);
}

//***************************************************************************************************************//

inline Cmpnts nAbs(Cmpnts v)
{
    Cmpnts vR;

    vR.x = fabs(v.x);
    vR.y = fabs(v.y);
    vR.z = fabs(v.z);

    return(vR);
}

//***************************************************************************************************************//

inline Cmpnts nInv(Cmpnts v)
{
    Cmpnts vR;

    vR.x = 1.0 / v.x;
    vR.y = 1.0 / v.y;
    vR.z = 1.0 / v.z;

    return(vR);
}

//***************************************************************************************************************//

inline Cmpnts nPow(Cmpnts v, PetscReal c)
{
    Cmpnts vR;

    vR.x = pow(v.x, c);
    vR.y = pow(v.y, c);
    vR.z = pow(v.z, c);

    return(vR);
}

//***************************************************************************************************************//

inline Cmpnts nScale(PetscReal c, Cmpnts v)
{
    Cmpnts vScaled;

    vScaled.x = c * v.x;
    vScaled.y = c * v.y;
    vScaled.z = c * v.z;

    return(vScaled);
}

//***************************************************************************************************************//

inline Cmpnts nMultElWise(Cmpnts v1, Cmpnts v2)
{
    Cmpnts vR;

    vR.x = v1.x * v2.x;
    vR.y = v1.y * v2.y;
    vR.z = v1.z * v2.z;

    return(vR);
}

//***************************************************************************************************************//

inline Cmpnts nDivElWise(Cmpnts v1, Cmpnts v2)
{
    Cmpnts vR;

    vR.x = v1.x / v2.x;
    vR.y = v1.y / v2.y;
    vR.z = v1.z / v2.z;

    return(vR);
}

//***************************************************************************************************************//

inline Cmpnts nScaleX(PetscReal c, Cmpnts v)
{
    Cmpnts vScaled;

    vScaled.x = c * v.x;
    vScaled.y = v.y;
    vScaled.z = v.z;

    return(vScaled);
}

//***************************************************************************************************************//

inline Cmpnts nScaleY(PetscReal c, Cmpnts v)
{
    Cmpnts vScaled;

    vScaled.x = v.x;
    vScaled.y = c * v.y;
    vScaled.z = v.z;

    return(vScaled);
}

//***************************************************************************************************************//

inline Cmpnts nScaleZ(PetscReal c, Cmpnts v)
{
    Cmpnts vScaled;

    vScaled.x = v.x;
    vScaled.y = v.y;
    vScaled.z = c * v.z;

    return(vScaled);
}

//***************************************************************************************************************//

inline PetscReal nDot(Cmpnts v1, Cmpnts v2)
{
    return
    (
        v1.x*v2.x +
        v1.y*v2.y +
        v1.z*v2.z
    );
}

//***************************************************************************************************************//

inline PetscReal nMag(Cmpnts v)
{
    return
    (
        sqrt
        (
            v.x * v.x +
            v.y * v.y +
            v.z * v.z
        )
    );
}

//***************************************************************************************************************//

inline Cmpnts nCross(Cmpnts v1, Cmpnts v2)
{
    Cmpnts vCross;

    vCross.x = v1.y * v2.z - v1.z * v2.y;
    vCross.y = v1.z * v2.x - v1.x * v2.z;
    vCross.z = v1.x * v2.y - v1.y * v2.x;

    // if(nMag(vCross)<1e-20)
    // {
    //     char error[512];
    //     sprintf(error, "vectors are parallel...was told to exit\n");
    //     fatalErrorInFunction("nCross",  error);
    // }

    return(vCross);
}

//***************************************************************************************************************//

inline Cmpnts nRot(Cmpnts axis, Cmpnts vStart, PetscReal theta)
{
    Cmpnts vRotated = {0.0,0.0,0.0};
    PetscReal costheta,sintheta;

    // normalize rotation axis vector (add a small value for division)
    PetscReal axisMag
    =
    sqrt
    (
        axis.x * axis.x +
        axis.y * axis.y +
        axis.z * axis.z
    );

    if(axisMag < 1.0e-20)
    {
        char error[512];
        sprintf(error, "provided rotation axis too small to working precision\n");
        fatalErrorInFunction("nRot",  error);
    }

    axis.x = axis.x / axisMag;
    axis.y = axis.y / axisMag;
    axis.z = axis.z / axisMag;

    // compute sine and cosine
    costheta = cos(theta);
    sintheta = sin(theta);

    // explicitly compute rotation
    vRotated.x += (costheta + (1 - costheta) * axis.x * axis.x         ) * vStart.x;
    vRotated.x += ((1 - costheta) * axis.x * axis.y - axis.z * sintheta) * vStart.y;
    vRotated.x += ((1 - costheta) * axis.x * axis.z + axis.y * sintheta) * vStart.z;

    vRotated.y += ((1 - costheta) * axis.x * axis.y + axis.z * sintheta) * vStart.x;
    vRotated.y += (costheta + (1 - costheta) * axis.y * axis.y         ) * vStart.y;
    vRotated.y += ((1 - costheta) * axis.y * axis.z - axis.x * sintheta) * vStart.z;

    vRotated.z += ((1 - costheta) * axis.x * axis.z - axis.y * sintheta) * vStart.x;
    vRotated.z += ((1 - costheta) * axis.y * axis.z + axis.x * sintheta) * vStart.y;
    vRotated.z += (costheta + (1 - costheta) * axis.z * axis.z         ) * vStart.z;

    return(vRotated);

}

//***************************************************************************************************************//

inline Cmpnts nTra(Cmpnts point, Cmpnts translation)
{
    Cmpnts vTranslated;

    vTranslated.x = point.x + translation.x;
    vTranslated.y = point.y + translation.y;
    vTranslated.z = point.z + translation.z;

    return(vTranslated);
}

//***************************************************************************************************************//

inline Cmpnts nUnit(Cmpnts v)
{
    PetscReal eps=1.e-12;

    PetscReal vMag
    =
    sqrt
    (
        v.x * v.x +
        v.y * v.y +
        v.z * v.z
    );

    Cmpnts vUnit;

    vUnit.x = v.x / (vMag + eps);
    vUnit.y = v.y / (vMag + eps);
    vUnit.z = v.z / (vMag + eps);

    return(vUnit);
}

//***************************************************************************************************************//

inline Cmpnts nSum(Cmpnts v1, Cmpnts v2)
{
    Cmpnts vSum;

    vSum.x = v1.x + v2.x;
    vSum.y = v1.y + v2.y;
    vSum.z = v1.z + v2.z;

    return(vSum);
}

//***************************************************************************************************************//

inline Cmpnts nSub(Cmpnts v1, Cmpnts v2)
{
    Cmpnts vSub;

    vSub.x = v1.x - v2.x;
    vSub.y = v1.y - v2.y;
    vSub.z = v1.z - v2.z;

    return(vSub);
}

//***************************************************************************************************************//

inline Cmpnts nSet(Cmpnts value)
{
    Cmpnts newVector;

    newVector.x = value.x;
    newVector.y = value.y;
    newVector.z = value.z;

    return(newVector);
}

//***************************************************************************************************************//

inline Cmpnts nSetFromComponents(PetscReal vx, PetscReal vy, PetscReal vz)
{
    Cmpnts newVector;

    newVector.x = vx;
    newVector.y = vy;
    newVector.z = vz;

    return(newVector);
}

//***************************************************************************************************************//

inline Cmpnts nSetZero()
{
    Cmpnts newVector;

    newVector.x = 0.0;
    newVector.y = 0.0;
    newVector.z = 0.0;

    return(newVector);
}

//***************************************************************************************************************//

inline void mUnit(Cmpnts &base)
{
    PetscReal vMag
    =
    sqrt
    (
        base.x * base.x +
        base.y * base.y +
        base.z * base.z
    );

    base.x = base.x / vMag;
    base.y = base.y / vMag;
    base.z = base.z / vMag;

    return;
}

//***************************************************************************************************************//

inline void mSum(Cmpnts &base, Cmpnts add)
{
    base.x = base.x + add.x;
    base.y = base.y + add.y;
    base.z = base.z + add.z;

    return;
}

//***************************************************************************************************************//

inline void mSub(Cmpnts &base, Cmpnts sub)
{
    base.x = base.x - sub.x;
    base.y = base.y - sub.y;
    base.z = base.z - sub.z;

    return;
}

//***************************************************************************************************************//

inline void mScale(PetscReal c, Cmpnts &v)
{
    v.x = c * v.x;
    v.y = c * v.y;
    v.z = c * v.z;

    return;
}

//***************************************************************************************************************//

inline void mScaleX(PetscReal c, Cmpnts &v)
{
    v.x = c * v.x;

    return;
}

//***************************************************************************************************************//

inline void mScaleY(PetscReal c, Cmpnts &v)
{
    v.y = c * v.y;

    return;
}

//***************************************************************************************************************//

inline void mScaleZ(PetscReal c, Cmpnts &v)
{
    v.z = c * v.z;

    return;
}

//***************************************************************************************************************//

inline void mSetScale(PetscReal scale, Cmpnts &a, Cmpnts b)
{
    a.x = scale*b.x;
    a.y = scale*b.y;
    a.z = scale*b.z;
}

//***************************************************************************************************************//

inline void mRot(Cmpnts axis, Cmpnts &v, PetscReal theta)
{
    Cmpnts vRotated = {0.0,0.0,0.0};
    PetscReal costheta,sintheta;

    // normalize rotation axis vector (add a small value for division)
    PetscReal axisMag
    =
    sqrt
    (
        axis.x * axis.x +
        axis.y * axis.y +
        axis.z * axis.z
    );

    if(axisMag < 1.0e-20)
    {
        char error[512];
        sprintf(error, "provided rotation axis too small to working precision\n");
        fatalErrorInFunction("mRot",  error);
    }

    axis.x = axis.x / axisMag;
    axis.y = axis.y / axisMag;
    axis.z = axis.z / axisMag;

    // compute sine and cosine
    costheta = cos(theta);
    sintheta = sin(theta);

    // explicitly compute rotation
    vRotated.x += (costheta + (1 - costheta) * axis.x * axis.x         ) * v.x;
    vRotated.x += ((1 - costheta) * axis.x * axis.y - axis.z * sintheta) * v.y;
    vRotated.x += ((1 - costheta) * axis.x * axis.z + axis.y * sintheta) * v.z;

    vRotated.y += ((1 - costheta) * axis.x * axis.y + axis.z * sintheta) * v.x;
    vRotated.y += (costheta + (1 - costheta) * axis.y * axis.y         ) * v.y;
    vRotated.y += ((1 - costheta) * axis.y * axis.z - axis.x * sintheta) * v.z;

    vRotated.z += ((1 - costheta) * axis.x * axis.z - axis.y * sintheta) * v.x;
    vRotated.z += ((1 - costheta) * axis.y * axis.z + axis.x * sintheta) * v.y;
    vRotated.z += (costheta + (1 - costheta) * axis.z * axis.z         ) * v.z;

    v.x = vRotated.x;
    v.y = vRotated.y;
    v.z = vRotated.z;

    return;
}

//***************************************************************************************************************//

inline void mTra(Cmpnts &point, Cmpnts translation)
{
    point.x = point.x + translation.x;
    point.y = point.y + translation.y;
    point.z = point.z + translation.z;

    return;
}

//***************************************************************************************************************//

inline void mSet(Cmpnts &base, Cmpnts value)
{
    base.x = value.x;
    base.y = value.y;
    base.z = value.z;

    return;
}

//***************************************************************************************************************//

inline void mSetValue(Cmpnts &base, PetscReal value)
{
    base.x = value;
    base.y = value;
    base.z = value;

    return;
}

//! C=aX+bY
inline void AxByC ( PetscReal a, Cmpnts &X, PetscReal b, Cmpnts &Y, Cmpnts *C)
{
    (*C).x = a*X.x + b*Y.x;
    (*C).y = a*X.y + b*Y.y;
    (*C).z = a*X.z + b*Y.z;

    return;
}
//***************************************************************************************************************//

// domain checks
inline bool isInsideBoundingBox ( Cmpnts pt, const boundingBox &simBox)
{
    return (pt.x > simBox.xmin
       &&   pt.x < simBox.xmax
       &&   pt.y > simBox.ymin
       &&   pt.y < simBox.ymax
       &&   pt.z > simBox.zmin
       &&   pt.z < simBox.zmax ) ;

}

// FIELD AVERAGES 
// ============================================================================================================= //

inline std::vector<PetscReal> jPlaneScalarMean(mesh_  *mesh, PetscReal ***s, PetscInt nLevels)
{
    DM            da = mesh->da, fda = mesh->fda;
    DMDALocalInfo info = mesh->info;
    PetscInt      xs = info.xs, xe = info.xs + info.xm;
    PetscInt      ys = info.ys, ye = info.ys + info.ym;
    PetscInt      zs = info.zs, ze = info.zs + info.zm;
    PetscInt      mx = info.mx, my = info.my, mz = info.mz;

    PetscInt      lxs, lxe, lys, lye, lzs, lze;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    PetscInt      i, j, k;
    PetscReal     ***aj;

    std::vector<PetscReal> lMean(nLevels, 0.0);
    std::vector<PetscReal> gMean(nLevels, 0.0);
    std::vector<PetscReal> lVol(nLevels, 0.0);
    std::vector<PetscReal> gVol(nLevels, 0.0);

    DMDAVecGetArray(da,  mesh->lAj, &aj);

    for (k=lzs; k<lze; k++)
    {
        for (j=lys; j<lye; j++)
        {
            for (i=lxs; i<lxe; i++)
            {
                lMean[j-1] += s[k][j][i] / aj[k][j][i];
                lVol[j-1]  += 1.0 / aj[k][j][i];
            }
        }
    }

    MPI_Allreduce(&lMean[0], &gMean[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);
    MPI_Allreduce(&lVol[0], &gVol[0], nLevels, MPIU_REAL, MPIU_SUM, mesh->MESH_COMM);

    for(j=0; j<nLevels; j++)
    {
        gMean[j] = gMean[j]/gVol[j];
    }
    
    DMDAVecRestoreArray(da,  mesh->lAj, &aj);

    return gMean;
}

// MESH OPERATIONS
// ============================================================================================================= //

inline void Calculate_Covariant_metrics(PetscReal g[3][3], PetscReal G[3][3])
{
    /*
        | csi.x  csi.y csi.z |-1    | x.csi   x.eta  x.zet |
        | eta.x  eta.y eta.z |    = | y.csi   y.eta  y.zet |
        | zet.x  zet.y zet.z |      | z.csi   z.eta  z.zet |

    */
    const PetscReal a11=g[0][0], a12=g[0][1], a13=g[0][2];
    const PetscReal a21=g[1][0], a22=g[1][1], a23=g[1][2];
    const PetscReal a31=g[2][0], a32=g[2][1], a33=g[2][2];

    PetscReal det= a11*(a33*a22-a32*a23)-a21*(a33*a12-a32*a13)+a31*(a23*a12-a22*a13);

    G[0][0] = (a33*a22-a32*a23)/det,    G[0][1] = - (a33*a12-a32*a13)/det,   G[0][2] = (a23*a12-a22*a13)/det;
    G[1][0] = -(a33*a21-a31*a23)/det,   G[1][1] = (a33*a11-a31*a13)/det,     G[1][2] = - (a23*a11-a21*a13)/det;
    G[2][0] = (a32*a21-a31*a22)/det,    G[2][1] = - (a32*a11-a31*a12)/det,   G[2][2] = (a22*a11-a21*a12)/det;

    return;
};

inline void calculateNormal(Cmpnts &csi, Cmpnts &eta, Cmpnts &zet, PetscReal ni[3], PetscReal nj[3], PetscReal nk[3])
{
    PetscReal g[3][3];
    PetscReal G[3][3];

    g[0][0]=csi.x, g[0][1]=csi.y, g[0][2]=csi.z;
    g[1][0]=eta.x, g[1][1]=eta.y, g[1][2]=eta.z;
    g[2][0]=zet.x, g[2][1]=zet.y, g[2][2]=zet.z;

    Calculate_Covariant_metrics(g, G);
    PetscReal xcsi=G[0][0], ycsi=G[1][0], zcsi=G[2][0];
    PetscReal xeta=G[0][1], yeta=G[1][1], zeta=G[2][1];
    PetscReal xzet=G[0][2], yzet=G[1][2], zzet=G[2][2];

    PetscReal nx_i = xcsi, ny_i = ycsi, nz_i = zcsi;
    PetscReal nx_j = xeta, ny_j = yeta, nz_j = zeta;
    PetscReal nx_k = xzet, ny_k = yzet, nz_k = zzet;

    PetscReal sum_i=sqrt(nx_i*nx_i+ny_i*ny_i+nz_i*nz_i);
    PetscReal sum_j=sqrt(nx_j*nx_j+ny_j*ny_j+nz_j*nz_j);
    PetscReal sum_k=sqrt(nx_k*nx_k+ny_k*ny_k+nz_k*nz_k);

    nx_i /= sum_i, ny_i /= sum_i, nz_i /= sum_i;
    nx_j /= sum_j, ny_j /= sum_j, nz_j /= sum_j;
    nx_k /= sum_k, ny_k /= sum_k, nz_k /= sum_k;

    ni[0] = nx_i, ni[1] = ny_i, ni[2] = nz_i;
    nj[0] = nx_j, nj[1] = ny_j, nj[2] = nz_j;
    nk[0] = nx_k, nk[1] = ny_k, nk[2] = nz_k;

    return;
}

// CELL DEFINITION TEST
// ============================================================================================================= //

//***************************************************************************************************************//

//! check if ibm cell within the neighbouring box of i,j,k
inline bool isBoxIBMCell (int k, int j, int i, PetscReal ***nvert)
{
    return (nvert[k  ][j  ][i  ] +
            nvert[k+1][j  ][i  ] +
            nvert[k-1][j  ][i  ] +
            nvert[k  ][j+1][i  ] +
            nvert[k  ][j-1][i  ] +
            nvert[k  ][j  ][i+1] +
            nvert[k  ][j  ][i-1] > 0.1);
}

//***************************************************************************************************************//

//! \brief check if purely a fluid cell
inline bool isFluidCell (PetscInt k, PetscInt j, PetscInt i, PetscReal ***nvert)
{
    return (nvert[k][j][i] < 0.1);
}

//***************************************************************************************************************//

//! \brief check if IBM cell - this includes IBM fluid cells (cell next to IBM body) and IBM solid cells
inline bool isIBMCell (PetscInt k, PetscInt j, PetscInt i, PetscReal ***nvert)
{
    return (nvert[k][j][i] > 0.1);
}

//***************************************************************************************************************//

//! \brief check if IBM Fluid cell (cell next to IBM body)
inline bool isIBMFluidCell (PetscInt k, PetscInt j, PetscInt i, PetscReal ***nvert)
{
    return ((nvert[k][j][i] > 0.1) && (nvert[k][j][i] < 1.5));
}

//***************************************************************************************************************//

//! \brief check if IBM solid cell
inline bool isIBMSolidCell (PetscInt k, PetscInt j, PetscInt i, PetscReal ***nvert)
{
    return (nvert[k][j][i] > 1.5);
}

//***************************************************************************************************************//

//! \brief check is a fluid k face - shared by 2 fluid cells
inline bool isFluidKFace (PetscInt kL, PetscInt j, PetscInt i, PetscInt kR, PetscReal ***nvert)
{
    return ((nvert[kL][j][i]+nvert[kR][j][i]) < 0.1);
}

//***************************************************************************************************************//

//! \brief check is a fluid j face - shared by 2 fluid cells
inline bool isFluidJFace (PetscInt k, PetscInt jL, PetscInt i, PetscInt jR, PetscReal ***nvert)
{
    return ((nvert[k][jL][i]+nvert[k][jR][i]) < 0.1);
}

//***************************************************************************************************************//

//! \brief check is a fluid i face - shared by 2 fluid cells
inline bool isFluidIFace (PetscInt k, PetscInt j, PetscInt iL, PetscInt iR, PetscReal ***nvert)
{
    return ((nvert[k][j][iL]+nvert[k][j][iR]) < 0.1);
}

//***************************************************************************************************************//

//! \brief check is a IBM k face - shared by atleast one IBM cell
inline bool isIBMKFace (PetscInt kL, PetscInt j, PetscInt i, PetscInt kR, PetscReal ***nvert)
{
    return ((nvert[kL][j][i]+nvert[kR][j][i]) > 0.1);
}

//***************************************************************************************************************//

//! \brief check is a IBM j face - shared by atleast one IBM cell
inline bool isIBMJFace (PetscInt k, PetscInt jL, PetscInt i, PetscInt jR, PetscReal ***nvert)
{
    return ((nvert[k][jL][i]+nvert[k][jR][i]) > 0.1);
}

//***************************************************************************************************************//

//! \brief check is a IBM i face - shared by atleast one IBM cell
inline bool isIBMIFace (PetscInt k, PetscInt j, PetscInt iL, PetscInt iR, PetscReal ***nvert)
{
    return ((nvert[k][j][iL]+nvert[k][j][iR]) > 0.1);
}

//***************************************************************************************************************//

//! \brief check if IBM fluid k face - atleast one of the shared cells is IBM fluid
inline bool isIBMFluidKFace (PetscInt kL, PetscInt j, PetscInt i, PetscInt kR, PetscReal ***nvert)
{
    return ((nvert[kL][j][i] + nvert[kR][j][i] > 0.1) && (nvert[kL][j][i] + nvert[kR][j][i] < 2.5));
}

//***************************************************************************************************************//

//! \brief check if IBM fluid j face - at least one of the shared cells is IBM fluid
inline bool isIBMFluidJFace (PetscInt k, PetscInt jL, PetscInt i, PetscInt jR, PetscReal ***nvert)
{
    return ((nvert[k][jL][i] + nvert[k][jR][i] > 0.1) && (nvert[k][jL][i] + nvert[k][jR][i] < 2.5));
}

//***************************************************************************************************************//

//! \brief check if IBM fluid i face - at least one of the shared cells is IBM fluid
inline bool isIBMFluidIFace (PetscInt k, PetscInt j, PetscInt iL, PetscInt iR, PetscReal ***nvert)
{
    return ((nvert[k][j][iL] + nvert[k][j][iR] > 0.1) && (nvert[k][j][iL] + nvert[k][j][iR] < 2.5));
}

//***************************************************************************************************************//

// IBM solid faces - between ibm solid cells or between ibm solid and ibm fluid cell
inline bool isIBMSolidKFace (PetscInt kL, PetscInt j, PetscInt i, PetscInt kR, PetscReal ***nvert)
{
    return ((nvert[kL][j][i] + nvert[kR][j][i] > 3.5));
}

//***************************************************************************************************************//

//! \brief check if IBM fluid j face - at least one of the shared cells is IBM fluid
inline bool isIBMSolidJFace (PetscInt k, PetscInt jL, PetscInt i, PetscInt jR, PetscReal ***nvert)
{
    return ((nvert[k][jL][i] + nvert[k][jR][i] > 3.5));
}

//***************************************************************************************************************//

//! \brief check if IBM fluid i face - at least one of the shared cells is IBM fluid
inline bool isIBMSolidIFace (PetscInt k, PetscInt j, PetscInt iL, PetscInt iR, PetscReal ***nvert)
{
    return ((nvert[k][j][iL] + nvert[k][j][iR] > 3.5));
}
//***************************************************************************************************************//

// check if iverset cell within the neighbouring box of i,j,k
inline bool isBoxOversetCell (int k, int j, int i, PetscReal ***meshTag)
{
    return (meshTag[k  ][j  ][i  ] +
            meshTag[k+1][j  ][i  ] +
            meshTag[k-1][j  ][i  ] +
            meshTag[k  ][j+1][i  ] +
            meshTag[k  ][j-1][i  ] +
            meshTag[k  ][j  ][i+1] +
            meshTag[k  ][j  ][i-1] > 0.1);
}
//***************************************************************************************************************//

//! \brief check if the given cell is interpolated
inline bool isInterpolatedCell (int k, int j, int i, PetscReal ***meshTag)
{
    return ((meshTag[k][j][i] > 0.1) && (meshTag[k][j][i] < 1.5));
}
//***************************************************************************************************************//

//! \brief check if the given cell is set to 0 - interior of overset mesh
inline bool isZeroedCell (int k, int j, int i, PetscReal ***meshTag)
{
    return (meshTag[k][j][i] > 1.5);
}
//***************************************************************************************************************//

inline bool isOversetCell (int k, int j, int i, PetscReal ***meshTag)
{
    return (meshTag[k][j][i] > 0.1);
}
//***************************************************************************************************************//

inline bool isCalculatedCell (int k, int j, int i, PetscReal ***meshTag)
{
    return (meshTag[k][j][i] < 0.1);
}
//***************************************************************************************************************//

//! \brief check if the given face is interpolated
inline bool isInterpolatedKFace (int kL, int j, int i, int kR, PetscReal ***meshTag)
{
    return ((meshTag[kL][j][i] + meshTag[kR][j][i] > 0.1) && (meshTag[kL][j][i] + meshTag[kR][j][i] < 2.5));
}
//***************************************************************************************************************//

inline bool isInterpolatedJFace (int k, int jL, int i, int jR, PetscReal ***meshTag)
{
    return ((meshTag[k][jL][i] + meshTag[k][jR][i] > 0.1) && (meshTag[k][jL][i] + meshTag[k][jR][i] < 2.5));
}
//***************************************************************************************************************//

inline bool isInterpolatedIFace (int k, int j, int iL, int iR, PetscReal ***meshTag)
{
    return ((meshTag[k][j][iL] + meshTag[k][j][iR] > 0.1) && (meshTag[k][j][iL] + meshTag[k][j][iR] < 2.5));
}
//***************************************************************************************************************//

//! \brief check if the given face is shared by a zeroed cell
inline bool isZeroedKFace (int kL, int j, int i, int kR, PetscReal ***meshTag)
{
    return ((meshTag[kL][j][i]+meshTag[kR][j][i]) > 2.5);
}
//***************************************************************************************************************//

inline bool isZeroedJFace (int k, int jL, int i, int jR, PetscReal ***meshTag)
{
    return ((meshTag[k][jL][i]+meshTag[k][jR][i]) > 2.5);
}
//***************************************************************************************************************//

inline bool isZeroedIFace (int k, int j, int iL, int iR, PetscReal ***meshTag)
{
    return ((meshTag[k][j][iL]+meshTag[k][j][iR]) > 2.5);
}
//***************************************************************************************************************//

//! \brief check if the given face is shared by an overset cell
inline bool isOversetKFace (int kL, int j, int i, int kR, PetscReal ***meshTag)
{
    return ((meshTag[kL][j][i]+meshTag[kR][j][i]) > 0.1);
}
//***************************************************************************************************************//

inline bool isOversetJFace (int k, int jL, int i, int jR, PetscReal ***meshTag)
{
    return ((meshTag[k][jL][i]+meshTag[k][jR][i]) > 0.1);
}
//***************************************************************************************************************//

inline bool isOversetIFace (int k, int j, int iL, int iR, PetscReal ***meshTag)
{
    return ((meshTag[k][j][iL]+meshTag[k][j][iR]) > 0.1);
}
//***************************************************************************************************************//

//! \brief check if the given face is shared by only calculated cells
inline bool isCalculatedKFace (int kL, int j, int i, int kR, PetscReal ***meshTag)
{
    return ((meshTag[kL][j][i]+meshTag[kR][j][i]) < 0.1);
}
//***************************************************************************************************************//

inline bool isCalculatedJFace (int k, int jL, int i, int jR, PetscReal ***meshTag)
{
    return ((meshTag[k][jL][i]+meshTag[k][jR][i]) < 0.1);
}
//***************************************************************************************************************//

inline bool isCalculatedIFace (int k, int j, int iL, int iR, PetscReal ***meshTag)
{
    return ((meshTag[k][j][iL]+meshTag[k][j][iR]) < 0.1);
}
//***************************************************************************************************************//

//! \brief Resets the value at the non-solved centers to zero 
inline void resetNonResolvedCellCentersScalar(mesh_ *mesh,  Vec &V)
{
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;

    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         i, j, k;

    PetscReal        ***nvert, ***s, ***meshTag;

    DMDAVecGetArray(mesh->da, mesh->lNvert, &nvert);
    DMDAVecGetArray(mesh->da, mesh->lmeshTag, &meshTag);
    DMDAVecGetArray(mesh->da, V, &s);

    // Resets to zero the values of a scalar field at the non-resolved
    // faces of the mesh. Note: doesn't scatter to local.


    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
                if
                (
                    i==0 ||
                    i==mx-1 ||
                    j==0 ||
                    j==my-1 ||
                    k==0 ||
                    k==mz-1 ||
                    isIBMCell(k, j, i, nvert) || 
                    isOversetCell(k, j, i, meshTag)
                )
                {
                    s[k][j][i] = 0.0;
                }
            }
        }
    }

    DMDAVecRestoreArray(mesh->da, V, &s);
    DMDAVecRestoreArray(mesh->da, mesh->lNvert, &nvert);
    DMDAVecRestoreArray(mesh->da, mesh->lmeshTag, &meshTag);

    return;

}

//***************************************************************************************************************//

//! \brief Interpolates cartesian velocity at cell center from the contravariant flux around a cell
inline void ContravariantToCartesianPoint(Cmpnts &csi, Cmpnts &eta, Cmpnts &zet, Cmpnts &ucont, Cmpnts *ucat)
{
    // solving the linear system using Cramer method
    PetscReal det
           =
           csi.x * (eta.y * zet.z - eta.z * zet.y) -
           csi.y * (eta.x * zet.z - eta.z * zet.x) +
           csi.z * (eta.x * zet.y - eta.y * zet.x);

    PetscReal det0
           =
           ucont.x * (eta.y * zet.z - eta.z * zet.y) -
           ucont.y * (csi.y * zet.z - csi.z * zet.y) +
           ucont.z * (csi.y * eta.z - csi.z * eta.y);

    PetscReal det1
           =
          -ucont.x * (eta.x * zet.z - eta.z * zet.x) +
           ucont.y * (csi.x * zet.z - csi.z * zet.x) -
           ucont.z * (csi.x * eta.z - csi.z * eta.x);

    PetscReal det2
           =
           ucont.x * (eta.x * zet.y - eta.y * zet.x) -
           ucont.y * (csi.x * zet.y - csi.y * zet.x) +
           ucont.z * (csi.x * eta.y - csi.y * eta.x);

    (*ucat).x = det0 / det;
    (*ucat).y = det1 / det;
    (*ucat).z = det2 / det;

    return;
};

//***************************************************************************************************************//
//! \brief Resets the value at the non-solved faces to zero except for zeroGradient, overset and periodic
inline void resetNonResolvedCellFaces(mesh_ *mesh, Vec &V)
{
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;

    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         i, j, k;

    PetscReal        ***nvert, ***meshTag;
    Cmpnts           ***v;

    DMDAVecGetArray(mesh->fda, V, &v);
    DMDAVecGetArray(mesh->da,  mesh->lNvert, &nvert);
    DMDAVecGetArray(mesh->da,  mesh->lmeshTag, &meshTag);

    // Resets to zero the values of a vector field at the non-resolved
    // faces of the mesh. Note: doesn't scatter to local.


    for (k=zs; k<ze; k++)
    {
        for (j=ys; j<ye; j++)
        {
            for (i=xs; i<xe; i++)
            {
              //iface
              if(i==0)
              {
                if(mesh->boundaryU.iLeft=="zeroGradient")
                {
                    //
                }
                else
                {
                  v[k][j][i].x = 0.0;
                }
              }
              else if(i==mx-1)
              {
                  v[k][j][i].x = 0.0;
              }
              else if(i==mx-2 && !(mesh->i_periodic) && !(mesh->ii_periodic))
              {
                if(mesh->boundaryU.iRight=="zeroGradient")
                {
                  //
                }
                else
                {
                  v[k][j][i].x = 0.0;
                }
              }
              else if(isIBMIFace(k, j, i, i+1, nvert))
              {
                  v[k][j][i].x = 0.0;
              }
              else if(isOversetIFace(k, j, i, i+1, meshTag))
              {
                  v[k][j][i].x = 0.0;
              }

              //jface
              if(j==0)
              {
                if(mesh->boundaryU.jLeft=="zeroGradient")
                {
                    //
                }
                else
                {
                  v[k][j][i].y = 0.0;
                }
              }
              else if(j==my-1)
              {
                  v[k][j][i].y = 0.0;
              }
              else if(j==my-2 && !(mesh->j_periodic) && !(mesh->jj_periodic))
              {
                if(mesh->boundaryU.jRight=="zeroGradient")
                {
                  //
                }
                else
                {
                  v[k][j][i].y = 0.0;
                }
              }
              else if(isIBMJFace(k, j, i, j+1, nvert))
              {
                  v[k][j][i].y = 0.0;
              }
              else if(isOversetJFace(k, j, i, j+1, meshTag))
              {
                  v[k][j][i].y = 0.0;
              }

              //kface
              if(k==0)
              {
                if(mesh->boundaryU.kLeft=="zeroGradient")
                {
                    //
                }
                else
                {
                  v[k][j][i].z = 0.0;
                }
              }
              else if(k==mz-1)
              {
                  v[k][j][i].z = 0.0;
              }
              else if(k==mz-2 && !(mesh->k_periodic) && !(mesh->kk_periodic))
              {
                if(mesh->boundaryU.kRight=="zeroGradient")
                {
                  //
                }
                else
                {
                  v[k][j][i].z = 0.0;
                }
              }
              else if(isIBMKFace(k, j, i, k+1, nvert))
              {
                  v[k][j][i].z = 0.0;
              }
              else if(isOversetKFace(k, j, i, k+1, meshTag))
              {
                  v[k][j][i].z = 0.0;
              }

            }
        }
    }

    DMDAVecRestoreArray(mesh->fda, V, &v);
    DMDAVecRestoreArray(mesh->da,  mesh->lNvert, &nvert);
    DMDAVecRestoreArray(mesh->da,  mesh->lmeshTag, &meshTag);

    return;
}

//***************************************************************************************************************//

inline bool isOnNonResolvedCellCenters(PetscInt i, PetscInt j, PetscInt k, DMDALocalInfo info)
{
    PetscInt         mx = info.mx,
                     my = info.my,
                     mz = info.mz;

    bool             result;

    // i-left ghost plane
    if(i==0)
    {
        result = true;
    }
    // j-left ghost plane
    else if(j==0)
    {
        result = true;
    }
    // k-left ghost plane
    else if(k==0)
    {
        result = true;
    }
    // i-right ghost plane
    else if(i==mx-1)
    {
       result = true;
    }
    // j-right ghost plane
    else if(j==my-1)
    {
        result = true;
    }
    // k-right ghost plane
    else if(k==mz-1)
    {
        result = true;
    }
    // resolved cells
    else
    {
        result = false;
    }

    return(result);
}

//***************************************************************************************************************//

inline bool isOnCornerCellCenters(PetscInt i, PetscInt j, PetscInt k, DMDALocalInfo info)
{
    PetscInt         mx = info.mx,
                     my = info.my,
                     mz = info.mz;

    bool             result;

    // line at i=0, k=0
    if(i==0 && k==0)
    {
        result = true;
    }
    // line at i=mx-1, k=0
    else if(i==mx-1 && k==0)
    {
        result = true;
    }
    // line at j=0, k=0
    else if(j==0 && k==0)
    {
        result = true;
    }
    // line at j=my-1, k=0
    else if(j==my-1 && k==0)
    {
       result = true;
    }
    // line at i=0, k=mz-1
    else if(i==0 && k==mz-1)
    {
        result = true;
    }
    // line at i=mx-1, k=mz-1
    else if(i==mx-1 && k==mz-1)
    {
        result = true;
    }
    // line at j=0, k=mz-1
    else if(j==0 && k==mz-1)
    {
        result = true;
    }
    // line at j=my-1, k=mz-1
    else if(j==my-1 && k==mz-1)
    {
        result = true;
    }
    // line at i=0, j=0
    else if(i==0 && j==0)
    {
        result = true;
    }
    // line at i=mx-1, j=0
    else if(i==mx-1 && j==0)
    {
        result = true;
    }
    // line at i=0, j=my-1
    else if(i==0 && j==my-1)
    {
        result = true;
    }
    // line at i=mx-1, j=my-1
    else if(i==mx-1 && j==my-1)
    {
        result = true;
    }
    // non-corner cells
    else
    {
        result = false;
    }

    return(result);
}

// STENCILS
// ============================================================================================================= //

//! \brief Get cell's idxs of a 3 cell-sencil in csi direction given the cell idx. NEVER USE GHOSTS
inline void getCell2Cell3StencilCsiNoGhost(mesh_ *mesh, PetscInt i, PetscInt mx, PetscInt *iL, PetscInt *iR)
{
    // apply zero gradient at boundaries
    if(i == 1)
    {
        *iL = i;
        *iR = i + 1;
    }
    if(i == mx-2)
    {
        *iL = i - 1;
        *iR = i;
    }
    else
    {
        *iL = i-1, *iR = i+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get cell's idxs of a 3 cell-sencil in eta direction given the cell idx. NEVER USE GHOSTS
inline void getCell2Cell3StencilEtaNoGhost(mesh_ *mesh, PetscInt j, PetscInt my, PetscInt *jL, PetscInt *jR)
{
    // apply zero gradient at boundaries
    if(j == 1)
    {
        *jL = j;
        *jR = j + 1;
    }
    else if(j == my-2)
    {
        *jL = j-1;
        *jR = j;
    }
    else
    {
        *jL = j-1, *jR = j+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get cell's idxs of a 3 cell-sencil in zet direction given the cell idx. NEVER USE GHOSTS
inline void getCell2Cell3StencilZetNoGhost(mesh_ *mesh, PetscInt k, PetscInt mz, PetscInt *kL, PetscInt *kR)
{
    // apply zero gradient at boundaries
    if(k == 1)
    {
        *kL = k;
        *kR = k + 1;
    }
    else if(k == mz-2)
    {
        *kL = k-1;
        *kR = k;
    }
    else
    {
        *kL = k-1, *kR = k+1;
    }

    return;
}

//! \brief Get cell's idxs of a 3 cell-sencil in csi direction given the cell idx
inline void getCell2Cell3StencilCsi(mesh_ *mesh, PetscInt i, PetscInt mx, PetscInt *iL, PetscInt *iR)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(i==1 || i==mx-2)
    {
        if(mesh->i_periodic  && i==mx-2)
        {
            *iL = mx-3, *iR = 1;
        }
        if(mesh->i_periodic  && i==1)
        {
            *iL = mx-2, *iR = 2;
        }
        else if(mesh->ii_periodic && i==mx-2)
        {
            *iR = mx+1;
            *iL = i-1;
        }
        else if(mesh->ii_periodic && i==1)
        {
            *iL = -2;
            *iR =  2;
        }
        else
        {
            if(i == 1)
            {
                *iL = i;
                *iR = i + 1;
            }

            if(i == mx-2)
            {
                *iL = i - 1;
                *iR = i;
            }
        }
    }

    else
    {
        *iL = i-1, *iR = i+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get cell's idxs of a 3 cell-sencil in eta direction given the cell idx
inline void getCell2Cell3StencilEta(mesh_ *mesh, PetscInt j, PetscInt my, PetscInt *jL, PetscInt *jR)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(j==1 || j==my-2)
    {
        if(mesh->j_periodic  && j==my-2)
        {
            *jL = my-3, *jR = 1;
        }
        if(mesh->j_periodic  && j==1)
        {
            *jL = my-2, *jR = 2;
        }
        else if(mesh->jj_periodic && j==my-2)
        {
            *jR = my+1;
            *jL = j-1;
        }
        else if(mesh->jj_periodic && j==1)
        {
            *jL= -2;
            *jR = 2;
        }
        else
        {
            if(j == 1)
            {
                *jL = j;
                *jR = j + 1;
            }

            if(j == my-2)
            {
                *jL = j-1;
                *jR = j;
            }
        }
    }

    else
    {
        *jL = j-1, *jR = j+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get cell's idxs of a 3 cell-sencil in zet direction given the cell idx
inline void getCell2Cell3StencilZet(mesh_ *mesh, PetscInt k, PetscInt mz, PetscInt *kL, PetscInt *kR)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(k==1 || k==mz-2)
    {
        if(mesh->k_periodic && k==mz-2)
        {
            *kL = mz-3, *kR = 1;
        }
        if(mesh->k_periodic && k==1)
        {
            *kL = mz-2, *kR = 2;
        }
        else if(mesh->kk_periodic && k==mz-2)
        {
            *kR = mz+1;
            *kL = k-1;
        }
        else if(mesh->kk_periodic && k==1)
        {
            *kL = -2;
            *kR =  2;
        }
        else
        {
            if(k == 1)
            {
                *kL = k;
                *kR = k + 1;
            }

            if(k == mz-2)
            {
                *kL = k-1;
                *kR = k;
            }
        }
    }

    else
    {
        *kL = k-1, *kR = k+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get outmost cell's idxs of a 4 cell-sencil in csi direction given the face idx
inline void getFace2Cell4StencilCsi(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt mx, PetscInt *iL, PetscInt *iR, PetscReal *denom, PetscReal ***nvert, PetscReal ***meshTag)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(i==0 || i==mx-2)
    {
        if(mesh->i_periodic)
        {
            *iL = mx-3, *iR = 2;
            *denom = 3.;
        }
        else if(mesh->ii_periodic && i==mx-2)
        {
            *iR=mx+2;
            *iL = i-1;
            *denom = 3.;
        }
        else if(mesh->ii_periodic && i==0)
        {
            *iL=-3;
            *iR = i+2;
            *denom = 3.;
        }
        else
        {
          if(i == 0)
          {
              *iL = i;
              *iR = i + 2;
              *denom = 3.;
          }

          if(i == mx-2)
          {
              *iL = i-1;
              *iR = i + 1;
              *denom = 3.;
          }
        }
    }
    else if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        *iL = i;
        *iR = i + 2;
        *denom = 3.; 
    }
    else if(isIBMCell(k, j, i+1, nvert) || isOversetCell(k, j, i+1, meshTag))
    {
        *iL = i-1;
        *iR = i + 1;
        *denom = 3.;
    }
    else
    {
        *iL = i-1, *iR = i+2;
        *denom = 3.;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get outmost cell's idxs of a 4 cell-sencil in eta direction given the face idx
inline void getFace2Cell4StencilEta(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt my, PetscInt *jL, PetscInt *jR, PetscReal *denom, PetscReal ***nvert, PetscReal ***meshTag)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(j==0 || j==my-2)
    {
        if(mesh->j_periodic)
        {
            *jL = my-3, *jR = 2;
            *denom = 3.;
        }
        else if(mesh->jj_periodic && j==my-2)
        {
            *jR=my+2;
            *jL = j-1;
            *denom = 3.;
        }
        else if(mesh->jj_periodic && j==0)
        {
            *jL=-3;
            *jR = j+2;
            *denom = 3.;
        }
        else
        {
          if(j == 0)
          {
              *jL = j;
              *jR = j + 2;
              *denom = 3.;
          }

          if(j == my-2)
          {
              *jL = j-1;
              *jR = j + 1;
              *denom = 3.;
          }
        }
    }
    else if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        *jL = j;
        *jR = j + 2;
        *denom = 3.; 
    }
    else if(isIBMCell(k, j+1, i, nvert) || isOversetCell(k, j+1, i, meshTag))
    {
        *jL = j-1;
        *jR = j + 1;
        *denom = 3.;
    }
    else
    {
        *jL = j-1, *jR = j+2;
        *denom = 3.;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get outmost cell's idxs of a 4 cell-sencil in zet direction given the face idx
inline void getFace2Cell4StencilZet(mesh_ *mesh, PetscInt k, PetscInt j, PetscInt i, PetscInt mz, PetscInt *kL, PetscInt *kR, PetscReal *denom, PetscReal ***nvert, PetscReal ***meshTag)
{
    // if on non-periodic boundary return 2 cell-stencil extrema idxs, if periodic
    // get the right data.
    if(k==0 || k==mz-2)
    {
        if(mesh->k_periodic)
        {
            *kL = mz-3, *kR = 2;
            *denom = 3.;
        }
        else if(mesh->kk_periodic && k==mz-2)
        {
            *kR=mz+2;
            *kL = k-1;
            *denom = 3.;
        }
        else if(mesh->kk_periodic && k==0)
        {
            *kL=-3;
            *kR = k+2;
            *denom = 3.;
        }
        else
        {
          if(k == 0)
          {
              *kL = k;
              *kR = k + 2;
              *denom = 3.;
          }

          if(k == mz-2)
          {
              *kL = k-1;
              *kR = k + 1;
              *denom = 3.;
          }
        }
    }
    else if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        *kL = k;
        *kR = k + 2;
        *denom = 3.; 
    }
    else if(isIBMCell(k+1, j, i, nvert) || isOversetCell(k+1, j, i, meshTag))
    {
        *kL = k-1;
        *kR = k + 1;
        *denom = 3.;
    }
    else
    {
        *kL = k-1, *kR = k+2;
        *denom = 3.;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get neighboring face idxs of a 3 face-sencil in csi direction given the central face idx
inline void getFace2Face3StencilCsi(mesh_ *mesh, PetscInt i, PetscInt mx, PetscInt *iL, PetscInt *iR)
{
    // cell faces don't have ghosts, if on non-periodic boundary return the central
    // face label as the ghost face label. If periodic get neighboring face idxs.
    if(i==0 || i==mx-2)
    {
        if(mesh->i_periodic)
        {
            *iL = mx-3, *iR = 1;
        }
        else if(mesh->ii_periodic && i==mx-2)
        {
            *iR=mx+1;
            *iL = i-1;
        }
        else if(mesh->ii_periodic && i==0)
        {
            *iL=-3;
            *iR = i+1;
        }
        else if(i==0)
        {
            *iL = i;
            *iR = i + 1;
        }
        else if(i==mx-2)
        {
            *iL = i - 1;
            *iR = i;
        }
    }
    else
    {
        *iL = i-1, *iR = i+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get neighboring face idxs of a 3 face-sencil in eta direction given the central face idx
inline void getFace2Face3StencilEta(mesh_ *mesh, PetscInt j, PetscInt my, PetscInt *jL, PetscInt *jR)
{
    // cell faces don't have ghosts, if on non-periodic boundary return the central
    // face label as the ghost face label. If periodic get neighboring face idxs.
    if(j==0 || j==my-2)
    {
        if(mesh->j_periodic)
        {
            *jL = my-3, *jR = 1;
        }
        else if(mesh->jj_periodic && j==my-2)
        {
            *jR=my+1;
            *jL = j-1;
        }
        else if(mesh->jj_periodic && j==0)
        {
            *jL=-3;
            *jR = j+1;
        }
        else if(j==0)
        {
            *jL = j;
            *jR = j + 1;
        }
        else if(j==my-2)
        {
            *jL = j - 1;
            *jR = j;
        }
    }
    else
    {
        *jL = j-1, *jR = j+1;
    }

    return;
}

//***************************************************************************************************************//

//! \brief Get neighboring face idxs of a 3 face-sencil in zeta direction given the central face idx
inline void getFace2Face3StencilZet(mesh_ *mesh, PetscInt k, PetscInt mz, PetscInt *kL, PetscInt *kR)
{
    // cell faces don't have ghosts, if on non-periodic boundary return the central
    // face label as the ghost face label. If periodic get neighboring face idxs.
    if(k==0 || k==mz-2)
    {
        if(mesh->k_periodic)
        {
            *kL = mz-3, *kR = 1;
        }
        else if(mesh->kk_periodic && k==mz-2)
        {
            *kR=mz+1;
            *kL = k-1;
        }
        else if(mesh->kk_periodic && k==0)
        {
            *kL=-3;
            *kR = k+1;
        }
        else if(k==0)
        {
            *kL = k;
            *kR = k + 1;
        }
        else if(k==mz-2)
        {
            *kL = k - 1;
            *kR = k;
        }
    }
    else
    {
        *kL = k-1, *kR = k+1;
    }

    return;
}

// GRADIENT SCHEMES
// ============================================================================================================= //

inline void Compute_du_i
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)
{
    *dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
    *dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
    *dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;

    if (isIBMSolidIFace(k, j+1, i, i+1, nvert) || isZeroedIFace(k, j+1, i, i+1, meshTag) || 
            (j==my-2 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dude = (ucat[k][j  ][i+1].x + ucat[k][j  ][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
        *dvde = (ucat[k][j  ][i+1].y + ucat[k][j  ][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
        *dwde = (ucat[k][j  ][i+1].z + ucat[k][j  ][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
    }
    else if  (isIBMSolidIFace(k, j-1, i, i+1, nvert) || isZeroedIFace(k, j-1, i, i+1, meshTag) || 
            (j==1 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j  ][i+1].x - ucat[k][j  ][i].x) * 0.5;
        *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j  ][i+1].y - ucat[k][j  ][i].y) * 0.5;
        *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j  ][i+1].z - ucat[k][j  ][i].z) * 0.5;
    }
    else
    {
        *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
        *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
        *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
    }

    if (isIBMSolidIFace(k+1, j, i, i+1, nvert) || isZeroedIFace(k+1, j, i, i+1, meshTag) || 
            (k==mz-2 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dudz = (ucat[k  ][j][i+1].x + ucat[k  ][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
        *dvdz = (ucat[k  ][j][i+1].y + ucat[k  ][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
        *dwdz = (ucat[k  ][j][i+1].z + ucat[k  ][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
    }
    else if (isIBMSolidIFace(k-1, j, i, i+1, nvert) || isZeroedIFace(k-1, j, i, i+1, meshTag) || 
            (k==1 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k  ][j][i+1].x - ucat[k  ][j][i].x) * 0.5;
        *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k  ][j][i+1].y - ucat[k  ][j][i].y) * 0.5;
        *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k  ][j][i+1].z - ucat[k  ][j][i].z) * 0.5;
    }
    else
    {
        *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
        *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
        *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
    }

    return;
}

//***************************************************************************************************************//

inline void Compute_du_j
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)

{
    if (isIBMSolidJFace(k, j, i+1, j+1, nvert) || isZeroedJFace(k, j, i+1, j+1, meshTag) ||
          (i==mx-2 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudc = (ucat[k][j+1][i  ].x + ucat[k][j][i  ].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
        *dvdc = (ucat[k][j+1][i  ].y + ucat[k][j][i  ].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
        *dwdc = (ucat[k][j+1][i  ].z + ucat[k][j][i  ].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
    }
    else if (isIBMSolidJFace(k, j, i-1, j+1, nvert) || isZeroedJFace(k, j, i-1, j+1, meshTag) || 
            (i==1 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i  ].x - ucat[k][j][i  ].x) * 0.5;
        *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i  ].y - ucat[k][j][i  ].y) * 0.5;
        *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i  ].z - ucat[k][j][i  ].z) * 0.5;
    }
    else
    {
        *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
        *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
        *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
    }

    *dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
    *dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
    *dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;

    if (isIBMSolidJFace(k+1, j, i, j+1, nvert) || isZeroedJFace(k+1, j, i, j+1, meshTag) || 
            (k==mz-2 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudz = (ucat[k  ][j+1][i].x + ucat[k  ][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
        *dvdz = (ucat[k  ][j+1][i].y + ucat[k  ][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
        *dwdz = (ucat[k  ][j+1][i].z + ucat[k  ][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
    }
    else if (isIBMSolidJFace(k-1, j, i, j+1, nvert) || isZeroedJFace(k-1, j, i, j+1, meshTag) ||
            (k==1 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k  ][j+1][i].x - ucat[k  ][j][i].x) * 0.5;
        *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k  ][j+1][i].y - ucat[k  ][j][i].y) * 0.5;
        *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k  ][j+1][i].z - ucat[k  ][j][i].z) * 0.5;
    }
    else
    {
        *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
        *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
        *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
    }

    return;
}

//***************************************************************************************************************//

inline void Compute_du_k
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)

{
    if (isIBMSolidKFace(k, j, i+1, k+1, nvert) || isZeroedKFace(k, j, i+1, k+1, meshTag) || 
          (i==mx-2 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dudc = (ucat[k+1][j][i  ].x + ucat[k][j][i  ].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.5;
        *dvdc = (ucat[k+1][j][i  ].y + ucat[k][j][i  ].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.5;
        *dwdc = (ucat[k+1][j][i  ].z + ucat[k][j][i  ].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.5;
    }
    else if (isIBMSolidKFace(k, j, i-1, k+1, nvert) || isZeroedKFace(k, j, i-1, k+1, meshTag) ||
            (i==1 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i  ].x - ucat[k][j][i  ].x) * 0.5;
        *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i  ].y - ucat[k][j][i  ].y) * 0.5;
        *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i  ].z - ucat[k][j][i  ].z) * 0.5;
    }
    else
    {
        *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x) * 0.25;
        *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y) * 0.25;
        *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z) * 0.25;
    }

    if (isIBMSolidKFace(k, j+1, i, k+1, nvert) || isZeroedKFace(k, j+1, i, k+1, meshTag) || 
            (j==my-2 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dude = (ucat[k+1][j  ][i].x + ucat[k][j  ][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.5;
        *dvde = (ucat[k+1][j  ][i].y + ucat[k][j  ][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.5;
        *dwde = (ucat[k+1][j  ][i].z + ucat[k][j  ][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.5;
    }
    else if (isIBMSolidKFace(k, j-1, i, k+1, nvert) || isZeroedKFace(k, j-1, i, k+1, meshTag) || 
            (j==1 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j  ][i].x - ucat[k][j  ][i].x) * 0.5;
        *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j  ][i].y - ucat[k][j  ][i].y) * 0.5;
        *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j  ][i].z - ucat[k][j  ][i].z) * 0.5;
    }
    else
    {
        *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x) * 0.25;
        *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y) * 0.25;
        *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z) * 0.25;
    }

    *dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
    *dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
    *dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;

    return;
}

//***************************************************************************************************************//

inline void Compute_du_i_4
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)
{
    if(i >= 1 && i <= mx - 3) 
    {
        *dudc = (-ucat[k][j][i+2].x + 7.0 * ucat[k][j][i+1].x - 7.0 * ucat[k][j][i].x + ucat[k][j][i-1].x) / 12.0;
        *dvdc = (-ucat[k][j][i+2].y + 7.0 * ucat[k][j][i+1].y - 7.0 * ucat[k][j][i].y + ucat[k][j][i-1].y) / 12.0;
        *dwdc = (-ucat[k][j][i+2].z + 7.0 * ucat[k][j][i+1].z - 7.0 * ucat[k][j][i].z + ucat[k][j][i-1].z) / 12.0;
    } 

    else 
    {
        *dudc = ucat[k][j][i+1].x - ucat[k][j][i].x;
        *dvdc = ucat[k][j][i+1].y - ucat[k][j][i].y;
        *dwdc = ucat[k][j][i+1].z - ucat[k][j][i].z;
    }

    if(j >= 1 && j <= my - 3 
       && !(isIBMSolidIFace(k, j+1, i, i+1, nvert) || isZeroedIFace(k, j+1, i, i+1, meshTag))
       && !(isIBMSolidIFace(k, j-1, i, i+1, nvert) || isZeroedIFace(k, j-1, i, i+1, meshTag)) )
    {
        *dude = (- (ucat[k][j+2][i+1].x + ucat[k][j+2][i].x) + 7.0 * (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x) - 7.0 * (ucat[k][j][i+1].x + ucat[k][j][i].x) + (ucat[k][j-1][i+1].x + ucat[k][j-1][i].x)) / 24.0;
        *dvde = (- (ucat[k][j+2][i+1].y + ucat[k][j+2][i].y) + 7.0 * (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y) - 7.0 * (ucat[k][j][i+1].y + ucat[k][j][i].y) + (ucat[k][j-1][i+1].y + ucat[k][j-1][i].y)) / 24.0;
        *dwde = (- (ucat[k][j+2][i+1].z + ucat[k][j+2][i].z) + 7.0 * (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z) - 7.0 * (ucat[k][j][i+1].z + ucat[k][j][i].z) + (ucat[k][j-1][i+1].z + ucat[k][j-1][i].z)) / 24.0;
    } 
    
    else if (isIBMSolidIFace(k, j+1, i, i+1, nvert) || isZeroedIFace(k, j+1, i, i+1, meshTag) || 
    (j==my-2 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dude = (ucat[k][j][i+1].x + ucat[k][j][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.5;
        *dvde = (ucat[k][j][i+1].y + ucat[k][j][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.5;
        *dwde = (ucat[k][j][i+1].z + ucat[k][j][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.5;
    } 
    
    else if  (isIBMSolidIFace(k, j-1, i, i+1, nvert) || isZeroedIFace(k, j-1, i, i+1, meshTag) || 
            (j==1 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j][i+1].x - ucat[k][j][i].x) * 0.5;
        *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j][i+1].y - ucat[k][j][i].y) * 0.5;
        *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j][i+1].z - ucat[k][j][i].z) * 0.5;
    }

    else 
    {
        *dude = (ucat[k][j+1][i+1].x + ucat[k][j+1][i].x - ucat[k][j-1][i+1].x - ucat[k][j-1][i].x) * 0.25;
        *dvde = (ucat[k][j+1][i+1].y + ucat[k][j+1][i].y - ucat[k][j-1][i+1].y - ucat[k][j-1][i].y) * 0.25;
        *dwde = (ucat[k][j+1][i+1].z + ucat[k][j+1][i].z - ucat[k][j-1][i+1].z - ucat[k][j-1][i].z) * 0.25;
    }

    if(k >= 1 && k <= mz - 3 
        && !(isIBMSolidIFace(k+1, j, i, i+1, nvert) || isZeroedIFace(k+1, j, i, i+1, meshTag))
        && !(isIBMSolidIFace(k-1, j, i, i+1, nvert) || isZeroedIFace(k-1, j, i, i+1, meshTag)) )
    {
        *dudz = (- (ucat[k+2][j][i+1].x + ucat[k+2][j][i].x) + 7.0 * (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x) - 7.0 * (ucat[k][j][i+1].x + ucat[k][j][i].x) + (ucat[k-1][j][i+1].x + ucat[k-1][j][i].x)) / 24.0;
        *dvdz = (- (ucat[k+2][j][i+1].y + ucat[k+2][j][i].y) + 7.0 * (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y) - 7.0 * (ucat[k][j][i+1].y + ucat[k][j][i].y) + (ucat[k-1][j][i+1].y + ucat[k-1][j][i].y)) / 24.0;
        *dwdz = (- (ucat[k+2][j][i+1].z + ucat[k+2][j][i].z) + 7.0 * (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z) - 7.0 * (ucat[k][j][i+1].z + ucat[k][j][i].z) + (ucat[k-1][j][i+1].z + ucat[k-1][j][i].z)) / 24.0;
    } 
    
    else if(isIBMSolidIFace(k+1, j, i, i+1, nvert) || isZeroedIFace(k+1, j, i, i+1, meshTag) || 
    (k==mz-2 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dudz = (ucat[k][j][i+1].x + ucat[k][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.5;
        *dvdz = (ucat[k][j][i+1].y + ucat[k][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.5;
        *dwdz = (ucat[k][j][i+1].z + ucat[k][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.5;
    }

    else if (isIBMSolidIFace(k-1, j, i, i+1, nvert) || isZeroedIFace(k-1, j, i, i+1, meshTag) || 
            (k==1 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k][j][i+1].x - ucat[k][j][i].x) * 0.5;
        *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k][j][i+1].y - ucat[k][j][i].y) * 0.5;
        *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k][j][i+1].z - ucat[k][j][i].z) * 0.5;
    } 
    
    else 
    {
        *dudz = (ucat[k+1][j][i+1].x + ucat[k+1][j][i].x - ucat[k-1][j][i+1].x - ucat[k-1][j][i].x) * 0.25;
        *dvdz = (ucat[k+1][j][i+1].y + ucat[k+1][j][i].y - ucat[k-1][j][i+1].y - ucat[k-1][j][i].y) * 0.25;
        *dwdz = (ucat[k+1][j][i+1].z + ucat[k+1][j][i].z - ucat[k-1][j][i+1].z - ucat[k-1][j][i].z) * 0.25;
    }

}
//***************************************************************************************************************//
inline void Compute_du_j_4
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)

{
    if(i >= 1 && i <= mx - 3 
        && !(isIBMSolidJFace(k, j, i+1, j+1, nvert) || isZeroedJFace(k, j, i+1, j+1, meshTag))
        && !(isIBMSolidJFace(k, j, i-1, j+1, nvert) || isZeroedJFace(k, j, i-1, j+1, meshTag)) )
    {
        *dudc = (- (ucat[k][j][i+2].x + ucat[k][j+1][i+2].x) + 7.0*(ucat[k][j][i+1].x + ucat[k][j+1][i+1].x) - 7.0*(ucat[k][j][i].x + ucat[k][j+1][i].x) + (ucat[k][j][i-1].x + ucat[k][j+1][i-1].x)) / 24.0;
        *dvdc = (- (ucat[k][j][i+2].y + ucat[k][j+1][i+2].y) + 7.0*(ucat[k][j][i+1].y + ucat[k][j+1][i+1].y) - 7.0*(ucat[k][j][i].y + ucat[k][j+1][i].y) + (ucat[k][j][i-1].y + ucat[k][j+1][i-1].y)) / 24.0;
        *dwdc = (- (ucat[k][j][i+2].z + ucat[k][j+1][i+2].z) + 7.0*(ucat[k][j][i+1].z + ucat[k][j+1][i+1].z) - 7.0*(ucat[k][j][i].z + ucat[k][j+1][i].z) + (ucat[k][j][i-1].z + ucat[k][j+1][i-1].z)) / 24.0;
    }

    else if (isIBMSolidJFace(k, j, i+1, j+1, nvert) || isZeroedJFace(k, j, i+1, j+1, meshTag) ||
    (i==mx-2 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) )) 
    {
        *dudc = (ucat[k][j+1][i].x + ucat[k][j][i].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.5;
        *dvdc = (ucat[k][j+1][i].y + ucat[k][j][i].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.5;
        *dwdc = (ucat[k][j+1][i].z + ucat[k][j][i].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.5;
    }

    else if (isIBMSolidJFace(k, j, i-1, j+1, nvert) || isZeroedJFace(k, j, i-1, j+1, meshTag) || 
            (i==1 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) )) 
    {
        *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i].x - ucat[k][j][i].x) * 0.5;
        *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i].y - ucat[k][j][i].y) * 0.5;
        *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i].z - ucat[k][j][i].z) * 0.5;
    } 
    
    else 
    {
        *dudc = (ucat[k][j+1][i+1].x + ucat[k][j][i+1].x - ucat[k][j+1][i-1].x - ucat[k][j][i-1].x) * 0.25;
        *dvdc = (ucat[k][j+1][i+1].y + ucat[k][j][i+1].y - ucat[k][j+1][i-1].y - ucat[k][j][i-1].y) * 0.25;
        *dwdc = (ucat[k][j+1][i+1].z + ucat[k][j][i+1].z - ucat[k][j+1][i-1].z - ucat[k][j][i-1].z) * 0.25;
    }

    if(j >= 1 && j <= my - 3) 
    {
        *dude = (-ucat[k][j+2][i].x + 7.0*ucat[k][j+1][i].x - 7.0*ucat[k][j][i].x + ucat[k][j-1][i].x)/12.0;
        *dvde = (-ucat[k][j+2][i].y + 7.0*ucat[k][j+1][i].y - 7.0*ucat[k][j][i].y + ucat[k][j-1][i].y)/12.0;
        *dwde = (-ucat[k][j+2][i].z + 7.0*ucat[k][j+1][i].z - 7.0*ucat[k][j][i].z + ucat[k][j-1][i].z)/12.0;
    } 
    else 
    {
        *dude = ucat[k][j+1][i].x - ucat[k][j][i].x;
        *dvde = ucat[k][j+1][i].y - ucat[k][j][i].y;
        *dwde = ucat[k][j+1][i].z - ucat[k][j][i].z;
    }

    if(k >= 1 && k <= mz - 3 
        && !(isIBMSolidJFace(k+1, j, i, j+1, nvert) || isZeroedJFace(k+1, j, i, j+1, meshTag))
        && !(isIBMSolidJFace(k-1, j, i, j+1, nvert) || isZeroedJFace(k-1, j, i, j+1, meshTag)) ) 
    {
        *dudz = (- (ucat[k+2][j][i].x + ucat[k+2][j+1][i].x) + 7.0*(ucat[k+1][j][i].x + ucat[k+1][j+1][i].x) - 7.0*(ucat[k][j][i].x + ucat[k][j+1][i].x) + (ucat[k-1][j][i].x + ucat[k-1][j+1][i].x)) / 24.0;
        *dvdz = (- (ucat[k+2][j][i].y + ucat[k+2][j+1][i].y) + 7.0*(ucat[k+1][j][i].y + ucat[k+1][j+1][i].y) - 7.0*(ucat[k][j][i].y + ucat[k][j+1][i].y) + (ucat[k-1][j][i].y + ucat[k-1][j+1][i].y)) / 24.0;
        *dwdz = (- (ucat[k+2][j][i].z + ucat[k+2][j+1][i].z) + 7.0*(ucat[k+1][j][i].z + ucat[k+1][j+1][i].z) - 7.0*(ucat[k][j][i].z + ucat[k][j+1][i].z) + (ucat[k-1][j][i].z + ucat[k-1][j+1][i].z)) / 24.0;
    } 
    
    else if (isIBMSolidJFace(k+1, j, i, j+1, nvert) || isZeroedJFace(k+1, j, i, j+1, meshTag) || 
    (k==mz-2 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudz = (ucat[k][j+1][i].x + ucat[k][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.5;
        *dvdz = (ucat[k][j+1][i].y + ucat[k][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.5;
        *dwdz = (ucat[k][j+1][i].z + ucat[k][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.5;
    } 
    
    else if (isIBMSolidJFace(k-1, j, i, j+1, nvert) || isZeroedJFace(k-1, j, i, j+1, meshTag) ||
    (k==1 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k][j+1][i].x - ucat[k][j][i].x) * 0.5;
        *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k][j+1][i].y - ucat[k][j][i].y) * 0.5;
        *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k][j+1][i].z - ucat[k][j][i].z) * 0.5;
    }

    else 
    {
        *dudz = (ucat[k+1][j+1][i].x + ucat[k+1][j][i].x - ucat[k-1][j+1][i].x - ucat[k-1][j][i].x) * 0.25;
        *dvdz = (ucat[k+1][j+1][i].y + ucat[k+1][j][i].y - ucat[k-1][j+1][i].y - ucat[k-1][j][i].y) * 0.25;
        *dwdz = (ucat[k+1][j+1][i].z + ucat[k+1][j][i].z - ucat[k-1][j+1][i].z - ucat[k-1][j][i].z) * 0.25;
    }
    
    return;

}
//***************************************************************************************************************//
inline void Compute_du_k_4
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)

{
    if(i >= 1 && i <= mx - 3 
        && !(isIBMSolidKFace(k, j, i+1, k+1, nvert) || isZeroedKFace(k, j, i+1, k+1, meshTag))
        && !(isIBMSolidKFace(k, j, i-1, k+1, nvert) || isZeroedKFace(k, j, i-1, k+1, meshTag)) )
    {
        *dudc = (- (ucat[k+1][j][i+2].x + ucat[k][j][i+2].x) + 7.0*(ucat[k+1][j][i+1].x + ucat[k][j][i+1].x) - 7.0*(ucat[k+1][j][i].x + ucat[k][j][i].x) + (ucat[k+1][j][i-1].x + ucat[k][j][i-1].x)) / 24.0;
        *dvdc = (- (ucat[k+1][j][i+2].y + ucat[k][j][i+2].y) + 7.0*(ucat[k+1][j][i+1].y + ucat[k][j][i+1].y) - 7.0*(ucat[k+1][j][i].y + ucat[k][j][i].y) + (ucat[k+1][j][i-1].y + ucat[k][j][i-1].y)) / 24.0;
        *dwdc = (- (ucat[k+1][j][i+2].z + ucat[k][j][i+2].z) + 7.0*(ucat[k+1][j][i+1].z + ucat[k][j][i+1].z) - 7.0*(ucat[k+1][j][i].z + ucat[k][j][i].z) + (ucat[k+1][j][i-1].z + ucat[k][j][i-1].z)) / 24.0;
    }

    else if (isIBMSolidKFace(k, j, i+1, k+1, nvert) || isZeroedKFace(k, j, i+1, k+1, meshTag) || 
    (i==mx-2 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dudc = (ucat[k+1][j][i].x + ucat[k][j][i].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x)*0.5;
        *dvdc = (ucat[k+1][j][i].y + ucat[k][j][i].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y)*0.5;
        *dwdc = (ucat[k+1][j][i].z + ucat[k][j][i].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z)*0.5;
    } 
    
    else if (isIBMSolidKFace(k, j, i-1, k+1, nvert) || isZeroedKFace(k, j, i-1, k+1, meshTag) ||
    (i==1 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i].x - ucat[k][j][i].x)*0.5;
        *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i].y - ucat[k][j][i].y)*0.5;
        *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i].z - ucat[k][j][i].z)*0.5;
    } 
    
    else 
    {
        *dudc = (ucat[k+1][j][i+1].x + ucat[k][j][i+1].x - ucat[k+1][j][i-1].x - ucat[k][j][i-1].x)*0.25;
        *dvdc = (ucat[k+1][j][i+1].y + ucat[k][j][i+1].y - ucat[k+1][j][i-1].y - ucat[k][j][i-1].y)*0.25;
        *dwdc = (ucat[k+1][j][i+1].z + ucat[k][j][i+1].z - ucat[k+1][j][i-1].z - ucat[k][j][i-1].z)*0.25;
    }
    
    if(j >= 1 && j <= my - 3 
        && !(isIBMSolidKFace(k, j+1, i, k+1, nvert) || isZeroedKFace(k, j+1, i, k+1, meshTag))
        && !(isIBMSolidKFace(k, j-1, i, k+1, nvert) || isZeroedKFace(k, j-1, i, k+1, meshTag)) ) 
    {
        *dude = (- (ucat[k+1][j+2][i].x + ucat[k][j+2][i].x) + 7.0*(ucat[k+1][j+1][i].x + ucat[k][j+1][i].x) - 7.0*(ucat[k+1][j][i].x + ucat[k][j][i].x) + (ucat[k+1][j-1][i].x + ucat[k][j-1][i].x)) / 24.0;
        *dvde = (- (ucat[k+1][j+2][i].y + ucat[k][j+2][i].y) + 7.0*(ucat[k+1][j+1][i].y + ucat[k][j+1][i].y) - 7.0*(ucat[k+1][j][i].y + ucat[k][j][i].y) + (ucat[k+1][j-1][i].y + ucat[k][j-1][i].y)) / 24.0;
        *dwde = (- (ucat[k+1][j+2][i].z + ucat[k][j+2][i].z) + 7.0*(ucat[k+1][j+1][i].z + ucat[k][j+1][i].z) - 7.0*(ucat[k+1][j][i].z + ucat[k][j][i].z) + (ucat[k+1][j-1][i].z + ucat[k][j-1][i].z)) / 24.0;
    } 

    else if (isIBMSolidKFace(k, j+1, i, k+1, nvert) || isZeroedKFace(k, j+1, i, k+1, meshTag) || 
    (j==my-2 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) )) 
    {
        *dude = (ucat[k+1][j][i].x + ucat[k][j][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x)*0.5;
        *dvde = (ucat[k+1][j][i].y + ucat[k][j][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y)*0.5;
        *dwde = (ucat[k+1][j][i].z + ucat[k][j][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z)*0.5;
    }

    else if (isIBMSolidKFace(k, j-1, i, k+1, nvert) || isZeroedKFace(k, j-1, i, k+1, meshTag) || 
            (j==1 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) )) 
    {
        *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j][i].x - ucat[k][j][i].x)*0.5;
        *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j][i].y - ucat[k][j][i].y)*0.5;
        *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j][i].z - ucat[k][j][i].z)*0.5;
    } 
    
    else 
    {
        *dude = (ucat[k+1][j+1][i].x + ucat[k][j+1][i].x - ucat[k+1][j-1][i].x - ucat[k][j-1][i].x)*0.25;
        *dvde = (ucat[k+1][j+1][i].y + ucat[k][j+1][i].y - ucat[k+1][j-1][i].y - ucat[k][j-1][i].y)*0.25;
        *dwde = (ucat[k+1][j+1][i].z + ucat[k][j+1][i].z - ucat[k+1][j-1][i].z - ucat[k][j-1][i].z)*0.25;
    }

    if(k >= 1 && k <= mz-3) 
    {
        *dudz = (-ucat[k+2][j][i].x + 7.0*ucat[k+1][j][i].x - 7.0*ucat[k][j][i].x + ucat[k-1][j][i].x) / 12.0;
        *dvdz = (-ucat[k+2][j][i].y + 7.0*ucat[k+1][j][i].y - 7.0*ucat[k][j][i].y + ucat[k-1][j][i].y) / 12.0;
        *dwdz = (-ucat[k+2][j][i].z + 7.0*ucat[k+1][j][i].z - 7.0*ucat[k][j][i].z + ucat[k-1][j][i].z) / 12.0;
    } 
    
    else 
    {
        *dudz = ucat[k+1][j][i].x - ucat[k][j][i].x;
        *dvdz = ucat[k+1][j][i].y - ucat[k][j][i].y;
        *dwdz = ucat[k+1][j][i].z - ucat[k][j][i].z;
    }

    return;

}

//***************************************************************************************************************//

inline void Compute_dscalar_i
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    PetscReal ***K, PetscReal ***nvert, PetscReal ***meshTag, 
    PetscReal *dkdc, PetscReal *dkde, PetscReal *dkdz
)
{
    *dkdc = K[k][j][i+1] - K[k][j][i];

    if (isIBMSolidIFace(k, j+1, i, i+1, nvert) || isZeroedIFace(k, j+1, i, i+1, meshTag) || 
            (j==my-2 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dkde = (K[k][j  ][i+1] + K[k][j  ][i] - K[k][j-1][i+1] - K[k][j-1][i]) * 0.5;
    }
    else if  (isIBMSolidIFace(k, j-1, i, i+1, nvert) || isZeroedIFace(k, j-1, i, i+1, meshTag) || 
            (j==1 && (i==0 || i==mx-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))    
    {
        *dkde = (K[k][j+1][i+1] + K[k][j+1][i] - K[k][j  ][i+1] - K[k][j  ][i]) * 0.5;
    }
    else
    {
        *dkde = (K[k][j+1][i+1] + K[k][j+1][i] - K[k][j-1][i+1] - K[k][j-1][i]) * 0.25;
    }
    if (isIBMSolidIFace(k+1, j, i, i+1, nvert) || isZeroedIFace(k+1, j, i, i+1, meshTag) || 
            (k==mz-2 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dkdz = (K[k  ][j][i+1] + K[k  ][j][i] - K[k-1][j][i+1] - K[k-1][j][i]) * 0.5;
    }
    else if (isIBMSolidIFace(k-1, j, i, i+1, nvert) || isZeroedIFace(k-1, j, i, i+1, meshTag) || 
            (k==1 && (i==0 || i==mx-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->i_periodic && !mesh->ii_periodic)) ))
    {
        *dkdz = (K[k+1][j][i+1] + K[k+1][j][i] - K[k  ][j][i+1] - K[k  ][j][i]) * 0.5;
    }
    else
    {
        *dkdz = (K[k+1][j][i+1] + K[k+1][j][i] - K[k-1][j][i+1] - K[k-1][j][i]) * 0.25;
    }

    return;
}

//***************************************************************************************************************//

inline void Compute_dscalar_j
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    PetscReal ***K, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dkdc, PetscReal *dkde, PetscReal *dkdz
)
{
    if (isIBMSolidJFace(k, j, i+1, j+1, nvert) || isZeroedJFace(k, j, i+1, j+1, meshTag) ||
          (i==mx-2 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dkdc = (K[k][j+1][i  ] + K[k][j][i  ] - K[k][j+1][i-1] - K[k][j][i-1]) * 0.5;
    }
    else if (isIBMSolidJFace(k, j, i-1, j+1, nvert) || isZeroedJFace(k, j, i-1, j+1, meshTag) || 
            (i==1 && (j==0 || j==my-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dkdc = (K[k][j+1][i+1] + K[k][j][i+1] - K[k][j+1][i  ] - K[k][j][i  ]) * 0.5;
    }
    else
    {
        *dkdc = (K[k][j+1][i+1] + K[k][j][i+1] - K[k][j+1][i-1] - K[k][j][i-1]) * 0.25;
    }

    *dkde = K[k][j+1][i] - K[k][j][i];

    if (isIBMSolidJFace(k+1, j, i, j+1, nvert) || isZeroedJFace(k+1, j, i, j+1, meshTag) || 
            (k==mz-2 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dkdz = (K[k  ][j+1][i] + K[k  ][j][i] - K[k-1][j+1][i] - K[k-1][j][i]) * 0.5;
    }
    else if (isIBMSolidJFace(k-1, j, i, j+1, nvert) || isZeroedJFace(k-1, j, i, j+1, meshTag) ||
            (k==1 && (j==0 || j==my-2) && ((!mesh->k_periodic && !mesh->kk_periodic) || (!mesh->j_periodic && !mesh->jj_periodic)) ))
    {
        *dkdz = (K[k+1][j+1][i] + K[k+1][j][i] - K[k  ][j+1][i] - K[k  ][j][i]) * 0.5;
    }
    else
    {
        *dkdz = (K[k+1][j+1][i] + K[k+1][j][i] - K[k-1][j+1][i] - K[k-1][j][i]) * 0.25;
    }

    return;
}

//***************************************************************************************************************//

inline void Compute_dscalar_k
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    PetscReal ***K, PetscReal ***nvert, PetscReal ***meshTag,
    PetscReal *dkdc, PetscReal *dkde, PetscReal *dkdz
)
{
    if (isIBMSolidKFace(k, j, i+1, k+1, nvert) || isZeroedKFace(k, j, i+1, k+1, meshTag) || 
          (i==mx-2 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dkdc = (K[k+1][j][i  ] + K[k][j][i  ] - K[k+1][j][i-1] - K[k][j][i-1]) * 0.5;
    }
    else if (isIBMSolidKFace(k, j, i-1, k+1, nvert) || isZeroedKFace(k, j, i-1, k+1, meshTag) ||
            (i==1 && (k==0 || k==mz-2) && ((!mesh->i_periodic && !mesh->ii_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dkdc = (K[k+1][j][i+1] + K[k][j][i+1] - K[k+1][j][i  ] - K[k][j][i  ]) * 0.5;
    }
    else
    {
        *dkdc = (K[k+1][j][i+1] + K[k][j][i+1] - K[k+1][j][i-1] - K[k][j][i-1]) * 0.25;
    }

    if (isIBMSolidKFace(k, j+1, i, k+1, nvert) || isZeroedKFace(k, j+1, i, k+1, meshTag) || 
            (j==my-2 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dkde = (K[k+1][j  ][i] + K[k][j  ][i] - K[k+1][j-1][i] - K[k][j-1][i]) * 0.5;
    }
    else if (isIBMSolidKFace(k, j-1, i, k+1, nvert) || isZeroedKFace(k, j-1, i, k+1, meshTag) || 
            (j==1 && (k==0 || k==mz-2) && ((!mesh->j_periodic && !mesh->jj_periodic) || (!mesh->k_periodic && !mesh->kk_periodic)) ))
    {
        *dkde = (K[k+1][j+1][i] + K[k][j+1][i] - K[k+1][j  ][i] - K[k][j  ][i]) * 0.5;
    }
    else
    {
        *dkde = (K[k+1][j+1][i] + K[k][j+1][i] - K[k+1][j-1][i] - K[k][j-1][i]) * 0.25;
    }

    *dkdz = K[k+1][j][i] - K[k][j][i];

    return;
}

//***************************************************************************************************************//

inline void Compute_du_center
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k,  PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***ucat, PetscReal ***nvert, PetscReal ***meshTag, 
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc,
    PetscReal *dude, PetscReal *dvde, PetscReal *dwde,
    PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)
{
    if (isIBMSolidCell(k, j, i+1, nvert) || isZeroedCell(k, j, i+1, meshTag) || (!mesh->i_periodic &&  !mesh->ii_periodic && i==mx-2))
    {
        *dudc = ( ucat[k][j][i].x - ucat[k][j][i-1].x );
        *dvdc = ( ucat[k][j][i].y - ucat[k][j][i-1].y );
        *dwdc = ( ucat[k][j][i].z - ucat[k][j][i-1].z );
    }
    else if (isIBMSolidCell(k, j, i-1, nvert) || isZeroedCell(k, j, i-1, meshTag) || (!mesh->i_periodic &&  !mesh->ii_periodic && i==1))
    {
        *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i].x );
        *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i].y );
        *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i].z );
    }
    else
    {
        *dudc = ( ucat[k][j][i+1].x - ucat[k][j][i-1].x ) * 0.5;
        *dvdc = ( ucat[k][j][i+1].y - ucat[k][j][i-1].y ) * 0.5;
        *dwdc = ( ucat[k][j][i+1].z - ucat[k][j][i-1].z ) * 0.5;
    }

    if (isIBMSolidCell(k, j+1, i, nvert) || isZeroedCell(k, j+1, i, meshTag) || (!mesh->j_periodic &&  !mesh->jj_periodic && j==my-2))
    {
        *dude = ( ucat[k][j][i].x - ucat[k][j-1][i].x );
        *dvde = ( ucat[k][j][i].y - ucat[k][j-1][i].y );
        *dwde = ( ucat[k][j][i].z - ucat[k][j-1][i].z );
    }
    else if (isIBMSolidCell(k, j-1, i, nvert) || isZeroedCell(k, j-1, i, meshTag) || (!mesh->j_periodic &&  !mesh->jj_periodic && j==1))
    {
        *dude = ( ucat[k][j+1][i].x - ucat[k][j][i].x );
        *dvde = ( ucat[k][j+1][i].y - ucat[k][j][i].y );
        *dwde = ( ucat[k][j+1][i].z - ucat[k][j][i].z );

        // *dude = ( 4.0 * ucat[k][j+1][i].x - 3.0 * ucat[k][j][i].x - ucat[k][j+2][i].x) * 0.5;
        // *dvde = ( 4.0 * ucat[k][j+1][i].y - 3.0 * ucat[k][j][i].y - ucat[k][j+2][i].y) * 0.5;
        // *dwde = ( 4.0 * ucat[k][j+1][i].z - 3.0 * ucat[k][j][i].z - ucat[k][j+2][i].z) * 0.5;
    }
    else
    {
        *dude = ( ucat[k][j+1][i].x - ucat[k][j-1][i].x ) * 0.5;
        *dvde = ( ucat[k][j+1][i].y - ucat[k][j-1][i].y ) * 0.5;
        *dwde = ( ucat[k][j+1][i].z - ucat[k][j-1][i].z ) * 0.5;
    }

    if (isIBMSolidCell(k+1, j, i, nvert) || isZeroedCell(k+1, j, i, meshTag) || ( !mesh->k_periodic &&  !mesh->kk_periodic && k==mz-2))
    {
        *dudz = ( ucat[k][j][i].x - ucat[k-1][j][i].x );
        *dvdz = ( ucat[k][j][i].y - ucat[k-1][j][i].y );
        *dwdz = ( ucat[k][j][i].z - ucat[k-1][j][i].z );
    }
    else if (isIBMSolidCell(k-1, j, i, nvert) || isZeroedCell(k-1, j, i, meshTag) || (!mesh->k_periodic &&  !mesh->kk_periodic && k==1))
    {
        *dudz = ( ucat[k+1][j][i].x - ucat[k][j][i].x );
        *dvdz = ( ucat[k+1][j][i].y - ucat[k][j][i].y );
        *dwdz = ( ucat[k+1][j][i].z - ucat[k][j][i].z );
    }
    else
    {
        *dudz = ( ucat[k+1][j][i].x - ucat[k-1][j][i].x ) * 0.5;
        *dvdz = ( ucat[k+1][j][i].y - ucat[k-1][j][i].y ) * 0.5;
        *dwdz = ( ucat[k+1][j][i].z - ucat[k-1][j][i].z ) * 0.5;
    }

    return;
}

//***************************************************************************************************************//

inline void Compute_dscalar_center
(
    mesh_ *mesh,
    PetscInt i, PetscInt j, PetscInt k,  PetscInt mx, PetscInt my, PetscInt mz,
    PetscReal ***K, PetscReal ***nvert, PetscReal ***meshTag, 
    PetscReal *dkdc, PetscReal *dkde, PetscReal *dkdz
)
{

	if      (i==mx-1) *dkdc = ( K[k][j][i] - K[k][j][i-1] );
	else if (i==0)    *dkdc = ( K[k][j][i+1] - K[k][j][i] );
	else if (isIBMSolidCell(k, j, i+1, nvert) || isZeroedCell(k, j, i+1, meshTag))
    {
		*dkdc = ( K[k][j][i] - K[k][j][i-1] );
	}
	else if (isIBMSolidCell(k, j, i-1, nvert) || isZeroedCell(k, j, i-1, meshTag))
    {
		*dkdc = ( K[k][j][i+1] - K[k][j][i] );
	}
	else
    {
		*dkdc = ( K[k][j][i+1] - K[k][j][i-1] ) * 0.5;
	}

	if      (j==my-1) *dkde = ( K[k][j][i] - K[k][j-1][i] );
	else if (j==0)    *dkde = ( K[k][j+1][i] - K[k][j][i] );
	else if (isIBMSolidCell(k, j+1, i, nvert) || isZeroedCell(k, j+1, i, meshTag))
    {
		*dkde = ( K[k][j][i] - K[k][j-1][i] );
	}
	else if (isIBMSolidCell(k, j-1, i, nvert) || isZeroedCell(k, j-1, i, meshTag))
    {
		*dkde = ( K[k][j+1][i] - K[k][j][i] );
	}
	else
    {
		*dkde = ( K[k][j+1][i] - K[k][j-1][i] ) * 0.5;
	}

	if      (k==mz-1) *dkdz = ( K[k][j][i] - K[k-1][j][i] );
	else if (k==0)    *dkdz = ( K[k+1][j][i] - K[k][j][i] );
	else if (isIBMSolidCell(k+1, j, i, nvert) || isZeroedCell(k+1, j, i, meshTag) )
    {
		*dkdz = ( K[k][j][i] - K[k-1][j][i] );
	}
	else if (isIBMSolidCell(k-1, j, i, nvert) || isZeroedCell(k-1, j, i, meshTag))
    {
		*dkdz = ( K[k+1][j][i] - K[k][j][i] );
	}
	else
    {
		*dkdz = ( K[k+1][j][i] - K[k-1][j][i] ) * 0.5;
	}

    return;
}

//***************************************************************************************************************//

// Simpson integration algorithm

inline PetscReal integrateTestfilterSimpson(PetscReal val[3][3][3], PetscReal w[3][3][3])
{

    PetscReal wsum=0, valsum=0;

    for(PetscInt i=0; i<3; i++)
    for(PetscInt j=0; j<3; j++)
    for(PetscInt k=0; k<3; k++)
    {
        PetscReal simpson_w = 1.0;
        if(i==1) simpson_w *= 4.;
        if(j==1) simpson_w *= 4.;
        if(k==1) simpson_w *= 4.;

        wsum   += simpson_w * w[i][j][k];
        valsum += simpson_w * w[i][j][k] * val[i][j][k];
    }

    return valsum / wsum;
};

inline PetscReal integrateTestfilterSimpson5x5(PetscReal val[5][5][5], PetscReal w[5][5][5])
{
    PetscReal wsum=0, valsum=0;

    for(PetscInt i=0; i<5; i++)
    for(PetscInt j=0; j<5; j++)
    for(PetscInt k=0; k<5; k++)
    {
        PetscReal simpson_w = 1.0;
        simpson_w *= (i == 0 || i == 4) ? 1 : (i % 2 == 0 ? 2 : 4);
        simpson_w *= (j == 0 || j == 4) ? 1 : (j % 2 == 0 ? 2 : 4);
        simpson_w *= (k == 0 || k == 4) ? 1 : (k % 2 == 0 ? 2 : 4);

        wsum   += simpson_w * w[i][j][k];
        valsum += simpson_w * w[i][j][k] * val[i][j][k];
    }

    return valsum / wsum;
};

//***************************************************************************************************************//

inline void Compute_du_dxyz
(
    mesh_ *mesh,
    PetscReal csi0, PetscReal csi1, PetscReal csi2,
    PetscReal eta0, PetscReal eta1, PetscReal eta2,
    PetscReal zet0, PetscReal zet1, PetscReal zet2,
    PetscReal ajc,
	PetscReal dudc, PetscReal dvdc, PetscReal dwdc,
    PetscReal dude, PetscReal dvde, PetscReal dwde,
    PetscReal dudz, PetscReal dvdz, PetscReal dwdz,
	PetscReal *du_dx, PetscReal *dv_dx, PetscReal *dw_dx,
    PetscReal *du_dy, PetscReal *dv_dy, PetscReal *dw_dy,
    PetscReal *du_dz, PetscReal *dv_dz, PetscReal *dw_dz
)
{
	*du_dx = (dudc * csi0 + dude * eta0 + dudz * zet0) * ajc;
	*du_dy = (dudc * csi1 + dude * eta1 + dudz * zet1) * ajc;
	*du_dz = (dudc * csi2 + dude * eta2 + dudz * zet2) * ajc;
	*dv_dx = (dvdc * csi0 + dvde * eta0 + dvdz * zet0) * ajc;
	*dv_dy = (dvdc * csi1 + dvde * eta1 + dvdz * zet1) * ajc;
	*dv_dz = (dvdc * csi2 + dvde * eta2 + dvdz * zet2) * ajc;
	*dw_dx = (dwdc * csi0 + dwde * eta0 + dwdz * zet0) * ajc;
	*dw_dy = (dwdc * csi1 + dwde * eta1 + dwdz * zet1) * ajc;
	*dw_dz = (dwdc * csi2 + dwde * eta2 + dwdz * zet2) * ajc;

    return;
}

//***************************************************************************************************************//

inline void Compute_dscalar_dxyz
(
    mesh_ *mesh,
    PetscReal csi0, PetscReal csi1, PetscReal csi2,
    PetscReal eta0, PetscReal eta1, PetscReal eta2,
    PetscReal zet0, PetscReal zet1, PetscReal zet2,
    PetscReal ajc,
	PetscReal dkdc, PetscReal dkde, PetscReal dkdz,
    PetscReal *dk_dx, PetscReal *dk_dy, PetscReal *dk_dz
)
{
	*dk_dx = (dkdc * csi0 + dkde * eta0 + dkdz * zet0) * ajc;
	*dk_dy = (dkdc * csi1 + dkde * eta1 + dkdz * zet1) * ajc;
	*dk_dz = (dkdc * csi2 + dkde * eta2 + dkdz * zet2) * ajc;

    return;
}

//***************************************************************************************************************//

inline void Compute_du_wmLocal
(
    mesh_ *mesh,
    Cmpnts eN, Cmpnts eT1, Cmpnts eT2,
    PetscReal du_dx, PetscReal dv_dx, PetscReal dw_dx,
    PetscReal du_dy, PetscReal dv_dy, PetscReal dw_dy,
    PetscReal du_dz, PetscReal dv_dz, PetscReal dw_dz,
    PetscReal *dut1dn, PetscReal *dut2dn, PetscReal *dundn,
    PetscReal *dut1dt1, PetscReal *dut2dt1, PetscReal *dundt1,
    PetscReal *dut1dt2, PetscReal *dut2dt2, PetscReal *dundt2
)
{
    PetscReal dudn = du_dx * eN.x + du_dy * eN.y + du_dz * eN.z;
    PetscReal dvdn = dv_dx * eN.x + dv_dy * eN.y + dv_dz * eN.z;
    PetscReal dwdn = dw_dx * eN.x + dw_dy * eN.y + dw_dz * eN.z;

    PetscReal dudt1 = du_dx * eT1.x + du_dy * eT1.y + du_dz * eT1.z;
    PetscReal dvdt1 = dv_dx * eT1.x + dv_dy * eT1.y + dv_dz * eT1.z;
    PetscReal dwdt1 = dw_dx * eT1.x + dw_dy * eT1.y + dw_dz * eT1.z;

    PetscReal dudt2 = du_dx * eT2.x + du_dy * eT2.y + du_dz * eT2.z;
    PetscReal dvdt2 = dv_dx * eT2.x + dv_dy * eT2.y + dv_dz * eT2.z;
    PetscReal dwdt2 = dw_dx * eT2.x + dw_dy * eT2.y + dw_dz * eT2.z;

    *dut1dn = dudn * eT1.x + dvdn * eT1.y + dwdn * eT1.z;
    *dut2dn = dudn * eT2.x + dvdn * eT2.y + dwdn * eT2.z;
    *dundn = dudn * eN.x + dvdn * eN.y + dwdn * eN.z;

    *dut1dt1 = dudt1 * eT1.x + dvdt1 * eT1.y + dwdt1 * eT1.z;
    *dut2dt1 = dudt1 * eT2.x + dvdt1 * eT2.y + dwdt1 * eT2.z;
    *dundt1 = dudt1 * eN.x + dvdt1 * eN.y + dwdt1 * eN.z;

    *dut1dt2 = dudt2 * eT1.x + dvdt2 * eT1.y + dwdt2 * eT1.z;
    *dut2dt2 = dudt2 * eT2.x + dvdt2 * eT2.y + dwdt2 * eT2.z;
    *dundt2 = dudt2 * eN.x + dvdt2 * eN.y + dwdt2 * eN.z;

    return;
}

//***************************************************************************************************************//

inline void Comput_JacobTensor_i
(
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***coor, PetscReal *dxdc, PetscReal *dxde, PetscReal *dxdz,
    PetscReal *dydc, PetscReal *dyde, PetscReal *dydz, PetscReal *dzdc, PetscReal *dzde, PetscReal *dzdz
)
{

    PetscReal centx, centy, centz;
    PetscReal centx_ip1, centy_ip1, centz_ip1;
    PetscReal centx_im1, centy_im1, centz_im1;
    PetscReal centx_jp1, centy_jp1, centz_jp1;
    PetscReal centx_jm1, centy_jm1, centz_jm1;
    PetscReal centx_kp1, centy_kp1, centz_kp1;
    PetscReal centx_km1, centy_km1, centz_km1;


	PetscInt i1=i,j1=j,k1=k;

    centx = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
                 coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
    centy = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
                 coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
    centz = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
                 coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;

	if (i!=mx-2)
    {
        i1=i+1,j1=j,k1=k;

    	centx_ip1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
            	     coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_ip1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
    	             coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
        centz_ip1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
    	             coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (i!=0)
    {
        i1=i-1,j1=j,k1=k;

        centx_im1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
    	             coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_im1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
    	             coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
        centz_im1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
    	             coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (j!=my-2)
    {
    	i1=i,j1=j+1,k1=k;

        centx_jp1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
    	             coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_jp1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
    	             coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
        centz_jp1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
    	             coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (j!=1)
    {
        i1=i,j1=j-1,k1=k;

        centx_jm1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
    	             coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_jm1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
            	     coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
    	centz_jm1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
                     coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (k!=mz-2)
    {
    	i1=i,j1=j,k1=k+1;

    	centx_kp1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
            	     coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_kp1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
    	             coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
        centz_kp1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
    	             coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (k!=1)
    {
    	i1=i,j1=j,k1=k-1;

        centx_km1 = (coor[k1  ][j1  ][i1].x + coor[k1-1][j1  ][i1].x +
    	             coor[k1  ][j1-1][i1].x + coor[k1-1][j1-1][i1].x) * 0.25;
        centy_km1 = (coor[k1  ][j1  ][i1].y + coor[k1-1][j1  ][i1].y +
            	     coor[k1  ][j1-1][i1].y + coor[k1-1][j1-1][i1].y) * 0.25;
        centz_km1 = (coor[k1  ][j1  ][i1].z + coor[k1-1][j1  ][i1].z +
                     coor[k1  ][j1-1][i1].z + coor[k1-1][j1-1][i1].z) * 0.25;
	}

	if (i==0)
    {
	  *dxdc = centx_ip1 - centx;
	  *dydc = centy_ip1 - centy;
	  *dzdc = centz_ip1 - centz;
	}
	else if (i==mx-2)
    {
	  *dxdc = centx - centx_im1;
	  *dydc = centy - centy_im1;
	  *dzdc = centz - centz_im1;
	}
	else
    {
	  *dxdc = (centx_ip1 - centx_im1) * 0.5;
	  *dydc = (centy_ip1 - centy_im1) * 0.5;
	  *dzdc = (centz_ip1 - centz_im1) * 0.5;
	}


	if (j==1)
    {
	  *dxde = centx_jp1 - centx;
	  *dyde = centy_jp1 - centy;
	  *dzde = centz_jp1 - centz;
	}
	else if (j==my-2)
    {
	  *dxde = centx - centx_jm1;
	  *dyde = centy - centy_jm1;
	  *dzde = centz - centz_jm1;
	}
	else
    {
	  *dxde = (centx_jp1 - centx_jm1) * 0.5;
	  *dyde = (centy_jp1 - centy_jm1) * 0.5;
	  *dzde = (centz_jp1 - centz_jm1) * 0.5;
	}

	if (k==1)
    {
	  *dxdz = (centx_kp1 - centx);
	  *dydz = (centy_kp1 - centy);
	  *dzdz = (centz_kp1 - centz);
	}
	else if (k==mz-2)
    {
	  *dxdz = (centx - centx_km1);
	  *dydz = (centy - centy_km1);
	  *dzdz = (centz - centz_km1);
	}
	else
    {
	  *dxdz = (centx_kp1 - centx_km1) * 0.5;
	  *dydz = (centy_kp1 - centy_km1) * 0.5;
	  *dzdz = (centz_kp1 - centz_km1) * 0.5;
	}

    return;
}

//***************************************************************************************************************//

inline void Comput_JacobTensor_j
(
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***coor, PetscReal *dxdc, PetscReal *dxde, PetscReal *dxdz,
    PetscReal *dydc, PetscReal *dyde, PetscReal *dydz, PetscReal *dzdc, PetscReal *dzde, PetscReal *dzdz
)
{
	PetscReal centx, centy, centz;
	PetscReal centx_ip1, centy_ip1, centz_ip1;
	PetscReal centx_im1, centy_im1, centz_im1;
	PetscReal centx_jp1, centy_jp1, centz_jp1;
	PetscReal centx_jm1, centy_jm1, centz_jm1;
	PetscReal centx_kp1, centy_kp1, centz_kp1;
	PetscReal centx_km1, centy_km1, centz_km1;

	PetscInt i1=i,j1=j,k1=k;
    centx = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
             coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
    centy = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
             coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
    centz = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
             coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	if (i!=mx-2)
    {
	        i1=i+1,j1=j,k1=k;

        	centx_ip1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                        coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
	        centy_ip1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                        coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
	        centz_ip1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                        coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (i!=1)
    {
        i1=i-1,j1=j,k1=k;

        centx_im1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
        centy_im1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
        centz_im1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (j!=my-2)
    {
        i1=i,j1=j+1,k1=k;

        centx_jp1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
        centy_jp1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
        centz_jp1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (j!=0)
    {
        i1=i,j1=j-1,k1=k;

        centx_jm1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
        centy_jm1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
        centz_jm1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (k!=mz-2)
    {
        i1=i,j1=j,k1=k+1;

        centx_kp1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
        centy_kp1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
        centz_kp1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (k!=1)
    {
        i1=i,j1=j,k1=k-1;

        centx_km1 = (coor[k1  ][j1][i1  ].x + coor[k1-1][j1][i1  ].x +
                coor[k1  ][j1][i1-1].x + coor[k1-1][j1][i1-1].x) * 0.25;
        centy_km1 = (coor[k1  ][j1][i1  ].y + coor[k1-1][j1][i1  ].y +
                coor[k1  ][j1][i1-1].y + coor[k1-1][j1][i1-1].y) * 0.25;
        centz_km1 = (coor[k1  ][j1][i1  ].z + coor[k1-1][j1][i1  ].z +
                coor[k1  ][j1][i1-1].z + coor[k1-1][j1][i1-1].z) * 0.25;

	}

	if (i==1)
    {
	  *dxdc = centx_ip1 - centx;
	  *dydc = centy_ip1 - centy;
	  *dzdc = centz_ip1 - centz;
	}
	else if (i==mx-2)
    {
	  *dxdc = centx - centx_im1;
	  *dydc = centy - centy_im1;
	  *dzdc = centz - centz_im1;
	}
	else {
	  *dxdc = (centx_ip1 - centx_im1) * 0.5;
	  *dydc = (centy_ip1 - centy_im1) * 0.5;
	  *dzdc = (centz_ip1 - centz_im1) * 0.5;
	}


	if (j==0)
    {
	  *dxde = centx_jp1 - centx;
	  *dyde = centy_jp1 - centy;
	  *dzde = centz_jp1 - centz;
	}
	else if (j==my-2)
    {
	  *dxde = centx - centx_jm1;
	  *dyde = centy - centy_jm1;
	  *dzde = centz - centz_jm1;
	}
	else
    {
	  *dxde = (centx_jp1 - centx_jm1) * 0.5;
	  *dyde = (centy_jp1 - centy_jm1) * 0.5;
	  *dzde = (centz_jp1 - centz_jm1) * 0.5;
	}

	if (k==1)
    {
	  *dxdz = (centx_kp1 - centx);
	  *dydz = (centy_kp1 - centy);
	  *dzdz = (centz_kp1 - centz);
	}
	else if (k==mz-2)
    {
	  *dxdz = (centx - centx_km1);
	  *dydz = (centy - centy_km1);
	  *dzdz = (centz - centz_km1);
	}
	else
    {
	  *dxdz = (centx_kp1 - centx_km1) * 0.5;
	  *dydz = (centy_kp1 - centy_km1) * 0.5;
	  *dzdz = (centz_kp1 - centz_km1) * 0.5;
	}

}

//***************************************************************************************************************//

inline void Comput_JacobTensor_k
(
    PetscInt i, PetscInt j, PetscInt k, PetscInt mx, PetscInt my, PetscInt mz,
    Cmpnts ***coor, PetscReal *dxdc, PetscReal *dxde, PetscReal *dxdz,
    PetscReal *dydc, PetscReal *dyde, PetscReal *dydz, PetscReal *dzdc, PetscReal *dzde, PetscReal *dzdz
)
{
	PetscReal centx, centy, centz;
	PetscReal centx_ip1, centy_ip1, centz_ip1;
	PetscReal centx_im1, centy_im1, centz_im1;
	PetscReal centx_jp1, centy_jp1, centz_jp1;
	PetscReal centx_jm1, centy_jm1, centz_jm1;
	PetscReal centx_kp1, centy_kp1, centz_kp1;
	PetscReal centx_km1, centy_km1, centz_km1;

	PetscInt i1=i,j1=j,k1=k;

    centx = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
             coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
    centy = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
             coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
    centz = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
             coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	if (i!=mx-2)
    {
        i1=i+1,j1=j,k1=k;

        centx_ip1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
             	     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
    	centy_ip1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
             	     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
    	centz_ip1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
             	     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (i!=1)
    {
        i1=i-1,j1=j,k1=k;

        centx_im1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
        centy_im1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
                     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
        centz_im1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
                     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (j!=my-2)
    {
        i1=i,j1=j+1,k1=k;

        centx_jp1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
        centy_jp1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
                     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
        centz_jp1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
                     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (j!=1)
    {
        i1=i,j1=j-1,k1=k;

        centx_jm1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
        centy_jm1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
                     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
        centz_jm1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
                     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (k!=mz-2)
    {
        i1=i,j1=j,k1=k+1;

        centx_kp1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
        centy_kp1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
                     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
        centz_kp1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
                     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (k!=0)
    {
        i1=i,j1=j,k1=k-1;

        centx_km1 = (coor[k1  ][j1][i1  ].x + coor[k1][j1-1][i1  ].x +
                     coor[k1  ][j1][i1-1].x + coor[k1][j1-1][i1-1].x) * 0.25;
        centy_km1 = (coor[k1  ][j1][i1  ].y + coor[k1][j1-1][i1  ].y +
                     coor[k1  ][j1][i1-1].y + coor[k1][j1-1][i1-1].y) * 0.25;
        centz_km1 = (coor[k1  ][j1][i1  ].z + coor[k1][j1-1][i1  ].z +
                     coor[k1  ][j1][i1-1].z + coor[k1][j1-1][i1-1].z) * 0.25;

	}

	if (i==1)
    {
	  *dxdc = centx_ip1 - centx;
	  *dydc = centy_ip1 - centy;
	  *dzdc = centz_ip1 - centz;
	}
	else if (i==mx-2)
    {
	  *dxdc = centx - centx_im1;
	  *dydc = centy - centy_im1;
	  *dzdc = centz - centz_im1;
	}
	else
    {
	  *dxdc = (centx_ip1 - centx_im1) * 0.5;
	  *dydc = (centy_ip1 - centy_im1) * 0.5;
	  *dzdc = (centz_ip1 - centz_im1) * 0.5;
	}


	if (j==1)
    {
	  *dxde = centx_jp1 - centx;
	  *dyde = centy_jp1 - centy;
	  *dzde = centz_jp1 - centz;
	}
	else if (j==my-2)
    {
	  *dxde = centx - centx_jm1;
	  *dyde = centy - centy_jm1;
	  *dzde = centz - centz_jm1;
	}
	else
    {
	  *dxde = (centx_jp1 - centx_jm1) * 0.5;
	  *dyde = (centy_jp1 - centy_jm1) * 0.5;
	  *dzde = (centz_jp1 - centz_jm1) * 0.5;
	}

	if (k==0)
    {
	  *dxdz = (centx_kp1 - centx);
	  *dydz = (centy_kp1 - centy);
	  *dzdz = (centz_kp1 - centz);
	}
	else if (k==mz-2)
    {
	  *dxdz = (centx - centx_km1);
	  *dydz = (centy - centy_km1);
	  *dzdz = (centz - centz_km1);
	}
	else
    {
	  *dxdz = (centx_kp1 - centx_km1) * 0.5;
	  *dydz = (centy_kp1 - centy_km1) * 0.5;
	  *dzdz = (centz_kp1 - centz_km1) * 0.5;
	}

}

//***************************************************************************************************************//

// From local wall model grid to computational grid
inline void Compute_du_Compgrid
(
    PetscReal dxdc, PetscReal dxde, PetscReal dxdz, PetscReal dydc, PetscReal dyde, PetscReal dydz, PetscReal dzdc, PetscReal dzde, PetscReal dzdz,
    PetscReal nx, PetscReal ny, PetscReal nz, PetscReal t1x, PetscReal t1y, PetscReal t1z, PetscReal t2x, PetscReal t2y, PetscReal t2z,
    PetscReal dut1dn, PetscReal dut2dn, PetscReal dundn, PetscReal dut1dt1, PetscReal dut2dt1, PetscReal dundt1, PetscReal dut1dt2, PetscReal dut2dt2, PetscReal dundt2,
    PetscReal *dudc, PetscReal *dvdc, PetscReal *dwdc, PetscReal *dude, PetscReal *dvde, PetscReal *dwde, PetscReal *dudz, PetscReal *dvdz, PetscReal *dwdz
)
{

	PetscReal dxdn=nx, dydn=ny, dzdn=nz;
	PetscReal dxdt1=t1x, dydt1=t1y, dzdt1=t1z;
	PetscReal dxdt2=t2x, dydt2=t2y, dzdt2=t2z;

	PetscReal dndx = dydt1*dzdt2-dydt2*dzdt1;
	PetscReal dt1dx = dydt2*dzdn-dydn*dzdt2;
	PetscReal dt2dx = dydn*dzdt1-dydt1*dzdn;

    PetscReal dndy = dzdt1*dxdt2-dzdt2*dxdt1;
    PetscReal dt1dy = dzdt2*dxdn-dzdn*dxdt2;
    PetscReal dt2dy = dzdn*dxdt1-dzdt1*dxdn;

    PetscReal dndz = dxdt1*dydt2-dxdt2*dydt1;
    PetscReal dt1dz = dxdt2*dydn-dxdn*dydt2;
    PetscReal dt2dz = dxdn*dydt1-dxdt1*dydn;


	PetscReal dundx = dundn*dndx+dundt1*dt1dx+dundt2*dt2dx;
	PetscReal dundy = dundn*dndy+dundt1*dt1dy+dundt2*dt2dy;
	PetscReal dundz = dundn*dndz+dundt1*dt1dz+dundt2*dt2dz;

    PetscReal dut1dx = dut1dn*dndx+dut1dt1*dt1dx+dut1dt2*dt2dx;
    PetscReal dut1dy = dut1dn*dndy+dut1dt1*dt1dy+dut1dt2*dt2dy;
    PetscReal dut1dz = dut1dn*dndz+dut1dt1*dt1dz+dut1dt2*dt2dz;

    PetscReal dut2dx = dut2dn*dndx+dut2dt1*dt1dx+dut2dt2*dt2dx;
    PetscReal dut2dy = dut2dn*dndy+dut2dt1*dt1dy+dut2dt2*dt2dy;
    PetscReal dut2dz = dut2dn*dndz+dut2dt1*dt1dz+dut2dt2*dt2dz;


	PetscReal du_dx = dundx*nx+dut1dx*t1x+dut2dx*t2x;
	PetscReal du_dy = dundy*nx+dut1dy*t1x+dut2dy*t2x;
	PetscReal du_dz = dundz*nx+dut1dz*t1x+dut2dz*t2x;

    PetscReal dv_dx = dundx*ny+dut1dx*t1y+dut2dx*t2y;
    PetscReal dv_dy = dundy*ny+dut1dy*t1y+dut2dy*t2y;
    PetscReal dv_dz = dundz*ny+dut1dz*t1y+dut2dz*t2y;

    PetscReal dw_dx = dundx*nz+dut1dx*t1z+dut2dx*t2z;
    PetscReal dw_dy = dundy*nz+dut1dy*t1z+dut2dy*t2z;
    PetscReal dw_dz = dundz*nz+dut1dz*t1z+dut2dz*t2z;

	*dudc = du_dx*dxdc+du_dy*dydc+du_dz*dzdc;
	*dude = du_dx*dxde+du_dy*dyde+du_dz*dzde;
	*dudz = du_dx*dxdz+du_dy*dydz+du_dz*dzdz;

    *dvdc = dv_dx*dxdc+dv_dy*dydc+dv_dz*dzdc;
    *dvde = dv_dx*dxde+dv_dy*dyde+dv_dz*dzde;
    *dvdz = dv_dx*dxdz+dv_dy*dydz+dv_dz*dzdz;

    *dwdc = dw_dx*dxdc+dw_dy*dydc+dw_dz*dzdc;
    *dwde = dw_dx*dxde+dw_dy*dyde+dw_dz*dzde;
    *dwdz = dw_dx*dxdz+dw_dy*dydz+dw_dz*dzdz;
}

// DIVERGENCE SCHEMES
// ============================================================================================================= //

inline PetscReal minMod(PetscReal m1, PetscReal m2)
{
	return 0.5 * ( sign(m1)+sign(m2) ) * PetscMin ( fabs(m1), fabs(m2) );
}

//***************************************************************************************************************//

inline PetscReal vanLeer(PetscReal f0, PetscReal f1, PetscReal f2)
{
    PetscReal r, n, d;

    if(f1>0)
    {
        n = f1 - f0; d = PetscMax(f2 - f1, 1e-5); r = n/d;
	}
	else
    {
        n = f1 - f2; d = PetscMax(f0 - f1, 1e-5); r = n/d;
	}

    // Van Leer limiter function
    return 0.5 * (r + fabs(r)) / (1.0 + fabs(r));
}

//***************************************************************************************************************//

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
        }
        else if(mesh->ii_periodic && i==mx-2)
        {
            im1 = i-1, i0  = i, ip1  = i+1, ip2 = mx+2;
            
            result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
            result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
            result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;
        }
        else if(mesh->ii_periodic && i==0)
        {
            im1 = -3, i0  = i, ip1  = i+1, ip2 = i+2;

            result.x = (-ucat[k][j][im1].x + 7.0 * ucat[k][j][i0].x + 7.0 * ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 12.0;
            result.y = (-ucat[k][j][im1].y + 7.0 * ucat[k][j][i0].y + 7.0 * ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 12.0;
            result.z = (-ucat[k][j][im1].z + 7.0 * ucat[k][j][i0].z + 7.0 * ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 12.0;
        }
        else
        {
            if(i == 0)
            {
                im1 = i, i0 = i+1, ip1 = i+2, ip2 = i+3;
                
                result.x = (5.0 * ucat[k][j][im1].x + 5.0 * ucat[k][j][i0].x - ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 8.0;
                result.y = (5.0 * ucat[k][j][im1].y + 5.0 * ucat[k][j][i0].y - ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 8.0;
                result.z = (5.0 * ucat[k][j][im1].z + 5.0 * ucat[k][j][i0].z - ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 8.0;
            }

            if(i == mx-2)
            {
                im1 = i-2, i0 = i-1, ip1 = i, ip2 = i+1;

                result.x = (-ucat[k][j][im1].x - ucat[k][j][i0].x + 5.0 * ucat[k][j][ip1].x + 5.0 * ucat[k][j][ip2].x) / 8.0;
                result.y = (-ucat[k][j][im1].y - ucat[k][j][i0].y + 5.0 * ucat[k][j][ip1].y + 5.0 * ucat[k][j][ip2].y) / 8.0;
                result.z = (-ucat[k][j][im1].z - ucat[k][j][i0].z + 5.0 * ucat[k][j][ip1].z + 5.0 * ucat[k][j][ip2].z) / 8.0;
            }
        }
    }
    else if(isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        im1 = i, i0 = i+1, ip1 = i+2, ip2 = i+3;
        result.x = (5.0 * ucat[k][j][im1].x + 5.0 * ucat[k][j][i0].x - ucat[k][j][ip1].x - ucat[k][j][ip2].x) / 8.0;
        result.y = (5.0 * ucat[k][j][im1].y + 5.0 * ucat[k][j][i0].y - ucat[k][j][ip1].y - ucat[k][j][ip2].y) / 8.0;
        result.z = (5.0 * ucat[k][j][im1].z + 5.0 * ucat[k][j][i0].z - ucat[k][j][ip1].z - ucat[k][j][ip2].z) / 8.0;
         
    }
    else if(isIBMCell(k, j, i+1, nvert) || isOversetCell(k, j, i+1, meshTag))
    {
        im1 = i-2, i0 = i-1, ip1 = i, ip2 = i+1;

        result.x = (-ucat[k][j][im1].x - ucat[k][j][i0].x + 5.0 * ucat[k][j][ip1].x + 5.0 * ucat[k][j][ip2].x) / 8.0;
        result.y = (-ucat[k][j][im1].y - ucat[k][j][i0].y + 5.0 * ucat[k][j][ip1].y + 5.0 * ucat[k][j][ip2].y) / 8.0;
        result.z = (-ucat[k][j][im1].z - ucat[k][j][i0].z + 5.0 * ucat[k][j][ip1].z + 5.0 * ucat[k][j][ip2].z) / 8.0;
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
        }
        else if (mesh->jj_periodic && j == my-2) // Assuming jj_periodic for y-direction
        {
            jm1 = j-1, j0 = j, jp1 = j+1, jp2 = my+2;

            result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
            result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
            result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;
        }
        else if (mesh->jj_periodic && j == 0)
        {
            jm1 = -3, j0 = j, jp1 = j+1, jp2 = j+2;

            result.x = (-ucat[k][jm1][i].x + 7.0 * ucat[k][j0][i].x + 7.0 * ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 12.0;
            result.y = (-ucat[k][jm1][i].y + 7.0 * ucat[k][j0][i].y + 7.0 * ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 12.0;
            result.z = (-ucat[k][jm1][i].z + 7.0 * ucat[k][j0][i].z + 7.0 * ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 12.0;
        }
        else
        {
            if (j == 0) // Forward stencil at lower boundary
            {
                jm1 = j, j0 = j+1, jp1 = j+2, jp2 = j+3;

                result.x = (5.0 * ucat[k][jm1][i].x + 5.0 * ucat[k][j0][i].x - ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 8.0;
                result.y = (5.0 * ucat[k][jm1][i].y + 5.0 * ucat[k][j0][i].y - ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 8.0;
                result.z = (5.0 * ucat[k][jm1][i].z + 5.0 * ucat[k][j0][i].z - ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 8.0;
            }
            if (j == my-2) // Backward stencil at upper boundary
            {
                jm1 = j-2, j0 = j-1, jp1 = j, jp2 = j+1;

                result.x = (-ucat[k][jm1][i].x - ucat[k][j0][i].x + 5.0 * ucat[k][jp1][i].x + 5.0 * ucat[k][jp2][i].x) / 8.0;
                result.y = (-ucat[k][jm1][i].y - ucat[k][j0][i].y + 5.0 * ucat[k][jp1][i].y + 5.0 * ucat[k][jp2][i].y) / 8.0;
                result.z = (-ucat[k][jm1][i].z - ucat[k][j0][i].z + 5.0 * ucat[k][jp1][i].z + 5.0 * ucat[k][jp2][i].z) / 8.0;
            }
        }
    }
    else if (isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        jm1 = j, j0 = j+1, jp1 = j+2, jp2 = j+3;

        result.x = (5.0 * ucat[k][jm1][i].x + 5.0 * ucat[k][j0][i].x - ucat[k][jp1][i].x - ucat[k][jp2][i].x) / 8.0;
        result.y = (5.0 * ucat[k][jm1][i].y + 5.0 * ucat[k][j0][i].y - ucat[k][jp1][i].y - ucat[k][jp2][i].y) / 8.0;
        result.z = (5.0 * ucat[k][jm1][i].z + 5.0 * ucat[k][j0][i].z - ucat[k][jp1][i].z - ucat[k][jp2][i].z) / 8.0;
    }
    else if (isIBMCell(k, j+1, i, nvert) || isOversetCell(k, j+1, i, meshTag))
    {
        jm1 = j-2, j0 = j-1, jp1 = j, jp2 = j+1;

        result.x = (-ucat[k][jm1][i].x - ucat[k][j0][i].x + 5.0 * ucat[k][jp1][i].x + 5.0 * ucat[k][jp2][i].x) / 8.0;
        result.y = (-ucat[k][jm1][i].y - ucat[k][j0][i].y + 5.0 * ucat[k][jp1][i].y + 5.0 * ucat[k][jp2][i].y) / 8.0;
        result.z = (-ucat[k][jm1][i].z - ucat[k][j0][i].z + 5.0 * ucat[k][jp1][i].z + 5.0 * ucat[k][jp2][i].z) / 8.0;
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
        }
        else if (mesh->kk_periodic && k == mz-2) // Assuming kk_periodic for z-direction
        {
            km1 = k-1, k0 = k, kp1 = k+1, kp2 = mz+2;

            result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
            result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
            result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;
        }
        else if (mesh->kk_periodic && k == 0)
        {
            km1 = -3, k0 = k, kp1 = k+1, kp2 = k+2;

            result.x = (-ucat[km1][j][i].x + 7.0 * ucat[k0][j][i].x + 7.0 * ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 12.0;
            result.y = (-ucat[km1][j][i].y + 7.0 * ucat[k0][j][i].y + 7.0 * ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 12.0;
            result.z = (-ucat[km1][j][i].z + 7.0 * ucat[k0][j][i].z + 7.0 * ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 12.0;
        }
        else
        {
            if (k == 0) // Forward stencil at lower boundary
            {
                km1 = k, k0 = k+1, kp1 = k+2, kp2 = k+3;

                result.x = (5.0 * ucat[km1][j][i].x + 5.0 * ucat[k0][j][i].x - ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 8.0;
                result.y = (5.0 * ucat[km1][j][i].y + 5.0 * ucat[k0][j][i].y - ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 8.0;
                result.z = (5.0 * ucat[km1][j][i].z + 5.0 * ucat[k0][j][i].z - ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 8.0;
            }
            if (k == mz-2) // Backward stencil at upper boundary
            {
                km1 = k-2, k0 = k-1, kp1 = k, kp2 = k+1;

                result.x = (-ucat[km1][j][i].x - ucat[k0][j][i].x + 5.0 * ucat[kp1][j][i].x + 5.0 * ucat[kp2][j][i].x) / 8.0;
                result.y = (-ucat[km1][j][i].y - ucat[k0][j][i].y + 5.0 * ucat[kp1][j][i].y + 5.0 * ucat[kp2][j][i].y) / 8.0;
                result.z = (-ucat[km1][j][i].z - ucat[k0][j][i].z + 5.0 * ucat[kp1][j][i].z + 5.0 * ucat[kp2][j][i].z) / 8.0;
            }
        }
    }
    else if (isIBMCell(k, j, i, nvert) || isOversetCell(k, j, i, meshTag))
    {
        km1 = k, k0 = k+1, kp1 = k+2, kp2 = k+3;

        result.x = (5.0 * ucat[km1][j][i].x + 5.0 * ucat[k0][j][i].x - ucat[kp1][j][i].x - ucat[kp2][j][i].x) / 8.0;
        result.y = (5.0 * ucat[km1][j][i].y + 5.0 * ucat[k0][j][i].y - ucat[kp1][j][i].y - ucat[kp2][j][i].y) / 8.0;
        result.z = (5.0 * ucat[km1][j][i].z + 5.0 * ucat[k0][j][i].z - ucat[kp1][j][i].z - ucat[kp2][j][i].z) / 8.0;
    }
    else if (isIBMCell(k+1, j, i, nvert) || isOversetCell(k+1, j, i, meshTag))
    {
        km1 = k-2, k0 = k-1, kp1 = k, kp2 = k+1;

        result.x = (-ucat[km1][j][i].x - ucat[k0][j][i].x + 5.0 * ucat[kp1][j][i].x + 5.0 * ucat[kp2][j][i].x) / 8.0;
        result.y = (-ucat[km1][j][i].y - ucat[k0][j][i].y + 5.0 * ucat[kp1][j][i].y + 5.0 * ucat[kp2][j][i].y) / 8.0;
        result.z = (-ucat[km1][j][i].z - ucat[k0][j][i].z + 5.0 * ucat[kp1][j][i].z + 5.0 * ucat[kp2][j][i].z) / 8.0;
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

// FIELD OPERATIONS
// ============================================================================================================= //

inline void resetNoPenetrationFluxes(ueqn_ *ueqn)
{
    mesh_             *mesh = ueqn->access->mesh;
    DM               da = mesh->da, fda = mesh->fda;
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         lxs, lxe, lys, lye, lzs, lze;
    PetscInt         i, j, k;

    Cmpnts           ***ucont;

    lxs = xs; lxe = xe; if (xs==0) lxs = xs+1; if (xe==mx) lxe = xe-1;
    lys = ys; lye = ye; if (ys==0) lys = ys+1; if (ye==my) lye = ye-1;
    lzs = zs; lze = ze; if (zs==0) lzs = zs+1; if (ze==mz) lze = ze-1;

    // We go through all the boundary conditions implying no penetration and set
    // the wall normal flux to zero in order to exactly enforce mass conservation.
    // This is done because, especially the non-linear solver, could give fluxes not
    // exactly zero at the wall.

    DMDAVecGetArray(fda, ueqn->Ucont, &ucont);

    for (k=zs; k<lze; k++)
    {
        for (j=ys; j<lye; j++)
        {
            for (i=xs; i<lxe; i++)
            {
                // noslip BC
                if(i==0 && mesh->boundaryU.iLeft=="noSlip")     ucont[k][j][i].x = 0.0;
                if(i==mx-2 && mesh->boundaryU.iRight=="noSlip") ucont[k][j][i].x = 0.0;
                if(j==0 && mesh->boundaryU.jLeft=="noSlip")     ucont[k][j][i].y = 0.0;
                if(j==my-2 && mesh->boundaryU.jRight=="noSlip") ucont[k][j][i].y = 0.0;
                if(k==0 && mesh->boundaryU.kLeft=="noSlip")     ucont[k][j][i].z = 0.0;
                if(k==mz-2 && mesh->boundaryU.kRight=="noSlip") ucont[k][j][i].z = 0.0;

                // wall models
                if(i==0 && mesh->boundaryU.iLeft=="velocityWallFunction")     ucont[k][j][i].x = 0.0;
                if(i==mx-2 && mesh->boundaryU.iRight=="velocityWallFunction") ucont[k][j][i].x = 0.0;
                if(j==0 && mesh->boundaryU.jLeft=="velocityWallFunction")     ucont[k][j][i].y = 0.0;
                if(j==my-2 && mesh->boundaryU.jRight=="velocityWallFunction") ucont[k][j][i].y = 0.0;
                if(k==0 && mesh->boundaryU.kLeft=="velocityWallFunction")     ucont[k][j][i].z = 0.0;
                if(k==mz-2 && mesh->boundaryU.kRight=="velocityWallFunction") ucont[k][j][i].z = 0.0;

                // slip BC
                if (mesh->boundaryU.iLeft=="slip" && i==0)     ucont[k][j][i].x = 0.0;
                if (mesh->boundaryU.iRight=="slip" && i==mx-2) ucont[k][j][i].x = 0.0;
                if (mesh->boundaryU.jLeft=="slip" && j==0)     ucont[k][j][i].y = 0.0;
                if (mesh->boundaryU.jRight=="slip" && j==my-2) ucont[k][j][i].y = 0.0;
                if (mesh->boundaryU.kLeft=="slip" && k==0)     ucont[k][j][i].z = 0.0;
                if (mesh->boundaryU.kRight=="slip" && k==mz-2) ucont[k][j][i].z = 0.0;

                // zero gradient BC (if reverse flow)
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
            }
        }
    }

    DMDAVecRestoreArray (fda, ueqn->Ucont, &ucont);

    DMGlobalToLocalBegin(fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);
    DMGlobalToLocalEnd  (fda, ueqn->Ucont, INSERT_VALUES, ueqn->lUcont);

    return;
}

//***************************************************************************************************************//

inline void resetFacePeriodicFluxesVector(mesh_ *mesh, Vec V, Vec lV, const char *scatterType)
{
    DMDALocalInfo    info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         i, j, k;

    Cmpnts           ***v, ***lv;

    // If the flow is periodic in l-direction, the variable is solved also
    // at the right boundary faces. Its value must be copied at the left boundary
    // faces to ensure consistency in the momentum computation.

    if(scatterType=="globalToLocal")
    {
        DMDAVecGetArray(mesh->fda, V, &v);
        DMDAVecGetArray(mesh->fda, lV, &lv);
    }
    else if(scatterType=="localToLocal")
    {
        DMDAVecGetArray(mesh->fda, V, &v);
    }

    if(mesh->access->info->periodic)
    {
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                for (i=xs; i<xe; i++)
                {
                    PetscInt flag=0, a=i, b=j, c=k;

                    if(mesh->i_periodic && i==0) a=mx-2, flag=1;
                    //else if(i_periodic && i==mx-1) a=1, flag=1;

                    if(mesh->j_periodic && j==0) b=my-2, flag=1;
                    //else if(j_periodic && j==my-1) b=1, flag=1;

                    if(mesh->k_periodic && k==0) c=mz-2, flag=1;
                    //else if(k_periodic && k==mz-1) c=1, flag=1;

                    if(mesh->ii_periodic && i==0) a=-2, flag=1;
                    //else if(ii_periodic && i==mx-1) a=mx+1, flag=1;

                    if(mesh->jj_periodic && j==0) b=-2, flag=1;
                    //else if(jj_periodic && j==my-1) b=my+1, flag=1;

                    if(mesh->kk_periodic && k==0) c=-2, flag=1;
                    //else if(kk_periodic && k==mz-1) c=mz+1, flag=1;

                    if(flag)
                    {
                        if(scatterType=="globalToLocal")
                        {
                            v[k][j][i] = lv[c][b][a];
                        }
                        else if(scatterType=="localToLocal")
                        {
                            v[k][j][i] = v[c][b][a];
                        }
                    }
                }
            }
        }
    }

    if(scatterType=="globalToLocal")
    {
        DMDAVecRestoreArray (mesh->fda, V, &v);
        DMDAVecRestoreArray (mesh->fda, lV, &lv);
        DMGlobalToLocalBegin(mesh->fda, V, INSERT_VALUES, lV);
        DMGlobalToLocalEnd  (mesh->fda, V, INSERT_VALUES, lV);
    }
    else if(scatterType=="localToLocal")
    {
        DMDAVecRestoreArray (mesh->fda, V, &v);
        DMLocalToLocalBegin (mesh->fda, V, INSERT_VALUES, V);
        DMLocalToLocalEnd   (mesh->fda, V, INSERT_VALUES, V);
    }
    else
    {
        char error[512];
        sprintf(error, "wrong scatter type selected, available types are\n    - globalToLocal\n    - localToLocal\n");
        fatalErrorInFunction("resetCellPeriodicFluxes",  error);
    }

    return;
}

//***************************************************************************************************************//

//! \brief Reset periodic values and scatter
inline void resetCellPeriodicFluxes(mesh_ *mesh, Vec &V, Vec &lV, const char *type, const char *scatterType)
{
    DM            da    = mesh->da, fda = mesh->fda;
    DMDALocalInfo  info = mesh->info;
    PetscInt         xs = info.xs, xe = info.xs + info.xm;
    PetscInt         ys = info.ys, ye = info.ys + info.ym;
    PetscInt         zs = info.zs, ze = info.zs + info.zm;
    PetscInt         mx = info.mx, my = info.my, mz = info.mz;

    PetscInt         i, j, k;

    Cmpnts           ***v, ***lv;
    PetscScalar      ***s, ***ls;

    // If the flow is periodic in l-th direction, put the first internal value
    // at the ghost cell of the opposite patch.

    if(type=="scalar")
    {
        if(scatterType=="globalToLocal")
        {
            DMDAVecGetArray(da, V, &s);
            DMDAVecGetArray(da, lV, &ls);
        }
        else if(scatterType=="localToLocal")
        {
            DMDAVecGetArray(da, V, &s);
        }
    }
    else if(type=="vector")
    {
        if(scatterType=="globalToLocal")
        {
            DMDAVecGetArray(fda, V, &v);
            DMDAVecGetArray(fda, lV, &lv);
        }
        else if(scatterType=="localToLocal")
        {
            DMDAVecGetArray(fda, V, &v);
        }
    }
    else
    {
        char error[512];
        sprintf(error, "wrong type selected, available types are\n    - scalar\n    - vector\n");
        fatalErrorInFunction("resetCellPeriodicFluxes",  error);
    }

    if(mesh->access->info->periodic)
    {
        for (k=zs; k<ze; k++)
        {
            for (j=ys; j<ye; j++)
            {
                for (i=xs; i<xe; i++)
                {
                    PetscInt flag=0, a=i, b=j, c=k;

                    if(mesh->i_periodic && i==0) a=mx-2, flag=1;
                    else if(mesh->i_periodic && i==mx-1) a=1, flag=1;

                    if(mesh->j_periodic && j==0) b=my-2, flag=1;
                    else if(mesh->j_periodic && j==my-1) b=1, flag=1;

                    if(mesh->k_periodic && k==0) c=mz-2, flag=1;
                    else if(mesh->k_periodic && k==mz-1) c=1, flag=1;

                    if(mesh->ii_periodic && i==0) a=-2, flag=1;
                    else if(mesh->ii_periodic && i==mx-1) a=mx+1, flag=1;

                    if(mesh->jj_periodic && j==0) b=-2, flag=1;
                    else if(mesh->jj_periodic && j==my-1) b=my+1, flag=1;

                    if(mesh->kk_periodic && k==0) c=-2, flag=1;
                    else if(mesh->kk_periodic && k==mz-1) c=mz+1, flag=1;

                    if(flag)
                    {
                        if(type=="scalar")
                        {
                            if(scatterType=="globalToLocal")
                            {
                                s[k][j][i] = ls[c][b][a];
                            }
                            else if(scatterType=="localToLocal")
                            {
                                s[k][j][i] = s[c][b][a];
                            }
                        }
                        else if(type=="vector")
                        {
                            if(scatterType=="globalToLocal")
                            {
                                v[k][j][i] = lv[c][b][a];
                            }
                            else if(scatterType=="localToLocal")
                            {
                                v[k][j][i] = v[c][b][a];
                            }
                        }
                    }
                }
            }
        }
    }

    if(type=="scalar")
    {
        if(scatterType=="globalToLocal")
        {
            DMDAVecRestoreArray(da, V, &s);
            DMDAVecRestoreArray(da, lV, &ls);
            DMGlobalToLocalBegin(da, V, INSERT_VALUES, lV);
            DMGlobalToLocalEnd(da, V, INSERT_VALUES, lV);
        }
        else if(scatterType=="localToLocal")
        {
            DMDAVecRestoreArray(da, V, &s);
            DMLocalToLocalBegin(da, V, INSERT_VALUES, V);
            DMLocalToLocalEnd(da, V, INSERT_VALUES, V);
        }
        else
        {
            char error[512];
            sprintf(error, "wrong scatter type selected, available types are\n    - globalToLocal\n    - localToLocal\n");
            fatalErrorInFunction("resetCellPeriodicFluxes",  error);
        }
    }
    else if(type=="vector")
    {
        if(scatterType=="globalToLocal")
        {
            DMDAVecRestoreArray(fda, V, &v);
            DMDAVecRestoreArray(fda, lV, &lv);
            DMGlobalToLocalBegin(fda, V, INSERT_VALUES, lV);
            DMGlobalToLocalEnd(fda, V, INSERT_VALUES, lV);
        }
        else if(scatterType=="localToLocal")
        {
            DMDAVecRestoreArray(fda, V, &v);
            DMLocalToLocalBegin(fda, V, INSERT_VALUES, V);
            DMLocalToLocalEnd(fda, V, INSERT_VALUES, V);
        }
        else
        {
            char error[512];
            sprintf(error, "wrong scatter type selected, available types are\n    - globalToLocal\n    - localToLocal\n");
            fatalErrorInFunction("resetCellPeriodicFluxes",  error);
        }
    }

    return;
}

// DAMPING VISCOSITY & FRINGE REGION RELATED
// ============================================================================================================= //

inline PetscReal viscRayleigh(PetscReal &alpha, PetscReal &hS, PetscReal &hE, PetscReal &h)
{
    PetscReal h_hat = (h-hS)/(hE-hS);
    PetscReal viscosity;

    viscosity
    =
    alpha *
    (
        1 - std::cos
        (
            M_PI / 2.0 * PetscMax(h_hat, 0.0)
        )
    );

    return(viscosity);
}

//***************************************************************************************************************//

inline PetscReal viscCosAscending(PetscReal &alpha, PetscReal &hS, PetscReal &hE, PetscReal &h)
{
    PetscReal h_hat = M_PI*(h-hS)/(hE-hS);
    PetscReal viscosity = 0.0;

    if(h > hS && h < hE)
    {
        viscosity
        =
        alpha * 0.5 *
        (
            1 - std::cos(h_hat)
        );
    }

    return(viscosity);
}

//***************************************************************************************************************//

inline PetscReal viscCosDescending(PetscReal &alpha, PetscReal &hS, PetscReal &hE, PetscReal &h)
{
    PetscReal h_hat = M_PI*(h-hS)/(hE-hS);
    PetscReal viscosity = 0;

    if(h > hS && h < hE)
    {
        viscosity
        =
        alpha * 0.5 *
        (
            1 + std::cos(h_hat)
        );
    }

    return(viscosity);
}

//***************************************************************************************************************//
inline PetscReal viscNordstrom(PetscReal &alpha, PetscReal &hS, PetscReal &hE, PetscReal &delta, PetscReal &h)
{
    // make sure delta is positive and strictly greater than zero
    delta = PetscMax(fabs(delta), 1e-5);

    PetscReal h1_hat = (h-hS)/delta;
    PetscReal h2_hat = (h-hE)/delta + 1.0;

    PetscReal viscosity;

    PetscReal s1, s2;

    // compute smoothing function rise component
    if(h1_hat <= 0.0)
    {
        s1 = 0.0;
    }
    else if(h1_hat >= 1.0)
    {
        s1 = 1.0;
    }
    else
    {
        s1 = 1.0 / (1.0 + std::exp(1.0/(h1_hat - 1.0) + 1.0 / h1_hat));
    }


    // compute smoothing function falling component
    if(h2_hat <= 0.0)
    {
        s2 = 0.0;
    }
    else if(h2_hat >= 1.0)
    {
        s2 = 1.0;
    }
    else
    {
        s2 = 1.0 / (1.0 + std::exp(1.0/(h2_hat - 1.0) + 1.0 / h2_hat));
    }

    viscosity = alpha * (s1 - s2);

    return(viscosity);
}

//***************************************************************************************************************//

inline double viscNordstromNoVertFilter(double hS, double hE, double &delta, double &h)
{
    // make sure delta is positive and strictly greater than zero
    delta = std::max(fabs(delta), 1e-5);

    double h1_hat = (h-hS)/delta;
    double h2_hat = (h-hE)/delta + 1.0;

    double viscosity;

    double s1, s2;

    // compute smoothing function rise component
    if(h1_hat <= 0.0)
    {
        s1 = 0.0;
    }
    else if(h1_hat >= 1.0)
    {
        s1 = 1.0;
    }
    else
    {
        s1 = 1.0 / (1.0 + std::exp(1.0/(h1_hat - 1.0) + 1.0 / h1_hat));
    }

    // compute smoothing function falling component
    if(h2_hat <= 0.0)
    {
        s2 = 0.0;
    }
    else if(h2_hat >= 1.0)
    {
        s2 = 1.0;
    }
    else
    {
        s2 = 1.0 / (1.0 + std::exp(1.0/(h2_hat - 1.0) + 1.0 / h2_hat));
    }

	viscosity = 1.0 - (s1 - s2);

    return(viscosity);
}

//***************************************************************************************************************//

inline double viscStipa(double &hS, double &hE, double &delta, double &h, double &z, double &H)
{
    // make sure delta is positive and strictly greater than zero
    delta = std::max(fabs(delta), 1e-5);

    double h1_hat = (h-hS)/delta;
    double h2_hat = (h-hE)/delta + 1.0;

    double viscosity;

    double s1, s2;

    // compute smoothing function rise component
    if(h1_hat <= 0.0)
    {
        s1 = 0.0;
    }
    else if(h1_hat >= 1.0)
    {
        s1 = 1.0;
    }
    else
    {
        s1 = 1.0 / (1.0 + std::exp(1.0/(h1_hat - 1.0) + 1.0 / h1_hat));
    }

    // compute smoothing function falling component
    if(h2_hat <= 0.0)
    {
        s2 = 0.0;
    }
    else if(h2_hat >= 1.0)
    {
        s2 = 1.0;
    }
    else
    {
        s2 = 1.0 / (1.0 + std::exp(1.0/(h2_hat - 1.0) + 1.0 / h2_hat));
    }

	if(z >= H)
	{
		viscosity = 1.0 - (s1 - s2);
	}
	else
	{
		viscosity = 1.0;
	}

    return(viscosity);
}

//***************************************************************************************************************//

inline double viscStipaDelta(double &hS, double &hE, double &deltaS, double &deltaE, double &h, double &z, double &H)
{
    // make sure delta is positive and strictly greater than zero
    deltaS = std::max(fabs(deltaS), 1e-5);
    deltaE = std::max(fabs(deltaE), 1e-5);

    double h1_hat = (h-hS)/deltaS;
    double h2_hat = (h-hE)/deltaE + 1.0;

    double viscosity;

    double s1, s2;

    // compute smoothing function rise component
    if(h1_hat <= 0.0)
    {
        s1 = 0.0;
    }
    else if(h1_hat >= 1.0)
    {
        s1 = 1.0;
    }
    else
    {
        s1 = 1.0 / (1.0 + std::exp(1.0/(h1_hat - 1.0) + 1.0 / h1_hat));
    }

    // compute smoothing function falling component
    if(h2_hat <= 0.0)
    {
        s2 = 0.0;
    }
    else if(h2_hat >= 1.0)
    {
        s2 = 1.0;
    }
    else
    {
        s2 = 1.0 / (1.0 + std::exp(1.0/(h2_hat - 1.0) + 1.0 / h2_hat));
    }

    if(z >= H)
	{
		viscosity = 1.0 - (s1 - s2);
	}
	else
	{
		viscosity = 1.0;
	}

    return(viscosity);
}

//***************************************************************************************************************//

inline void findInterpolationWeights(PetscReal *weights, PetscInt *labels, PetscReal *pvec, PetscInt npts, PetscReal pval)
{
    // equality tolerance: if equal data points up to re-tol = 1e-5 are found,
    // that interpolation point is skipped because it does not provide information
    // and would result in a division by zero when computing the weights.
    PetscReal tol1, tol2, eqTol = 1e-5;

    // check that extrapolation is not necessary
    // (up to a minimum tol of (pEnd-pStart) / npts that it is small enough not to have extrapolation errors)
    PetscReal eps = (pvec[npts-1] - pvec[0]) / npts;
    // if(pval > pvec[npts-1]+eps || pval < pvec[0]-eps)
    // {
    //    char error[512];
    //     sprintf(error, "query point %lf outside of bounds [%lf - %lf]... was told not to extrapolate\n", pval, pvec[0]-eps, pvec[npts-1]+eps);
    //     fatalErrorInFunction("findInterpolationWeights",  error);
    // }

    PetscReal diff[npts];
    PetscInt    idx_1 = 0,
           idx_2 = 1;

    // find the difference
    for(PetscInt i=0; i<npts; i++)
    {
        diff[i] = fabs(pvec[i]-pval);
    }

    // find the two closest labels
    for(PetscInt i=0; i<npts; i++)
    {
        tol1 = fabs(diff[i] - diff[idx_2]) / diff[i];

        if(diff[i] < diff[idx_1] && tol1 > eqTol)
        {
            idx_2 = idx_1;
            idx_1 = i;
        }

        tol2 = fabs(diff[i] - diff[idx_1]) / diff[i];

        if(diff[i] < diff[idx_2] && i != idx_1 && tol2 > eqTol)
        {
            idx_2 = i;
        }
    }

    // always put the lower value at idx_1 and higher at idx_2
    if(pvec[idx_2] < pvec[idx_1])
    {
        PetscInt idx_tmp = idx_2;
        idx_2 = idx_1;
        idx_1 = idx_tmp;
    }

    // store the labels
    labels[0] = idx_1;
    labels[1] = idx_2;

    // should never divide by zero (increase equality tolerance if it happens)
    PetscReal idx = (idx_2 - idx_1) / (pvec[idx_2] - pvec[idx_1]) * (pval - pvec[idx_1]) + idx_1;

    // store the weights
    weights[0] = (idx_2 - idx) / (idx_2 - idx_1);
    weights[1] = (idx - idx_1) / (idx_2 - idx_1);

    return;
}

//***************************************************************************************************************//

inline void findInterpolationWeightsWithExtrap(PetscReal *weights, PetscInt *labels, PetscReal *pvec, PetscInt npts, PetscReal pval)
{
    if(pval >= pvec[npts-1])
    {
        weights[0] = 0.5;
        weights[1] = 0.5;

        labels[0]  = npts-1;
        labels[1]  = npts-1;
    }
    else if (pval <= pvec[0])
    {
        weights[0] = 0.5;
        weights[1] = 0.5;

        labels[0]  = 0;
        labels[1]  = 0;
    }
    else
    {
        findInterpolationWeights(weights, labels, pvec, npts, pval);
    }

    return;
}

// ============================================================================================================= //

inline PetscReal scaleHyperTangBot(PetscReal h, PetscReal H, PetscReal delta)
{
    // the input parameter H is the height at which the function is almost zero.
    // cb is the center of the smoothing region, i.e. H - 0.5*delta
    PetscReal cb = H - 0.5*delta;
    return
    (
        0.5 *
        (
            1.0 - std::tanh(7.0 * (h-cb)/delta)
        )
    );
}

// ============================================================================================================= //

inline PetscReal scaleHyperTangTop(PetscReal h, PetscReal H, PetscReal delta)
{
    // the input parameter H is the height at which the function is almost one.
    // cb is the center of the smoothing region, i.e. H - 0.5*delta
    PetscReal cb = H - 0.5*delta;
    return
    (
        0.5 *
        (
            1.0 + std::tanh(7.0 * (h-cb)/delta)
        )
    );
}


// FLUID QUANTITIES
// ============================================================================================================= //

//! \brief Computes mechanical energy (Em = <u><u> + <v><v> + <w><w> + <u'u'> + <v'v'> + <w'w'> + <p>/rho)
inline PetscReal computeEm(Cmpnts &avgU, symmTensor &avgUprimeUprime, PetscReal avgP)
{
    return
    (
        0.5 *
        (
            avgU.x*avgU.x + avgU.y*avgU.y + avgU.z*avgU.z +
            avgUprimeUprime.xx + avgUprimeUprime.yy + avgUprimeUprime.zz
        ) + avgP
    );
}

//! \brief Computes modified mechanical energy (Em = <U><U> + <V><V> + <W><W> + <U'U'> + <V'V'> + <W'W'>) using contrav. fluxes
inline PetscReal computeEmTilde(Cmpnts &avgU, symmTensor &avgUprimeUprime)
{
    return
    (
        0.5 *
        (
            avgU.x*avgU.x + avgU.y*avgU.y + avgU.z*avgU.z +
            avgUprimeUprime.xx + avgUprimeUprime.yy + avgUprimeUprime.zz
        )
    );
}

inline PetscReal computeMKE(Cmpnts &avgU)
{
    return
    (
        0.5 *
        (
            avgU.x*avgU.x + avgU.y*avgU.y + avgU.z*avgU.z
        )
    );
}

inline PetscReal computeTKE(symmTensor &avgUprimeUprime)
{
    return
    (
        0.5 *
        (
            avgUprimeUprime.xx + avgUprimeUprime.yy + avgUprimeUprime.zz
        )
    );
}

// INTERPOLATION
// ============================================================================================================= //

//! \brief Trilinear volume interpolation function for scalars
inline void scalarPointLocalVolumeInterpolation
(
    mesh_ *mesh,
    PetscReal px, PetscReal py, PetscReal pz,
    PetscInt ic, PetscInt jc, PetscInt kc,
    Cmpnts ***cent,
    PetscReal ***v,
    PetscReal &result
)
{
    // Given generic point coordinates px, py, pz and the indices of the closest
    // cell center, this functions finds the closest 8 cells to the given point which form
    // its surrounding box and interpolates the point value using the values at the 8
    // surrounding cell centers (tri-linear interpolation).
    // Note: cell indices passed to this function must be internal cells, not ghosts.
    // Throws an error otherwise.
    // Interpolation weights are bounded from 0 to 1 to handle query points lying outside
    // of the boundaries: this means that those points are projected on the boundaries and
    // the value there is returned (Prevents extrapolation ouside of boundaries).

    DMDALocalInfo    info = mesh->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt         gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt         gys  = info.gys, gye = info.gys + info.gym;
    PetscInt         gzs  = info.gzs, gze = info.gzs + info.gzm;

    PetscInt         i, j, k;

    PetscInt iL, iR;
    // do the search if surrounded by cells
    if(ic > 1 && ic < mx-2)
    {
        if (ic-1 < gxs)
        {
            iL = ic; iR = ic + 1;
        }
        else if(ic+1 >= gxe)
        {
            iL = ic - 1; iR = ic;
        }
        else
        {
            PetscReal dcsiRightSq
            =
            std::pow(cent[kc][jc][ic+1].x - px,2) +
            std::pow(cent[kc][jc][ic+1].y - py,2) +
            std::pow(cent[kc][jc][ic+1].z - pz,2);

            PetscReal dcsiLeftSq
            =
            std::pow(px - cent[kc][jc][ic-1].x,2) +
            std::pow(py - cent[kc][jc][ic-1].y,2) +
            std::pow(pz - cent[kc][jc][ic-1].z,2);

            if(dcsiRightSq <= dcsiLeftSq)
            {
                iR = ic + 1; iL = ic;
            }
            else
            {
                iR = ic; iL = ic-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(ic == 1)
    {
        iR = ic + 1; iL = ic;
    }
    // on boundary, already know where to pick data
    else if(ic == mx-2)
    {
        iR = ic; iL = ic-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "i-index out of bounds: ended up in the ghosts when looking for interpolation data\n");
        fatalErrorInFunction("pointLocalVolumeInterpolation",  error);
    }

    PetscInt jL, jR;
    // do the search if surrounded by cells
    if(jc > 1 && jc < my-2)
    {
        // check processor bounds
        if (jc-1 < gys)
        {
            jL = jc; jR = jc + 1;
        }
        else if(jc+1 >= gye)
        {
            jL = jc - 1; jR = jc;
        }
        else
        {
            PetscReal detaRightSq
            =
            std::pow(cent[kc][jc+1][ic].x - px,2) +
            std::pow(cent[kc][jc+1][ic].y - py,2) +
            std::pow(cent[kc][jc+1][ic].z - pz,2);

            PetscReal detaLeftSq
            =
            std::pow(px - cent[kc][jc-1][ic].x,2) +
            std::pow(py - cent[kc][jc-1][ic].y,2) +
            std::pow(pz - cent[kc][jc-1][ic].z,2);

            if(detaRightSq <= detaLeftSq)
            {
                jR = jc + 1; jL = jc;
            }
            else
            {
                jR = jc; jL = jc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(jc == 1)
    {
        jR = jc + 1; jL = jc;
    }
    // on boundary, already know where to pick data
    else if(jc == my-2)
    {
        jR = jc; jL = jc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "j-index out of bounds: ended up in the ghosts when looking for interpolation data\n");
        fatalErrorInFunction("pointLocalVolumeInterpolation",  error);
    }

    PetscInt kL, kR;
    // do the search if surrounded by cells
    if(kc > 1 && kc < mz-2)
    {
        // check processor bounds
        if (kc-1 < gzs)
        {
            kL = kc; kR = kc + 1;
        }
        else if(kc+1 >= gze)
        {
            kL = kc - 1; kR = kc;
        }
        else
        {
            PetscReal dzetRightSq
            =
            std::pow(cent[kc+1][jc][ic].x - px,2) +
            std::pow(cent[kc+1][jc][ic].y - py,2) +
            std::pow(cent[kc+1][jc][ic].z - pz,2);

            PetscReal dzetLeftSq
            =
            std::pow(px - cent[kc-1][jc][ic].x,2) +
            std::pow(py - cent[kc-1][jc][ic].y,2) +
            std::pow(pz - cent[kc-1][jc][ic].z,2);

            if(dzetRightSq <= dzetLeftSq)
            {
                kR = kc + 1; kL = kc;
            }
            else
            {
                kR = kc; kL = kc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(kc == 1)
    {
        kR = kc + 1; kL = kc;
    }
    // on boundary, already know where to pick data
    else if(kc == mz-2)
    {
        kR = kc; kL = kc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "k-index out of bounds: ended up in the ghosts when looking for interpolation data\n");
        fatalErrorInFunction("pointLocalVolumeInterpolation",  error);
    }

    /*
    // Debugging information (only i and j)
    PetscPrintf(PETSC_COMM_WORLD, "Neighboring cells\n");
    PetscPrintf(PETSC_COMM_WORLD, "(iL,jL) idxs: %d %d, C = %.2f %.2f\n", iL, jL, cent[kL][jL][iL].y, cent[kL][jL][iL].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iL,jR) idxs: %d %d, C = %.2f %.2f\n", iL, jR, cent[kL][jR][iL].y, cent[kL][jR][iL].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iR,jL) idxs: %d %d, C = %.2f %.2f\n", iR, jL, cent[kL][jL][iR].y, cent[kL][jL][iR].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iR,jR) idxs: %d %d, C = %.2f %.2f\n", iR, jR, cent[kL][jR][iR].y, cent[kL][jR][iR].z);
    PetscPrintf(PETSC_COMM_WORLD, "(ic,jc) idxs: %d %d\n", ic, jc);
    PetscPrintf(PETSC_COMM_WORLD, "P = %.2f %.2f\n\n", py, pz);
    */

    // find directional gradients
    Cmpnts csiHat, etaHat, zetHat;
    PetscReal csiNor, etaNor, zetNor;

    csiHat.x = cent[kL][jL][iR].x - cent[kL][jL][iL].x;
    csiHat.y = cent[kL][jL][iR].y - cent[kL][jL][iL].y;
    csiHat.z = cent[kL][jL][iR].z - cent[kL][jL][iL].z;

    etaHat.x = cent[kL][jR][iL].x - cent[kL][jL][iL].x;
    etaHat.y = cent[kL][jR][iL].y - cent[kL][jL][iL].y;
    etaHat.z = cent[kL][jR][iL].z - cent[kL][jL][iL].z;

    zetHat.x = cent[kR][jL][iL].x - cent[kL][jL][iL].x;
    zetHat.y = cent[kR][jL][iL].y - cent[kL][jL][iL].y;
    zetHat.z = cent[kR][jL][iL].z - cent[kL][jL][iL].z;

    csiNor = std::sqrt(csiHat.x*csiHat.x + csiHat.y*csiHat.y + csiHat.z*csiHat.z);
    etaNor = std::sqrt(etaHat.x*etaHat.x + etaHat.y*etaHat.y + etaHat.z*etaHat.z);
    zetNor = std::sqrt(zetHat.x*zetHat.x + zetHat.y*zetHat.y + zetHat.z*zetHat.z);

    csiHat.x = csiHat.x / csiNor; csiHat.y = csiHat.y / csiNor; csiHat.z = csiHat.z / csiNor;
    etaHat.x = etaHat.x / etaNor; etaHat.y = etaHat.y / etaNor; etaHat.z = etaHat.z / etaNor;
    zetHat.x = zetHat.x / zetNor; zetHat.y = zetHat.y / zetNor; zetHat.z = zetHat.z / zetNor;

    // bound deltas from 0 to 1: this means that if the query point is outside of boundaries (we have
    // no data there), we project that point on the boundary and give the interpolated value at that
    // point as a result.

    PetscReal deltaCsi
    =
    PetscMax
    (
        PetscMin
        (
            (
                csiHat.x * (px - cent[kL][jL][iL].x) +
                csiHat.y * (py - cent[kL][jL][iL].y) +
                csiHat.z * (pz - cent[kL][jL][iL].z)

            ) / csiNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaEta
    =
    PetscMax
    (
        PetscMin
        (
            (
                etaHat.x * (px - cent[kL][jL][iL].x) +
                etaHat.y * (py - cent[kL][jL][iL].y) +
                etaHat.z * (pz - cent[kL][jL][iL].z)

            ) / etaNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaZet
    =
    PetscMax
    (
        PetscMin
        (
            (
                zetHat.x * (px - cent[kL][jL][iL].x) +
                zetHat.y * (py - cent[kL][jL][iL].y) +
                zetHat.z * (pz - cent[kL][jL][iL].z)

            ) / zetNor,
            1.0
        ),
        0.0
    );

    // find interpolation weights
    PetscReal v000 = v[kL][jL][iL],
           v100 = v[kL][jL][iR],
           v110 = v[kL][jR][iR],
           v010 = v[kL][jR][iL],
           v001 = v[kR][jL][iL],
           v101 = v[kR][jL][iR],
           v111 = v[kR][jR][iR],
           v011 = v[kR][jR][iL];

    PetscReal c0 = v000,
           c1 = v100 - v000,
           c2 = v010 - v000,
           c3 = v001 - v000,
           c4 = v110 - v010 - v100 + v000,
           c5 = v011 - v001 - v010 + v000,
           c6 = v101 - v001 - v100 + v000,
           c7 = v111 - v011 - v101 - v110 + v100 + v001 + v010 - v000;

    // interpolate
    result
    =
    c0 +
    c1 * deltaCsi +
    c2 * deltaEta +
    c3 * deltaZet +
    c4 * deltaCsi * deltaEta +
    c5 * deltaEta * deltaZet +
    c6 * deltaZet * deltaCsi +
    c7 * deltaCsi * deltaEta * deltaZet;

    return;
}

//***************************************************************************************************************//

//! \brief Trilinear volume interpolation function for scalars to return interpolation Weights
inline void PointInterpolationWeights
(
    mesh_ *mesh,
    PetscReal px, PetscReal py, PetscReal pz,
    PetscInt ic, PetscInt jc, PetscInt kc,
    Cmpnts ***cent,
    PetscReal *intWts,
    PetscInt *intId
)
{
    // Given generic point coordinates px, py, pz and the indices of the closest
    // cell center, this functions finds the closest 8 cells to the given point which form
    // its surrounding box and interpolates the point value using the values at the 8
    // surrounding cell centers (tri-linear interpolation).
    // Note: cell indices passed to this function must be internal cells, not ghosts.
    // Throws an error otherwise.
    // Interpolation weights are bounded from 0 to 1 to handle query points lying outside
    // of the boundaries: this means that those points are projected on the boundaries and
    // the value there is returned (Prevents extrapolation ouside of boundaries).

    DMDALocalInfo    info = mesh->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt         gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt         gys  = info.gys, gye = info.gys + info.gym;
    PetscInt         gzs  = info.gzs, gze = info.gzs + info.gzm;

    PetscInt         i, j, k;

    PetscInt iL, iR;

    // do the search if surrounded by cells
    if(ic > 1 && ic < mx-2)
    {
        if (ic-1 < gxs)
        {
            iL = ic; iR = ic + 1;
        }
        else if(ic+1 >= gxe)
        {
            iL = ic - 1; iR = ic;
        }
        else
        {
            PetscReal dcsiRightSq
            =
            std::pow(cent[kc][jc][ic+1].x - px,2) +
            std::pow(cent[kc][jc][ic+1].y - py,2) +
            std::pow(cent[kc][jc][ic+1].z - pz,2);

            PetscReal dcsiLeftSq
            =
            std::pow(px - cent[kc][jc][ic-1].x,2) +
            std::pow(py - cent[kc][jc][ic-1].y,2) +
            std::pow(pz - cent[kc][jc][ic-1].z,2);

            if(dcsiRightSq <= dcsiLeftSq)
            {
                iR = ic + 1; iL = ic;
            }
            else
            {
                iR = ic; iL = ic-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(ic == 1)
    {
        iR = ic + 1; iL = ic;
    }
    // on boundary, already know where to pick data
    else if(ic == mx-2)
    {
        iR = ic; iL = ic-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "i-index out of bounds: ended up in the ghosts when looking for interpolation (i = %ld)\n", ic);
        fatalErrorInFunction("scalarPointInterpolationWeights",  error);
    }

    PetscInt jL, jR;
    // do the search if surrounded by cells
    if(jc > 1 && jc < my-2)
    {
        // check processor bounds
        if (jc-1 < gys)
        {
            jL = jc; jR = jc + 1;
        }
        else if(jc+1 >= gye)
        {
            jL = jc - 1; jR = jc;
        }
        else
        {
            PetscReal detaRightSq
            =
            std::pow(cent[kc][jc+1][ic].x - px,2) +
            std::pow(cent[kc][jc+1][ic].y - py,2) +
            std::pow(cent[kc][jc+1][ic].z - pz,2);

            PetscReal detaLeftSq
            =
            std::pow(px - cent[kc][jc-1][ic].x,2) +
            std::pow(py - cent[kc][jc-1][ic].y,2) +
            std::pow(pz - cent[kc][jc-1][ic].z,2);

            if(detaRightSq <= detaLeftSq)
            {
                jR = jc + 1; jL = jc;
            }
            else
            {
                jR = jc; jL = jc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(jc == 1)
    {
        jR = jc + 1; jL = jc;
    }
    // on boundary, already know where to pick data
    else if(jc == my-2)
    {
        jR = jc; jL = jc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "j-index out of bounds: ended up in the ghosts when looking for interpolation (j = %ld)\n", jc);
        fatalErrorInFunction("scalarPointInterpolationWeights",  error);
    }

    PetscInt kL, kR;
    // do the search if surrounded by cells
    if(kc > 1 && kc < mz-2)
    {
        // check processor bounds
        if (kc-1 < gzs)
        {
            kL = kc; kR = kc + 1;
        }
        else if(kc+1 >= gze)
        {
            kL = kc - 1; kR = kc;
        }
        else
        {
            PetscReal dzetRightSq
            =
            std::pow(cent[kc+1][jc][ic].x - px,2) +
            std::pow(cent[kc+1][jc][ic].y - py,2) +
            std::pow(cent[kc+1][jc][ic].z - pz,2);

            PetscReal dzetLeftSq
            =
            std::pow(px - cent[kc-1][jc][ic].x,2) +
            std::pow(py - cent[kc-1][jc][ic].y,2) +
            std::pow(pz - cent[kc-1][jc][ic].z,2);

            if(dzetRightSq <= dzetLeftSq)
            {
                kR = kc + 1; kL = kc;
            }
            else
            {
                kR = kc; kL = kc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(kc == 1)
    {
        kR = kc + 1; kL = kc;
    }
    // on boundary, already know where to pick data
    else if(kc == mz-2)
    {
        kR = kc; kL = kc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "k-index out of bounds: ended up in the ghosts when looking for interpolation (k = %ld)\n", kc);
        fatalErrorInFunction("scalarPointInterpolationWeights",  error);
    }

    intId[0] = kL;
    intId[1] = kR;
    intId[2] = jL;
    intId[3] = jR;
    intId[4] = iL;
    intId[5] = iR;

    //check that the point is inside the interpolation box, send warning if not
    // if
    // (
    //     ( ( (px < cent[kL][jL][iL].x) && (px < cent[kR][jR][iR].x) ) || ( (px > cent[kR][jR][iR].x) && (px > cent[kL][jL][iL].x) ) ) ||
    //     ( ( (py < cent[kL][jL][iL].y) && (py < cent[kR][jR][iR].y) ) || ( (py > cent[kR][jR][iR].y) && (py > cent[kL][jL][iL].y) ) )||
    //     ( ( (pz < cent[kL][jL][iL].z) && (pz < cent[kR][jR][iR].z) ) || ( (pz > cent[kR][jR][iR].z) && (pz > cent[kL][jL][iL].z) ) )
    //
    // )
    // {
    //   char warning[512];
    //   sprintf(warning, "\n\n\npoint %lf %lf %lf is outside the trilinear interpolation box of %ld %ld %ld, %ld %ld %ld - min bound = %lf %lf %lf, max bound %lf %lf %lf, results may be incorrect!!\n\n\n", px, py, pz, kL, jL, iL, kR, jR, iR, cent[kL][jL][iL].x, cent[kL][jL][iL].y, cent[kL][jL][iL].z, cent[kR][jR][iR].x, cent[kR][jR][iR].y, cent[kR][jR][iR].z);
    //   warningInFunction("PointInterpolationWeights",  warning);
    // }

    // find directional gradients
    Cmpnts csiHat, etaHat, zetHat;
    PetscReal csiNor, etaNor, zetNor;

    csiHat.x = cent[kL][jL][iR].x - cent[kL][jL][iL].x;
    csiHat.y = cent[kL][jL][iR].y - cent[kL][jL][iL].y;
    csiHat.z = cent[kL][jL][iR].z - cent[kL][jL][iL].z;

    etaHat.x = cent[kL][jR][iL].x - cent[kL][jL][iL].x;
    etaHat.y = cent[kL][jR][iL].y - cent[kL][jL][iL].y;
    etaHat.z = cent[kL][jR][iL].z - cent[kL][jL][iL].z;

    zetHat.x = cent[kR][jL][iL].x - cent[kL][jL][iL].x;
    zetHat.y = cent[kR][jL][iL].y - cent[kL][jL][iL].y;
    zetHat.z = cent[kR][jL][iL].z - cent[kL][jL][iL].z;

    csiNor = std::sqrt(csiHat.x*csiHat.x + csiHat.y*csiHat.y + csiHat.z*csiHat.z);
    etaNor = std::sqrt(etaHat.x*etaHat.x + etaHat.y*etaHat.y + etaHat.z*etaHat.z);
    zetNor = std::sqrt(zetHat.x*zetHat.x + zetHat.y*zetHat.y + zetHat.z*zetHat.z);

    csiHat.x = csiHat.x / csiNor; csiHat.y = csiHat.y / csiNor; csiHat.z = csiHat.z / csiNor;
    etaHat.x = etaHat.x / etaNor; etaHat.y = etaHat.y / etaNor; etaHat.z = etaHat.z / etaNor;
    zetHat.x = zetHat.x / zetNor; zetHat.y = zetHat.y / zetNor; zetHat.z = zetHat.z / zetNor;

    // bound deltas from 0 to 1: this means that if the query point is outside of boundaries (we have
    // no data there), we project that point on the boundary and give the interpolated value at that
    // point as a result.

    PetscReal deltaCsi
    =
    PetscMax
    (
        PetscMin
        (
            (
                csiHat.x * (px - cent[kL][jL][iL].x) +
                csiHat.y * (py - cent[kL][jL][iL].y) +
                csiHat.z * (pz - cent[kL][jL][iL].z)

            ) / csiNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaEta
    =
    PetscMax
    (
        PetscMin
        (
            (
                etaHat.x * (px - cent[kL][jL][iL].x) +
                etaHat.y * (py - cent[kL][jL][iL].y) +
                etaHat.z * (pz - cent[kL][jL][iL].z)

            ) / etaNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaZet
    =
    PetscMax
    (
        PetscMin
        (
            (
                zetHat.x * (px - cent[kL][jL][iL].x) +
                zetHat.y * (py - cent[kL][jL][iL].y) +
                zetHat.z * (pz - cent[kL][jL][iL].z)

            ) / zetNor,
            1.0
        ),
        0.0
    );

    // find interpolation weights
    intWts[0] = 1 - deltaCsi - deltaEta - deltaZet + (deltaCsi * deltaEta) + (deltaEta * deltaZet) + (deltaZet * deltaCsi) - (deltaCsi * deltaEta * deltaZet),
    intWts[1] = deltaCsi - (deltaCsi * deltaEta) - (deltaZet * deltaCsi) + (deltaCsi * deltaEta * deltaZet),
    intWts[2] = deltaEta - (deltaCsi * deltaEta) - (deltaEta * deltaZet) + (deltaCsi * deltaEta * deltaZet),
    intWts[3] = (deltaCsi * deltaEta) - (deltaCsi * deltaEta * deltaZet),
    intWts[4] = deltaZet - (deltaEta * deltaZet) - (deltaZet * deltaCsi) + (deltaCsi * deltaEta * deltaZet),
    intWts[5] = (deltaZet * deltaCsi) - (deltaCsi * deltaEta * deltaZet),
    intWts[6] = (deltaEta * deltaZet) - (deltaCsi * deltaEta * deltaZet),
    intWts[7] = (deltaCsi * deltaEta * deltaZet);

    return;
}
//***************************************************************************************************************//
//! \brief Trilinear volume interpolation function for scalars to return interpolation cells
inline void PointInterpolationCells
(
    mesh_ *mesh,
    PetscReal px, PetscReal py, PetscReal pz,
    PetscInt ic, PetscInt jc, PetscInt kc,
    Cmpnts ***cent,
    PetscInt *intId
)
{
    // Given generic point coordinates px, py, pz and the indices of the closest
    // cell center, this functions finds the closest 8 cells to the given point which form
    // its surrounding box and interpolates the point value using the values at the 8
    // surrounding cell centers (tri-linear interpolation).
    // Note: cell indices passed to this function must be internal cells, not ghosts.
    // Throws an error otherwise.
    // Interpolation weights are bounded from 0 to 1 to handle query points lying outside
    // of the boundaries: this means that those points are projected on the boundaries and
    // the value there is returned (Prevents extrapolation ouside of boundaries).

    DMDALocalInfo    info = mesh->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt         gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt         gys  = info.gys, gye = info.gys + info.gym;
    PetscInt         gzs  = info.gzs, gze = info.gzs + info.gzm;

    PetscInt         i, j, k;

    PetscInt iL, iR;

    // do the search if surrounded by cells
    if(ic > 1 && ic < mx-2)
    {
        if (ic-1 < gxs)
        {
            iL = ic; iR = ic + 1;
        }
        else if(ic+1 >= gxe)
        {
            iL = ic - 1; iR = ic;
        }
        else
        {
            PetscReal dcsiRightSq
            =
            std::pow(cent[kc][jc][ic+1].x - px,2) +
            std::pow(cent[kc][jc][ic+1].y - py,2) +
            std::pow(cent[kc][jc][ic+1].z - pz,2);

            PetscReal dcsiLeftSq
            =
            std::pow(px - cent[kc][jc][ic-1].x,2) +
            std::pow(py - cent[kc][jc][ic-1].y,2) +
            std::pow(pz - cent[kc][jc][ic-1].z,2);

            if(dcsiRightSq <= dcsiLeftSq)
            {
                iR = ic + 1; iL = ic;
            }
            else
            {
                iR = ic; iL = ic-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(ic == 1)
    {
        iR = ic + 1; iL = ic;
    }
    // on boundary, already know where to pick data
    else if(ic == mx-2)
    {
        iR = ic; iL = ic-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "i-index out of bounds: ended up in the ghosts when looking for interpolation (i = %ld)\n", ic);
        fatalErrorInFunction("scalarPointInterpolationCells",  error);
    }

    PetscInt jL, jR;
    // do the search if surrounded by cells
    if(jc > 1 && jc < my-2)
    {
        // check processor bounds
        if (jc-1 < gys)
        {
            jL = jc; jR = jc + 1;
        }
        else if(jc+1 >= gye)
        {
            jL = jc - 1; jR = jc;
        }
        else
        {
            PetscReal detaRightSq
            =
            std::pow(cent[kc][jc+1][ic].x - px,2) +
            std::pow(cent[kc][jc+1][ic].y - py,2) +
            std::pow(cent[kc][jc+1][ic].z - pz,2);

            PetscReal detaLeftSq
            =
            std::pow(px - cent[kc][jc-1][ic].x,2) +
            std::pow(py - cent[kc][jc-1][ic].y,2) +
            std::pow(pz - cent[kc][jc-1][ic].z,2);

            if(detaRightSq <= detaLeftSq)
            {
                jR = jc + 1; jL = jc;
            }
            else
            {
                jR = jc; jL = jc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(jc == 1)
    {
        jR = jc + 1; jL = jc;
    }
    // on boundary, already know where to pick data
    else if(jc == my-2)
    {
        jR = jc; jL = jc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "j-index out of bounds: ended up in the ghosts when looking for interpolation (j = %ld)\n", jc);
        fatalErrorInFunction("scalarPointInterpolationCells",  error);
    }

    PetscInt kL, kR;
    // do the search if surrounded by cells
    if(kc > 1 && kc < mz-2)
    {
        // check processor bounds
        if (kc-1 < gzs)
        {
            kL = kc; kR = kc + 1;
        }
        else if(kc+1 >= gze)
        {
            kL = kc - 1; kR = kc;
        }
        else
        {
            PetscReal dzetRightSq
            =
            std::pow(cent[kc+1][jc][ic].x - px,2) +
            std::pow(cent[kc+1][jc][ic].y - py,2) +
            std::pow(cent[kc+1][jc][ic].z - pz,2);

            PetscReal dzetLeftSq
            =
            std::pow(px - cent[kc-1][jc][ic].x,2) +
            std::pow(py - cent[kc-1][jc][ic].y,2) +
            std::pow(pz - cent[kc-1][jc][ic].z,2);

            if(dzetRightSq <= dzetLeftSq)
            {
                kR = kc + 1; kL = kc;
            }
            else
            {
                kR = kc; kL = kc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(kc == 1)
    {
        kR = kc + 1; kL = kc;
    }
    // on boundary, already know where to pick data
    else if(kc == mz-2)
    {
        kR = kc; kL = kc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "k-index out of bounds: ended up in the ghosts when looking for interpolation (k = %ld)\n", kc);
        fatalErrorInFunction("scalarPointInterpolationCells",  error);
    }

    intId[0] = kL;
    intId[1] = kR;
    intId[2] = jL;
    intId[3] = jR;
    intId[4] = iL;
    intId[5] = iR;

    return;
}
//***************************************************************************************************************//
//! \brief Trilinear volume interpolation function for vectors
inline void vectorPointLocalVolumeInterpolation
(
    mesh_ *mesh,
    PetscReal px, PetscReal py, PetscReal pz,
    PetscInt ic, PetscInt jc, PetscInt kc,
    Cmpnts ***cent,
    Cmpnts ***v,
    Cmpnts &result
)
{
    // Given generic point coordinates px, py, pz and the indices of the closest
    // cell center, this functions finds the closest 8 cells to the given point which form
    // its surrounding box and interpolates the point value using the values at the 8
    // surrounding cell centers (tri-linear interpolation).
    // Note: cell indices passed to this function must be internal cells, not ghosts.
    // Throws an error otherwise.
    // Interpolation weights are bounded from 0 to 1 to handle query points lying outside
    // of the boundaries: this means that those points are projected on the boundaries and
    // the value there is returned (Prevents extrapolation ouside of boundaries).

    DMDALocalInfo    info = mesh->info;
    PetscInt         xs   = info.xs, xe = info.xs + info.xm;
    PetscInt         ys   = info.ys, ye = info.ys + info.ym;
    PetscInt         zs   = info.zs, ze = info.zs + info.zm;
    PetscInt         mx   = info.mx, my = info.my, mz = info.mz;
    PetscInt         gxs  = info.gxs, gxe = info.gxs + info.gxm;
    PetscInt         gys  = info.gys, gye = info.gys + info.gym;
    PetscInt         gzs  = info.gzs, gze = info.gzs + info.gzm;

    PetscInt         i, j, k;

    PetscInt iL, iR;
    // do the search if surrounded by cells
    if(ic > 1 && ic < mx-2)
    {
        if (ic-1 < gxs)
        {
            iL = ic; iR = ic + 1;
        }
        else if(ic+1 >= gxe)
        {
            iL = ic - 1; iR = ic;
        }
        else
        {
            PetscReal dcsiRightSq
            =
            std::pow(cent[kc][jc][ic+1].x - px,2) +
            std::pow(cent[kc][jc][ic+1].y - py,2) +
            std::pow(cent[kc][jc][ic+1].z - pz,2);

            PetscReal dcsiLeftSq
            =
            std::pow(px - cent[kc][jc][ic-1].x,2) +
            std::pow(py - cent[kc][jc][ic-1].y,2) +
            std::pow(pz - cent[kc][jc][ic-1].z,2);

            if(dcsiRightSq <= dcsiLeftSq)
            {
                iR = ic + 1; iL = ic;
            }
            else
            {
                iR = ic; iL = ic-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(ic == 1)
    {
        iR = ic + 1; iL = ic;
    }
    // on boundary, already know where to pick data
    else if(ic == mx-2)
    {
        iR = ic; iL = ic-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "i-index out of bounds: ended up in the ghosts when looking for interpolation data (i = %ld)\n", ic);
        fatalErrorInFunction("vectorPointLocalVolumeInterpolation",  error);
    }

    PetscInt jL, jR;
    // do the search if surrounded by cells
    if(jc > 1 && jc < my-2)
    {
        // check processor bounds
        if (jc-1 < gys)
        {
            jL = jc; jR = jc + 1;
        }
        else if(jc+1 >= gye)
        {
            jL = jc - 1; jR = jc;
        }
        else
        {
            PetscReal detaRightSq
            =
            std::pow(cent[kc][jc+1][ic].x - px,2) +
            std::pow(cent[kc][jc+1][ic].y - py,2) +
            std::pow(cent[kc][jc+1][ic].z - pz,2);

            PetscReal detaLeftSq
            =
            std::pow(px - cent[kc][jc-1][ic].x,2) +
            std::pow(py - cent[kc][jc-1][ic].y,2) +
            std::pow(pz - cent[kc][jc-1][ic].z,2);

            if(detaRightSq <= detaLeftSq)
            {
                jR = jc + 1; jL = jc;
            }
            else
            {
                jR = jc; jL = jc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(jc == 1)
    {
        jR = jc + 1; jL = jc;
    }
    // on boundary, already know where to pick data
    else if(jc == my-2)
    {
        jR = jc; jL = jc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "j-index out of bounds: ended up in the ghosts when looking for interpolation data (j = %ld)\n", jc);
        fatalErrorInFunction("vectorPointLocalVolumeInterpolation",  error);
    }

    PetscInt kL, kR;
    // do the search if surrounded by cells
    if(kc > 1 && kc < mz-2)
    {
        // check processor bounds
        if (kc-1 < gzs)
        {
            kL = kc; kR = kc + 1;
        }
        else if(kc+1 >= gze)
        {
            kL = kc - 1; kR = kc;
        }
        else
        {
            PetscReal dzetRightSq
            =
            std::pow(cent[kc+1][jc][ic].x - px,2) +
            std::pow(cent[kc+1][jc][ic].y - py,2) +
            std::pow(cent[kc+1][jc][ic].z - pz,2);

            PetscReal dzetLeftSq
            =
            std::pow(px - cent[kc-1][jc][ic].x,2) +
            std::pow(py - cent[kc-1][jc][ic].y,2) +
            std::pow(pz - cent[kc-1][jc][ic].z,2);

            if(dzetRightSq <= dzetLeftSq)
            {
                kR = kc + 1; kL = kc;
            }
            else
            {
                kR = kc; kL = kc-1;
            }
        }

    }
    // on boundary, already know where to pick data
    else if(kc == 1)
    {
        kR = kc + 1; kL = kc;
    }
    // on boundary, already know where to pick data
    else if(kc == mz-2)
    {
        kR = kc; kL = kc-1;
    }
    // outside the boundary: the Vec V might not contain useful data.
    else
    {
        char error[512];
        sprintf(error, "k-index out of bounds: ended up in the ghosts when looking for interpolation data (k = %ld)\n", kc);
        fatalErrorInFunction("vectorPointLocalVolumeInterpolation",  error);
    }

    /*
    // Debugging information (only i and j)
    PetscPrintf(PETSC_COMM_WORLD, "Neighboring cells\n");
    PetscPrintf(PETSC_COMM_WORLD, "(iL,jL) idxs: %d %d, C = %.2f %.2f\n", iL, jL, cent[kL][jL][iL].y, cent[kL][jL][iL].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iL,jR) idxs: %d %d, C = %.2f %.2f\n", iL, jR, cent[kL][jR][iL].y, cent[kL][jR][iL].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iR,jL) idxs: %d %d, C = %.2f %.2f\n", iR, jL, cent[kL][jL][iR].y, cent[kL][jL][iR].z);
    PetscPrintf(PETSC_COMM_WORLD, "(iR,jR) idxs: %d %d, C = %.2f %.2f\n", iR, jR, cent[kL][jR][iR].y, cent[kL][jR][iR].z);
    PetscPrintf(PETSC_COMM_WORLD, "(ic,jc) idxs: %d %d\n", ic, jc);
    PetscPrintf(PETSC_COMM_WORLD, "P = %.2f %.2f\n\n", py, pz);
    */

    // find directional gradients
    Cmpnts csiHat, etaHat, zetHat;
    PetscReal csiNor, etaNor, zetNor;

    csiHat.x = cent[kL][jL][iR].x - cent[kL][jL][iL].x;
    csiHat.y = cent[kL][jL][iR].y - cent[kL][jL][iL].y;
    csiHat.z = cent[kL][jL][iR].z - cent[kL][jL][iL].z;

    etaHat.x = cent[kL][jR][iL].x - cent[kL][jL][iL].x;
    etaHat.y = cent[kL][jR][iL].y - cent[kL][jL][iL].y;
    etaHat.z = cent[kL][jR][iL].z - cent[kL][jL][iL].z;

    zetHat.x = cent[kR][jL][iL].x - cent[kL][jL][iL].x;
    zetHat.y = cent[kR][jL][iL].y - cent[kL][jL][iL].y;
    zetHat.z = cent[kR][jL][iL].z - cent[kL][jL][iL].z;

    csiNor = std::sqrt(csiHat.x*csiHat.x + csiHat.y*csiHat.y + csiHat.z*csiHat.z);
    etaNor = std::sqrt(etaHat.x*etaHat.x + etaHat.y*etaHat.y + etaHat.z*etaHat.z);
    zetNor = std::sqrt(zetHat.x*zetHat.x + zetHat.y*zetHat.y + zetHat.z*zetHat.z);

    csiHat.x = csiHat.x / csiNor; csiHat.y = csiHat.y / csiNor; csiHat.z = csiHat.z / csiNor;
    etaHat.x = etaHat.x / etaNor; etaHat.y = etaHat.y / etaNor; etaHat.z = etaHat.z / etaNor;
    zetHat.x = zetHat.x / zetNor; zetHat.y = zetHat.y / zetNor; zetHat.z = zetHat.z / zetNor;

    // bound deltas from 0 to 1: this means that if the query point is outside of boundaries (we have
    // no data there), we project that point on the boundary and give the interpolated value at that
    // point as a result.

    PetscReal deltaCsi
    =
    PetscMax
    (
        PetscMin
        (
            (
                csiHat.x * (px - cent[kL][jL][iL].x) +
                csiHat.y * (py - cent[kL][jL][iL].y) +
                csiHat.z * (pz - cent[kL][jL][iL].z)

            ) / csiNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaEta
    =
    PetscMax
    (
        PetscMin
        (
            (
                etaHat.x * (px - cent[kL][jL][iL].x) +
                etaHat.y * (py - cent[kL][jL][iL].y) +
                etaHat.z * (pz - cent[kL][jL][iL].z)

            ) / etaNor,
            1.0
        ),
        0.0
    );

    PetscReal deltaZet
    =
    PetscMax
    (
        PetscMin
        (
            (
                zetHat.x * (px - cent[kL][jL][iL].x) +
                zetHat.y * (py - cent[kL][jL][iL].y) +
                zetHat.z * (pz - cent[kL][jL][iL].z)

            ) / zetNor,
            1.0
        ),
        0.0
    );

    // find interpolation weights
    Cmpnts v000 = nSet(v[kL][jL][iL]),
           v100 = nSet(v[kL][jL][iR]),
           v110 = nSet(v[kL][jR][iR]),
           v010 = nSet(v[kL][jR][iL]),
           v001 = nSet(v[kR][jL][iL]),
           v101 = nSet(v[kR][jL][iR]),
           v111 = nSet(v[kR][jR][iR]),
           v011 = nSet(v[kR][jR][iL]);

    Cmpnts c0, c1, c2, c3, c4, c5, c6, c7;

    c0.x = v000.x,
    c1.x = v100.x - v000.x,
    c2.x = v010.x - v000.x,
    c3.x = v001.x - v000.x,
    c4.x = v110.x - v010.x - v100.x + v000.x,
    c5.x = v011.x - v001.x - v010.x + v000.x,
    c6.x = v101.x - v001.x - v100.x + v000.x,
    c7.x = v111.x - v011.x - v101.x - v110.x + v100.x + v001.x + v010.x - v000.x;

    c0.y = v000.y,
    c1.y = v100.y - v000.y,
    c2.y = v010.y - v000.y,
    c3.y = v001.y - v000.y,
    c4.y = v110.y - v010.y - v100.y + v000.y,
    c5.y = v011.y - v001.y - v010.y + v000.y,
    c6.y = v101.y - v001.y - v100.y + v000.y,
    c7.y = v111.y - v011.y - v101.y - v110.y + v100.y + v001.y + v010.y - v000.y;

    c0.z = v000.z,
    c1.z = v100.z - v000.z,
    c2.z = v010.z - v000.z,
    c3.z = v001.z - v000.z,
    c4.z = v110.z - v010.z - v100.z + v000.z,
    c5.z = v011.z - v001.z - v010.z + v000.z,
    c6.z = v101.z - v001.z - v100.z + v000.z,
    c7.z = v111.z - v011.z - v101.z - v110.z + v100.z + v001.z + v010.z - v000.z;


    // interpolate
    result.x
    =
    c0.x +
    c1.x * deltaCsi +
    c2.x * deltaEta +
    c3.x * deltaZet +
    c4.x * deltaCsi * deltaEta +
    c5.x * deltaEta * deltaZet +
    c6.x * deltaZet * deltaCsi +
    c7.x * deltaCsi * deltaEta * deltaZet;

    result.y
    =
    c0.y +
    c1.y * deltaCsi +
    c2.y * deltaEta +
    c3.y * deltaZet +
    c4.y * deltaCsi * deltaEta +
    c5.y * deltaEta * deltaZet +
    c6.y * deltaZet * deltaCsi +
    c7.y * deltaCsi * deltaEta * deltaZet;

    result.z
    =
    c0.z +
    c1.z * deltaCsi +
    c2.z * deltaEta +
    c3.z * deltaZet +
    c4.z * deltaCsi * deltaEta +
    c5.z * deltaEta * deltaZet +
    c6.z * deltaZet * deltaCsi +
    c7.z * deltaCsi * deltaEta * deltaZet;

    return;
}

//***************************************************************************************************************//

// Dynamic array allocation for increasing array size: linked list
inline void initlist(list *ilist)
{

  ilist->head = PETSC_NULL;

  return;
}

//***************************************************************************************************************//

// initialize a cell list to null pointer
inline void initCellList(cellList *ilist)
{

  ilist->head = PETSC_NULL;
  return;
}

//***************************************************************************************************************//

// insert a node into a list and check if previously added
inline bool insertnode(list *ilist, PetscInt Node)
{

    node *_new;
    node *current;

    current = ilist->head;

    PetscBool Exist = PETSC_FALSE;

    while(current)
    {
        if (Node == current->Node)
        {
            Exist = PETSC_TRUE;
        }

        if (Exist) break;
        current = current->next;
    }

    if (!Exist)
    {
        PetscMalloc(sizeof(node), &_new);
        _new->next = ilist->head;
        _new->Node = Node;
        ilist->head = _new;
    }

    if(Exist)
        return false;
    else
        return true;
}

//***************************************************************************************************************//

// insert a node into a list
inline void insertnode1(list *ilist, PetscInt Node)
{
    node *_new;

    PetscMalloc(sizeof(node), &_new);

    _new->next = ilist->head;

    _new->Node = Node;

    ilist->head = _new;

    return;
}

//***************************************************************************************************************//

// insert a node into a list and check for previous existing node
inline bool insertCellNode(cellList *ilist, cellIds Node)
{

  cellNode *_new;
  cellNode *current;
  current = ilist->head;

  PetscBool Exist = PETSC_FALSE;
  while(current) {
    if (Node.i == current->Node.i
      && Node.j == current->Node.j
      && Node.k == current->Node.k)
       {
         Exist = PETSC_TRUE;
       }

    if (Exist) break;
    current = current->next;
  }

  if (!Exist) {
    PetscMalloc(sizeof(cellNode), &_new);
    _new->next = ilist->head;

    _new->Node.i = Node.i;
    _new->Node.j = Node.j;
    _new->Node.k = Node.k;

    ilist->head = _new;
  }

  if(Exist)
    return false;
  else
    return true;
}

//***************************************************************************************************************//

// insert a node into a list and dont check if it previously exists
inline void insertCellNode1(cellList *ilist, cellIds Node)
{
  cellNode *_new;

  PetscMalloc(sizeof(cellNode), &_new);

  _new->next = ilist->head;

  _new->Node.i = Node.i;
  _new->Node.j = Node.j;
  _new->Node.k = Node.k;

  ilist->head = _new;

  return;
}

// detroy a list
inline void destroy(list *ilist)
{

  node *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }

  return;
}

//***************************************************************************************************************//

// detroy a cell list
inline void destroyCellList(cellList *ilist)
{

  cellNode *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }

  return;
}

//***************************************************************************************************************//
// inverse matrix functions

//***************************************************************************************************************//

inline void reorderMatrix(PetscReal **A, PetscReal **B, PetscInt N, PetscInt pivot)
{
    PetscInt  maxrow=pivot;
    PetscInt  row,col;
    PetscReal temp;

    for (row=pivot+1;row<N;row++)
    {
        if ( fabs(A[row][pivot]) > fabs(A[maxrow][pivot]) ) maxrow = row;
    }
        

   if (fabs(A[maxrow][pivot])<1.e-11)
   {
        char error[512];
        sprintf(error, "matrix inverse reorder algorithm failed\n");
        fatalErrorInFunction("reorderMatrix",  error);
   }
       
    if (maxrow != pivot)
    {
        for (col=0; col<N; col++)
        {
            temp           = A[maxrow][col];
            A[maxrow][col] = A[pivot][col];
            A[pivot][col]  = temp;

            temp           = B[maxrow][col];
            B[maxrow][col] = B[pivot][col];
            B[pivot][col]  = temp;
        }
    }

    return;
}

//***************************************************************************************************************//

inline void elimJordanAlg(PetscReal **A, PetscReal **B, PetscInt N, PetscInt pivot)
{
    
    PetscInt row,col;

    for (col=pivot+1; col<N; col++)
    {
        A[pivot][col] = A[pivot][col]/A[pivot][pivot];
    } 

    for (col=0; col<N; col++)
    {
        B[pivot][col] = B[pivot][col]/A[pivot][pivot];
    }
        
    A[pivot][pivot]=1.;

    for (row=0; row<N; row++)
    {
        if (row == pivot) continue;

        for (col=pivot+1; col<N; col++) A[row][col] -= A[row][pivot]*A[pivot][col];
        for (col=0; col<N; col++)       B[row][col] -= A[row][pivot]*B[pivot][col];
        A[row][pivot]=0.;
    }

    return;
}

//***************************************************************************************************************//
// inverse of generic matrix of any order 
inline void inverseMatrix(PetscReal **A, PetscReal **inv_A, PetscInt size)
{
    PetscInt i,j,k;

    for (i=0; i<size; i++)
    {
        for (j=0; j<size; j++)
        {
            inv_A[i][j]=0.;
        }
    }
        
    for (i=0; i<size; i++)
    {
        inv_A[i][i]=1.0;
    } 

    PetscInt pivot=0;

    while (pivot < size)
    {
        reorderMatrix(A, inv_A, size, pivot);
        elimJordanAlg(A, inv_A, size, pivot);

        pivot++;
    }

    return;
}

//***************************************************************************************************************//
//invert a 4 by 4 matrix  -- faster
inline void inv_4by4(PetscReal **A, PetscReal **inv_A, PetscInt size)
{

    if( size!= 4 )
    {
        char error[512];
        sprintf(error, "expected array of size 4\n");
        fatalErrorInFunction("inv_4by4",  error);
    }

    PetscReal A2323 = A[2][2] * A[3][3] - A[2][3] * A[3][2] ;
    PetscReal A1323 = A[2][1] * A[3][3] - A[2][3] * A[3][1] ;
    PetscReal A1223 = A[2][1] * A[3][2] - A[2][2] * A[3][1] ;
    PetscReal A0323 = A[2][0] * A[3][3] - A[2][3] * A[3][0] ;
    PetscReal A0223 = A[2][0] * A[3][2] - A[2][2] * A[3][0] ;
    PetscReal A0123 = A[2][0] * A[3][1] - A[2][1] * A[3][0] ;
    PetscReal A2313 = A[1][2] * A[3][3] - A[1][3] * A[3][2] ;
    PetscReal A1313 = A[1][1] * A[3][3] - A[1][3] * A[3][1] ;
    PetscReal A1213 = A[1][1] * A[3][2] - A[1][2] * A[3][1] ;
    PetscReal A2312 = A[1][2] * A[2][3] - A[1][3] * A[2][2] ;
    PetscReal A1312 = A[1][1] * A[2][3] - A[1][3] * A[2][1] ;
    PetscReal A1212 = A[1][1] * A[2][2] - A[1][2] * A[2][1] ;
    PetscReal A0313 = A[1][0] * A[3][3] - A[1][3] * A[3][0] ;
    PetscReal A0213 = A[1][0] * A[3][2] - A[1][2] * A[3][0] ;
    PetscReal A0312 = A[1][0] * A[2][3] - A[1][3] * A[2][0] ;
    PetscReal A0212 = A[1][0] * A[2][2] - A[1][2] * A[2][0] ;
    PetscReal A0113 = A[1][0] * A[3][1] - A[1][1] * A[3][0] ;
    PetscReal A0112 = A[1][0] * A[2][1] - A[1][1] * A[2][0] ;

    PetscReal det = A[0][0] * ( A[1][1] * A2323 - A[1][2] * A1323 + A[1][3] * A1223 )
        - A[0][1] * ( A[1][0] * A2323 - A[1][2] * A0323 + A[1][3] * A0223 )
        + A[0][2] * ( A[1][0] * A1323 - A[1][1] * A0323 + A[1][3] * A0123 )
        - A[0][3] * ( A[1][0] * A1223 - A[1][1] * A0223 + A[1][2] * A0123 ) ;

    det = 1 / det;

    inv_A[0][0] = det *   ( A[1][1] * A2323 - A[1][2] * A1323 + A[1][3] * A1223 );
    inv_A[0][1] = det * - ( A[0][1] * A2323 - A[0][2] * A1323 + A[0][3] * A1223 );
    inv_A[0][2] = det *   ( A[0][1] * A2313 - A[0][2] * A1313 + A[0][3] * A1213 );
    inv_A[0][3] = det * - ( A[0][1] * A2312 - A[0][2] * A1312 + A[0][3] * A1212 );
    inv_A[1][0] = det * - ( A[1][0] * A2323 - A[1][2] * A0323 + A[1][3] * A0223 );
    inv_A[1][1] = det *   ( A[0][0] * A2323 - A[0][2] * A0323 + A[0][3] * A0223 );
    inv_A[1][2] = det * - ( A[0][0] * A2313 - A[0][2] * A0313 + A[0][3] * A0213 );
    inv_A[1][3] = det *   ( A[0][0] * A2312 - A[0][2] * A0312 + A[0][3] * A0212 );
    inv_A[2][0] = det *   ( A[1][0] * A1323 - A[1][1] * A0323 + A[1][3] * A0123 );
    inv_A[2][1] = det * - ( A[0][0] * A1323 - A[0][1] * A0323 + A[0][3] * A0123 );
    inv_A[2][2] = det *   ( A[0][0] * A1313 - A[0][1] * A0313 + A[0][3] * A0113 );
    inv_A[2][3] = det * - ( A[0][0] * A1312 - A[0][1] * A0312 + A[0][3] * A0112 );
    inv_A[3][0] = det * - ( A[1][0] * A1223 - A[1][1] * A0223 + A[1][2] * A0123 );
    inv_A[3][1] = det *   ( A[0][0] * A1223 - A[0][1] * A0223 + A[0][2] * A0123 );
    inv_A[3][2] = det * - ( A[0][0] * A1213 - A[0][1] * A0213 + A[0][2] * A0113 );
    inv_A[3][3] = det *   ( A[0][0] * A1212 - A[0][1] * A0212 + A[0][2] * A0112 );

    return;
}

//***************************************************************************************************************//
//invert a 3 by 3 matrix
inline void inv_3by3(PetscReal **A, PetscReal **inv_A, PetscInt size)
{
    if( size!= 3 )
    {
        char error[512];
        sprintf(error, "expected array of size 3\n");
        fatalErrorInFunction("inv_3by3",  error);
    }

    PetscReal det = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
                    A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
                    A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);

    det = 1 / det;

    // Compute elements of the inverse matrix
    inv_A[0][0] = (A[1][1] * A[2][2] - A[2][1] * A[1][2]) * det;
    inv_A[0][1] = (A[0][2] * A[2][1] - A[0][1] * A[2][2]) * det;
    inv_A[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) * det;
    inv_A[1][0] = (A[1][2] * A[2][0] - A[1][0] * A[2][2]) * det;
    inv_A[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) * det;
    inv_A[1][2] = (A[1][0] * A[0][2] - A[0][0] * A[1][2]) * det;
    inv_A[2][0] = (A[1][0] * A[2][1] - A[2][0] * A[1][1]) * det;
    inv_A[2][1] = (A[2][0] * A[0][1] - A[0][0] * A[2][1]) * det;
    inv_A[2][2] = (A[0][0] * A[1][1] - A[1][0] * A[0][1]) * det;

    return;
}

//***************************************************************************************************************//

inline void reorder_20(PetscReal A[][20], PetscReal B[][20],PetscInt N,PetscInt pivot)
{
    PetscInt maxrow=pivot;
    PetscInt row,col;
    PetscReal temp;
    for (row=pivot+1;row<N;row++)
        if (fabs(A[row][pivot])>fabs(A[maxrow][pivot])) maxrow=row;

//    if (fabs(A[maxrow][pivot])<1.e-11)
//        exit(0);
    if (maxrow!=pivot)
        for (col=0;col<N;col++){
            temp= A[maxrow][col];
            A[maxrow][col]=A[pivot][col];
            A[pivot][col]=temp;
            temp= B[maxrow][col];
            B[maxrow][col]=B[pivot][col];
            B[pivot][col]=temp;
        }
    return;
}

//***************************************************************************************************************//

inline void elim_jordan_20(PetscReal A[][20], PetscReal B[][20],PetscInt N,PetscInt pivot){
    PetscInt row,col;
    for (col=pivot+1;col<N;col++) A[pivot][col]=A[pivot][col]/A[pivot][pivot];
    for (col=0;col<N;col++)
        B[pivot][col]=B[pivot][col]/A[pivot][pivot];
    A[pivot][pivot]=1.;
    for (row=0;row<N;row++){
        if (row==pivot) continue;
        for (col=pivot+1;col<N;col++) A[row][col] -= A[row][pivot]*A[pivot][col];
        for (col=0;col<N;col++)       B[row][col] -= A[row][pivot]*B[pivot][col];
        A[row][pivot]=0.;
    }

    return;
}

//***************************************************************************************************************//

inline void inv_20(PetscReal A[20][20], PetscReal inv_A[][20],PetscInt N){
    PetscInt i,j,k;
    for (i=0;i<N;i++)
        for (j=0;j<N;j++)
            inv_A[i][j]=0.;
    for (i=0;i<N;i++) inv_A[i][i]=1.;
    PetscInt pivot=0;
    while (pivot<N){
        reorder_20(A, inv_A, N, pivot);
        elim_jordan_20(A,inv_A,N,pivot);
        pivot++;
    }

    return;
}

//***************************************************************************************************************//

inline void reorder_10(PetscReal A[][10], PetscReal B[][10],PetscInt N,PetscInt pivot){
    PetscInt maxrow=pivot;
    PetscInt row,col;
    PetscReal temp;
    for (row=pivot+1;row<N;row++)
        if (fabs(A[row][pivot])>fabs(A[maxrow][pivot])) maxrow=row;

//    if (fabs(A[maxrow][pivot])<1.e-11)
//        exit(0);
    if (maxrow!=pivot)
        for (col=0;col<N;col++){
            temp= A[maxrow][col];
            A[maxrow][col]=A[pivot][col];
            A[pivot][col]=temp;
            temp= B[maxrow][col];
            B[maxrow][col]=B[pivot][col];
            B[pivot][col]=temp;
        }
    return;
}

//***************************************************************************************************************//

inline void elim_jordan_10(PetscReal A[][10], PetscReal B[][10],PetscInt N,PetscInt pivot){
    PetscInt row,col;
    for (col=pivot+1;col<N;col++) A[pivot][col]=A[pivot][col]/A[pivot][pivot];
    for (col=0;col<N;col++)
        B[pivot][col]=B[pivot][col]/A[pivot][pivot];
    A[pivot][pivot]=1.;
    for (row=0;row<N;row++){
        if (row==pivot) continue;
        for (col=pivot+1;col<N;col++) A[row][col] -= A[row][pivot]*A[pivot][col];
        for (col=0;col<N;col++)       B[row][col] -= A[row][pivot]*B[pivot][col];
        A[row][pivot]=0.;
    }

    return;
}

//***************************************************************************************************************//

inline void inv_10(PetscReal A[10][10], PetscReal inv_A[][10],PetscInt N){
    PetscInt i,j,k;
    for (i=0;i<N;i++)
        for (j=0;j<N;j++)
            inv_A[i][j]=0.;
    for (i=0;i<N;i++) inv_A[i][i]=1.;
    PetscInt pivot=0;
    while (pivot<N){
        reorder_10(A, inv_A, N, pivot);
        elim_jordan_10(A,inv_A,N,pivot);
        pivot++;
    }

    return;
}

//***************************************************************************************************************//

inline void mult_mats3_lin(PetscReal **inv_A,PetscReal **B,PetscInt nsupport,PetscReal *PHI){

  for (PetscInt i=0;i<nsupport;i++) {
    PHI[i]=0.;
    PetscReal temp=0.;
    for (PetscInt k=0;k<4;k++) temp += inv_A[0][k]*B[k][i];
    PHI[i] = temp;
  }

  return;
}


#endif
