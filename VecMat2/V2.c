/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/*****************************/
/* single vector with scalar */
/*****************************/

/* unary vector */

/* a[] := 0; return a[] */
double *V2zero (double a[2])
{
    a[0] = 0;
    a[1] = 0;
    return(a);
} /* end V2zero() */


/* generate a unit vector: a[i] := 1, a[j!=i] := 0; return a[] */
double *V2unit (double a[2], int i)
{
    a[0] = 0;
    a[1] = 0;
    a[i] = 1;
    return(a);
} /* end V2unit() */


/* generate a[] with 2 independent random components on (0.,1.); return a[] */
double *V2frandom (double a[2])
{
    a[0] = Frandom();
    a[1] = Frandom();
    return(a);
} /* end V2frandom() */


/* generate a unit vector of spherically uniform random orientation */
double *V2randomunit (double a[2])
{
    double length;
    while(1)
    {
        a[0] = FRANDOM();
        a[1] = FRANDOM();
        length = V2LENGTH(a);
        if ( (length > SMALL) && (length < 0.5 - SMALL) )
        { /* ensure numerical accuracy */
            a[0] /= length;
            a[1] /= length;
            return(a);
        }
    }
} /* end V2randomunit() */


/* b[] := a[]; then return b[] */
double *V2eqv (double a[2], double b[2])
{
    b[0] = a[0];
    b[1] = a[1];
    return(b);
} /* end V2eqv() */


/* b[] := -a[]; then return b[] */
double *V2neg (double a[2], double b[2])
{
    b[0] = -a[0];
    b[1] = -a[1];
    return(b);
} /* end V2neg() */


/* a[] := -a[]; then return a[] */
double *V2Neg (double a[2])
{
    a[0] = -a[0];
    a[1] = -a[1];
    return(a);
} /* end V2Neg() */


/* change a[] to its own image in [0,1)^2; then return a[] */
double *V2Trim (double a[2])
{
    a[0] = TRIM(a[0]);
    a[1] = TRIM(a[1]);
    return(a);
} /* end V2Trim() */


/* make b[] an image of a[] in [0,1)^2; then return b[] */
double *V2trim (double a[2], double b[2])
{
    b[0] = TRIM(a[0]);
    b[1] = TRIM(a[1]);
    return(b);
} /* end V2trim() */


/* make b[] image of a[]'s in [-0.5,0.5)^2; then return b[] */
double *V2image (double a[2], double b[2])
{
    b[0] = IMAGE(a[0]);
    b[1] = IMAGE(a[1]);
    return(b);
} /* end V2image() */


/* change a[] to its own image in [-0.5,0.5)^2; then return a[] */
double *V2Image (double a[2])
{
    a[0] = IMAGE(a[0]);
    a[1] = IMAGE(a[1]);
    return(a);
} /* end V2Image() */


/* scalar & vector */


/* b[] := multiplier * a[]; then return b[] */
double *V2mul (double multiplier, double a[2], double b[2])
{
    b[0] = multiplier * a[0];
    b[1] = multiplier * a[1];
    return(b);
} /* end V2mul() */


/* a[] := multiplier * a[]; then return a[] */
double *V2Mul (double multiplier, double a[2])
{
    a[0] *= multiplier;
    a[1] *= multiplier;
    return(a);
} /* end V2Mul() */


/* b[] := a[] / divisor; then return b[] */
double *V2div (double a[2], double divisor, double b[2])
{
    if (divisor == 0.) pe ("V2div: divisor = %e\n", divisor);
    b[0] = a[0] / divisor;
    b[1] = a[1] / divisor;
    return(b);
} /* end V2div() */


/* a[] := a[] / divisor; then return a[] */
double *V2Div (double a[2], double divisor)
{
    if (divisor == 0.) pe("V2DIV: divisor = %e\n", divisor);
    a[0] /= divisor;
    a[1] /= divisor;
    return(a);
} /* end V2DIV() */


/* length of a[] */
double V2length (double a[2])
{
    return (sqrt(a[0]*a[0]+a[1]*a[1]));
} /* end V2length() */


/* squared length of a[] */
double V2length2 (double a[2])
{
    return (a[0]*a[0]+a[1]*a[1]);
} /* end V2length2() */


/* arbitrary exponent-th norm of a[] */
double V2norm (double a[2], double exponent)
{
    if (exponent <= 0.)
        pe ("V2norm: norm exponent = %lf <= 0 is illegal\n", exponent);
    else if (exponent >= EXPONENT_INFINITY)
        return (MAX(fabs(a[0]),fabs(a[1])));
    else
        return (pow(pow(fabs(a[0]),exponent)+pow(fabs(a[1]),exponent),
                    1./exponent));
    return (0.);
} /* end V2norm() */


/* normalize a[] to unit length; then return a[] */
double *V2normalize (double a[2])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]);
    if (length == 0.)
    {
        pr ("warning: V2normalize: length = 0,\n"
            "will not normalize input vector\n");
        return (a);
    }
    a[0] /= length;
    a[1] /= length;
    return (a);
} /* end V2normalize() */


/* normalize a[] to unit length; then return its original length */
double V2Normalize (double a[2])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]);
    if (length == 0.)
    {
        pr ("warning: V2NORMALIZE: length = 0,"
            "will not normalize input vector\n");
        return (length);
    }
    a[0] /= length;
    a[1] /= length;
    return (length);
} /* end V2Normalize() */


/* b[] := a[] / |a[]|, return b[] */
double *V2direction (double a[2], double b[2])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]);
    if (length == 0.)
    {
        pr ("warning: V2direction: length = 0,\n"
            "will not normalize input vector\n");
        return (b);
    }
    b[0] = a[0] / length;
    b[1] = a[1] / length;
    return (b);
} /* end V2direction() */


/* b[] := a[] / |a[]|, return |a[]| */
double V2Direction (double a[2], double b[2])
{
    double length = sqrt(a[0]*a[0]+a[1]*a[1]);
    if (length == 0.)
    {
        pr ("warning: V2Direction: length = 0,\n"
            "will not normalize input vector\n");
        return (0);
    }
    b[0] = a[0] / length;
    b[1] = a[1] / length;
    return (length);
} /* end V2Direction() */


/* return 0..1 permutation idx[] such that length[idx[]] is ascending */
int *V2sort (double length[2], int idx[2])
{
    if (length[0] <= length[1]) SEQ2(idx);
    else 
    {
        idx[0] = 1;
        idx[1] = 0;
    }
    return (idx);
} /* end V2sort() */


/* return 0..1 random permutation index idx[] */
int *I2randompermutation (int idx[2])
{
    idx[0] = Fran(0,1);
    idx[1] = 1 - idx[0];
    return(idx);
} /* end I2randompermutation() */


/* return 0..1 permutation idx[] such that length[idx[]]  */
/* is non-decreasing. Equal length indices are randomized */
int *V2randomsort (double length[2], int idx[2])
{
    int tmp;
    I2RANDOMPERMUTATION (idx);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    return (idx);
} /* end V2randomsort() */


/* return 0..1 permutation idx[] such that length[idx[]] is     */
/* approximately increasing within epsilon. idx[] is randomized */
int *V2Randomsort (double length[2], int idx[2], double epsilon)
{
    int tmp;
    I2RANDOMPERMUTATION(idx);
    if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmp);
    BEPOSITIVE(epsilon);
    if ( SMALLSEPARATION(length[idx[0]],length[idx[1]],epsilon)
         && HALF_DECISION() )  SWAP(idx[0],idx[1],tmp);
    return (idx);
} /* end V2Randomsort() */
