/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/************************/
/* row space operations */
/************************/


/* rowlength[i] := | A[i][] |; returns rowlength[] */
double *M2rowlengths (double A[2][2], double rowlength[2])
{
    rowlength[0] = V2LENGTH (A[0]);
    rowlength[1] = V2LENGTH (A[1]);
    return (rowlength);
} /* end M2rowlengths() */


/* returns the maximum Euclidean length of A[0][], A[1][], A[2][] */
double M2maxrowlength (double A[2][2])
{
    double length0, length1;
    length0 = V2LENGTH (A[0]);
    length1 = V2LENGTH (A[1]);
    return (MAX(length0,length1));
} /* end M2maxrowlength() */


/* c[] := a[] * B[][]; then return c[] */
double *V2mulM2 (double a[2], double B[2][2], double c[2])
{
    c[0] = a[0]*B[0][0] + a[1]*B[1][0];
    c[1] = a[0]*B[0][1] + a[1]*B[1][1];
    return(c);
} /* end V2mulM2() */


/* a[] := a[] * B[][]; then return a[] */
double *V2MULM2 (double a[2], double B[2][2])
{
    double c[2];
    c[0] = a[0];
    c[1] = a[1];
    a[0] = c[0]*B[0][0] + c[1]*B[1][0];
    a[1] = c[0]*B[0][1] + c[1]*B[1][1];
    return(a);
} /* end V2MULM2() */


/* d[] := a[] + b[] * C[][]; then return d[] */
double *V2addmulM2 (double a[2], double b[2], double C[2][2], double d[2])
{
    d[0] = a[0] + b[0]*C[0][0] + b[1]*C[1][0];
    d[1] = a[1] + b[0]*C[0][1] + b[1]*C[1][1];
    return(d);
} /* end V2addmulM2() */


/* a[] := a[] + b[] * C[][]; then return a[] */
double *V2ADDmulM2 (double a[2], double b[2], double C[2][2])
{
    a[0] += b[0]*C[0][0] + b[1]*C[1][0];
    a[1] += b[0]*C[0][1] + b[1]*C[1][1];
    return(a);
} /* end V2ADDmulM2() */


/* d[] := a[] - b[] * C[][]; then return d[] */
double *V2submulM2 (double a[2], double b[2], double C[2][2], double d[2])
{
    d[0] = a[0] - b[0]*C[0][0] - b[1]*C[1][0];
    d[1] = a[1] - b[0]*C[0][1] - b[1]*C[1][1];
    return(d);
} /* end V2submulM2() */


/* d[] := b[] * C[][] - a[]; then return d[] */
double *V2mulsubM2 (double b[2], double C[2][2], double a[2], double d[2])
{
    d[0] = b[0]*C[0][0] + b[1]*C[1][0] - a[0];
    d[1] = b[0]*C[0][1] + b[1]*C[1][1] - a[1];
    return(d);
} /* end V2mulsubM2() */


/* a[] := a[] - b[] * C[][]; then return a[] */
double *V2SUBmulM2 (double a[2], double b[2], double C[2][2])
{
    a[0] -= b[0]*C[0][0] + b[1]*C[1][0];
    a[1] -= b[0]*C[0][1] + b[1]*C[1][1];
    return(a);
} /* end V2SUBmulM2() */
