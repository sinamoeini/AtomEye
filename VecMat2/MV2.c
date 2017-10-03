/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/***************************/
/* column space operations */
/***************************/


/* columnlength[i] := | A[][i] |; returns columnlength[] */
double *M2columnlengths (double A[2][2], double columnlength[2])
{
    columnlength[0] = DISTANCE2D (A[0][0],A[1][0]);
    columnlength[1] = DISTANCE2D (A[0][1],A[1][1]);
    return (columnlength);
} /* end M2columnlengths() */


/* returns the maximum Euclidean length of A[][0], A[][1] */
double M2maxcolumnlength (double A[2][2])
{
    double length0, length1;
    length0 = DISTANCE2D (A[0][0],A[1][0]);
    length1 = DISTANCE2D (A[0][1],A[1][1]);
    return (MAX(length0,length1));
} /* end M2maxcolumnlength() */


/* column[] := A[][i]; return column[] */
double *M2column (double A[2][2], int i, double column[2])
{
    column[0] = A[0][i];
    column[1] = A[1][i];
    return(column);
} /* end M2column() */


/* c[] := A[][] * b[]; then return c[] */
double *M2mulV2 (double A[2][2], double b[2], double c[2])
{
    c[0] = A[0][0]*b[0] + A[0][1]*b[1];
    c[1] = A[1][0]*b[0] + A[1][1]*b[1];
    return(c);
} /* end M2mulV2() */


/* b[] := A[][] * b[]; then return b[] */
double *M2MULV2 (double A[2][2], double b[2])
{
    double a[2];
    a[0] = b[0];
    a[1] = b[1];
    b[0] = A[0][0]*a[0] + A[0][1]*a[1];
    b[1] = A[1][0]*a[0] + A[1][1]*a[1];
    return(b);
} /* end M2MULV2() */


/* d[] := A[][] * b[] + c[]; then return d[] */
double *M2muladdV2 (double A[2][2], double b[2], double c[2], double d[2])
{
    d[0] = A[0][0]*b[0] + A[0][1]*b[1] + c[0];
    d[1] = A[1][0]*b[0] + A[1][1]*b[1] + c[1];
    return(d);
} /* end M2muladdV2() */


/* c[] := c[] + A[][] * b[]; then return c[] */
double *M2mulADDV2 (double A[2][2], double b[2], double c[2])
{
    c[0] += A[0][0]*b[0] + A[0][1]*b[1];
    c[1] += A[1][0]*b[0] + A[1][1]*b[1];
    return(c);
} /* end M2mulADDV2() */


/* d[] := A[][] * b[] - c[]; then return d[] */
double *M2mulsubV2 (double A[2][2], double b[2], double c[2], double d[2])
{
    d[0] = A[0][0]*b[0] + A[0][1]*b[1] - c[0];
    d[1] = A[1][0]*b[0] + A[1][1]*b[1] - c[1];
    return(d);
} /* end M2mulsubV2() */


/* d[] := c[] - A[][] * b[]; then return d[] */
double *M2submulV2 (double c[2], double A[2][2], double b[2], double d[2])
{
    d[0] = c[0] - A[0][0]*b[0] - A[0][1]*b[1];
    d[1] = c[1] - A[1][0]*b[0] - A[1][1]*b[1];
    return(d);
} /* end M2submulV2() */


/* c[] := c[] - A[][] * b[]; then return c[] */
double *M2mulSUBV2 (double A[2][2], double b[2], double c[2])
{
    c[0] -= A[0][0]*b[0] + A[0][1]*b[1];
    c[1] -= A[1][0]*b[0] + A[1][1]*b[1];
    return (c);
} /* end M2mulSUBV2() */
