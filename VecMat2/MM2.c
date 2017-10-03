/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/********************************/
/* matrix & matrix with scalars */
/********************************/


/* matrix & matrix */


/* C[][] := A[][] + B[][] */
void M2add (double A[2][2], double B[2][2], double C[2][2])
{
    C[0][0] = A[0][0] + B[0][0];
    C[0][1] = A[0][1] + B[0][1];
    C[1][0] = A[1][0] + B[1][0];
    C[1][1] = A[1][1] + B[1][1];       
    return;
} /* end M2add() */


/* C[][] := A[][] - B[][] */
void M2sub (double A[2][2], double B[2][2], double C[2][2])
{
    C[0][0] = A[0][0] - B[0][0];
    C[0][1] = A[0][1] - B[0][1];
    C[1][0] = A[1][0] - B[1][0];         
    C[1][1] = A[1][1] - B[1][1];       
    return;
} /* end M2sub() */


/* C[][] := A[][] * B[][] */
void M2mul (double A[2][2], double B[2][2], double C[2][2])
{
    C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    return;
} /* end M2mul() */


/* D[][] := A[][] * B[][] * C[][] */
void M2mul2 (double A[2][2], double B[2][2], double C[2][2], double D[2][2])
{
    static double E[2][2];
    M2MUL (A, B, E);
    M2MUL (E, C, D);
    return;
} /* end M2mul2() */


/* matrix & matrix & scalar */
