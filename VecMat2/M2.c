/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/**********************************************/
/* single matrix & scalar & generating vector */
/**********************************************/


/* matrix */


/* generate a zero matrix A[][] := 0 */
void M2zero (double A[2][2])
{
    A[0][0] = 0.;
    A[0][1] = 0.;
    A[1][0] = 0.;
    A[1][1] = 0.;
    return;
} /* end M2zero() */


/* generate an identity matrix A[][] := I[][] */
void M2identity (double A[2][2])
{
    A[0][0] = 1.;
    A[0][1] = 0.;
    A[1][0] = 0.;
    A[1][1] = 1.;
    return;
} /* end M2identity() */


/* generate A[][] with 4 independent random components on (0.,1.) */
void M2frandom (double A[2][2])
{
    A[0][0] = Frandom();
    A[0][1] = Frandom();
    A[1][0] = Frandom();
    A[1][1] = Frandom();
    return;
} /* end M2frandom() */


/* B[][] := A[][] */
void M2eqv (double A[2][2], double B[2][2])
{
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
    return;
} /* end M2eqv() */


/* arbitrary exponent-th norm of A[][] (all components) */
double M2norm (double A[2][2], double exponent)
{
    double tmp;
    if (exponent <= 0.)
        pe("M2norm: norm exponent = %lf <= 0 is illegal\n", exponent);
    else if ( exponent >= EXPONENT_INFINITY )
    {
        tmp = ABS( A[0][0] );
        if ( A[0][1] > tmp ) tmp = A[0][1];
        else if ( -A[0][1] > tmp ) tmp = -A[0][1];
	if ( A[1][0] > tmp ) tmp = A[1][0];
	else if ( -A[1][0] > tmp ) tmp = -A[1][0];
	if ( A[1][1] > tmp ) tmp = A[1][1];
	else if ( -A[1][1] > tmp ) tmp = -A[1][1];
	return (tmp);
    }
    else
    {
        tmp = pow(ABS(A[0][0]), exponent) +
            pow(ABS(A[0][1]), exponent) +
            pow(ABS(A[1][0]), exponent) +
            pow(ABS(A[1][1]), exponent);
        return (pow(tmp, 1./exponent));
    }
    return (0.);
} /* end M2norm() */


/* return the largest component of A[][] in absolute value */
double M2infnorm (double A[2][2])
{
    register double tmp = ABS( A[0][0] );
    if ( A[0][1] > tmp ) tmp = A[0][1];
    else if ( -A[0][1] > tmp ) tmp = -A[0][1];
    if ( A[1][0] > tmp ) tmp = A[1][0];
    else if ( -A[1][0] > tmp ) tmp = -A[1][0];
    if ( A[1][1] > tmp ) tmp = A[1][1];
    else if ( -A[1][1] > tmp ) tmp = -A[1][1];
    return (tmp);
} /* end M2infnorm() */


/* return sqrt{ \sum_{i=0..1} \sum_{j=0..1} |A_ij|^2 } */
double M22norm (double A[2][2])
{
    return (sqrt(V2LENGTH2(A[0])+V2LENGTH2(A[1])));
} /* end M22norm() */


/* B[][] := -A[][] */
void M2neg (double A[2][2], double B[2][2])
{
    B[0][0] = -A[0][0];
    B[0][1] = -A[0][1];
    B[1][0] = -A[1][0];
    B[1][1] = -A[1][1];
    return;
} /* end M2neg() */


/* A[][] := -A[][] */
void M2Neg (double A[2][2])
{
    A[0][0] = -A[0][0];
    A[0][1] = -A[0][1];
    A[1][0] = -A[1][0];
    A[1][1] = -A[1][1];
    return;
} /* end M2Neg() */


/* B[][] := A[][]' */
void M2transpose (double A[2][2], double B[2][2])
{
    B[0][0] = A[0][0];
    B[0][1] = A[1][0];
    B[1][0] = A[0][1];
    B[1][1] = A[1][1];
    return;
} /* end M2transpose() */


/* A[][] := A[][]' */
void M2Transpose (double A[2][2])
{
    double tmp;
    tmp = A[0][1];
    A[0][1] = A[1][0];
    A[1][0] = tmp;
    return;
} /* end M2Transpose() */


/* B[][] := (A[][]+A[][]')/2 */
void M2symmetrize (double A[2][2], double B[2][2])
{
    B[0][0] = A[0][0];
    B[1][1] = A[1][1];
    B[0][1] = B[1][0] = (A[1][0]+A[0][1])/2.;
    return;
} /* end M2symmetrize() */


/* return the trace of A[][] */
double M2Tr (double A[2][2])
{
    return (A[0][0]+A[1][1]);
} /* end M2Tr() */


/* B[][] := trace(A[][]) / 2 * I[][]; return trace(A) */
double M2trace (double A[2][2], double B[2][2])
{
    double trace = A[0][0]+A[1][1];
    B[0][0] = trace/2.;
    B[1][0] = B[0][1] = 0.;
    B[1][1] = trace/2.;
    return (trace);
} /* return M2trace() */


/* A[][] := trace(A[][])/2 * I[][]; return original trace(A) */
double M2Trace (double A[2][2])
{
    double trace = A[0][0]+A[1][1];
    A[0][0] = trace/2.;
    A[1][0] = A[0][1] = 0.;
    A[1][1] = trace/2.;
    return (trace);
} /* return M2Trace() */


/* B[][] := A[][] - trace(A[][])/2 * I[][]; return trace(A) */
double M2traceless (double A[2][2], double B[2][2])
{
    double trace = A[0][0]+A[1][1];
    B[0][0] = A[0][0] - trace/2.;
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] - trace/2.;
    return (trace);
} /* return M2traceless() */


/* A[][] := A[][] - trace(A[][])/2 * I[][]; return original trace(A) */
double M2Traceless (double A[2][2])
{
    double trace = A[0][0]+A[1][1];
    A[0][0] -= trace/2.;
    A[1][1] -= trace/2.;
    return (trace);
} /* return M2Traceless() */


/* decompose A[][] to b*I[][] + C[][], where b := trace(A)/2; return b */
double M2tracedecompose (double A[2][2], double B[2][2], double C[2][2])
{
    double trace = A[0][0] + A[1][1];
    B[0][0] = trace/2.;
    B[1][0] = B[0][1] = 0.;
    B[1][1] = trace/2.;
    C[0][0] = A[0][0] - trace/2.;
    C[0][1] = A[0][1];
    C[1][0] = A[1][0];
    C[1][1] = A[1][1] - trace/2.;
    return (trace);
} /* end M2tracedecompose() */



/* matrix & scalar */



/* generate an identity matrix A[][] := a x I[][] */
void M2Identity (double a, double A[2][2])
{
    A[0][0] = a;
    A[1][0] = A[0][1] = 0.;
    A[1][1] = a;
    return;
} /* end M2Identity() */


/* B[][] := multiplier * A[][] */
void M2multiply (double multiplier, double A[2][2], double B[2][2])
{
    B[0][0] = multiplier * A[0][0];
    B[0][1] = multiplier * A[0][1];
    B[1][0] = multiplier * A[1][0];
    B[1][1] = multiplier * A[1][1];
    return;
} /* end M2multiply() */


/* A[][] := multiplier * A[][] */
void M2Multiply (double multiplier, double A[2][2])
{
    A[0][0] *= multiplier;
    A[0][1] *= multiplier;
    A[1][0] *= multiplier;
    A[1][1] *= multiplier;
    return;
} /* end M2Multiply() */


/* B[][] := A[][] / divisor */
void M2divide (double A[2][2], double divisor, double B[2][2])
{
    if (divisor == 0.) pe ("M2divide: divisor = %e\n", divisor);
    B[0][0] = A[0][0] / divisor;
    B[0][1] = A[0][1] / divisor;
    B[1][0] = A[1][0] / divisor;
    B[1][1] = A[1][1] / divisor;
    return;
} /* end M2divide() */


/* A[][] := A[][] / divisor */
void M2Divide (double A[2][2], double divisor)
{
    if (divisor == 0.) pe ("M2Divide: divisor = %e\n", divisor);
    A[0][0] /= divisor;
    A[0][1] /= divisor;
    A[1][0] /= divisor;
    A[1][1] /= divisor;
    return;
} /* end M2Divide() */


/* B[][] := A[][] + a * I */
void M2adddiag (double A[2][2], double a, double B[2][2])
{
    B[0][0] = A[0][0] + a;
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] + a;
    return;
} /* end M2adddiag() */


/* A[][] := A[][] + a * I */
void M2Adddiag (double A[2][2], double a)
{
    A[0][0] += a;
    A[1][1] += a;
    return;
} /* end M2Adddiag() */


/* B[][] := A[][] - a * I */
void M2subdiag (double A[2][2], double a, double B[2][2])
{
    B[0][0] = A[0][0] - a;
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1] - a;
    return;
} /* end M2subdiag() */


/* A[][] := A[][] - a * I */
void M2Subdiag (double A[2][2], double a)
{
    A[0][0] -= a;
    A[1][1] -= a;
    return;
} /* end M2Subdiag() */
