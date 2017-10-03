/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/*******************************/
/* vector & vector with scalar */
/*******************************/

/* vector & vector */

/* c[] := a[] + b[]; then return c[] */
double *V2add (double a[2], double b[2], double c[2])
{
    c[0] = a[0] + b[0];
    c[1] = a[1] + b[1];
    return (c);
} /* end V2add() */


/* b[] := a[] + b[]; then return b[] */
double *V2Add (double a[2], double b[2])
{
    b[0] += a[0];
    b[1] += a[1];
    return (b);
} /* end V2Add() */


/* c[] := a[] - b[]; then return c[] */
double *V2sub (double a[2], double b[2], double c[2])
{
    c[0] = a[0] - b[0];
    c[1] = a[1] - b[1];
    return (c);
} /* end V2sub() */


/* a[] := a[] - b[]; then return a[] */
double *V2Sub (double a[2], double b[2])
{
    a[0] -= b[0];
    a[1] -= b[1];
    return (a);
} /* end V2SUB() */


/* dot product of a[] and b[] */
double V2dot (double a[2], double b[2])
{
    return (a[0]*b[0]+a[1]*b[1]);
} /* end V2dot() */


/* vector & vector & scalar */


/* c[] := a[] + multiplier * b[]; then return c[] */
double *V2addmul (double a[2], double multiplier, double b[2], double c[2])
{
    c[0] = a[0] + multiplier * b[0];
    c[1] = a[1] + multiplier * b[1];
    return (c);
} /* end V2addmul() */


/* c[] := a[] + b[] / divisor; then return c[] */
double *V2adddiv (double a[2], double b[2], double divisor, double c[2])
{
    if (divisor == 0.) pe ("error: V2adddiv: divisor = %e\n", divisor);
    c[0] = a[0] + b[0] / divisor;
    c[1] = a[1] + b[1] / divisor;
    return (c);
} /* end V2adddiv() */


/* b[] := b[] + multiplier * a[]; then return b[] */
double *V2ADDmul (double multiplier, double a[2], double b[2])
{
    b[0] += multiplier * a[0];
    b[1] += multiplier * a[1];
    return (b);
} /* end V2ADDmul() */


/* b[] := b[] + a[] / divisor; then return b[] */
double *V2ADDdiv (double a[2], double divisor, double b[2])
{
    if (divisor == 0.) pe ("V2ADDdiv: divisor = %e\n", divisor);
    b[0] += a[0] / divisor;
    b[1] += a[1] / divisor;
    return (b);
} /* end V2ADDdiv() */


/* c[] := a[] - multiplier * b[]; then return c[] */
double *V2submul (double a[2], double multiplier, double b[2], double c[2])
{
    c[0] = a[0] - multiplier * b[0];
    c[1] = a[1] - multiplier * b[1];
    return (c);
} /* end V2submul() */


/* c[] := a[] - b[] / divisor; then return c[] */
double *V2subdiv (double a[2], double b[2], double divisor, double c[2])
{
    if (divisor == 0.) pe ("V2subdiv: divisor = %e\n", divisor);
    c[0] = a[0] - b[0] / divisor;
    c[1] = a[1] - b[1] / divisor;
    return (c);
} /* end V2subdiv() */


/* a[] := a[] - multiplier * b[]; then return a[] */
double *V2SUBmul (double a[2], double multiplier, double b[2])
{
    a[0] -= multiplier * b[0];
    a[1] -= multiplier * b[1];
    return (a);
} /* end V2SUBmul() */


/* a[] := a[] - b[] / divisor; then return a[] */
double *V2SUBdiv (double a[2], double b[2], double divisor)
{
    if (divisor == 0.) pe ("V2SUBdiv: divisor = %e\n", divisor);
    a[0] -= b[0] / divisor;
    a[1] -= b[1] / divisor;
    return (a);
} /* end V2SUBdiv() */


/* c[] := part of a[] that is perpendicular to b[]; return c[] */
double *V2perpendicular (double a[2], double b[2], double c[2])
{
    double b2 = V2LENGTH2 (b);
    if ( b2 == 0 ) pe ("V2perpendicular: b[] = 0\n");
    return ( V2submul(a, V2DOT(a,b)/b2, b, c) );
} /* end V2perpendicular() */


/* a[] := part of a[] that is perpendicular to b[]; return modified a[] */
double *V2Perpendicular (double a[2], double b[2])
{
    double b2 = V2LENGTH2 (b);
    if ( b2 == 0 ) pe ("V2Perpendicular: b[] = 0\n");
    return ( V2SUBmul(a, V2DOT(a,b)/b2, b) );
} /* end V2Perpendicular() */
