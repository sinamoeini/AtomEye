/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/****************************************/
/* rotation generator & representations */
/****************************************/


/* b[] := right-handed rotation of a[] by theta [RAD] */
double *V2rotate (double a[2], double theta, double b[2])
{
    double ra[2], cosine, sine;
    V2rotate90 (a,ra);
    cosine = cos(theta);
    sine = sin(theta);
    V2ADDMULMUL (cosine,a,sine,ra,b);
    return (b);
} /* end V2rotate() */


/* If a rotation makes a[] -> b[], what happens to v[]?    */
/* Assume a[] and b[] are already NORMALIZED; returns v[]. */
double *V2geodesic (double a[2], double b[2], double v[2])
{
    double ra[2], rb[2], va, vra;
    V2rotate90(a,ra);
    V2rotate90(b,rb);
    va  = V2DOT(v,a);
    vra = V2DOT(v,ra);
    V2ADDMULMUL(va,b,vra,rb,v);
    return(v);
} /* end V2geodesic() */


/* Compute the rotational matrix R[][] corresponding to */
/* a rotation that makes a[]->b[] as v[] := v[]*R[][]   */
/* Assume a[] and b[] are already NORMALIZED.           */
void M2geodesic (double a[2], double b[2], double R[2][2])
{
    R[0][0] = 1;
    R[0][1] = 0;
    V2geodesic (a, b, R[0]);
    R[1][0] = -R[0][1];
    R[1][1] = R[0][0];
    return;
} /* end M2geodesic() */


#ifdef _M2geodesic_TEST
int main (int argc, char *argv[])
{
    double a[2]={1,0}, b[2]={1/SQRT2,1/SQRT2}, R[2][2];
    M2geodesic (a, b, R);
    S2PR ("%M",R);
    return (0);
}
#endif /* _M2geodesic_TEST */
