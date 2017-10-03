/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"

/****************************/
/* Properties of H matrices */
/****************************/

#ifdef _Crystallographic2D_to_H_TEST
int main (int argc, char *argv[])
{
    double H[2][2];
    Crystallographic2D X;
    Crystallographic2DAssign(1.2,1.2,60,X);
    Crystallographic2D_to_H(X,H);
    S2PR("H = %M\n", H);
    H_to_Crystallographic2D(H,X);
    printf ("%f %f %f\n", X.a, X.b, X.gamma);
    return (0);
}
#endif /* _Crystallographic2D_to_H_TEST */


/* determine which of the three vertices of a col-parallelepiped */
/* formed by A[][0], A[][1] is the farthest from origin.         */
double M2maxcolumnradius (double H[2][2])
{
    int i;
    double dx[2], r2, r2max;
    double ds[3][2] = {{1,0},{0,1},{1,1}};
    for (r2max=i=0; i<3; i++)
    {
        M2mV2 (H, ds[i], dx);
        r2 = V2DOT (dx,dx);
        if (r2 > r2max) r2max = r2;
    }
    return (sqrt(r2max));
} /* end M2maxcolumnradius() */


/* determine which of the three vertices of a row-parallelepiped */
/* formed by A[0][], A[1][] is the farthest from origin.         */
double M2maxrowradius (double H[2][2])
{
    int i;
    double dx[2], r2, r2max;
    double ds[3][2] = {{1,0},{0,1},{1,1}};
    for (r2max=i=0; i<3; i++)
    {
        V2mM2 (ds[i], H, dx);
        r2 = V2DOT (dx,dx);
        if (r2 > r2max) r2max = r2;
    }
    return (sqrt(r2max));
} /* end M2maxrowradius() */


/* returns the thickness (>=0) of the parallelepiped formed */
/* by (H[0][], H[1][] in the row i direction.               */
double M2rowthickness (double H[2][2], int i)
{
    return (M2VOLUME(H) / V2LENGTH(H[1-i]));
} /* end M2rowthickness() */


/* returns the two thicknesses (>=0) of the parallelepiped */
/* formed by (H[0][], H[1][]), in thickness[].             */
double *M2rowthicknesses (double H[2][2], double thickness[2])
{
    double volume = M2VOLUME(H);
    thickness[0] = volume / V2LENGTH(H[1]);
    thickness[1] = volume / V2LENGTH(H[0]);
    return (thickness);
} /* end M2rowthicknesses() */


/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[0..2][] is guaranteed */
/* to cover sphere of radius R; returns total number of replicas */
int M2rows_to_cover_sphere (double H[2][2], double R, int nc[2])
{
    int i;
    static double thickness[2];
    M2rowthicknesses (H, thickness);
    for (i=0; i<2; i++)
        if (thickness[i] == 0.)
            pe ("M2rows_to_cover_sphere: thickness[%d] = 0.\n", i);
        else nc[i] =(int)ceil(ABS(R)/thickness[i]);
    return ((2*nc[0]+1)*(2*nc[1]+1));
} /* end M2rows_to_cover_sphere() */


/* return the thickness (>=0) of the parallelepiped formed */
/* by (H[][0], H[][1], H[][2]) in the column i direction.  */
double M2columnthickness (double H[2][2], int i)
{
    static double HT[2][2];
    M2TRANSPOSE (H, HT);
    return (M2rowthickness(HT,i));
} /* end M2columnthickness() */


/* return the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[][0], H[][1], H[][2]), in thickness[].      */
double *M2columnthicknesses (double H[2][2], double thickness[2])
{
    static double HT[2][2];
    M2TRANSPOSE (H, HT);
    return (M2rowthicknesses(HT,thickness));
} /* end M2columnthicknesses() */


/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[][0..2] is guaranteed */
/* to cover sphere of radius R; return total number of replicas. */
int M2columns_to_cover_sphere (double H[2][2], double R, int nc[2])
{
    static double HT[2][2];
    M2TRANSPOSE (H, HT);
    return (M2rows_to_cover_sphere(HT,R,nc));
} /* end M2columns_to_cover_sphere() */


/* dl^2 = dx_i * (1 + 2 * eta_{ij}) * dx_j, where dx meshes the undeformed */
/* because dx = ds * H0, dx' = ds * H, dl^2 = dx' * (dx')^T => above       */
void Lagrangian_strain2D (double H0[2][2], double H[2][2], double eta[2][2])
{
    double determinant;
    M2 H0I,J;
    M2INV (H0,H0I,determinant);
    if (determinant == 0)
        pe ("Lagrangian_strain2D: H0 is singular.\n");
    M2MUL (H0I,H,J);
    M2MULT (J,eta);
    M2SubdiaG (eta, 1);
    M2DividE (eta, 2);
    return;
} /* end Lagrangian_strain2D() */


/* achieve eta without rotation. M := sqrt(1 + 2*eta), H := H0 * M */
void pure_deform2D (double H0[2][2], double eta[2][2], double H[2][2])
{
    int i;
    static double eigval[2],L[2][2]={{0}},Q[2][2],QT[2][2],M[2][2],tmp[2][2];
    M2diag (eta, eigval, Q);
    for (i=0; i<2; i++)
    {
        eigval[i] = 1. + 2 * eigval[i];
        if (eigval[i] < 0)
        {
            Mfprintf (stderr, "\neta = %2M|| %15.8le %15.8le |\n ", eta[0]);
            pe ("pure_deform2D: Lagrangian strain infeasible.\n");
        }
        L[i][i] = sqrt(eigval[i]);
    }
    M2TRANSPOSE (Q, QT);
    M2MUL2 (QT, L, Q, M, tmp);
    if (H[0] == H0[0])
    { /* H and H0 are the same matrix */
        M2MUL (H0, M, tmp);
        M2EQV (tmp, H);
    }
    else M2MUL (H0, M, H);
    return;
} /* end pure_deform2D() */


#ifdef _pure_deform2D_TEST
int main (int argc, char *argv[])
{
    double H0[2][2], H[2][2], eta[2][2];
    M2FRANDOM(H0);
    S2PR("H0 = %M\n", H0);
    M2FRANDOM(eta);
    M2Symmetrize(eta);
    M2MultiplY(0.2,eta);
    S2PR("eta = %M\n", eta);
    pure_deform2D(H0,eta,H);
    S2PR("H = %M\n", H);
    Lagrangian_strain2D(H0,H,eta);
    S2PR("eta = %M\n", eta);
    return (0);
}
#endif /* _pure_deform2D_TEST */
