/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#include "VecMat2.h"


/*********************************************************/
/* Matrix inversion and symmetric matrix diagonalization */
/*********************************************************/


/* B[][] := A[][]^-1; return det(A) */
double M2inv (double A[2][2], double B[2][2])
{
    double determinant;
    determinant = M2DETERMINANT(A);
    if (determinant == 0.)
        pe ("M2inv: determinant = %e, matrix is singular\n", determinant);
    B[0][0] =  A[1][1] / determinant;
    B[1][1] =  A[0][0] / determinant;
    B[1][0] = -A[1][0] / determinant;
    B[0][1] = -A[0][1] / determinant;
    return (determinant);
} /* end M2inv() */


#ifdef _M2inv_TEST
int main()
{
    double A[2][2], B[2][2], C[2][2];
    M2FRANDOM(A);
    S2PR("\nA = %M\n ", A);
    printf ("determinant = %f\n", M2inv(A,B));
    S2PR("\nB = %M\n ", B);
    M2mul(A,B,C);
    S2PR("A * B = %M\n ", C);
    return(0);
}
#endif  /* _M2inv_TEST */


/* A[][] := A[][]^-1; return original det(A) */
double M2Inv (double A[2][2])
{
    double tmp, determinant = M2DETERMINANT(A);
    if (determinant == 0.)
        pe ("M2Inv: determinant = %e, matrix is singular\n", determinant);
    tmp = A[0][0];
    A[0][0] = A[1][1] / determinant;
    A[1][1] =     tmp / determinant;
    A[1][0] /= -determinant;
    A[0][1] /= -determinant;
    return (determinant);
} /* end M2Inv() */


#ifdef _M2Inv_TEST
int main()
{
    double A[2][2], B[2][2], C[2][2];
    M2FRANDOM(A);
    M2EQV(A,B);
    S2PR("\nA = %M\n ", A);
    printf ("determinant = %f\n", M2Inv(A));
    S2PR("\nAfter M2Inv, A = %M\n ", A);
    M2mul(A,B,C);
    S2PR("A * A^-1 = %M\n ", C);
    return(0);
}
#endif  /* _M2Inv_TEST */


/******************************************************************/
/* Diagonalize 2x2 real symmetric matrix: A = V^T*Diag(eigval)*V, */
/* eigval[i] will be stored in ascending order, i = 0..1; the     */
/* corresponding eigenvectors V[i][] form a right-handed system.  */
/* Return index i of the eigenvalue of largest absolute value.    */
/******************************************************************/
int M2diag (double A[2][2], double eigval[2], double V[2][2])
{
    double B[2][2];
    double scale, b, c, d, tmp;
    /* find the largest component of A[][] in absolute value */
    M2INFNORM (A, scale);
    if ( scale == 0. )
    { /* this is a perfectly zero matrix */
        M2IDENTITY (V);
        V2ZERO (eigval);
        return (0);
    }
    /* set up the zoom: safe to do so */
    M2DIVIDE (A, scale, B);
    /* we really go after numerical accuracy */
    if ( M2ISSYMMETRIC(B) ) M2Symmetrize(B);
    else
    {
        Mfprintf (stderr, "\n%2M|| %15.8le %15.8le |\n ", A[0]);
        pe ("M2diag: matrix is not symmetric\n");
    }
    /* eigenvalue equation: x^2 - b*x + c = 0 */
    b = M2TR(B);
    c = M2DETERMINANT(B);
    d = sqrt(SQUARE(b)-4*c);
    eigval[0] = ( b - d ) / 2;
    eigval[1] = eigval[0] + d;
    M2SubdiaG (B, eigval[0]);
    V[0][0] = V2LENGTH(B[1]);
    V[0][1] = V2LENGTH(B[0]) * ( (V2DOT(B[0],B[1]) > 0) ? -1 : 1 );
    if (V2EQZERO(V[0])) V[0][0] = 1;
    else V2NORMALIZE(V[0], tmp);
    V2rotate90 (V[0], V[1]);
    V2MuL (scale,eigval);
    if (fabs(eigval[0]) > fabs(eigval[1])) return(0);
    else return(1);
} /* end M2diag() */

#ifdef _M2diag_TEST
int main()
{
    double A[2][2],eigval[2],Q[2][2],QT[2][2],L[2][2],QI[2][2];
    TimeRandomize();
    /* M2FRANDOM (A); */
    /* M2Symmetrize (A); */
    M2IDENTITY(A);
    A[1][1] = 0;
    M2diag (A,eigval,Q);
    S2PR ("\nA = %M\n ", A);
    V2pr ("eigval[] = %M\n ", eigval);
    printf ("%25.18e %25.18e\n",eigval[0], eigval[1]);
    S2PR ("\nQ = %M\n ", Q);
    M2TRANSPOSE (Q,QT);
    M2MUL2 (Q,A,QT,L,QI);
    S2PR ("Q * A * Q^T = %M\n ", L);
    M2MUL (Q,QT,L);
    S2PR("Q * Q^T = %M\n ", L);
    printf ("det|Q| = %f\n\n", M2DETERMINANT(Q));
    return(0);
}
#endif  /* _M2diag_TEST */


/* A = M*R, where M is symmetric, R (do not return) is orthogonal */
void M2RightPolarDecompose (M2 A, M2 M)
{
    V2 eigval;
    M2 A2,V,E,VT,TMP;
    if (M2EQZERO(A))
    {
        M2ZERO(M);
        return;
    }
    M2MULT(A,A2);
    M2diag(A2,eigval,V);
    M2diagonal(sqrt(eigval[0]),sqrt(eigval[1]),E);
    M2TRANSPOSE(V,VT);
    M2MUL2 (VT,E,V,M,TMP);
    return;
} /* end M2RightPolarDecompose() */


/* A = M*R, where M is symmetric, R is orthogonal */
void M2RightPolarDECOMPOSE (M2 A, M2 M, M2 R)
{
    M2 MI;
    M2RightPolarDecompose (A, M);
    if (M2DETERMINANT(M) == 0.)
        pe("M2RightPolarDECOMPOSE: A is singular.\n");
    M2inv (M,MI);
    M2MUL (MI,A,R);
    return;
} /* end M2RightPolarDECOMPOSE() */

#ifdef _M2RightPolarDECOMPOSE_TEST
int main (int argc, char *argv[])
{
    M2 A = GuineaPigM2, M,R,RT,RRT, RESIDUE;
    M2RightPolarDECOMPOSE (A, M, R);
    M2MUL (M, R, RESIDUE);
    M2SuB (RESIDUE, A);
    S2PR("residue = %M", RESIDUE);
    S2PR("R = %M", R);
    M2TRANSPOSE (R,RT);
    M2MUL (R,RT,RRT);
    S2PR("R * R^T = %M", RRT);
    return (0);
}
#endif /* _M2RightPolarDECOMPOSE_TEST */


/* A = L*M, where L (do not return) is orthogonal, M is symmetric */
void M2LeftPolarDecompose (M2 A, M2 M)
{
    V2 eigval;
    M2 A2,V,E,VT,TMP;
    if (M2EQZERO(A))
    {
        M2ZERO(M);
        return;
    }
    M2TMUL(A,A2);
    M2diag(A2,eigval,V);
    M2diagonal(sqrt(eigval[0]),sqrt(eigval[1]),E);
    M2TRANSPOSE(V,VT);
    M2MUL2 (VT,E,V,M,TMP);
    return;
} /* end M2LeftPolarDecompose() */


/* A = L*M, where L is orthogonal, M is symmetric */
void M2LeftPolarDECOMPOSE (M2 A, M2 L, M2 M)
{
    M2 MI;
    M2LeftPolarDecompose (A, M);
    if (M2DETERMINANT(M) == 0.)
        pe("M2LeftPolarDECOMPOSE: A is singular.\n");
    M2inv (M,MI);
    M2MUL (A,MI,L);
    return;
} /* end M2LeftPolarDECOMPOSE() */

#ifdef _M2LeftPolarDECOMPOSE_TEST
int main (int argc, char *argv[])
{
    M2 A = GuineaPigM2, L,M,LT,LLT, RESIDUE;
    M2LeftPolarDECOMPOSE (A, L, M);
    M2MUL (L, M, RESIDUE);
    M2SuB (RESIDUE, A);
    S2PR("residue = %M", RESIDUE);
    S2PR("L = %M", L);
    M2TRANSPOSE (L,LT);
    M2MUL (L,LT,LLT);
    S2PR("L * L^T = %M", LLT);
    return (0);
}
#endif /* _M2LeftPolarDECOMPOSE_TEST */
