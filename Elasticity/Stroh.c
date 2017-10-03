#include "Elasticity.h"

/*****************/
/* Stroh tensors */
/*****************/

#ifdef _ElasticityStrohContraction_TEST
/* Hirth & Lothe: pp. 467: Eqn. 13-162 */
int main (int argc, char *argv[])
{
    int i,j,k,l;
    V3 a,b;
    M3 ab,AB;
    T4 C;
    TimeRandomize();
    for (i=0; i<6; i++)
        for (j=0; j<=i; j++)
            C[j][i] = C[i][j] = Frandom();
    V3Frandom(a);
    V3Frandom(b);
    ElasticityStrohContraction(a,b,C,ab);
    for (j=0; j<3; j++)
        for (k=0; k<3; k++)
        {
            AB[j][k] = 0;
            for (i=0; i<3; i++)
                for (l=0; l<3; l++)
                    AB[j][k] += a[i] * T4Element(C,i,j,k,l) * b[l];
            printf ("%e\n", AB[j][k]-ab[j][k]);
        }
    return (0);
}
#endif /* _ElasticityStrohContraction_TEST */


/* Subtract momentum derivatives from ElasticityStrohContraction() for */
/* uniformly moving system: D.J. Bacon, D.M. Barnett, R.O. Scattergood */
/* Progress in Materials Science 23 (1978) 53-262. pp. 136.            */
void ElasticityStrohCONTRACTION
(V3 a, V3 b, T4 C, double mass_density, V3 velocity, M3 ab)
{
    double cc;
    ElasticityStrohContraction(a,b,C,ab);
    cc = mass_density * V3DOT(a,velocity) * V3DOT(b,velocity);
    M3SubdiaG (ab, cc);
    return;
} /* end ElasticityStrohCONTRACTION() */


/* C[][] := A[][] * B[][] */
void ElasticityStrohM6mul (M6 A, M6 B, M6 C)
{
    int i,j,k;
    for (i=0; i<6; i++)
        for (j=0; j<6; j++)
        {
            C[i][j] = 0;
            for (k=0; k<6; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
    return;
} /* end ElasticityStrohM6mul() */


/* Diagonalize 6x6 real nonsymmetric matrix */
#define RANK  6
#define LWORK 256
void ElasticityStrohM6Diag (M6 A, V6 wr, V6 wi, M6 VR, M6 VI)
{
    int i,j,rank=RANK,lwork=LWORK,info;
    double ap[RANK][RANK],vp[RANK][RANK],work[LWORK];
    for (i=0; i<RANK; i++)
        for (j=0; j<RANK; j++)
            ap[i][j] = A[j][i];
    dgeev_ ("N", "V", &rank, ap[0], &rank, wr, wi,
            NULL, &rank, vp[0], &rank, work, &lwork, &info);
    /* printf ("optimal LWORK = %d\n", (int)work[0]); */
    if (info > 0)
        pe("ElasticityStrohM6Diag: QR algorithm failed to compute\n"
           "eigenvalues %d to %d, no eigenvectors were computed.\n",
           0, info-1);
    for (i=0; i<RANK; i++)
        if (wi[i] == 0)
            for (j=0; j<RANK; j++)
            {
                VR[j][i] = vp[i][j];
                VI[j][i] = 0;
            }
        else
        {
            for (j=0; j<RANK; j++)
            {
                VR[j][i] = vp[i][j];
                VI[j][i] = vp[i+1][j];
                VR[j][i+1] =  VR[j][i];
                VI[j][i+1] = -VI[j][i];
            }
            i++;
            continue;
        }
    return;
} /* end ElasticityStrohM6Diag() */
#undef LWORK
#undef RANK

#ifdef _ElasticityStrohM6Diag_TEST
int main (int argc, char *argv[])
{
    int i,j;
    M6 A,VR,VI,AVR,AVI;
    V6 wr,wi;
    TimeRandomize();
    for (i=0; i<6; i++)
        for (j=0; j<6; j++)
            A[i][j] = Frandom();
    ElasticityStrohM6Diag (A,wr,wi,VR,VI);
    S6PR("\n A = %M\n ", A);
    V6pr("wr = %M", wr);
    V6pr("wi = %M\n ", wi);
    S6PR("VR = %M", VR);
    S6PR("VI = %M\n ", VI);
    ElasticityStrohM6mul (A, VR, AVR);
    for (i=0; i<6; i++)
        for (j=0; j<6; j++)
            printf ("%e\n", AVR[j][i]-wr[i]*VR[j][i]+wi[i]*VI[j][i]);
    ElasticityStrohM6mul (A, VI, AVI);
    for (i=0; i<6; i++)
        for (j=0; j<6; j++)
            printf ("%e\n", AVI[j][i]-wr[i]*VI[j][i]-wi[i]*VR[j][i]);
    return (0);
}
#endif /* _ElasticityStrohM6Diag_TEST */
