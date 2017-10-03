#include "Elasticity.h"

/*************************************/
/* Nonlinear homogeneous deformation */
/*************************************/

static M3P HomogeneousWork_H0, HomogeneousWork_sout, HomogeneousWork_H;
static M3 HomogeneousWork_DeltaH;

/* How much work is done by changing U x dc */
double HomogeneousWork_Integrand (double c)
{
    M3 H, HI, DeltaH, FS;
    double volumec;
    M3ADDMULMUL (1-c, HomogeneousWork_H0, c, HomogeneousWork_H, H);
    M3InV (H, HI, volumec);
    M3MUL (HomogeneousWork_sout, HI, FS);
    return (volumec*M3TRPROD(HomogeneousWork_DeltaH,FS));
} /* end HomogeneousWork_Integrand() */


/* Calculate the work done by constant external stress sout[][] */
/* from H0[][] to H[][] by straight-path Romberg integration.   */
double HomogeneousWork (M3 H0, M3 sout, M3 H)
{
    double work;
    /* M3EQV ( H0,   HomogeneousWork_H0 ); */
    /* M3EQV ( sout, HomogeneousWork_sout ); */
    /* M3EQV ( H,    HomogeneousWork_H ); */
    HomogeneousWork_H0   = H0;
    HomogeneousWork_sout = sout;
    HomogeneousWork_H    = H;
    M3SUB (HomogeneousWork_H, HomogeneousWork_H0, HomogeneousWork_DeltaH);
    work = qromb(&HomogeneousWork_Integrand, 0., 1.);
    return (work);
} /* end HomogeneousWork() */


#ifdef _HomogeneousWork_TEST
int main (int argc, char *argv[])
{
    double pressure = 1;
    M3 H0, sout, H;
    TimeRandomize();
    M3Frandom (H0);
    M3AdddiaG (H0, 2);
    M3Frandom (H);
    M3AdddiaG (H, 2);
    M3Diagonal (-pressure, sout);
    printf ("%.16g %.16g\n", HomogeneousWork(H0, sout, H),
            -pressure*(M3VOLUME(H)-M3VOLUME(H0)));
    return (0);
}
#endif /* _HomogeneousWork_TEST */
