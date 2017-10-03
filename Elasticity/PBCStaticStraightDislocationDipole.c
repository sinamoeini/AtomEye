#include "Elasticity.h"

/**********************************************************/
/* Straight/static dislocation dipole along H[2][] of PBC */
/**********************************************************/

const int
ElasticityPBCStaticStraightDislocationDipoleSufficientImages[2]={30,30};


/***************************************************************/
/* Stroh solution in PBC: Hfree[][] is the PBC supercell in a  */
/* stress-free crystal in the observer frame. C is the elastic */
/* constant of the crystal in observer frame. The dislocation  */
/* at s0[] has xi // H[2][] and Burgers vector b[]. The cut is */
/* generated from s0[] to s1[] in the convex supercell. The    */
/* Green's function images will be summed                      */
/* -images[0]:images[0] by -images[1]:images[1] times.         */
/***************************************************************/
void ElasticityPBCStaticStraightDislocationDipoleAssign
(M3 Hfree, T4 C, double s0[2], double s1[2], V3 B, int images[2],
 ElasticityPBCStaticStraightDislocationDipole *dipole)
{
    int i, j, k;
    V3 R, tmp, u;
    T2 TMP;
    M3 dxref;

    M3EQV(Hfree, dipole->Hfree);
    dipole->linelength = V3LENGTH(dipole->Hfree[2]);
    dipole->volumefree = M3VOLUME(dipole->Hfree);
    ElasticityT4Eqv(C, dipole->C);
    ElasticityT4Inverse (dipole->C, dipole->S);
    V2EQV(s0, dipole->s0);
    V2Trim(dipole->s0);
    V3ADDMULMUL(dipole->s0[0], dipole->Hfree[0],
                dipole->s0[1], dipole->Hfree[1],
                dipole->x0);
    V2EQV(s1, dipole->s1);
    V2Trim(dipole->s1);
    V3ADDMULMUL(dipole->s1[0]-dipole->s0[0], dipole->Hfree[0],
                dipole->s1[1]-dipole->s0[1], dipole->Hfree[1],
                dipole->dx);
    V3ADDMULMUL(1-dipole->s0[0], dipole->Hfree[0],
                0-dipole->s0[1], dipole->Hfree[1], dxref[0]); /* s=(1,0) */
    V3ADDMULMUL(0-dipole->s0[0], dipole->Hfree[0],
                1-dipole->s0[1], dipole->Hfree[1], dxref[1]); /* s=(0,1) */
    V3ADDMULMUL(0-dipole->s0[0], dipole->Hfree[0],
                0-dipole->s0[1], dipole->Hfree[1], dxref[2]); /* s=(0,0) */

    ElasticityStraightStaticDislocation
        (dipole->dx, dipole->Hfree[2], B, dipole->C, dipole->d);
    V3mul(- ElasticityStraightDislocationCartesianR(dipole->d,dipole->dx) *
          dipole->linelength / dipole->volumefree,
          dipole->d->n, tmp);
    M3ASSIGNSYMMETRIZEDV3V3(tmp, dipole->d->b, dipole->plasticStrain);
    V2EQV(images, dipole->images);

    M3ZERO(dipole->greenu);
    dipole->green->Ebench = 0;
    M3ZERO(dipole->green->Stress0);
    M3ZERO(dipole->green->Stress1);
    for (i=-dipole->images[0]; i<=dipole->images[0]; i++)
        for (j=-dipole->images[1]; j<=dipole->images[1]; j++)
        {
            V3ADDMULMUL(i, dipole->Hfree[0], j, dipole->Hfree[1], R);
            for (k=0; k<3; k++)
            { /* supercell floats */
                V3SUB(dxref[k], R, tmp);
                ElasticityStraightDislocation_u(dipole->d, tmp, u);
                V3AdD(u, dipole->greenu[k]);
                V3SuB(tmp, dipole->dx);
                ElasticityStraightDislocation_u(dipole->d, tmp, u);
                V3SuB (dipole->greenu[k], u);
            }
            if ( (i==0) && (j==0) )
            { /* energy and stress */
                dipole->green->Ebench +=
                    ElasticityStraightDislocationStaticDipoleSelfEnergy
                    (dipole->d, dipole->dx);
                ElasticityStraightDislocation_stress
                    (dipole->d, dipole->dx, TMP);
                M3AdD(TMP, dipole->green->Stress1);
                V3NEG(dipole->dx, tmp);
                ElasticityStraightDislocation_stress (dipole->d, tmp, TMP);
                M3SuB(dipole->green->Stress0, TMP);
            }
            else
            {
                dipole->green->Ebench +=
                    ElasticityStraightDislocationStaticDipoleInteractionEnergy
                    (dipole->d, dipole->dx, R) / 2;
                V3NeG(R);
                ElasticityStraightDislocation_stress (dipole->d, R, TMP);
                M3AdD(TMP, dipole->green->Stress0);
                M3SuB(dipole->green->Stress1, TMP);
                V3ADD(R, dipole->dx, tmp);
                ElasticityStraightDislocation_stress (dipole->d, tmp, TMP);
                M3AdD(TMP, dipole->green->Stress1);
                V3SUB(R, dipole->dx, tmp);
                ElasticityStraightDislocation_stress (dipole->d, tmp, TMP);
                M3SuB(dipole->green->Stress0, TMP);
            }
        }
    dipole->green->Ebench *= dipole->linelength;

    M3EQV(dipole->Hfree, dipole->green->H);
    V3SuB(dipole->greenu[0], dipole->greenu[2]);
    V3SuB(dipole->greenu[1], dipole->greenu[2]);
    V3AdD(dipole->greenu[0], dipole->green->H[0]);
    V3AdD(dipole->greenu[1], dipole->green->H[1]);
    dipole->green->volume = fabs(M3inv(dipole->green->H,dipole->green->HI));
    simple_D (dipole->Hfree, dipole->green->H, dipole->green->D);
    M3SYMMETRIZE(dipole->green->D, dipole->green->simpleStrain);

    M3SUB(dipole->green->simpleStrain, dipole->plasticStrain, TMP);
    ElasticityT4MulT2 (dipole->C, TMP, dipole->green->Stresscell);
    dipole->Egreenwork = ElasticityT2Contraction
        (dipole->green->Stresscell, dipole->green->simpleStrain) / 2 *
        dipole->volumefree;
    dipole->green->Ebench += dipole->Egreenwork;
    /* to obtain P-K Force0, Force1 in Green's supercell */
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, dipole->green);

    ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown(dipole->stay);
    M3EQV(dipole->Hfree, dipole->stay->H);
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, dipole->stay);

    ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown(dipole->z0);
    M3ZERO(dipole->z0->Stress0);
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, dipole->z0);

    ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown(dipole->z1);
    M3ZERO(dipole->z1->Stress1);
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, dipole->z1);

    ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown
        (dipole->natural);
    M3ZERO(dipole->natural->Stresscell);
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, dipole->natural);

    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleAssign() */


/* Reduced coordinate transformation */
void ElasticityPBCStaticStraightDislocationDipoleTransform
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, V3 snew)
{
    register int i, j;
    V3 xfree, xnew;
    V3 tmp, u, R;
    V3mM3(sfree, dipole->Hfree, xfree);
    V3EQV(xfree, xnew);
    for (i=-dipole->images[0]; i<=dipole->images[0]; i++)
        for (j=-dipole->images[1]; j<=dipole->images[1]; j++)
        {
            V3ADDMULMUL(i, dipole->Hfree[0], j, dipole->Hfree[1], R);
            V3SUBSUB(xfree, dipole->x0, R, tmp);
            ElasticityStraightDislocation_u (dipole->d, tmp, u);
            V3AdD (u, xnew);
            V3SuB (tmp, dipole->dx);
            ElasticityStraightDislocation_u (dipole->d, tmp, u);
            V3SuB (xnew, u);
        }
    V3SuB(xnew, dipole->greenu[2]);
    V3mM3(xnew, dipole->green->HI, snew);
    return;
} /* ElasticityPBCStaticStraightDislocationDipoleTransform() */


/* Local strain in the Green's supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleGreenStrain
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, T2 Strain)
{
    register int i, j;
    V3 xfree;
    V3 tmp, R;
    T2 tmpstrain;
    V3mM3(sfree, dipole->Hfree, xfree);
    M3ZERO(Strain);
    for (i=-dipole->images[0]; i<=dipole->images[0]; i++)
        for (j=-dipole->images[1]; j<=dipole->images[1]; j++)
        {
            V3ADDMULMUL(i, dipole->Hfree[0], j, dipole->Hfree[1], R);
            V3SUBSUB(xfree, dipole->x0, R, tmp);
            ElasticityStraightDislocation_strain (dipole->d, tmp, tmpstrain);
            M3AdD (tmpstrain, Strain);
            V3SuB (tmp, dipole->dx);
            ElasticityStraightDislocation_strain (dipole->d, tmp, tmpstrain);
            M3SuB (Strain, tmpstrain);
        }
    return;
} /* ElasticityPBCStaticStraightDislocationDipoleGreenStrain() */


/* Local stress in the Green's supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleGreenStress
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, T2 Stress)
{
    register int i, j;
    V3 xfree;
    V3 tmp, R;
    T2 tmpstress;
    V3mM3(sfree, dipole->Hfree, xfree);
    M3ZERO(Stress);
    for (i=-dipole->images[0]; i<=dipole->images[0]; i++)
        for (j=-dipole->images[1]; j<=dipole->images[1]; j++)
        {
            V3ADDMULMUL(i, dipole->Hfree[0], j, dipole->Hfree[1], R);
            V3SUBSUB(xfree, dipole->x0, R, tmp);
            ElasticityStraightDislocation_stress (dipole->d, tmp, tmpstress);
            M3AdD (tmpstress, Stress);
            V3SuB (tmp, dipole->dx);
            ElasticityStraightDislocation_stress (dipole->d, tmp, tmpstress);
            M3SuB (Stress, tmpstress);
        }
    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleGreenStress() */


/* Set all user-setup inputs as unknown */
void ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown
(ElasticityPBCStaticStraightDislocationDipoleUserSetup *user)
{
    ElasticityPBCStaticStraightDislocationDipoleUserSetup unknown[1] =
    {UnknownElasticityPBCStaticStraightDislocationDipoleUserSetup};
    MEMCpy(user, unknown,
           ElasticityPBCStaticStraightDislocationDipoleUserSetup);
    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown() */


/* Autocomplete user-setup */
void ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
(ElasticityPBCStaticStraightDislocationDipole *dipole, 
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user)
{
    M3 TMP,STRESS;
    if ( ElasticityM3Known(user->H) &&
         ElasticityM3Known(user->D) &&
         ElasticityM3Known(user->Stress0) &&
         ElasticityM3Known(user->Stress1) &&
         ElasticityM3Known(user->Stresscell) )
    {
        ElasticityStraightDislocation_pkforce
            (user->Stress0, dipole->d, user->Force0);
        ElasticityStraightDislocation_pkforce
            (user->Stress1, dipole->d, user->Force1);
        V3NeG(user->Force1);
        return;
    }

    if (ElasticityM3Known(user->H) ||
        ElasticityM3Known(user->D) )
    {
        if (ElasticityM3Known(user->H))
            simple_D (dipole->Hfree, user->H, user->D);
        else simple_deform (dipole->Hfree, user->D, user->H);
        M3SYMMETRIZE (user->D, user->simpleStrain);
        M3SUB(user->simpleStrain, dipole->green->simpleStrain, TMP);
        ElasticityT4MulT2 (dipole->C, TMP, STRESS);
        M3ADD(dipole->green->Stress0, STRESS, user->Stress0);
        M3ADD(dipole->green->Stress1, STRESS, user->Stress1);
        M3ADD(dipole->green->Stresscell, STRESS, user->Stresscell);
    }
    else if ( ElasticityM3Known(user->Stress0) ||
              ElasticityM3Known(user->Stress1) ||
              ElasticityM3Known(user->Stresscell) )
    {
        if (ElasticityM3Known(user->Stress0))
        {
            M3SUB(user->Stress0, dipole->green->Stress0, STRESS);
            M3ADD(dipole->green->Stress1, STRESS, user->Stress1);
            M3ADD(dipole->green->Stresscell, STRESS, user->Stresscell);
        }
        else if (ElasticityM3Known(user->Stress1))
        {
            M3SUB(user->Stress1, dipole->green->Stress1, STRESS);
            M3ADD(dipole->green->Stress0, STRESS, user->Stress0);
            M3ADD(dipole->green->Stresscell, STRESS, user->Stresscell);
        }
        else
        {
            M3SUB(user->Stresscell, dipole->green->Stresscell, STRESS);
            M3ADD(dipole->green->Stress0, STRESS, user->Stress0);
            M3ADD(dipole->green->Stress1, STRESS, user->Stress1);
        }
        ElasticityT4MulT2 (dipole->S, STRESS, TMP);
        M3ADD(dipole->green->D, TMP, user->D);
        /* we strive to maintain the dislocation-axis unchanged */
        user->D[0][1] += user->D[1][0]; user->D[1][0] = 0;
        user->D[0][2] += user->D[2][0]; user->D[2][0] = 0;
        user->D[1][2] += user->D[2][1]; user->D[2][1] = 0;
        simple_deform(dipole->Hfree, user->D, user->H);
        M3SYMMETRIZE (user->D, user->simpleStrain);
    }
    else pe ("ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete:\n"
             "unable to complete because not enough is Known.\n");

    user->volume = fabs(M3inv(user->H,user->HI));
    M3SUB(user->simpleStrain, dipole->green->simpleStrain, TMP);
    user->Ebench = dipole->green->Ebench +
        ElasticityT2ContractionTRAPEZOIDAL
        (dipole->green->Stresscell, user->Stresscell, TMP, STRESS) *
        dipole->volumefree;
    /* to obtain P-K Force0, Force1 */
    ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
        (dipole, user);
    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete() */


/* Strain field in the user-setup supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleUserStrain
(ElasticityPBCStaticStraightDislocationDipole *dipole,
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user,
 V3 sfree, T2 Strain)
{
    M3 STRAIN;
    M3SUB(user->simpleStrain, dipole->green->simpleStrain, STRAIN);
    ElasticityPBCStaticStraightDislocationDipoleGreenStrain
        (dipole, sfree, Strain);
    M3AdD (STRAIN, Strain);
    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleUserStrain() */


/* Stress field in the user-setup supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleUserStress
(ElasticityPBCStaticStraightDislocationDipole *dipole,
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user,
 V3 sfree, T2 Stress)
{
    M3 STRESS;
    M3SUB(user->Stresscell, dipole->green->Stresscell, STRESS);
    ElasticityPBCStaticStraightDislocationDipoleGreenStress
        (dipole, sfree, Stress);
    M3AdD (STRESS, Stress);
    return;
} /* end ElasticityPBCStaticStraightDislocationDipoleUserStress() */


#ifdef _PBCStaticStraightDislocationDipoleEgreenwork_TEST
/* This proves the boundary work integral _is_ Egreenwork */
static ElasticityPBCStaticStraightDislocationDipole *energyfluxdipole;
static double energyflux (double s)
{
    int i;
    double flux=0, sfree[3]={0};
    V3 surface, force;
    T2 stress;
    for (i=0; i<2; i++)
    {
        V3CROSS(energyfluxdipole->Hfree[i],
                energyfluxdipole->Hfree[2], surface);
        if (V3DOT(surface, energyfluxdipole->Hfree[1-i]) < 0) V3NeG(surface);
        sfree[i] = s; sfree[1-i] = 1;
        ElasticityPBCStaticStraightDislocationDipoleGreenStress
            (energyfluxdipole, sfree, stress);
        V3mM3(surface, stress, force);
        flux += V3DOT(force, energyfluxdipole->greenu[1-i])/2;
    }
    return(flux);
} /* end energyflux() */
int main (int argc, char *argv[])
{
    int i,images[][2]={{30,30},{30,100}};
    M3 H;
    ElasticityPBCStaticStraightDislocationDipole dipole[1];
    Elasticity e = UnknownElasticity;
    double dipole_s0[2], dipole_s1[2], dipole_b[3];

    M3diagonal(53.2122231030966,56.4400856973166,7.68052283318781,H);
    ElasticityCubicT4Assign(
        1.009057365175595e+00, 5.092577892287198e-01,
        3.761627302304579e-01, e.c);
    /* ElasticityIsotropicT4Assign( */
    /* 1.009057365175595e+00, 2.092577892287198e-01, e.c); */
    V3ASSIGN(1,1,-2,e.o[0]);
    V3ASSIGN(1,1, 1,e.o[1]);
    ElasticityObserverComplete(e.o);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    V2ASSIGN(0.25, 0.5,  dipole_s0);
    V2ASSIGN(0.75, 0.65, dipole_s1);
    V3ASSIGN(0, 0, H[2][2]/2, dipole_b);
    for (i=0; i<sizeof(images)/sizeof(int)/2; i++)
    {
        printf ("%d x %d Green's function images:\n",
                1+2*images[i][0], 1+2*images[i][1]);
        ElasticityPBCStaticStraightDislocationDipoleAssign
            (H, e.C, dipole_s0, dipole_s1, dipole_b, images[i], dipole);
        energyfluxdipole = dipole;
        printf ("Egreenwork = %g, Eboundarywork = %g\n",
                dipole->Egreenwork, qromb(energyflux,0.,1.));
    }
    return (0);
}
#endif /* _PBCStaticStraightDislocationDipoleEgreenwork_TEST */


#ifdef _PBCStaticStraightDislocationDipoleInvariant_TEST
/* Displacement field in the PBC supercell */
void ElasticityPBCStaticStraightDislocationDipole_U
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, V3 U)
{
    int i, j;
    V3 xfree;
    V3 tmp, u, R;
    V3mM3(sfree, dipole->Hfree, xfree);
    V3ZERO(U);
    for (i=-dipole->images[0]; i<=dipole->images[0]; i++)
        for (j=-dipole->images[1]; j<=dipole->images[1]; j++)
        {
            V3ADDMULMUL(i, dipole->Hfree[0], j, dipole->Hfree[1], R);
            V3SUBSUB(xfree, dipole->x0, R, tmp);
            ElasticityStraightDislocation_u (dipole->d, tmp, u);
            V3AdD (u, U);
            V3SuB (tmp, dipole->dx);
            ElasticityStraightDislocation_u (dipole->d, tmp, u);
            V3SuB (U, u);
        }
    return;
} /* ElasticityPBCStaticStraightDislocationDipole_U() */
#define MESH 10
int main (int argc, char *argv[])
{
    int i,j,images[][2]={{100,100},{100,300}};
    M3 H;
    ElasticityPBCStaticStraightDislocationDipole dipole[1];
    Elasticity e = UnknownElasticity;
    double dipole_s0[2], dipole_s1[2], dipole_b[3];
    V3 s, U0, U1;
    M3 stress;
    
    M3diagonal(53.2122231030966,56.4400856973166,7.68052283318781,H);
    ElasticityCubicT4Assign(
        1.009057365175595e+00, 5.092577892287198e-01,
        3.761627302304579e-01, e.c);
    /* ElasticityIsotropicT4Assign( */
    /* 1.009057365175595e+00, 2.092577892287198e-01, e.c); */
    V3ASSIGN(1,1,-2,e.o[0]);
    V3ASSIGN(1,1, 1,e.o[1]);
    ElasticityObserverComplete(e.o);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    S6PR("C = %M\n",e.C);
    V2ASSIGN(0.25, 0.5, dipole_s0);
    V2ASSIGN(0.75, 0.65, dipole_s1);
    V3ASSIGN(0, 0, H[2][2]/2, dipole_b);
    for (i=0; i<sizeof(images)/sizeof(int)/2; i++)
    {
        printf ("%d x %d Green's function images:\n",
                1+2*images[i][0], 1+2*images[i][1]);
        ElasticityPBCStaticStraightDislocationDipoleAssign
            (H, e.C, dipole_s0, dipole_s1, dipole_b, images[i], dipole);
        printf ("Egreen = %g, Estay = %g, Ez0 = %g, Emin = %g\n",
                dipole->green->Ebench,
                dipole->stay->Ebench,
                dipole->z0->Ebench,
                dipole->natural->Ebench);
        V3ASSIGN(0., 0.3, 0.3, s);
        ElasticityPBCStaticStraightDislocationDipoleUserStress
            (dipole, dipole->z0, s, stress);
        S3PR("%M\n", stress);
        s[0] = 1;
        ElasticityPBCStaticStraightDislocationDipoleUserStress
            (dipole, dipole->z0, s, stress);
        S3PR("%M\n", stress);
        for (j=0; j<MESH; j++)
        {
            s[0] = (j+0.5) / MESH; s[1] = 0.; s[2] = 0.34;
            ElasticityPBCStaticStraightDislocationDipole_U(dipole,s,U0);
            s[1]++; s[2] = 0.74;
            ElasticityPBCStaticStraightDislocationDipole_U(dipole,s,U1);
            V3SuB(U1,U0);
            V3pr("dU = %M\n",U1);
        }
        cr();
    }
    return (0);
}
#endif /* _PBCStaticStraightDislocationDipoleInvariant_TEST */
