#include "Elasticity.h"

/**********************************************************/
/* Straight dislocation in infinite linear elastic medium */
/**********************************************************/

/* Stroh solution of an infinite, straight, uniformly moving dislocation */
void ElasticityStraightMovingDislocation
(V3 m, V3 xi, V3 b, double mass_density, V3 velocity, T4 c,
 ElasticityStraightDislocation *d)
{
    int a,i,j,k,l;
    double ccr,cci,cc,ddr,ddi;
    M3 nn,nni,nm,mn,mm,A1,A2,A3,A4;
    M6 N,VR,VI;
    V6 wr,wi;

    V3direction ( xi, d->xi );
    V3CROSS ( d->xi, m, d->n );
    V3normalize ( d->n );
    V3CROSS ( d->n, d->xi, d->m );
    V3EQV (b, d->b);
    d->mass_density = mass_density;
    V3EQV (velocity, d->velocity);
    /* the velocity component along xi has no meaning */
    ccr = V3DOT (d->velocity, d->xi);
    V3SUBmuL (d->velocity, ccr, d->xi);
    ElasticityT4Eqv (c, d->c);
    ElasticityT4Complete (d->c);

    /* To break the isotropic degeneracy */
    ElasticityENSUREAnisotropy
        (d->c,ElasticityStraightDislocationMinAnisotropy);
    
    M3Assign (d->m, d->n, d->xi, d->observer);
    V3AbsoluteToObserver(d->observer,d->b,d->B);
    V3AbsoluteToObserver(d->observer,d->velocity,d->Velocity);
    ElasticityT4AbsoluteToObserver (d->observer, d->c, d->C);

    ElasticityStrohCONTRACTION
        (d->n,d->n, d->c, d->mass_density,d->velocity, nn);
    M3inv (nn, nni);
    ElasticityStrohCONTRACTION
        (d->n,d->m, d->c, d->mass_density,d->velocity, nm);
    /* velocity term breaks full elastic constant symmetry */
    M3TRANSPOSE (nm, mn);
    /* but this symmetry maintains */
    ElasticityStrohCONTRACTION
        (d->m,d->m, d->c, d->mass_density,d->velocity, mm);
    M3NEG (nni, A2);          /* upper right block */
    M3MUL (A2, nm, A1);       /* upper left block  */
    M3MUL2 (mn, A2, nm, A3, A4);
    M3AdD (mm, A3);           /* lower left block  */
    /* Hirth & Lothe pp. 468 contains a typo, should be -(mn)(nn)^-1 */
    M3MUL (mn, A2, A4);       /* lower right block */
    /* see D.J. Bacon, D.M. Barnett, R.O. Scattergood, */
    /* Progress in Materials Science 23 (1978) 53-262. pp. 123. */
    for (a=0; a<3; a++)
    { /* assemble 6x6 matrix N */
        V3EQV(A1[a], N[a]);
        V3EQV(A2[a], N[a]+3);
        V3EQV(A3[a], N[a+3]);
        V3EQV(A4[a], N[a+3]+3);
    }
    ElasticityStrohM6Diag(N, wr, wi, VR, VI);
    for (a=0; a<6; a++)
        if (wi[a]==0)
            pe("ElasticityStraightMovingDislocation: prescribed\n"
               "dislocation is supersonic or transonic, Baby!\n");
    for (a=0; a<6; a+=2)
        for (j=0; j<3; j++)
        {
            d->pr[a/2] = wr[a];
            d->pi[a/2] = wi[a];
            d->ar[a/2][j] = VR[j][a];
            d->ai[a/2][j] = VI[j][a];
            d->lr[a/2][j] = VR[j+3][a];
            d->li[a/2][j] = VI[j+3][a];
        }
    for (a=0; a<3; a++)
    {
        ccr = V3DOT(d->ar[a],d->lr[a]) - V3DOT(d->ai[a],d->li[a]);
        cci = V3DOT(d->ar[a],d->li[a]) + V3DOT(d->ai[a],d->lr[a]);
        cc = sqrt( SQUARE(ccr) + SQUARE(cci) );
        cci = atan2(cci, ccr);
        ccr = sqrt(.5 / cc) * cos(- cci / 2);
        cci = sqrt(.5 / cc) * sin(- cci / 2);
        for (j=0; j<3; j++)
        {
            ddr = d->ar[a][j];
            ddi = d->ai[a][j];
            d->ar[a][j] = ddr * ccr - ddi * cci;
            d->ai[a][j] = ddi * ccr + ddr * cci;
            ddr = d->lr[a][j];
            ddi = d->li[a][j];
            d->lr[a][j] = ddr * ccr - ddi * cci;
            d->li[a][j] = ddi * ccr + ddr * cci;
        }
    }

    /* Hirth & Lothe: pp. 470: Eqn. 13-177. Verify skew-orthogonality */
    /* for eigenvectors with DIFFERENT eigenvalues: */
    /* for (a=0; a<3; a++) */
    /* for (i=0; i<=a; i++) */
    /* printf ("%e %e\n", */
    /* V3DOT(d->ar[a],d->lr[i])-V3DOT(d->ai[a],d->li[i])+ */
    /* V3DOT(d->ar[i],d->lr[a])-V3DOT(d->ai[i],d->li[a]), */
    /* V3DOT(d->ar[a],d->li[i])+V3DOT(d->ai[a],d->lr[i])+ */
    /* V3DOT(d->ar[i],d->li[a])+V3DOT(d->ai[i],d->lr[a])); */

    M3ZERO(d->k);
    for (a=0; a<3; a++)
    {   /* Hirth & Lothe: pp. 471: Eqn. 13-189 */
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
            {
                d->kmode[a][i][j] = - 2 * (d->lr[a][i] * d->li[a][j] +
                                           d->li[a][i] * d->lr[a][j]);
                d->k[i][j] += d->kmode[a][i][j];
            }
        /* H&L pp. 470: Eqn. 13-182 */
        d->dr[a] = V3DOT(d->lr[a], d->b);
        d->di[a] = V3DOT(d->li[a], d->b);
        /* displacement field constant: H&L pp. 471: Eqn. 13-183 */
        for (k=0; k<3; k++)
        {
            d->ur[a][k] =
                ( d->ar[a][k] * d->dr[a] - d->ai[a][k] * d->di[a] ) / PI;
            d->ui[a][k] =
                ( d->ar[a][k] * d->di[a] + d->ai[a][k] * d->dr[a] ) / PI;
        }
        V3AbsoluteToObserver(d->observer,d->ur[a],d->UR[a]);
        V3AbsoluteToObserver(d->observer,d->ui[a],d->UI[a]);
        /* displacement gradient field constant: H&L pp. 471: Eqn. 13-184 */
        for (l=0; l<3; l++)
        {
            ccr = d->m[l] + d->pr[a] * d->n[l];
            cci =           d->pi[a] * d->n[l];
            for (k=0; k<3; k++)
            {
                d->uxr[a][k][l] = d->ur[a][k] * ccr - d->ui[a][k] * cci;
                d->uxi[a][k][l] = d->ui[a][k] * ccr + d->ur[a][k] * cci;
            }
        }
        ElasticityT2AbsoluteToObserver(d->observer,d->uxr[a],d->UXR[a]);
        ElasticityT2AbsoluteToObserver(d->observer,d->uxi[a],d->UXI[a]);
        /* velocity field constant: u = u(x-vt) */
        for (k=0; k<3; k++)
        {
            M3mV3 (d->uxr[a], d->velocity, d->udotr[a]);
            V3NeG (d->udotr[a]);
            M3mV3 (d->uxi[a], d->velocity, d->udoti[a]);
            V3NeG (d->udoti[a]);
        }
        V3AbsoluteToObserver(d->observer,d->udotr[a],d->UdotR[a]);
        V3AbsoluteToObserver(d->observer,d->udoti[a],d->UdotI[a]);
        /* strain field constant: */
        for (k=0; k<3; k++)
            for (l=0; l<=k; l++)
            {
                d->strainr[a][l][k] = d->strainr[a][k][l]
                    = 0.5 * (d->uxr[a][k][l] + d->uxr[a][l][k]);
                d->straini[a][l][k] = d->straini[a][k][l]
                    = 0.5 * (d->uxi[a][k][l] + d->uxi[a][l][k]);
            }
        ElasticityT2AbsoluteToObserver
            (d->observer,d->strainr[a],d->StrainR[a]);
        ElasticityT2AbsoluteToObserver
            (d->observer,d->straini[a],d->StrainI[a]);
        /* stress field constant: H&L pp. 471: Eqn. 13-185 */
        ElasticityT4MulT2(d->c,d->strainr[a],d->stressr[a]);
        ElasticityT4MulT2(d->c,d->straini[a],d->stressi[a]);
        ElasticityT2AbsoluteToObserver
            (d->observer,d->stressr[a],d->StressR[a]);
        ElasticityT2AbsoluteToObserver
            (d->observer,d->stressi[a],d->StressI[a]);
    }
    ElasticityT2AbsoluteToObserver(d->observer,d->k,d->K);
    V3ZERO(d->PEmode);
    for (a=0; a<3; a++)
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                d->PEmode[a] += d->b[i] * d->kmode[a][i][j] * d->b[j] / 4 / PI;
    d->PE = V3SUM(d->PEmode);

    d->KE = 0;
    for (a=0; a<3; a++)
        for (l=0; l<3; l++)
        {
            ccr = d->pr[a] - d->pr[l];
            cci = d->pi[a] + d->pi[l];
            cc = SQUARE(ccr) + SQUARE(cci);
            ccr /= cc;
            cci /= -cc;
            ddr = V3DOT(d->udotr[a],d->udotr[l]) +
                V3DOT(d->udoti[a],d->udoti[l]);
            ddi = V3DOT(d->udotr[a],d->udoti[l]) -
                V3DOT(d->udoti[a],d->udotr[l]);
            d->KE -= PI * d->mass_density * (ccr * ddi + cci * ddr);
        }
    d->PE += d->KE;
    d->E = d->PE + d->KE;
    return;
} /* end ElasticityStraightMovingDislocation() */


/* Stroh solution of an infinite, straight, static dislocation */
void ElasticityStraightStaticDislocation
(V3 m, V3 xi, V3 b, T4 c, ElasticityStraightDislocation *d)
{
    V3 velocity = {0.,0.,0.};
    ElasticityStraightMovingDislocation (m, xi, b, 0., velocity, c, d);
    return;
} /* end ElasticityStraightStaticDislocation() */


/* Among np atoms s[0..3*np-1], pick out the atom w/ the smallest radial */
/* distance to dislocation at (s0,s1,0) which is also parallel to H[2][] */
int ElasticityStraightDislocationSmallestCartesianR2Atom
(M3 H, double s0, double s1, int np, double *s)
{
    register int i,j;
    V3 xi, dx;
    double cc, dmin=DOUBLE_PRECISION_INFINITY;
    V3DIRECTION(H[2], xi, cc);
    for (i=np; i--;)
    {
        V3ADDMULMUL( s[3*i]-s0, H[0], s[3*i+1]-s1, H[1], dx );
        cc = V3LENGTH2(dx) - SQUARE(V3DOT(xi,dx));
        if (cc < dmin)
        {
            dmin = cc;
            j = i;
        }
    }
    return(j);
} /* end ElasticityStraightDislocationSmallestCartesianR2Atom() */


/* Scaled distance squared for a pair of Stroh modes */
/* like the denominator in H&L pp. 455: Eqn. 13-129. */
double ElasticityStraightDislocationStrohR2
(ElasticityStraightDislocation *d, int a, V3 dx)
{
    double X,Y,etar,etai;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    etar = X + d->pr[a] * Y;
    etai =     d->pi[a] * Y;
    return(SQUARE(etar) + SQUARE(etai));
} /* end ElasticityStraightDislocationStrohR2() */


/* Peach-Koehler force[] per unit-length of dislocation in */
/* absolute frame at external stress[][] in absolute frame */
void ElasticityStraightDislocation_pkforce
(T2 stress, ElasticityStraightDislocation *d, V3 force)
{
    V3 tmp;
    ElasticityStraightDislocation_pkFORCE(stress,d,force,tmp);
    return;
} /* end ElasticityStraightDislocation_pkforce() */


/* P-K force per unit-length of dislocation resolved in normalized  */
/* absolute direction[] at an external stress[][] in absolute frame */
double ElasticityStraightDislocation_pkresolved
(T2 stress, ElasticityStraightDislocation *d, V3 direction)
{
    V3 tmp, force;
    ElasticityStraightDislocation_pkFORCE(stress,d,force,tmp);
    return(V3DOT(force,direction));
} /* end ElasticityStraightDislocation_pkresolved() */


/* Peach-Koehler Force[] per unit-length of dislocation in */
/* observer frame at external stress[][] in absolute frame */
void ElasticityStraightDislocation_pkForce
(T2 stress, ElasticityStraightDislocation *d, V3 Force)
{
    V3 tmp, force;
    ElasticityStraightDislocation_pkFORCE(stress,d,force,tmp);
    V3AbsoluteToObserver(d->observer,force,Force);
    return;
} /* end ElasticityStraightDislocation_pkForce() */


/* P-K force per unit-length of dislocation resolved in normalized  */
/* observer Direction[] at an external stress[][] in absolute frame */
double ElasticityStraightDislocation_pkResolved
(T2 stress, ElasticityStraightDislocation *d, V3 Direction)
{
    V3 tmp, force, Force;
    ElasticityStraightDislocation_pkFORCE(stress,d,force,tmp);
    V3AbsoluteToObserver(d->observer,force,Force);
    return(V3DOT(Force,Direction));
} /* end ElasticityStraightDislocation_pkResolved() */


/* Peach-Koehler Force[] per unit-length of dislocation in */
/* observer frame at external Stress[][] in observer frame */
void ElasticityStraightDislocation_PKForce
(T2 Stress, ElasticityStraightDislocation *d, V3 Force)
{
    ElasticityStraightDislocation_PKFORCE(Stress,d,Force);
    return;
} /* end ElasticityStraightDislocation_PKForce() */


/* P-K force per unit-length of dislocation resolved in normalized  */
/* observer Direction[] at an external Stress[][] in observer frame */
double ElasticityStraightDislocation_PKResolved
(T2 Stress, ElasticityStraightDislocation *d, V3 Direction)
{
    V3 Force;
    ElasticityStraightDislocation_PKFORCE(Stress,d,Force);
    return(V3DOT(Force,Direction));
} /* end ElasticityStraightDislocation_PKResolved() */


/* Peach-Koehler force[] per unit-length of dislocation in */
/* absolute frame at external Stress[][] in observer frame */
void ElasticityStraightDislocation_PKforce
(T2 Stress, ElasticityStraightDislocation *d, V3 force)
{
    V3 Force;
    ElasticityStraightDislocation_PKFORCE(Stress,d,Force);
    V3ObserverToAbsolute(d->observer,Force,force);
    return;
} /* end ElasticityStraightDislocation_PKforce() */


/* P-K force per unit-length of dislocation resolved in normalized  */
/* absolute direction[] at an external Stress[][] in observer frame */
double ElasticityStraightDislocation_PKresolved
(T2 Stress, ElasticityStraightDislocation *d, V3 direction)
{
    V3 Force, force;
    ElasticityStraightDislocation_PKFORCE(Stress,d,Force);
    V3ObserverToAbsolute(d->observer,Force,force);
    return(V3DOT(force,direction));
} /* end ElasticityStraightDislocation_PKresolved() */


/* Self-energy of a straight static dislocation dipole */
/* separated by dx[] in infinite linear elastic medium */
double ElasticityStraightDislocationStaticDipoleSelfEnergy
(ElasticityStraightDislocation *d, V3 dx)
{
    return( d->PE * log(ElasticityStraightDislocationCartesianR2(d,dx)) );
    /************************************************************/
    /* According to H&L pp. 474, Axiom 13-1, PE is invariant    */
    /* to m; and there is Eself = 2*Ecore + 2*PE*log(R/r0),     */
    /* where PE=b'*k*b/4/PI. For screw dislocation in isotropic */
    /* medium, k=mu. Here we compute only the 2*PE*log(R) part. */
    /************************************************************/
} /* end ElasticityStraightDislocationStaticDipoleSelfEnergy() */


/* Eself(dipole) = 2*Ecore + 2*PE*log(R/r0) + 2*A(theta): return A(theta); */
/* theta is with respect to input dx[], and xi of d, so A(theta=0) = 0.    */
double ElasticityStraightDislocationStaticDipoleSelfAngularEnergy
(ElasticityStraightDislocation *d, V3 dx, double theta)
{
    register int a;
    double AngularEnergy = 0, dxrotated[3];
    V3axialrotate (d->xi, theta, dx, dxrotated);
    for (a=0; a<3; a++)
        AngularEnergy += d->PEmode[a] *
            log( ElasticityStraightDislocationStrohR2(d,a,dxrotated) /
                 ElasticityStraightDislocationStrohR2(d,a,dx) );
    return (AngularEnergy / 2);
} /* end ElasticityStraightDislocationStaticDipoleSelfAngularEnergy() */


/* Total coupling energy between two identical straight */
/* static dislocation dipoles of dx[] separated by R[]. */
double ElasticityStraightDislocationStaticDipoleInteractionEnergy
(ElasticityStraightDislocation *d, V3 dx, V3 R)
{
    register int a;
    double pote=0, R2;
    V3 Rp,Rm;
    for (a=0; a<3; a++)
    {
        V3ADD(R, dx, Rp);
        V3SUB(R, dx, Rm);
        R2 = ElasticityStraightDislocationStrohR2(d,a,R);
        pote += d->PEmode[a] *
            log( ElasticityStraightDislocationStrohR2(d,a,Rp) *
                 ElasticityStraightDislocationStrohR2(d,a,Rm) / R2 / R2);
    }
    return(pote);
} /* end ElasticityStraightDislocationStaticDipoleInteractionEnergy() */


/* displacement field: dx[] (input) and u[] are in absolute frame */
void ElasticityStraightDislocation_u
(ElasticityStraightDislocation *d, V3 dx, V3 u)
{
    int a,k;
    double X,Y,etar,etai,lnetar,lnetai;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    V3ZERO(u);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        lnetar = log(SQUARE(etar) + SQUARE(etai)) / 2;
        /* -pi at y=0+, pi at y=0- */
        lnetai = atan2(-etai,-etar);
        for (k=0; k<3; k++)
            u[k] += d->ur[a][k] * lnetai + d->ui[a][k] * lnetar;
    }
    return;
} /* end ElasticityStraightDislocation_u() */


/* Displacement field: (X,Y,0) input and U[] are in (m,n,xi) frame */
void ElasticityStraightDislocation_U
(ElasticityStraightDislocation *d, double X, double Y, V3 U)
{
    int a,k;
    double etar,etai,lnetar,lnetai;
    V3ZERO(U);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        lnetar = log(SQUARE(etar) + SQUARE(etai)) / 2;
        /* -pi at y=0+, pi at y=0- */
        lnetai = atan2(-etai,-etar);
        for (k=0; k<3; k++)
            U[k] += d->UR[a][k] * lnetai + d->UI[a][k] * lnetar;
    }
    return;
} /* end ElasticityStraightDislocation_U() */


/* displacement gradient field: dx[] (input) and */
/* ux[k][l] := du[k]/dx[l] are in absolute frame */
void ElasticityStraightDislocation_ux
(ElasticityStraightDislocation *d, V3 dx, double ux[3][3])
{
    int a,k,l;
    double X,Y,etar,etai,eta2;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    M3ZERO(ux);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<3; l++)
                ux[k][l] += d->uxr[a][k][l] * etai +
                    d->uxi[a][k][l] * etar;
    }
    return;
} /* end ElasticityStraightDislocation_ux */


/* Displacement Gradient field: (X,Y,0) input and */
/* UX[k][l] := dU[k]/dX[l] are in (m,n,xi) frame. */
void ElasticityStraightDislocation_UX
(ElasticityStraightDislocation *d, double X, double Y, double UX[3][3])
{
    int a,k,l;
    double etar,etai,eta2;
    M3ZERO(UX);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<3; l++)
                UX[k][l] += d->UXR[a][k][l] * etai +
                    d->UXI[a][k][l] * etar;
    }
    return;
} /* end ElasticityStraightDislocation_UX */


/* velocity field: dx[] (input) and udot[] are in absolute frame */
void ElasticityStraightDislocation_udot
(ElasticityStraightDislocation *d, V3 dx, V3 udot)
{
    int a,k;
    double X,Y,etar,etai,eta2;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    V3ZERO(udot);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            udot[k] += d->udotr[a][k] * etai + d->udoti[a][k] * etar;
    }
    return;
} /* end ElasticityStraightDislocation_udot() */


/* Velocity field: (X,Y,0) input and Udot[] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Udot
(ElasticityStraightDislocation *d, double X, double Y, V3 Udot)
{
    int a,k;
    double etar,etai,eta2;
    V3ZERO(Udot);
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            Udot[k] += d->UdotR[a][k] * etai + d->UdotI[a][k] * etar;
    }
    return;
} /* end ElasticityStraightDislocation_Udot() */


/* strain field: dx[] (input) and strain[][] are in absolute frame */
void ElasticityStraightDislocation_strain
(ElasticityStraightDislocation *d, V3 dx, T2 strain)
{
    int a,k,l;
    double X,Y,etar,etai,eta2;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    for (k=0; k<3; k++)
        for (l=0; l<=k; l++)
            strain[k][l] = 0;
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<=k; l++)
                strain[k][l] += d->strainr[a][k][l] * etai +
                    d->straini[a][k][l] * etar;
    }
    for (k=0; k<3; k++)
        for (l=0; l<k; l++)
            strain[l][k] = strain[k][l];
    return;
} /* end ElasticityStraightDislocation_strain */


/* Strain field: (X,Y,0) input and Strain[][] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Strain
(ElasticityStraightDislocation *d, double X, double Y, T2 Strain)
{
    int a,k,l;
    double etar,etai,eta2;
    for (k=0; k<3; k++)
        for (l=0; l<=k; l++)
            Strain[k][l] = 0;
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<=k; l++)
                Strain[k][l] += d->StrainR[a][k][l] * etai +
                    d->StrainI[a][k][l] * etar;
    }
    for (k=0; k<3; k++)
        for (l=0; l<k; l++)
            Strain[l][k] = Strain[k][l];
    return;
} /* end ElasticityStraightDislocation_Strain */


/* stress field: dx[] (input) and stress[][] are in absolute frame */
void ElasticityStraightDislocation_stress
(ElasticityStraightDislocation *d, V3 dx, T2 stress)
{
    int a,k,l;
    double X,Y,etar,etai,eta2;
    X = V3DOT(d->m,dx);
    Y = V3DOT(d->n,dx);
    for (k=0; k<3; k++)
        for (l=0; l<=k; l++)
            stress[k][l] = 0;
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<=k; l++)
                stress[k][l] += d->stressr[a][k][l] * etai +
                    d->stressi[a][k][l] * etar;
    }
    for (k=0; k<3; k++)
        for (l=0; l<k; l++)
            stress[l][k] = stress[k][l];
    return;
} /* end ElasticityStraightDislocation_stress */


/* Stress field: (X,Y,0) input and Stress[][] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Stress
(ElasticityStraightDislocation *d, double X, double Y, T2 Stress)
{
    int a,k,l;
    double etar,etai,eta2;
    for (k=0; k<3; k++)
        for (l=0; l<=k; l++)
            Stress[k][l] = 0;
    for (a=0; a<3; a++)
    {
        etar = X + d->pr[a] * Y;
        etai =     d->pi[a] * Y;
        eta2 = SQUARE(etar) + SQUARE(etai);
        etar /= eta2;
        etai /= -eta2;
        for (k=0; k<3; k++)
            for (l=0; l<=k; l++)
                Stress[k][l] += d->StressR[a][k][l] * etai +
                    d->StressI[a][k][l] * etar;
    }
    for (k=0; k<3; k++)
        for (l=0; l<k; l++)
            Stress[l][k] = Stress[k][l];
    return;
} /* end ElasticityStraightDislocation_Stress */


/* strain energy per unit volume: dx[] (input) is in absolute frame */
double ElasticityStraightDislocation_pote_density
(ElasticityStraightDislocation *d, V3 dx)
{
    T2 strain;
    ElasticityStraightDislocation_strain(d, dx, strain);
    return(ElasticityStrainEnergyDensity(d->c, strain));
} /* end ElasticityStraightDislocation_pote_density() */


/* Strain energy per unit volume: (X,Y,0) input is in (m,n,xi) frame */
double ElasticityStraightDislocation_Pote_Density
(ElasticityStraightDislocation *d, double X, double Y)
{
    T2 Strain;
    ElasticityStraightDislocation_Strain(d, X, Y, Strain);
    return(ElasticityStrainEnergyDensity(d->C, Strain));
} /* end ElasticityStraightDislocation_Pote_Density() */


/* kinetic energy per unit volume: dx[] (input) is in absolute frame */
double ElasticityStraightDislocation_kine_density
(ElasticityStraightDislocation *d, V3 dx)
{
    V3 udot;
    ElasticityStraightDislocation_udot(d, dx, udot);
    return(d->mass_density*V3LENGTH2(udot)/2);
} /* end ElasticityStraightDislocation_kine_density() */


/* Kinetic energy per unit volume: (X,Y,0) input is in (m,n,xi) frame */
double ElasticityStraightDislocation_Kine_Density
(ElasticityStraightDislocation *d, double X, double Y)
{
    V3 Udot;
    ElasticityStraightDislocation_Udot(d, X, Y, Udot);
    return(d->mass_density*V3LENGTH2(Udot)/2);
} /* end ElasticityStraightDislocation_Kine_Density() */


#ifdef _ElasticityStraightStaticDislocation_TEST
int main (int argc, char *argv[])
{
    int i;
    double c11=10.21,c12=3.34,c44=5.78,delta;
    Elasticity e = UnknownElasticity;
    ElasticityStraightDislocation d[1];
    V3 m,xi,U,x,xtest,utest;
    M3 dx,u,J;
    T2 stress,strain;
    
    V3ASSIGN(1,1,-2,e.o[0]);
    V3ASSIGN(1,1,1,e.o[1]);
    ElasticityObserverComplete(e.o);
    ElasticityCubicT4Assign(c11,c12,c44,e.c);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);

    V3ASSIGN(1,0,0,m);
    V3ASSIGN(0,0,1,xi);
    ElasticityStraightStaticDislocation(m,xi,xi,e.C,d);

    V3randomunit(x); delta = SMALL;
    for (i=0; i<3; i++)
    {
        V3randomunit(dx[i]);
        V3ADDMUL(x,delta,dx[i],xtest);
        ElasticityStraightDislocation_u(d,xtest,utest);
        V3EQV(utest,u[i]);
        V3SUBMUL(x,delta,dx[i],xtest);
        ElasticityStraightDislocation_u(d,xtest,utest);
        V3SuB(u[i],utest);
    }
    M3DividE(u,2*delta);
    M3Inv(dx);
    M3MUL(dx,u,J);
    M3Symmetrize(J);
    S3PR("strain (dif) = %M\n ", J);
    ElasticityStraightDislocation_strain(d,x,strain);
    S3PR("strain (ana) = %M\n ", strain);
    ElasticityT4MulT2 (e.C, J, stress);
    S3PR("stress (dif) = %M\n ", stress); 
    ElasticityStraightDislocation_stress(d,x,stress);
    S3PR("stress (ana) = %M\n ", stress);
    return (0);
}
#endif /* _ElasticityStraightStaticDislocation_TEST */


#ifdef _ElasticityStraightStaticIsotropicScrewDislocation_TEST
int main (int argc, char *argv[])
{
    T4 c;
    double delta,lambda=10.2,mu=3.5;
    double b=2,X=1.2,Y=0.4,EPSILONXZ,EPSILONYZ,SIGMAXZ,SIGMAYZ;
    ElasticityStraightDislocation d[1];
    V3 m,xi,U;
    T2 Stress,Strain;
    V3ASSIGN(1,0,0,m);
    V3ASSIGN(0,0,1,xi);
    V3NORMALIZE(xi,SIGMAXZ);
    V3MuL(b,xi);
    ElasticityIsotropicT4Assign(lambda,mu,c);
    ElasticityStraightStaticDislocation(m,xi,xi,c,d);
    ElasticityStraightDislocation_U(d,X,Y,U);
    printf("UZ = %e\n", b/PI/2*atan2(-Y,-X));
    V3pr("U = %M\n ", U);
    EPSILONXZ = -b/4/PI*Y/(SQUARE(X)+SQUARE(Y));
    EPSILONYZ =  b/4/PI*X/(SQUARE(X)+SQUARE(Y));
    printf("EPSILONXZ = %e,  EPSILONYZ = %e\n", EPSILONXZ, EPSILONYZ);
    ElasticityStraightDislocation_Strain(d,X,Y,Strain);
    S3PR("Strain = %M\n ", Strain);
    SIGMAXZ = -mu*b/2/PI*Y/(SQUARE(X)+SQUARE(Y));
    SIGMAYZ =  mu*b/2/PI*X/(SQUARE(X)+SQUARE(Y));
    printf("SIGMAXZ = %e,  SIGMAYZ = %e\n", SIGMAXZ, SIGMAYZ);
    ElasticityStraightDislocation_Stress(d,X,Y,Stress);
    S3PR("Stress = %M\n ", Stress);
    printf("pote density = %.15e = %.15e\n\n",
           EPSILONXZ*SIGMAXZ + EPSILONYZ*SIGMAYZ,
           ElasticityStraightDislocation_Pote_Density(d,X,Y));
    printf("energy prefactor = %.15e = %.15e\n",
           mu * SQUARE(b) / 4 / PI, d->E);
    return (0);
}
#endif /* _ElasticityStraightStaticIsotropicScrewDislocation_TEST */


#ifdef _ElasticityStraightMovingIsotropicScrewDislocation_TEST
int main (int argc, char *argv[])
{
    T4 c;
    double lambda=10.2,mu=3.5,mass_density=1.95435,v=1.03123;
    double b=2,X=1.2,Y=0.4,EPSILONXZ,EPSILONYZ,SIGMAXZ,SIGMAYZ;
    double Ct,gamma;
    ElasticityStraightDislocation d[1];
    V3 m,xi,U,velocity;
    T2 Stress,Strain;
    V3ASSIGN(1,0,0,m);
    V3ASSIGN(0,0,1,xi);
    V3NORMALIZE(xi,SIGMAXZ);
    V3MuL(b,xi);
    V3MUL(v,m,velocity);
    ElasticityIsotropicT4Assign(lambda,mu,c);
    ElasticityStraightMovingDislocation(m,xi,xi,mass_density,velocity,c,d);
    Ct = sqrt(mu/mass_density);
    gamma = sqrt(1-SQUARE(v/Ct));
    printf("Ct = %g,   gamma = %g\n\n", Ct,gamma);
    ElasticityStraightDislocation_U(d,X,Y,U);
    printf("UZ = %e\n", b/PI/2*atan2(-gamma*Y,-X));
    V3pr("U = %M\n ", U);
    SIGMAXZ = -mu*b/2/PI*gamma*Y/(SQUARE(X)+SQUARE(gamma*Y));
    SIGMAYZ =  mu*b/2/PI*gamma*X/(SQUARE(X)+SQUARE(gamma*Y));
    printf("SIGMAXZ = %e,  SIGMAYZ = %e\n", SIGMAXZ, SIGMAYZ);
    ElasticityStraightDislocation_Stress(d,X,Y,Stress);
    S3PR("Stress = %M\n ", Stress);
    printf("pote density = %.15e = %.15e\n\n",
           (SQUARE(SIGMAXZ)+SQUARE(SIGMAYZ))/2/mu,
           ElasticityStraightDislocation_Pote_Density(d,X,Y));
    printf("kinetic energy prefactor = %.15e\n", d->KE);
    printf("energy prefactor = %.15e = %.15e\n",
           mu * SQUARE(b) / 4 / PI / gamma, d->E);
    return (0);
}
#endif /* _ElasticityStraightMovingIsotropicScrewDislocation_TEST */


#ifdef _ElasticityStraightMovingDislocation_TEST
int main (int argc, char *argv[])
{
    T4 c;
    V3 xi,m,b,velocity;
    double mass_density=1.;
    ElasticityStraightDislocation d;
    double C11,C12,C13,C16,C22,C36,C44,C45,C55,C66,C11avg,H,Km,Kn,Kxi;
    /* values for aluminum */
    double c11=10.82,c12=6.13,c44=2.85;
    ElasticityCubicT4Assign(c11,c12,c44,c);
    /* ElasticityIsotropicT4Assign(5.93,2.65,C); */
    /* Hirth & Lothe pp. 461: Fig. 13-10 */
    V3ASSIGN(0,1,0,m);
    V3ASSIGN(1,0,1,xi);
    V3ASSIGN(1,1,-1,b);
    V3ASSIGN(0,0,0,velocity);
    ElasticityStraightMovingDislocation
        (m, xi, b, mass_density, velocity, c, &d);
    S3PR("observer(m,n,xi) = %M\n ", d.observer);
    printf ("Warning: only in this frame do c16,c26 vanish (H&L 13-99)\n"
            "upon which H&L 13-100 to 13-124 (pp. 449-454) are based,\n"
            "but this is not the frame we will likely use in practice!\n\n");
    S6PR("C(observer frame) = %M\n ", d.C);
    S3PR("K = %M\n ", d.K);
    C11 = d.C[0][0];
    C12 = d.C[0][1];
    C22 = d.C[1][1];
    C44 = d.C[3][3];
    C45 = d.C[3][4];
    C55 = d.C[4][4];
    C66 = d.C[5][5];
    C11avg = sqrt(C11*C22);
    Km = (C11avg + C12) *
        sqrt(C66 * (C11avg-C12) / (C11avg + C12 + 2*C66) / C22);
    Kn = (C11avg + C12) *
        sqrt(C66 * (C11avg-C12) / (C11avg + C12 + 2*C66) / C11);
    Kxi = sqrt(C44*C55 - C45*C45);
    printf ("Km=%g Kn=%g Kxi=%g\n\n", Km, Kn, Kxi);
    return (0);
}
#endif /* _ElasticityStraightMovingDislocation_TEST */
