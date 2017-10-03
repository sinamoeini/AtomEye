#include "Elasticity.h"

/*********************************/
/* Linear anisotropic elasticity */
/*********************************/

/* Voigt notation: 11->1, 22->2, 33->3, 23->4, 13->5, 12->6 */
const int ElasticityVoigtIndices[3][3] = {{0,5,4},{5,1,3},{4,3,2}};


/* Autocomplete orthonormal observer directions matrix */
void ElasticityObserverComplete (M3 observer)
{
    int i;
    double cc;
    bool known[3];
    for (i=0; i<3; i++)
        if (known[i]=ElasticityV3Known(observer[i]))
        {
            cc = V3LENGTH2(observer[i]);
            if (cc == 0.)
                pe ("ElasticityObserverComplete: declared\n"
                    "Known but is a zero-vector.\n");
            else if (cc != 1.)
            {
                cc = sqrt(cc);
                V3DiV(observer[i], cc);
            }
        }
    if ( known[0] && known[1] && known[2] )
        if ( M3NONORTHOGONAL(observer) )
        {
            S3PR("observer directions = %M -->", observer);
            pe ("ElasticityObserverComplete: declared\n"
                "Known but is not orthogonal\n");
        }
        else return;
    else if ( known[0] && known[1] )
    {
        V3CROSS(observer[0], observer[1], observer[2]);
        V3normalize(observer[2]);
        V3CROSS(observer[2], observer[0], observer[1]);
    }
    else if ( known[0] && known[2] )
    {
        V3CROSS(observer[2], observer[0], observer[1]);
        V3normalize(observer[1]);
        V3CROSS(observer[0], observer[1], observer[2]);
    }
    else if ( known[1] && known[2] )
    {
        V3CROSS(observer[1], observer[2], observer[0]);
        V3normalize(observer[0]);
        V3CROSS(observer[0], observer[1], observer[2]);
    }
    else
        pe ("ElasticityObserverComplete: unable\n"
            "to complete because not enough is Known.\n");
    return;
} /* end ElasticityObserverComplete() */

#ifdef _ElasticityObserverComplete_TEST
int main (int argc, char *argv[])
{
    Elasticity e = UnknownElasticity;
    /* V3ASSIGN(2,0,0,e.o[0]); */
    V3ASSIGN(0,4,0,e.o[1]);
    V3ASSIGN(0,0,1,e.o[2]);
    ElasticityObserverComplete(e.o);
    S3PR("observer directions = %M.", e.o);
    return (0);
}
#endif /* _ElasticityObserverComplete_TEST */


/* Autocomplete a rank-2 symmetric tensor */
void ElasticityT2Complete (T2 tensor)
{
    if ( ElasticityValueUnknown(tensor[0][0]) ||
         ElasticityValueUnknown(tensor[1][1]) ||
         ElasticityValueUnknown(tensor[2][2]) )
        pe ("ElasticityT2Complete: cannot complete\n"
            "because diagonal element is unknown.\n");
    if ( ElasticityValueKnown(tensor[0][1]) )
        tensor[1][0] = tensor[0][1];
    else if ( ElasticityValueKnown(tensor[1][0]) )
        tensor[0][1] = tensor[1][0];
    else pe ("ElasticityT2Complete: cannot complete\n"
             "because 12/21 element is unknown.\n");
    if ( ElasticityValueKnown(tensor[0][2]) )
        tensor[2][0] = tensor[0][2];
    else if ( ElasticityValueKnown(tensor[2][0]) )
        tensor[0][2] = tensor[2][0];
    else pe ("ElasticityT2Complete: cannot complete\n"
             "because 13/31 element is unknown.\n");
    if ( ElasticityValueKnown(tensor[1][2]) )
        tensor[2][1] = tensor[1][2];
    else if ( ElasticityValueKnown(tensor[2][1]) )
        tensor[1][2] = tensor[2][1];
    else pe ("ElasticityT2Complete: cannot complete\n"
             "because 23/32 element is unknown.\n");
    return;
} /* end ElasticityT2Complete() */


/* Convert a rank-2 tensor in Observer notation to Absolute notation */
void ElasticityT2ObserverToAbsolute
(M3 observer, T2 observer_notation, T2 absolute_notation)
{
    M3 OT,TMP;
    M3TRANSPOSE (observer, OT);
    M3MUL2 (OT, observer_notation, observer, absolute_notation, TMP);
    return;
} /* end ElasticityT2ObserverToAbsolute() */


/* Convert a rank-2 tensor in Absolute notation to Observer notation */
void ElasticityT2AbsoluteToObserver
(M3 observer, T2 absolute_notation, T2 observer_notation)
{
    M3 OT,TMP;
    M3TRANSPOSE (observer, OT);
    M3MUL2 (observer, absolute_notation, OT, observer_notation, TMP);
    return;
} /* end ElasticityT2AbsoluteToObserver() */

#ifdef _T2_TEST
int main (int argc, char *argv[])
{
    Elasticity e = UnknownElasticity;
    V3ASSIGN(1,1,0,e.o[0]);
    V3ASSIGN(-1,1,0,e.o[1]);
    ElasticityObserverComplete(e.o);
    S3PR("observer = %M.\n ", e.o);
    e.Stress[0][0] = 1;
    e.Stress[1][1] = -1;
    e.Stress[2][2] = 0;
    e.Stress[0][1] = e.Stress[0][2] = e.Stress[1][2] = 0;
    ElasticityT2Complete(e.Stress);
    S3PR("observer stress = %M,", e.Stress);
    printf ("Hydro Invariant = %f,\n", ElasticityT2HydroInvariant(e.Stress));
    printf ("Mises Invariant = %f.\n\n", ElasticityT2MisesInvariant(e.Stress));
    ElasticityT2ObserverToAbsolute(e.o, e.Stress, e.stress);
    S3PR("absolute stress = %M,", e.stress);
    printf ("Hydro Invariant = %f,\n", ElasticityT2HydroInvariant(e.stress));
    printf ("Mises Invariant = %f.\n", ElasticityT2MisesInvariant(e.stress));
    return (0);
}
#endif /* _T2_TEST */


#ifdef _T2AbsoluteToObserver
int main (int argc, char *argv[])
{
    int i;
    M3 I=IdentityM3, t, T;
    Elasticity e = UnknownElasticity;
    char buf[FGETS_MAXSIZE]={0};

    printf ("\nPlease input observer frame's\n");
    for (i=0; i<3; i++)
    {
        printf ("%s edge: ", word_for_order(i+1)); fflush(stdout);
        Gets(buf);
        if (*blank_advance(buf)!=EOS) sscanf(buf, "%lf %lf %lf", V3e(e.o[i]));
    }
    ElasticityObserverComplete(e.o);
    printf ("\nPlease input tensor in absolute frame:\n"); fflush(stdout);
    for (i=0; i<3; i++)
    {
        Gets(buf);
        sscanf(buf, "%lf %lf %lf", V3e(t[i]));
    }
    S3PR("\nIn the %M frame,\n ", I);
    S3PR("tensor = %M (input).\n ", t);
    ElasticityT2AbsoluteToObserver(e.o, t, T);
    S3PR("In the %M observer frame,\n ", e.o);
    S3PR("tensor = %M.\n ", T);
    return (0);
}
#endif /* _T2AbsoluteToObserver */


#ifdef _T2ObserverToAbsolute
int main (int argc, char *argv[])
{
    int i;
    M3 I=IdentityM3, t, T;
    Elasticity e = UnknownElasticity;
    char buf[FGETS_MAXSIZE]={0};

    printf ("\nPlease input observer frame's\n");
    for (i=0; i<3; i++)
    {
        printf ("%s edge: ", word_for_order(i+1)); fflush(stdout);
        Gets(buf);
        if (*blank_advance(buf)!=EOS) sscanf(buf, "%lf %lf %lf", V3e(e.o[i]));
    }
    ElasticityObserverComplete(e.o);
    printf ("\nPlease input tensor in observer frame:\n"); fflush(stdout);
    for (i=0; i<3; i++)
    {
        Gets(buf);
        sscanf(buf, "%lf %lf %lf", V3e(T[i]));
    }
    S3PR("\nIn the %M observer frame,\n ", e.o);
    S3PR("tensor = %M (input).\n ", T);
    ElasticityT2ObserverToAbsolute(e.o, T, t);
    S3PR("In the %M frame,\n ", I);
    S3PR("tensor = %M.\n ", t);
    return (0);
}
#endif /* _T2ObserverToAbsolute */


/* ElasticityT2Contraction((A+B)/2,C) */
double ElasticityT2ContractionTrapezoidal(M3 A, M3 B, M3 C)
{
    M3 TMP;
    return(ElasticityT2ContractionTRAPEZOIDAL(A,B,C,TMP));
} /* end ElasticityT2ContractionTrapezoidal() */


/* Autocomplete a rank-4 symmetric tensor */
void ElasticityT4Complete (T4 tensor)
{
    int i,j;
    if ( ElasticityValueUnknown(tensor[0][0]) ||
         ElasticityValueUnknown(tensor[1][1]) ||
         ElasticityValueUnknown(tensor[2][2]) ||
         ElasticityValueUnknown(tensor[3][3]) ||
         ElasticityValueUnknown(tensor[4][4]) ||
         ElasticityValueUnknown(tensor[5][5]) )
        pe ("ElasticityT4Complete: cannot complete\n"
            "because diagonal element is unknown.\n");
    for (i=0; i<5; i++)
        for (j=i+1; j<6; j++)
            if (ElasticityValueKnown(tensor[i][j]))
                tensor[j][i] = tensor[i][j];
            else if (ElasticityValueKnown(tensor[j][i]))
                tensor[i][j] = tensor[j][i];
            else pe ("ElasticityT4Complete: cannot complete\n"
                     "because %1d%1d/%1d%1d element is unknown.\n",
                     i+1,j+1, j+1,i+1);
    return;
} /* end ElasticityT4Complete() */


#define RANK 6
/* T_{ijkl} * TI_{klmn} = delta_{im} * delta_{jn}, but both */
/* T and TI are expressed in Voigt notation due to symmetry */
void ElasticityT4Inverse (T4 tensor, T4 inverse)
{
    int i,j,k,rank=RANK,ipiv[RANK],info;
    double work[RANK],ap[RANK*(RANK+1)/2];
    for (k=j=0; j<RANK; j++)
        for (i=0; i<=j; i++,k++)
            ap[k] = tensor[i][j];
    dsptrf_ ( "U", &rank, ap, ipiv, &info );
    dsptri_ ( "U", &rank, ap, ipiv, work, &info );
    if (info > 0)
        pe("ElasticityT4Inverse: the matrix is singular\n"
           "and its inverse could not be computed.\n");
    for (k=j=0; j<RANK; j++)
        for (i=0; i<=j; i++,k++)
        {
            inverse[i][j] = ap[k];
            if (i>2) inverse[i][j] /= 2;
            if (j>2) inverse[i][j] /= 2;
            inverse[j][i] = inverse[i][j];
        }
    return;
} /* end ElasticityT4Inverse() */
#undef RANK

#ifdef _ElasticityT4Inverse_TEST
int main (int argc, char *argv[])
{
    int i,j;
    T2 Stress=GuineaPigSymmetricM3, Strain;
    Elasticity e = UnknownElasticity;
    V3ASSIGN(1,-2,1,e.o[0]);
    V3ASSIGN(1,1,1,e.o[1]);
    ElasticityObserverComplete(e.o);
    ElasticityCubicT4Assign(1.68486797165, 1.21325117263,
                            0.75413697334, e.c);
    S6PR("c = %M\n ",e.c);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    S6PR("C = %M\n ",e.C);
    ElasticityT4Inverse (e.C,e.S);
    S6PR("S = %M\n ",e.S);
    S3PR("Stress = %M", Stress);
    ElasticityT4MulT2 (e.S,Stress,Strain);
    ElasticityT4MulT2 (e.C,Strain,Stress);
    S3PR("Stress = %M", Stress);
    return (0);
}
#endif /* _ElasticityT4Inverse_TEST */


#define RANK 6
/* Element-by-element reciprocal */
void ElasticityT4ElementalReciprocal (T4 tensor, T4 reciprocal)
{
    int i,j;
    for (j=0; j<RANK; j++)
        for (i=0; i<=j; i++)
        {
            reciprocal[i][j] = 1 / tensor[i][j];
            reciprocal[j][i] = reciprocal[i][j];
        }
    return;
} /* end ElasticityT4ElementalReciprocal() */
#undef RANK


#define RANK 6
/* Transform a direct-6x6-inverse, which is a non-tensor, to a tensor */
void ElasticityT4tildeToT4 (T4 tilde, T4 tensor)
{
    int i,j;
    for (j=0; j<RANK; j++)
        for (i=0; i<=j; i++)
        {
            tensor[i][j] = tilde[i][j];
            if (i>2) tensor[i][j] /= 2;
            if (j>2) tensor[i][j] /= 2;
            tensor[j][i] = tensor[i][j];
        }
    return;
} /* end ElasticityT4tildeToT4() */
#undef RANK


#ifdef _cubiC
int main (int argc, char *argv[])
{
    int i,j;
    M3 I=IdentityM3; 
    Elasticity e = UnknownElasticity;
    char buf[FGETS_MAXSIZE]={0};

    if (argc != 4)
    {
        printf ("\nPurpose: suppose the crystal has cubic symmetry with\n"
                "elastic constants c11,c12,c44 in the current frame,\n"
                "print out the 6x6 elastic constant and elastic compliance\n"
                "matrices in any new frame.\n\n");
        printf ("Usage: %s c11 c12 c44\n\n", argv[0]);
        return (1);
    }
    printf ("\nPlease input new frame's\n");
    for (i=0; i<3; i++)
    {
        printf ("%s edge: ", word_for_order(i+1)); fflush(stdout);
        Gets(buf);
        if (*blank_advance(buf)!=EOS) sscanf(buf, "%lf %lf %lf", V3e(e.o[i]));
    }
    ElasticityObserverComplete(e.o);
    ElasticityCubicT4Assign(atof(argv[1]), atof(argv[2]), atof(argv[3]), e.c);
    S3PR("\nIn the current %M frame,\n ", I);
    S6PR("c = %M.\n ",e.c);
    S3PR("In the new %M observation frame,\n ", e.o);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    S6PR("C = %M,\n ",e.C);
    ElasticityT4Inverse (e.C,e.S);
    S6PR("S = %M.\n ",e.S);
    ElasticityT4RelaxedModuli(e.S,e.M);
    S6PR("M = %M.\n ",e.M);
    return (0);
}
#endif /* _cubiC */


#ifdef _hCp
int main (int argc, char *argv[])
{
    int i,j;
    M3 I=IdentityM3; 
    Elasticity e = UnknownElasticity;
    char buf[FGETS_MAXSIZE]={0};
    if (argc != 6)
    {
        printf ("\nPurpose: a hexagonal-close-packed crystal has elastic\n"
                "constants c11,c12,c13,c33,c44 (basal plane 1-2 invariant\n"
                "to axis-3 rotation). Print out the 6x6 elastic constant\n"
                "and elastic compliance matrices in any new frame.\n\n");
        printf ("Usage: %s c11 c12 c13 c33 c44\n\n", argv[0]);
        return (1);
    }
    printf ("\nPlease input new frame's\n");
    for (i=0; i<3; i++)
    {
        printf ("%s edge: ", word_for_order(i+1)); fflush(stdout);
        Gets(buf);
        if (*blank_advance(buf)!=EOS) sscanf(buf, "%lf %lf %lf", V3e(e.o[i]));
    }
    ElasticityObserverComplete(e.o);
    ElasticityPlanarIsotropyT4Assign
        (atof(argv[1]), atof(argv[2]), atof(argv[3]),
         atof(argv[4]), atof(argv[5]), e.c);
    S3PR("\nIn the current %M frame,\n ", I);
    S6PR("c = %M.\n ",e.c);
    S3PR("In the new %M observation frame,\n ", e.o);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    S6PR("C = %M,\n ",e.C);
    ElasticityT4Inverse (e.C,e.S);
    S6PR("S = %M.\n ",e.S);
    ElasticityT4RelaxedModuli(e.S,e.M);
    S6PR("M = %M.\n ",e.M);
    return (0);
}
#endif /* _hCp */


/* Convert a rank-4 tensor in Observer notation to Absolute notation */
void ElasticityT4ObserverToAbsolute
(M3 observer, T4 observer_notation, T4 absolute_notation)
{
    int i,j,k,l, ip,jp,kp,lp, a,b;
    for (i=0; i<3; i++)
        for (j=i; j<3; j++)
        {
            a = VoigtIndex(i,j);
            for (k=0; k<3; k++)
                for (l=k; l<3; l++)
                {
                    b = VoigtIndex(k,l);
                    absolute_notation[a][b] = 0;
                    /* pgcc bug here if "-funroll-loops" flag is turned on */
                    for (ip=0; ip<3; ip++)
                        for (jp=0; jp<3; jp++)
                            for (kp=0; kp<3; kp++)
                                for (lp=0; lp<3; lp++)
                                    absolute_notation[a][b] +=
                                        T4Element
                                        (observer_notation,ip,jp,kp,lp) *
                                        observer[ip][i] *
                                        observer[jp][j] *
                                        observer[kp][k] *
                                        observer[lp][l];
                    absolute_notation[b][a] = absolute_notation[a][b];
                }
        }                    
    return;
} /* end ElasticityT4ObserverToAbsolute() */


/* Convert a rank-4 tensor in Absolute notation to Observer notation */
void ElasticityT4AbsoluteToObserver
(M3 observer, T4 absolute_notation, T4 observer_notation)
{
    int i,j,k,l, ip,jp,kp,lp, a,b;
    for (i=0; i<3; i++)
        for (j=i; j<3; j++)
        {
            a = VoigtIndex(i,j);
            for (k=0; k<3; k++)
                for (l=k; l<3; l++)
                {
                    b = VoigtIndex(k,l);
                    observer_notation[a][b] = 0;
                    /* pgcc bug here if "-funroll-loops" flag is turned on */
                    for (ip=0; ip<3; ip++)
                        for (jp=0; jp<3; jp++)
                            for (kp=0; kp<3; kp++)
                                for (lp=0; lp<3; lp++)
                                    observer_notation[a][b] +=
                                        T4Element
                                        (absolute_notation,ip,jp,kp,lp) *
                                        observer[i][ip] *
                                        observer[j][jp] *
                                        observer[k][kp] *
                                        observer[l][lp];
                    observer_notation[b][a] = observer_notation[a][b];
                }
        }
    return;
} /* end ElasticityT4AbsoluteToObserver() */


/* What happens to a tensor, previously denoted "prev_notation" */
/* while an axis is labeled "prev_label", has the label changed */
/* to "current_label". The unrelated axis keeps its label.      */
void ElasticityT4Relabel
(T4 prev_notation, int prev_label, int current_label, T4 current_notation)
{
    int invariant_label;
    M3 observer = UnknownElasticityM3;
    T4 absolute_notation;
    if (prev_label == current_label)
    { /* nothing is switched */
        ElasticityT4Eqv (prev_notation, current_notation);
        return;
    }
    ElasticityT4Eqv (prev_notation, absolute_notation);
    invariant_label = 0 + 1 + 2 - prev_label - current_label;
    V3UNIT (observer[invariant_label],invariant_label);
    V3UNIT (observer[current_label],prev_label);
    ElasticityObserverComplete (observer);
    ElasticityT4AbsoluteToObserver
        (observer,absolute_notation,current_notation);
    return;
} /* end ElasticityT4Relabel() */


#ifdef _ElasticityT4_TEST
/* fcc example in Hirth & Lothe, pp. 434 */
int main (int argc, char *argv[])
{
    double t11=10.21,t12=3.34,t44=5.78;
    double T11,T12,T13,T16,T22,T44,T55,H;
    Elasticity e = UnknownElasticity;
    H = 2*t44 + t12 - t11;
    T11 = t11 + H / 2;
    T12 = t12 - H / 3;
    T13 = t12 - H / 6;
    T16 = H * SQRT2 / 6;
    T22 = t11 + H / 3 * 2;
    T44 = t44 - H / 3;
    T55 = t44 - H / 6;
    V3ASSIGN(1,-2,1,e.o[0]);
    V3ASSIGN(1,1,1,e.o[1]);
    ElasticityObserverComplete(e.o);
    S3PR("observer = %M\n ", e.o);
    ElasticityCubicT4Assign(t11,t12,t44,e.c);
    S6PR("c = %M\n ",e.c);
    ElasticityT4AbsoluteToObserver(e.o,e.c,e.C);
    S6PR("C = %M\n ",e.C);
    printf ("T11 = %f\n", T11);
    printf ("T12 = %f\n", T12);
    printf ("T13 = %f\n", T13);
    printf ("T16 = %f\n", T16);
    printf ("T22 = %f\n", T22);
    printf ("T44 = %f\n", T44);
    printf ("T55 = %f\n\n", T55);
    ElasticityT4ObserverToAbsolute(e.o,e.C,e.c);
    S6PR("c = %M\n ",e.c);
    ElasticityT4Relabel(e.C, 1, 2, e.C);
    /* 1 labeled the triple axis, convert to 2 */
    ElasticityT4Relabel(e.C, 1, 0, e.C);
    /* now: 1 labeled the reflection plane, convert to 0 */
    S6PR("C = %M\n ",e.C);
    /* this is actually Trigonal6 system */
    ElasticityTrigonal6T4Assign(T11,T13,T12,T16,T22,T44,e.C);
    S6PR("C = %M\n ",e.C);
    return (0);
}
#endif /* _ElasticityT4_TEST */


/* Return strain energy per unit volume given the */
/* elastic constant (c) and strain; one also can  */
/* pass in the elastic compliance (s) and stress. */
double ElasticityStrainEnergyDensity(T4 c, T2 strain)
{
    T2 stress;
    ElasticityT4MulT2(c,strain,stress);
    return (ElasticityT2Contraction(stress,strain)/2);
} /* end ElasticityStrainEnergyDensity() */


/* Young's modulus in "pull_axis" direction */
double ElasticityYoungsModulus (T4 c, int pull_axis)
{
    T4 s;
    ElasticityT4Inverse(c,s);
    return (1./s[pull_axis][pull_axis]);
} /* end ElasticityYoungsModulus() */

/* Autocomplete c, then return Young's modulus in "pull_axis" direction */
double ElasticityYOUNGSModulus (T4 c, int pull_axis)
{
    ElasticityT4Complete(c);
    return (ElasticityYoungsModulus(c,pull_axis));
} /* end ElasticityYOUNGSModulus() */


/* Negative ratio of the transverse to the longitudinal strain */
double ElasticityPoissonsRatio (T4 c, int pull_axis, int shrink_axis)
{
    T4 s;
    ElasticityT4Inverse(c,s);
    if (pull_axis == shrink_axis)
        pe ("ElasticityPoissonsRatio: pull_axis is shrink_axis.\n");
    return (-s[shrink_axis][pull_axis]/s[pull_axis][pull_axis]);
} /* end ElasticityPoissonsRatio() */

/* Autocomplete c, then return ElasticityPoissonsRatio() */
double ElasticityPOISSONSRatio (T4 c, int pull_axis, int shrink_axis)
{
    ElasticityT4Complete(c);
    return (ElasticityPoissonsRatio(c,pull_axis,shrink_axis));
} /* end ElasticityPOISSONSRatio() */

/* Average negative ratio of the transverse to the longitudinal strain */
double ElasticityPOISSONSRATIO (T4 c, int pull_axis)
{
    int i,j;
    T4 s;
    switch (pull_axis)
    {
        case 0:
            i = 1;
            j = 2;
            break;
        case 1:
            i = 0;
            j = 2;
            break;
        case 2:
            i = 0;
            j = 1;
            break;
        default:
            pe ("ElasticityPOISSONSRATIO: pull_axis = %d.\n", pull_axis);
    }
    ElasticityT4Complete(c);
    ElasticityT4Inverse(c,s);
    return (-0.5*(s[i][pull_axis]+s[j][pull_axis])/s[pull_axis][pull_axis]);
} /* end ElasticityPOISSONSRATIO() */


/* From any two of lambda,mu,E,nu,B, infer the rest */
void IsotropicStiffnessComplete (IsotropicStiffness *s)
{
    if (ElasticityValueKnown(s->lambda) &&
        ElasticityValueKnown(s->mu))
    { /* lambda,mu -> E,nu,B */
        s->E = s->mu * (3 * s->lambda + 2 * s->mu) / (s->lambda + s->mu);
        s->nu = s->lambda / 2 / (s->lambda + s->mu);
        s->B = s->lambda + 2 * s->mu / 3;
    }
    else if (ElasticityValueKnown(s->E) &&
             ElasticityValueKnown(s->nu))
    { /* E,nu -> lambda,mu,B */
        s->lambda = s->E * s->nu / (1 + s->nu) / (1 - 2 * s->nu);
        s->mu = s->E / 2 / (1 + s->nu);
        s->B = s->E / 3 / (1 - 2 * s->nu);
    }
    else if (ElasticityValueKnown(s->E) &&
             ElasticityValueKnown(s->B))
    { /* E,B -> lambda,mu,nu */
        s->lambda = 3 * s->B * (3 * s->B - s->E) / (9 * s->B - s->E);
        s->mu = 3 * s->B * s->E / (9 * s->B - s->E);
        s->nu = (3 * s->B - s->E) / 6 / s->B;
    }
    else if (ElasticityValueKnown(s->B) &&
             ElasticityValueKnown(s->nu))
    { /* B,nu -> lambda,mu,E */
        s->lambda = 3 * s->B * s->nu / (1 + s->nu);
        s->mu = 3 * s->B * (1 - 2 * s->nu) / 2 / (1 + s->nu);
        s->E = 3 * s->B * (1 - 2 * s->nu);
    }
    else if (ElasticityValueKnown(s->lambda) &&
             ElasticityValueKnown(s->nu))
    { /* lambda,nu -> mu,E,B */
        s->mu = s->lambda * (1 - 2 * s->nu) / 2 / s->nu;
        s->E = s->lambda * (1 + s->nu) * (1 - 2 * s->nu) / s->nu;
        s->B = s->lambda * (1 + s->nu) / 3 / s->nu;
    }
    else if (ElasticityValueKnown(s->mu) &&
             ElasticityValueKnown(s->nu))
    { /* mu,nu -> lambda,E,B */
        s->lambda = 2 * s->mu * s->nu / (1 - 2 * s->nu);
        s->E = 2 * s->mu * (1 + s->nu);
        s->B = 2 * s->mu * (1 + s->nu) / 3 / (1 - 2 * s->nu);
    }
    else if (ElasticityValueKnown(s->E) &&
             ElasticityValueKnown(s->mu))
    { /* E,mu -> lambda,nu,B */
        s->lambda = s->mu * (s->E - 2 * s->mu) / (3 * s->mu - s->E);
        s->nu = s->E / 2 / s->mu - 1;
        s->B = s->mu * s->E / 3 / (3 * s->mu - s->E);
    }
    else if (ElasticityValueKnown(s->B) &&
             ElasticityValueKnown(s->mu))
    { /* B,mu -> lambda,nu,E */
        s->lambda = s->B - 2 * s->mu / 3;
        s->nu = (3 * s->B - 2 * s->mu) / 2 / (3 * s->B + s->mu);
        s->E = 9 * s->B * s->mu / (3 * s->B + s->mu);
    }
    else if (ElasticityValueKnown(s->B) &&
             ElasticityValueKnown(s->lambda))
    { /* B,lambda -> mu,nu,E */
        s->mu = 3 * (s->B - s->lambda) / 2;
        s->nu = s->lambda / (3 * s->B - s->lambda);
        s->E = 9 * s->B * (s->B - s->lambda) / (3 * s->B - s->lambda);
    }
    else if (ElasticityValueKnown(s->E) &&
             ElasticityValueKnown(s->lambda))
    { /* E,lambda -> nu,mu,B */
        s->nu = (sqrt(SQUARE(s->E/s->lambda+1)+8) -
                 (s->E/s->lambda+1)) / 4;
        /* above is the +root: -root is feasible but unlikely */
        s->mu = s->lambda * (1 - 2 * s->nu) / 2 / s->nu;
        s->B = s->lambda * (1 + s->nu) / 3 / s->nu;
    }
    else pe ("IsotropicStiffnessComplete: Not enough \n"
             "is known to complete.\n");
    return;
} /* end IsotropicStiffnessComplete() */


#ifdef _IsotropicStiffnessComplete
int main (int argc, char *argv[])
{
    int i,j,known=0;
    IsotropicStiffness input=UnknownIsotropicStiffness,output;
    TermString s;
    printf ("\n============== Isotropic Elasticity ===========\n");
    printf ("From any 2 of lambda,mu,E,nu,B, infer the rest:\n");
    printf ("--------------------- Input -------------------\n");
    printf ("lambda = ");
    for (i=0; i<TERMCHAR; i++)
        if ((s[i]=getchar()) == '\n') break;
    s[i] = EOS;
    if (i>0)
    {
        input.lambda = atof(s);
        if ((++known)==2) goto finish;
    }
    printf ("    mu = ");
    for (i=0; i<TERMCHAR; i++)
        if ((s[i]=getchar()) == '\n') break;
    s[i] = EOS;
    if (i>0)
    {
        input.mu = atof(s);
        if ((++known)==2) goto finish;
    }
    printf ("     E = ");
    for (i=0; i<TERMCHAR; i++)
        if ((s[i]=getchar()) == '\n') break;
    s[i] = EOS;
    if (i>0)
    {
        input.E = atof(s);
        if ((++known)==2) goto finish;
    }
    printf ("    nu = ");
    for (i=0; i<TERMCHAR; i++)
        if ((s[i]=getchar()) == '\n') break;
    s[i] = EOS;
    if (i>0)
    {
        input.nu = atof(s);
        if ((++known)==2) goto finish;
    }
    printf ("     B = ");
    for (i=0; i<TERMCHAR; i++)
        if ((s[i]=getchar()) == '\n') break;
    s[i] = EOS;
    if (i>0)
    {
        input.B = atof(s);
        if ((++known)==2) goto finish;
    }
    pe ("Please input at least two.\n");
  finish:
    printf ("--------------------- Output ------------------\n");
    output = input;
    IsotropicStiffnessComplete (&output);
    if (ElasticityValueUnknown(input.lambda))
        printf ("lambda = %.16g\n", output.lambda);
    if (ElasticityValueUnknown(input.mu))
        printf ("    mu = %.16g\n", output.mu);
    if (ElasticityValueUnknown(input.E))
        printf ("     E = %.16g\n", output.E);
    if (ElasticityValueUnknown(input.nu))
        printf ("    nu = %.16g\n", output.nu);
    if (ElasticityValueUnknown(input.B))
        printf ("     B = %.16g\n", output.B);
    printf ("===============================================\n\n");
    return (0);
}
#endif /* _IsotropicStiffnessComplete */


/* Polycrystals under uniform strain approx. (Hirth & Lothe, pp. 425) */
void VoigtAveragedIsotropicStiffness(T4 c, IsotropicStiffness *stiffness)
{
    double i1, i2;
    i1 = c[0][0] + c[1][1] + c[2][2] + 2 * (c[3][3] + c[4][4] + c[5][5]);
    i2 = c[0][0] + c[1][1] + c[2][2] + 2 * (c[0][1] + c[0][2] + c[1][2]);
    stiffness->lambda = (2 * i2 - i1) / 15;
    stiffness->mu = (3 * i1 - i2) / 30;
    stiffness->E = stiffness->mu * (3 * stiffness->lambda + 2 * stiffness->mu)
        / (stiffness->lambda + stiffness->mu);
    stiffness->nu = stiffness->lambda / 2 /
        (stiffness->lambda + stiffness->mu);
    stiffness->B = stiffness->lambda + 2 * stiffness->mu / 3;
    return;
} /* end VoigtAveragedIsotropicStiffness() */

#ifdef _VoigtAveragedIsotropicStiffness_TEST
int main (int argc, char *argv[])
{
    T4 c;
    IsotropicStiffness stiffness;
    double lambda=10.2, mu=3.5;
    ElasticityIsotropicT4Assign(lambda,mu,c);
    VoigtAveragedIsotropicStiffness(c, &stiffness);
    printf ("%g %g\n", stiffness.lambda, stiffness.mu);
    return (0);
}
#endif /* _VoigtAveragedIsotropicStiffness_TEST */


/* Frame-invariant quantification of elastic anisotropy */
double ElasticityVoigtAnisotropy (T4 c)
{
    int i,j,k,l;
    double cc,dd;
    IsotropicStiffness stiffness[1];
    T4 cavg;
    VoigtAveragedIsotropicStiffness(c,stiffness);
    IsotropicStiffnessT4Assign(stiffness,cavg);
    cc = 0;
    dd = 0;
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            for (k=0; k<3; k++)
                for (l=0; l<3; l++)
                {
                    cc += SQUARE( T4Element(c,i,j,k,l) -
                                  T4Element(cavg,i,j,k,l) );
                    dd += SQUARE( T4Element(cavg,i,j,k,l) );
                }
    return(sqrt(cc/dd));
} /* end ElasticityVoigtAnisotropy() */


/* Ensure anisotropy by adding cubic symmetry perturbation: */
/* if anisotropy(cprev) > prev_tolerance, then caft=cprev;  */
/* else caft = VoigtAverage(cprev) + cubic(aft_amplitude).  */
void ElasticityEnsureAnisotropy
(T4 cprev, double prev_tolerance, double aft_amplitude, T4 caft)
{
    int i,j,k,l;
    double cc,dd;
    IsotropicStiffness stiffness[1];
    T4 cavg;
    VoigtAveragedIsotropicStiffness(cprev,stiffness);
    IsotropicStiffnessT4Assign(stiffness,cavg);
    cc = 0;
    dd = 0;
    for (i=0; i<3; i++)
        for (j=0; j<3; j++)
            for (k=0; k<3; k++)
                for (l=0; l<3; l++)
                {
                    cc += SQUARE( T4Element(cprev,i,j,k,l) -
                                  T4Element(cavg,i,j,k,l) );
                    dd += SQUARE( T4Element(cavg,i,j,k,l) );
                }
    if (sqrt(cc/dd) > prev_tolerance)
        ElasticityT4Eqv(cprev,caft);
    else
    {
        ElasticityT4Eqv(cavg,caft);
        dd = sqrt(dd/30) * aft_amplitude;
        caft[0][0] -= 2 * dd;
        caft[1][1] -= 2 * dd;
        caft[2][2] -= 2 * dd;
        caft[0][1] += dd;
        caft[0][2] += dd;
        caft[1][2] += dd;
        caft[1][0] += dd;
        caft[2][0] += dd;
        caft[2][1] += dd;
        caft[3][3] += dd;
        caft[4][4] += dd;
        caft[5][5] += dd;
    }
    return;
} /* end ElasticityEnsureAnisotropy() */

#ifdef _ElasticityEnsureAnisotropy_TEST
int main (int argc, char *argv[])
{
    T4 c,C;
    double lambda=10.2,mu=3.5,tolerance=0.01;
    M3 observer=GuineaPigOrthonormalRightHandedM3;
    ElasticityIsotropicT4Assign(lambda,mu,c);
    ElasticityENSUREAnisotropy(c,tolerance);
    S6PR("c = %M\n ", c);
    printf ("anisotropy = %e\n\n", ElasticityVoigtAnisotropy(c));
    ElasticityT4AbsoluteToObserver(observer,c,C);
    S6PR("C = %M\n ", C);
    printf ("anisotropy = %e\n", ElasticityVoigtAnisotropy(C));
    return (0);
}
#endif /* _ElasticityEnsureAnisotropy_TEST */


/* Polycrystals under uniform stress approx. (Hirth & Lothe, pp. 426) */
void ReussAveragedIsotropicStiffness(T4 c, IsotropicStiffness *stiffness)
{
    T4 s;
    double i1,i2;
    ElasticityT4Inverse(c,s);
    i1 = s[0][0] + s[1][1] + s[2][2] + 2 * (s[3][3] + s[4][4] + s[5][5]);
    i2 = s[0][0] + s[1][1] + s[2][2] + 2 * (s[0][1] + s[0][2] + s[1][2]);
    stiffness->E  = 15 / (2 * i1 + i2);
    stiffness->mu = 15 / (6 * i1 - 2 * i2);
    stiffness->lambda = stiffness->mu * (stiffness->E - 2 * stiffness->mu)
        / (3 * stiffness->mu - stiffness->E);
    stiffness->nu = stiffness->E / 2 / stiffness->mu - 1;
    stiffness->B = stiffness->mu * stiffness->E / 3
        / (3 * stiffness->mu - stiffness->E);
    return;
} /* end ReussAveragedIsotropicStiffness() */

#ifdef _ReussAveragedIsotropicStiffness_TEST
int main (int argc, char *argv[])
{
    T4 c;
    IsotropicStiffness stiffness;
    double lambda=10.2, mu=3.5;
    ElasticityIsotropicT4Assign(lambda,mu,c);
    ReussAveragedIsotropicStiffness(c, &stiffness);
    printf ("%g %g\n", stiffness.lambda, stiffness.mu);
    return (0);
}
#endif /* _ReussAveragedIsotropicStiffness_TEST */
