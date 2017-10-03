#ifndef _Elasticity_h
#define _Elasticity_h

#include <IO.h>
#include <VecMat.h>
#include <VecMat2.h>
#include <VecMat3.h>

typedef double M6[6][6];
typedef double (*M6P)[6];
typedef M3 T2; /* rank-2 symmetric tensor */
typedef M6 T4; /* rank-4 symmetric tensor */

/* memory expansion */
#define V6e(a) (a),(a)+1,(a)+2,(a)+3,(a)+4,(a)+5
#define V6E(a) (a)[0],(a)[1],(a)[2],(a)[3],(a)[4],(a)[5]
#define M6e(A) V6e(A[0]),V6e(A[1]),V6e(A[2]),V6e(A[3]),V6e(A[4]),V6e(A[5])
#define M6E(A) V6E(A[0]),V6E(A[1]),V6E(A[2]),V6E(A[3]),V6E(A[4]),V6E(A[5])
#define V9e(a) (a),(a)+1,(a)+2,(a)+3,(a)+4,(a)+5,(a)+6,(a)+7,(a)+8
#define V9E(a) (a)[0],(a)[1],(a)[2],(a)[3],(a)[4],(a)[5],(a)[6],(a)[7],(a)[8]

/* Elasticity.c: */

extern const int ElasticityVoigtIndices[3][3];
#define VoigtIndex(i,j) ElasticityVoigtIndices[i][j]
#define VoigtINDEX(i,j) (ElasticityVoigtIndices[(i)-1][(j)-1]+1)
#define T4Element(A,i,j,k,l) ((A)[VoigtIndex(i,j)][VoigtIndex(k,l)])
#define T4ELEMENT(A,i,j,k,l) \
  ((A)[VoigtIndex((i)-1,(j)-1)][VoigtIndex((k)-1,(l)-1)])

/*************************************************************/
/* All variables evaluated in absolute frame are lower-case  */
/* All variables evaluated in observer frame are capitalized */
/* --------------------------------------------------------- */
/* To understand the notation, take fcc lattice as example:  */
/* Absolute frame: the cubic vantage point                   */
/* c[][] would correspond to c11,c12,c44                     */
/* --------------------------------------------------------- */
/* Observer frame: 0=[1 1 -2], 1=[1 1 1], 2=[1 -1 0]         */
/* H[][] of MD is often in view of the observer frame        */
/* C[][] would correspond to resolved shear modulus etc      */
/*************************************************************/
typedef struct
{
    M3 o;      /* local observer orthogonal frame in absolute frame */
    M3 h0;     /* undeformed RVE edges in absolute frame */
    M3 H0;     /* undeformed RVE edges in observer frame */
    M3 h;      /* current RVE edges in absolute frame */
    M3 H;      /* current RVE edges in observer frame */
    T2 strain; /* strain tensor in absolute frame */
    T2 Strain; /* strain tensor in observer frame */
    T2 stress; /* stress tensor in absolute frame */
    T2 Stress; /* stress tensor in observer frame */
    T4 c;      /* elastic constant in absolute frame */
    T4 C;      /* elastic constant in observer frame */
    T4 s;      /* elastic compliance in absolute frame */
    T4 S;      /* elastic compliance in observer frame */
    T4 m;      /* "relaxed" elastic modulus in absolute frame */
    T4 M;      /* "relaxed" elastic modulus in observer frame */
} Elasticity;

#define UnknownElasticityValue  DOUBLE_PRECISION_INFINITY
#define UnknownElasticityV3 \
  {UnknownElasticityValue,UnknownElasticityValue,UnknownElasticityValue}
#define UnknownElasticityM3 \
  {UnknownElasticityV3,UnknownElasticityV3,UnknownElasticityV3}
#define UnknownElasticityV6 \
  {UnknownElasticityValue,UnknownElasticityValue,UnknownElasticityValue, \
  UnknownElasticityValue,UnknownElasticityValue,UnknownElasticityValue}
#define UnknownElasticityM6 \
  {UnknownElasticityV6,UnknownElasticityV6,UnknownElasticityV6, \
  UnknownElasticityV6,UnknownElasticityV6,UnknownElasticityV6}
#define UnknownElasticity \
  {UnknownElasticityM3,UnknownElasticityM3,UnknownElasticityM3, \
  UnknownElasticityM3,UnknownElasticityM3,UnknownElasticityM3, \
  UnknownElasticityM3,UnknownElasticityM3,UnknownElasticityM3, \
  UnknownElasticityM6,UnknownElasticityM6,UnknownElasticityM6, \
  UnknownElasticityM6}

#define ElasticityValueUnknown(x) ((x)==UnknownElasticityValue)
#define ElasticityV3Unknown(V)    ( ElasticityValueUnknown((V)[0]) || \
  ElasticityValueUnknown((V)[1]) || ElasticityValueUnknown((V)[2]) )
#define ElasticityM3Unknown(M) ( ElasticityV3Unknown((M)[0]) || \
  ElasticityV3Unknown((M)[1]) || ElasticityV3Unknown((M)[2]) )
#define ElasticityV6Unknown(V)    ( ElasticityValueUnknown((V)[0]) || \
  ElasticityValueUnknown((V)[1]) || ElasticityValueUnknown((V)[2]) || \
  ElasticityValueUnknown((V)[3]) || ElasticityValueUnknown((V)[4]) || \
  ElasticityValueUnknown((V)[5]) )
#define ElasticityM6Unknown(M) ( ElasticityV6Unknown((M)[0]) || \
  ElasticityV6Unknown((M)[1]) || ElasticityV6Unknown((M)[2]) || \
  ElasticityV6Unknown((M)[3]) || ElasticityV6Unknown((M)[4]) || \
  ElasticityV6Unknown((M)[5]) )

#define ElasticityValueKnown(x) ((x)!=UnknownElasticityValue)
#define ElasticityV3Known(V)    ( ElasticityValueKnown((V)[0]) && \
  ElasticityValueKnown((V)[1]) && ElasticityValueKnown((V)[2]) )
#define ElasticityM3Known(M) ( ElasticityV3Known((M)[0]) && \
  ElasticityV3Known((M)[1]) && ElasticityV3Known((M)[2]) )
#define ElasticityV6Known(V)    ( ElasticityValueKnown((V)[0]) && \
  ElasticityValueKnown((V)[1]) && ElasticityValueKnown((V)[2]) && \
  ElasticityValueKnown((V)[3]) && ElasticityValueKnown((V)[4]) && \
  ElasticityValueKnown((V)[5]) )
#define ElasticityM6Known(M) ( ElasticityV6Known((M)[0]) && \
  ElasticityV6Known((M)[1]) && ElasticityV6Known((M)[2]) && \
  ElasticityV6Known((M)[3]) && ElasticityV6Known((M)[4]) && \
  ElasticityV6Known((M)[5]) )

#define ElasticityV3Written(V)    ( ElasticityValueKnown((V)[0]) || \
  ElasticityValueKnown((V)[1]) || ElasticityValueKnown((V)[2]) )
#define ElasticityM3Written(M) ( ElasticityV3Written((M)[0]) || \
  ElasticityV3Written((M)[1]) || ElasticityV3Written((M)[2]) )
#define ElasticityV6Written(V)    ( ElasticityValueKnown((V)[0]) || \
  ElasticityValueKnown((V)[1]) || ElasticityValueKnown((V)[2]) || \
  ElasticityValueKnown((V)[3]) || ElasticityValueKnown((V)[4]) || \
  ElasticityValueKnown((V)[5]) )
#define ElasticityM6Written(M) ( ElasticityV6Written((M)[0]) || \
  ElasticityV6Written((M)[1]) || ElasticityV6Written((M)[2]) || \
  ElasticityV6Written((M)[3]) || ElasticityV6Written((M)[4]) || \
  ElasticityV6Written((M)[5]) )

/* Autocomplete orthonormal observer directions matrix */
void ElasticityObserverComplete (M3 observer);

/* Convert a vector in Observer notation to Absolute notation */
#define V3ObserverToAbsolute(observer,observer_notation,absolute_notation) \
  V3mM3(observer_notation,observer,absolute_notation)
/* Convert a vector in Absolute notation to Observe notation */
#define V3AbsoluteToObserver(observer,absolute_notation,observer_notation) \
  M3mV3(observer,absolute_notation,observer_notation)

/* Autocomplete a rank-2 symmetric tensor */
void ElasticityT2Complete (T2 tensor);

/* Hydrostatic invariant of the tensor */
#define ElasticityT2HydroInvariant(tensor) (M3TR(tensor)/3.)

/* von Mises invariant of the tensor */
#define ElasticityT2MisesInvariant(t) \
  sqrt( SQUARE((t)[0][1]) + SQUARE((t)[0][2]) + SQUARE((t)[1][2]) + \
  (SQUARE((t)[0][0]-(t)[1][1]) + SQUARE((t)[0][0]-(t)[2][2]) + \
  SQUARE((t)[1][1]-(t)[2][2]))/6. )

/* Convert a rank-2 tensor in Observer notation to Absolute notation */
void ElasticityT2ObserverToAbsolute
(M3 observer, T2 observer_notation, T2 absolute_notation);

/* Convert a rank-2 tensor in Absolute notation to Observer notation */
void ElasticityT2AbsoluteToObserver
(M3 observer, T2 absolute_notation, T2 observer_notation);

/* \sum_{i=0}^{2}\sum_{j=0}^{2} A_{ij}*B_{ij}. A,B both are symmetric */
#define ElasticityT2Contraction(A,B) ( \
  A[0][0]*B[0][0] + A[1][1]*B[1][1] + A[2][2]*B[2][2] + 2 * ( \
  A[0][1]*B[0][1] + A[0][2]*B[0][2] + A[1][2]*B[1][2] ) )

/* return ElasticityT2Contraction((A+B)/2,C) */
double ElasticityT2ContractionTrapezoidal(M3 A, M3 B, M3 C);
#define ElasticityT2ContractionTRAPEZOIDAL(A,B,C,TMP) ( \
  M3ADD(A,B,TMP), ElasticityT2Contraction(TMP,C)/2. )

/* Autocomplete a rank-4 symmetric tensor */
void ElasticityT4Complete (T4 tensor);

#define ElasticityT4Zero(tensor) bzero((void *)(tensor[0]),36*sizeof(double))
/* b[][] = a[][] */
#define ElasticityT4Eqv(a,b) \
  memcpy((void *)(b[0]),(void *)(a[0]),36*sizeof(double))

/* Cubic: 3 independents */
#define ElasticityCubicT4Assign(T11,T12,T44,tensor) \
  ( ElasticityT4Zero(tensor), tensor[0][0]=tensor[1][1]=tensor[2][2]=T11, \
  tensor[0][1]=tensor[1][0]=tensor[0][2]=tensor[2][0]= \
  tensor[1][2]=tensor[2][1]=T12, tensor[3][3]=tensor[4][4]=tensor[5][5]=T44 )

/* Isotropic: 2 independents */
/* T_{ijkl} = Lambda * delta_{ij} * delta_{kl} + Mu *  */
/* (delta_{ik} * delta_{jl} + delta_{il} * delta_{jk}) */
#define ElasticityIsotropicT4Assign(Lambda,Mu,tensor) \
  ElasticityCubicT4Assign((Lambda)+2*(Mu),Lambda,Mu,tensor)

/* 1-2 basal plane isotropic, 3 sticking out */
#define ElasticityPlanarIsotropyT4Assign(T11,T12,T13,T33,T44,tensor) \
  ( ElasticityT4Zero(tensor), tensor[0][0]=tensor[1][1]=T11, \
  tensor[0][1]=tensor[1][0]=T12, tensor[2][2]=T33, \
  tensor[0][2]=tensor[2][0]=tensor[1][2]=tensor[2][1]=T13, \
  tensor[3][3]=tensor[4][4]=T44, tensor[5][5]=((T11)-(T12))/2 )

/* Trigonal (3 is triple axis): */
/* Trigonal7 (3,\bar{3}): 7 independents (T14=-T24=T56,T15=-T25=-T46) */
#define ElasticityTrigonal7T4Assign(T11,T12,T13,T14,T15,T33,T44,tensor) \
  ( tensor[0][0]=T11, \
  tensor[0][1]=tensor[1][0]=T12, \
  tensor[0][2]=tensor[2][0]=T13, \
  tensor[0][3]=tensor[3][0]=T14, \
  tensor[0][4]=tensor[4][0]=T15, \
  tensor[0][5]=tensor[5][0]=0, \
  tensor[1][1]=T11, \
  tensor[1][2]=tensor[2][1]=T13, \
  tensor[1][3]=tensor[3][1]=-(T14), \
  tensor[1][4]=tensor[4][1]=-(T15), \
  tensor[1][5]=tensor[5][1]=0, \
  tensor[2][2]=T33, \
  tensor[2][3]=tensor[3][2]=0, \
  tensor[2][4]=tensor[4][2]=0, \
  tensor[2][5]=tensor[5][2]=0, \
  tensor[3][3]=T44, \
  tensor[3][4]=tensor[4][3]=0, \
  tensor[3][5]=tensor[5][3]=-(T15), \
  tensor[4][4]=T44, \
  tensor[4][5]=tensor[5][4]=T14, \
  tensor[5][5]=((T11)-(T12))/2 )
/* Rhombohedral is just another name for Trigonal */
#define ElasticityRhombohedral7T4Assign(T11,T12,T13,T14,T15,T33,T44,tensor) \
            ElasticityTrigonal7T4Assign(T11,T12,T13,T14,T15,T33,T44,tensor)

/* If additionally, 1 or 2 is mirror plane or two-fold axes: */
/* Trigonal6 (32,\bar{3}m,3m): 6 independents (T15=-T25=-T46=0) */
#define ElasticityTrigonal6T4Assign(T11,T12,T13,T14,T33,T44,tensor) \
        ElasticityTrigonal7T4Assign(T11,T12,T13,T14,0,T33,T44,tensor)
#define ElasticityRhombohedral6T4Assign(T11,T12,T13,T14,T33,T44,tensor) \
            ElasticityTrigonal6T4Assign(T11,T12,T13,T14,T33,T44,tensor)

/* Tetragonal is orthorhombic with equivalent 1,2 axes (3 is 4-fold axis): */
/* Tetragonal7 (4,\bar{4},4/m): 7 independents (T26=-T16) */
#define ElasticityTetragonal7T4Assign(T11,T12,T13,T16,T33,T44,T66,tensor) \
  ( tensor[0][0]=T11, \
  tensor[0][1]=tensor[1][0]=T12, \
  tensor[0][2]=tensor[2][0]=T13, \
  tensor[0][3]=tensor[3][0]=0, \
  tensor[0][4]=tensor[4][0]=0, \
  tensor[0][5]=tensor[5][0]=T16, \
  tensor[1][1]=T11, \
  tensor[1][2]=tensor[2][1]=T13, \
  tensor[1][3]=tensor[3][1]=0, \
  tensor[1][4]=tensor[4][1]=0, \
  tensor[1][5]=tensor[5][1]=-(T16), \
  tensor[2][2]=T33, \
  tensor[2][3]=tensor[3][2]=0, \
  tensor[2][4]=tensor[4][2]=0, \
  tensor[2][5]=tensor[5][2]=0, \
  tensor[3][3]=T44, \
  tensor[3][4]=tensor[4][3]=0, \
  tensor[3][5]=tensor[5][3]=0, \
  tensor[4][4]=T44, \
  tensor[4][5]=tensor[5][4]=0, \
  tensor[5][5]=T66 )

/* If additionally, 1,2 are mirror planes or two-fold axes: */
/* Tetragonal6 (4mm,\bar{4}2m,422,4/mmm): 6 independents (T26=-T16=0) */
#define ElasticityTetragonal6T4Assign(T11,T12,T13,T33,T44,T66,tensor) \
        ElasticityTetragonal7T4Assign(T11,T12,T13,0,T33,T44,T66,tensor)

/* Orthorhombic (3 orthogonal planes): 9 independents */
#define ElasticityOrthorhombicT4Assign(T11,T12,T13,T22,T23,T33,T44,T55,T66,tensor) \
  ( tensor[0][0]=T11, \
  tensor[0][1]=tensor[1][0]=T12, \
  tensor[0][2]=tensor[2][0]=T13, \
  tensor[0][3]=tensor[3][0]=0, \
  tensor[0][4]=tensor[4][0]=0, \
  tensor[0][5]=tensor[5][0]=0, \
  tensor[1][1]=T22, \
  tensor[1][2]=tensor[2][1]=T23, \
  tensor[1][3]=tensor[3][1]=0, \
  tensor[1][4]=tensor[4][1]=0, \
  tensor[1][5]=tensor[5][1]=0, \
  tensor[2][2]=T33, \
  tensor[2][3]=tensor[3][2]=0, \
  tensor[2][4]=tensor[4][2]=0, \
  tensor[2][5]=tensor[5][2]=0, \
  tensor[3][3]=T44, \
  tensor[3][4]=tensor[4][3]=0, \
  tensor[3][5]=tensor[5][3]=0, \
  tensor[4][4]=T55, \
  tensor[4][5]=tensor[5][4]=0, \
  tensor[5][5]=T66 )

/* Monoclinic (3 is a two-fold axis): 13 independents */
#define ElasticityMonoclinicT4Assign(T11,T12,T13,T16,T22,T23,T26,T33,T36,T44,T45,T55,T66,tensor) \
  ( tensor[0][0]=T11, \
  tensor[0][1]=tensor[1][0]=T12, \
  tensor[0][2]=tensor[2][0]=T13, \
  tensor[0][3]=tensor[3][0]=0, \
  tensor[0][4]=tensor[4][0]=0, \
  tensor[0][5]=tensor[5][0]=T16, \
  tensor[1][1]=T22, \
  tensor[1][2]=tensor[2][1]=T23, \
  tensor[1][3]=tensor[3][1]=0, \
  tensor[1][4]=tensor[4][1]=0, \
  tensor[1][5]=tensor[5][1]=T26, \
  tensor[2][2]=T33, \
  tensor[2][3]=tensor[3][2]=0, \
  tensor[2][4]=tensor[4][2]=0, \
  tensor[2][5]=tensor[5][2]=T36, \
  tensor[3][3]=T44, \
  tensor[3][4]=tensor[4][3]=T45, \
  tensor[3][5]=tensor[5][3]=0, \
  tensor[4][4]=T55, \
  tensor[4][5]=tensor[5][4]=0, \
  tensor[5][5]=T66 )

/* Triclinic (no symmetry at all): 21 independents */
#define ElasticityTriclinicT4Assign(T11,T12,T13,T14,T15,T16,T22,T23,T24,T25,T26,T33,T34,T35,T36,T44,T45,T46,T55,T56,T66,tensor) \
  ( tensor[0][0]=T11, \
  tensor[0][1]=tensor[1][0]=T12, \
  tensor[0][2]=tensor[2][0]=T13, \
  tensor[0][3]=tensor[3][0]=T14, \
  tensor[0][4]=tensor[4][0]=T15, \
  tensor[0][5]=tensor[5][0]=T16, \
  tensor[1][1]=T22, \
  tensor[1][2]=tensor[2][1]=T23, \
  tensor[1][3]=tensor[3][1]=T24, \
  tensor[1][4]=tensor[4][1]=T25, \
  tensor[1][5]=tensor[5][1]=T26, \
  tensor[2][2]=T33, \
  tensor[2][3]=tensor[3][2]=T34, \
  tensor[2][4]=tensor[4][2]=T35, \
  tensor[2][5]=tensor[5][2]=T36, \
  tensor[3][3]=T44, \
  tensor[3][4]=tensor[4][3]=T45, \
  tensor[3][5]=tensor[5][3]=T46, \
  tensor[4][4]=T55, \
  tensor[4][5]=tensor[5][4]=T56, \
  tensor[5][5]=T66 )

/* c_{ij} = \sum_{i=0}^{2}\sum_{j=0}^{2} A_{ijkl} * b_{kl} */
#define ElasticityT4MulT2(A,b,c) \
  ( c[0][0] = A[0][0]*b[0][0] + A[0][1]*b[1][1] + A[0][2]*b[2][2] + \
  2 * (A[0][3]*b[1][2] + A[0][4]*b[0][2] + A[0][5]*b[0][1]), \
  c[1][1] = A[1][0]*b[0][0] + A[1][1]*b[1][1] + A[1][2]*b[2][2] + \
  2 * (A[1][3]*b[1][2] + A[1][4]*b[0][2] + A[1][5]*b[0][1]), \
  c[2][2] = A[2][0]*b[0][0] + A[2][1]*b[1][1] + A[2][2]*b[2][2] + \
  2 * (A[2][3]*b[1][2] + A[2][4]*b[0][2] + A[2][5]*b[0][1]), \
  c[1][2] = c[2][1] = A[3][0]*b[0][0] + A[3][1]*b[1][1] + A[3][2]*b[2][2] + \
  2 * (A[3][3]*b[1][2] + A[3][4]*b[0][2] + A[3][5]*b[0][1]), \
  c[0][2] = c[2][0] = A[4][0]*b[0][0] + A[4][1]*b[1][1] + A[4][2]*b[2][2] + \
  2 * (A[4][3]*b[1][2] + A[4][4]*b[0][2] + A[4][5]*b[0][1]), \
  c[0][1] = c[1][0] = A[5][0]*b[0][0] + A[5][1]*b[1][1] + A[5][2]*b[2][2] + \
  2 * (A[5][3]*b[1][2] + A[5][4]*b[0][2] + A[5][5]*b[0][1]) )

/* T_{ijkl} * TI_{klmn} = delta_{im} * delta_{jn}, but both */
/* T and TI are expressed in Voigt notation due to symmetry */
void ElasticityT4Inverse (T4 tensor, T4 inverse);

/* Element-by-element reciprocal */
void ElasticityT4ElementalReciprocal (T4 tensor, T4 reciprocal);

/* Transform a direct-6x6-inverse, which is a non-tensor, to a tensor */
void ElasticityT4tildeToT4 (T4 tilde, T4 tensor);
#define ElasticityT4tildeTOT4(tensor) ElasticityT4tildeToT4(tensor,tensor)

#define ElasticityT4RelaxedModuli(S,M) ( \
  ElasticityT4ElementalReciprocal(S,M), ElasticityT4tildeTOT4(M) )

/* Convert a rank-4 tensor in Observer notation to Absolute notation */
void ElasticityT4ObserverToAbsolute
(M3 observer, T4 observer_notation, T4 absolute_notation);

/* Convert a rank-4 tensor in Absolute notation to Observer notation */
void ElasticityT4AbsoluteToObserver
(M3 observer, T4 absolute_notation, T4 observer_notation);

/* What happens to a tensor, previously denoted "prev_notation" */
/* while an axis is labeled "prev_label", has the label changed */
/* to "current_label". The unrelated axis keeps its label.      */
void ElasticityT4Relabel
(T4 prev_notation, int prev_label, int current_label, T4 current_notation);

/* Return strain energy per unit volume given the */
/* elastic constant (c) and strain; one also can  */
/* pass in the elastic compliance (s) and stress. */
double ElasticityStrainEnergyDensity(T4 c, T2 strain);

/* B = -dP / dlnV */
#define ElasticityBulkModulus(c) ((c[0][0] + c[1][1] + c[2][2] + \
  2 * (c[0][1] + c[0][2] + c[1][2])) / 9)

/* Young's modulus in "pull_axis" direction */
double ElasticityYoungsModulus (T4 c, int pull_axis);
/* Autocomplete c, then return Young's modulus in "pull_axis" direction */
double ElasticityYOUNGSModulus (T4 c, int pull_axis);

/* Negative ratio of the transverse to the longitudinal strain */
double ElasticityPoissonsRatio (T4 c, int pull_axis, int shrink_axis);
/* Autocomplete c, then return ElasticityPoissonsRatio() */
double ElasticityPOISSONSRatio (T4 c, int pull_axis, int shrink_axis);
/* Average negative ratio of the transverse to the longitudinal strain */
double ElasticityPOISSONSRATIO (T4 c, int pull_axis);

typedef struct
{ /* 5 elastic properties of isotropic material */
    double lambda;  /* c12 */
    double mu;      /* shear modulus (c44) */
    double E;       /* Young's modulus */
    double nu;      /* Poisson's ratio */
    double B;       /* bulk modulus */
} IsotropicStiffness;

#define UnknownIsotropicStiffness \
  { UnknownElasticityValue, UnknownElasticityValue, UnknownElasticityValue, \
  UnknownElasticityValue, UnknownElasticityValue }

/* From any two of lambda,mu,E,nu,B, infer the rest */
void IsotropicStiffnessComplete (IsotropicStiffness *s);

/* Convert IsotropicStiffness with filled lambda,mu to T4 */
#define IsotropicStiffnessT4Assign(stiffness,c) \
  ElasticityIsotropicT4Assign((stiffness)->lambda,(stiffness)->mu,c)
/* Convert any valid IsotropicStiffness to T4 */
#define IsotropicStiffnessT4ASSIGN(stiffness,c) \
  (IsotropicStiffnessComplete(stiffness), \
  IsotropicStiffnessT4Assign(stiffness,c))

/* Polycrystals under uniform strain approx. (Hirth & Lothe, pp. 425) */
void VoigtAveragedIsotropicStiffness(T4 c, IsotropicStiffness *stiffness);

/* Frame-invariant quantification of elastic anisotropy */
double ElasticityVoigtAnisotropy (T4 c);

/* Ensure anisotropy by adding cubic symmetry perturbation: */
/* if anisotropy(cprev) > prev_tolerance, then caft=cprev;  */
/* else caft = VoigtAverage(cprev) + cubic(aft_amplitude).  */
void ElasticityEnsureAnisotropy
(T4 cprev, double prev_tolerance, double aft_amplitude, T4 caft);
#define ElasticityENSUREAnisotropy(c,tolerance) \
  ElasticityEnsureAnisotropy(c,tolerance,tolerance,c);

/* Polycrystals under uniform stress approx. (Hirth & Lothe, pp. 426) */
void ReussAveragedIsotropicStiffness(T4 c, IsotropicStiffness *stiffness);
/* Autocomplete c, then return ReussAveragedIsotropicStiffness() */
#define ReussAVERAGEDIsotropicStiffness(c,stiffness) \
  (ElasticityT4Complete(c), ReussAveragedIsotropicStiffness(c,stiffness))


/* Stroh.c: */

typedef double V6[6];

/* Hirth & Lothe: pp. 467: Eqn. 13-162 */
#define ElasticityStrohContraction(a,b,C,ab) ( \
  ab[0][0] = \
  (a)[0] * ( (b)[0]*C[0][0] + (b)[1]*C[0][5] + (b)[2]*C[0][4]) + \
  (a)[1] * ( (b)[0]*C[5][0] + (b)[1]*C[5][5] + (b)[2]*C[5][4]) + \
  (a)[2] * ( (b)[0]*C[4][0] + (b)[1]*C[4][5] + (b)[2]*C[4][4]), \
  ab[0][1] = \
  (a)[0] * ( (b)[0]*C[0][5] + (b)[1]*C[0][1] + (b)[2]*C[0][3]) + \
  (a)[1] * ( (b)[0]*C[5][5] + (b)[1]*C[5][1] + (b)[2]*C[5][3]) + \
  (a)[2] * ( (b)[0]*C[4][5] + (b)[1]*C[4][1] + (b)[2]*C[4][3]), \
  ab[0][2] = \
  (a)[0] * ( (b)[0]*C[0][4] + (b)[1]*C[0][3] + (b)[2]*C[0][2]) + \
  (a)[1] * ( (b)[0]*C[5][4] + (b)[1]*C[5][3] + (b)[2]*C[5][2]) + \
  (a)[2] * ( (b)[0]*C[4][4] + (b)[1]*C[4][3] + (b)[2]*C[4][2]), \
  ab[1][0] = \
  (a)[0] * ( (b)[0]*C[5][0] + (b)[1]*C[5][5] + (b)[2]*C[5][4]) + \
  (a)[1] * ( (b)[0]*C[1][0] + (b)[1]*C[1][5] + (b)[2]*C[1][4]) + \
  (a)[2] * ( (b)[0]*C[3][0] + (b)[1]*C[3][5] + (b)[2]*C[3][4]), \
  ab[1][1] = \
  (a)[0] * ( (b)[0]*C[5][5] + (b)[1]*C[5][1] + (b)[2]*C[5][3]) + \
  (a)[1] * ( (b)[0]*C[1][5] + (b)[1]*C[1][1] + (b)[2]*C[1][3]) + \
  (a)[2] * ( (b)[0]*C[3][5] + (b)[1]*C[3][1] + (b)[2]*C[3][3]), \
  ab[1][2] = \
  (a)[0] * ( (b)[0]*C[5][4] + (b)[1]*C[5][3] + (b)[2]*C[5][2]) + \
  (a)[1] * ( (b)[0]*C[1][4] + (b)[1]*C[1][3] + (b)[2]*C[1][2]) + \
  (a)[2] * ( (b)[0]*C[3][4] + (b)[1]*C[3][3] + (b)[2]*C[3][2]), \
  ab[2][0] = \
  (a)[0] * ( (b)[0]*C[4][0] + (b)[1]*C[4][5] + (b)[2]*C[4][4]) + \
  (a)[1] * ( (b)[0]*C[3][0] + (b)[1]*C[3][5] + (b)[2]*C[3][4]) + \
  (a)[2] * ( (b)[0]*C[2][0] + (b)[1]*C[2][5] + (b)[2]*C[2][4]), \
  ab[2][1] = \
  (a)[0] * ( (b)[0]*C[4][5] + (b)[1]*C[4][1] + (b)[2]*C[4][3]) + \
  (a)[1] * ( (b)[0]*C[3][5] + (b)[1]*C[3][1] + (b)[2]*C[3][3]) + \
  (a)[2] * ( (b)[0]*C[2][5] + (b)[1]*C[2][1] + (b)[2]*C[2][3]), \
  ab[2][2] = \
  (a)[0] * ( (b)[0]*C[4][4] + (b)[1]*C[4][3] + (b)[2]*C[4][2]) + \
  (a)[1] * ( (b)[0]*C[3][4] + (b)[1]*C[3][3] + (b)[2]*C[3][2]) + \
  (a)[2] * ( (b)[0]*C[2][4] + (b)[1]*C[2][3] + (b)[2]*C[2][2]) )

/* Subtract momentum derivatives from ElasticityStrohContraction() for */
/* uniformly moving system: D.J. Bacon, D.M. Barnett, R.O. Scattergood */
/* Progress in Materials Science 23 (1978) 53-262. pp. 136.            */
void ElasticityStrohCONTRACTION
(V3 a, V3 b, T4 C, double mass_density, V3 velocity, M3 ab);

/* C[][] := A[][] * B[][] */
void ElasticityStrohM6mul (M6 A, M6 B, M6 C);

/* Diagonalize 6x6 real nonsymmetric matrix */
void ElasticityStrohM6Diag (M6 A, V6 wr, V6 wi, M6 VR, M6 VI);


/* StraightDislocation.c: */

#define MODES 3
/* total 6 modes, but we use just the 3 modes with positive imag(eigenvalue) */
typedef struct
{
    /** Input: **/
    V3 xi;  /* dislocation line direction in absolute frame */
    V3 m;   /* cut direction in absolute frame */
    V3 b;   /* Burgers vector in absolute frame */
    double mass_density;  /* mass density of the crystal */
    V3 velocity;  /* velocity of the dislocation in absolute frame */
    T4 c;   /* crystal elastic constant in absolute frame */

    /** Results: **/
    V3 n;   /* m x n = xi: Hirth & Lothe: pp. 468: Fig. 13-12 */
    M3 observer;  /* (m,n,xi) frame */

    V3 B;  /* Burgers vector in observer frame */
    V3 Velocity;  /* translation velocity of the dislocation in o.f. */
    T4 C;   /* crystal elastic constant in observer frame */

    double pr[MODES];  /* real( modal eigenvalue ) */
    double pi[MODES];  /* imag( modal eigenvalue ) */

    double ar[MODES][3];  /* real( displacement eigenvector ) */
    double ai[MODES][3];  /* imag( displacement eigenvector ) */

    double lr[MODES][3];  /* real( stress function eigenvector ) */
    double li[MODES][3];  /* imag( stress function eigenvector ) */

    double dr[MODES];     /* real( combination coefficient ) */
    double di[MODES];     /* imag( combination coefficient ) */

    double ur[MODES][3];  /* real( displacement field ) */
    double ui[MODES][3];  /* imag( displacement field ) */
    double UR[MODES][3];  /* real( observer frame displacement field ) */
    double UI[MODES][3];  /* imag( observer frame displacement field ) */

    double uxr[MODES][3][3]; /* real( displacement gradient field ) */
    double uxi[MODES][3][3]; /* imag( displacement gradient field ) */
    double UXR[MODES][3][3]; /* real( o.f. displacement gradient field ) */
    double UXI[MODES][3][3]; /* imag( o.f. displacement gradient field ) */

    double udotr[MODES][3];  /* real( velocity field ) */
    double udoti[MODES][3];  /* imag( velocity field ) */
    double UdotR[MODES][3];  /* real( o.f. velocity field ) */
    double UdotI[MODES][3];  /* imag( o.f. velocity field ) */

    double strainr[MODES][3][3]; /* real( strain field ) */
    double straini[MODES][3][3]; /* imag( strain field ) */
    double StrainR[MODES][3][3]; /* real( o.f. strain field ) */
    double StrainI[MODES][3][3]; /* imag( o.f. strain field ) */

    double stressr[MODES][3][3]; /* real( stress field ) */
    double stressi[MODES][3][3]; /* imag( stress field ) */
    double StressR[MODES][3][3]; /* real( o.f. stress field ) */
    double StressI[MODES][3][3]; /* imag( o.f. stress field ) */

    /* "if-static" energy coefficient tensor (H&L p.471, Eqn. 13-189): */
    double kmode[MODES][3][3];
    M3 k;
    /* see also, A.N. Stroh, J. Math. Phys. 41 (1962) 77-103: p.91-92. */
    M3 K; /* observer frame "if-static" energy coefficient tensor */

    /***********************************************************/
    /* If we create a static dipole of distance R by making a  */
    /* cut along m, the leading behavior in R of the work done */
    /* is 2*PE*log(R) as H&L p.471 13-187, w/ PE=b'*k*b/4/PI.  */
    /* There, the dislocation(s) is presumed to be NOT MOVING. */
    /* Also, if we choose another m', creating a different     */
    /* dipole along a different cut, this leading behavior     */
    /* remains the same, so PE is invariant to m (H&L pp. 474, */
    /* Axiom 13-1). But for the next-order term which's zeroth */
    /* order of R, things can be tricky. This term may depend  */
    /* on m (angular force between dislocations) as H&L p.117, */
    /* 5-16 indicates just due to the asymptotic double loop   */
    /* setup, and in reality one may also observe this effect  */
    /* due to a variable integration cutoff radius r0(m).      */
    /*                                                         */
    /* For a single, MOVING dislocation one can study the R-   */
    /* dependence of the potential and kinetic energies INSIDE */
    /* a circle of radius R, see A.N. Stroh, J. Math. Phys. 41 */
    /* (1962) 77-103: p.91, Eqns. 7.3 and 7.4. One finds that  */
    /* kine=KE*log(R), pote=PE*log(R), and PE=KE+b'*k*b/4/PI.  */
    /* When it's not moving, this single dislocation's strain  */
    /* energy integral PE is the SAME as the PE from the work  */
    /* done by the (or any) dipole cut, providing convenience. */
    /***********************************************************/

    double PEmode[MODES];
    double PE; /* prelogarithmic potential energy factor */
    double KE; /* prelogarithmic kinetic energy factor */
    double E;  /* prelogarithmic total energy factor */
} ElasticityStraightDislocation;
#undef MODES

/* To break the isotropic degeneracy */
#define ElasticityStraightDislocationMinAnisotropy (SMALL/1000)

/* Stroh solution of an infinite, straight, uniformly moving dislocation */
void ElasticityStraightMovingDislocation
(V3 m, V3 xi, V3 b, double mass_density, V3 velocity, T4 c,
 ElasticityStraightDislocation *d);

/* Stroh solution of an infinite, straight, static dislocation */
void ElasticityStraightStaticDislocation
(V3 m, V3 xi, V3 b, T4 c, ElasticityStraightDislocation *d);

/* Cartesian radial distance */
#define ElasticityStraightDislocationCartesianR2(d,dx) \
  ( V3LENGTH2(dx) - SQUARE(V3DOT((d)->xi,dx)) )
#define ElasticityStraightDislocationCartesianR(d,dx) \
  sqrt( ElasticityStraightDislocationCartesianR2(d,dx) )

/* Among np atoms s[0..3*np-1], pick out the atom w/ the smallest radial */
/* distance to dislocation at (s0,s1,0) which is also parallel to H[2][] */
int ElasticityStraightDislocationSmallestCartesianR2Atom
(M3 H, double s0, double s1, int np, double *s);

/* Scaled distance squared for a pair of Stroh modes */
/* like the denominator in H&L pp. 455: Eqn. 13-129. */
double ElasticityStraightDislocationStrohR2
(ElasticityStraightDislocation *d, int a, V3 dx);
#define ElasticityStraightDislocationStrohR(d,a,dx) \
  sqrt( ElasticityStraightDislocationStrohR2(d,a,dx) )

/* Peach-Koehler force[] per unit-length of dislocation in */
/* absolute frame at external stress[][] in absolute frame */
void ElasticityStraightDislocation_pkforce
(T2 stress, ElasticityStraightDislocation *d, V3 force);
#define ElasticityStraightDislocation_pkFORCE(stress,d,force,tmp) \
  ( V3mM3(d->b,stress,tmp), V3CROSS(tmp,d->xi,force) )

/* P-K force per unit-length of dislocation resolved in normalized  */
/* absolute direction[] at an external stress[][] in absolute frame */
double ElasticityStraightDislocation_pkresolved
(T2 stress, ElasticityStraightDislocation *d, V3 direction);

/* Peach-Koehler Force[] per unit-length of dislocation in */
/* observer frame at external stress[][] in absolute frame */
void ElasticityStraightDislocation_pkForce
(T2 stress, ElasticityStraightDislocation *d, V3 Force);

/* P-K force per unit-length of dislocation resolved in normalized  */
/* observer Direction[] at an external stress[][] in absolute frame */
double ElasticityStraightDislocation_pkResolved
(T2 stress, ElasticityStraightDislocation *d, V3 Direction);

/* Peach-Koehler Force[] per unit-length of dislocation in */
/* observer frame at external Stress[][] in observer frame */
void ElasticityStraightDislocation_PKForce
(T2 Stress, ElasticityStraightDislocation *d, V3 Force);
#define ElasticityStraightDislocation_PKFORCE(Stress,d,Force) ( \
  (Force)[0] =   (d)->B[0] * Stress[0][1] + \
  (d)->B[1] * Stress[1][1] + (d)->B[2] * Stress[2][1], \
  (Force)[1] = - (d)->B[0] * Stress[0][0] - \
  (d)->B[1] * Stress[1][0] - (d)->B[2] * Stress[2][0], \
  (Force)[2] = 0 )

/* P-K force per unit-length of dislocation resolved in normalized  */
/* observer Direction[] at an external Stress[][] in observer frame */
double ElasticityStraightDislocation_PKResolved
(T2 Stress, ElasticityStraightDislocation *d, V3 Direction);

/* Peach-Koehler force[] per unit-length of dislocation in */
/* absolute frame at external Stress[][] in observer frame */
void ElasticityStraightDislocation_PKforce
(T2 Stress, ElasticityStraightDislocation *d, V3 force);

/* P-K force per unit-length of dislocation resolved in normalized  */
/* absolute direction[] at an external Stress[][] in observer frame */
double ElasticityStraightDislocation_PKresolved
(T2 Stress, ElasticityStraightDislocation *d, V3 direction);

/* Self-energy of a straight static dislocation dipole */
/* separated by dx[] in infinite linear elastic medium */
double ElasticityStraightDislocationStaticDipoleSelfEnergy
(ElasticityStraightDislocation *d, V3 dx);

/* Eself(dipole) = 2*Ecore + 2*PE*log(R/r0) + 2*A(theta): return A(theta); */
/* theta is with respect to input dx[], and xi of d, so A(theta=0) = 0.    */
double ElasticityStraightDislocationStaticDipoleSelfAngularEnergy
(ElasticityStraightDislocation *d, V3 dx, double theta);

/* Total coupling energy between two identical straight */
/* static dislocation dipoles of dx[] separated by R[]. */
double ElasticityStraightDislocationStaticDipoleInteractionEnergy
(ElasticityStraightDislocation *d, V3 dx, V3 R);

/* displacement field: dx[] (input) and u[] are in absolute frame */
void ElasticityStraightDislocation_u
(ElasticityStraightDislocation *d, V3 dx, V3 u);

/* Displacement field: (X,Y,0) input and U[] are in (m,n,xi) frame */
void ElasticityStraightDislocation_U
(ElasticityStraightDislocation *d, double X, double Y, V3 U);

/* displacement gradient field: dx[] (input) and */
/* ux[k][l] := du[k]/dx[l] are in absolute frame */
void ElasticityStraightDislocation_ux
(ElasticityStraightDislocation *d, V3 dx, double ux[3][3]);

/* Displacement Gradient field: (X,Y,0) input and */
/* UX[k][l] := dU[k]/dX[l] are in (m,n,xi) frame. */
void ElasticityStraightDislocation_UX
(ElasticityStraightDislocation *d, double X, double Y, double UX[3][3]);

/* velocity field: dx[] (input) and udot[] are in absolute frame */
void ElasticityStraightDislocation_udot
(ElasticityStraightDislocation *d, V3 dx, V3 udot);

/* Velocity field: (X,Y,0) input and Udot[] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Udot
(ElasticityStraightDislocation *d, double X, double Y, V3 Udot);

/* strain field: dx[] (input) and strain[][] are in absolute frame */
void ElasticityStraightDislocation_strain
(ElasticityStraightDislocation *d, V3 dx, T2 strain);

/* Strain field: (X,Y,0) input and Strain[][] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Strain
(ElasticityStraightDislocation *d, double X, double Y, T2 Strain);

/* stress field: dx[] (input) and stress[][] are in absolute frame */
void ElasticityStraightDislocation_stress
(ElasticityStraightDislocation *d, V3 dx, T2 stress);

/* Stress field: (X,Y,0) input and Stress[][] are in (m,n,xi) frame */
void ElasticityStraightDislocation_Stress
(ElasticityStraightDislocation *d, double X, double Y, T2 Stress);

/* strain energy per unit volume: dx[] (input) is in absolute frame */
double ElasticityStraightDislocation_pote_density
(ElasticityStraightDislocation *d, V3 dx);

/* Strain energy per unit volume: (X,Y,0) input is in (m,n,xi) frame */
double ElasticityStraightDislocation_Pote_Density
(ElasticityStraightDislocation *d, double X, double Y);

/* kinetic energy per unit volume: dx[] (input) is in absolute frame */
double ElasticityStraightDislocation_kine_density
(ElasticityStraightDislocation *d, V3 dx);

/* Kinetic energy per unit volume: (X,Y,0) input is in (m,n,xi) frame */
double ElasticityStraightDislocation_Kine_Density
(ElasticityStraightDislocation *d, double X, double Y);


/* PBCStaticStraightDislocationDipole.c: */

/************************************************************/
/* Note that the observer frame for observing a slip system */
/* is usually the "absolute" frame for later operations.    */
/************************************************************/

/* User-setup */
typedef struct
{
    /* possible inputs */
    M3 H; /* new PBC dimensions in observer/"absolute" frame */
    T2 D; /* H = H0 * (1+D) */
    T2 Stress0; /* local stress at the first dislocation */
    T2 Stress1; /* local stress at the second dislocation */
    T2 Stresscell; /* supercell Virial stress average */

    /* outputs */
    T2 simpleStrain; /* simple strain compared to Hfree[][]: = (D+D^T)/2 */
    double volume;
    M3 HI; /* H^-1[][] */
    double Ebench; /* elastic benchmark energy */
    /* Eatomistic = Ebench + 2*(Ecore-PE*log(r0))*linelength */
    V3 Force0; /* P-K force/length at the first dislocation */
    V3 Force1; /* P-K force/length at the second dislocation */
} ElasticityPBCStaticStraightDislocationDipoleUserSetup;

#define UnknownElasticityPBCStaticStraightDislocationDipoleUserSetup { \
  UnknownElasticityM3, UnknownElasticityM3, UnknownElasticityM3, \
  UnknownElasticityM3, UnknownElasticityM3 }

typedef struct
{
    /** Input: **/
    M3 Hfree;
    /* PBC supercell in a stress-free crystal in observer frame */
    T4 C;
    /* elastic constant of the stress-free crystal in observer frame */
    double s0[2]; /* in-plane reduced coordinates of the first dislocation */
    double s1[2]; /* in-plane reduced coordinates of the second dislocation */
    int images[2]; /* number of Green's function images */

    /** Results: **/
    double linelength; /* (single) dislocation length in supercell */
    double volumefree;
    V3 x0, dx; /* dislocation sites in stress-free crystal */
    ElasticityStraightDislocation d[1];
    /* the dislocations are parallel to H[2][] */
    T2 plasticStrain; /* plastic strain embedded in the Green's supercell */
    T4 S; /* elastic compliance of the stress-free crystal in observer frame */

    /* Green's supercell: the initially stress-free  */
    /* crystal after summing image Green's functions */
    ElasticityPBCStaticStraightDislocationDipoleUserSetup green[1];
    M3 greenu; /* Green's displacement at sfree=(1,0),(0,1) and (0,0) */
    double Egreenwork; /* this affine work breaks conditional convergence */

    /* if H[][] remains Hfree[][] */
    ElasticityPBCStaticStraightDislocationDipoleUserSetup stay[1];
    /* if local stress on the first dislocation is zero */
    ElasticityPBCStaticStraightDislocationDipoleUserSetup z0[1];
    /* if local stress on the second dislocation is zero */
    ElasticityPBCStaticStraightDislocationDipoleUserSetup z1[1];
    /* if supercell Virial stress average is zero: */
    /* this is also the minimum total energy point */
    ElasticityPBCStaticStraightDislocationDipoleUserSetup natural[1];
} ElasticityPBCStaticStraightDislocationDipole;

extern const int
ElasticityPBCStaticStraightDislocationDipoleSufficientImages[2];

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
 ElasticityPBCStaticStraightDislocationDipole *dipole);
#define \
  ElasticityPBCStaticStraightDislocationDipoleASSIGN(Hfree,C,s0,s1,B,dipole) \
  ElasticityPBCStaticStraightDislocationDipoleAssign(Hfree,C,s0,s1,B,\
  INTP(ElasticityPBCStaticStraightDislocationDipoleSufficientImages),dipole)

/* Reduced coordinate transformation */
void ElasticityPBCStaticStraightDislocationDipoleTransform
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, V3 snew);
#define ElasticityPBCStaticStraightDislocationDipoleTRANSFORM(dipole,sfree) \
  ElasticityPBCStaticStraightDislocationDipoleTransform(dipole,sfree,sfree)

/* Local strain in the Green's supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleGreenStrain
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, T2 Strain);

/* Local stress in the Green's supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleGreenStress
(ElasticityPBCStaticStraightDislocationDipole *dipole, V3 sfree, T2 Stress);

/* Set all user-setup inputs as unknown */
void ElasticityPBCStaticStraightDislocationDipoleUserSetupUnknown
(ElasticityPBCStaticStraightDislocationDipoleUserSetup *user);

/* Autocomplete user-setup */
void ElasticityPBCStaticStraightDislocationDipoleUserSetupComplete
(ElasticityPBCStaticStraightDislocationDipole *dipole, 
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user);

/* Strain field in the user-setup supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleUserStrain
(ElasticityPBCStaticStraightDislocationDipole *dipole,
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user,
 V3 sfree, T2 Strain);

/* Stress field in the user-setup supercell: the input is s[] of Hfree[][] */
void ElasticityPBCStaticStraightDislocationDipoleUserStress
(ElasticityPBCStaticStraightDislocationDipole *dipole,
 ElasticityPBCStaticStraightDislocationDipoleUserSetup *user,
 V3 sfree, T2 Stress);

#define \
  ElasticityPBCStaticStraightDislocationDipoleAngularEnergy(dipole,theta) \
  ElasticityStraightDislocationStaticDipoleSelfAngularEnergy((dipole)->d, \
  (dipole)->dx, theta)


/* Homogeneous.c: */

/*************************************/
/* Nonlinear homogeneous deformation */
/*************************************/

/* Calculate the work done by constant external stress sout[][] */
/* from H0[][] to H[][] by straight-path Romberg integration.   */
double HomogeneousWork(M3 H0, M3 sout, M3 H);

#endif  /* _Elasticity_h */
