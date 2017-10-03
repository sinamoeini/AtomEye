/******************************************/
/* libVecMat2: -lScalar -lIO -lm          */
/*                                        */
/* Two-dimensional Euclidean space vector */
/* and matrix library.                    */
/*                                        */
/* May 24, 2001 Ju Li <liju99@mit.edu>    */
/******************************************/

#ifndef _VecMat2_h
#define _VecMat2_h

#include <IO.h>
#include <Scalar.h>

#define ABS2(a,b)  ((b)[0]= ABS((a)[0]),(b)[1]= ABS((a)[1]))
#define FABS2(a,b) ((b)[0]=fabs((a)[0]),(b)[1]=fabs((a)[1]))

typedef double V2[2];
typedef double M2[2][2];
typedef double (*M2P)[2];

/* memory expansion */
#define V2e(a) (a),(a)+1
#define V2E(a) (a)[0],(a)[1]


/* V2.c: */
/*****************************/
/* single vector with scalar */
/*****************************/

/* unary vector */

/* a[] := 0; return a[] */
double *V2zero (double a[2]);
#define V2ZERO(a) ((a)[0]=0, (a)[1]=0)
#define V2EYE(a,x) ((x)[0]=(a), (x)[1]=(a))

/* generate a unit vector: a[i] := 1, a[j!=i] := 0; return a[] */
double *V2unit (double a[2], int i);
#define V2UNIT(a,i) ((a)[0]=0,(a)[1]=0,(a)[i]=1)

#define V2ASSIGN(x0,x1,a) ((a)[0]=(x0),(a)[1]=(x1))

/* generate a[] with 2 independent random components on (0.,1.); return a[] */
double *V2frandom (double a[2]);
#define V2Frandom(a) ((a)[0]=Frandom(),(a)[1]=Frandom())

/* generate a with 2 independent random components on (-0.5,0.5) */
#define V2FRANDOM(a) ((a)[0]=FRANDOM(),(a)[1]=FRANDOM())

/* generate a unit vector of spherically uniform random orientation */
double *V2randomunit (double a[2]);

/* b[] := a[]; then return b[] */
double *V2eqv (double a[2], double b[2]);
#define V2EQV(a,b) ((b)[0]=(a)[0], (b)[1]=(a)[1])

#define V2EQ(a,b) ( ((b)[0]==(a)[0]) && ((b)[1]==(a)[1]) )
#define V2NE(a,b) ( ((b)[0]!=(a)[0]) || ((b)[1]!=(a)[1]) )
#define V2EQZERO(a) ( (0 == (a)[0]) && (0 == (a)[1]) )
#define V2NEZERO(a) ( (0 != (a)[0]) || (0 != (a)[1]) )

#define V2ISTINY(a)  ( V2LENGTH2(a) < TINY*TINY )
#define V2ISSMALL(a) ( V2LENGTH2(a) < SMALL*SMALL )
#define V2ISDELTA(a) ( V2LENGTH2(a) < DELTA*DELTA )

/* b[] := -a[]; then return b[] */
double *V2neg (double a[2], double b[2]);
#define V2NEG(a,b) ((b)[0]=-(a)[0], (b)[1]=-(a)[1])
/* a[] := -a[]; then return a[] */
double *V2Neg (double a[2]);
#define V2NeG(a) V2NEG(a,a)

#define V2SWAP(a,b,tmp) ( SWAP((a)[0],(b)[0],tmp), SWAP((a)[1],(b)[1],tmp) )

#define V2XIN(v,lower_bound,upper_bound) (XIN((v)[0],lower_bound,upper_bound) \
  && XIN((v)[1],lower_bound,upper_bound))

#define V2XOU(v,lower_bound,upper_bound) (XOU((v)[0],lower_bound,upper_bound) \
  || XOU((v)[1],lower_bound,upper_bound))

#define V2NONEED_TRIM(x)  ( (((x)[0]) >= 0.) && (((x)[0]) < 1.) && \
  (((x)[1]) >= 0.) && (((x)[1]) < 1.) )
#define V2NEED_TRIM(x)    ( (((x)[0]) <  0.) || (((x)[0]) >= 1.) || \
  (((x)[1]) <  0.) || (((x)[1]) >= 1.) )
/* change a[] to its own image in [0,1)^2; then return a[] */
#define V2TriM(a) { Trim((a)[0]); Trim((a)[1]); }
/* faster than V2TRIM() when a[] is around 1 */
#define V2TRIM(a,b) ( (b)[0]=TRIM((a)[0]), (b)[1]=TRIM((a)[1]) )
/* faster than V2TriM() when a[] is very large */

/* change a[] to its own image in [0,1)^2; then return a[] */
double *V2Trim (double a[2]);
/* make b[] an image of a[] in [0,1)^2; then return b[] */
double *V2trim (double a[2], double b[2]);

#define V2NONEED_IMAGE(x)  ( ((x[0]) >= -0.5) && ((x[0]) < 0.5) && \
  ((x[1]) >= -0.5) && ((x[1]) < 0.5) )
#define V2NEED_IMAGE(x)    ( ((x[0]) <  -0.5) || ((x[0]) >= 0.5) || \
  ((x[1]) <  -0.5) || ((x[1]) >= 0.5) )

/* change a[] to its own image in [-0.5,0.5)^2; then return a[] */
#define V2ImagE(a) { Image((a)[0]); Image((a)[1]); }
/* faster than V2IMAGE() when a[] is around 1 */
#define V2IMAGE(a,b) (  (b)[0]=IMAGE((a)[0]), (b)[1]=IMAGE((a)[1]) )
/* faster than V2ImagE() when a[] is very large */

/* change a[] to its own image in [-0.5,0.5)^2; then return a[] */
double *V2Image (double a[2]);
/* make b[] image of a[]'s in [-0.5,0.5)^2; then return b[] */
double *V2image (double a[2], double b[2]);

/* scalar & vector */

/* b[i] := a[i] + x, i=0..1 */
#define V2aDD(a,x,b) ((b)[0]=(a)[0]+(x),(b)[1]=(a)[1]+(x))

/* a[i] := a[i] + x, i=0..1 */
#define V2aDd(a,x)   ((a)[0]+=(x),(a)[1]+=(x))

/* b[] := multiplier * a[]; then return b[] */
double *V2mul (double multiplier, double a[2], double b[2]);
#define V2MUL(multiplier,a,b) ( \
  (b)[0] = (multiplier)*(a)[0], (b)[1] = (multiplier)*(a)[1] )

#define V2muL(multiplier0,multiplier1,multiplier2,a,b)  ( \
  (b)[0] = (multiplier0) * (a)[0], (b)[1] = (multiplier1) * (a)[1] )

/* a[] := multiplier * a[]; then return a[] */
double *V2Mul (double multiplier, double a[2]);
#define V2MuL(multiplier,a) ( \
  (a)[0] *= (multiplier), (a)[1] *= (multiplier) )

#define V2mUL(multiplier0,multiplier1,a) ( \
  (a)[0] *= (multiplier0), (a)[1] *= (multiplier1) )

/* b[] := a[] / divisor; then return b[] */
double *V2div (double a[2], double divisor, double b[2]);
#define V2DIV(a,divisor,b) ( \
  (b)[0] = (a)[0] / (divisor), (b)[1] = (a)[1] / (divisor) )

#define V2diV(a,divisor0,divisor1,b) ( \
  (b)[0] = (a)[0] / (divisor0), (b)[1] = (a)[1] / (divisor1) )
/* c[] := a[] ./ b[] */
#define V2ddiV(a,b,c) V2diV(a,b[0],b[1],b[2],c)
/* a[] := a[] ./ b[] */
#define V2ddiv(a,b) ( (a)[0]/=(b)[0], (a)[1]/=(b)[1] )

/* a[] := a[] / divisor; then return a[] */
double *V2Div (double a[2], double divisor);
#define V2DiV(a,divisor) \
  ((a)[0]/=(divisor), (a)[1]/=(divisor))

#define V2dIV(a,divisor0,divisor1,divisor2)  ( \
  (a)[0]/=(divisor0), (a)[1]/=(divisor1) )

/* length of a[] */
double V2length (double a[2]);
#define V2LENGTH(a) (sqrt((a)[0]*(a)[0]+(a)[1]*(a)[1]))

/* squared length of a[] */
double V2length2 (double a[2]);
#define V2LENGTH2(a) ((a)[0]*(a)[0]+(a)[1]*(a)[1])

/* arbitrary exponent-th norm of a[] */
double V2norm (double a[2], double exponent);

/* the largest component of a[] in absolute value */
#define V2INFNORM(a,infnorm) { infnorm = ABS(a[0]); \
  if ( a[1] > infnorm ) infnorm = a[1]; \
  else if ( -a[1] > infnorm ) infnorm = -a[1]; }

/* normalize a[] to unit length; then return a[] */
double *V2normalize (double a[2]);
#define V2NORMALIZE(a,r) (r=V2LENGTH(a),a[0]/=r,a[1]/=r)

/* normalize a[] to unit length; then return its original length */
double V2Normalize (double a[2]);

/* b[] := a[] / |a[]|, return b[] */
double *V2direction (double a[2], double b[2]);
#define V2DIRECTION(a,b,r) (r=V2LENGTH(a),V2DIV(a,r,b))

/* b[] := a[] / |a[]|, return |a[]| */
double V2Direction (double a[2], double b[2]);
#define V2DirectioN(a,b,r) (r=V2LENGTH(a),V2DIV(a,r,b),r)

/* assign index array idx[] := 0..1 */
#define SEQ2(idx) ( idx[0]=0, idx[1]=1 )

/* return 0..1 permutation idx[] such that length[idx[]] is ascending */
int *V2sort (double length[2], int idx[2]);
#define V2SORT(length,idx,tmpi) { SEQ2(idx); \
  if (length[idx[0]]>length[idx[1]]) SWAP(idx[0],idx[1],tmpi); \
  if (length[idx[1]]>length[idx[2]]) SWAP(idx[1],idx[2],tmpi); }

/* return 0..1 random permutation index idx[] */
int *I2randompermutation (int idx[2]);
#define I2RANDOMPERMUTATION(idx) ( idx[0] = Fran(0,1), idx[1] = 1 - idx[0] )

/* return 0..1 permutation idx[] such that length[idx[]]  */
/* is non-decreasing. Equal length indices are randomized */
int *V2randomsort (double length[2], int idx[2]);

/* return 0..1 permutation idx[] such that length[idx[]] is     */
/* approximately increasing within epsilon. idx[] is randomized */
int *V2Randomsort (double length[2], int idx[2], double epsilon);



/* VV2.c: */

/*******************************/
/* vector & vector with scalar */
/*******************************/


/* vector & vector */

/* c[] := a[] + b[]; then return c[] */
double *V2add (double a[2], double b[2], double c[2]);
#define V2ADD(a,b,c) ( (c)[0] = (a)[0] + (b)[0], (c)[1] = (a)[1] + (b)[1] )

/* b[] := a[] + b[]; then return b[] */
double *V2Add (double a[2], double b[2]);
#define V2AdD(a,b) ( (b)[0]+=(a)[0], (b)[1]+=(a)[1] )

/* c[] := a[] - b[]; then return c[] */
double *V2sub (double a[2], double b[2], double c[2]);
#define V2SUB(a,b,c) ( (c)[0] = (a)[0] - (b)[0], (c)[1] = (a)[1] - (b)[1] )

/* a[] := a[] - b[]; then return a[] */
double *V2Sub (double a[2], double b[2]);
#define V2SuB(a,b) ( (a)[0]-=(b)[0], (a)[1]-=(b)[1] )

/* dot product of a[] and b[] */
double V2dot (double a[2], double b[2]);
#define V2DOT(a,b) ((a)[0]*(b)[0]+(a)[1]*(b)[1])

/* determine if a[] and b[] are collinear: b[] = c*a[], c=(-inf,inf) */
#define V2COLLINEAR(a,b,tmp) ( tmp = V2DOT(a,b), \
  tmp *= tmp/V2LENGTH2(a)/V2LENGTH2(b), tmp--, ((tmp>-TINY)&&(tmp<TINY)) )
/* assuming a[] and b[] have non-zero length */

/* vector & vector & scalar */

/* c[] := a[] + multiplier * b[]; then return c[] */
double *V2addmul (double a[2], double multiplier, double b[2], double c[2]);
#define V2ADDMUL(a,multiplier,b,c) ( \
  (c)[0] = (a)[0]+(multiplier)*(b)[0], (c)[1] = (a)[1]+(multiplier)*(b)[1] )

/* c[] := aa * a[] + bb * b[] */
#define V2ADDMULMUL(aa,a,bb,b,c) ( \
  (c)[0] = (aa)*(a)[0]+(bb)*(b)[0], (c)[1] = (aa)*(a)[1]+(bb)*(b)[1] )

/* c[] := a[] + b[] / divisor; then return c[] */
double *V2adddiv (double a[2], double b[2], double divisor, double c[2]);
#define V2ADDDIV(a,b,divisor,c) ( \
  (c)[0] = (a)[0] + (b)[0] / (divisor), (c)[1] = (a)[1] + (b)[1] / (divisor) )

/* b[] := b[] + multiplier * a[]; then return b[] */
double *V2ADDmul (double multiplier, double a[2], double b[2]);
#define V2ADDmuL(multiplier,a,b) ( \
  (b)[0] += (multiplier) * (a)[0], (b)[1] += (multiplier) * (a)[1] )

/* b[] := b[] + a[] / divisor; then return b[] */
double *V2ADDdiv (double a[2], double divisor, double b[2]);
#define V2ADDdiV(a,divisor,b) ( \
  (b)[0] += (a)[0] / (divisor), (b)[1] += (a)[1] / (divisor) )

/* c[] := a[] - multiplier * b[]; then return c[] */
double *V2submul (double a[2], double multiplier, double b[2], double c[2]);
#define V2SUBMUL(a,multiplier,b,c) ( \
  (c)[0] = (a)[0]-(multiplier)*(b)[0], (c)[1] = (a)[1]-(multiplier)*(b)[1] )

/* c[] := multiplier * a[] - b[]; */
#define V2MULSUB(multiplier,a,b,c) ( \
  (c)[0]=(multiplier)*(a)[0]-(b)[0], (c)[1] = (multiplier)*(a)[1]-(b)[1] )

/* c[] := a[] - b[] / divisor; then return c[] */
double *V2subdiv (double a[2], double b[2], double divisor, double c[2]);
#define V2SUBDIV(a,b,divisor,c) ( \
  (c)[0] = (a)[0] - (b)[0] / (divisor), (c)[1] = (a)[1] - (b)[1] / (divisor) )

/* a[] := a[] - multiplier * b[]; then return a[] */
double *V2SUBmul (double a[2], double multiplier, double b[2]);
#define V2SUBmuL(a,multiplier,b) ( \
  (a)[0] -= (multiplier) * (b)[0], (a)[1] -= (multiplier) * (b)[1] )

/* a[] := a[] - b[] / divisor; then return a[] */
double *V2SUBdiv (double a[2], double b[2], double divisor);
#define V2SUBdiV(a,b,divisor) ( \
  (a)[0] -= (b)[0] / (divisor), (a)[1] -= (b)[1] / (divisor) )

/* c[] := part of a[] that is perpendicular to b[]; return c[] */
double *V2perpendicular (double a[2], double b[2], double c[2]);
#define V2PERPENDICULAR(a,b,c,tmp) \
  ( tmp=V2DOT(a,b)/V2LENGTH2(b), V2SUBMUL(a, tmp, b, c) )

/* a[] := part of a[] that is perpendicular to b[]; return modified a[] */
double *V2Perpendicular (double a[2], double b[2]);
#define V2PerpendiculaR(a,b,tmp) \
  ( tmp=V2DOT(a,b)/V2LENGTH2(b), V2SUBmuL(a, tmp, b) )


/* vector & vector & vector */

/* d[] := a[] + b[] - c[] */
#define V2ADDSUB(a,b,c,d) ( (d)[0] = (a)[0] + (b)[0] - (c)[0], \
  (d)[1] = (a)[1] + (b)[1] - (c)[1] )
/* d[] := a[] - b[] - c[] */
#define V2SUBSUB(a,b,c,d) ( (d)[0] = (a)[0] - (b)[0] - (c)[0], \
  (d)[1] = (a)[1] - (b)[1] - (c)[1] )


/* M2.c: */

#define GuineaPigM2 \
  {{-0.43256481152822,0.12533230647483},{-1.66558437823810,0.28767642035855}}
#define GuineaPigSymmetricM2 \
  {{-1.14647135068146,1.19003983364755},{1.19003983364755,-0.03763327659332}}
#define GuineaPigSymmetricPositiveDefiniteM2 \
  {{0.14198038281628,-0.07835326351823},{-0.07835326351823,0.55727075019712}}
/* O * O' = O' * O = I */
#define GuineaPigOrthonormalM2 \
  {{-0.58768503408637,0.80908979768064},{0.80908979768064,0.58768503408637}}
/* det|O| = 1 */
#define GuineaPigOrthonormalRightHandedM2 \
  {{-0.59419785867637,0.80431890736474},{-0.80431890736474,-0.59419785867637}}
/* det|O| = -1 */
#define GuineaPigOrthonormalLeftHandedM2 \
  {{0.71497723335124,-0.69914773531022},{-0.69914773531022,-0.71497723335124}}

/* memory expansion */
#define M2e(A) A[0][0],A[0][1],A[1][0],A[1][1]

/**********************************************/
/* single matrix & scalar & generating vector */
/**********************************************/

/* matrix */

/* symmetric 2x2 -> 0-2 index conversion */
#define M2OD(i,j) ((3-MIN(i,j))*MIN(i,j)/2+MAX(i,j))

#define M2ISSYMMETRIC(A) EQUAL((A)[0][1],(A)[1][0])

#define M2ASSIGN(x00,x01,x10,x11,A) ( A[0][0]=(x00),A[0][1]=(x01), \
  A[1][0]=(x10),A[1][1]=(x11) )

/* A[0][] := x0[], A[1][] := x1[] */
#define M2Assign(x0,x1,x2,A) M2ASSIGN(x0[0],x0[1], x1[0],x1[1], A)

/* A[][0] := y0[], A[][1] := y1[] */
#define M2assign(y0,y1,y2,A) M2ASSIGN(y0[0],y1[0], y0[1],y1[1], A)

/* C[][] := a[]' * b[] */
#define M2ASSIGNV2V2(a,b,C) ( \
  (C)[0][0] = (a)[0] * (b)[0], \
  (C)[0][1] = (a)[0] * (b)[1], \
  (C)[1][0] = (a)[1] * (b)[0], \
  (C)[1][1] = (a)[1] * (b)[1] )

/* generate a zero matrix A[][] := 0 */
void M2zero (double A[2][2]);
#define M2ZERO(A) ( A[0][0]=0, A[0][1]=0, A[1][0]=0, A[1][1]=0 )

/* generate an identity matrix A[][] := I[][] */
void M2identity (double A[2][2]);
#define M2IDENTITY(A) ( A[0][0]=1, A[0][1]=0, A[1][0]=0, A[1][1]=1 )

/* generate A[][] with 4 independent random components on (0.,1.) */
void M2frandom (double A[2][2]);
#define M2Frandom(A) ( A[0][0] = Frandom(), A[0][1] = Frandom(), \
  A[1][0] = Frandom(), A[1][1] = Frandom() )

/* generate A[][] with 4 independent random components on (-0.5,0.5) */
#define M2FRANDOM(A) ( A[0][0] = FRANDOM(), A[0][1] = FRANDOM(), \
  A[1][0] = FRANDOM(), A[1][1] = FRANDOM() )

/* B[][] := A[][] */
void M2eqv (double A[2][2], double B[2][2]);
#define M2EQV(A,B) ( B[0][0] = A[0][0], B[0][1] = A[0][1], \
  B[1][0] = A[1][0], B[1][1] = A[1][1] )

#define M2EQ(A,B) ( (B[0][0] == A[0][0]) && (B[0][1] == A[0][1]) && \
  (B[1][0] == A[1][0]) && (B[1][1] == A[1][1]) )
#define M2EQZERO(A) ( (0 == A[0][0]) && (0 == A[0][1]) && \
  (0 == A[1][0]) && (0 == A[1][1]) )

#define M2NE(A,B) ( (B[0][0] != A[0][0]) || (B[0][1] != A[0][1]) || \
  (B[1][0] != A[1][0]) || (B[1][1] != A[1][1]) )
#define M2NEZERO(A) ( (0 != A[0][0]) || (0 != A[0][1]) || \
  (0 != A[1][0]) || (0 != A[1][1]) )

#define M2ORTHOGONAL(A) ( \
  ISTINY(V2DOT(A[0],A[1])) && ISTINY(V2DOT(A[0],A[2])) )
#define M2NONORTHOGONAL(A) ( \
  NOTTINY(V2DOT(A[0],A[1])) || NOTTINY(V2DOT(A[0],A[2])) )

/* arbitrary exponent-th norm of A[][] (all components) */
double M2norm (double A[2][2], double exponent);

/* return the largest component of A[][] in absolute value */
double M2infnorm (double A[2][2]);
#define M2INFNORM(A,infnorm) { infnorm = ABS( A[0][0] ); \
  if ( A[0][1] > infnorm ) infnorm = A[0][1]; \
  else if ( -A[0][1] > infnorm ) infnorm = -A[0][1]; \
  if ( A[1][0] > infnorm ) infnorm = A[1][0]; \
  else if ( -A[1][0] > infnorm ) infnorm = -A[1][0]; \
  if ( A[1][1] > infnorm ) infnorm = A[1][1]; \
  else if ( -A[1][1] > infnorm ) infnorm = -A[1][1]; }

/* return sqrt{ \sum_{i=0..1} \sum_{j=0..1} |A_ij|^2 } */
double M22norm (double A[2][2]);
#define M22NORM2(A) (V2LENGTH2(A[0])+V2LENGTH2(A[1]))
#define M22NORM(A) sqrt(V2LENGTH2(A[0])+V2LENGTH2(A[1]))

/* B[][] := -A[][] */
void M2neg (double A[2][2], double B[2][2]);
#define M2NEG(A,B) (  B[0][0] = -A[0][0], B[0][1] = -A[0][1], \
  B[1][0] = -A[1][0], B[1][1] = -A[1][1] )

/* A[][] := -A[][] */
void M2Neg (double A[2][2]);
#define M2NeG(A) M2NEG(A,A)

/* B[][] := A[][]' */
void M2transpose (double A[2][2], double B[2][2]);
/* !! please use M2Transpose() for A[][] := A[][]' !! */
#define M2TRANSPOSE(A,B) ( B[0][0] = A[0][0], B[0][1] = A[1][0], \
  B[1][0] = A[0][1], B[1][1] = A[1][1] )
/* !! please use M2Transpose() for A[][] := A[][]' !! */

/* A[][] := A[][]' */
void M2Transpose (double A[2][2]);
#define M2TransposE(A,tmp) ( tmp = A[0][1], A[0][1] = A[1][0], A[1][0] = tmp )

/* B[][] := (A[][]+A[][]')/2 */
void M2symmetrize (double A[2][2], double B[2][2]);
/* please use M2Symmetrize() for self-application */
#define M2SYMMETRIZE(A,B) ( B[0][0] = A[0][0], B[1][1] = A[1][1], \
  B[0][1] = B[1][0] = ( A[1][0] + A[0][1] ) / 2. )

/* A[][] := (A[][]+A[][]')/2 */
#define M2Symmetrize(A) ( A[0][1] = A[1][0] = (A[1][0]+A[0][1])/2. )

/* A[][] := A[][]+A[][]' */
#define M2Symmetrize2(A) ( A[0][0]*=2, A[1][1]*=2, A[2][2]*=2, \
  A[0][1] = A[1][0] = A[1][0]+A[0][1] )

/* return the trace of A[][] */
double M2Tr (double A[2][2]);
#define M2TR(A) (A[0][0]+A[1][1])

/* Tr(AB) */
#define M2TRPROD(A,B) (B[0][0]*A[0][0] + B[0][1]*A[1][0] + \
  B[1][0]*A[0][1] + B[1][1]*A[1][1] )

/* B[][] := trace(A[][]) / 2 * I[][] */
double M2trace (double A[2][2], double B[2][2]);
#define M2TRACE(A,B,trace) ( (trace) = A[0][0]+A[1][1], \
  B[0][0] = (trace)/2., B[0][1] = 0., \
  B[1][0] = 0., B[1][1] = (trace)/2. )

/* A[][] := trace(A[][])/2 * I[][] */
double M2Trace (double A[2][2]);
#define M2TracE(A,trace) ( (trace) = A[0][0]+A[1][1], \
  A[0][0] = (trace)/2., A[0][1] = 0., \
  A[1][0] = 0., A[1][1] = (trace)/2. )

/* B[][] := A[][] - trace(A[][])/2 * I[][]; return trace(A) */
double M2traceless (double A[2][2], double B[2][2]);
#define M2TRACELESS(A,B,trace) ( (trace) = A[0][0]+A[1][1], \
  B[0][0] = A[0][0] - (trace)/2., B[0][1] = A[0][1], \
  B[1][0] = A[1][0], B[1][1] = A[1][1] - (trace)/2. )

/* A[][] := A[][] - trace(A[][])/2 * I[][]; return original trace(A) */
double M2Traceless(double A[2][2]);
#define M2TracelesS(A,trace) ( (trace) = A[0][0]+A[1][1], \
  A[0][0] -= (trace)/2., A[1][1] -= (trace)/2. )

/* decompose A[][] to b*I[][] + C[][], where b := trace(A)/2; return b */
double M2tracedecompose (double A[2][2], double B[2][2], double C[2][2]);
#define M2TRACEDECOMPOSE(A,B,C,trace) ( (trace) = A[0][0]+A[1][1], \
  B[0][0] = (trace)/2., B[0][1] = 0., B[1][0] = 0., \
  B[1][1] = (trace)/2., C[0][0] = A[0][0] - (trace)/2., C[0][1] = A[0][1], \
  C[1][0] = A[1][0], C[1][1] = A[1][1] - (trace)/2. )

/* matrix & scalar */

/* generate an identity matrix A[][] := a x I[][] */
void M2Identity (double a, double A[2][2]);
#define M2IdentitY(a,A) ( A[0][0]=(a), A[0][1]=0., \
  A[1][0] = 0., A[1][1] = (a) )

/* B[][] := multiplier * A[][] */
void M2multiply (double multiplier, double A[2][2], double B[2][2]);
#define M2MULTIPLY(multiplier,A,B) (B[0][0] = (multiplier) * A[0][0], \
  B[0][1] = (multiplier) * A[0][1], B[1][0] = (multiplier) * A[1][0], \
  B[1][1] = (multiplier) * A[1][1])

/* A[][] := multiplier * A[][] */
void M2Multiply (double multiplier, double A[2][2]);
#define M2MultiplY(multiplier,A) (A[0][0] *= (multiplier), \
  A[0][1] *= (multiplier), A[1][0] *= (multiplier), A[1][1] *= (multiplier) )

/* B[][] := A[][] / divisor */
void M2divide (double A[2][2], double divisor, double B[2][2]);
#define M2DIVIDE(A,divisor,B) (  B[0][0] = A[0][0] / (divisor), \
  B[0][1] = A[0][1] / (divisor), B[1][0] = A[1][0] / (divisor), \
  B[1][1] = A[1][1] / (divisor) )

/* A[][] := A[][] / divisor */
void M2Divide (double A[2][2], double divisor);
#define M2DividE(A,divisor) ( A[0][0] /= (divisor), A[0][1] /= (divisor), \
  A[1][0] /= (divisor), A[1][1] /= (divisor) )

/* B[][] := A[][] + a * I */
void M2adddiag (double A[2][2], double a, double B[2][2]);
#define M2ADDDIAG(A,a,B) ( B[0][0] = A[0][0] + (a), B[0][1] = A[0][1], \
  B[1][0] = A[1][0], B[1][1] = A[1][1] + (a) )

/* B[][] := a * A[][] + x * I */
#define M2MULADDDIAG(a,A,x,B) ( \
  B[0][0] = (a)*A[0][0]+(x), B[0][1] = (a)*A[0][1], \
  B[1][0] = (a)*A[1][0], B[1][1] = (a)*A[1][1]+(x) )

/* A[][] := A[][] + a * I */
void M2Adddiag (double A[2][2], double a);
#define M2AdddiaG(A,a) ( A[0][0] += (a), A[1][1] += (a) )

/* B[][] := A[][] - a * I */
void M2subdiag (double A[2][2], double a, double B[2][2]);
#define M2SUBDIAG(A,a,B) ( B[0][0] = A[0][0] - (a), \
  B[0][1] = A[0][1], B[1][0] = A[1][0], B[1][1] = A[1][1] - (a) )

/* A[][] := A[][] - a * I */
void M2Subdiag (double A[2][2], double a);
#define M2SubdiaG(A,a) ( A[0][0] -= (a), A[1][1] -= (a) )


/* matrix & generating vector */

/* generate a diagonal matrix A[i][i] := a_i, others = 0 */
#define M2diagonal(a_0,a_1,A) ( \
  A[0][0] = (a_0), A[0][1] = 0., A[1][0] = 0., A[1][1] = (a_1) )
#define M2Diagonal(a,A) M2diagonal(a,a,A)

/* generate a diagonal matrix A[i][i] := a[i], others = 0 */
#define M2DIAGONAL(a,A) ( \
  A[0][0] = (a)[0], A[0][1] = 0., A[1][0] = 0., A[1][1] = (a)[1] )


/* VM2.c: */

/************************/
/* row space operations */
/************************/

/* rowlength[i] := | A[i][] |; returns rowlength[] */
double *M2rowlengths (double A[2][2], double rowlength[2]);
#define M2ROWLENGTHS(A,rowlength) ( \
  rowlength[0] = V2LENGTH(A[0]), rowlength[1] = V2LENGTH(A[1]) )

/* returns the maximum Euclidean length of A[0][], A[1][] */
double M2maxrowlength (double A[2][2]);

/* c[] := a[] * B[][]; then return c[] */
double *V2mulM2 (double a[2], double B[2][2], double c[2]);
#define V2mM2(a,B,c) ( \
  (c)[0] = (a)[0]*B[0][0] + (a)[1]*B[1][0], \
  (c)[1] = (a)[0]*B[0][1] + (a)[1]*B[1][1] )

#define V2M2LENGTH2(ds,H,dx) ( V2mM2(ds,H,dx), (dx)[2]=V2LENGTH2(dx) )
#define V2M2LENGTH(ds,H,dx)  ( V2mM2(ds,H,dx), (dx)[2]=V2LENGTH(dx) )

/* x*A*y' */
#define V2ADOT(x,A,y,tmp) (V2mM2(x,A,tmp),V2DOT(tmp,y))

/* a[] := a[] * B[][]; then return a[] */
double *V2MULM2 (double a[2], double B[2][2]);
#define V2MM2(a,B,tmp) (V2mM2(a,B,tmp), V2EQV(tmp,a))

/* a[] := a[] * B[][] / divisor */
#define V2MM2DIV(a,B,divisor,tmp) (V2mM2(a,B,tmp), V2DIV(tmp,divisor,a))

/* d[] := a[] + b[] * C[][]; then return d[] */
double *V2addmulM2 (double a[2], double b[2], double C[2][2], double d[2]);
#define V2addmM2(a,b,C,d) ( \
  (d)[0] = (a)[0] + (b)[0]*C[0][0] + (b)[1]*C[1][0], \
  (d)[1] = (a)[1] + (b)[0]*C[0][1] + (b)[1]*C[1][1] )

/* a[] := a[] + b[] * C[][]; then return a[] */
double *V2ADDmulM2 (double a[2], double b[2], double C[2][2]);
#define V2ADDmM2(a,b,C,d) ( \
  (a)[0] += (b)[0]*C[0][0] + (b)[1]*C[1][0], \
  (a)[1] += (b)[0]*C[0][1] + (b)[1]*C[1][1] )

/* d[] := a[] - b[] * C[][]; then return d[] */
double *V2submulM2 (double a[2], double b[2], double C[2][2], double d[2]);
#define V2submM2(a,b,C,d) ( \
  (d)[0] = (a)[0] - (b)[0]*C[0][0] - (b)[1]*C[1][0], \
  (d)[1] = (a)[1] - (b)[0]*C[0][1] - (b)[1]*C[1][1] )

/* d[] := b[] * C[][] - a[]; then return d[] */
double *V2mulsubM2 (double b[2], double C[2][2], double a[2], double d[2]);
#define V2msubM2(b,C,a,d) ( \
  (d)[0] = (b)[0]*C[0][0] + (b)[1]*C[1][0], \
  (d)[1] = (b)[0]*C[0][1] + (b)[1]*C[1][1] )

/* a[] := a[] - b[] * C[][]; then return a[] */
double *V2SUBmulM2 (double a[2], double b[2], double C[2][2]);
#define V2SUBmM2(a,b,C) ( \
  (a)[0] -= (b)[0]*C[0][0] + (b)[1]*C[1][0], \
  (a)[1] -= (b)[0]*C[0][1] + (b)[1]*C[1][1] )


/* MV2.c: */

/***************************/
/* column space operations */
/***************************/

/* columnlength[i] := |A[][i]|; returns columnlength[] */
double *M2columnlengths (double A[2][2], double columnlength[2]);
#define M2COLUMNLENGTHS(A,columnlength) ( \
  (columnlength)[0] = DISTANCE2D(A[0][0],A[1][0]), \
  (columnlength)[1] = DISTANCE2D(A[0][1],A[1][1]) )

/* returns the maximum Euclidean length of A[][0], A[][1], A[][2] */
double M2maxcolumnlength (double A[2][2]);

/* column[] := A[][i]; return column[] */
double *M2column (double A[2][2], int i, double column[2]);
#define M2COLUMN(A,i,column) \
  ( (column)[0] = A[0][i], (column)[1] = A[1][i] )

/* c[] := A[][] * b[]; then return c[] */
double *M2mulV2 (double A[2][2], double b[2], double c[2]);
#define M2mV2(A,b,c) ( \
  (c)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1], \
  (c)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] )

/* b[] := A[][] * b[]; then return b[] */
double *M2MULV2 (double A[2][2], double b[2]);
#define M2MV2(A,b,tmp) ( M2mV2(A,b,tmp), V2EQV(tmp,b) )

/* d[] := A[][] * b[] + c[]; then return d[] */
double *M2muladdV2 (double A[2][2], double b[2], double c[2], double d[2]);
#define M2maddV2(A,b,c,d) ( \
  (d)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1], \
  (d)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] )

/* c[] := c[] + A[][] * b[]; then return c[] */
double *M2mulADDV2 (double A[2][2], double b[2], double c[2]);
#define M2mADDV2(A,b,c) ( \
  (c)[0] += A[0][0]*(b)[0] + A[0][1]*(b)[1], \
  (c)[1] += A[1][0]*(b)[0] + A[1][1]*(b)[1] )

/* d[] := A[][] * b[] - c[]; then return d[] */
double *M2mulsubV2 (double A[2][2], double b[2], double c[2], double d[2]);
#define M2msubV2(A,b,c,d) ( \
  (d)[0] = A[0][0]*(b)[0] + A[0][1]*(b)[1] - (c)[0], \
  (d)[1] = A[1][0]*(b)[0] + A[1][1]*(b)[1] - (c)[1] )

/* d[] := c[] - A[][] * b[]; then return d[] */
double *M2submulV2 (double c[2], double A[2][2], double b[2], double d[2]);
#define M2submV2(c,A,b,d) ( \
  (d)[0] = (c)[0] - A[0][0]*(b)[0] - A[0][1]*(b)[1], \
  (d)[1] = (c)[1] - A[1][0]*(b)[0] - A[1][1]*(b)[1] )

/* c[] := c[] - A[][] * b[]; then return c[] */
double *M2mulSUBV2 (double A[2][2], double b[2], double c[2]);
#define M2mSUBV2(A,b,c) ( \
  (c)[0] -= A[0][0]*(b)[0] + A[0][1]*(b)[1], \
  (c)[1] -= A[1][0]*(b)[0] + A[1][1]*(b)[1] )


/* M2diag.c */

/*********************************************************/
/* Matrix inversion and symmetric matrix diagonalization */
/*********************************************************/

/* B[][] := A[][]^-1; return det(A) */
double M2inv (double A[2][2], double B[2][2]);

/* A[][] := A[][]^-1; return original det(A) */
double M2Inv (double A[2][2]);

#define M2DETERMINANT(A) ( A[0][0]*A[1][1] - A[1][0]*A[0][1] )
#define M2VOLUME(A) fabs(M2DETERMINANT(A))

/* B[][] := A[][]^-1; determinant := det(A) */
#define M2INV(A,B,determinant) ( \
  (determinant) = M2DETERMINANT(A), \
  B[0][0] =  A[1][1] / (determinant), \
  B[1][1] =  A[0][0] / (determinant), \
  B[1][0] = -A[1][0] / (determinant), \
  B[0][1] = -A[0][1] / (determinant) )

#define M2InV(A,B,volume) { M2INV(A,B,volume); \
  if ((volume) < 0) (volume) = -(volume); }

/* B[][] := A[][]^-T; determinant := det(A) */
#define M2INVTRANSPOSE(A,B,determinant) ( \
  (determinant) = M2DETERMINANT(A), \
  B[0][0] =  A[1][1] / (determinant), \
  B[1][1] =  A[0][0] / (determinant), \
  B[1][0] = -A[0][1] / (determinant), \
  B[0][1] = -A[1][0] / (determinant) )

/******************************************************************/
/* Diagonalize 2x2 real symmetric matrix: A = V^T*Diag(eigval)*V, */
/* eigval[i] will be stored in ascending order, i = 0..1; the     */
/* corresponding eigenvectors V[i][] form a right-handed system.  */
/* Return index i of the eigenvalue of largest absolute value.    */
/******************************************************************/
int M2diag (double A[2][2], double eigval[2], double Q[2][2]);

/***********************************************************************/
/* Find a linear combination of columns A[][0], A[][1], A[][2] so that */
/* A[][]*c[]=0. c[] is normalized and "random" if possible. return c[] */
/***********************************************************************/
double *M2nullvector (double A[2][2], double c[2]);

/* A = M*R, where M is symmetric, R (do not return) is orthogonal */
void M2RightPolarDecompose (M2 A, M2 M);
/* A = M*R, where M is symmetric, R is orthogonal */
void M2RightPolarDECOMPOSE (M2 A, M2 M, M2 R);

/* A = L*M, where L (do not return) is orthogonal, M is symmetric*/
void M2LeftPolarDecompose (M2 A, M2 M);
/* A = L*M, where L is orthogonal, M is symmetric */
void M2LeftPolarDECOMPOSE (M2 A, M2 L, M2 M);

/* MM2.c */

/********************************/
/* matrix & matrix with scalars */
/********************************/

/* matrix & matrix */

/* C[][] := A[][] + B[][] */
void M2add (double A[2][2], double B[2][2], double C[2][2]);
#define M2ADD(A,B,C) (C[0][0] = A[0][0] + B[0][0], \
  C[0][1] = A[0][1] + B[0][1], C[1][0] = A[1][0] + B[1][0], \
  C[1][1] = A[1][1] + B[1][1])

/* B[][] := A[][] + B[][] */
#define M2AdD(A,B) ( B[0][0] += A[0][0], B[0][1] += A[0][1], \
  B[1][0] += A[1][0], B[1][1] += A[1][1] )

/* D[][] := A[][] + b[]' * c[] */
#define M2ADDV2V2(A,b,c,D) ( \
  (D)[0][0] = (A)[0][0] + (b)[0] * (c)[0], \
  (D)[0][1] = (A)[0][1] + (b)[0] * (c)[1], \
  (D)[1][0] = (A)[1][0] + (b)[1] * (c)[0], \
  (D)[1][1] = (A)[1][1] + (b)[1] * (c)[1] )

/* C[][] := C[][] + a[]' * b[] */
#define M2AdDV2V2(a,b,C) ( \
  (C)[0][0] += (a)[0] * (b)[0], \
  (C)[0][1] += (a)[0] * (b)[1], \
  (C)[1][0] += (a)[1] * (b)[0], \
  (C)[1][1] += (a)[1] * (b)[1] )

/* C[][] := A[][] - B[][] */
void M2sub (double A[2][2], double B[2][2], double C[2][2]);
#define M2SUB(A,B,C) ( C[0][0] = A[0][0] - B[0][0], \
  C[0][1] = A[0][1] - B[0][1], C[1][0] = A[1][0] - B[1][0], \
  C[1][1] = A[1][1] - B[1][1] )

/* A[][] := A[][] - B[][] */
#define M2SuB(A,B) ( A[0][0] -= B[0][0], A[0][1] -= B[0][1], \
  A[1][0] -= B[1][0], A[1][1] -= B[1][1] )

/* C[][] := A[][] * B[][] */
void M2mul (double A[2][2], double B[2][2], double C[2][2]);
#define M2MUL(A,B,C) ( \
  C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0], \
  C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1], \
  C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0], \
  C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] )

/* C[][] := A[][] * A'[][] */
#define M2MULT(A,C) ( \
  C[0][0] = A[0][0]*A[0][0] + A[0][1]*A[0][1], \
  C[1][0] = C[0][1] = A[0][0]*A[1][0] + A[0][1]*A[1][1], \
  C[1][1] = A[1][0]*A[1][0] + A[1][1]*A[1][1] )

/* C[][] := A'[][] * A[][] */
#define M2TMUL(A,C) ( \
  C[0][0] = A[0][0]*A[0][0] + A[1][0]*A[1][0], \
  C[1][0] = C[0][1] = A[0][0]*A[0][1] + A[1][0]*A[1][1], \
  C[1][1] = A[0][1]*A[0][1] + A[1][1]*A[1][1] )

/* D[][] := A[][] * B[][] * C[][] */
void M2mul2 (double A[2][2], double B[2][2], double C[2][2], double D[2][2]);
#define M2MUL2(A,B,C,D,TMP) (M2MUL(A,B,TMP),M2MUL(TMP,C,D))

/* matrix & matrix & scalar */


/* strain.c: */

/****************************/
/* Properties of H matrices */
/****************************/

/* calibrated to pure pressure */
#define SymmetricM2HydroInvariant(A) (M2TR(A)/2)
/* calibrated to single shear */
#define SymmetricM2MisesInvariant(A) \
  sqrt( SQUARE(A[0][1]) + SQUARE(A[0][0]-A[1][1])/4 )

/* Crystallographic Notation: http://spot.colorado.edu/~smyth/G20102.html */
typedef struct
{
    double a;
    double b;
    double gamma;  /* in degrees */
} Crystallographic2D;

#define Crystallographic2DAssign(A,B,GAMMA,X) ( \
  (X).a=(A), (X).b=(B), (X).gamma=(GAMMA) )

/* construct crystallographic notation from H[][] */
#define H_to_Crystallographic2D(H,X) ( \
  (X).a = V2LENGTH(H[0]), (X).b = V2LENGTH(H[1]), \
  (X).gamma = RADIAN_TO_DEGREE( acos(V2DOT(H[0],H[1])/(X).a/(X).b) ) )

/* construct H[][] from crystallographic notation */
#define Crystallographic2D_to_H(X,H) ( \
  H[0][0] = (X).a, H[0][1] = 0, \
  H[1][0] = cos( DEGREE_TO_RADIAN((X).gamma) ) * (X).b, \
  H[1][1] = sin( DEGREE_TO_RADIAN((X).gamma) ) * (X).b )

/* determine which of the three vertices of a col-parallelepiped */
/* formed by A[][0], A[][1] is the farthest from origin.         */
double M2maxcolumnradius (double H[2][2]);

/* determine which of the three vertices of a row-parallelepiped */
/* formed by A[0][], A[1][] is the farthest from origin.         */
double M2maxrowradius (double H[2][2]);

/* returns the thickness (>=0) of the parallelepiped formed */
/* by (H[0][], H[1][], H[2][]) in the row i direction.      */
double M2rowthickness (double H[2][2], int i);

/* returns the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[0][], H[1][], H[2][]), in thickness[].       */
double *M2rowthicknesses (double H[2][2], double thickness[2]);

/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[0..1][] is guaranteed */
/* to cover sphere of radius R; returns total number of replicas */
int M2rows_to_cover_sphere (double H[2][2], double R, int nc[2]);
#define total_replica_2D(nc) ( (2*(nc)[0]+1) * (2*(nc)[1]+1) )

/* return the thickness (>=0) of the parallelepiped formed */
/* by (H[][0], H[][1], H[][2]) in the column i direction.  */
double M2columnthickness (double H[2][2], int i);

/* return the three thicknesses (>=0) of the parallelepiped */
/* formed by (H[][0], H[][1], H[][2]), in thickness[].      */
double *M2columnthicknesses (double H[2][2], double thickness[2]);

/* Calculate multiplication factors nc[] whereby -nc[0]:nc[0],   */
/* -nc[1]:nc[1], -nc[2]:nc[2] replica of H[][0..1] is guaranteed */
/* to cover sphere of radius R; return total number of replicas. */
int M2columns_to_cover_sphere (double H[2][2], double R, int nc[2]);

/* dl^2 = dx_i * (1 + 2 * eta_{ij}) * dx_j, where dx meshes the undeformed */
/* because dx = ds * H0, dx' = ds * H, dl^2 = dx' * (dx')^T => above       */
void Lagrangian_strain2D (double H0[2][2], double H[2][2], double eta[2][2]);

/* achieve eta without rotation. M := sqrt(1 + 2*eta), H := H0 * M */
void pure_deform2D (double H0[2][2], double eta[2][2], double H[2][2]);
/* H := H * M */
#define pure_DEFORM2D(H,eta) pure_deform2D(H,eta,H)

/* rotation.c: */

/****************************************/
/* rotation generator & representations */
/****************************************/

/* b[] := right-handed rotation of a[] by theta [RAD] */
double *V2rotate (double a[2], double theta, double b[2]);
#define V2Rotate(theta,b) V2rotate(b,theta,b)

/* theta is now in degrees */
#define V2rotate_in_degrees(a,theta,b) V2rotate(a,DEGREE_TO_RADIAN(theta),b)
#define V2Rotate_in_degrees(theta,b) V2rotate_in_degrees(b,theta,b)

/* right-hand rotation by 90 degrees */
#define V2rotate90(a,b)  ((b)[0]=-(a)[1],(b)[1]=(a)[0])
/* left-hand rotation by 90 degrees */
#define V2rotate270(a,b) ((b)[0]=(a)[1],(b)[1]=-(a)[0])

/* If a rotation makes a[] -> b[], what happens to v[]?    */
/* Assume a[] and b[] are already NORMALIZED; returns v[]. */
double *V2geodesic (double a[2], double b[2], double v[2]);

/* Compute the rotational matrix R[][] corresponding to a geodesic */
/* ("big-circle") rotation that makes a[]->b[] as v[] := v[]*R[][] */
/* Assume a[] and b[] are already NORMALIZED and NOT collinear.    */
void M2geodesic (double a[2], double b[2], double R[2][2]);

/* Generate rotational matrix in z-axis of angle "theta". */
#define M2axialrmat(theta,R) \
  M2ASSIGN(cos(theta),sin(theta),-sin(theta),cos(theta),R)

#define M2RANDOMROTATION(R,theta) ((theta)=Frandom()*2*PI,M2axialrmat(theta,R))

#endif  /* _VecMat2_h */
