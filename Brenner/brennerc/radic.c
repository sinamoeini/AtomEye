/* bren.c - Copyright (c) 1998 Zyvex LLC.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    or its derived works must display the following acknowledgement:
 * 	This product includes software developed by Zyvex LLC.
 * 
 * This software is provided "as is" and any express or implied warranties,
 * including, but not limited to, the implied warranties of merchantability
 * or fitness for any particular purpose are disclaimed. In no event shall
 * Zyvex LLC be liable for any direct, indirect, incidental, special,
 * exemplary, or consequential damages (including, but not limited to,
 * procurement of substitute goods or services; loss of use, data, or
 * profits; or business interruption) however caused and on any theory of
 * liability, whether in contract, strict liability, or tort (including
 * negligence or otherwise) arising in any way out of the use of this
 * software, even if advised of the possibility of such damage.
 */

#include "brenner.h"
#include <math.h>
#include <stdio.h>
#include "myassert.h"
#include "expand.h"

/* The first coordinate of in3 is an index into CLMN or TLMN.
   The second coordinate is either 1, 2, or 3 for L (XINT1), M (XINT2), or N
   (CONJUG).
   The value is the power to raise XINT1, XINT2, or CONJUG to before
   multiplying by the corresponding element of CLMN or TLMN and then adding to
   evaluate the cubic spline. */
static int IN3[64+1][3+1];
/* Index is the inverse of IN3.  
   For II between 1 and 64 inclusive,
   index[IN3[II][1],IN3[II][2], IN3[II][3]] == II. */
static int index [4][4][4];
/* dIN3, strangely enough, has exactly the same values as IN3.  It's used to
   evaluate derivatives, although they could have used IN3 insetad. */
static Float dIN3[64+1][3+1];
/* Cubic spline coefficients for RADIC.  These have to be Double so we can have
 enough precision to get the right answer.*/
static Double CLMN[3+1][10+1][10+1][10+1][64+1];
/* The single precision version.  The difference is that the input to the
   single precision version of the polynomial is the fractional parts
   of XNT1, XNT2, and CONJUG, instead of the entire number.  This will
   hopefully allow us to get accurate answers with single precision arithmetic.
   The "S" at the end is for "Single". */
static Float CLMNS[3+1][10+1][10+1][10+1][64+1];
/* Double precision cubic spline coefficients for TOR. */
static Double TLMN[10+1][10+1][10+1][64+1];
/* Single precision.  Once again, the input to the single precision version is
   the fractional part of the input to the double precision version. */
static Float TLMNS[10+1][10+1][10+1][64+1];
/* in3_initialized is nonzero if we've copied the data files into CLMN and
   TLMN. */
static int in3_initialized = 0;

static void translate_poly (Float dest[10+1][10+1][10+1][64+1],
                            Double src[10+1][10+1][10+1][64+1]) {
  int L, M, N, II;
  for (L=1; L < 11; L++) 
    for (M=1; M < 11; M++) 
      for (N=1; N < 11; N++) {
        Double copy [64+1];
        for (II = 1; II < 65; II++) 
          copy [II] = 0.0;
        for (II = 1; II < 65; II++) {
          /* i, j, and k are the exponents of XNT1, XNT2, and CONJUG for the
             term at position II in array CLMN[L][M][N]. */
          int i, j, k;
          /* ti [z] is the coefficient of F1**z in the expansion of
             (L+F1)**i. Analogously for tj and tk. */
          int ti[4], tj[4], tk[4];
          /* ie will be the term of the expansion of (L+F1)**i that we'll be
             looking at.  Analogously je and ke. */
          int ie, je, ke;
          i = IN3[II][1];
          j = IN3[II][2];
          k = IN3[II][3];
          expand (L, ti, i);
          expand (M, tj, j);
          expand (N, tk, k);
          for (ie = 0; ie < 4; ie++)
            for (je = 0; je < 4; je++)
              for (ke = 0; ke < 4; ke++) {
                copy [index [ie][je][ke]] +=
                  src[L][M][N][II]*ti[ie]*tj[je]*tk[ke];
              }
        }
        for (II = 1; II < 65; II++) 
          dest [L][M][N][II] = copy [II];
      }
}

static void init_in3()
{
  int IC=0;
  int I3D,ITD;
  int i, I, J, K, L, M, N;
  FILE *inter3Dtors;
  for (I=1; I<5; I++) {
    for (J=1; J<5; J++) {
      for (K=1; K<5; K++) {
        IC=IC+1;
        // It's weird that IN3 and dIN3 are the same as each other.
        IN3[IC][1]=I-1;
        IN3[IC][2]=J-1;
        IN3[IC][3]=K-1;
        dIN3[IC][1]=I-1;
        dIN3[IC][2]=J-1;
        dIN3[IC][3]=K-1;
        index[IN3[IC][1]][IN3[IC][2]][IN3[IC][3]]=IC;
      }
    }
  }

    /* Read the tricubic spline coefficients for RADIC (that is, F). */

	for(i = 0; i < 3; ++i)
	{
	  static const char *inter_file[] = {"inter3d_iv_new.d",
					     "inter3d_ch.d",
					     "inter3d_h.d" };
	  FILE *inter3Divnew;
	  if ((inter3Divnew = fopen(inter_file[i],"r"))==NULL)
	  {
	    fprintf(stderr,"opening %s failed", inter_file[i]);
	    my_exit(-1);
	  }
	  fscanf(inter3Divnew, "%d\n",&I3D);
	  while(1) {
		if(fscanf(inter3Divnew, "%d %d %d\n",&L,&M,&N) != 3)
		  break;

		for (I=1; I<=64; ++I) {
		  double val;
		  if(!fscanf(inter3Divnew,"%lf\n", &val))
		    {
		      fprintf(stderr, "parse error in file %s, %d %d %d %d\n",
			      inter_file[i], L, M, N, I);
		      my_exit(-1);
		    }
		  CLMN[i+1][L][M][N][I] = val;
		}
	  }
	  fclose(inter3Divnew);
	}
	
      
	for (L=1; L<11; L++) {
		for (M=1; M<11; M++) {
			for (I=1; I<65; I++) {
				CLMN[1][L][M][10][I]=CLMN[1][L][M][9][I];
				CLMN[2][L][M][10][I]=CLMN[2][L][M][9][I];
				for (N=6; N<11; N++)
					CLMN[3][L][M][N][I]=CLMN[3][L][M][5][I];
			}
		}
	}
    {
      int littlei;
      for (littlei = 1; littlei <= 3; littlei++)
        translate_poly (CLMNS[littlei], CLMN[littlei]);
    }
    /* Read tricubic spline coefficients for torsional potential*/
	if ((inter3Dtors = fopen("inter3dtors.d","r"))==NULL)
		printf("opening inter3dtors failed");
	fscanf(inter3Dtors, "%d\n",&ITD);
	for (J=1; J<109; J++) {
		fscanf(inter3Dtors, "%d %d %d\n",&L,&M,&N);
		for (I=1; I<64; I+=3) {
		  double i0, i1, i2;
		  fscanf(inter3Dtors, "%lf %lf %lf\n", &i0, &i1, &i2);
		  TLMN[L][M][N][I] = i0;
		  TLMN[L][M][N][I+1] = i1;
		  TLMN[L][M][N][I+2] = i2;
		}
		{
		  double val;
		  fscanf(inter3Dtors, "%lf\n",&val);
		  TLMN[L][M][N][64] = val;
		}
	}
	fclose(inter3Dtors);

      
	for (L=1; L<11; L++)
		for (M=1; M<11; M++)
			for (N=4; N<11; N++)
				for (I=1; I<65; I++) {
					TLMN[L][M][N][I]=TLMN[L][M][3][I];
				}
    translate_poly (TLMNS, TLMN);
	in3_initialized = 1;
}

/* Fij(NiT,NjT,Nijconj) */
/* See equation 12 of HCnewpot1.ps. */
/* This version requires double-precision math. */
static inline Float RADICD
(int KI, int KJ, Float XNT1, Float XNT2, Float CONJUG,
 Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr)
{
	int J,L,M,N,KIKJ;
	Double X, RAD;
	Double drdl = 0.0;
	Double drdm = 0.0;
	Double drdn = 0.0;
	const Double *clmn_ptr;
	Double xnt1_pow[4];
	Double xnt2_pow[4];
	Double conj_pow[4];

	if(!in3_initialized)
	  init_in3();
	/*tricubic spline
	Conjug = 1 Noncongugated
	    >1 Conjugated*/

	L=(int)floor(XNT1);
	M=(int)floor(XNT2);
	N=(int)floor(CONJUG);
	RAD = 0;
	KIKJ = KI + KJ - 1;

	if (L>=4) {
		L =4;
		XNT1 = 4;
	}
	if (M>=4){
		M=4;
		XNT2=4;
	}
	/* I observed xnt2 to be 0.99999999999999989, which left m equal to 0
	 *     and then we referenced CLMN out of bounds.  Tim Freeman 15 Feb
	 *     2000.
	 */
	if (M < 1)
	{
	  M = 1;
	  XNT2 = 1;
	}
	/*     Likewise L and xnt1.  Tim Freeman 16 Feb 2000. */
	if (L < 1)
	{
	  L = 1;
	  XNT1 = 1;
	}

	if (N>=9){
		N=9;
		CONJUG=9;
	}

	xnt1_pow[0] = 1.0;
	xnt1_pow[1] = XNT1;
	xnt1_pow[2] = XNT1*XNT1;
	xnt1_pow[3] = XNT1*XNT1*XNT1;
	xnt2_pow[0] = 1.0;
	xnt2_pow[1] = XNT2;
	xnt2_pow[2] = XNT2*XNT2;
	xnt2_pow[3] = XNT2*XNT2*XNT2;
	conj_pow[0] = 1.0;
	conj_pow[1] = CONJUG; 
	conj_pow[2] = CONJUG*CONJUG;
	conj_pow[3] = CONJUG*CONJUG*CONJUG;
	clmn_ptr = &CLMN[KIKJ][L][M][N][0];
	for (J=1; J<65; ++J) {
		X=*++clmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
        // X has been observed to have values around 1.5 million when conj is
        // around 9.  The values from inter3d_iv_new.d can get really large.
        // Thus, double precision seems to really be required here.  We could
        // restructure things so the inputs to the tricubic spline are (XNT1 -
        // L), (XNT2 - M), and (CONJUG - N) instead of XNT1 and XNT2 and
        // CONJUG.  This would reduce the absolute values of the numbers we're
        // adding up and might therefore reduce the precision required enough
        // for things to fit into a single precision floating point number.
        // I'm thinking about this because Pentium 3's have these really fast
        // single precision floating point instructions.  Tim Freeman 30 Jul
        // 2000.
		RAD += X;
		drdl += X*dIN3[J][1]/XNT1;
		drdm += X*dIN3[J][2]/XNT2;
		drdn += X*dIN3[J][3]/CONJUG;
	}
	if(drdl_ptr) *drdl_ptr = drdl;
	if(drdm_ptr) *drdm_ptr = drdm;
	if(drdn_ptr) *drdn_ptr = drdn;
	return RAD;
}

/* Single-precision version of RADIC. */
static inline Float RADICS
(int KI, int KJ, Float XNT1, Float XNT2, Float CONJUG,
 Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr)
{
  int j,L,M,N,KIKJ;
  if(!in3_initialized)
    init_in3();
  /*tricubic spline
    Conjug = 1 Noncongugated
    >1 Conjugated*/
  
  L=(int)floor(XNT1);
  M=(int)floor(XNT2);
  N=(int)floor(CONJUG);
  KIKJ = KI + KJ - 1;
  
  if (L>=4) {
    L =4;
    XNT1 = 4;
  }
  if (M>=4){
    M=4;
    XNT2=4;
  }
  /* I observed xnt2 to be 0.99999999999999989, which left m equal to 0
   *     and then we referenced CLMN out of bounds.  Tim Freeman 15 Feb
   *     2000.
   */
  if (M < 1) {
    M = 1;
    XNT2 = 1;
  }
  /*     Likewise L and xnt1.  Tim Freeman 16 Feb 2000. */
  if (L < 1) {
    L = 1;
    XNT1 = 1;
  }
  if (N>=9){
    N=9;
    CONJUG=9;
  }
  {
    Float RAD = 0;
    Float drdl = 0.0;
    Float drdm = 0.0;
    Float drdn = 0.0;
    const Float *clmns_ptr;
    /* The next three are actually the powers of the *fractional* parts of
       xnt1, xnt2, and conj. */
    Float xnt1_pow[4];
    Float xnt2_pow[4];
    Float conj_pow[4];
    Float xnt1_d[4];
    Float xnt2_d[4];
    Float conj_d[4];
    const Float XNT1F = XNT1 - L;
    const Float XNT2F = XNT2 - M;
    const Float CONJUGF = CONJUG - N;
    clmns_ptr = &CLMNS[KIKJ][L][M][N][0];
    if (XNT1F < 0.000001 &&
        XNT2F < 0.000001 &&
        CONJUGF < 0.000001) {
          RAD=clmns_ptr [index[0][0][0]];
          drdl=clmns_ptr [index[1][0][0]];
          drdm=clmns_ptr [index[0][1][0]];
          drdn=clmns_ptr [index[0][0][1]];
    } else {
      xnt1_pow[0] = xnt2_pow[0] = conj_pow[0] = 1.0;
      xnt1_d[0] = xnt2_d[0] = conj_d[0] = 0;
      for (j = 1; j < 4; j++) {
        xnt1_pow[j] = xnt1_pow[j-1]*XNT1F;
        xnt2_pow[j] = xnt2_pow[j-1]*XNT2F;
        conj_pow[j] = conj_pow[j-1]*CONJUGF;
        xnt1_d[j] = xnt1_pow[j-1]*j;
        xnt2_d[j] = xnt2_pow[j-1]*j;
        conj_d[j] = conj_pow[j-1]*j;
      }
      for (j=1; j<65; ++j) {
        Float c = *++clmns_ptr;
        /* Trust the compiler to do common subexpression elimination well. */
        RAD += c*xnt1_pow[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_pow[IN3[j][3]];
        drdl += c*xnt1_d[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_pow[IN3[j][3]];
        drdm += c*xnt1_pow[IN3[j][1]]*xnt2_d[IN3[j][2]]*conj_pow[IN3[j][3]];
        drdn += c*xnt1_pow[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_d[IN3[j][3]];
      }
    }
    if(drdl_ptr) *drdl_ptr = drdl;
    if(drdm_ptr) *drdm_ptr = drdm;
    if(drdn_ptr) *drdn_ptr = drdn;
    return RAD;
  }
}

#ifdef NDEBUG
#define RADIC_SINGLE
#else
#define RADIC_BOTH
#endif
Float RADIC (int KI, int KJ, Float XNT1, Float XNT2, Float CONJUG,
             Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr)
{
#ifdef RADIC_BOTH
  Float drdl, drdm, drdn;
  Float R1 = RADICS (KI, KJ, XNT1, XNT2, CONJUG, drdl_ptr, drdm_ptr, drdn_ptr);
  Float R2 = RADICD (KI, KJ, XNT1, XNT2, CONJUG, &drdl, &drdm, &drdn);
  assert (fabs (R1 - R2) < 0.001);
  if (drdl_ptr) assert (fabs (*drdl_ptr - drdl) < 0.001);
  if (drdm_ptr) assert (fabs (*drdm_ptr - drdm) < 0.001);
  if (drdn_ptr) assert (fabs (*drdn_ptr - drdn) < 0.001);
  return R1;
#else
   #ifdef RADIC_SINGLE
  return RADICS (KI, KJ, XNT1, XNT2, CONJUG, drdl_ptr, drdm_ptr, drdn_ptr);
   #else
      #ifdef RADIC_DOUBLE  
         return RADICD (KI, KJ, XNT1, XNT2, CONJUG,
                        drdl_ptr, drdm_ptr, drdn_ptr);
      #else
         #error RADIC_DOUBLE, RADIC_SINGLE, or RADIC_BOTH must be set.
      #endif
   #endif
#endif
}

/* Tricubic spline for the torsional rotation.  It's Tij on page 15 of
   HCnewpot1.ps, equation 16. */
static inline Float TORD( Float XNT1, Float XNT2, Float CONJUG,
                          Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr)
{
	int J,L,M,N;
	Double X;
	Double ATOR;
	Double drdl = 0.0;
	Double drdm = 0.0;
	Double drdn = 0.0;

	if(!in3_initialized)
	  init_in3();
	/*tricubic spline for torsional interaction
	CONJUG = 1: NONCONJUGATED
		   > 1: CONJUGATED*/

	ATOR = 0;
	if (XNT1<4 || XNT2<4) {
		const Double *tlmn_ptr;
		Double xnt1_pow[4];
        Double xnt2_pow[4];
		Double conj_pow[4];
		xnt1_pow[0] = 1.0;
		xnt1_pow[1] = XNT1;
		xnt1_pow[2] = XNT1*XNT1;
		xnt1_pow[3] = XNT1*XNT1*XNT1;
		xnt2_pow[0] = 1.0;
		xnt2_pow[1] = XNT2;
		xnt2_pow[2] = XNT2*XNT2;
		xnt2_pow[3] = XNT2*XNT2*XNT2;
		conj_pow[0] = 1.0;
		conj_pow[1] = CONJUG;
		conj_pow[2] = CONJUG*CONJUG;
		conj_pow[3] = CONJUG*CONJUG*CONJUG;
		L = (int)floor(XNT1);
		M = (int)floor(XNT2);
		N = (int)floor(CONJUG);
		tlmn_ptr = &TLMN[L][M][N][0];
		for (J=1; J<65; J++) {
			X = *++tlmn_ptr*xnt1_pow[IN3[J][1]]*xnt2_pow[IN3[J][2]]*conj_pow[IN3[J][3]];
			ATOR += X;
			drdl += X*dIN3[J][1]/XNT1;
			drdm += X*dIN3[J][2]/XNT2;
			drdn += X*dIN3[J][3]/CONJUG;
		}
	}
	if(drdl_ptr) *drdl_ptr = drdl;
	if(drdm_ptr) *drdm_ptr = drdm;
	if(drdn_ptr) *drdn_ptr = drdn;
	return ATOR;
}
 
/* Tricubic spline for the torsional rotation.  It's Tij on page 15 of
   HCnewpot1.ps, equation 16. */
static inline Float TORS ( Float XNT1, Float XNT2, Float CONJUG,
                           Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr)
{
	Float ATOR;
	Float drdl = 0.0;
	Float drdm = 0.0;
	Float drdn = 0.0;

	if(!in3_initialized)
	  init_in3();
	/*tricubic spline for torsional interaction
	CONJUG = 1: NONCONJUGATED
		   > 1: CONJUGATED*/

	ATOR = 0;
	if (XNT1<4 || XNT2<4) {
		const Float *tlmns_ptr;
		Float xnt1_pow[4];
		Float xnt2_pow[4];
		Float conj_pow[4];
        Float xnt1_d[4];
        Float xnt2_d[4];
        Float conj_d[4];
		const int L = (int)floor(XNT1);
		const int M = (int)floor(XNT2);
		const int N = (int)floor(CONJUG);
        const Float XNT1F = XNT1 - L;
        const Float XNT2F = XNT2 - M;
        const Float CONJUGF = CONJUG - N;
        tlmns_ptr = &TLMNS[L][M][N][0];
        if (XNT1F < 0.000001 &&
            XNT2F < 0.000001 &&
            CONJUGF < 0.000001) {
          ATOR=tlmns_ptr [index[0][0][0]];
          drdl=tlmns_ptr [index[1][0][0]];
          drdm=tlmns_ptr [index[0][1][0]];
          drdn=tlmns_ptr [index[0][0][1]];
        } else {
          int j;
          xnt1_pow[0] = xnt2_pow[0] = conj_pow[0] = 1.0;
          xnt1_d[0] = xnt2_d[0] = conj_d[0] = 0;
          for (j = 1; j < 4; j++) {
            xnt1_pow[j] = xnt1_pow[j-1]*XNT1F;
            xnt1_d[j] = xnt1_pow[j-1]*j;
            xnt2_pow[j] = xnt2_pow[j-1]*XNT2F;
            xnt2_d[j] = xnt2_pow[j-1]*j;
            conj_pow[j] = conj_pow[j-1]*CONJUGF;
            conj_d[j] = conj_pow[j-1]*j;
          }
          for (j=1; j<65; j++) {
            const Float c = *++tlmns_ptr;
            ATOR +=
              c*xnt1_pow[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_pow[IN3[j][3]];
            drdl +=
              c*xnt1_d[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_pow[IN3[j][3]];
            drdm +=
              c*xnt1_pow[IN3[j][1]]*xnt2_d[IN3[j][2]]*conj_pow[IN3[j][3]];
            drdn +=
              c*xnt1_pow[IN3[j][1]]*xnt2_pow[IN3[j][2]]*conj_d[IN3[j][3]];
          }
		}
	}
	if(drdl_ptr) *drdl_ptr = drdl;
	if(drdm_ptr) *drdm_ptr = drdm;
	if(drdn_ptr) *drdn_ptr = drdn;
	return ATOR;
}
 
#ifdef NDEBUG
#define TOR_SINGLE
#else
#define TOR_BOTH
#endif

Float TOR ( Float XNT1, Float XNT2, Float CONJUG,
            Float *drdl_ptr, Float *drdm_ptr, Float *drdn_ptr) {
#ifdef TOR_BOTH
  Float drdl, drdm, drdn;
  Float T1 = TORS (XNT1, XNT2, CONJUG, drdl_ptr, drdm_ptr, drdn_ptr);
  Float T2 = TORD (XNT1, XNT2, CONJUG, &drdl, &drdm, &drdn);
  assert (fabs (T1 - T2) < 0.001);
  if (drdl_ptr) assert (fabs (*drdl_ptr - drdl) < 0.001);
  if (drdm_ptr) assert (fabs (*drdm_ptr - drdm) < 0.001);
  if (drdn_ptr) assert (fabs (*drdn_ptr - drdn) < 0.001);
  return T1;
#else
#ifdef TOR_SINGLE
  return TORS (XNT1, XNT2, CONJUG, drdl_ptr, drdm_ptr, drdn_ptr);
#else
  #ifdef TOR_DOUBLE
     return TORD (XNT1, XNT2, CONJUG, drdl_ptr, drdm_ptr, drdn_ptr);
  #else
     #error TOR_DOUBLE, TOR_SINGLE, or TOR_BOTH must be set.
  #endif
#endif
#endif
}
