/* bcuint.c - Copyright (c) 1998 Zyvex LLC.
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
#include "expand.h"
#include "myassert.h"

static const char* ktype_name[4+1] = {NULL, "Carbon", "Hydrogen", "Silicon", "Germanium"};

/* Bicubic spline */

/* Hij(NiH,NiC)*/

#define MAX_BC_NEIGHBORS 32

static int IN2[16+1][2+1];
/* Index is the inverse of IN2.
   For II between 1 and 16 inclusive,
   index[IN3[II][1],IN3[II][2]] == II. */
static int index [4][4];
static Double CLM[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];
static Float CLMS[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];

void
init_in2(FILE *inter2D)
{
	int IC=0;
	int I, J, K, L, M;
	for (I=1; I<5; I++) {
		for (J=1; J<5; J++) {
			IC=IC+1;
			IN2[IC][1]=I-1;
			IN2[IC][2]=J-1;
            index[IN2[IC][1]][IN2[IC][2]]=IC;
		}
	}

    /* Zero bicubic spline coefficients*/
	for (I=1; I<3; I++)
		for (L=1; L <= MAX_BC_NEIGHBORS; L++)
			for (M=1; M <= MAX_BC_NEIGHBORS; M++)
				for (J=1; J<17; J++)
					CLM[I][L][M][J]=0.0;
	for (K=1; K<73; K++) {
		fscanf(inter2D, "%d %d %d\n",&I,&L,&M);
		for (J=1; J<16; J+=4) {
		  double j0, j1, j2, j3;
		  fscanf(inter2D, "%lf %lf %lf %lf\n", &j0, &j1, &j2, &j3);
		  CLM[I][L][M][J] = (Float) j0;
		  CLM[I][L][M][J+1] = (Float) j1;
		  CLM[I][L][M][J+2] = (Float) j2;
		  CLM[I][L][M][J+3] = (Float) j3;
		}
	}
    {
      int littlei;
      for (littlei = 1; littlei < 3; littlei++) 
        for (L=1; L <= MAX_BC_NEIGHBORS; L++) 
          for (M=1; M <= MAX_BC_NEIGHBORS; M++) {
            int II;
            Double copy [16+1];
            for (II = 1; II < 17; II++) 
              copy [II] = 0.0;
            for (II = 1; II < 17; II++) {
              /* i and j are the exponents of XX1, XX2 for the
                 term at position II in array CLM[L][M]. */
              int i, j;
              /* ti [z] is the coefficient of F1**z in the expansion of
                 (L+F1)**i. Analogously for tj. */
              int ti[4], tj[4];
              /* ie will be the term of the expansion of (L+F1)**i that we'll
                 be looking at.  Analogously je. */
              int ie, je, ke;
              i = IN2[II][1];
              j = IN2[II][2];
              expand (L, ti, i);
              expand (M, tj, j);
              for (ie = 0; ie < 4; ie++)
                for (je = 0; je < 4; je++) {
                  copy [index [ie][je]] +=
                    CLM[littlei][L][M][II]*ti[ie]*tj[je];
                }
            }
            for (II = 1; II < 65; II++) 
              CLMS [littlei][L][M][II] = copy [II];
          }
    }
}

inline static Float
BCUINTD(int KJ, Float XX1, Float XX2, Float *ansy1_ptr, Float *ansy2_ptr)
{
	int J, NH, NC;
	Double X, ANSY;
	Double ansy1 = 0;
	Double ansy2 = 0;
	Double xx1_pow[4];
	Double xx2_pow[4];
	const Double *clm_ptr;

	/*bicubic spline*/

	/* FIXME Having 1.0e-12 here stops working when we're 
	 * single-precision.  The Fortran code didn't have this.  Why does the C
	 * code need it?  Rip it out until I understand why.  Tim Freeman 30 Jul
	 * 2000. */
#if 0
	if (XX1+1.0e-12 < 1)  /*can't have zero indices*/
	  XX1++; /*treat atom as bonded to itself*/
	if (XX2+1.0e-12 < 1)
	  XX2++;
#endif

    if (XX1 < 1) XX1 = 1;
    if (XX2 < 1) XX2 = 1;
    /* Took out adjustment adding 1.0e-8 to XX1 and XX2 before taking the
     * floor.  We still get comparable results to the Fortran version, so I
     * suppose it's okay.  Tim Freeman 30 Jul 2000.
     */
	NH=(int)floor(XX1);
	NC=(int)floor(XX2);

#if 0
	/* I had added the 1.0e-8 adjustment to handle suspicious behavior
	 * I saw using http://www.rahul.net/pcm/mmtk/PDB/crambin_diamond.pdb.
	 * Closer investigation with the single-precision code shows that
	 * while adding 1.0e-6 might be a good idea, the differences I have
	 * measured in the resulting atom positions are about 1 part in 1.0e6,
	 * which doesn't appear to matter. Enable this block of code to
	 * print the places where the differences show up. - pcm 2000-08-11 */
	{
	  int NH2=(int)floor(XX1+1.0e-6);
	  int NC2=(int)floor(XX2+1.0e-6);
	  if(NC2 != NC || NH2 != NH)
	    printf("NC2 %d NC %d NH2 %d NH %d XX1 %15.12f XX2 %15.12f\n",
		   NC2,NC,NH2,NH,XX1,XX2);
	}
#endif

    assert (0 != NH);
    assert (0 != NC);
    if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS)
    {
      fprintf(stderr,
	      "A %s atom has too many neighbors for the bicubic spline.\n"
	      "It has %d Hydrogen and %d Carbon neighbors; the maximum is %d (MAX_BC_NEIGHBORS).\n",
	      KJ >= 1 && KJ <= 4 ? ktype_name[KJ] : "invalid type",
	      NH, NC, MAX_BC_NEIGHBORS);
      my_exit(-1);
    }
	if (KJ==0) {
		fprintf(stderr,
                "error BCUINT unexpected zero: "
                "KJ %d NH %d NC %d XX1 %f XX2 %f\n",
                KJ, NH, NC, XX1, XX2);
		my_exit(-1);
	}

	xx1_pow[0] = 1.0;
	xx1_pow[1] = XX1;
	xx1_pow[2] = XX1*XX1;
	xx1_pow[3] = XX1*XX1*XX1;
	xx2_pow[0] = 1.0;
	xx2_pow[1] = XX2;
	xx2_pow[2] = XX2*XX2;
	xx2_pow[3] = XX2*XX2*XX2;
	clm_ptr = &CLM[KJ][NH][NC][0];
	ANSY=0;
	for (J=1; J<17; J++) {
		X = *++clm_ptr*xx1_pow[IN2[J][1]]*xx2_pow[IN2[J][2]];
		ANSY += X;
		ansy1 += X*IN2[J][1]/XX1;
		ansy2 += X*IN2[J][2]/XX2;
	}
	if(ansy1_ptr) *ansy1_ptr = ansy1;
	if(ansy2_ptr) *ansy2_ptr = ansy2;
	return ANSY;
}

inline static Float
BCUINTS(int KJ, Float XX1, Float XX2, Float *ansy1_ptr, Float *ansy2_ptr)
{
  int NH, NC;
  /*bicubic spline*/
  
  if (XX1 < 1) XX1 = 1;
  if (XX2 < 1) XX2 = 1;
  NH=(int)floor(XX1);
  NC=(int)floor(XX2);
  
  assert (0 != NH);
  assert (0 != NC);
  if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS)
  {
    fprintf(stderr,
	    "A %s atom has too many neighbors for the bicubic spline.\n"
	    "It has %d Hydrogen and %d Carbon neighbors; the maximum is %d (MAX_BC_NEIGHBORS).\n",
	    KJ >= 1 && KJ <= 4 ? ktype_name[KJ] : "invalid type",
	    NH, NC, MAX_BC_NEIGHBORS);
    my_exit(-1);
  }
  if (KJ==0) {
    fprintf(stderr,
            "error BCUINT unexpected zero: "
            "KJ %d NH %d NC %d XX1 %f XX2 %f\n",
            KJ, NH, NC, XX1, XX2);
    my_exit(-1);
  }
  {
    Float ANSY;
    Float ansy1 = 0;
    Float ansy2 = 0;
    const Float *clms_ptr = &CLMS[KJ][NH][NC][0];
    Float XX1F = XX1 - NH;
    Float XX2F = XX2 - NC;
    if (XX1F < 0.000001 && XX2F < 0.000001) {
      ANSY=clms_ptr [index[0][0]];
      ansy1=clms_ptr [index[1][0]];
      ansy2=clms_ptr [index[0][1]];
    } else {
      int j;
      Float xx1_pow[4];
      Float xx2_pow[4];
      Float xx1_d[4];
      Float xx2_d[4];
      xx1_pow[0] = xx2_pow[0] = 1;
      xx1_d[0] = xx2_d[0] = 0;
      for (j = 1; j < 4; j++) {
        xx1_pow[j] = xx1_pow[j-1]*XX1F;
        xx1_d[j] = xx1_pow[j-1]*j;
        xx2_pow[j] = xx2_pow[j-1]*XX2F;
        xx2_d[j] = xx2_pow[j-1]*j;
      }
      ANSY=0;
      for (j=1; j<17; j++) {
        const Float c = *++clms_ptr;
        ANSY += c*xx1_pow[IN2[j][1]]*xx2_pow[IN2[j][2]];
        ansy1 += c*xx1_d[IN2[j][1]]*xx2_pow[IN2[j][2]];
        ansy2 += c*xx1_pow[IN2[j][1]]*xx2_d[IN2[j][2]];
      }
    }
    if(ansy1_ptr) *ansy1_ptr = ansy1;
    if(ansy2_ptr) *ansy2_ptr = ansy2;
    return ANSY;
  }
}
  
#ifdef NDEBUG
#define BCUINT_SINGLE
#else
#define BCUINT_BOTH
#endif

/* Compute P sub i j, which is referenced in equation 7 on page 10.  
   KJ is probably the identity of the "j" atom in P sub i j.  I don't know why
   we don't have a KI argument.
   Judging by the variable names, XNT1 is probably N super H sub i and XNT2 is
   N super C sub i.  Note that this differs from the order of arguments of P
   sub i j in the paper.  
   ansy1_ptr is an output, the derivative of P sub i j with respect to
   XNT1.  
   ansy2_ptr is an output, the derivative of P sub i j with respect to
   XNT2.  */
Float BCUINT (int KJ, Float XNT1, Float XNT2,
              Float *ansy1_ptr, Float *ansy2_ptr) {
#ifdef BCUINT_BOTH
  Float ansy1, ansy2;
  Float T1 = BCUINTS (KJ, XNT1, XNT2, ansy1_ptr, ansy2_ptr);
  Float T2 = BCUINTD (KJ, XNT1, XNT2, &ansy1, &ansy2);
  assert (fabs (T1 - T2) < 0.001);
  if (ansy1_ptr) assert (fabs (*ansy1_ptr - ansy1) < 0.001);
  if (ansy2_ptr) assert (fabs (*ansy2_ptr - ansy2) < 0.001);
  return T1;
#else
#ifdef BCUINT_SINGLE
  return BCUINTS (KJ, XNT1, XNT2, ansy1_ptr, ansy2_ptr);
#else
  #ifdef BCUINT_DOUBLE
     return BCUINTD (KJ, XNT1, XNT2, ansy1_ptr, ansy2_ptr);
  #else
     #error BCUINT_DOUBLE, BCUINT_SINGLE, or BCUINT_BOTH must be set.
  #endif
#endif
#endif
}
