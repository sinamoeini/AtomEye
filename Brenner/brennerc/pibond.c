/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/
/* Next one defines pibond.c. */
#include "brenc.h"

#ifndef __XALLOC_H__
#include "xalloc.h"
#endif

#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

#ifndef __ATOM_PAIR_INFO_STATE_H__
#include "AtomPairInfoState.h"
#endif

#include "brenner.h"
#include "spgch.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>

#define ncc 10000
/* xqm is the constant 3.7 in equation 10 on page 13. */
#define xqm 3.7
/* att is the constant 3.2 in equation 10 on page 13. */
#define att 3.2
#define pq (M_PI/(xqm - att))
#define ndihed 2
#define PIDT (M_PI/0.3)


/* The next table implements the special cases for G sub C (cos (theta))
   described on pages 11 through 13.  The index to igc is cos (theta) suitably
   quantized, and the output is an index to SPGC, so it says which sixth order
   polynomial spline we will use. */
/* 16*4,2*3,2*2,5*1 */
static const int igc[] = { -1000000, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
			   4, 4, 4, 3, 3, 2, 2, 1, 1, 1, 1, 1 };

/* The next table implements the special cases for G sub H (cos (theta))
   described in the middle of page 17.  The index to igh is cos (theta)
   suitably quantized, and the output is an index to SPGH, so it says which
   sixth order polynomial spline we will use. */
/* 18*3,4*2,3*1 */
static const int igh[] = { -1000000, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
			   3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1 };

static Float xh[3][5][5];
static Float xh1[3][4][4];
static Float xh2[3][4][4];

static Float XDB[4][4][4];
static Float REG[4][4][4];


void init_xh(FILE *inter2D)
{
  /* bicubic spline*/
	int i, j, k;
	const Float RHH = 0.7415886997;
	const Float RCH = 1.09;

	while (1) {
	  /* fscanf expects a double pointer, so xhh must be a double, not a
	     Float. */ 
	  double xhh;
	  fscanf(inter2D, "%d %d %d %lf\n",&i,&j,&k,&xhh);
	  if (i <= 0) 
	    break;
	  xh[i][j][k] = (Float) xhh;
	}

	xh1[2][3][1] = (xh[2][4][1]-xh[2][2][1])/2.0;
	xh1[2][2][2] = (xh[2][3][2]-xh[2][1][2])/2.0;

	xh2[2][2][2] = (xh[2][2][3]-xh[2][2][1])/2.0;
	xh2[2][1][3] = (xh[2][1][4]-xh[2][1][2])/2.0;



	for(i = 1; i < 3; i++)
		for(j = 1; j < 3; j++)
			for (k = 1; k < 3; k++) {
				XDB[i][j][k] = 0.0;
				REG[i][j][k] = 1;
			}


	XDB[2][2][2]=4.0;

	XDB[2][1][2]=4.0;
	XDB[2][2][1]=4.0;

	XDB[2][1][1]=4.0;
	XDB[1][2][1]=0.0;
	XDB[1][2][2]=0.0;

	REG[2][1][2]=exp(XDB[2][1][2]*(RHH-RCH));
	REG[2][2][1]=exp(XDB[2][2][1]*(RCH-RHH));

}

/* k_this is the ktype of atom number index.  See the comments near the word
   "ktype" in brenner.h. */
/* For what it's worth, declaring all these inline measurably improves
   performance, from 5.20 seconds per timing run to 5.18. I didn't say it was
   worth much... */
/* Definition: if something is a per-neighbor input, then it's an array of 250
   somethings, and the first atm_num[index1].num_neighbors elements are inputs,
   and the rest are ignored. 
   
   If something is a per-neighbor output, then it's an array of 250 somethings,
   and the first atm_num[index1].num_neighbors elements are set, and the rest
   are ignored. */
/* This routine does the bulk of the work for p super delta pi sub i j defined
   in equation 7 on page 10.  This is used in equation 2 twice, where i and j
   are swapped the second time.  Not surprisingly, this routine is called
   exactly twice as well.  

   When I refer to "forumla-i" and "formula-j" and "formula-k"
   in the comments below, I mean the symbols in the *formula in the paper*.
   There are variables in this program with similar names but not necessarily
   similar meanings.

   The atom formula-i is in the middle of the angle, and we're dealing with the
   angle between the formula-i formula-j and formula-i formula-k line segments.
*/
static inline void
calc1side(const struct State *info,
          /* k_this is the ktype for atom formula-i.  This could have come from
             atm_num.  I put an assert in below to verify this. */
          const int k_this,
          /* kother is the ktype for atom formula-j. */
          const int kother,
          /* index1 is an index into atm_num, the index of atom formula-i. */
          const int index1,
          /* We call this routine twice.  In the first call, index2 is the
             variable j, which is an index into atom_pair.  In the second call,
             index2 is the variable i, which is an index into atm_num.  index2
             is only used in one place in the code below, to compute a boolean
             value.  It's the atom we want to skip to ensure that (in the
             notation in the paper) k is not j, that is, we don't do physics on
             the angle between a bond and itself. */
          const int index2,
          /* j is the index into atom_pair of the bond we're looking at.
             c_sign tells us which end of that bond we're looking at. */
          const int j,
          /* c_sign is -1 when we call calc1side the first time, 1 when we call
             it the second time.  It tells us the orientation of edge cj, and
             (due to what looks like a mis-fixed bug) the meaning of the index2
             parameter. */
          const int c_sign,
          /* xni points to an array of 2 Float's, either xni or xnj */
          Float *xni,
          /* atom_pairj is the neighbor pairs that j is an index into. */
          const struct AtomPairInfo *const atom_pairj,
          /* atom_pairk is the neighbor pairs for atom index1. */
          const struct AtomPairInfo *const atom_pairk,
          /* cj is the vector between formula-i and formula-j.  It's always a
             copy of atom_pairj[j].cor.  Whether the vector goes from formula-i
             to formula-j or the other way depends on c_sign. */
          const vector cj,
          Float *xsij, Float *xhc[2], Float *cfuni, Float *dcfuni,
          Float* conk, Float* dctjk, Float* dctij, Float* dctik,
          /* ssumk is the middle term under the square root in equation 7 on
             page 10. */
          Float* ssumk, Float* sdalik, Float* xsjk,
          /* xsik is a per-neighbor output. */
          Float* xsik,
          /* sij is the distance between atoms i and j, always equal to
             atom_pairj[j].rcor. */
          const Float sij,
          /* rsqij is the square of sij. */
          const Float rsqij,
          vector *xk /* xl */, Float *cosk, Float *sink,
          /* exnij is an output set to P sub i j, which is reference in
             equation 7 on page 10 and defined in the paper only as a bicubic
             spline. */ 
          Float *exnij,
          /* dexni is a pointer to two Float's.  One is the partial of exnij
             with respect to N super C sub i, and the other is the partial of
             exnij with respect to N super H sub i. */
          Float *dexni)
{
  int k, jj;
  Float qi;
  /* nk will be the number of the neighbor by the time we're done. */
  int nk = 0; /* nl */
  const struct NeighborState *cans = caNeighborsconst (info);
  /* stop_index is one past the index of the last neighbor. */
  const int stop_index = numNeighbors (index1, cans);
  const int *const iNeighbors = neighbors (index1, cans);
  Float gangle, dgdthet;
  Float gangle1, dgdthet1;
  assert (getKtype (index1, info) == k_this);
  assert (sij == atom_pairj[j].rcor);
  /* Apparently we don't get exact equality when we should on the next line. */
  assert (fabs (rsqij - sij * sij) < 0.00001);
  *xsij = 0.0; /* xsji */
  *ssumk = 0.0; /* ssuml */
  *conk = 0.0; /* conl */
  xni[0] = xhc[0][index1]; /* xnj */
  xni[1] = xhc[1][index1];
  /* The next assert says that other is either a hydrogen or carbon. */
  assert(kother-1 < 2);
  xni[kother-1] -= atom_pairj[j].ww;
  /* qi is N super t sub i, in equation 11 on page 13, . */
  qi = xni[0] + xni[1] - 2.0; /* qj */
  *sdalik = 0.0; /* SDALJL */
  for(k = 0; k < stop_index; ++k) /* for(l) */
  {
    /* k is the index of a bond between formula-i and formula-k.  c_sign
       determines which end of the bond that is formula-i; the other end is
       formula-k. */
    Float ali = 0.0; /* alj */
    Float dali = 0.0;  /* dalj */
    Float daldik = 0.0; /* DALDJL */
    Float s3, rr, ss, rsq3, rsq2;
    /* costh will be the cosine of the angle between ij and ik. */
    Float costh;
    int kn, kk;
    kn = iNeighbors[k]; /* ln */
    if((c_sign > 0 ? kn : k) == index2 || atom_pairk[k].lcheck != 1) continue;
    kk = getKtype (kn, info); /* kl */
    ++nk; /* nl */
    /* s3 is the length of the formula-i to formula-k bond. */
    s3 = atom_pairk[k].rcor;
    rsq3 = s3*s3;
    /* xk is the vector from formula-j to formula-k, or perhaps the other way
       around. */
    xk[nk].x = atom_pairk[k].cor.x + c_sign*cj.x;
    xk[nk].y = atom_pairk[k].cor.y + c_sign*cj.y;
    xk[nk].z = atom_pairk[k].cor.z + c_sign*cj.z;
    /* rsq2 is the square of the distance from formula-j to formula-k. */
    rsq2 = xk[nk].x*xk[nk].x + xk[nk].y*xk[nk].y + xk[nk].z*xk[nk].z;
    ss = 2.0*sij*s3;
    rr = rsqij-rsq3;
    /* If you do the analytic geometry, the next line does indeed set costh to
       the cosine of the formula-j formula-i formula-k angle. */
    costh = (rsqij+rsq3-rsq2)/ss;
    if(costh > 1.0) costh = 1.0;
    else if(costh < -1.0) costh = -1.0;
    cosk[nk] = costh;
    sink[nk] = sqrt(1.0-costh*costh); /* sinl */
    // According to the online docs, acos never returns a value greater
    // than pi.  What was the intent here?
    // Sent email to Brenner 29 Jul 2000.  Tim Freeman.
    // if(acos(costh) > M_PI) sink[nk] = -sink[nk];
    assert (acos(costh) <= M_PI);
    if(k_this == 1) {
      /* Atom number index is a carbon. */
      int ig = igc[(int)floor(-costh*12.0)+13];
      /* gangle is G(cos(theta sub j i k)), mentioned in equation 7 on page
         10.  I hope the compiler discovers the shared subexpressions. */
      gangle = SPGC[ig][1] + SPGC[ig][2]*costh + SPGC[ig][3]*costh*costh
        + SPGC[ig][4]*costh*costh*costh
        + SPGC[ig][5]*costh*costh*costh*costh
        + SPGC[ig][6]*costh*costh*costh*costh*costh;
      /* dgdthet is the derivative of gangle with respect to cosine theta. */
      dgdthet = SPGC[ig][2] + SPGC[ig][3]*2*costh
        + SPGC[ig][4]*3*costh*costh
        + SPGC[ig][5]*4*costh*costh*costh
        + SPGC[ig][6]*5*costh*costh*costh*costh;
      if(ig == 4) {
        /* Now we're revising g sub C, per instructions in the middle of page
           13. */
        int ig1;
        /* ali will be Q sub i (N super t sub i) as defined in equation 10 on
           page 13. */
        ali = 0.0;
        /* dali will be the partial of ali with respect to N super t sub i. */
        dali = 0.0;
        if(qi < xqm) {
          /* xqm is 3.7, so if qi >= xqm, we should leave ali and dali at 0. */
          ali = 1.0;
          if(qi > att) {
            Float dtemp = pq*(qi-att);
            ali = (1.0+cos(dtemp))/2.0;
            dali= -pq/2.0*sin(dtemp);
          }
        }
        ig1 = ig+1;
        /* ig was 4, so ig1 is 5.  Presumably this is the row of SPGC that
           defines gamma sub C. */
        gangle1 = SPGC[ig1][1] + SPGC[ig1][2]*costh + SPGC[ig1][3]*costh*costh
          + SPGC[ig1][4]*costh*costh*costh
          + SPGC[ig1][5]*costh*costh*costh*costh
          + SPGC[ig1][6]*costh*costh*costh*costh*costh;
        /* daldik is the partial of g sub C per change in coordination (that
           is, per N super T sub i in the paper, or qi above). */
        daldik = dali*(gangle1-gangle);
        /* Next line is equation 9 on page 13. */
        gangle += ali*(gangle1-gangle);
        /* dgdthet1 is the partial of gamma sub C (cos (theta)) with respect to
           cos (theta). */
        dgdthet1 = SPGC[ig1][2] + SPGC[ig1][3]*2*costh
          + SPGC[ig1][4]*3*costh*costh
          + SPGC[ig1][5]*4*costh*costh*costh
          + SPGC[ig1][6]*5*costh*costh*costh*costh;
        dgdthet += ali*(dgdthet1-dgdthet);
      }
    } else {
      /* k_this != 1 */
      /* Formula-i is not a carbon.  We only get here for hydrogens and
         carbons, so it must be a hydrogen. */
      int ig = igh[(int)floor(-costh*12.0)+13];
      gangle = SPGH[ig][1] + SPGH[ig][2]*costh + SPGH[ig][3]*costh*costh
        + SPGH[ig][4]*costh*costh*costh
        + SPGH[ig][5]*costh*costh*costh*costh
        + SPGH[ig][6]*costh*costh*costh*costh*costh;
      dgdthet = SPGH[ig][2] + SPGH[ig][3]*2*costh
        + SPGH[ig][4]*3*costh*costh
        + SPGH[ig][5]*4*costh*costh*costh
        + SPGH[ig][6]*5*costh*costh*costh*costh;
    }
    {
      Float xtemp, gs, exx, gfx;
      Float dctdjk, dctdij, dctdik;
      /* fc is the cutoff function for the distance between formula-i and
         formula-k. */
      Float fc = atom_pairk[k].ww;
      /* dfc is the derivative of the cutoff with respect to changes in the
         distance between the pair. */
      Float dfc = atom_pairk[k].dww;
      cfuni[nk] = 0.0; /* cfunj */
      dcfuni[nk] = 0.0; /* dcfunj */
      if(kk == 1) { /* kl */
        /* The neighbor of k is a carbon. */
        /* This if statement puts one term the sum in equation 13 on page 14
           into cfuni[nk], and the derivative with respect to the distance
           between the pair into dcfuni[nk]. */
        /* xx is x sub i k, used in equation 13 on page 14 and defined in
           equation 15 on page 14. */
        Float xx = xhc[0][kn]+xhc[1][kn]-fc-2.0;
        if(xx < 3.0) {
          if(xx <= 2.0)
            cfuni[nk] = 1.0;
          else {
            Float px = M_PI*(xx-2.0);
            cfuni[nk] = (1.0+cos(px))/2.0;
            dcfuni[nk] = -fc*sin(px)*M_PI/2.0;
          }
        }
      }
      /* *conk will be the middle sum on the right side of equation 13 on page
         14.  It looks like it would have been slightly better to move this
         inside the if (kk == 1) {...}.*/
      *conk += fc*cfuni[nk];

      /* exx is the exponential term in the middle sum under the square root in
         equation 7 on page 10. */
      if(XDB[k_this][kother][kk] != 0.0)
        /* sij is the distance between atoms formula-i and formula-j.  */
        /* The value inside the exp must be lambda sub i j k, although this is
           undefined in the document and it looks to be a constant in the
           document (and it isn't here). */
        exx = REG[k_this][kother][kk] * exp(XDB[k_this][kother][kk]*(sij-s3));
      else
        exx = 1.0;
      
      /* dctdjk is apparently partial of cosine theta per change in distance
         between formula-j and formula-k, divided by the distance between
         formula-j and formula-k. */
      dctdjk = -2.0/ss; /* DCTDIL */
      /* dctdij is apparently the partial of cosine theta sub j i k per change
         in distance between formula-i and formula-j, divided by the distance
         between formula-i and formula-j. */
      dctdij = (rr+rsq2)/(ss*rsqij); /* DCTDJI */
      /* dctdik is analogous to dctdij; it is the partial of cosine theta sub j
         i k per change in distance between formula-i and formula-k, divided by
         the distance between formula-i and formula-k. */
      dctdik = (-rr+rsq2)/(ss*rsq3); /* DCTDJL */
      dctjk[nk] = dctdjk; /* DCTIL */
      dctij[nk] = dctdij; /* DCTJI */
      dctik[nk] = dctdik; /* DCTJL */
      /* gangle is the G sub i (cos (theta sub j i k)) in the middle term under
         the square root in equation 7 on page 10, which is defined in equation
         9 on page 13. */
      /* exx is the exponential term in the middle sum under the square root in
         equation 7 on page 10. */
      gs = gangle*exx;
      /* ssumk is the middle term under the square root in equation 7 on page
         10. */
      *ssumk += fc*gs;
      xtemp = fc*exx*dgdthet;
      gfx = gs*fc*XDB[k_this][kother][kk];
      *xsij += xtemp*dctdij+gfx/sij;
      xsik[nk] = (gs*dfc-gfx)/s3+xtemp*dctdik; /* XSJL */
      *sdalik += exx*fc*daldik; /* SDALJL */
      xsjk[nk] = xtemp*dctdjk; /* XSIL */
    }
  }
  
  /* check side depends */
  *exnij = 0.0; /* EXNJI */
  dexni[0] = 0.0; /* DEXNJ */
  dexni[1] = 0.0;
  if(k_this == 1) {
    // The code immediately below is trying to test whether xni[1] and xni[0]
    // are approximately integers, so we can avoid doing a spline in that
    // case.  The constants 0.5 in the next two lines of code used to be
    // 1.0e-12, which doesn't work very well when we're doing single
    // precision.  I think 0.5 is right anyway, since in that case floor
    // returns the nearest integer to xni[whatever].  I don't know why it's
    // important to avoid doing a spline when they're integers.  Tim
    // Freeman 29 Jul 2000.
#define ALWAYS_SPLINE
#ifndef ALWAYS_SPLINE
    int nh = (int)floor(xni[1]+0.5);
    int nc = (int)floor(xni[0]+0.5);
    if((fabs((Float)nh-xni[1]) > 1.0e-08)
       || (fabs((Float)nc-xni[0]) > 1.0e-08)) {
#endif
      /* If xni[1] and xni[0] aren't integers, then do BCUINT, which is a
         bicubic spline. */
      *exnij = BCUINT(kother, xni[1], xni[0], &dexni[1], &dexni[0]);
#ifndef ALWAYS_SPLINE
    } else {
      /* xni[1] and xni[0] are approximately integers, so we don't need to
         spline. */
      *exnij = xh[kother][nh][nc];
      dexni[1] = xh1[kother][nh][nc];
      dexni[0] = xh2[kother][nh][nc];
    }
#endif
  }
}

/* This routine implements equation 16 on page 15, computing 
   pi super dh sub i j.  (dh is for dihedral.)
   We update the forces in info, so it can't be const.
*/
static inline void
calc_dihedral_terms(struct State *info, int i, int j, int jn,
                    const struct AtomPairInfoState *apis,
                    /* xnt1 is N super t sub i. */
                    Float xnt1,
                    /* xnt2 is N super t sub j. */
                    Float xnt2,
                    /* conjug is N sub conj. */
                    Float conjug, Float sij,
                    Float *sink, Float *sinl, Float *dctik, Float *dctjk,
                    Float *dctil, Float *dctij, Float *dctjl, Float *dctji,
                    Float *cosk, Float *cosl, vector cj,
                    Float vatt,
                    const vector *xk, const vector *xl,
                    /* btot will be the total energy, that is,
                       pi super dh sub i j. */ 
                    Float *btot,
                    /* dradi is returned as the partial of total energy with
                       respect to xnt1. */
                    Float *dradi,
                    /* dradj is returned as the partial of total energy with
                       respect to xnt2. */
                    Float *dradj,
                    /* drdc is returned as the partial of total energy with
                       respect to conjug. */
                    Float *drdc)
{
  // k will be an index into the list of neighbors of i.
  int k, nk;
  const struct NeighborState *cans = caNeighbors (info);
  int stop_index = numNeighbors (i, cans);
  const int *neighbor_start = neighbors (i, cans);
  const struct AtomPairInfo *iatom_pair = pairsconst (i, apis);
  Float dbtori = 0.0;
  Float dbtorj = 0.0;
  Float dbtorc = 0.0;
  Float btor = 0.0;
  Float datori,datorj ,datorc;
  Float ator = TOR(xnt1,xnt2,conjug,&datori,&datorj,&datorc);
  Float vatt_ator = vatt*ator;

  if(fabs(ator) <= 1.0e-08) return;

  nk = 0;
  for(k = 0; k < stop_index; ++k) {
    Float sink2, rck, fck, dfck;
    vector ck, cl, dt2djl;
    const vector xk_nk = xk[nk+1];
    Float crkx, crky, crkz;
    int nl, l;
    int kn;
    int stop_index2;
    Float cosk_nk;
    const int *jneighbor_start;
    const struct AtomPairInfo *jatom_pair;
    if(k == j || iatom_pair[k].lcheck != 1) continue;
    ++nk;
    /* Used to have a bug here, cosk was getting set to cosk[nk] before we
       incremented nk. */
    cosk_nk = cosk[nk];
    if(fabs(sink[nk]) < 1.0e-01) continue;
    sink2 = sink[nk]*sink[nk];
    kn = neighbor_start[k];
    ck.x = iatom_pair[k].cor.x;
    ck.y = iatom_pair[k].cor.y;
    ck.z = iatom_pair[k].cor.z;
    crkx = ck.y*cj.z-cj.y*ck.z;
    crky = ck.z*cj.x-cj.z*ck.x;
    crkz = ck.x*cj.y-cj.x*ck.y;
    dt2djl.x = -cj.y*crkz+cj.z*crky;
    dt2djl.y = -cj.z*crkx+cj.x*crkz;
    dt2djl.z = -cj.x*crky+cj.y*crkx;

    rck = iatom_pair[k].rcor;

    if(getKtype (kn, info) == 2) {
      /* It's a hydrogen.  This is computing f super c sub i j from equation 18
       on page 16, but I don't know why we aren't using the table so
       laboriously precomputed for this purpose in mtable.c.  Tim Freeman  5
       Aug 2000. */
      fck = 1.0;
      dfck = 0.0;
      if(rck >= 1.60) continue;
      if(rck >= 1.30) {
        Float dtemp= PIDT*(rck-1.30);
        fck = (1.0+cos(dtemp))/2.0;
        dfck = -PIDT/2.0*sin(dtemp);
      }
    } else {
      fck = iatom_pair[k].ww;
      dfck = iatom_pair[k].dww;
    }
    jneighbor_start = neighbors (jn, cans);
    stop_index2 = numNeighbors (jn, cans);
    jatom_pair = pairsconst (jn, apis);
    nl = 0;
    for(l = 0; l < stop_index2; ++l) {
      int ln = jneighbor_start[l];
      Float sinl2, rcl, t1, t2, cw, bt, fcl, dfcl;
      Float aa, aaa1, at2, rp1, rp2, rp3, rp4, rp5, rep;
      Float dt1dik, dt1djk, dt1djl, dt1dil, dt1dij;
      Float crlx, crly, crlz;
      vector dt2dik, dt2dij;
      if(ln == i || jatom_pair[l].lcheck != 1) continue;
      ++nl;
      if(fabs(sinl[nl]) < 1.0e-01) continue;
      sinl2 = sinl[nl]*sinl[nl];
      cl.x = jatom_pair[l].cor.x;
      cl.y = jatom_pair[l].cor.y;
      cl.z = jatom_pair[l].cor.z;
      rcl = jatom_pair[l].rcor;

      if(getKtype (ln, info) == 2) {
        fcl = 1.0;
        dfcl = 0.0;
        // FIXME Another trig computation that is probably in some lookup
        // table somewhere.  Tim Freeman 11 Sep 2000
        if(rcl >= 1.60) continue;
        if(rcl >= 1.30) {
          Float dtemp = PIDT*(rcl-1.30);
          fcl = (1.0+cos(dtemp))/2.0;
          dfcl = -PIDT/2.0*sin(dtemp);
        }
      } else {
        fcl = jatom_pair[l].ww;
        dfcl = jatom_pair[l].dww;
      }
      t1 = rck*rcl*sij*sij*sink[nk]*sinl[nl];
      
      dt1dik = 1.0/rck/rck - dctik[nk]/sink2*cosk_nk;
      
      dt1djk = -dctjk[nk]/sink2*cosk_nk;
      
      dt1djl = 1.0/rcl/rcl-dctjl[nl]/sinl2*cosl[nl];
      
      dt1dil = -dctil[nl]/sinl2*cosl[nl];

      dt1dij = 2.0/sij/sij-dctij[nk]/sink2*cosk_nk-dctji[nl]/sinl2*cosl[nl];

      crlx = cj.y*cl.z-cl.y*cj.z;
      crly = cj.z*cl.x-cl.z*cj.x;
      crlz = cj.x*cl.y-cl.x*cj.y;

      t2 = crkx*crlx+crky*crly+crkz*crlz;

      cw = t2/t1;
      bt = (1.0-cw*cw);
      btor += bt*fck*fcl;

      dt2dik.x = -cj.z*crly+cj.y*crlz;
      dt2dik.y = -cj.x*crlz+cj.z*crlx;
      dt2dik.z = -cj.y*crlx+cj.x*crly;

      dt2dij.x = ck.z*crly-cl.z*crky-ck.y*crlz+cl.y*crkz;
      dt2dij.y = ck.x*crlz-cl.x*crkz-ck.z*crlx+cl.z*crkx;
      dt2dij.z = ck.y*crlx-cl.y*crkx-ck.x*crly+cl.x*crky;

      aaa1 = vatt_ator*bt;
      aa = -vatt_ator*2.0*cw/t1*fcl*fck;
      at2 = aa*t2;

      rp1 = -dt1dij*at2;
      rp2 = -dt1dik*at2+aaa1*fcl*dfck/rck;
      rp3 = -dt1djl*at2+aaa1*fck*dfcl/rcl;
      rp4 = -dt1djk*at2;
      rp5 = -dt1dil*at2;

      transForcex (i, jn, info, rp1*cj.x + aa*dt2dij.x);
      transForcex (i, kn, info, rp2*ck.x + aa*dt2dik.x);
      transForcex (jn, ln, info, rp3*cl.x + aa*dt2djl.x);
      transForcex (jn, kn, info, rp4*xk_nk.x);
      transForcex (i, ln, info, rp5*xl[nl].x);

      transForcey (i, jn, info, rp1*cj.y + aa*dt2dij.y);
      transForcey (i, kn, info, rp2*ck.y + aa*dt2dik.y);
      transForcey (jn, ln, info, rp3*cl.y + aa*dt2djl.y);
      transForcey (jn, kn, info, rp4*xk_nk.y);
      transForcey (i, ln, info, rp5*xl[nl].y);

      transForcez (i, jn, info, rp1*cj.z + aa*dt2dij.z);
      transForcez (i, kn, info, rp2*ck.z + aa*dt2dik.z);
      transForcez (jn, ln, info, rp3*cl.z + aa*dt2djl.z);
      transForcez (jn, kn, info, rp4*xk_nk.z);
      transForcez (i,ln, info, rp5*xl[nl].z);
    }
  }

  *btot += btor*ator;
  *dradi += datori*btor;
  *dradj += datorj*btor;
  *drdc += datorc*btor;
}

static inline void
calc_many_body(struct State *info, const struct AtomPairInfoState *apis,
               /* i_side is really a boolean, true the first time and false the
                  second time. */
               int i_side, int index1, int index2,
               /* What does indexj mean now that things are reorganized?  It
                used to be the index of the neighbor pair in the entire
                neighbor array the second time, and -1 the first time.  We
                ignore it the first time, because i_side is true.  The second
                time, we use it to skip repeated examination of something we
                already looked at if k is equal to indexj.

                After the reorganization, it is a position in the list of
                neighbors of index1.
               */
               int indexj,
               Float vdbdi,
               /* xsik is a per-neighbor input. */
               Float *xsik,
               Float *dexni,
               Float vdrdi, Float vdrdc, Float *cfuni, Float sdalik,
               Float *xsjk, vector *xk, Float *dcfuni, Float conk)
{
  // The code that increments nk here has to match up with the code that
  // increments nk in calc1side, since we use nk to index into cfuni which is
  // computed in calc1side.
  int nk = 0;
  int n; /* nl */
  int k; /* l */
  int m;
  Float dwr;
  const struct NeighborState *cans = caNeighbors (info);
  const int stop_index = numNeighbors (index1, cans);
  const int *const kNeighbors = neighbors (index1, cans);
  const struct AtomPairInfo *const kPairs = pairsconst (index1, apis);
  for(k = 0; k < stop_index; ++k) {
    Float rp, rp1, rp2, rep, ddr;
    int kn = kNeighbors [k];
    int kk = getKtype (kn, info);
    if(i_side && k == indexj) continue;
    if(!i_side && kn == index2) continue;
    if(kPairs[k].lcheck != 1) continue;
    dwr = kPairs[k].dww/kPairs[k].rcor;
    ++nk; // nk only increases when we get past the continue's above.
    /*
     * First Neighbors
     */
    rp1 = vdbdi*(xsik[nk]+dwr*dexni[kk-1]) + dwr*(vdrdi+vdrdc*cfuni[nk])
      + vdbdi*dwr*sdalik;
    rp2 = vdbdi*xsjk[nk];
    transForcex (index1, kn, info, rp1*kPairs[k].cor.x);
    transForcey (index1, kn, info, rp1*kPairs[k].cor.y);
    transForcez (index1, kn, info, rp1*kPairs[k].cor.z);
    /*
     * Angular Forces
     */
    transForcex (index2, kn, info, rp2*xk[nk].x);
    transForcey (index2, kn, info, rp2*xk[nk].y);
    transForcez (index2, kn, info, rp2*xk[nk].z);
    /*
     * Second Neighbors via RADIC
     */
    ddr = vdrdc*dcfuni[nk]*2.0*conk;
    if(ddr == 0.0) continue;
    {
      int stop_index2 = numNeighbors (kn, cans);
      const int *const mNeighbors = neighbors (kn, cans);
      const struct AtomPairInfo *const mPairs = pairsconst (kn, apis);
      for(m = 0; m < stop_index2; ++m) {
        int mn;
        if(mPairs[m].lcheck != 1) continue;
        mn = mNeighbors[m];
        if(mn == kn) continue;
        rp = ddr*mPairs[m].dww/mPairs[m].rcor;
        transForcex (kn, mn, info, rp*mPairs[m].cor.x);
        transForcey (kn, mn, info, rp*mPairs[m].cor.y);
        transForcez (kn, mn, info, rp*mPairs[m].cor.z);
      }
    }
  }
}

/* The return value is the total energy due to the pi bonds. */
/* This routine should be "pure", in the sense that it only accesses the
   physical state of the system by using the api defined in api.h.
*/
/* We update the forces in info, so it can't be const. */
Float
pibond(struct State *const info, struct AtomPairInfoState *const apis)
{
  const int num_atms = pastLastIndex (info);
  /* When we're looking at a bond from atom i to atom j, and we're considering
     all neighbors of j, then xk will be the vector from atom i to each
     neighbor of j. */
  vector xk[250+1];
  /* Similarly, when we're considering all neighbors of i, xl will be the
     vector from atom j to each neighbor of i. */
  vector xl[250+1],cj,ck,cl,rk,rl,dt2dik,dt2djl,dt2dij;
  Float xsik[250+1],xsjk[250+1],xsil[250+1],xsjl[250+1]
    ,dexni[2],dexnj[2],xni[2],xnj[2]
    ,cfuni[ncc],cfunj[ncc],dcfuni[ncc],dcfunj[ncc]
    ,cosk[250+1],cosl[250+1],sink[250+1],sinl[250+1]
    ,dctjk[250+1],dctij[250+1],dctik[250+1],dctil[250+1],dctji[250+1],
    dctjl[250+1];
  Float tote = 0;
  int i;
  Float *hydrogens_connected = (Float*)xmalloc(num_atms*sizeof(Float));
  Float *carbons_connected = (Float*)xmalloc(num_atms*sizeof(Float));
  Float *xhc[2];
  const struct NeighborState *const cans = caNeighborsconst (info);
  xhc[0] = hydrogens_connected;
  xhc[1] = carbons_connected;
  for(i = 0; i < num_atms; ++i)
    xhc[0][i] = xhc[1][i] = 0;
  for(i = 0; i < ncc; ++i)
    cfuni[i] = cfunj[i] = dcfuni[i] = dcfunj[i] = 0;
  dexni[0]=dexnj[0]=xni[0]=xnj[0] = 0;
  dexni[1]=dexnj[1]=xni[1]=xnj[1] = 0;
  /*
   * Find the number of hydrogens and carbons connected to each atom.
   */
  for(i = firstIndex (info); i < num_atms; ++i) {
    const int *const iNeighbors = neighbors (i, cans);
    const struct AtomPairInfo *const pairs = pairsconst (i, apis);
    int nCount = numNeighbors (i, cans);
    int jn;
    /* Why do we start counting at 1 instead of 0 here?  We might not have any
       hydrogens at all in the figure, and in that case we'll set some field of
       hydrogens_connected to 1.  Weird. Tim Freeman 6 Sep 2000. */
    hydrogens_connected[i] = 1;
    carbons_connected[i] = 1;
    for (jn = 0; jn < nCount; jn++) {
      if(1 == pairs[jn].lcheck) {
        const int j = iNeighbors [jn];
        const int atype = getKtype (j, info);
        assert(atype > 0 && atype <= 2);
        xhc[atype-1][i] += pairs[jn].ww;
      }
    }
  }
  /*
   * Sum over bonds between atoms I and J
   */
  for(i = 0; i < num_atms; ++i) {
    /* j is the index of this neighbor pair in the neighbors array.  Note that
       j is an entirely different beast from i, which is an index into the atom
       array. */
    int j;
    const int ki = getKtype (i, info);
    const int * const iNeighbors = neighbors (i, cans);
    const int nCount = numNeighbors (i, cans);
    const struct AtomPairInfo *const atom_pair = pairsconst (i, apis);
    for (j = 0; j < nCount; j++) {
      /* jn will be the atom number of the current neighbor. */
      const int jn = iNeighbors [j];
      /* calc_dihedral_terms adds a value to btot but does not otherwise read
         its value. */
      Float conjug, xnt1, xnt2, btot;
      /* vdbdi, vdbdj, vdrdi, vdrdj, and vdrdc are variables that are
         consequences of what we got from calc_dihedral_terms and then inputs
         into calc_many_body. */
      Float vdbdi, vdbdj, vdrdi, vdrdj, vdrdc;
      /* sdalik and sdaljl are outputs from calc1side I don't understand
         yet and inputs to calc_many_body. */
      Float sdalik, sdaljl;
      /* conk and conl are outputs from calc1side I don't understand yet and
         inputs into calc_many_body. */
      Float conk, conl;
      /* cj is the vector from i to j. */
      vector cj;
      /* If one of them isn't hydrogen or carbon, then look at the 
         next pair. */
      if(atom_pair[j].lcheck != 1) continue;
      if(i >= jn) continue;
      /* Now we're considering each pair only once. */
      cj.x = atom_pair[j].cor.x;
      cj.y = atom_pair[j].cor.y;
      cj.z = atom_pair[j].cor.z;
      {
        /* kj will be the ktype for the current neighbor. */
        const int kj = getKtype (jn, info);
        /* sij will be the distance from one atom to the other. */
        const Float sij = atom_pair[j].rcor;
        const Float rsqij = sij*sij;
        /* xsij and xsji are an output from calc1side I don't understand
           yet. */
        Float xsij, xsji;
        /* ssumk is the middle term under the square root in equation 7 on page
           10. */
        Float ssumk, ssuml;
        /* exnij is pi super rc sub i j, mentioned in equation 3 on page 6 and
           defined in equation 12 on page 14.  exnji is pi super rc sub j i. */
        Float exnij, exnji;
        calc1side(info, ki, kj, i,  j, j, -1,
                  xni, atom_pair, atom_pair,
                  cj, &xsij, xhc, cfuni, dcfuni,
                  &conk, dctjk, dctij, dctik, &ssumk, &sdalik, xsjk, xsik,
                  sij, rsqij,
                  xk, cosk, sink, &exnij, dexni);
        calc1side(info, kj, ki, jn, i, j, 1,
                  xnj, atom_pair,
                  pairs (jn, apis),
                  cj, &xsji, xhc, cfunj, dcfunj,
                  &conl, dctil, dctji, dctjl, &ssuml, &sdaljl, xsil, xsjl,
                  sij, rsqij,
                  xl, cosl, sinl, &exnji, dexnj);
        {
          Float dbdzi, dbdzj;
          Float dradi = 0.0;
          Float dradj = 0.0;
          Float drdc = 0.0;
          {
            Float bji, bij, rad;
            {
              Float dij;
              /* dij is the expression inside the square root of equation 7 on page
                 10. */
              dij = (1.0+exnij+ssumk);
              /* The next line is equation 7 on page 10.  
                 bij is p super sigma pi sub i j. */
              bij = 1/sqrt(dij);
              dbdzi = -0.50*bij/dij;
            } /* Forget dij. */
            {
              Float dji = (1.0+exnji+ssuml);
              bji = 1/sqrt(dji);
              dbdzj = -0.50*bji/dji;
            } /* Forget dji. */
            /* Next line is equation 13 on page 14.
               conjug is N super conj sub i j. */
            conjug = 1.0 + (conk*conk) + (conl*conl);
            /* Next line is equation 8a on page 11.
               xnt1 is N super C sub i. */
            xnt1 = xni[0]+xni[1]-1.0;
            /* Next line is equation 8b on page 11.
               xnt2 is N super H sub i. */
            xnt2 = xnj[0]+xnj[1]-1.0;
            /* Next line is equation 12 on page 14.
               rad is pi super rc sub i j. 
               rc probably stands for radical. */
            rad = RADIC(ki,kj,xnt1,xnt2,conjug,&dradi,&dradj,&drdc);
            btot = (bji+bij+rad);
          } /* Forget bji, bij, rad. */
          {
            /* vatt is the energy due to attraction between the pair. */
            const Float vatt = atom_pair[j].exx1;
            const int kikj = ki+kj;
            /*
             * Dihedral terms
             */
            /* ndihed is 2, and kikj is the sum of the ktypes ki and kj, so
               "kikj == ndihed" is a very obscure way to say that both i and j
               are carbon (ktype 1), since ktype's are integers and are never
               zero or negative. */
            if(kikj == ndihed)
              calc_dihedral_terms(info, i, j, jn, apis, xnt1, xnt2,
                                  conjug, sij, sink, sinl, dctik, dctjk,
                                  dctil, dctij, dctjl, dctji, cosk, cosl, cj,
                                  vatt, xk, xl, &btot, &dradi, &dradj, &drdc);
            /* btot is the sum of the b's, specifically the b sub i j in
               equation 1 on page 4.  vatt is V super A sub i j in equation 1
               on page 4. */
            tote -= btot*vatt;
            /* dbdzi is the derivative of b with respect to i-dont-know-what.
               vdbdi is the derivative of the total energy with respect to the
               same thing, due to changes in b. */
            vdbdi = vatt*dbdzi;
            vdbdj = vatt*dbdzj;
            vdrdc = vatt*drdc;
            vdrdi = vatt*dradi;
            vdrdj = vatt*dradj;
          } /* Forget vatt. */
        } /* Forget dbdzi, dbdzj, drdc, dradi, dradj */
        {
          Float rp = vdbdi*xsij + vdbdj*xsji + btot*atom_pair[j].dexx1;
#ifndef NDEBUG
          if(!(rp < 1.e99))
            printf("rp %d %d %d %11.8f %11.8f %11.8f %11.8f %11.8f %11.8f\n",
                   j, i, jn, vdbdi,xsij, vdbdj,xsji, btot, atom_pair[j].dexx1);
          
#endif
          transForcex (i, jn, info, rp*cj.x);
          transForcey (i, jn, info, rp*cj.y);
          transForcez (i, jn, info, rp*cj.z);
        } /* Forget rp. */
      } /* Forget sij, rsqij, xsij, xsji. */
      /*
       * Add many-body forces
       */
      calc_many_body(info, apis, 1, i, jn, j, /* i side */
                     vdbdi, xsik, dexni, vdrdi, vdrdc, cfuni, sdalik, xsjk, xk,
                     dcfuni, conk);
      calc_many_body(info, apis, 0, jn, i, -1, /* j side */
                     vdbdj, xsjl, dexnj, vdrdj, vdrdc, cfunj, sdaljl, xsil, xl,
                     dcfunj, conl);
    }
  }
  free(hydrogens_connected);
  free(carbons_connected);
  return tote;
}
