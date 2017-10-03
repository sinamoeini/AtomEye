/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/

/* The code in this file should be "pure", in the sense that it should only
   access the physical state of the system through the API defined in api.h and
   neighborstate.h.  Then, when I want to make the code access Fungimol data
   structures, all I have to do is change the next two includes.  Tim Freeman
   5 Sep 2000. */
#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

#ifndef __ATOM_PAIR_INFO_STATE_H__
#include "AtomPairInfoState.h"
#endif

#ifndef __RMAX_H__
#include "rmax.h"
#endif

#ifndef __myassert_h__
#include "myassert.h"
#endif

#ifndef __CALC_DIST_H__
#include "calc_dist.h"
#endif

#include "brenner.h"
#include "sili_germ.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "spgch.h"

#undef min
#define min(x, y) ((x) < (y) ? (x) : (y))


static Float drtable[4][4][NTAB];
/* ddtab is the amount we multiply a radius by to get an index into
   tabfc and tabdfc.  Note that it used to be an amount we *divide*
   by, but I changed this  5 Aug 2000 because multiplication is
   faster.  Tim Freeman */
static Float ddtab[4][4];
/* rtable holds values for V super R sub i j, defined in equation 4 on page
   7. */
static Float rtable[4][4][NTAB];
/* atable holds values for V super A sub i j, defined in equation 5 on page
   7. */
static Float atable[4][4][NTAB];
/* datable holds values for the derivative of V super A sub i j with respect to
   r. */
static Float datable[4][4][NTAB];
/* tabfc encodes f super c sub ij (r sub ij), mentioned in equations 4 and 5 on
   page 7, defined in equation 18 on page 16. The first coordinates are the
   atomic numbers of the two atoms, and the last coordinate is the radius,
   divided by a constant chosen to make the range right, then truncated to an
   integer. */ 
static Float tabfc[4][4][NTAB];
/* tabdfc has the derivative of f super c sub i j with respect to r. */
static Float tabdfc[4][4][NTAB];

void init_c(BrennerMainInfo *info)
{
  mtable(drtable, ddtab, rtable, atable, datable, tabfc, tabdfc, info->rb2);
}

/* kend is the number of neighbors in the "neighbors" variable.  We've already
   figured out neighbors.  Do the physics for Brenner's bond order potential.
   This routine and all code that it calls should be "pure", in the sense that
   it should only access the state through the api defined in api.h, because
   that gives hope of making this into a Fungimol plugin.

   I can't remember what "ca" stands for; replacing this sentence
   with positive information would be of value.  Tim Freeman  5 Sep 2000*/
Float
caguts(struct State *info)
{
  const int num_atms = pastLastIndex (info);
  // k is the position of the current neighbor in apis.
  int k = 0;
  int l;
  struct NeighborState *cans = caNeighbors (info);
  static struct AtomPairInfoState *apis = 0;
  /* tote is the total energy, for the return value of caguts. */
  Float tote = 0;
  
  if (0 == apis) {
    apis = newAtomPairInfoState ();
  }
  allocateAtomPairAtom (apis, num_atms);
  /* Loop over each bond.  Actually we'll see each bond twice, once from each
     end. */
  {
    int i;
    for (i = firstIndex (info); i < num_atms; i++) {
      const int iNeighborCount = numNeighbors (i, cans);
      const int *const iNeighbors = neighbors (i, cans);
      const int ki = getKtype (i, info) - 1;
      struct AtomPairInfo *pair_i;
      int jn;
      apis->ai_start [i] = k;
      k += iNeighborCount;
      allocate_AtomPairInfo (apis, k);
      pair_i = pairs (i, apis);
      for (jn = 0; jn < iNeighborCount; jn++) {
        const int j = iNeighbors [jn];
        const int kj = getKtype (j, info) - 1;
        int it, floor_rt;
        /* rsq will be the squared length of the bond. */
        Float rsq;
        /* rc will be the length of the bond. */
        Float rc, rt, vv, rp;
        /* We're going to fill in *pair_k. */
        struct AtomPairInfo *thisPair = &(pair_i[jn]);
        thisPair->lcheck = 0;
        rsq = CALC_DIST(info, i, j, &thisPair->cor);
        /* Now rsq is the square of the distance between the two neighbors, and
           thisPair->cor is the vector from j to i. */
        /* Hmm, many pairs of atoms will be far apart by any measure, so we might
           want to have one number that's the maximum rmax to avoid some array
           indexing in this common case.  Small deal. */
        if(rsq > rmax[ki][kj]) continue;
        /* If both are hydrogen or carbon, set lcheck to 1. */
        if(kj <= 1 && ki <= 1) thisPair->lcheck = 1;
        /* If neither are hydrogen or carbon, set lcheck to 2. */
        if(kj >= 2 && ki >= 2) thisPair->lcheck = 2;
        /* otherwise lcheck remains 0. */
        /* set rc to the distance between them. */
        rc = sqrt(rsq);
        /* We used to divide by the value in ddtab here, but I changed the
           definition of ddtab so we multiply now.  Tim Freeman  5 Aug 2000. */
        rt = rc * ddtab[ki][kj];
        floor_rt = (int)floor(rt);
        it = min(floor_rt, NTAB-2);
        thisPair->rcor = rc;
        /* Linear interpolation to compute f super c sub i j (rc), defined in
           equation 18 on page 16. */
        thisPair->ww = tabfc[ki][kj][it]
          + (tabfc[ki][kj][it+1] - tabfc[ki][kj][it])*(rt-it);
        thisPair->dww = tabdfc[ki][kj][it]
          + (tabdfc[ki][kj][it+1] - tabdfc[ki][kj][it])*(rt-it);
        thisPair->exx1 = atable[ki][kj][it]
          + (atable[ki][kj][it+1] - atable[ki][kj][it])*(rt-it);
        thisPair->dexx1 = datable[ki][kj][it] +
          (datable[ki][kj][it+1] - datable[ki][kj][it])*(rt-it);
        if(i >= j) continue;
        /* For the rest of the loop we're only seeing each bond once. */
        vv = rtable[ki][kj][it]
          + (rtable[ki][kj][it+1] - rtable[ki][kj][it])*(rt-it);
        /* Now vv is the energy for the pair due to their pairwise repulsion. */
        rp = drtable[ki][kj][it]
          + (drtable[ki][kj][it+1] - drtable[ki][kj][it])*(rt-it);
        /* Now rp is the magnitude of the force from pairwise repulsion. */
        tote = tote + vv;
        
        /* rpp is the repulsive force between the pair. */
        thisPair->rpp.x = rp*thisPair->cor.x;
        thisPair->rpp.y = rp*thisPair->cor.y;
        thisPair->rpp.z = rp*thisPair->cor.z;
      } /* End loop over neighors. */
    } /* End loop over atoms. */
    assert (num_atms == i);
    apis->ai_start [i] = k;
  } /* Forget i. */
  
  { // FIXME Is there any reason not to merge this loop with the previous one?
    int i;
    for (i = firstIndex (info); i < num_atms; i++) {
      const int iNeighborCount = numNeighbors (i, cans);
      const int *const iNeighbors = neighbors (i, cans);
      const struct AtomPairInfo *const apiNeighbors = pairs (i, apis);
      int jn;
      assert (iNeighborCount == numPairs(i, apis));
      for (jn = 0; jn < iNeighborCount; jn++) {
        const int j = iNeighbors [jn];
        if(apiNeighbors[jn].lcheck == 0) continue;
        if(i >= j) continue;
        /* FIXME The compiler might need some help seeing the common
           subexpressions and the advantages of reordering here. */
        transForcex (i, j, info, apiNeighbors[jn].rpp.x);
        transForcey (i, j, info, apiNeighbors[jn].rpp.y);
        transForcez (i, j, info, apiNeighbors[jn].rpp.z);
      }
    }
  }
  /* If there are hydrogens or carbons, then do pibond. */
  if(getNoa (0, info) + getNoa (1, info) != 0)
    tote += pibond(info, apis);
  /* If there are silicons or germaniums, then do sili_germ. */
  if(getNoa (3, info) + getNoa (4, info) != 0)
    tote += sili_germ(info, apis);
  return tote;
}
