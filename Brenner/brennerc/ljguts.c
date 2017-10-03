/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/
#ifndef __BRENNER_H__
#include "brenner.h"
#endif

#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

/* calc_dist needs State to be defined, so it has to appear after ccapi.h. */
#ifndef __CALC_DIST_H__
#include "calc_dist.h"
#endif

#ifndef __XALLOC_H__
#include "xalloc.h"
#endif

#include <math.h>
/* stdlib.h declares qsort, among other things. */
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

#define MAX_LOOK 10000

static Float rmaxlj[NTYPES+1][NTYPES+1], rslj[NTYPES+1][NTYPES+1]
     ,xmm[NTYPES+1][NTYPES+1],xmms[NTYPES+1][NTYPES+1] 
     ,xm[NTYPES+1][NTYPES+1]
     ,epss[NTYPES][NTYPES],sigs[NTYPES][NTYPES]
     ,surf ,tau[NTYPES]
     ,ndir 
      ,dellj,vlook[MAX_LOOK+1][NTYPES+1][NTYPES+1]
     ,dlook[MAX_LOOK+1][NTYPES+1][NTYPES+1];

static Float eps[NTYPES+1][NTYPES+1];
static Float sig[NTYPES+1][NTYPES+1];

void
ljparam_init_kt(int kt1, Float sigt, Float epst)
{
    sig[kt1][kt1] = sigt;
    eps[kt1][kt1] = epst;
}

/* Set up the spline lookup tables for the lennard-jones potential. */
void ljparam(BrennerMainInfo *info)
{
  int i, j;
  Float rspl[NTYPES+1][NTYPES+1],rspls[NTYPES+1][NTYPES+1];
  Float c2[NTYPES+1][NTYPES+1], c3[NTYPES+1][NTYPES+1];
  /* Evaluate the LJ constants: units are in Ang. and eV */
  
  /* convert K to eV */
  for(i = 1; i <= info->ktmax; ++i) {
    /* FIXME Instead of 11605, this should use the EPSI parameter 
       in bren2.c. */
    eps[i][i] = 4.0*eps[i][i]/11605.0;
  }
  
  for(i = 1; i <= info->ktmax; ++i) {
    for(j = 1; j <= info->ktmax; ++j)
      {
        eps[i][j] = sqrt(eps[i][i]*eps[j][j]);
        sig[i][j] = (sig[i][i]+sig[j][j])/2.0;
        xmm[i][j] = 0.0;
        rspl[i][j] = 0.0;
      }
  }
  
  /*
   * inner and outer points for spline for C-C,H-H,C-H
   */
  for(i = 1; i <= 2; ++i) {
    for(j = 1; j <= 2; ++j) {
      if(eps[i][j] != 0.0) {
        xmm[i][j] = info->rb2[i-1][j-1];
        rspl[i][j] = 0.95*sig[i][j];
      }
    }
  }
  /*
   * inner and outer points for spline for Si-Si,Ge-Ge,Si-Ge
   */
  for(i = 3; i <= 4; ++i) {
    for(j = 3; j <= 4; ++j) {
      if(eps[i][j] != 0.0) {
        xmm[i][j] = info->rb2[i-1][j-1];
        
        rspl[i][j] = 0.95*sig[i][j];
      }
    }
  }
  
  for(i = 1; i <= info->ktmax; ++i) {
    for(j = 1; j <= info->ktmax; ++j) {
      if(eps[i][j] != 0.0) {
        xmms[i][j] = xmm[i][j] - info->RLL;
        xm[i][j] = xmm[i][j];
        rmaxlj[i][j] = 2.50*sig[i][j];
        rslj[i][j] = rmaxlj[i][j]+info->RLL;
        if(0) printf("rslj[%d][%d] %6.3f xmms %6.3f\n",i,j, rslj[i][j],
                     xmms[i][j]);
        rmaxlj[i][j] *= rmaxlj[i][j];
        rslj[i][j] *= rslj[i][j];
        xmm[i][j] *= xmm[i][j];
        xmms[i][j] *= xmms[i][j];
        rspls[i][j] = rspl[i][j]*rspl[i][j];
        sig[i][j] *= sig[i][j];
      } else {
        rmaxlj[i][j] = 0.0;
        rslj[i][j] = 0.0;
      }
    }
  }
  /*
   * find spline parameters for C-C, H-H, and C-H
   */
  for(i = 1; i <= 2; ++i) {
    for(j = 1; j <= 2; ++j) {
      if(eps[i][j] != 0.0) {
        Float dr = rspl[i][j] - info->rb2[i-1][j-1];
        Float rtmp = sig[i][j]/rspls[i][j];
        Float r6 = rtmp*rtmp*rtmp;
        Float vlj = eps[i][j]*r6*(r6-1.0);
        Float dvlj = -eps[i][j]/rspl[i][j]*r6*(12.0*r6 - 6.0);
        c2[i][j] = (3.0/dr*vlj - dvlj)/dr;
        c3[i][j] = (vlj/(dr*dr) - c2[i][j])/dr;
      }
    }
  }
  /*
   * find spline parameters for Si-Si,Ge-Ge,Si-Ge
   */
  for(i = 3; i <= 4; ++i) {
    for(j = 3; j <= 4; ++j) {
      if(eps[i][j] != 0.0) {
        Float dr = rspl[i][j] - info->rb2[i-1][j];
        Float rtmp = sig[i][j]/rspls[i][j];
        Float r6 = rtmp*rtmp*rtmp;
        Float vlj = eps[i][j]*r6*(r6-1.0);
        Float dvlj = -eps[i][j]/rspl[i][j]*r6*(12.0*r6 - 6.0);
        c2[i][j] = (3.0/dr*vlj - dvlj)/dr;
        c3[i][j] = (vlj/(dr*dr) - c2[i][j])/dr;
      }
    }
  }
  
  /*
   * generate table look up for pair interactions
   */
  dellj = 0.001;
  for(i = 1; i <= info->ktmax; ++i) {
    for(j = 1; j <+ info->ktmax; ++j) {
      int k;
      int kmax = (int)floor((sqrt(rmaxlj[i][j]) - xm[i][j])/dellj);
      if(kmax > MAX_LOOK) {
        fprintf(stderr, "kmax overflow %d %d %d\n", kmax,i,j);
        my_exit(-1);
      }
      
      for(k = kmax+1; k <= MAX_LOOK; ++k) {
        vlook[k][i][j] = 0.0;
        dlook[k][i][j] = 0.0;
      }
      
      for(k = kmax; k >= 1; --k) {
        Float r = (k-1)*dellj + xm[i][j];
        Float rsqs;
        Float vdw, dvdw;
        if(((i > 2) || (j > 2)) && (r == 0.0))
          r = dellj;
        rsqs = r*r;
        if(rsqs < rmaxlj[i][j]) {
          if(rsqs >= rspls[i][j]) {
            Float rtmp = sig[i][j]/rsqs;
            Float r6 = rtmp*rtmp*rtmp;
            vdw = eps[i][j] * r6 * (r6 - 1.0);
            dvdw = eps[i][j]/rsqs * r6*(12.0*r6 - 6.0);
          } else {
            Float dr = r - info->rb2[i-1][j-1];
            vdw = dr*dr*(dr*c3[i][j] + c2[i][j]);
            dvdw = -dr*(3.0*dr*c3[i][j] + 2.0*c2[i][j])/r;
          }
        } else {
          vdw = 0.0;
          dvdw = 0.0;
        }
        vlook[k][i][j] =  vdw;
        dlook[k][i][j] = dvdw;
      }
    }
  }
}

/* This is copied from random logic at the end of pred_corr.f that has the
comment "Update Neighbor List?".  

If nobody moved less than RLL/2 from prev_coord, return 2, and leave the
prev_coord's untouched.

If somebody moved more than RLL/2 from prev_coord, return 1, and set the
prev_coords to the present coordinates. */
int
choose_lj(BrennerMainInfo *info)
{
  BrenAtom *atm_num = info->atm_num;
  int num_atms = info->num_atms;
  const Float RLL = info->RLL;
#ifndef INFINITE_CUBE
  const Float *cube = info->cube;
#endif
  Float rsqcrt = (0.5*RLL) * (0.5*RLL);
  int i;
  /* lchk is the return value. */
  int lchk = 0;
  for(i = 0; i < num_atms; ++i) {
    Float rsq = 0.0;
    Float r;
    r = atm_num[i].coord.x - atm_num[i].prev_coord.x;
#ifndef INFINITE_CUBE
    r -= cube[0]*floor(r/cube[0] + 0.5);
#endif
    rsq += r*r;
    r = atm_num[i].coord.y - atm_num[i].prev_coord.y;
#ifndef INFINITE_CUBE
    r -= cube[1]*floor(r/cube[1] + 0.5);
#endif
    rsq += r*r;
    r = atm_num[i].coord.z - atm_num[i].prev_coord.z;
#ifndef INFINITE_CUBE
    r -= cube[2]*floor(r/cube[2] + 0.5);
#endif
    rsq += r*r;
    if(rsq > rsqcrt){ lchk = 1; break; }
  }
  if(lchk == 0) return 2;
  for(i = 0; i < num_atms; ++i)
  {
    atm_num[i].prev_coord.x = atm_num[i].coord.x;
    atm_num[i].prev_coord.y = atm_num[i].coord.y;
    atm_num[i].prev_coord.z = atm_num[i].coord.z;
  }
  assert (1 == lchk);
  return lchk;
}

static void findNeighbors (BrennerMainInfo *state) {
  BrenAtom *atm_num = state->atm_num;
  const int num_atms = state->num_atms;
  const int lchk = state->lchk;
  struct NeighborState *ljns = state->ljNeighborState;
  /*
   * Find the LJ neighbor list
   * Hmm, an n^2 algorithm, and we're doing van-der-waals forces, which have a
   * fairly long interaction distance.  I bet we spend a lot of time here.  Why
   * do we do neighbors twice, once for lj and once for the hydrocarbon
   * potential?  Tim Freeman 24 Jul 2000.
   * We spend most of our time here for molecules with at least a few hundred
   * atoms, and there's no way around it with the standard rslj cutoffs.
   * The very different cuttoffs for the lj and hydrocarbon neighbors
   * makes it hard to efficiently use the same approach for both.
   * - pcm 2000-08-10.
   */
  if(lchk == 1) {
    if (num_atms >= state->min_lj_range_tree) {
      make_neighbor_lists_tree(state, ljns, rslj, xmms);
    } else {
      int *neighbor_start;
      int kli = 0;
      int i, j;
      allocate_atom (ljns, num_atms);
      neighbor_start = ljns->neighbor_start;
      /* Next loop used to go from 0 to num_atms - 2, but the code for the last
         atom is cleaner if it isn't a special case. */
      for(i = 0; i < num_atms; ++i) {
        const BrenAtom *atm_num_i = &atm_num[i];
        const BrenAtom *atm_num_j = atm_num_i + 1;
        int ki = atm_num_i->ktype;
        const Float *rslj_ki = rslj[ki];
        const Float *xmms_ki = xmms[ki];
        neighbor_start [i] = kli;
        allocate_neighbor (ljns, kli + num_atms);
        for(j = i+1; j < num_atms; ++j, ++atm_num_j) /* O(N^2) */ {
          int kj = atm_num_j->ktype;
          Float rsqs;
          if(rslj_ki[kj] == 0.0) continue;
          rsqs = CALC_DIST_NO_RR(toState (state), i, j);
          if(rsqs > rslj_ki[kj]) continue;
          if(rsqs < xmms_ki[kj]) continue;
          /* if((rsqs > rslj_ki[kj]) | (rsqs < xmms_ki[kj])) continue; */
          assert (kli < ljns->neighbors_allocated);
          ljns->neighbor_list[kli] = j;
          ++kli;
        }
      } /* for each atom */
      /* Fill in the start of the atom-after-last, so we get the right number
         of neighbors for the last atom. */
      assert (i == num_atms);
      neighbor_start [num_atms] = kli;
    } /* whether to use n^2 algorithm or n log^3 n */
  } /* if lchk == 1 */
  if(0) printf("lchk %d \n", lchk);
}

/* This routine does the physics for the Lennard-Jones potential.  It only
   accesses the given State via the functions defined in api.h. */
/* FIXME Transplant this to some other file, to sort out the pure stuff. 
   "Pure" meaning code that does physics using only the API.
*/
static Float ljPure (struct State *s) {
  /* pvdw accumulates the total energy so it can be returned as the value of
     ljguts. */
  Float pvdw = 0;
  struct NeighborState *ljns = ljNeighbors (s);
  const int num_atms = pastLastIndex (s);
  int i;
  for (i = 0; i < num_atms; i++) {
    const int numNabr = numNeighbors (i, ljns);
    const int *const nabrs = neighbors (i, ljns);
    int np;
    const int ki = getKtype (i, s);
    for (np = 0; np < numNabr; np++) {
      const int j = nabrs [np];
      const int kj = getKtype (j, s);
      int ii, l;
      
      Float vdw, dvdw, r, dr, rt;
      /* rrs will be the vector from j to i. */
      vector rrs;
      /* rsqs will be the square of the length from i to j. */
      Float rsqs = CALC_DIST(s, i, j, &rrs);
      if(rsqs > rmaxlj[ki][kj]) continue;
      
      if(rsqs < xmm[ki][kj]) continue;
      /* table look up with linear interpolation */
      /* We could go whole hog with the table lookup and have different tables
         for each ki/kj combination.  In this case we could incorporate the
         square root into the function computed by the table, thus eliminating
         the need to compute the square root. */
      r = sqrt(rsqs) - xm[ki][kj];
      rt = r/dellj;
      /* Why ceiling here, and then we look at elements ii and ii+1?  Weird. */
      ii = (int)ceil(rt);
      vdw  = vlook[ii][ki][kj]
        + (vlook[ii+1][ki][kj] - vlook[ii][ki][kj])*(rt-ii+1);
      dvdw = dlook[ii][ki][kj]
        + (dlook[ii+1][ki][kj] - dlook[ii][ki][kj])*(rt-ii+1);
      
      pvdw += vdw;
      {      
        Float dr;
        dr = dvdw*rrs.x;
        setForcex (i, s, getForcex (i, s) + dr);
        setForcex (j, s, getForcex (j, s) - dr);
        dr = dvdw*rrs.y;
        setForcey (i, s, getForcey (i, s) + dr);
        setForcey (j, s, getForcey (j, s) - dr);
        dr = dvdw*rrs.z;
        setForcez (i, s, getForcez (i, s) + dr);
        setForcez (j, s, getForcez (j, s) - dr);
      }
    }
  }
  return pvdw;
}

/* Do the lennard-jones potential. */
Float ljguts(BrennerMainInfo *state)
{
  findNeighbors (state);
  return ljPure (toState (state));
}
