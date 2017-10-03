/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/
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
#include "sili_germ.h"
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>


typedef struct sili_germ
{
  Float xtn2;
  Float xtn1;
  Float adb;
  Float cdb;
  Float cdb2;
  Float ddb;
  Float ddb2;
  Float hdb;
} Sili_Germ;

static Sili_Germ si_ge_parms[4+1];

static int si_ge_initialized;

static void
si_ge_init()
{
	si_ge_parms[3].xtn2 = 0.78734;
	si_ge_parms[3].xtn1 = 1/(2.0*si_ge_parms[3].xtn2);
	si_ge_parms[3].adb = 1.0999e-6;
	si_ge_parms[3].cdb = 1.0039e+5;
	si_ge_parms[3].cdb2 = si_ge_parms[3].cdb*si_ge_parms[3].cdb;
	si_ge_parms[3].ddb = 16.218;
	si_ge_parms[3].ddb2 = si_ge_parms[3].ddb*si_ge_parms[3].ddb;
	si_ge_parms[3].hdb = 0.59826;

	si_ge_parms[4].xtn2 = 0.75627;
	si_ge_parms[4].xtn1 = 1/(2.0*si_ge_parms[4].xtn2);
	si_ge_parms[4].adb = 9.0166e-07;
	si_ge_parms[4].cdb = 1.0643e+5;
	si_ge_parms[4].cdb2 = si_ge_parms[4].cdb*si_ge_parms[4].cdb;
	si_ge_parms[4].ddb = 15.652;
	si_ge_parms[4].ddb2 = si_ge_parms[4].ddb*si_ge_parms[4].ddb;
	si_ge_parms[4].hdb = -0.43884;
	si_ge_initialized = 1;
}

/* I don't have a test case for this code, so I probably broke it.  Tim Freeman
   13 Sep 2000. */
Float
sili_germ(struct State *const info, struct AtomPairInfoState *const apis)
{
  const int num_atms = pastLastIndex (info);
  const struct NeighborState* const cans = caNeighborsconst (info);
  /*
   *     silicon: e = 4.6297255 eV/atom
   *              rnn = 2.3521697 A
   */
  vector xt[250+1];
  Float xslj[250+1],xsij[250+1];
  Float s3, ss, rr, rsq2, rsq3, costh, aarg, aarg2, dgdthet;
  Float rep, rp, vv, vatt, dbdz, gangle;
  Float dctdij, dctdli, dctdlj, dli, bli;
  Float tote = 0;
  int l;
  if(!si_ge_initialized)
    si_ge_init();
  for(l = firstIndex (info); l < num_atms; ++l) {
    const struct AtomPairInfo *const lpairs = pairsconst (l, apis);
    const int *const lNeighbors = neighbors (l, cans);
    const int stop_index = numNeighbors (l, cans);
    const int kl = getKtype (l, info);
    /* i is the index of a neighbor of l. */
    int i;
    for(i = 0; i < stop_index; ++i) {
      Float xsli, ssum;
      int in, ki;
      if (lpairs[i].lcheck != 2) continue;
      in = lNeighbors[i];
      /*  THIS IS THE IL TERM */
      /* There used to be a bug here -- we used i, which was an index into the
         neighbors list, as an index into the atom list to get the ktype of
         something.  They must have meant the ktype of the atom we're talking
         about, which has the index in in the atom list.  Tim Freeman 13 Sep
         2000. */
      ki = getKtype (in, info);
      xsli = 0.0;
      ssum = 0.0;
      /*  RSQ1 IS THE LI TERM */
      if(stop_index > 1) {
        Float s1 = lpairs[i].rcor;
        Float rsq1 = s1*s1;
        /* j is also an index into the list of neighbors of l. */
        int j;
        for(j = 0; j < stop_index; ++j) {
          if(j == i) continue;
          if(lpairs[j].lcheck != 2) continue;
          /*  RSQ3 IS THE LJ TERM */
          s3 = lpairs[j].rcor;
          rsq3 = s3*s3;
          /*   RSQ2 IS THE IJ TERM */
          rsq2 = 0.0;
          
          xt[j].x = lpairs[j].cor.x - lpairs[i].cor.x;
          rsq2 += xt[j].x*xt[j].x;
          xt[j].y = lpairs[j].cor.y - lpairs[i].cor.y;
          rsq2 += xt[j].y*xt[j].y;
          xt[j].z = lpairs[j].cor.z - lpairs[i].cor.z;
          rsq2 += xt[j].z*xt[j].z;
          ss = 2.0*s1*s3;
          rr = rsq1-rsq3;
          costh = (rsq1+rsq3-rsq2)/ss;
          aarg = si_ge_parms[kl].hdb + costh;
          aarg2 = si_ge_parms[kl].ddb2 + aarg*aarg;
          gangle = 1 + si_ge_parms[kl].cdb2*(1/si_ge_parms[kl].ddb2 - 1/aarg2);
          /*  DG / DCOSTHETA */
          dgdthet = 2.*si_ge_parms[kl].cdb2*aarg/(aarg2*aarg2);
          /*  DCOSTHETA / DRIJ  * (1./RIJ) */
          dctdij = -2.0/ss;
          dctdli = (rr+rsq2)/(ss*rsq1);
          dctdlj = (-rr+rsq2)/(ss*rsq3);
          {
            const Float exx = si_ge_parms[kl].adb * lpairs[j].ww;
            ssum += gangle*exx;
            /*  THE XSIJ  TERMS ARE GOING TO BE  DZ/DRIJ */
            {
              Float xtemp = exx*dgdthet;
              xsli += xtemp*dctdli;
              xslj[j] = exx * gangle * lpairs[j].dww / s3 + xtemp*dctdlj;
              xsij[j] = xtemp*dctdij;
            }
          }            
        }
      }
      dli = (1.0 + pow(ssum, si_ge_parms[kl].xtn2));
      bli = pow(dli, -si_ge_parms[kl].xtn1);
      
      dbdz = 0.0;
      if(ssum != 0.0) {
        dbdz = -si_ge_parms[kl].xtn1*bli/dli*si_ge_parms[kl].xtn2
          * pow(ssum, si_ge_parms[kl].xtn2 - 1.0);
      }
      vatt = lpairs[i].exx1;
      vv = -vatt*bli;
      tote += vv;
      if(0) printf("si ge tote %2d %10.7f %10.7f %10.7f %10.7f\n",
                   i, tote, vv, ssum, vatt);
      
      /*  LI: */
      rp = vatt* dbdz*xsli + bli * lpairs[i].dexx1;
      
      transForcex (l, in, info, rp*lpairs[i].cor.x);
      transForcey (l, in, info, rp*lpairs[i].cor.y);
      transForcez (l, in, info, rp*lpairs[i].cor.z);
      
      if(stop_index > 1) {
        int j;
        for(j = 0; j < stop_index; ++j) {
          int jn;
          if(i == j || lpairs[j].lcheck != 2) continue;
          jn = lNeighbors[j];
          /*  LJ: */
          rp = vatt*dbdz*xslj[j];
          
          transForcex (l, jn, info, rp*lpairs[j].cor.x);
          transForcey (l, jn, info, rp*lpairs[j].cor.y);
          transForcez (l, jn, info, rp*lpairs[j].cor.z);
          
          /*  IJ: */
          rp = vatt*dbdz*xsij[j];
          
          transForcex (in, jn, info, rp*xt[j].x);
          transForcey (in, jn, info, rp*xt[j].y);
          transForcez (in, jn, info, rp*xt[j].z);
        }
      }
    }
  }
  return tote;
}
