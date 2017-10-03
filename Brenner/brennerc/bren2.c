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

/* Reads in coordinates of atoms in Hyperchem, determines their
   energies, steps them in space to lower their energies, and updates
   Hyperchem coordinates.  Also maybe accept keyboard input for velocities
   of and forces on atoms.  Uses Brenner's potential: Phys. Rev. B (1990) v42 n15 p9458-70.
   Based on new version from http://www.engr.ncsu.edu/mat/CompMatSci/projects.html
   File HCnewpot1.ps is unpublished Brenner paper describing the new version. */

#ifndef __RMAX_H__
#include "rmax.h"
#endif

#include "brenner.h"
#include "sili_germ.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include "spgch.h"
typedef unsigned char BYTE;
typedef long DWORD;
typedef char *LPSTR;

#ifndef __BOOL_H__
#include "bool.h"
#endif

#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

#ifndef __CALC_DIST_H__
#include "calc_dist.h"
#endif

#ifndef __CCNEIGHBORSTATE_H__
#include "CCNeighborState.h"
#endif

typedef void VOID;
#define RBOXES 64
/*#define BOXWIDTH 1.7 */
#define BOXWIDTH 4
#define BOXATMS 80
#define JITTER 1.e-3

#define BOLZ 1.380662
#define EPSI0 11604.5 /* EPSI as initialized in setin.f */
#define EPSI 11605.0
#define AVO 6.02205
#define ECONV ((1/(AVO*BOLZ*EPSI0))*1e7)

/*Debye temperature converted to fs-1 X 2Pi*/
#define WD (2230*(2.08365e-5)*2*M_PI)
/*LANGEVIN PARAMETERS*/

/* A box has a list of atoms in it.  There are RBOXES^3 of them. */

struct BOX { int list[BOXATMS]; int numAtoms;};
typedef struct BOX SB;
static SB box[RBOXES][RBOXES][RBOXES];

static Float RLIST[NTYPES+1][NTYPES+1];

/* The code below reads this from the third line of input.d.  (It used to read
   it from the second line, but that makes things incompatible with the sample
   input.d that came with the Fortran version, so I revised it to read from the
   third line.)  The comment at this place in the sample input.d from
   Brennermd's Execute directory says this is "=1 REBO (C,H,Si,Ge), =2
   tight-binding for C".  I disagree with Peter's comment (which I have
   deleted) saying that it's the periodic box length.  The dimensions of the
   periodic box are the fourth line of coord.d.  Tim Freeman 20 Jul 2000.*/

static int IPOT;

static int most_atoms_in_a_box;

static int IGC[26], IGH[26];

static void initkt(BrennerMainInfo *info, int *kflag);
void append_file(FILE *fp, Float TTIME, Float DELTA, Float *cube);
static void initialize_atoms(BrennerMainInfo *info, int kflag);
static void apply_thermostat(BrennerMainInfo *);
static void third_order_corrector(BrennerMainInfo *);
static void third_order_predictor(BrennerMainInfo *);
static void make_neighbor_lists(BrennerMainInfo *);
static int look_in_this(SB *box1, BrenAtom *pAtm, BrennerMainInfo *state,
                        int kend);
static void find_boxes(BrennerMainInfo *);
static void PARAM(BrennerMainInfo * info);

/* Get coefficients for H and F and torsion. */

static void PARAM(BrennerMainInfo * info)
{
  /*Parameters for hydrocarbons*/
  int I,IC,I2D,J,K,L,M,N;
  int i,j;
  Float XHH;
  FILE *inter2D, *inter3Dh, *inter3Dch;
  
  for (I=1; I<17; I++)
    IGC[I] = 4;
  for (I=17; I<19; I++)
    IGC[I] = 3;
  for (I=19; I<21; I++)
    IGC[I] = 2;
  for (I=21; I<26; I++) 
    IGC[I] = 1;
  
  for (I=1; I<19; I++)
    IGH[I] = 3;
  for (I=19; I<23; I++)
    IGH[I] = 2;
  for (I=23; I<26; I++)
    IGH[I] = 1;
  
  if ((inter2D = fopen("inter2d_iv.d","r"))==NULL)
	{
      fprintf(stderr, "opening inter2d_iv failed\n");
      my_exit(-1);
	}
  fscanf(inter2D, "%d\n",&I2D);
  
  /*Read bicubic spline coefficients*/
  init_xh(inter2D);
  init_in2(inter2D);
  fclose(inter2D);
  for(i = 0; i < 4; ++i)
    {
      for(j = 0; j < 4; ++j)
        {
          Float r2 = info->rb2[i][j];
          RLIST[i+1][j+1] = (r2+info->RLL)*(r2+info->RLL);
#ifndef NDEBUG
          if(0) printf("rlist[%d][%d] %9.4f %9.4f %9.4f rmax %7.3f\n",
                       i,j,RLIST[i+1][j+1],r2,info->RLL, r2*r2);
#endif
          rmax[i][j] = r2*r2;
        }
    }
  ljparam(info);
}

static int cmp_n(const void *v1, const void *v2)
{
  const IntPair *i1 = (const IntPair*)v1;
  const IntPair *i2 = (const IntPair*)v2;
  return i1->v1 == i2->v1 ? (i1->v2 - i2->v2) : (i1->v1 - i2->v1);
}
 
/* Compute neighbor lists if there are less than 200 atoms.  Straightforward
   n-squared algorithm, where we look at each atom and then look at all
   other atoms as possible neighbors.  The return value will be the number of
   neighbors that we found. */
static int
n2_neighbor_list(BrennerMainInfo *info)
{
  struct NeighborState *cans = info->caNeighborState;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  /* We will return k. */
  int k = 0;
  int i, j;
  allocate_atom (cans, num_atms);
  for(i = 0; i < num_atms; ++i) {
    const Float *rlist_ki;
    int ki = atm_num[i].ktype;
    allocate_neighbor (cans, k+num_atms);
    cans->neighbor_start[i] = k;
    
    /*
     * cuts out all but C,H,Si, and Ge
     */
    if(ki >= 5) continue;
    rlist_ki = RLIST[ki];
    
    /* We're not making use of the fact that if atom i is a neighbor of
       j, then j is a neighbor of i.  This loop could probably start at j=i+1
       and run twice as fast. Tim Freeman 23 Jul 2000. */
    /* Added code to do that if N2_ASYM defined; the speed difference is
       too small to measure for typical molecules. pcm 2000-08-11 */
    /* And the new neighbor data structure has to be created in order, so the
       code Peter put in doesn't work any more and it isn't obvious how to fix
       it.  Thus I'm ripping out the code I badgered Peter into putting in
       here.  Sorry.  Tim Freeman 28 Aug 2000. */
    for(j = 0; j < num_atms; ++j) {
      int kj;
      Float rlis, rsq;
      if(i == j) continue;
      kj = atm_num[j].ktype;
      /*
       * cuts out all but C,H,Si, and Ge
       */
      if(kj >= 5) continue;
      rlis = rlist_ki[kj];
      rsq = CALC_DIST_NO_RR(toState (info), i, j);
      if(rsq > rlis) continue;
      cans->neighbor_list[k] = j;
      ++k;
    }
  }
  cans->neighbor_start[i] = k;
  return k;
}

/* Sort the atoms into boxes.  The boxes are in the file-scope array "box"
   declared above.
*/
static void find_boxes(BrennerMainInfo *info)
{
  int i, j, k;
  int x,y,z;
  int num_atms = info->num_atms;
  const BrenAtom *const atm_num = info->atm_num;
  static SB *boxes_used[MAX_ATOMS];
  static int num_boxes_used;
#ifndef INFINITE_CUBE
  const Float *cube = info->cube;
#endif
  
  /* initialize boxes to empty */
  for(i = 0; i < num_boxes_used; ++i)
    boxes_used[i]->numAtoms = 0;
  /* assign atoms to boxes */ 
  num_boxes_used = 0;
  for (i=0; i<num_atms; i++) {
    Double t;
    /* move the atoms +RBOXES/2 in every coordinate to put the atoms
       in the middle of the box group */
    t = atm_num[i].coord.x;
#ifndef INFINITE_CUBE
    t -= cube[0]*floor(t/cube[0] + 0.5);
#endif
    x = (int)floor(t/BOXWIDTH)+RBOXES/2;
    t = atm_num[i].coord.y;
#ifndef INFINITE_CUBE
    t -= cube[1]*floor(t/cube[1] + 0.5);
#endif
    y = (int)floor(t/BOXWIDTH)+RBOXES/2;
    t = atm_num[i].coord.z;
#ifndef INFINITE_CUBE
    t -= cube[2]*floor(t/cube[2] + 0.5);
#endif
    z = (int)floor(t/BOXWIDTH)+RBOXES/2;
    if (box[x][y][z].numAtoms < BOXATMS) {
      box[x][y][z].list[box[x][y][z].numAtoms] = i;
      if(!box[x][y][z].numAtoms)
        boxes_used[num_boxes_used++] = &box[x][y][z];
      if(++box[x][y][z].numAtoms > most_atoms_in_a_box) {
        most_atoms_in_a_box = box[x][y][z].numAtoms;
        if(most_atoms_in_a_box > 48)
          printf("most_atoms_in_a_box %d\n", most_atoms_in_a_box);
      }
    } else {
      fprintf(stderr,
              "Too many atoms in box %d at %d, %d, %d\n",
              box[x][y][z].numAtoms, x, y, z);
      fprintf(stderr, "atom %d coords %.2f %.2f %.2f\n",
              i, atm_num[i].coord.x, atm_num[i].coord.y,
              atm_num[i].coord.z);
#ifdef INFINITE_CUBE
      fprintf(stderr,"This may indicate a need to increase RBOXES\n");
#endif
      my_exit(-1);
    }
  }
}

/* Given an atom and a box, read box.list[] into the neighbor structure.
   The caller should have allocated enough room in the neighbor structure for
   this already.
   Return the new number of neighbors.
*/
static inline int
look_in_this(SB *box1, BrenAtom *pAtm, BrennerMainInfo *state, int kend)
{
	int i;
	int ki = pAtm->ktype;
	const Float *rlist_ki = RLIST[ki];
    const BrenAtom *const atm_num = state->atm_num;
	for (i = 0; i < box1->numAtoms; i++) {
		int index2 = box1->list[i];
		const BrenAtom *pAtm2 = &atm_num[index2];
		int kj = pAtm2->ktype;
		Float rlis, rsq;
		/*
		 * cuts out all but C,H,Si, and Ge
		 */
		if(kj >= 5) continue;
        /* FIXME We could put them in the neighbor list regardless of distance,
           and filter out the bogus ones when we scan the neighbor list later.
           This would make for longer neighbor lists, with the benefit of one
           less distance computation.  I don't know whether this is a net win.
           Peter says he'll replace this with a range tree anyway, so I'm not
           interested in trying this optimization soon.  Tim Freeman 28 Aug
           2000.
           But the code that would filter out the neighbors that are too far
           apart is already present in caguts.c, so the only cost of the
           experiment should be ifdef'ing out the code here.  Maybe worth
           trying anyway once I've transformed everything to use the API.  Tim
           Freeman 28 Aug 2000. */
		rlis = rlist_ki[kj];
		rsq = CALC_DIST_NO_RR(toState (state), pAtm->number, index2);
        /* If it's too far away, leave it out. */
		if(rsq > rlis) continue;
        /* The following if statement ensures that no atom is counted ase its
           own neighbor. */
		if (index2 != pAtm->number) {
          /* FIXME I hope the compiler figures out that
             state->caNeighborState->neighbor_list is constant and hoists it
             out of this inline code.  Seems dubious.  Might want to check the
             assembly language.
             Tim Freeman 28 Aug 2000 */
          state->caNeighborState->neighbor_list [kend] = index2;
          ++kend;
		}
	}
	return kend;
}

int
out_of_box(double x)
{
  return fabs(x) > RBOXES*BOXWIDTH/2;
}

/* Fill in neighbor_list to reflect who is a neighbor of who.  At the end, for
   each atom (suppose it's number is i), info->atom_num[i].num_neighbors will
   be the number of neighbors of atom i.  The neighbors will be represented as
   consecutive elements of neighbor_list.  The first one will be neighbor_list
   [info->atom_num[i].neighbor_start].  If j is a neighbor of i, then one
   element of neighbor_list will have v1 equal to i and v2 equal to j.  In the
   likely event that i is also a neighbor of j, then another element of
   neighbor_list will apparently have v1 equal to j and v2 equal to i. 
   Tim Freeman 24 Jul 2000. */
static void make_neighbor_lists(BrennerMainInfo *info)
{
	int i,j,k,n,x,y,z;
	int kend = 0;
	BrenAtom *pAtom;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
    struct NeighborState *cans = info->caNeighborState;
#ifndef INFINITE_CUBE
	const Float *cube = info->cube;
	const int cubex = (int)(info->cube[0]/BOXWIDTH) + RBOXES/2;
	const int cubey = (int)(info->cube[1]/BOXWIDTH) + RBOXES/2;
	const int cubez = (int)(info->cube[2]/BOXWIDTH) + RBOXES/2;	
#endif
	if(num_atms >= info->min_range_tree
#ifndef INFINITE_CUBE 
       || out_of_box(info->cube[0])
	   || out_of_box(info->cube[1])
       || out_of_box(info->cube[2])
#endif
	   )
	{
	  static Float xmms[NTYPES+1][NTYPES+1]; /* initialized by the compiler to
                                                all zeros */
	  n2_neighbor_list(info);
	} else if(num_atms < 50) {
      n2_neighbor_list(info);
    } else {
      find_boxes(info);
      /*find the neighbor list for each atom by going through 
        the 27 surrouding boxes. */
      allocate_atom (cans, num_atms);
      for (n=0; n<num_atms; n++) {
        int xplus1, yplus1, zplus1;
        int xminus1, yminus1, zminus1;
        double t;
        allocate_neighbor (cans, kend + num_atms);
        cans->neighbor_start [n] = kend;
        pAtom = &atm_num[n];
        if(pAtom->ktype >= 5) continue;
        t = atm_num[n].coord.x;
#ifndef INFINITE_CUBE
        t -= cube[0]*floor(t/cube[0] + 0.5);
#endif
        x = (int)floor(t/BOXWIDTH)+RBOXES/2;
        /* Peter's code here used to care what the cube was, even in the infinite
           cube case.  I don't know the right thing to do in that case.  I think
           Peter's code was wrong in the infinite cube case because xminus1 might
           be negative if the x coordinate has become sufficiently negative, thus
           causing us to reference our array out of bounds.  Instead I'll use
           simpler code that I absolutely know is wrong, by just assuming the
           entire scene will fit within the grid of boxes.  This code will simply
           do the wrong thing if atoms in the scene escape from the region that
           we've made boxes for.  Tim Freeman 26 Aug 2000. */
#ifdef INFINITE_CUBE
        xminus1 = x - 1;
        xplus1 = x + 1;
#else
        xminus1 = x ? x-1 : cubex - 1;
        xplus1 = x+1 >= cubex ? 0 : x+1;
#endif
        t = atm_num[n].coord.y;
#ifndef INFINITE_CUBE
        t -= cube[1]*floor(t/cube[1] + 0.5);
#endif
        y = (int)floor(t/BOXWIDTH)+RBOXES/2;
#ifdef INFINITE_CUBE
        yminus1 = y - 1;
        yplus1 = y + 1;
#else
        yminus1 = y ? y-1 : cubey - 1;
        yplus1 = y+1 >= cubey ? 0 : y+1;
#endif
        t = atm_num[n].coord.z;
#ifndef INFINITE_CUBE
        t -= cube[2]*floor(t/cube[2] + 0.5);
#endif
        z = (int)floor(t/BOXWIDTH)+RBOXES/2;
#ifdef INFINITE_CUBE
        zminus1 = z - 1;
        zplus1 = z + 1;
#else
        zminus1 = z ? z-1 : cubez - 1;
        zplus1 = z+1 >= cubez ? 0 : z+1;
#endif
        kend = look_in_this(&box[xminus1][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xminus1][yplus1][zplus1],	pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[x][yplus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yminus1][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][y][zplus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][zminus1], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][z], pAtom,
                            info, kend);
        kend = look_in_this(&box[xplus1][yplus1][zplus1], pAtom,
                            info, kend);
      }
      cans->neighbor_start [n] = kend;
      // The new neighbor data structure makes the old sorting code
      // unusable.  Take it out, and put in something new if we really need
      // sorting.  Tim Freeman  6 Sep 2000.
      // qsort(info->neighbor_list,kend, sizeof(info->neighbor_list[0]),
      // cmp_n); 
    }
}

static void third_order_predictor(BrennerMainInfo *info)
{
  int i;
  int num_atms = info->num_atms;
  BrenAtom *atm_num = info->atm_num;
  
  for (i=0; i<num_atms; i++) 
    if (atm_num[i].movable) {
      atm_num[i].coord =
        dplus(dplus(dplus(atm_num[i].coord,atm_num[i].velocity),
                    atm_num[i].accel),
              atm_num[i].dx3dt3);
      atm_num[i].velocity =
        plus(plus(atm_num[i].velocity,product(2,atm_num[i].accel)),
             product(3,atm_num[i].dx3dt3));
      atm_num[i].accel = plus(atm_num[i].accel,product(3,atm_num[i].dx3dt3));
    }						
}

static void third_order_corrector(BrennerMainInfo *info)
{
	int i;
	Float DE;
	vector RI;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
	
	DE=(info->timestep*info->timestep/2)/ECONV;

	for (i=0; i<num_atms; i++) 
	  if (atm_num[i].movable) {
	    RI = plus(atm_num[i].accel,product(-DE/atm_num[i].mass,atm_num[i].force));
	    atm_num[i].coord = dplus(atm_num[i].coord, product(-1.0/6,RI));
#ifndef NDEBUG
	    if(!(atm_num[i].coord.x > -1.e9))
	    {
	      fprintf(stderr, "absurd coord %d %e xforce %e a %e mass %e\n",
		      i, RI.x, atm_num[i].force.x, atm_num[i].accel.x, atm_num[i].mass);
	      my_exit(-1);
	    }
#endif
	    atm_num[i].velocity = plus(atm_num[i].velocity, product(-5.0/6,RI));
	    atm_num[i].accel = plus(atm_num[i].accel, product(-1.0,RI));
	    atm_num[i].dx3dt3 = plus(atm_num[i].dx3dt3, product(-1/3.0,RI));
	  }
}

/* FRICTION AND RANDOM FORCE */
/* I'm uncertain whether this works as well as the Fortran version; the
 * overall behavior looks suspiscious, but looking at individual changes
 * I don't see anything that can be explained by randomness. pcm 2000-08-30.
 */
static void
apply_gleq_thermostat(BrennerMainInfo *info)
{
  int i;
  int ii;
  int nta = 0;
  int nlr;
  const double PI2 = M_PI*2;
  const Float TR = info->temperature/EPSI/ECONV;
  const Float BET = WD*M_PI*ECONV/6.0/info->timestep;
  const Float GSIG = sqrt(2.0*TR*ECONV*BET);
  Float *gl;
  int *nlist = (int *)malloc(info->num_atms * sizeof(*nlist));

  for(i = 0; i < info->num_atms; ++i) 
    if(info->atm_num[i].thermostated)
    {
      nlist[nta++] = i;
    }
  nlr = 3*nta/2;
  gl = (Float *)malloc(2*nlr*sizeof(gl[0]));

  for(i = 0; i < nlr; ++i)
  {
    double rr = rand()/(double)RAND_MAX;
    if(rr >= 1.0e-6)
    {
      float pre = sqrt(-2.0*log(rr));
      gl[i] = pre*GSIG*(Float)cos(PI2*rand()/RAND_MAX);
      gl[i+nlr] = pre*GSIG*(Float)sin(PI2*rand()/RAND_MAX);
    }
  }

  for(ii = 0; ii < nta; ++ii)
  {
    i = nlist[ii];
    if(info->atm_num[i].thermostated)
    {
      Float mass = info->xmass[info->atm_num[i].ktype];
      Float bm = BET*mass;
      Float sm = sqrt(mass);
      info->atm_num[i].force.x -= bm*info->atm_num[i].velocity.x + sm*gl[ii];
      info->atm_num[i].force.y -= bm*info->atm_num[i].velocity.y + sm*gl[ii+nta];
      info->atm_num[i].force.z -= bm*info->atm_num[i].velocity.z + sm*gl[ii+2*nta];
    }
  }
  free(gl);
  free(nlist);
}

/* USE EVANS-HOOVER SCHEME */

static void
apply_hoov_thermostat(BrennerMainInfo *info)
{
  /* this used for all atoms */
  int i;
  Float sc;
  Float ff = 0.0;
  Float df = 0.0;
  for(i = 0; i < info->num_atms; ++i)
  {
    Float mass = info->xmass[info->atm_num[i].ktype];
    Float velx = info->atm_num[i].velocity.x;
    Float vely = info->atm_num[i].velocity.y;
    Float velz = info->atm_num[i].velocity.z;
    ff += info->atm_num[i].force.x*velx;
    df += velx*velx*mass;
    ff += info->atm_num[i].force.y*vely;
    df += vely*vely*mass;
    ff += info->atm_num[i].force.z*velz;
    df += velz*velz*mass;
  }
  if(df == 0)
  {
    fprintf(stderr, "Error in hoov thermostat - atom velocities are all zero?\n");
    my_exit(-1);
  }
  sc = ff/df;
  for(i = 0; i < info->num_atms; ++i)
  {
    Float mass = info->xmass[info->atm_num[i].ktype];
    info->atm_num[i].force.x -= sc*info->atm_num[i].velocity.x*mass;
    info->atm_num[i].force.y -= sc*info->atm_num[i].velocity.y*mass;
    info->atm_num[i].force.z -= sc*info->atm_num[i].velocity.z*mass;
  }
}

static void apply_thermostat(BrennerMainInfo *info)	/*USE BERENDSEN SCHEME*/
{
	int i;
	Float XX, SC, SM;
	Float TR = info->temperature/EPSI/ECONV;
	Float BET = (WD*M_PI*ECONV/6/info->timestep);
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;

	XX=0;
	for (i=0; i<num_atms; i++) 
		if (atm_num[i].thermostated)
			XX += dot(atm_num[i].velocity,atm_num[i].velocity)*atm_num[i].mass;

	if (XX<1e-7) 
	{
		fprintf(stderr, "T=0, Reset Thermostat to other than 1, fatal");
		exit(-1);
	}

	SC=BET*(TR*info->dnla/XX-1);

	for (i=0; i<num_atms; i++) 
		if (atm_num[i].thermostated) {
			SM=atm_num[i].mass*SC;
			atm_num[i].force = plus(atm_num[i].force,product(SM,atm_num[i].velocity));
			if(0)
			  printf("therm %2d %9.6f V %10.8f\n",
			       i, atm_num[i].force.x, atm_num[i].velocity.x);
		}
}

/* Despite the name, this finds forces too.  This is not static because it is
   also called from minimize.c. */
void
find_atom_energies(BrennerMainInfo *info)
{
	int i;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
	static int cnt1;

	if(info->lchk == 1) {
	   make_neighbor_lists(info);
#ifndef NDEBUG
       sortNeighbors (info->caNeighborState, num_atms);
#endif
    }
	info->system_energy = 0;
	for (i=0; i<num_atms; i++) {
		atm_num[i].force.x = 0;
		atm_num[i].force.y = 0;
		atm_num[i].force.z = 0;
	}

	info->system_energy += caguts(toState (info));
	info->system_energy += ljguts(info);
}

int
write_coord_file(FILE *fp, BrennerMainInfo *info)
{
  int i;
  const BrenAtom *atm_num = info->atm_num;
  fprintf(fp, "%s\n", info->head);
  fprintf(fp, "%6d\n", info->num_atms);
  fprintf(fp, "%20.11e %20.11e\n", info->starttime, info->timestep);
#ifdef INFINITE_CUBE
  fprintf(fp, "%20.11e %20.11e %20.11e\n", 1000.0, 1000.0, 1000.0);
#else
  fprintf(fp, "%20.11e %20.11e %20.11e\n", info->cube[0], info->cube[1],
	  info->cube[2]);
#endif
  for(i = 0; i < info->num_atms; ++i)
  {
    fprintf(fp,"%5d %5d %20.11e %20.11e %20.11e %3d\n",
	    i+1, info->kt2[atm_num[i].ktype], atm_num[i].coord.x, 
	    atm_num[i].coord.y, atm_num[i].coord.z,
	    atm_num[i].movable ? atm_num[i].thermostated : 2);
  }

  for(i = 0; i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].velocity.x, 
	    atm_num[i].velocity.y, atm_num[i].velocity.z);
  }

  for(i = 0; i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].accel.x, 
	    atm_num[i].accel.y, atm_num[i].accel.z);
  }

  for(i = 0; i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, atm_num[i].dx3dt3.x, 
	    atm_num[i].dx3dt3.y, atm_num[i].dx3dt3.z);
  }
  for(i = 0; i < info->num_atms; ++i)
  {
    fprintf(fp, "%5d %20.11e %20.11e %20.11e\n", i+1, 0.0, 0.0, 0.0);
  }
  return 1;
}

static void initialize_atoms(BrennerMainInfo *info, int kflag)
{
	int i;
	Float small;
	int num_atms = info->num_atms;
	BrenAtom *atm_num = info->atm_num;
    
    /* FIXME This value for small won't work in single-precision. */
	small=1e-12;
	for (i=0; i<num_atms; i++) 
      if (atm_num[i].movable==1){
        if(kflag != 4) {
          /* Add some noise, so that when we try to do the thermostat we
             never have zero temperature, since that leads to dividing by
             zero. */
          atm_num[i].velocity.x = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].velocity.y = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].velocity.z = (rand()/(Float)RAND_MAX)*2*JITTER-JITTER;
          atm_num[i].accel.x = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].accel.y = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].accel.z = (rand()/(Float)RAND_MAX)*2*small-small;
          
          atm_num[i].dx3dt3.x = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].dx3dt3.y = (rand()/(Float)RAND_MAX)*2*small-small;
          atm_num[i].dx3dt3.z = (rand()/(Float)RAND_MAX)*2*small-small;
        }
        atm_num[i].prev_coord.x = atm_num[i].coord.x;
        atm_num[i].prev_coord.y = atm_num[i].coord.y;
        atm_num[i].prev_coord.z = atm_num[i].coord.z;
      }
}

#ifndef __SAFE_FGETS_H__
#include "safe_fgets.h"
#endif

static void
initkt (BrennerMainInfo *info, int *kflag)
{
  char buf[128];
  int i;
  int kuc, maxkb;
  FILE *fp = fopen("input.d", "r");
  if(!fp)
  {
    fprintf(stderr, "can't read input.d\n");
    my_exit(-1);
  }
  safe_fgets(buf, sizeof(buf), fp, "input.d line 1");
  if(sscanf(buf, "%d %d %d %d", &kuc, &maxkb, kflag, &info->nxmol) != 4)
  {
    fprintf(stderr, "parse error in input.d line 1\n");
    my_exit(-1);
  }
  /* Random # seed, neighbor list, temperature(K) */
  safe_fgets(buf, sizeof(buf), fp, "input.d line 2");
  /* I want to be able to feed the same files to the distributed Fortran
     brennermd program and to this program, so we need compatibility.  This
     code used to demand that ipot be on the same line as the rest.  Note that
     we never use IPOT.  Tim Freeman 24 Jul 2000*/
  {
    double pseed, rll, temp;
    if(sscanf(buf, "%lf %lf %lf", &pseed, &rll, &temp) != 3) {
      fprintf(stderr, "parse error in input.d line 2\n");
      my_exit(-1);
    }
    info->PSEED = pseed;
    info->RLL = rll;
    info->temperature = temp;
  }
  safe_fgets(buf, sizeof(buf), fp, "unused potential flag");
  if(sscanf(buf, "%d", &IPOT) != 1)
    {
      fprintf(stderr, "parse error in input.d line 3\n");
      my_exit(-1);
    }
  for(i = 0; i <= NTYPES; ++i)
    info->kt[i] = info->kt2[i] = 0;
  info->kt[HYDROGEN] = 2;
  info->kt[CARBON] = 1;
  /* The next two used to write to kt out of bounds, back before I made it have
     the size MAX_ATOMNO+1 instead of NTYPES+1.  Tim Freeman  5 Oct 2000. */
  info->kt[SILICON] = 3;
  info->kt[GERMANIUM] = 4;
  info->kt2[2] = HYDROGEN;
  info->kt2[1] = CARBON;
  info->kt2[3] = SILICON;
  info->kt2[4] = GERMANIUM;
  info->ktmax = 4;
  while(fgets(buf, sizeof(buf), fp)) {
      int natom;
      Float xma, epst, sigt;
      {
        double xmad, epstd, sigtd;
        if (sscanf(buf, "%d %lf %lf %lf", &natom, &xmad, &epstd , &sigtd) != 4) {
          fprintf (stderr, "Parse error in Lennard-Jones paramters for this "
                   "line:\n%s\n", buf);
          my_exit (-1);
        }
        xma = xmad;
        epst = epstd;
        sigt = sigtd;
      }
      if(natom < 0)
        continue;
      if(info->kt[natom] == 0) {
        if(info->ktmax > NTYPES) {
          fprintf(stderr, "Error - number of types (%d) exceeds limit\n", info->ktmax);
          exit(-1);
        }
        info->kt[natom] = info->ktmax;
        info->kt2[info->ktmax] = natom;
      }
      info->xmass[info->kt[natom]] = xma;
      ljparam_init_kt(info->kt[natom], sigt, epst);
      ++info->ktmax;
  }
  fclose(fp);
}

static void
remove_com_velocity(BrennerMainInfo *info)
{
  Float xmt = 0.0;
  int i;
  vector com;
/*
****IMPORTANT********************************
*                                           *
* IF NO RIGID ATOMS, SUBTRACT COM VELOCITY  *
*                                           *
* REMOVE THIS CODE FOR MOLECULAR COLLISION  *
*                                           *
*********************************************
*/
  for(i = 0; i < info->num_atms; ++i)
  {
    if(!info->atm_num[i].movable)
      return;
    xmt += info->xmass[info->atm_num[i].ktype];
  }
  com.x = com.y = com.z = 0.0;
  for(i = 0; i < info->num_atms; ++i)
  {
    com.x += info->atm_num[i].velocity.x*info->xmass[info->atm_num[i].ktype];
    com.y += info->atm_num[i].velocity.y*info->xmass[info->atm_num[i].ktype];
    com.z += info->atm_num[i].velocity.z*info->xmass[info->atm_num[i].ktype];
  }
  com.x /= xmt;
  com.y /= xmt;
  com.z /= xmt;
  for(i = 0; i < info->num_atms; ++i)
  {
    info->atm_num[i].velocity.x -= com.x;
    info->atm_num[i].velocity.y -= com.y;
    info->atm_num[i].velocity.z -= com.z;
  }
}

void
bren_1_step(BrennerMainInfo *info, int kflag)
{
	third_order_predictor(info);
	find_atom_energies(info);

	if(kflag == -1)
	  apply_gleq_thermostat(info);
	else if(kflag == 2)
	{
	  int i;
	  for(i = 0; i < info->num_atms; ++i)
	    info->atm_num[i].velocity.x = info->atm_num[i].velocity.y = info->atm_num[i].velocity.z = 0.0;
	}
	else if(kflag == 3)
	  apply_hoov_thermostat(info);
	else if (kflag != 4) {
	  if(info->num_thermostated_atoms > 0)
	    apply_thermostat(info);
	}
	third_order_corrector(info);
	info->lchk = choose_lj(info);
	if(info->zero_com_velocity)
	  remove_com_velocity(info);
#ifndef INFINITE_CUBE
	if(info->volume_scale_dir != -1)
	  vscale(info);
#endif
}

void
update_info(BrennerMainInfo *info)
{
}

BrennerMainInfo *
alloc_bren(int *kflag)
{
	BrennerMainInfo *info = (BrennerMainInfo*)malloc(sizeof(BrennerMainInfo));
	memset(info, 0, sizeof(*info));
    info->ljNeighborState = newNeighborState ();
    info->caNeighborState = newNeighborState ();
	info->temperature = 300;
	info->num_atms = 0;
	info->steps = 100000;
	info->timestep = 0.5;
	info->starttime = 0;
	info->lchk = 1;
	info->num_thermostated_atoms = 0;
	info->head[0] = 0;
	info->xmass[1] = 12.0;
	info->xmass[2] = 1.0;
	info->xmass[3] = 28.0;
	info->xmass[4] = 72.0;
	info->volume_scale_dir = -1;
	info->min_range_tree = 100;
	info->min_lj_range_tree = 100;
	initkt(info, kflag);
	init_c(info);
	PARAM(info);
	return info;
}

void
init_bren(BrennerMainInfo *info, int kflag, int tight_binding)
{
	int i;
	int nonzero_vel = 0;
	for (i=0; i < info->num_atms; i++)
	{
		if(info->atm_num[i].thermostated)
		  ++info->num_thermostated_atoms;
		if(info->atm_num[i].velocity.x != 0 || 
		   info->atm_num[i].velocity.y != 0 ||
		   info->atm_num[i].velocity.z != 0)
		  nonzero_vel = 1;
	}
	if(!nonzero_vel) initialize_atoms(info, kflag);
	info->dnla=3*info->timestep*info->timestep*info->num_thermostated_atoms;
	if(kflag == 6)
	  minimize(info);
}
