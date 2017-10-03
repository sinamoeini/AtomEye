/*
 This code was released into the public domain by Tim Freeman on 29 Aug 2000.
*/

#ifndef __ATOM_PAIR_INFO_STATE_H__
#define __ATOM_PAIR_INFO_STATE_H__

struct AtomPairInfo;

/* This data structure holds miscellaneous information about the bonds.  It's
   similar to the CCNeighborState structureq.  The difference is that I intend
   to replace CCNeighborState with something else when this becomes a Fungimol
   plugin, but AtomPairInfoState will be permanently a part of the bond order
   force field.  Thus the force field can manipulate AtomPairInfoState
   directly, but it has to go through an API to deal with CCNeighborState.  Tim
   Freeman 28 Aug 2000 */
struct AtomPairInfoState {
  /* Put a pointer to the AtomPairInfo list first, since we probably access
     this field most often. */ 
  struct AtomPairInfo *ai_list;
  /* ai_allocated is the number of AtomPairInfo's we can fit into
     ai_list. */
  int ai_allocated;
  /* Atoms_allocated is the number of integers we can fit into
     neighbor_start. */
  int atoms_allocated;
  /* For atom number i, neighbor_list [neighbor_start [i]] is the first
     neighbor of atom i.  neighbor_list [neighbor_start [i+1] - 1] is the last
     neighbor of atom i.  Looking at the last atom, we'll need room for one
     more value in neighbor_start than we have atoms. */
  int *ai_start;
};

/* Call this every time after a new neighbor is added.  AtomPairInfo_end should
   be the position of where we want to add the *next* neighbor, and we'll make
   room for it.  By using conservative high estimates for AtomPairInfo_end, we
   can call this only once per atom instead of once per neighbor.  Thus this
   routine doesn't need to be inline. */
void allocate_AtomPairInfo (struct AtomPairInfoState *ais,
                            const int AtomPairInfo_end);
int numPairs (int ai, const struct AtomPairInfoState *s);

struct AtomPairInfo *pairs (int ai, struct AtomPairInfoState *s);

const struct AtomPairInfo *pairsconst (int ai,
                                       const struct AtomPairInfoState *s);

void allocateAtomPairAtom (struct AtomPairInfoState *ljns, int max_atoms);

struct AtomPairInfoState *newAtomPairInfoState ();

#endif

