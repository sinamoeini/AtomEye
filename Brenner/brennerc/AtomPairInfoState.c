/*
 This code was released into the public domain by Tim Freeman on 8/23/2000.
*/

#include "AtomPairInfoState.h"

#ifndef __XMALLOC_H__
#include "xalloc.h"
#endif

#ifndef __BRENC_H__
#include "brenc.h"
#endif

/* Call this to ensure we have room for the given number of atoms.  We actually
   make room for one more, since the starting neighbor of the next atom is the
   ending neighbor of the previous atom. */
void allocateAtomPairAtom (struct AtomPairInfoState *apis, int max_atoms) {
  /* Here's where we make the one extra. */
  max_atoms++;
  if(max_atoms >= apis->atoms_allocated) {
    apis->atoms_allocated = 2 * apis->atoms_allocated;
    if (apis->atoms_allocated <= max_atoms) {
      apis->atoms_allocated = max_atoms + 1;
    }
    apis->ai_start =
      (int *) xrealloc (apis->ai_start,
                        sizeof(int) * apis->atoms_allocated);
  }
}

int numPairs (int ai, const struct AtomPairInfoState *s) {
  return s->ai_start[ai+1] - s->ai_start [ai];
}

struct AtomPairInfo *pairs (int ai, struct AtomPairInfoState *s) {
  return &s->ai_list[s->ai_start[ai]];
}

const struct AtomPairInfo *pairsconst
(int ai, const struct AtomPairInfoState *s) {
  return &s->ai_list[s->ai_start[ai]];
}

struct AtomPairInfoState *newAtomPairInfoState () {
  const int initial_allocation =
#ifdef NDEBUG
    1000
#else
    /* For debugging, it's good to ensure that reallocation happens. */
    3
#endif
    ;
  struct AtomPairInfoState *apis = xmalloc (sizeof (struct AtomPairInfoState));
  apis->ai_allocated = initial_allocation;
  apis->atoms_allocated = initial_allocation;
  apis->ai_start = xmalloc (sizeof (int) * initial_allocation);
  apis->ai_list = xmalloc (sizeof (struct AtomPairInfo) * initial_allocation);
  return apis;
}

/* Call this every time after a new neighbor is added.  AtomPairInfo_end should
   be the position of where we want to add the *next* neighbor, and we'll make
   room for it.  By using conservative high estimates for AtomPairInfo_end, we
   can call this only once per atom instead of once per neighbor.  Thus this
   routine probably doesn't need to be inline. */
void allocate_AtomPairInfo (struct AtomPairInfoState *ais,
                            const int AtomPairInfo_end) {
  if(AtomPairInfo_end >= ais->ai_allocated) {
    ais->ai_allocated = 2 * ais->ai_allocated;
    if (AtomPairInfo_end >= ais->ai_allocated) {
      ais->ai_allocated = AtomPairInfo_end + 1;
    }
    ais->ai_list =
      (struct AtomPairInfo *) xrealloc (ais->ai_list,
                                 sizeof (struct AtomPairInfo) *
                                 ais->ai_allocated);
  }
}

