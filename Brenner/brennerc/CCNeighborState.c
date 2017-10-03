/*
 This code was released into the public domain by Tim Freeman on 8/23/2000.
*/

#include "CCNeighborState.h"

#ifndef __myassert_h__
#include "myassert.h"
#endif

/* stdlib.h is needed because it declares qsort. */
#ifndef __stdlib_h__
#include <stdlib.h>
#define __stdlib_h__
#endif

/* Call this to ensure we have room for the given number of atoms.  We actually
   make room for one more, since the starting neighbor of the next atom is the
   ending neighbor of the previous atom. */
void allocate_atom (struct NeighborState *ljns, int max_atoms) {
  /* We'll have one more neighbor than atom. */
  max_atoms++;
  if(max_atoms >= ljns->atoms_allocated) {
    ljns->atoms_allocated = 2 * ljns->atoms_allocated;
    if (ljns->atoms_allocated <= max_atoms) {
      ljns->atoms_allocated = max_atoms + 1;
    }
    ljns->neighbor_start =
      (int *) xrealloc (ljns->neighbor_start,
                        sizeof(int) * ljns->atoms_allocated);
  }
}

int numNeighbors (int ai, const struct NeighborState *s) {
  return s->neighbor_start[ai+1] - s->neighbor_start [ai];
}

const int *neighbors (int ai, const struct NeighborState *s) {
  return &s->neighbor_list[s->neighbor_start[ai]];
}

struct NeighborState *newNeighborState () {
  const int initial_allocation =
#ifdef NDEBUG
    1000
#else
    /* For debugging, it's good to ensure that reallocation happens. */
    3
#endif
    ;
  struct NeighborState *ljns = xmalloc (sizeof (struct NeighborState) *
                                        initial_allocation);
  ljns->neighbors_allocated = initial_allocation;
  ljns->atoms_allocated = initial_allocation;
  ljns->neighbor_start = xmalloc (sizeof (int) * initial_allocation);
  ljns->neighbor_list = xmalloc (sizeof (int) * initial_allocation);
  return ljns;
}

/* cmp_int is only used for qsort.  Return -1 if i1 < i2, 0 if they're equal,
   or 1 if i1 > i2.  i1 and i2 are really int's.  We're taking advantage of the
   fact that ints and pointers are the same size at the moment. */
static int cmp_int (const void *i1, const void *i2) {
  const int n1 = *(const int *) i1;
  const int n2 = *(const int *) i2;
  return (n1>n2)-(n1<n2);
}

void sortNeighbors (struct NeighborState *const ns, const int natoms) {
  int i;
  assert (sizeof(int) == sizeof (const void *));
  for (i = 0; i < natoms; i++) {
    qsort (&(ns->neighbor_list[ns->neighbor_start[i]]),
           ns->neighbor_start[i+1] - ns->neighbor_start[i],
           sizeof (int),
           cmp_int);
  }
}
