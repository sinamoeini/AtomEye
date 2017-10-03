/*
 This code was released into the public domain by Tim Freeman on 8/23/2000.
*/

#ifndef __NEIGHBORSTATE_H__
#include "NeighborState.h"
#endif

/* We can't make this abstract.  That is, we can't just have "struct
   NeighborState;" here and put the real structure declaration in
   CCNeighborState.c.  The problem is that ljguts.c manipulated this structure
   directly when adding neighbors.  Tim Freeman  4 Sep 2000. */
struct NeighborState {
  /* Put a pointer to the actual neighbor list first, since we probably access
     this field most often. */ 
  int *neighbor_list;
  /* Neighbors_allocated is the number of integers we can fit into
     nneighbor_list.*/
  int neighbors_allocated;
  /* Atoms_allocated is the number of integers we can fit into
     neighbor_start. */
  int atoms_allocated;
  /* For atom number i, neighbor_list [neighbor_start [i]] is the first
     neighbor of atom i.  neighbor_list [neighbor_start [i+1] - 1] is the last
     neighbor of atom i.  Looking at the last atom, we'll need room for one
     more value in neighbor_start than we have atoms.*/
  int *neighbor_start;
};

/* FIXME Only need to include xalloc.h here because of the xrealloc call in the
   inline allocate_neighbor below. */

#ifndef __XMALLOC_H__
#include "xalloc.h"
#endif

/* Call this every time before we add a bunch of neighbors.  lj_neighbor_end
   should be the end position of the last neighbor we might be adding.

   FIXME Try making this a non-inline subroutine and measure the performace.
   If there's no perceptible drop, then make this permanently non-inline.  Tim
   Freeman 29 Aug 2000.  */
static inline void allocate_neighbor (struct NeighborState *ljns,
                                      const int lj_neighbor_end) {
  if(lj_neighbor_end >= ljns->neighbors_allocated) {
    ljns->neighbors_allocated = 2 * ljns->neighbors_allocated;
    if (lj_neighbor_end >= ljns->neighbors_allocated) {
      ljns->neighbors_allocated = lj_neighbor_end + 1;
    }
    ljns->neighbor_list =
      (int *) xrealloc (ljns->neighbor_list,
                        sizeof (int) * ljns->neighbors_allocated);
  }
}

// Sort the neighbors so inessential differences are limited, which makes
// debugging easier.  natoms must be the total number of atoms in the
// NeighborState.
void sortNeighbors (struct NeighborState* ns, int natoms); 
