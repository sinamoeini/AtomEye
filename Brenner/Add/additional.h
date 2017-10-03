/* additional.h: Ju Li */
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
