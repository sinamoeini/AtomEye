/*
 This code was released into the public domain by Tim Freeman on 11 Aug 2000.
*/

#ifndef __API_H__
#define __API_H__

/* 

This is an API that allows C code to access atoms and neighbors of atoms.  The
goal is to allow one version of the source code for the C port of Brenner's
molecular dynamics code to either become a Fungimol plugin or to be part of a
purely C program that does batch manipulation of text files without any
significant loss of efficiency.  We get the required flexibility by having
multiple implementations of this api, and we get efficiency by having lots of
inline functions and being willing to recompile the library to use it in
multiple contexts rather than just relink it.

For starters, we just have the version of the API that supports the batch
program.

*/

/* Need vector.h because some of the API functions declared below use
   vectors. */
#ifndef __VECTOR_H__
#include "vector.h"
#endif

#ifndef __BOOL_H__
#include "bool.h"
#endif

struct NeighborState;
struct State;



/* gcc doesn't optimize inline functions that return structures very well, so
   we have to provide inline functions that get the individual fields too. */












// The next three add the specified amount to the given force.  FIXME Strangely
// enough, I don't use the next three yet as of 11 Sep 2000 and may never use
// them.  All forces are between pairs of atoms, and in that case transForcex
// &c. are better.



// The next three add the specified amount to the force for the first index
// and decrement for the second index.





#ifndef INFINITE_CUBE

#endif
// There isn't much point in the following be inline, since we call them only
// a constant number of times per timestep. 
struct NeighborState *ljNeighbors (struct State *s);
struct NeighborState *caNeighbors (struct State *s);
const struct NeighborState *caNeighborsconst (const struct State *s);

#endif
