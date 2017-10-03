/*
 This code was released into the public domain by Tim Freeman on 15 Aug 2000.
*/

#ifndef __CCAPI_H__
#include "ccapi.h"
#endif

struct NeighborState *ljNeighbors (struct State * s) {
  return s->bmi.ljNeighborState;
}

struct NeighborState *caNeighbors (struct State * s) {
  return s->bmi.caNeighborState;
}

const struct NeighborState *caNeighborsconst
(const struct State * s)
{
  return s->bmi.caNeighborState;
}

