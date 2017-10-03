/*
 This code was released into the public domain by Tim Freeman on 15 Aug 2000.
*/

/* This is the C-to-C version of the api. It exists so I can replace it with a
   C-to-C++ version when I do the Fungimol integration. 
   Tim Freeman 15 Aug 2000.
*/

#ifndef __API_H__
#include "api.h"
#endif

/* Need brenner.h for BrennerMainInfo. */
#ifndef __BRENNER_H__
#include "brenner.h"
#endif

/* A State is a state of the atoms in the system.  For the C to C interface,
   it's a BrennerMainInfo. */
struct State {
  BrennerMainInfo bmi;
};

static inline struct State *toState (BrennerMainInfo *bmi) {
  return (struct State *) bmi;
}

static inline struct State *fromState (BrennerMainInfo *s) {
  return (struct State *) s;
}

/* Returns the first Atom. */
static inline int firstIndex (const struct State *s) {
  return 0;
}

static inline int pastLastIndex (const struct State * s) {
  return s->bmi.num_atms;
}

static inline dvector pos (int ai, struct State * s) {
  return s->bmi.atm_num[ai].coord;
}

static inline Double posx (int ai, struct State * s) {
  return s->bmi.atm_num[ai].coord.x;
}

static inline Double posy (int ai, struct State * s) {
  return s->bmi.atm_num[ai].coord.y;
}

static inline Double posz (int ai, struct State * s) {
  return s->bmi.atm_num[ai].coord.z;
}

static inline vector vel (int ai, struct State * s) {
  return s->bmi.atm_num[ai].velocity;
}

static inline vector getForce (int ai, struct State * s) {
  return s->bmi.atm_num[ai].force;
}

static inline Float getForcex (int ai, struct State *s) {
  return s->bmi.atm_num[ai].force.x;
}

static inline Float getForcey (int ai, struct State *s) {
  return s->bmi.atm_num[ai].force.y;
}

static inline Float getForcez (int ai, struct State *s) {
  return s->bmi.atm_num[ai].force.z;
}

static inline void setForce (int ai, struct State * s, vector *v) {
  s->bmi.atm_num[ai].force = *v;
}

static inline void setForcex (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.x = f;
}

static inline void setForcey (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.y = f;
}

static inline void setForcez (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.z = f;
}

static inline void incForcex (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.x += f;
}

static inline void incForcey (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.y += f;
}

static inline void incForcez (int ai, struct State * s, Float f) {
  s->bmi.atm_num[ai].force.z += f;
}

static inline void transForcex (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.x += f;
  s->bmi.atm_num[j].force.x -= f;
}

static inline void transForcey (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.y += f;
  s->bmi.atm_num[j].force.y -= f;
}

static inline void transForcez (int i, int j, struct State * s, Float f) {
  s->bmi.atm_num[i].force.z += f;
  s->bmi.atm_num[j].force.z -= f;
}

static inline int getKtype (int ai, const struct State * s) {
  return s->bmi.atm_num[ai].ktype;
}

static inline int getNoa (int ktype, const struct State *s) {
  return s->bmi.noa[ktype];
}

#ifndef INFINITE_CUBE
static inline Float *getCube (struct State * s) {
  return s->bmi.cube;
}
#endif

