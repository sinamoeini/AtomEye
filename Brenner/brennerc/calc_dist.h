/* brenner.h - Copyright (c) 1998 Zyvex LLC.
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

#ifndef __CALC_DIST_H__
#define __CALC_DIST_H__

/* You'd better define the type State before including calc_dist.h. */
#ifndef __API_H__
#include "api.h"
#endif

#ifndef __VECTOR_H__
#include "vector.h"
#endif

/* Takes a State and two atom numbers.  Puts the vector from atom j to atom i
   into rr, and returns the squared distance between them.  I tried having one
   definition for CALC_DIST with conditionals to see if rr was zero, but it
   seemed that gcc 2.95.2 didn't always optimize away the conditional test. 
*/

static inline Float CALC_DIST_NO_RR (struct State *s, int atm_i, int atm_j)
{
  Float result;
  Float d;
#ifndef INFINITE_CUBE
  const Float *const cube = getCube (s);
#endif
  d = posx (atm_i, s) - posx (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[0]*floor(d/cube[0] + 0.5);
#endif
  result = d * d;
  d = posy (atm_i, s) - posy (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[1]*floor(d/cube[1] + 0.5);
#endif
  result += d * d;
  d = posz (atm_i, s) - posz (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[2]*floor(d/cube[2] + 0.5);
#endif
  result += d * d;
  return result;
}

static inline Float CALC_DIST (struct State *s, int atm_i, int atm_j,
                               vector *rr)
{
  Float result;
  Float d;
#ifndef INFINITE_CUBE
  const Float *const cube = getCube (s);
#endif
  d = posx (atm_i, s) - posx (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[0]*floor(d/cube[0] + 0.5);
#endif
  result = d * d;
  rr->x = d;
  d = posy (atm_i, s) - posy (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[1]*floor(d/cube[1] + 0.5);
#endif
  result += d * d;
  rr->y = d;
  d = posz (atm_i, s) - posz (atm_j, s);
#ifndef INFINITE_CUBE
  d -= cube[2]*floor(d/cube[2] + 0.5);
#endif
  result += d * d;
  rr->z = d;
  return result;
}

#endif
