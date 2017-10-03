#ifndef __expand_h__
#define __expand_h__
/* Set values in ti such that ti[j] is the coefficient of x**j in the expansion
   of (x+(L-1))**i.  Only works for i < 4.
*/
void expand (int L, int ti [4], int i);

#endif
