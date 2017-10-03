/*
 This code was released into the public domain by Tim Freeman on 20 Aug 2000.

 Memory allocation with error checking.
*/

#ifndef __XALLOC_H__
#define __XALLOC_H__
/* Next one is for size_t. */
#include <stdlib.h>

/* Stop if reallocation failed. */
void *xrealloc (void *ptr, size_t size);
/* Stop if allocation failed. */
void *xmalloc (size_t size);

#endif
