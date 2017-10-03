/*
 This code was released into the public domain by Tim Freeman on 20 Aug 2000.

 Memory allocation with error checking.
*/

#include "xalloc.h"

/* Next one for fprintf, stderr. */
#include "stdio.h"
/* Next one for exit. */
#include <unistd.h>

/* Next one for exit. */
#include <unistd.h>

void *xrealloc (void *ptr, size_t size) {
  void *result = realloc (ptr, size);
  if (0 == result && 0 != size) {
    fprintf (stderr, "Failed to reallocate to get %d bytes.\n", size);
    exit (1);
  }
  return result;
}

void *xmalloc (size_t size) {
  void *result = malloc (size);
  if (0 == result) {
    fprintf (stderr, "Failed to allocate %d bytes.\n", size);
    exit (1);
  }
  return result;
}
