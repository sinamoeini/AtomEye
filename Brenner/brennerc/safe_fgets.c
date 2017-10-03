#include "safe_fgets.h"
/* Need next one for my_exit. */
#include "brenner.h"
char *safe_fgets (char *s, int size, FILE *stream, const char *what) {
  char *result = fgets (s, size, stream);
  if (0 == result) {
    fprintf (stderr, "Can't read %s", what);
    my_exit (-1);
  }
  return result;
}

