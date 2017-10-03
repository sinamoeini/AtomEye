#ifndef __SAFE_FGETS_H__
#define __SAFE_FGETS_H__
/* Next one requires for FILE. */
#include <stdio.h>
char *safe_fgets (char *s, int size, FILE *stream, const char *what);
#endif
