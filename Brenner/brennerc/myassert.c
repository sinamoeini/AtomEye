// Fungimol - an extensible system for designing atomic-scale objects.
// Copyright (C) 2000 Tim Freeman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
// 
// You should have received a copy of the GNU Library General Public
// License along with this library in the file COPYING.txt; if not,
// write to the Free Software Foundation, Inc., 59 Temple Place -
// Suite 330, Boston, MA 02111-1307, USA
//
// The author can be reached by email at tim@infoscreen.com, or by
// paper mail at:
//
// Tim Freeman
// 655 S. FairOaks Ave., Apt B-316
// Sunnyvale, CA 94086
//

#include "myassert.h"

// For abort.
#ifndef __stdlib_h__
#include <stdlib.h>
#define __stdlib_h__
#endif

// For printf.
#ifndef __stdio_h__
#include <stdio.h>
#define __stdio_h__
#endif

int myassertimplementation
(const char *string, const char *file, const int line,
 const char *pretty_function) {
  fprintf (stderr, "Assertion failed at line %d of %s in %s:\n%s\n",
           line, file, pretty_function, string);
  abort ();
  return 1;
}
