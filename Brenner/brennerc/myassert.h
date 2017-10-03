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

#ifndef __myassert_h__
#define __myassert_h__

#if 1
   int myassertimplementation
   (const char *string, const char *file, const int line,
    const char *pretty_function);
   // GCC's assert in assert.h garbages the caller's stack frame when it
   // triggers, making debugging difficult.
   // If they included <assert.h>, get rid of the old assert macro.
   #ifdef assert
      #error You included the systems assert.h!
   #endif
   #undef assert
   // Next line should cause the external include guards not to
   // redundantly include <assert.h>  
   #define __assert_h__
   #ifdef NDEBUG
      #define assert(x) ((void) 0)
   #else
      #define assert(x) ((void) ((x)?0: \
		      myassertimplementation(#x, __FILE__, __LINE__, \
			       __PRETTY_FUNCTION__)))
   #endif
#else
   #ifndef __assert_h__
   #include <assert.h>
   #define __assert_h__
   #endif
#endif
#endif
