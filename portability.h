/* Portability header file for the CADO project
 
Copyright 2013 Pierrick Gaudry, Alexander Kruppa,
               Emmanuel Thome, Paul Zimmermann

This file is part of the CADO project.

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

*/

/* This header file defines macros (and perhaps static functions) to improve 
   portability of the CADO code. They aim to provide a wrapper for some C99
   and POSIX functionality for systems that lack those. */

#ifndef CADO_PORTABILITY_H_
#define CADO_PORTABILITY_H_

#ifndef CADO_VERSION_MAJOR
#error cado_config.h must be included before portability.h
#endif

#ifndef HAVE_STRDUP
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline char * strdup(const char * const s)
{
    const size_t size = strlen(s) + 1;
    char * const r = (char *) malloc(size * sizeof(char));
    if (r != NULL)
        memcpy(r, s, size);
    return r;
}
#ifdef __cplusplus
}
#endif
#endif /* HAVE_STRDUP */

#ifndef HAVE_STRNDUP
/* Not every libc has this, and providing a workalike is very easy */
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline char * strndup(const char * const a, const size_t n)
{
    const size_t l = strlen(a);
    const size_t size = (l < n ? l : n) + 1;
    char * const r = (char *) malloc(size * sizeof(char));
    if (r != NULL) {
        memcpy(r, a, size);
        r[size] = '\0';
    }
    return r;
}
#ifdef __cplusplus
}
#endif
#endif /* HAVE_STRNDUP */

/* MS VS and MinGW use the MS RTL (called MSVCRT for MinGW) which does not
   know the "%zu" format, they use "%Iu" instead. On MinGW, we use wrapper 
   functions that rewrite the %zu format accordingly, so the bulk of the
   code can continue to use C99 syntax.
   We do these renames only if stdio.h has been parsed before this file.
   Header files that need a certain include order are ugly, but that never
   stopped us and renaming printf() before parsing stdio.h would be "bad." 
   This way, when the renames are needed but don't happen, with any luck 
   gcc will complain about not understanding "%zu". Note that C++ does not 
   define _STDIO_H, so we test instead for a constant that all stdio.h 
   headers should define, SEEK_SET, and hope that this is in fact a 
   preprocessor macro (at least under MinGW which is the case that matters).
   If NO_PRINTF_RENAME is defined, no renames happen. This is meant to allow 
   the code that implements the format substitutions to refer to the plain
   libc functions. */

#ifdef HAVE_MINGW

#ifndef SEEK_SET
#error stdio.h must be included before portability.h
#endif

#include <stdarg.h>
#include <stdlib.h>
#include "macros.h"

static inline const char *
subst_zu(const char * const s)
{
    char * const r = strdup(s);
    const char * const prisiz = "Iu";
    const char * const priptrdiff = "Id";
    const size_t l = strlen(r);
    size_t i;
    
    ASSERT_ALWAYS(strlen(prisiz) == 2);
    ASSERT_ALWAYS(strlen(priptrdiff) == 2);
    ASSERT_ALWAYS(r != NULL);
    for (i = 0; i + 2 < l; i++)
        if (r[i] == '%' && r[i+1] == 'z' && r[i+2] == 'u') {
            r[i+1] = prisiz[0];
            r[i+2] = prisiz[1];
        } else if (r[i] == '%' && r[i+1] == 't' && r[i+2] == 'd') {
            r[i+1] = priptrdiff[0];
            r[i+2] = priptrdiff[1];
        }
    return r;
}

static inline int
scanf_subst_zu (const char * const format, ...)
{
  va_list ap;
  const char * const subst_format = subst_zu (format);
  int r;
  
  va_start (ap, format);
  r = vscanf (subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

static inline int
fscanf_subst_zu (FILE * const stream, const char * const format, ...)
{
  va_list ap;
  const char * const subst_format = subst_zu (format);
  int r;
  
  va_start (ap, format);
  r = vfscanf (stream, subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

static inline int
sscanf_subst_zu (const char * const str, const char * const format, ...)
{
  va_list ap;
  const char * const subst_format = subst_zu (format);
  int r;
  
  va_start (ap, format);
  r = vsscanf (str, subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

#define scanf scanf_subst_zu
#define fscanf fscanf_subst_zu
#define sscanf sscanf_subst_zu

#endif /* ifdef HAVE_MINGW */

#ifndef HAVE_ASPRINTF
/* Copied and improved from
 * http://mingw-users.1079350.n2.nabble.com/Query-regarding-offered-alternative-to-asprintf-td6329481.html
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __cplusplus
extern "C" {
#endif
static inline int 
vasprintf( char ** const sptr, const char *const fmt, va_list argv )
{
    int wanted = vsnprintf( *sptr = NULL, 0, fmt, argv );
    if (wanted<0)
        return -1;
    *sptr = (char *) malloc(1 + wanted);
    if (!*sptr)
        return -1;
    int rc = vsnprintf(*sptr, 1+wanted, fmt, argv );
    return rc;
}

static inline int 
asprintf( char ** const sptr, const char * const fmt, ... )
{
    int retval;
    va_list argv;
    va_start(argv, fmt);
    retval = vasprintf(sptr, fmt, argv);
    va_end(argv);
    return retval;
}
#ifdef __cplusplus
}
#endif
#endif  /* HAVE_ASPRINTF */

#ifndef HAVE_GETC_UNLOCKED
#define getc_unlocked getc
#endif
#ifndef HAVE_LRAND48
#define lrand48 rand
#endif

#if defined(HAVE_MINGW) && !defined(HAVE_REALPATH)
#include <io.h>
#include <stdlib.h>
#include <errno.h>
/* realpath() function copied from 
 * http://sourceforge.net/p/mingw/patches/256/?page=0
 * Its copyright notice is:
 * Written by Keith Marshall <keithmarshall@users.sourceforge.net>
 *
 * This is free software.  You may redistribute and/or modify it as you
 * see fit, without restriction of copyright. */
    
static inline char __cdecl
*realpath( const char *__restrict__ name, char *__restrict__ resolved )
{
  char *retname = NULL;

  if( name == NULL )
    errno = EINVAL;
  else if( access( name, 4 ) == 0 )
  {
    if( (retname = resolved) == NULL )
    {
      retname = (char *) malloc( _MAX_PATH );
    }
    if( retname == NULL )
      errno = ENOMEM;
    else if( (retname = _fullpath( retname, name, _MAX_PATH )) == NULL )
      errno = ENAMETOOLONG;
  }
  return retname;
}
#endif

#endif /* ifndef CADO_PORTABILITY_H_ */
