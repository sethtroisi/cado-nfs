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
   portability of the CADO code. They aim to provide a wrappers for some C99 
   and POSIX functionality for systems that lack those. */

#ifndef CADO_PORTABILITY_H_
#define CADO_PORTABILITY_H_

#ifndef CADO_VERSION_MAJOR
#error cado_config.h must be included before portability.h
#endif

#ifndef HAVE_STRDUP
#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
#include <string.h>
static inline char *
strdup(const char * const s)
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
#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>
#include <string.h>
static inline char *
strndup(const char * const a, const size_t n)
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

#ifndef HAVE_ASPRINTF
/* Copied and improved from
 * http://mingw-users.1079350.n2.nabble.com/Query-regarding-offered-alternative-to-asprintf-td6329481.html
 */
#ifdef __cplusplus
extern "C" {
#endif
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
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

#endif /* ifndef CADO_PORTABILITY_H_ */
