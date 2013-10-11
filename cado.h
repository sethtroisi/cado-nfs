/* Common header file for the CADO project

Copyright 2007, 2008, 2009, 2010, 2011 Pierrick Gaudry, Alexander Kruppa,
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

#ifndef CADO_H
#define CADO_H

/* The purpose of this header is to define some feature macros, which
 * tweak the behaviour of include files. Our intent is to define here in
 * a unique place the required macros exposing the functions we like to
 * have.
 *
 * It is necessary that this file appears only on top of the compilation
 * units, and before any other header. We promise to never include
 * another header file as a side-effect of this one (except
 * cado_config.h, which makes sense to include here as well).
 */

#if defined(__CYGWIN__) && defined(__STRICT_ANSI__)
/* The C library which comes with cygwin has no feature test macros. We do
   something ugly then.  */
#undef __STRICT_ANSI__
#endif

/* OpenBSD and FreeBSD expose *all* functions by default, and feature
 * macros are (apparently) used the other way around to restrict the
 * exposed interfaces.
 * FIXME: It's not entirely clear. Maybe it has been so in some version,
 * but it could also well be that I wholly misunderstood the thing.
 * OpenBSD 4.9 appears at least to grok _BSD_SOURCE as I expect (i.e.
 * _do_ expose BSD prototypes as an _addition_ to the rest).
 */
#ifdef __OpenBSD__
#define _BSD_SOURCE
#elif defined(__FreeBSD__)
/* XXX should check whether my former fear mentioned about turns out to
 * be true.
 */
#else
#define _POSIX_C_SOURCE 200112L /* strtoumax */
/* POSIX: popen/pclose with -std=c99, -pedantic or -ansi (requires
 * _POSIX_C_SOURCE==2 ?) fileno */
#define _XOPEN_SOURCE   700     /* posix_memalign lrand48 strndup */
#define _BSD_SOURCE     /* M_LN2 gethostname strdup random */
#define _ISOC99_SOURCE  /* Sometimes there's link trickery which causes fscanf to be linked in *only* when this is defined */
#ifndef __cplusplus
#define _GNU_SOURCE         /* asprintf vasprintf */
#endif
#define _DARWIN_C_SOURCE    /* asprintf ; _ANSI_SOURCE must be undefined */
#define _NETBSD_SOURCE      /* asprintf vasprintf */
#endif

#ifdef __cplusplus
#define __STDC_LIMIT_MACROS
#define __STDC_FORMAT_MACROS    /* PRIu32 in lingen_mat_types.hpp */
#endif

#if	defined(_AIX) && !defined(__cplusplus)
/* AIX seems to have weird rules. For example I don't understand exactly
 * how strndup is exposed. It seems that our choice of enabling feature
 * test macros above has the unfortunate effect of switching from the
 * default "expose-everything" mode for the libc headers to the more
 * restrictice "expose-only-what-has-been asked" mode. Putting
 * _ALL_SOURCE, which is obviously a catch-all, seems to maximize the
 * number of exported prototypes (it's clearly by lack of a better
 * solution). */
/* (for C++, _ALL_SOURCE is forcibly enabled by default) */
#define _ALL_SOURCE
#endif  /* _AIX */


#ifdef  __GNUC__
#pragma GCC diagnostic error "-W"
#pragma GCC diagnostic error "-Wall"
#endif

#include "cado_config.h"

#endif  /* CADO_H_ */
