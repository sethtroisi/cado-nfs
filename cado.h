/* Common header file for the CADO project
 
Copyright 2007, 2008, 2009, 2010 Pierrick Gaudry, Alexander Kruppa,
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

#include <gmp.h>

#include "cado_config.h"

#include "macros.h"

/* merge has some pressure on I/O, so having hex here speeds up the process a
 * bit */
/* TODO: This has to go. Some functions write to .purged, or parse it.
 * They should all be merged into purgedfile.c, and this #define (if
 * kept) would go there instead. */
#define PURGE_INT_FORMAT "%x"

#endif  /* CADO_H_ */
