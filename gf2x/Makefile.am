# This file is part of the gf2x library.
#
# Copyright 2007, 2008, 2009, 2010
# Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2.1 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with CADO-NFS; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA 02110-1301, USA.

ACLOCAL_AMFLAGS=-I config --install
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)


lib_LTLIBRARIES=libgf2x.la
libgf2x_la_SOURCES=gf2x.c toom.c fft.c
nobase_include_HEADERS=gf2x.h	\
		gf2x/gf2x-impl.h		\
		gf2x/gf2x-tfft.h		\
		gf2x/gf2x-small.h

# All these are generated.
nobase_nodist_include_HEADERS=	\
		gf2x/gf2x-config.h		\
		gf2x/gf2x-thresholds.h		\
		gf2x/gf2x_mul1.h		\
		gf2x/gf2x_mul2.h		\
		gf2x/gf2x_mul3.h		\
		gf2x/gf2x_mul4.h		\
		gf2x/gf2x_mul5.h		\
		gf2x/gf2x_mul6.h		\
	        gf2x/gf2x_mul7.h		\
		gf2x/gf2x_mul8.h		\
		gf2x/gf2x_mul9.h


libgf2x_la_LDFLAGS=-version-info @gf2x_lib_version@

#if ENABLE_NTL_CHECKS
#    MAYBE_APPS=apps
#endif

SUBDIRS=. tests src 

DIST_SUBDIRS=tests src 

DISTCLEANFILES=gf2x/gf2x-thresholds.h	\
	       gf2x/gf2x_mul1.h		\
	       gf2x/gf2x_mul2.h		\
	       gf2x/gf2x_mul3.h		\
	       gf2x/gf2x_mul4.h		\
	       gf2x/gf2x_mul5.h		\
	       gf2x/gf2x_mul6.h		\
	       gf2x/gf2x_mul7.h		\
	       gf2x/gf2x_mul8.h		\
	       gf2x/gf2x_mul9.h

EXTRA_DIST=
EXTRA_DIST+=version.sh.in
EXTRA_DIST+=BUGS
EXTRA_DIST+=patch-wrt-ntl-5.4.2 inject-ntl-5.4.2
EXTRA_DIST+=already_tuned
dist-hook:
	-find $(distdir) \( -name .svn -type d \) -o -name "*~" -o -name ".#*" \
                | xargs rm -rf

tune-lowlevel tune-toom tune-fft: all
	cd src; $(MAKE) $@

distclean-local:
	-find already_tuned/tuned -maxdepth 1 -name "gf2x*" | xargs -r rm -f
