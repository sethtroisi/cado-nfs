#  This file is part of the gf2x library.
#
#  Copyright 2007, 2008, 2009, 2010, 2013
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of either:
#   - If the archive contains a file named toom-gpl.c (not a trivial
#     placeholder), the GNU General Public License as published by the
#     Free Software Foundation; either version 3 of the License, or (at
#     your option) any later version.
#   - If the archive contains a file named toom-gpl.c which is a trivial
#     placeholder, the GNU Lesser General Public License as published by
#     the Free Software Foundation; either version 2.1 of the License, or
#     (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
#  
#  You should have received a copy of the GNU General Public License as
#  well as the GNU Lesser General Public License along with this program;
#  see the files COPYING and COPYING.LIB.  If not, write to the Free
#  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301, USA.


# Be aware that ``make check'' from here cannot ensure that the
# upper-level library is up-to-date.

AUTOMAKE_OPTIONS=subdir-objects
AM_CPPFLAGS = -I$(top_srcdir) -I$(srcdir)/.. -I$(top_srcdir)/fft -I$(top_srcdir)/fft/mpfq -I$(top_builddir)

SYNTHESIZED_TESTS =
check_PROGRAMS =

general_fftcheck_sources=check.c	\
		      check-and-bench-common.c		\
		      ../test-tools.c

# This must be exactly the same list of configurations as in fft/Makefile.am
# the flags WITHOUT_CANTOR_TRUNCATION, CANTOR_GM, and CANTOR_GM_TRUNCATE
# are irrelevant here, but we lazily copy the config block from
# fft/Makefile.am
CCOMMON=-DENGINE_CANTOR -DESCAPE_CANTOR_SYMBOLS
C64=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=64
C64nt=$(C64) -DWITHOUT_CANTOR_TRUNCATION
C64gm=$(C64) -DCANTOR_GM
C64gmt=$(C64gm) -DCANTOR_GM_TRUNCATE
C128=$(CCOMMON) -DCANTOR_BASE_FIELD_SIZE=128
C128nt=$(C128) -DWITHOUT_CANTOR_TRUNCATION
C128gm=$(C128) -DCANTOR_GM
C128gmt=$(C128gm) -DCANTOR_GM_TRUNCATE

# TODO: auto-generate this for all variants !
check_PROGRAMS+=check_cantor64
check_cantor64_SOURCES=$(general_fftcheck_sources)
check_cantor64_CFLAGS=$(AM_CFLAGS) $(C64)
check_cantor64_LDADD=$(top_builddir)/fft/libcantor64.la
SYNTHESIZED_TESTS+=	\
    check_cantor64_1000_1_c28bc2e3.fft		\
    check_cantor64_1000_2_77db95cc.fft		\
    check_cantor64_1000_3_6dc06b2f.fft		\
    check_cantor64_10000_1_2706ca53.fft		\
    check_cantor64_10000_2_3f32a516.fft		\
    check_cantor64_10000_3_fd7bba07.fft		\
    check_cantor64_100000_1_fca70ced.fft		\
    check_cantor64_100000_2_cd1293a5.fft		\
    check_cantor64_100000_3_91e7d307.fft

$(SYNTHESIZED_TESTS) : ; @touch $@

TESTS = $(SYNTHESIZED_TESTS)
TEST_EXTENSIONS = .fft
# We'll run the *.fft tests above as:
# 	./do-check-fft.sh -m check_cantor64_100_1_00000000.fft
FFT_LOG_COMPILER = $(srcdir)/do-check-fft.sh
AM_FFT_LOG_FLAGS=-m

EXTRA_DIST=do-check-fft.sh
