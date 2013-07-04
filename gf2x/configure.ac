#  This file is part of the gf2x library.
# 
#  Copyright 2007, 2008, 2009, 2010
#  Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann
# 
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 2.1 of the License, or (at
#  your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#  License for more details.
#  
#  You should have received a copy of the GNU Lesser General Public
#  License along with CADO-NFS; see the file COPYING.  If not, write to
#  the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
#  Boston, MA 02110-1301, USA.


AC_INIT([gf2x], [1.0])
: ${CFLAGS=-g -O4}
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CANONICAL_TARGET

AC_SUBST([gf2x_lib_version], [1:0:0])

AM_INIT_AUTOMAKE

AM_MAINTAINER_MODE([disable])

LT_INIT()

AC_ARG_ENABLE([ntl-checks],
              [AS_HELP_STRING([--enable-ntl-checks],
                              [Turn on ntl checks])])
AC_ARG_ENABLE([sse2],
              [AS_HELP_STRING([--enable-sse2],
                              [Turn on sse-2 code (default is yes)])])

AC_ARG_ENABLE([pclmul],
              [AS_HELP_STRING([--enable-pclmul],
                              [Turn on pclmul code (default is yes)])])

AM_CONDITIONAL([ENABLE_NTL_CHECKS],[test x$enable_ntl_checks = xyes])

AC_PROG_CC
AC_PROG_CXX
AC_COMPILE_WARNINGS
AC_PROG_CC_C99

GF2X_PROG_CC_FOR_BUILD
GF2X_PROG_EXEEXT_FOR_BUILD

# This macro is used for tuning
AM_PROG_CC_C_O

# A fallback for generic builds. Otherwise unused.
AC_CHECK_SIZEOF([unsigned long])
gf2x_wordsize=`expr 8 \* $ac_cv_sizeof_unsigned_long`
AC_SUBST([GF2X_WORDSIZE], [$gf2x_wordsize])


CHECK_SSE2_SUPPORT()

CHECK_PCLMUL_SUPPORT()

: ${ABI='default'}

if test x$ABI != xdefault ; then
 VERIFY_WORDSIZE([$ABI],[whether $CC and ABI=$ABI fit together])
fi

if test x$hwdir = x ; then
 # Our _default_ setting for hwdir is to use 8*sizeof(ulong).
 hwdir=generic$gf2x_wordsize

 # Yet there are several situations in which we select another directory.
 # This refined selection mechanism is impacted by the ABI= setting.
 case "$target_cpu" in
  # All 32-bit pentium patterns
  i?86|pentium*|athlon|prescott)
   if test "$gf2x_cv_cc_supports_sse2" = "no" ; then
    hwdir=x86_nosse2
   else
    hwdir=x86_sse2
    # Maybe on some funny platform from hell, we could get i386 yet have a
    # 64-bit cpu ?
    if test x$ABI = x64 ; then
     hwdir=x86_64
    fi
   fi
   ;;
  core2|opteron|x86_64|nocona|k10)
   # Note: could set to something else if we have separate tuning files.
   if test "$gf2x_cv_cc_supports_pclmul" = "no" ; then
     hwdir=x86_64
   else
     hwdir=x86_64_pclmul
   fi
   if test x$ABI = x32 ; then
    hwdir=x86_sse2
   elif test x$ABI = xdefault ; then
    RUNTIME_ULONG_BITS()
    case x$gf2x_cv_ulongbits in
     x64|xdontknow) true ;;
     *)
       AC_MSG_WARN([The default ABI for this compiler has $gf2x_cv_ulongbits-bit unsigned longs, not 64-bit as the CPU supports. Using $gf2x_cv_ulongbits-bit])
       hwdir=x86_sse2;;
    esac
   fi
   ;;
 esac
else
 AC_MSG_NOTICE([Using supplied hwdir value $hwdir])
fi

AC_MSG_NOTICE([using ABI="$ABI"])
AC_MSG_NOTICE([      CC="$CC"])
AC_MSG_NOTICE([      CFLAGS="$CFLAGS"])
AC_MSG_NOTICE([      CPPFLAGS="$CPPFLAGS"])
AC_MSG_NOTICE([      hwdir="$hwdir"])

for f in gf2x-thresholds.h      \
        gf2x_mul1.h gf2x_mul2.h gf2x_mul3.h gf2x_mul4.h \
        gf2x_mul5.h gf2x_mul6.h gf2x_mul7.h gf2x_mul8.h \
        gf2x_mul9.h ; do
 who=tuned
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=$hwdir; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=generic$gf2x_wordsize; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then who=generic; fi
 if ! test -f $srcdir/already_tuned/$who/$f ; then AC_MSG_ERROR([$f not found]); fi
 AC_CONFIG_LINKS([gf2x/$f:already_tuned/$who/$f])
 if test "$f" = "gf2x-thresholds.h" ; then
    tuned_nbits=[`sed -n 's/^#define[ 	][ 	]*GF2X_WORDSIZE[ 	][ 	]*\([0-9][0-9]*\).*$/\1/p' $srcdir/already_tuned/$who/$f`]
 fi
done

if test x$tuned_nbits = x ; then
 tuned_nbits=$gf2x_wordsize
else
 if test x$ABI != xdefault ; then
  if test x$tuned_nbits != x$ABI ; then
   AC_MSG_ERROR([already_tuned/$hwdir/ assumes $tuned_nbits-bits unsigned longs, which conflicts with ABI=$ABI])
  fi
 else
  # At this point we haven't run the standard check.
  VERIFY_WORDSIZE([$tuned_nbits],[whether already_tuned/$hwdir/ is right assuming $tuned_nbits-bits unsigned longs])
 fi
fi

# This is used in src/Makefile.am
AM_CONDITIONAL([GF2X_32BIT_SOURCES],[test "x$tuned_nbits" = x32])
AM_CONDITIONAL([GF2X_64BIT_SOURCES],[test "x$tuned_nbits" = x64])
AM_CONDITIONAL([GF2X_SSE2_AVAILABLE],[test "x$gf2x_cv_cc_supports_sse2" != xno])
AM_CONDITIONAL([GF2X_PCLMUL_AVAILABLE],[test "x$gf2x_cv_cc_supports_pclmul" != xno])

AC_CONFIG_HEADERS([gf2x/gf2x-config.h])

AC_CONFIG_FILES([version.sh Makefile tests/Makefile src/Makefile gf2x-fft/Makefile])
AC_OUTPUT

# vim: set sw=1:
