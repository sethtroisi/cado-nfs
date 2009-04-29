
AC_INIT([gf2x], [0.9.5])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
AC_CANONICAL_TARGET

AC_SUBST([gf2x_lib_version], [1:0:0])

AM_INIT_AUTOMAKE

LT_INIT

AC_ARG_ENABLE([ntl-checks],
              [AS_HELP_STRING([--enable-ntl-checks],
                              [Turn on ntl checks])])
AC_ARG_ENABLE([sse2],
              [AS_HELP_STRING([--enable-sse2],
                              [Turn on sse-2 code (default is yes)])])

AM_CONDITIONAL([ENABLE_NTL_CHECKS],[test x$enable_ntl_checks = xyes])

AC_PROG_CC
AC_PROG_CXX
AC_COMPILE_WARNINGS
AC_PROG_CC_C99

# This macro is used for tuning
AM_PROG_CC_C_O

# A fallback for generic builds. Otherwise unused.
AC_CHECK_SIZEOF([unsigned long])
gf2x_wordsize=`expr 8 \* $ac_cv_sizeof_unsigned_long`
AC_SUBST([GF2X_WORDSIZE], [$gf2x_wordsize])


CHECK_SSE2_SUPPORT()

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
   hwdir=x86_64
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

for f in gf2x-thresholds.h gf2x_mul1.h gf2x_mul2.h gf2x_mul3.h gf2x_mul4.h ; do
 who=tuned
 if ! test -f $srcdir/hardware/$who/$f ; then who=$hwdir; fi
 if ! test -f $srcdir/hardware/$who/$f ; then who=generic$gf2x_wordsize; fi
 if ! test -f $srcdir/hardware/$who/$f ; then who=generic; fi
 if ! test -f $srcdir/hardware/$who/$f ; then AC_MSG_ERROR([$f not found]); fi
 AC_CONFIG_LINKS([gf2x/$f:hardware/$who/$f])
 if test "$f" = "gf2x-thresholds.h" ; then
    tuned_nbits=[`sed -n 's/^#define[ 	][ 	]*GF2X_WORDSIZE[ 	][ 	]*\([0-9][0-9]*\).*$/\1/p' $srcdir/hardware/$who/$f`]
 fi
done

if test x$tuned_nbits = x ; then
 tuned_nbits=$gf2x_wordsize
else
 if test x$ABI != xdefault ; then
  if test x$tuned_nbits != x$ABI ; then
   AC_MSG_ERROR([hardware/$hwdir/ assumes $tuned_nbits-bits unsigned longs, which conflicts with ABI=$ABI])
  fi
 else
  # At this point we haven't run the standard check.
  VERIFY_WORDSIZE([$tuned_nbits],[whether hardware/$hwdir/ is right assuming $tuned_nbits-bits unsigned longs])
 fi
fi

# This is used in tuning/Makefile.am
AM_CONDITIONAL([GF2X_32BIT_SOURCES],[test "x$tuned_nbits" = x32])
AM_CONDITIONAL([GF2X_64BIT_SOURCES],[test "x$tuned_nbits" = x64])
AM_CONDITIONAL([GF2X_SSE2_AVAILABLE],[test "x$gf2x_cv_cc_supports_sse2" != xno])

AC_CONFIG_HEADERS([gf2x/gf2x-config.h])

AC_CONFIG_FILES([version.sh Makefile tests/Makefile apps/Makefile tuning/Makefile])
AC_OUTPUT

# vim: set sw=1:
