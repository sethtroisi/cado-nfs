
# Be aware that ``make check'' from here cannot ensure that the
# upper-level library is up-to-date.

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_builddir)

check_PROGRAMS=check-mul
check_mul_SOURCES=check-mul.c
check_mul_LDADD=../libgf2x.la

TESTS_ENVIRONMENT=$(SHELL) 2>/dev/null
TESTS=do-check-mul.sh

EXTRA_DIST=$(TESTS) check-mul.res
