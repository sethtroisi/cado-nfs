
AM_CPPFLAGS = -I$(top_srcdir) -T$(top_builddir)

check_PROGRAMS=factor
factor_SOURCES=factor.cpp halfgcd.cpp
noinst_HEADERS=halfgcd.hpp
factor_LDADD=../libgf2x.la -lntl -lgmp

TESTS_ENVIRONMENT=$(SHELL) 2>/dev/null
TESTS=check-factor-000.sh check-factor-001.sh	\
      check-factor-002.sh check-factor-003.sh

EXTRA_DIST=$(TESTS)
