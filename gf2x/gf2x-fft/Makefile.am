AM_CFLAGS = -I$(top_srcdir) -I$(top_srcdir)/.. -I$(top_srcdir)/../utils -I$(top_builddir) -I$(top_builddir)/.. -I$(top_builddir)/../utils
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/.. -I$(top_srcdir)/../utils -I$(top_builddir) -I$(top_builddir)/.. -I$(top_builddir)/../utils

lib_LTLIBRARIES = libgf2x-fft.la libcantor128.la
libgf2x_fft_la_SOURCES = fake_fft.c
libcantor128_la_SOURCES = cantor128.c

check_PROGRAMS=cantor_bench cantor
cantor_SOURCES=cantor_main.c
cantor_LDADD=../libgf2x.la libcantor128.la
cantor_bench_SOURCES=cantor_bench.c
cantor_bench_LDADD=../libgf2x.la libcantor128.la
