AM_CFLAGS = -I$(top_srcdir) -I$(top_srcdir)/.. -I$(top_srcdir)/../utils -I$(top_builddir) -I$(top_builddir)/.. -I$(top_builddir)/../utils
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/.. -I$(top_srcdir)/../utils -I$(top_builddir) -I$(top_builddir)/.. -I$(top_builddir)/../utils

lib_LTLIBRARIES = libgf2x-fft.la libcantor128.la


check_LTLIBRARIES =libcantor128gm.la libcantor128gmt.la
libgf2x_fft_la_SOURCES = fake_fft.c
libcantor128_la_SOURCES = cantor128.c
libcantor128gm_la_SOURCES = cantor128.c
libcantor128gmt_la_SOURCES = cantor128.c
libcantor128gm_la_CFLAGS = $(AM_CFLAGS) -DCANTOR_GM
libcantor128gmt_la_CFLAGS = $(AM_CFLAGS) -DCANTOR_GM -DCANTOR_GM_TRUNCATE

# check_PROGRAMS=cantor cantor_gm cantor_gmt
# cantor_SOURCES=cantor_main.c
# cantor_LDADD=../libgf2x.la libcantor128.la -lgmp
# cantor_gm_SOURCES=cantor_main.c
# cantor_gm_LDADD=../libgf2x.la libcantor128gm.la -lgmp
# cantor_gmt_SOURCES=cantor_main.c
# cantor_gmt_LDADD=../libgf2x.la libcantor128gmt.la -lgmp
 
#check_PROGRAMS=cantor_bench cantor_gm_bench cantor_gmt_bench cantor cantor_gm cantor_gmt
#cantor_bench_SOURCES=cantor_bench.c
#cantor_bench_LDADD=../libgf2x.la libcantor128.la -lgmp
#cantor_gm_bench_SOURCES=cantor_bench.c
#cantor_gm_bench_LDADD=../libgf2x.la libcantor128gm.la -lgmp
#cantor_gmt_bench_SOURCES=cantor_bench.c
#cantor_gmt_bench_LDADD=../libgf2x.la libcantor128gmt.la -lgmp
