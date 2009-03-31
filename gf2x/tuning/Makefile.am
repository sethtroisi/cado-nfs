
# This impacts the storage class of the functions in the inlines files.
AM_CPPFLAGS = -I$(top_builddir) -I$(top_srcdir) -DTUNING
EXTRA_DIST=

# We've got a lot of per-target cfiles here, but it's rather silly to
# have automake produce separate .o files for the files which fo _not_
# change in this regard. So we resort to having lots of silly libraries,
# as a workaround.

EXTRA_LTLIBRARIES=libtiming.la
libtiming_la_SOURCES=timing.c

EXTRA_LTLIBRARIES+=libtuneup-s1.la
libtuneup_s1_la_SOURCES=tuneup.c
libtuneup_s1_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=1
EXTRA_LTLIBRARIES+=libtuneup-s2.la
libtuneup_s2_la_SOURCES=tuneup.c
libtuneup_s2_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=2
EXTRA_LTLIBRARIES+=libtuneup-s3.la
libtuneup_s3_la_SOURCES=tuneup.c
libtuneup_s3_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=3
EXTRA_LTLIBRARIES+=libtuneup-s4.la
libtuneup_s4_la_SOURCES=tuneup.c
libtuneup_s4_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=4
#EXTRA_LTLIBRARIES+=libtuneup-s5.la
#libtuneup_s5_la_SOURCES=tuneup.c
#libtuneup_s5_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=5
#EXTRA_LTLIBRARIES+=libtuneup-s6.la
#libtuneup_s6_la_SOURCES=tuneup.c
#libtuneup_s6_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=6
#EXTRA_LTLIBRARIES+=libtuneup-s7.la
#libtuneup_s7_la_SOURCES=tuneup.c
#libtuneup_s7_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=7
#EXTRA_LTLIBRARIES+=libtuneup-s8.la
#libtuneup_s8_la_SOURCES=tuneup.c
#libtuneup_s8_la_CPPFLAGS=-I$(top_builddir) -I$(top_srcdir) -DTUNING -DTSIZE=8

# noinst_PROGRAMS=tuneup_mul1 tuneup_mul2 tuneup_mul4

EXTRA_PROGRAMS=gen_bb_mul_code$(EXEEXT)
gen_bb_mul_code_SOURCES=gen_bb_mul_code.c

LOWLEVEL=
CLEANFILES=

# The list below is the result of variations around the following vim macro:
# s/.*/EXTRA_PROGRAMS+=tune_&^Mtune_&_SOURCES=tuneup.c
# &.c^Mtune_&_CPPFLAGS=-DTSIZE=&/g
#
#  which is fed the list of all functions to be tested.
######################################################################
if GF2X_32BIT_SOURCES
LOWLEVEL+=tune_gen_mul1_w32k1$(EXEEXT)
nodist_tune_gen_mul1_w32k1_SOURCES=gen_mul1_w32k1.c
tune_gen_mul1_w32k1_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k1.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 1 > $@
CLEANFILES+=gen_mul1_w32k1.c

LOWLEVEL+=tune_gen_mul1_w32k2$(EXEEXT)
nodist_tune_gen_mul1_w32k2_SOURCES=gen_mul1_w32k2.c
tune_gen_mul1_w32k2_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k2.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 2 > $@
CLEANFILES+=gen_mul1_w32k2.c

LOWLEVEL+=tune_gen_mul1_w32k3$(EXEEXT)
nodist_tune_gen_mul1_w32k3_SOURCES=gen_mul1_w32k3.c
tune_gen_mul1_w32k3_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k3.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 3 > $@
CLEANFILES+=gen_mul1_w32k3.c

LOWLEVEL+=tune_gen_mul1_w32k4$(EXEEXT)
nodist_tune_gen_mul1_w32k4_SOURCES=gen_mul1_w32k4.c
tune_gen_mul1_w32k4_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k4.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 4 > $@
CLEANFILES+=gen_mul1_w32k4.c

LOWLEVEL+=tune_gen_mul1_w32k5$(EXEEXT)
nodist_tune_gen_mul1_w32k5_SOURCES=gen_mul1_w32k5.c
tune_gen_mul1_w32k5_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k5.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 5 > $@
CLEANFILES+=gen_mul1_w32k5.c

LOWLEVEL+=tune_gen_mul1_w32k6$(EXEEXT)
nodist_tune_gen_mul1_w32k6_SOURCES=gen_mul1_w32k6.c
tune_gen_mul1_w32k6_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k6.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 6 > $@
CLEANFILES+=gen_mul1_w32k6.c

LOWLEVEL+=tune_gen_mul1_w32k7$(EXEEXT)
nodist_tune_gen_mul1_w32k7_SOURCES=gen_mul1_w32k7.c
tune_gen_mul1_w32k7_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k7.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 7 > $@
CLEANFILES+=gen_mul1_w32k7.c

LOWLEVEL+=tune_gen_mul1_w32k8$(EXEEXT)
nodist_tune_gen_mul1_w32k8_SOURCES=gen_mul1_w32k8.c
tune_gen_mul1_w32k8_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w32k8.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 32 8 > $@
CLEANFILES+=gen_mul1_w32k8.c
endif

######################################################################

if GF2X_64BIT_SOURCES
LOWLEVEL+=tune_gen_mul1_w64k1$(EXEEXT)
nodist_tune_gen_mul1_w64k1_SOURCES=gen_mul1_w64k1.c
tune_gen_mul1_w64k1_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k1.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 1 > $@
CLEANFILES+=gen_mul1_w64k1.c

LOWLEVEL+=tune_gen_mul1_w64k2$(EXEEXT)
nodist_tune_gen_mul1_w64k2_SOURCES=gen_mul1_w64k2.c
tune_gen_mul1_w64k2_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k2.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 2 > $@
CLEANFILES+=gen_mul1_w64k2.c

LOWLEVEL+=tune_gen_mul1_w64k3$(EXEEXT)
nodist_tune_gen_mul1_w64k3_SOURCES=gen_mul1_w64k3.c
tune_gen_mul1_w64k3_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k3.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 3 > $@
CLEANFILES+=gen_mul1_w64k3.c

LOWLEVEL+=tune_gen_mul1_w64k4$(EXEEXT)
nodist_tune_gen_mul1_w64k4_SOURCES=gen_mul1_w64k4.c
tune_gen_mul1_w64k4_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k4.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 4 > $@
CLEANFILES+=gen_mul1_w64k4.c

LOWLEVEL+=tune_gen_mul1_w64k5$(EXEEXT)
nodist_tune_gen_mul1_w64k5_SOURCES=gen_mul1_w64k5.c
tune_gen_mul1_w64k5_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k5.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 5 > $@
CLEANFILES+=gen_mul1_w64k5.c

LOWLEVEL+=tune_gen_mul1_w64k6$(EXEEXT)
nodist_tune_gen_mul1_w64k6_SOURCES=gen_mul1_w64k6.c
tune_gen_mul1_w64k6_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k6.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 6 > $@
CLEANFILES+=gen_mul1_w64k6.c

LOWLEVEL+=tune_gen_mul1_w64k7$(EXEEXT)
nodist_tune_gen_mul1_w64k7_SOURCES=gen_mul1_w64k7.c
tune_gen_mul1_w64k7_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k7.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 7 > $@
CLEANFILES+=gen_mul1_w64k7.c

LOWLEVEL+=tune_gen_mul1_w64k8$(EXEEXT)
nodist_tune_gen_mul1_w64k8_SOURCES=gen_mul1_w64k8.c
tune_gen_mul1_w64k8_LDADD=libtuneup-s1.la libtiming.la ../libgf2x.la
gen_mul1_w64k8.c: gen_bb_mul_code$(EXEEXT) ; ./gen_bb_mul_code 64 8 > $@
CLEANFILES+=gen_mul1_w64k8.c
endif

######################################################################
LOWLEVEL+=tune_mul2k$(EXEEXT)
tune_mul2k_SOURCES=mul2k.c
tune_mul2k_LDADD=libtuneup-s2.la libtiming.la ../libgf2x.la

LOWLEVEL+=tune_mul3k$(EXEEXT)
tune_mul3k_SOURCES=mul3k.c
tune_mul3k_LDADD=libtuneup-s3.la libtiming.la ../libgf2x.la

LOWLEVEL+=tune_mul4k$(EXEEXT)
tune_mul4k_SOURCES=mul4k.c
tune_mul4k_LDADD=libtuneup-s4.la libtiming.la ../libgf2x.la

######################################################################
if GF2X_SSE2_AVAILABLE
if GF2X_64BIT_SOURCES
LOWLEVEL+=tune_mul2t$(EXEEXT)
tune_mul2t_SOURCES=mul2t.c
tune_mul2t_LDADD=libtuneup-s2.la libtiming.la ../libgf2x.la
endif
if GF2X_32BIT_SOURCES
LOWLEVEL+=tune_mul4t$(EXEEXT)
tune_mul4t_SOURCES=mul4t.c
tune_mul4t_LDADD=libtuneup-s4.la libtiming.la ../libgf2x.la

LOWLEVEL+=tune_mul3t$(EXEEXT)
tune_mul3t_SOURCES=mul3t.c
tune_mul3t_LDADD=libtuneup-s3.la libtiming.la ../libgf2x.la
endif
endif

######################################################################
EXTRA_PROGRAMS+=$(LOWLEVEL)
tune-lowlevel: $(LOWLEVEL)
	cd $(top_builddir) ; $(MAKE)
	MAKE="$(MAKE)" srcdir=$(srcdir)	\
	perl $(srcdir)/tune-lowlevel.pl $(LOWLEVEL)
EXTRA_DIST+=tune-lowlevel.pl

######################################################################
noinst_HEADERS=timing.h tuning-common.h replace.h

EXTRA_PROGRAMS+=update-thresholds$(EXEEXT)
update_thresholds_SOURCES=update-thresholds.c replace.c

######################################################################
EXTRA_PROGRAMS+=tunetoom$(EXEEXT)
tunetoom_SOURCES=tunetoom.c tuning-common.c
tunetoom_LDADD=libtiming.la ../libgf2x.la

TOOM_TUNING_LIMIT=2048

tunetoom.res: tunetoom$(EXEEXT)
	./tunetoom -s 1.05 $(TOOM_TUNING_LIMIT) -o tunetoom.res

tune-toom: tunetoom.res update-thresholds$(EXEEXT)
	cd $(top_builddir) ; $(MAKE)
	cat ../gf2x/gf2x-thresholds.h > tuned_thresholds.h
	./update-thresholds -o tuned_thresholds.h < tunetoom.res
	mkdir -p ../hardware/tuned/
	mv -f tuned_thresholds.h ../hardware/tuned/gf2x-thresholds.h
	rm -f ../gf2x/gf2x-thresholds.h
	ln -sf ../hardware/tuned/gf2x-thresholds.h ../gf2x/
	cd $(top_builddir) ; $(MAKE)

######################################################################
EXTRA_PROGRAMS+=tunefft$(EXEEXT)
tunefft_SOURCES=tunefft.c tuning-common.c
tunefft_LDADD=libtiming.la ../libgf2x.la

FFT_TUNING_LIMIT=8000000

tunefft.res: tunefft$(EXEEXT)
	./tunefft -s 1.02 $(FFT_TUNING_LIMIT) -o tunefft.res

tune-fft: tunefft.res update-thresholds$(EXEEXT)
	cd $(top_builddir) ; $(MAKE)
	cat ../gf2x/gf2x-thresholds.h > tuned_thresholds.h
	./update-thresholds -o tuned_thresholds.h < tunefft.res
	mkdir -p ../hardware/tuned/
	mv -f tuned_thresholds.h ../hardware/tuned/gf2x-thresholds.h
	rm -f ../gf2x/gf2x-thresholds.h
	ln -sf ../hardware/tuned/gf2x-thresholds.h ../gf2x/
	cd $(top_builddir) ; $(MAKE)

CLEANFILES+=$(EXTRA_PROGRAMS) $(EXTRA_LTLIBRARIES) tunetoom.res tunefft.res
