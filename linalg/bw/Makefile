
# Default flags ; can be overriden by methods specified in main.mkf,
# like:
#
# STATIC=yes --> adds -static
# GMP=<path> --> uses GMP as in path
# MPFR=<path> --> uses MPFR as in path

TOPPATH:=.

CC:=gcc
CXX:=g++
ARQCV:=ar qcv
OFLAGS:=-O2 -DNDEBUG
PFLAGS:=
IFLAGS:=-I ../.. -I ../../utils -I ../../gf2x -I ../../gf2x/cantor -I $(TOPPATH)/include -I . -I old
LFLAGS:=-L ../../utils -L ../../gf2x -L ../../gf2x/cantor

-include Makefile.local
#ifdef GMP
#GMPINCLUDE?=$(GMP)/include
#GMPLIB?=-L$(GMP)/lib
#IFLAGS+=$(GMPINCLUDE)
#LFLAGS+=$(GMPLIB)
#endif

ifdef PATCHED_GMP
IFLAGS+=-DPATCHED_GMP
endif

# No more variable changes after this include.
include $(TOPPATH)/util/main.mki

######################################################################

MKBITSTRINGS_SOURCES:=mkbitstrings.c
-include .mkbitstrings.mkf

#old/structure_automatic.h: old/structure.m4
#	m4 -Doutputfile=MAIN < $< > $@
#
#old/structure_inlines_automatic.h: old/structure.m4
#	m4 -Doutputfile=COMPANION < $< > $@
#
#old/structure_debug_automatic.h:   old/structure_debug.m4
#	m4 -Doutputfile=DEBUG < $< > $@
#
BW_PREP_SOURCES:=					\
		auxfuncs.c				\
		cxx/prep.cpp cxx/files.cpp		\
		cxx/matrix_header.cpp cxx/matrix.cpp	\
		cxx/preconditioner.cpp			\
		cxx/matrix_line.cpp
BW_PREP_LIBS:=-lgmpxx -lgmp
-include .bw-prep.mkf

BW_KRYLOV_SOURCES:=					\
		auxfuncs.c barrier.c			\
		cxx/config_file.cpp			\
		cxx/ticks.cpp				\
		cxx/krylov.cpp cxx/state.cpp		\
		cxx/threads.cpp	cxx/corners.cpp		\
		cxx/matrix.cpp cxx/matrix_line.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp	\
		cxx/preconditioner.cpp
BW_KRYLOV_LIBS:=-lgmpxx -lgmp
BW_KRYLOV_EXTRA_MODELS:=mt
BW_KRYLOV_MT_LIBS:=-lgmpxx -lgmp
-include .bw-krylov.mkf

# BW_SLAVE_OLD_SOURCES:=				\
# 		auxfuncs.c barrier.c			\
# 		old/bw_scalar.c old/old-endian.c		\
# 		old/filenames.c old/options.c		\
# 		old/tagfile.c old/timer.c		\
# 		old/variables.c old/addmul.c		\
# 		old/corners.c old/matrix.c		\
# 		old/recovery.c	old/slave.c		\
# 		old/state.c old/threaded.c		\
# 		old/bw_lvblock.c old/certificate.c
# BW_SLAVE_OLD_LIBS:=-lgmp
# BW_SLAVE_OLD_EXTRA_MODELS:=mt
# BW_SLAVE_OLD_MT_LIBS:=-lgmp -lrt
# -include .bw-slave-old.mkf
 
BW_PRINTMAGMA_SOURCES:=					\
		cxx/config_file.cpp			\
		auxfuncs.c cxx/matrix_header.cpp	\
		cxx/files.cpp				\
		cxx/matrix.cpp				\
		cxx/matrix_line.cpp cxx/printmagma.cpp
BW_PRINTMAGMA_LIBS:=-lgmpxx -lgmp
-include .bw-printmagma.mkf


# These are the sources of the new lingen code thing in C++. So far, it's
# such a mess that it doesn't even compile. It's also terribly slower
# than it should (when it actually compiles...).
# ifeq ($(BROKEN),also)
# BW_LINGEN_SOURCES:=					\
# 		auxfuncs.c				\
# 		cxx/lingen.cpp cxx/config_file.cpp	\
# 		cxx/files.cpp cxx/matrix_header.cpp	\
# 		cxx/ticks.cpp cxx/matrices_printing.cpp
# BW_LINGEN_LIBS:=-lgmpxx -lgmp
# -include .bw-lingen.mkf
# endif

BW_LINGEN_BINARY_SOURCES:=		\
		auxfuncs.c		\
		cxx/files.cpp		\
		cxx/lingen_binary.cpp	\
		cxx/matrix.cpp		\
		cxx/matrix_header.cpp		\
		cxx/matrix_line.cpp		\
		cxx/config_file.cpp	\
		cxx/fake_fft.c
BW_LINGEN_BINARY_LIBS:=-lgmpxx -lgmp -lutils -lgf2x ../../gf2x/cantor/cantor128.o
-include .bw-lingen-binary.mkf

ifeq (0,0)
BW_LINGEN2_SOURCES:=						\
		auxfuncs.c					\
		old/bw_scalar.c old/old-endian.c			\
		old/timer.c					\
		old/variables.c	old/fft_on_matrices.cpp		\
		old/ops_poly_fft.cpp				\
		old/ops_poly_schoolbook.cpp			\
		old/ops_poly_wrap.cpp				\
		old/fft_core.c					\
		old/field_common.c old/field_prime.c		\
		old/field_quad.c old/field_complex.c		\
		old/field_usage.c old/modulus_hacks.c		\
		old/e_polynomial.c old/structure.c		\
		old/structure_debug.c old/lingen2.cpp		\
		old/twisting_polynomials.c			\
		old/lingen_common.c
ifdef	PATCHED_GMP
BW_LINGEN2_SOURCES+=old/ops_poly_integerfft.cpp
endif
BW_LINGEN2_LIBS:=-lgmp
-include .bw-lingen2.mkf
endif

BW_LINGEN_OLD_SOURCES:=						\
		auxfuncs.c barrier.c				\
		old/bw_scalar.c old/old-endian.c			\
		old/timer.c			\
		old/variables.c old/field_common.c		\
		old/field_prime.c old/field_quad.c		\
		old/field_complex.c old/field_usage.c		\
		old/modulus_hacks.c old/structure.c		\
		old/structure_debug.c old/lingen.c		\
		old/lingen_common.c
BW_LINGEN_OLD_LIBS:=-lgmp
-include .bw-lingen-old.mkf


BENCH_POLMATMUL_SOURCES:=					\
		cxx/bench_polmatmul.cpp				\
		cxx/fake_fft.c
BENCH_POLMATMUL_LIBS:=-lgmp -lutils -lgf2x ../../gf2x/cantor/cantor128.o
-include .bench_polmatmul.mkf

# BW_FFT_BASICTEST_SOURCES:=					\
# 		auxfuncs.c					\
# 		old/bw_scalar.c					\
# 		old/variables.c 				\
# 		old/ops_poly_fft.cpp				\
# 		old/ops_poly_schoolbook.cpp			\
# 		old/ops_poly_wrap.cpp				\
# 		old/fft_core.c					\
# 		old/field_common.c old/field_prime.c		\
# 		old/field_quad.c old/field_complex.c		\
# 		old/field_usage.c old/modulus_hacks.c		\
# 		old/structure.c					\
# 		old/structure_debug.c				\
# 		old/fft-basictest.cpp
# ifdef	PATCHED_GMP
# BW_FFT_BASICTEST_SOURCES+=old/ops_poly_integerfft.cpp
# endif
# BW_FFT_BASICTEST_LIBS:=-lgmp
# -include .bw-fft-basictest.mkf
# 
# BW_BENCHFFT_SOURCES:=						\
# 		auxfuncs.c					\
# 		old/bw_scalar.c old/old-endian.c			\
# 		old/options.c					\
# 		old/timer.c					\
# 		old/variables.c old/fft_on_matrices.cpp		\
# 		old/ops_poly_fft.cpp				\
# 		old/ops_poly_schoolbook.cpp			\
# 		old/ops_poly_wrap.cpp			\
# 		old/fft_core.c					\
# 		old/field_common.c old/field_prime.c		\
# 		old/field_quad.c old/field_complex.c		\
# 		old/field_usage.c old/modulus_hacks.c		\
# 		old/structure.c old/structure_debug.c		\
# 		old/testfft.cpp
# ifdef	PATCHED_GMP
# BW_BENCHFFT_SOURCES+=old/ops_poly_integerfft.cpp
# endif
# BW_BENCHFFT_LIBS:=-lgmp
# -include .bw-benchfft.mkf
# 
BW_GATHER_SOURCES:=					\
		auxfuncs.c				\
		cxx/config_file.cpp			\
		cxx/gather.cpp				\
		cxx/matrix.cpp	cxx/matrix_line.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp	\
		cxx/preconditioner.cpp
BW_GATHER_LIBS:=-lgmpxx -lgmp
-include .bw-gather.mkf

BW_MUL_SOURCES:=					\
		auxfuncs.c	\
		cxx/mul.cpp				\
		cxx/matrix.cpp	cxx/matrix_line.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp
BW_MUL_LIBS:=-lgmpxx -lgmp
-include .bw-mul.mkf

BW_RANDOM_SOURCES:=			\
		auxfuncs.c		\
		cxx/matrix_header.cpp	\
		cxx/matrix_line.cpp	\
		cxx/random.cpp
BW_RANDOM_LIBS:=-lgmpxx -lgmp
-include .bw-random.mkf

BW_DUMP_SOURCES:=auxfuncs.c old/old-endian.c old/dump.c
BW_DUMP_LIBS:=
-include .bw-dump.mkf

BW_TRANSPOSE_SOURCES:=				\
		auxfuncs.c			\
		cxx/transpose.cpp cxx/files.cpp	\
		cxx/matrix_header.cpp		\
		cxx/matrix_line.cpp
BW_TRANSPOSE_LIBS:=
-include .bw-transpose.mkf

NFSFORGE_SOURCES:=				\
		auxfuncs.c			\
		cxx/nfsforge.cpp cxx/files.cpp	\
		cxx/matrix_header.cpp		\
		cxx/matrix_line.cpp
NFSFORGE_LIBS:=-lgmpxx -lgmp
-include .nfsforge.mkf

BW_BALANCE_SOURCES:=				\
		auxfuncs.c			\
		cxx/balance.cpp cxx/files.cpp	\
		cxx/matrix_header.cpp		\
		cxx/matrix_line.cpp
BW_BALANCE_LIBS:=
-include .bw-balance.mkf

BW_SECURE_SOURCES:=						\
		auxfuncs.c cxx/files.cpp cxx/matrix_header.cpp	\
		cxx/matrix_line.cpp cxx/secure.cpp
BW_SECURE_LIBS:=-lgmpxx -lgmp
-include .bw-secure.mkf

BW_PCOLUMN_SOURCES:=auxfuncs.c old/old-endian.c old/pcolumn.c
BW_PCOLUMN_LIBS:=
-include .bw-pcolumn.mkf

BW_MATRIX_WRITER_SOURCES:=auxfuncs.c old/old-endian.c old/matrix_writer.c
BW_MATRIX_WRITER_LIBS:=
-include .bw-matrix-writer.mkf

matlingen_dist:
	D=matlingen-`date '+%Y%m%d'` ;		\
	cp -a matlingen.skel $$D ;		\
	for f in `cat matlingen_files` ; do	\
		./cp-license.sh $$f $$D/ ;	\
	done ;					\
	tar cvvzf $$D.tgz $$D ;			\
	rm -rf $$D

######################################################################

$(eval $(MAKEFILE_TRAILER))
