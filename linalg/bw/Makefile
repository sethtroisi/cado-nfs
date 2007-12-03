
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
OFLAGS:=
PFLAGS:=
IFLAGS:=-iquote $(TOPPATH)/include -iquote . -iquote old
LFLAGS:=

# No more variable changes after this include.
include $(TOPPATH)/util/main.mki

######################################################################


old/structure_automatic.h: old/structure.m4
	m4 -Doutputfile=MAIN < $< > $@

old/structure_inlines_automatic.h: old/structure.m4
	m4 -Doutputfile=COMPANION < $< > $@

old/structure_debug_automatic.h:   old/structure_debug.m4
	m4 -Doutputfile=DEBUG < $< > $@

BW_PREP_SOURCES:=					\
		auxfuncs.c				\
		cxx/prep.cpp cxx/files.cpp		\
		cxx/matrix_header.cpp
BW_PREP_LIBS:=-lgmp -lgmpxx
-include .bw-prep.mkf

BW_SLAVE_SOURCES:=					\
		auxfuncs.c barrier.c			\
		cxx/detect_params.cpp			\
		cxx/ticks.cpp				\
		cxx/slave.cpp cxx/state.cpp		\
		cxx/threads.cpp	cxx/corners.cpp		\
		cxx/matrix.cpp cxx/matrix_line.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp
BW_SLAVE_LIBS:=-lgmp -lgmpxx
BW_SLAVE_EXTRA_MODELS:=mt
BW_SLAVE_MT_LIBS:=-lgmp -lgmpxx
-include .bw-slave.mkf

BW_SLAVE_OLD_SOURCES:=					\
		auxfuncs.c barrier.c			\
		old/bw_scalar.c old/endian.c		\
		old/filenames.c old/options.c		\
		old/tagfile.c old/timer.c		\
		old/variables.c old/addmul.c		\
		old/corners.c old/matrix.c		\
		old/recovery.c	old/slave.c		\
		old/state.c old/threaded.c		\
		old/bw_lvblock.c old/certificate.c
BW_SLAVE_OLD_LIBS:=-lgmp
BW_SLAVE_OLD_EXTRA_MODELS:=mt
BW_SLAVE_OLD_MT_LIBS:=-lgmp -lrt
-include .bw-slave-old.mkf

BW_PRINTMAGMA_SOURCES:=					\
		auxfuncs.c cxx/matrix_header.cpp	\
		cxx/files.cpp				\
		cxx/matrix_line.cpp cxx/printmagma.cpp
-include .bw-printmagma.mkf

ifeq ($(BROKEN),also)
BW_MASTER_SOURCES:=					\
		auxfuncs.c				\
		cxx/master.cpp cxx/detect_params.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp	\
		cxx/ticks.cpp cxx/matrices_printing.cpp
BW_MASTER_LIBS:=-lgmp -lgmpxx
-include .bw-master.mkf
endif

ifeq (0,0)
BW_MASTER2_SOURCES:=						\
		auxfuncs.c					\
		old/bw_scalar.c old/endian.c			\
		old/filenames.c	old/options.c			\
		old/tagfile.c old/timer.c			\
		old/variables.c	old/fft_on_matrices.c		\
		old/field_common.c old/field_prime.c		\
		old/field_quad.c old/field_complex.c		\
		old/field_usage.c old/modulus_hacks.c		\
		old/e_polynomial.c old/structure.c		\
		old/structure_debug.c old/master2.c		\
		old/twisting_polynomials.c
BW_MASTER2_LIBS:=-lgmp
-include .bw-master2.mkf
endif

BW_MASTER_OLD_SOURCES:=						\
		auxfuncs.c barrier.c				\
		old/bw_scalar.c old/endian.c			\
		old/filenames.c old/options.c			\
		old/tagfile.c old/timer.c			\
		old/variables.c old/field_common.c		\
		old/field_prime.c old/field_quad.c		\
		old/field_complex.c old/field_usage.c		\
		old/modulus_hacks.c old/structure.c		\
		old/structure_debug.c old/master.c
BW_MASTER_OLD_LIBS:=-lgmp
-include .bw-master-old.mkf

BW_BENCHFFT_SOURCES:=						\
		auxfuncs.c					\
		old/bw_scalar.c old/endian.c			\
		old/filenames.c old/options.c			\
		old/tagfile.c old/timer.c			\
		old/variables.c old/fft_on_matrices.c		\
		old/field_common.c old/field_prime.c		\
		old/field_quad.c old/field_complex.c		\
		old/field_usage.c old/modulus_hacks.c		\
		old/structure.c old/structure_debug.c		\
		old/testfft.c
BW_BENCHFFT_LIBS:=-lgmp
-include .bw-benchfft.mkf

BW_FINDSOL_SOURCES:=						\
		auxfuncs.c					\
		old/bw_scalar.c old/endian.c			\
		old/filenames.c old/options.c			\
		old/tagfile.c old/variables.c			\
		old/field_common.c old/field_prime.c		\
		old/field_quad.c old/field_complex.c		\
		old/field_usage.c old/modulus_hacks.c		\
		old/structure.c old/structure_debug.c		\
		old/findsol.c
BW_FINDSOL_LIBS:=-lgmp
-include .bw-findsol.mkf

BW_REDRESS_SOURCES:=					\
		auxfuncs.c				\
		old/bw_scalar.c old/endian.c		\
		old/filenames.c old/options.c		\
		old/tagfile.c old/timer.c		\
		old/variables.c old/addmul.c		\
		old/bw_lvblock.c old/matrix.c		\
		old/redress.c old/threaded.c
BW_REDRESS_LIBS:=-lgmp
BW_REDRESS_EXTRA_MODELS:=mt
BW_REDRESS_MT_LIBS:=-lgmp -lrt
-include .bw-redress.mkf

BW_GATHER_SOURCES:=					\
		auxfuncs.c barrier.c			\
		cxx/detect_params.cpp			\
		cxx/gather.cpp				\
		cxx/matrix.cpp	cxx/matrix_line.cpp	\
		cxx/files.cpp cxx/matrix_header.cpp
BW_GATHER_LIBS:=-lgmp -lgmpxx
-include .bw-gather.mkf

BW_GATHER_OLD_SOURCES:=					\
		auxfuncs.c old/bw_scalar.c		\
		old/endian.c old/filenames.c		\
		old/options.c old/tagfile.c		\
		old/timer.c old/variables.c		\
		old/addmul.c old/bw_lvblock.c		\
		old/matrix.c old/gather.c		\
		old/threaded.c
BW_GATHER_OLD_LIBS:=-lgmp
BW_GATHER_OLD_EXTRA_MODELS:=mt
BW_GATHER_OLD_MT_LIBS:=-lgmp -lrt
-include .bw-gather-old.mkf

BW_TRYMUL_SOURCES:=					\
		auxfuncs.c				\
		old/bw_scalar.c old/endian.c		\
		old/filenames.c	old/options.c		\
		old/tagfile.c old/timer.c		\
		old/variables.c old/addmul.c 		\
		old/bw_lvblock.c old/matrix.c		\
		old/trymul.c old/threaded.c
BW_TRYMUL_LIBS:=-lgmp
BW_TRYMUL_EXTRA_MODELS:=mt
BW_TRYMUL_MT_LIBS:=-lgmp -lrt
-include .bw-trymul.mkf

BW_ADDVEC_SOURCES:=					\
		auxfuncs.c				\
		old/bw_scalar.c old/endian.c		\
		old/filenames.c old/tagfile.c		\
		old/variables.c old/addvec.c
BW_ADDVEC_LIBS:=-lgmp
-include .bw-addvec.mkf

BW_RANDOM_SOURCES:=			\
		auxfuncs.c		\
		cxx/matrix_header.cpp	\
		cxx/matrix_line.cpp	\
		cxx/random.cpp
-include .bw-random.mkf

BW_DUMP_SOURCES:=auxfuncs.c old/endian.c old/dump.c
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
NFSFORGE_LIBS:=-lgmp -lgmpxx
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
BW_SECURE_LIBS:=-lgmp -lgmpxx
-include .bw-secure.mkf

BW_PCOLUMN_SOURCES:=auxfuncs.c old/endian.c old/pcolumn.c
BW_PCOLUMN_LIBS:=
-include .bw-pcolumn.mkf

BW_MATRIX_WRITER_SOURCES:=auxfuncs.c old/endian.c old/matrix_writer.c
BW_MATRIX_WRITER_LIBS:=
-include .bw-matrix-writer.mkf
######################################################################

$(eval $(MAKEFILE_TRAILER))
