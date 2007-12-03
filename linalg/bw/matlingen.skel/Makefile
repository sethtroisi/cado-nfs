
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
IFLAGS:=-iquote $(TOPPATH)/include -iquote .
LFLAGS:=

# No more variable changes after this include.
include $(TOPPATH)/util/main.mki

TARNAME:=matlingen-0.5
DISTFILES:=$(shell cat 02FILES)

######################################################################


structure_automatic.h: structure.m4
	m4 -Doutputfile=MAIN < $< > $@

structure_inlines_automatic.h: structure.m4
	m4 -Doutputfile=COMPANION < $< > $@

structure_debug_automatic.h:   structure_debug.m4
	m4 -Doutputfile=DEBUG < $< > $@

BW_MASTER_SOURCES:=						\
		auxfuncs.c \
		bw_scalar.c \
		timer.c			\
		master_common.c	\
		variables.c field_common.c		\
		field_prime.c field_quad.c		\
		field_complex.c field_usage.c		\
		modulus_hacks.c structure.c		\
		structure_debug.c master.c
BW_MASTER_LIBS:=-lgmp
-include .bw-master.mkf

BW_MASTER2_SOURCES:=						\
		auxfuncs.c					\
		bw_scalar.c	\
		timer.c			\
		variables.c	fft_on_matrices.c		\
		field_common.c field_prime.c		\
		field_quad.c field_complex.c		\
		field_usage.c modulus_hacks.c		\
		e_polynomial.c structure.c		\
		structure_debug.c master2.c		\
		master_common.c	\
		twisting_polynomials.c
BW_MASTER2_LIBS:=-lgmp
-include .bw-master2.mkf

######################################################################

$(eval $(MAKEFILE_TRAILER))
