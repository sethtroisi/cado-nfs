TOP:=../..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_CFLAGS+=-I$(TOP)/utils -I..
MY_CFLAGS+=-I ../.. -I ../../utils
MY_CFLAGS+=-I ../../gf2x -I ../../gf2x/cantor
MY_CFLAGS+=-I $(TOPPATH)/include -I .

MY_LFLAGS+=-L$(TOP)/utils
MY_LFLAGS+=-L ../../utils -L ../../gf2x -L ../../gf2x/cantor


-include Makefile.local

ifeq ($(WITH_PTHREADS),1)
MY_CFLAGS+=-DWITH_PTHREADS
endif

ifeq ($(WITH_MPI),1)
MY_CFLAGS+=-DWITH_MPI
CC:=mpicc
CXX:=mpic++

# Work around bugs that clutter output.
CC:=$(shell $(CC) -showme | sed -e 's/-I/-isystem /g')
CXX:=$(shell $(CXX) -showme | sed -e 's/-I/-isystem /g')

CC:=$(CC) -Wl,-rpath,$(shell cat /etc/alternatives/mpilibs64)
CXX:=$(CXX) -Wl,-rpath,$(shell cat /etc/alternatives/mpilibs64)
endif


## names:=build.c matmul.cpp
## src:=$(patsubst %,%.c,$(names))
## obj:=$(patsubst %,%.o,$(names))

bins:=build prep

all:$(bins)

build_objs:=build.o matmul.o ../readbuffer.o ../readmat.o

prep_objs:=prep.o intersections.o matmul.o	\
	matmul_top.o parallelizing_info.o	\
	info_file.o ../readbuffer.o ../readmat.o	\
	auxfuncs.o config_file.o


build: $(build_objs)
	$(CXX) $(MY_LFLAGS) $(build_objs) -o $@ -lutils -lgmp

prep: $(prep_objs)
	$(CXX) $(MY_LFLAGS) $(prep_objs) -o $@ -lutils -lgmp -lgmpxx

clean:
	-rm -f *.bak *.o $(bins) *~
	-rm -rf .deps
	

# ALL source files must come here. Using $(wildcard *.c) is discouraged
# for obvious reasons.
sources:=	\
    args.cpp auxfuncs.c build.c config_file.cpp info_file.cpp		\
    intersections.cpp matmul.cpp matmul_top.cpp parallelizing_info.c	\
    prep.c preparegrid.c
dep: $(patsubst %,.deps/%.d,$(sources))
-include $(patsubst %,.deps/%.d,$(sources))

