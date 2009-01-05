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

# offer the possibility of changing the mpi library by specifying other
# compilers on the make command line.
#
# As an alternative to specifying both MPICC and MPICXX, it's possible to
# give the MPI= parameter, pointing to the bin/ directory of the MPI
# installed tree.
ifneq ($(MPI),)
MPICC:=$(MPI)/mpicc
MPICXX:=$(MPI)/mpicxx
else
MPICC:=mpicc
MPICXX:=mpicxx
endif

ifeq ($(WITH_MPI),1)
MY_CFLAGS+=-DWITH_MPI
CC:=$(MPICC)
CXX:=$(MPICXX)

# Work around bugs that clutter output. Avoid forced optimization level.
CC:=$(shell $(CC) -show | sed -e 's/-I/-isystem /g' -e 's/-O[0-9]//g')
CXX:=$(shell $(CXX) -show | sed -e 's/-I/-isystem /g' -e 's/-O[0-9]//g')

ifeq ($(CC),)
$(error "You must install an mpicc compiler")
endif

mpicclibs:=$(shell $(MPICC) -show | sed -e 's/^.*-L\(\S*\)\>.*$$/\1/')
mpicxxlibs:=$(shell $(MPICXX) -show | sed -e 's/^.*-L\(\S*\)\>.*$$/\1/')

# CC:=$(CC) -Wl,-rpath,$(shell cat /etc/alternatives/mpilibs64)
# CXX:=$(CXX) -Wl,-rpath,$(shell cat /etc/alternatives/mpilibs64)
CC:=$(CC) -Wl,-rpath,$(mpicclibs)
CXX:=$(CXX) -Wl,-rpath,$(mpicxxlibs)
endif


## names:=build.c matmul-sliced.cpp
## src:=$(patsubst %,%.c,$(names))
## obj:=$(patsubst %,%.o,$(names))

bins:=bench build prep test_hello test_pthread_reduce

all:$(bins)

build_objs:=build.o ../readbuffer.o ../readmat.o libmatmul.a
bench_objs:=bench.o ../readbuffer.o ../readmat.o libmatmul.a

libmatmul_objs:=matmul-sliced.o matmul-basic.o matmul.o

libmatmul.a: $(libmatmul_objs)
	-rm -f $@
	ar r $@ $(libmatmul_objs)

test_hello_objs:=test_hello.o parallelizing_info.o

prep_objs:=prep.o intersections.o 			\
	matmul_top.o parallelizing_info.o		\
	info_file.o ../readbuffer.o ../readmat.o	\
	auxfuncs.o xvectors.o				\
	libmatmul.a


# We need to make sure that pthreads are always included.
common_libs:=-lutils -lgmp
ifeq ($(WITH_PTHREADS),1)
# As a matter of cact, it is fairly common for mpi libs to include
# pthreads automatically. Hence in some cases, we would need this ifndef
# to avoid very redundant link lines. But strictly specaking, this is a
# behaviour we can't rely on.
# ifneq ($(WITH_MPI),1)
common_libs += -lpthread
# endif
endif

test_pthread_reduce_objs:=test_pthread_reduce.o

test_hello: $(test_hello_objs)
	$(CXX) $(MY_LFLAGS) $(test_hello_objs) -o $@ $(common_libs)

build: $(build_objs)
	$(CXX) $(MY_LFLAGS) $(build_objs) -o $@ $(common_libs)

bench: $(bench_objs)
	$(CXX) $(MY_LFLAGS) $(bench_objs) -o $@ $(common_libs)

prep: $(prep_objs)
	$(CXX) $(MY_LFLAGS) $(prep_objs) -o $@ $(common_libs) -lgmpxx

test_pthread_reduce: $(test_pthread_reduce_objs)
	$(CC) $(MY_LFLAGS) $(test_pthread_reduce_objs) -o $@ -lpthread -O4


ifeq ($(WITH_MPI),1)
ifeq ($(WITH_PTHREADS),1)
cregexp:=.
else
cregexp:='thr=1x1'
endif
else
ifeq ($(WITH_PTHREADS),1)
cregexp:='mpi=1x1'
else
cregexp:='mpi=1x1.*thr=1x1'
endif
endif

check: test_hello
	grep $(cregexp) checks.list.light | ./checks.sh

test: test_hello
	grep $(cregexp) checks.list | ./checks.sh

clean:
	-rm -f sketch.aux sketch.dvi sketch.log sketch.ps
	-rm -f *.bak *.o $(bins) *~
	-rm -rf .deps
	

# ALL source files must come here. Using $(wildcard *.c) is discouraged
# for obvious reasons.
sources:=	\
    auxfuncs.c bench.c build.c info_file.c		\
    intersections.cpp matmul-sliced.cpp matmul_top.c	\
    parallelizing_info.c prep.c test_hello.c xvectors.c	\
    matmul.c matmul-basic.c

dep: $(patsubst %,.deps/%.d,$(sources))
-include $(patsubst %,.deps/%.d,$(sources))

