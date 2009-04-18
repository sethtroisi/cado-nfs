TOP:=../..

# Our default setting has pthreads enabled, but not mpi
WITH_PTHREADS=1

include $(TOP)/Makefile.common

# We require the cado utils library
MY_CFLAGS+=-I..
MY_CFLAGS+=-I $(TOP) -I $(TOP)/utils
MY_CFLAGS+=-I $(TOP)/gf2x -I $(TOP)/cantor
MY_CFLAGS+=-I.

# Use "-L<path>", not "-L <path>" !
MY_LFLAGS+=-L$(TOP)/utils -L$(TOP)/gf2x -L$(TOP)/cantor

# -include Makefile.local

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

HAVE_POSIX_MEMALIGN ?= 1
ifneq ($(HAVE_POSIX_MEMALIGN), 0)
    MY_CFLAGS+=-DHAVE_POSIX_MEMALIGN=1
else
    MY_CFLAGS+=-DHAVE_POSIX_MEMALIGN=0
endif

ifdef DISABLE_ASM
    MY_CFLAGS+=-DDISABLE_ASM
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

# That's really a quirk, but otherwise gcc output is awfully cluttered
# with gcc 4.1.x
ifneq ($(TRY_TO_REMOVE_GCC_NOLINK_WARNINGS),)
# MY_LFLAGS has late expansion
lflags_extracted:=$(filter -Wl% -l%,$(CC))
MY_LFLAGS+=$(lflags_extracted)
CC:=$(filter-out -Wl% -l%,$(CC))
CXX:=$(filter-out -Wl% -l%,$(CXX))
endif


## names:=build.c matmul-sliced.cpp
## src:=$(patsubst %,%.c,$(names))
## obj:=$(patsubst %,%.o,$(names))

bins:=bench build prep krylov mksol

default-target: all
allbins:=test_hello split acollect random shell

common_libs:=-lutils -lgmp
common_link_deps:=$(TOP)/utils/libutils.a

# We need to make sure that pthreads are always included.
ifeq ($(WITH_PTHREADS),1)
# As a matter of fact, it is fairly common for mpi libs to include
# pthreads automatically. Hence in some cases, we would need this ifndef
# to avoid very redundant link lines. But strictly speaking, this is a
# behaviour we can't rely on.
# ifneq ($(WITH_MPI),1)
common_libs += -lpthread
# endif
endif
common_libs += $(EXTRA_LIBS)


test_hello_objs:=test_hello.o parallelizing_info.o
test_hello: $(test_hello_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) $(test_hello_objs) -o $@ $(common_libs)
-include $(patsubst %,.deps/%.d,$(test_hello_objs))

random_objs:=random.o
random: $(random_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) $(random_objs) -o $@ $(common_libs) -lm
-include $(patsubst %,.deps/%.d,$(random_objs))

split_objs:=split.o bw-common.o
split: $(split_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) $(split_objs) -o $@ $(common_libs)
-include $(patsubst %,.deps/%.d,$(split_objs))

acollect_objs:=acollect.o bw-common.o
acollect: $(acollect_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) $(acollect_objs) -o $@ $(common_libs)
-include $(patsubst %,.deps/%.d,$(acollect_objs))

shell_objs:=shell.o parallelizing_info.o bw-common.o bw-common-mpi.o
shell: $(shell_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) $(shell_objs) -o $@ $(common_libs)
-include $(patsubst %,.deps/%.d,$(shell_objs))

all_objlists:=libmatmul build bench prep gather secure krylov mksol shell

libmatmul_objs:=					\
    matmul-basic.o matmul-sliced.o matmul-threaded.o	\
    matmul.o

ifndef	DISABLE_ASM
ifeq ($(shell uname -m),x86_64)
libmatmul_objs+=matmul-sliced-asm.o
endif
endif

build_objs:=build.o libmatmul.a
bench_objs:=bench.o libmatmul.a
prep_objs:=prep.o intersections.o 			\
	matmul_top.o parallelizing_info.o info_file.o	\
	xvectors.o xymats.o debug.o bw-common.o		\
	rusage.o async.o bw-common-mpi.o		\
	libmatmul.a
krylov_objs:=krylov.o intersections.o 			\
	matmul_top.o parallelizing_info.o info_file.o	\
	xvectors.o xymats.o debug.o bw-common.o		\
	rusage.o async.o bw-common-mpi.o		\
	xdotprod.o					\
	libmatmul.a
mksol_objs:=mksol.o intersections.o 			\
	matmul_top.o parallelizing_info.o info_file.o	\
	xvectors.o xymats.o debug.o bw-common.o		\
	rusage.o async.o bw-common-mpi.o		\
	xdotprod.o					\
	libmatmul.a
gather_objs:=gather.o intersections.o 			\
	matmul_top.o parallelizing_info.o info_file.o	\
	xvectors.o debug.o bw-common.o			\
	rusage.o async.o bw-common-mpi.o		\
	libmatmul.a
secure_objs:=secure.o intersections.o 			\
	matmul_top.o parallelizing_info.o info_file.o	\
	xvectors.o debug.o bw-common.o			\
	rusage.o async.o bw-common-mpi.o		\
	libmatmul.a

mat_io_objs:=../readbuffer.o ../readmat.o readmat-easy.o
build_normalobjs:=$(mat_io_objs)
bench_normalobjs:=$(mat_io_objs) debug.o
prep_normalobjs:=$(mat_io_objs) ../gauss.o
krylov_normalobjs:=$(mat_io_objs)
mksol_normalobjs:=$(mat_io_objs)
secure_normalobjs:=$(mat_io_objs)
gather_normalobjs:=$(mat_io_objs)


# This list contains the object files which do not depend on the abase
# selection.
constant_objects:=parallelizing_info.o intersections.o xvectors.o bw-common.o bw-common-mpi.o rusage.o async.o debug.o

# hardlinking like this over nfs is prone to exhibit bugs with time
# difference w.r.t the nfs server.
define hardlink_object
$(1)/$(2): $(2)
	ln -f $$< $$@
endef

define	multiple_binary
$(1)/$(2): $$($(1)_$(2)_objs) $$($(2)_normalobjs) $(common_link_deps)
	$$(CXX) $$(MY_LFLAGS) $$($(1)_$(2)_objs) $$($(2)_normalobjs) -o $$@ $$(common_libs)
-include $$(patsubst %,.deps/%.d,$$($(1)_$(2)_objs))
allbins+=$(1)/$(2)
endef

define configure_abase
$$(foreach x,$$(constant_objects),$$(eval $$(call hardlink_object,$(1),$$(x))))

$(1)/%.o: %.c
	@mkdir -p $(1)
	$$(CC) -DSELECT_ABASE_$(1) $$(MY_CFLAGS) $$< -c -o $$@
$(1)/%.o: %.S
	@mkdir -p $(1)
	$$(CC) -DSELECT_ABASE_$(1) $$(MY_CFLAGS) $$< -c -o $$@
$(1)/%.o: %.cpp
	@mkdir -p $(1)
	$$(CXX) -DSELECT_ABASE_$(1) $$(MY_CXXFLAGS) $$< -c -o $$@
.deps/$(1)/%.o.d: %.c
	@mkdir -p .deps/$(1)
	$$(CC) -DSELECT_ABASE_$(1) $$(MY_CFLAGS) -M -MM -MG -MP -MF $$@ -MT '$$@ $(1)/$$*.o' $$<
.deps/$(1)/%.o.d: %.S
	@mkdir -p .deps/$(1)
	$$(CC) -DSELECT_ABASE_$(1) $$(MY_CFLAGS) -M -MM -MG -MP -MF $$@ -MT '$$@ $(1)/$$*.o' $$<
.deps/$(1)/%.o.d: %.cpp
	@mkdir -p .deps/$(1)
	$$(CXX) -DSELECT_ABASE_$(1) $$(MY_CXXFLAGS) -M -MM -MG -MP -MF $$@ -MT '$$@ $(1)/$$*.o' $$<

$$(foreach x,$$(all_objlists),$$(eval $(1)_$$(x)_objs:=$(1)/abase-$(1).o $$(patsubst %,$(1)/%,$$($$(x)_objs))))


$(1)/libmatmul.a: $$($(1)_libmatmul_objs)
	-rm -f $$@
	ar r $$@ $$($(1)_libmatmul_objs)

# These targets make sense for all binary instances -- although nobody sane
# would ever consider running them seriously for u64n
$(call multiple_binary,$(1),build)
$(call multiple_binary,$(1),bench)
$(call multiple_binary,$(1),krylov)
$(call multiple_binary,$(1),mksol)
endef

# u64 is the prime-time version. The others can be dropped to accelerate
# compilation.
$(eval $(call configure_abase,u64))
$(eval $(call configure_abase,u64k))
$(eval $(call configure_abase,u64n))
$(eval $(call configure_abase,u128))

# Only one version of prep and secure exist.
$(eval $(call multiple_binary,u64n,prep))
$(eval $(call multiple_binary,u64,secure))
$(eval $(call multiple_binary,u64n,gather))

all: $(allbins)


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
	-rm -f *.bak *.[oa] $(bins) *~ u*/*.[oa] $(allbins)
	-rm -rf .deps
