
TOP:=../../..
include $(TOP)/Makefile.common

# We require the cado utils library
# We require the cado utils library
MY_CFLAGS+=-I$(TOP)/utils -I..
MY_CFLAGS+=-I $(TOP) -I $(TOP)/utils
MY_CFLAGS+=-I $(TOP)/gf2x -I $(TOP)/cantor
MY_CFLAGS+=-I.

MY_LFLAGS+=-L$(TOP)/utils
MY_LFLAGS+=-L$(TOP)/utils -L$(TOP)/gf2x/.libs -L$(TOP)/cantor

-include Makefile.local

all: lingen

lingen_objs:=fake_fft.o lingen_binary.o
lingen_libs:=../bw-common.o -lgmpxx -lgmp -lutils -lgf2x $(TOP)/cantor/cantor128.o

common_link_deps:=$(TOP)/utils/libutils.a

lingen: $(lingen_objs) $(common_link_deps)
	$(CXX) $(MY_LFLAGS) -o $@ $(lingen_objs) $(lingen_libs)

clean:
	-rm -f $(lingen_objs) lingen

lingen_binary.o: ../bw-common.h ../filenames.h
lingen_binary.o: fft_adapter.hpp lingen_mat_types.hpp
lingen_binary.o: $(TOP)/cado.h $(TOP)/utils/utils.h $(TOP)/utils/manu.h
lingen_binary.o: $(TOP)/gf2x/gf2x.h

fake_fft.o: fake_fft.h
fake_fft.o: $(TOP)/gf2x/gf2x.h $(TOP)/utils/utils.h $(TOP)/utils/manu.h

lingen_binary.o fake_fft.o: alloc_proxy.h
