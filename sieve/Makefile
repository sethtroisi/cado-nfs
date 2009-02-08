
TOP:=..
include $(TOP)/Makefile.common

# Some compilation defaults specific to this directory
# The following flags are specific to gcc, and fix a regression in gcc 4.3.0
CFLAGS+=-fno-ivopts -fno-regmove -fno-tree-sra -fno-tree-ter

# if you compile with GCC 4.2 or later, you can also try -march/mtune=native
# instead of the default -march/mtune=generic.
# Cf http://gcc.gnu.org/ml/gcc-help/2008-07/msg00060.html

# using -DSSE_NORM_INIT seems to be slower (about 3% with gcc 4.3.0)
# CFLAGS+=-DSSE_NORM_INIT
# -fkeep-inline-functions does not yield any speedup, is it useful?
# CFLAGS+=-fkeep-inline-functions
# posix_memalign is not available everywhere... For the time being there is no
# automatic detection. HAVE_POSIX_MEMALIGN is set to one by default. If this
# fails, try to set it to zero.
HAVE_POSIX_MEMALIGN ?= 1
ifneq ($(HAVE_POSIX_MEMALIGN), 0)
    CFLAGS+=-DHAVE_POSIX_MEMALIGN=1
else
    CFLAGS+=-DHAVE_POSIX_MEMALIGN=0
endif

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils -Lecm
MY_CFLAGS+=-I$(TOP)/utils -Iecm

all: sieve makefb las

sieve: sieve.o fb.o ; 
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
makefb: makefb.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
las: las.o fb.o trialdiv.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lfacul -lutils -lgmp -lm

# Some dependencies on headers.
las.o: bucket.h ecm/facul.h
fb.o : fb.c fb.h basicnt.h config.h
sieve.o : sieve.c bucket.h config.h basicnt.h $(TOP)/cado.h
makefb.o : makefb.c
trialdiv.o : trialdiv.h fb.h

# Also give some info on extra dependencies.
sieve makefb las: $(TOP)/utils/libutils.a
las: ecm/libfacul.a


config.h: config_h.in
	cp config_h.in config.h

clean:
	rm -f sieve makefb las *.o
