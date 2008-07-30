
TOP:=..
include $(TOP)/Makefile.common

# Some compilation defaults specific to this directory
CFLAGS+=-DSSE_NORM_INIT
CFLAGS+=-fkeep-inline-functions
# posix_memalign is not available everywhere... For the time being there is no
# automatic detection. HAVE_POSIX_MEMALIGN is set to one by default. If this
# fails, try to set it to zero.
HAVE_POSIX_MEMALIGN ?= 1
ifneq ($(HAVE_POSIX_MEMALIGN), 0)
    CFLAGS+=-DHAVE_POSIX_MEMALIGN=1
else
    CFLAGS+=-DHAVE_POSIX_MEMALIGN=0
endif

include $(TOP)/Makefile.tifa-support

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils -Lecm
MY_CFLAGS+=-I$(TOP)/utils -Iecm

all: sieve makefb las

sieve: sieve.o fb.o ; 
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -ltifa -lgmp -lm
makefb: makefb.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -ltifa -lgmp -lm
las: las.o fb.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lfacul -lutils -ltifa -lgmp -lm

# Some dependencies on headers.
las.o: bucket.h ecm/facul.h
fb.o : fb.c fb.h config.h basicnt.h
sieve.o : sieve.c bucket.h config.h basicnt.h $(TOP)/cado.h
makefb.o : makefb.c


# Also give some info on extra dependencies.
sieve makefb las: $(TOP)/utils/libutils.a
las: $(TIFA_LIBDIR)/libtifa.a


config.h: config_h.in
	cp config_h.in config.h

clean:
	rm -f sieve makefb las *.o
