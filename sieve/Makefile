# Usage: make <target> 

# Paths and anything specific to a host go in Makefile.local. It can specify
# the full path/filename to the GMP library in GMPLIB, and full path to GMP
# header file in GMPINCLUDE. If not specified, defaults are derived by the
# variable assignments below.

include Makefile.local


# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

# Set USE_TIFA_FROM_SVN to one to use TIFA headers/library from local SVN
# repository. Of course TIFA has to be compiled.
USE_TIFA_FROM_SVN ?= 1

# Set default filenames/paths for TIFA if not already set.
ifeq ($(USE_TIFA_FROM_SVN), 1)
    TIFA = ../postsieve/tifa
    TIFALIB ?= $(TIFA)/build/lib/libtifa.a
    TIFAINCFLAGS = -I$(TIFA) -I$(TIFA)/lib -I$(TIFA)/lib/algo/include \
                   -I$(TIFA)/lib/data/include -I$(TIFA)/lib/utils/include
else
    TIFA ?= /usr/local
    TIFALIB ?= $(TIFA)/lib/libtifa.a
    TIFAINCFLAGS ?= $(TIFA)/include/tifa
endif

OPTCFLAGS ?= -O2 -DNDEBUG -g -Wall -W -DSSE_NORM_INIT
INCFLAGS=-I$(GMPINCLUDE) $(TIFAINCFLAGS)
# -std=c99 is needed for TIFA's inline functions
CFLAGS=$(INCFLAGS) $(OPTCFLAGS) -std=c99 -fkeep-inline-functions -DREV=\"`svnversion .`\"
LDFLAGS=$(TIFALIB) $(GMPLIB)

# posix_memalign is not available everywhere... For the time being there is no
# automatic detection. HAVE_POSIX_MEMALIGN is set to one by default. If this
# fails, try to set it to zero.

HAVE_POSIX_MEMALIGN ?= 1

ifneq ($(HAVE_POSIX_MEMALIGN), 0)
    EXTRA_CPPFLAGS = -DHAVE_POSIX_MEMALIGN=1
else
    EXTRA_CPPFLAGS = -DHAVE_POSIX_MEMALIGN=0
endif

all: sieve makefb las

Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

fb.o : fb.c fb.h config.h basicnt.h Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) fb.c

sieve.o : sieve.c bucket.h config.h basicnt.h ../cado.h Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) $(EXTRA_CPPFLAGS) sieve.c

makefb.o : makefb.c Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) makefb.c

sieve: sieve.o fb.o ../utils/libutils.a
	$(CC) $(CFLAGS) -o $@ sieve.o fb.o ../utils/libutils.a $(GMPLIB) -lm

makefb: makefb.o ../utils/libutils.a
	$(CC) $(CFLAGS) -o $@ makefb.o ../utils/libutils.a $(GMPLIB) -lm

config.h: config_h.in
	cp config_h.in config.h

las: las.c fb.o bucket.h sse_macros.h $(TIFALIB)
	$(CC) $(CFLAGS) $(EXTRA_CPPFLAGS) -o las -I.. las.c fb.o ../utils/libutils.a $(LDFLAGS) -lm

clean:
	rm -f sieve makefb las *.o *~
