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

OPTCFLAGS ?= -O2 -DNDEBUG -g -Wall -W
INCFLAGS=-I$(GMPINCLUDE) $(TIFAINCFLAGS)
CFLAGS=$(INCFLAGS) $(OPTCFLAGS) -fkeep-inline-functions -DREV=\"`svnversion .`\"
LDFLAGS=$(TIFALIB) $(GMPLIB)

all: sieve makefb las


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

fb.o : fb.c fb.h config.h basicnt.h sieve_aux.h Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) fb.c

sieve.o : sieve.c bucket.h config.h basicnt.h sieve_aux.h ../cado.h Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) sieve.c

makefb.o : makefb.c Makefile.local
	$(CC) -c -I. -I.. $(CFLAGS) makefb.c

sieve: sieve.o fb.o ../utils/libutils.a
	$(CC) $(CFLAGS) -o $@ sieve.o fb.o ../utils/libutils.a $(GMPLIB) -lm

makefb: makefb.o ../utils/libutils.a
	$(CC) $(CFLAGS) -o $@ makefb.o ../utils/libutils.a $(GMPLIB) -lm

config.h: config_h.in
	cp config_h.in config.h

las: las.c fb.o bucket.h $(TIFALIB)
	$(CC) $(CFLAGS) -o las -I.. las.c fb.o ../utils/libutils.a $(LDFLAGS) -lm

clean:
	rm -f sieve makefb las *.o *~
