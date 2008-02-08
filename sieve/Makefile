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

OPTCFLAGS ?= -O2 -g -Wall -W
CFLAGS=$(OPTCFLAGS) -fkeep-inline-functions -DREV=\"`svnversion .`\"
LDFLAGS=

all: sieve makefb las


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

fb.o : fb.c fb.h config.h basicnt.h sieve_aux.h Makefile.local
	gcc -c -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) fb.c

sieve.o : sieve.c config.h basicnt.h sieve_aux.h ../cado.h Makefile.local
	gcc -c -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) sieve.c

makefb.o : makefb.c Makefile.local
	gcc -c -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) makefb.c

sieve: sieve.o fb.o ../utils/libutils.a
	gcc $(CFLAGS) -o $@ sieve.o fb.o ../utils/libutils.a $(GMPLIB) -lm

makefb: makefb.o ../utils/libutils.a
	gcc $(CFLAGS) -o $@ makefb.o ../utils/libutils.a $(GMPLIB) -lm

config.h: config_h.in
	cp config_h.in config.h

las: las.c fb.o
	gcc -O4 -o las -I.. las.c fb.o ../utils/libutils.a $(GMPLIB) -lm

clean:
	rm -f sieve makefb las *.o *~
