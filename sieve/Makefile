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

all: sieve makefb


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

sieve: sieve.c fb.c ../utils/libutils.a config.h sieve_aux.h ../cado.h Makefile.local
	gcc -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) -o $@ -Wall -W sieve.c fb.c ../utils/libutils.a $(GMPLIB) -lm

makefb: makefb.c Makefile.local
	gcc -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) -o $@ $< ../utils/libutils.a $(GMPLIB) -lm

config.h: config_h.in
	cp config_h.in config.h

clean:
	rm -f sieve makefb
