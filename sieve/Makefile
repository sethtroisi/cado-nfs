# Usage: make <target> 

# Paths and anything specific to a host go in Makefile.local. It can specify
# the full path/filename to the ECM libarary in ECMLIB, and full path to ECM 
# header file in ECMINCLUDE. If not specified, defaults are derived by the
# variable assignments below.
# Same thing for GMPLIB and GMPINCLUDE for the GMP library.

include Makefile.local

CFLAGS=-O2 -g -Wall -W -DREV=\"`svnversion`\"
LDFLAGS=

# Set default filenames/paths if not already set.
ECM ?= /usr/local/
GMP ?= /usr/local/
ECMLIB ?= $(ECM)/lib/libecm.a
ECMINCLUDE ?= $(ECM)/include/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

all: sieve checknorms

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

sieve: sieve.c fb.c ../polyfile.c config.h sieve_aux.h
	gcc -I. -I.. -I$(GMPINCLUDE) $(CFLAGS) -o $@ -Wall -Wextra sieve.c fb.c ../polyfile.c $(GMPLIB) -lm

checknorms: checknorms.c
	gcc -I. -I.. -I$(ECMINCLUDE) $(CFLAGS) -o $@ $< $(ECMLIB) $(GMPLIB) -lm

config.h:
	echo '#define HAVE_MSRH 1'   >  $@
	echo '#define WANT_ASSERT 1' >> $@
	echo '#define L1SIZE 32768' >>  $@

clean:
	rm -f sieve checknorms
