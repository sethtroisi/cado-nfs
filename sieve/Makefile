# Usage: make <target> 

# Paths and anything specific to a host go in Makefile.local
# It must define ECMLIB which is the full path to and filename of the
# libecm.a library, and GMPLIB which is full path/name or just "-lgmp"

include Makefile.local

CFLAGS=-O2 -g -Wall -W
LDFLAGS=

all: sieve checknorms

sieve: sieve.c fb.c ../polyfile.c config.h
	gcc -I. -I.. $(CFLAGS) -o $@ -Wall -Wextra sieve.c fb.c ../polyfile.c $(GMPLIB) -lm

checknorms: checknorms.c
	gcc -I.. -I$(ECM)/include $(CFLAGS) -o $@ $< $(ECMLIB) $(GMPLIB) -lm

config.h:
	echo '#define HAVE_MSRH 1'   >  $@
	echo '#define WANT_ASSERT 1' >> $@
	echo '#define L1SIZE 32768' >>  $@

clean:
	rm -f sieve checknorms
