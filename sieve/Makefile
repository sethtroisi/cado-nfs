# Usage: make <target> GMP=$GMP ECM=$ECM where $GMP and $ECM are the
# installation directories of GMP and GMP-ECM

GMP=
ECM=

CFLAGS=-O2 -g
LDFLAGS=-lm $(GMP)/lib/libgmp.a

sieve: sieve.c fb.c config.h
	gcc -I.. $(CFLAGS) -o $@ -Wall -Wextra sieve.c fb.c $(LDFLAGS)

checknorms: checknorms.c
	gcc -I.. -I$(ECM)/include $(CFLAGS) -o $@ $< $(ECM)/lib/libecm.a $(LDFLAGS)

config.h:
	echo '#define HAVE_MSRH 1'   >  $@
	echo '#define WANT_ASSERT 1' >> $@

clean:
	\rm sieve checknorms


