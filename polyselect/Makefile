include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= -L $(GMP)/lib -lgmp
GMPINCLUDE ?= $(GMP)/include/


CFLAGS=-O2 -g -W -Wall -std=c99
EXTRAFLAGS=-I$(GMPINCLUDE) -DREV=\"`svnversion .`\"
LDFLAGS=../utils/libutils.a $(GMPLIB) -lm

polyselect: polyselect.c ../cado.h
	gcc $(CFLAGS) $(EXTRAFLAGS) -I.. $< $(LDFLAGS) -o $@

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

clean:
	-rm -f polyselect *~ *.o
