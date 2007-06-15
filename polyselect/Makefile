include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/


CFLAGS=-O2 -g -I$(GMPINCLUDE) -DREV=\"`svnversion`\"
LDFLAGS=$(GMPLIB) -lm

polyselect: polyselect.c ../cado.h
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

clean:
	rm polyselect *~
