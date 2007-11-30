include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/


CFLAGS=-g -Wall -I $(GMPINCLUDE) -DREV=\"`svnversion .`\"
LDFLAGS=../../utils/libutils.a $(GMPLIB)

all: fmalgsqrt

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

algsqrt: algsqrt.c poly.o
	gcc $(CFLAGS) -I../.. -o $@ $< poly.o $(LDFLAGS) 

fmalgsqrt: fmalgsqrt.c poly.o
	gcc $(CFLAGS) -I../.. -o $@ $< poly.o $(LDFLAGS) 

poly.o: poly.c
	gcc $(CLFAGS) -c $< -o $@

ratsqrt: ratsqrt.c
	gcc $(CFLAGS) -I../.. -o $@ $< $(LDFLAGS)


clean:
	rm -f poly.o algsqrt ratsqrt

