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

.c.o:
	gcc $(CFLAGS) -I../.. -o $@ -c $<

algsqrt: algsqrt.o poly.o
	gcc -o $@ $^ $(LDFLAGS) 

fmalgsqrt: fmalgsqrt.o poly.o
	gcc -o $@ $^ $(LDFLAGS)

ratsqrt: ratsqrt.o
	gcc -o $@ $^ $(LDFLAGS)

clean:
	-rm -f poly algsqrt ratsqrt *.o *~

