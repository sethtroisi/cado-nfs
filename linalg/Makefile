include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

CFLAGS=-O2 -g -DREV=\"`svnversion`\"
LDFLAGS=../utils/libutils.a $(GMPLIB)
INC=-I $(GMPINCLUDE) -I..

all: procrels linalg characters matsort

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

procrels: procrels.c ../utils/libutils.a
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

linalg: linalg.c gauss.o
	gcc $(CFLAGS) gauss.o $(INC) $< $(LDFLAGS) -o $@

gauss.o: gauss.c
	gcc -c $(CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $(INC) $< -o $@

characters: characters.c gauss.o
	gcc $(CFLAGS) $(INC) $< gauss.o $(LDFLAGS) -o $@

matsort: matsort.c gauss.o
	gcc $(CFLAGS) $(INC) $< gauss.o $(LDFLAGS) -o $@

clean:
	\rm procrels gauss.o linalg characters matsort
