include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

SVNVER:=$(shell svnversion .)
CFLAGS=-g -Wall -DREV=\"$(SVNVER)\" -O2 
LDFLAGS=../utils/libutils.a $(GMPLIB)
INC=-I $(GMPINCLUDE) -I.. -I../utils

binaries:=allsqrt replay merge purge procrels linalg	\
	characters matsort matprune transpose

all:	$(binaries)

Makefile.local:
	@echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

procrels: procrels.c ../utils/libutils.a
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

hashpair: hashpair.c
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

purge: purge.c ../utils/libutils.a hashpair.o
	gcc $(CFLAGS) $(INC) $< hashpair.o $(LDFLAGS) -o $@

merge: merge.c ../utils/libutils.a
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

replay: replay.c ../utils/libutils.a
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

linalg: linalg.o gauss.o readmat.o
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

gauss.o: gauss.c
	gcc -c $(CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $(INC) $< -o $@

.c.o:
	gcc -c $(CFLAGS) $(INC) $< -o $@

transpose: transpose.o readmat.o
	gcc $(CFLAGS) $(INC) $^ $(LDFLAGS) -o $@

characters: characters.c files.o gauss.o
	gcc $(CFLAGS) $(INC) $< files.o gauss.o $(LDFLAGS) -o $@

allsqrt: allsqrt.c files.o hashpair.o
	gcc $(CFLAGS) $(INC) $< files.o hashpair.o $(LDFLAGS) -o $@

matsort: matsort.o gauss.o readmat.o
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

matprune: matprune.c
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@



clean:
	-rm -f $(binaries) *.o *~ *.bak
