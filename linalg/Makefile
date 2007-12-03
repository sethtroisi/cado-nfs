include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

CFLAGS=-g -Wall -DREV=\"`svnversion .`\" -O2 
LDFLAGS=../utils/libutils.a $(GMPLIB)
INC=-I $(GMPINCLUDE) -I..

all: allsqrt replay merge purge procrels linalg characters matsort matprune

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
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

linalg: linalg.c gauss.o
	gcc $(CFLAGS) gauss.o $(INC) $< $(LDFLAGS) -o $@

gauss.o: gauss.c
	gcc -c $(CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $(INC) $< -o $@

files.o: files.c
	gcc -c $(CFLAGS) $(INC) $< -o $@

characters: characters.c files.o gauss.o
	gcc $(CFLAGS) $(INC) $< files.o gauss.o $(LDFLAGS) -o $@

allsqrt: allsqrt.c files.o hashpair.o
	gcc $(CFLAGS) $(INC) $< files.o hashpair.o $(LDFLAGS) -o $@

matsort: matsort.c gauss.o
	gcc $(CFLAGS) $(INC) $< gauss.o $(LDFLAGS) -o $@

matprune: matprune.c
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

clean:
	\rm replay merge purge procrels gauss.o linalg characters matsort matprune
