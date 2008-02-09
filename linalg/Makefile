include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/
OPTCFLAGS ?= -g -W -Wall -O2

SVNVER:=$(shell svnversion .)
CFLAGS=$(OPTCFLAGS) -DREV=\"$(SVNVER)\"

LDFLAGS=../utils/libutils.a $(GMPLIB) -lm# -lefence
INC=-I $(GMPINCLUDE) -I.. -I../utils

binaries:=allsqrt replay merge purge freerel linalg characters duplicates \
	  transpose checkdep

all:	$(binaries)

Makefile.local:
	@echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local

procrels: procrels.c ../utils/libutils.a
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

hashpair: hashpair.c hashpair.h
	gcc $(CFLAGS) $(INC) $< $(LDFLAGS) -o $@

prune.o: prune.c prune.h
	gcc -c $(CFLAGS) $(INC) $< -o $@

freerel: freerel.c ../utils/libutils.a ../sieve/fb.o
	gcc $(CFLAGS) $(INC) $< ../sieve/fb.o $(LDFLAGS) -o $@

duplicates: duplicates.c ../utils/libutils.a hashpair.o files.o
	gcc $(CFLAGS) $(INC) $< hashpair.o files.o $(LDFLAGS) -o $@

purge: purge.c ../utils/libutils.a hashpair.o files.o
	gcc $(CFLAGS) $(INC) $< hashpair.o files.o $(LDFLAGS) -o $@

merge: merge.c merge.h ../utils/libutils.a sparse.o prune.o
	gcc $(CFLAGS) $(INC) $< sparse.o prune.o $(LDFLAGS) -o $@

replay: replay.c ../utils/libutils.a sparse.o
	gcc $(CFLAGS) $(INC) $< sparse.o $(LDFLAGS) -o $@

linalg: linalg.o gauss.o readmat.o
	gcc $(CFLAGS) $^ $(LDFLAGS) -o $@

gauss.o: gauss.c gauss.h
	gcc -c $(CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $(INC) $< -o $@

.c.o: merge.h sparse.h readmat.h ../utils/utils.h
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
