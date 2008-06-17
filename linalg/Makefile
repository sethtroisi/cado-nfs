TOP:=..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

binaries:=allsqrt replay merge purge freerel linalg characters duplicates \
	  transpose checkdep balance
all:	$(binaries)

procrels: procrels.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
hashpair: hashpair.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

hashpair.o: hashpair.h
prune.o: prune.h
merge.o: merge.h

freerel: freerel.o ../sieve/fb.o hashpair.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
duplicates: duplicates.o hashpair.o files.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
purge: purge.o hashpair.o files.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
merge: merge.o sparse.o prune.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lm
replay: replay.o sparse.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
checkdep: checkdep.o
	$(CC) $(MY_LFLAGS) -o $@ $^
linalg: linalg.o gauss.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
balance: balance.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lm

gauss.o: gauss.c
	$(CC) $(MY_CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $< -c -o $@
gauss.o: gauss.h

transpose: transpose.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^

characters: characters.o files.o gauss.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

allsqrt: allsqrt.o files.o hashpair.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp

matsort: matsort.o gauss.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^

matprune: matprune.c
	$(CC) $(MY_LFLAGS) -o $@ $^

######################################################################
## MPI section
######################################################################
MPI_LIB=$(MPI_DIR)/lib
MPI_INC=-I$(MPI_DIR)/include
MPI_CC=mpicc

mpimerge.o: merge.c
	$(MPI_CC) -DUSE_MPI $(MPI_INC) $(MY_CFLAGS) $< -c -o $@

mpimerge: mpimerge.o sparse.o prune.o gzip.o
	$(MPI_CC) -DUSE_MPI $(MPI_INC) $(MY_LFLAGS) -L$(MPI_LIB) -o $@ $^ -lutils -lm

hismerge: hismerge.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

######################################################################

$(binaries): $(TOP)/utils/libutils.a

clean:
	-rm -f $(binaries) *.o *~ *.bak
