TOP:=..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_LXXFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

binaries:=allsqrt replay merge purge freerel linalg characters duplicates \
	  transpose checkdep balance apply_perm
all:	$(binaries)

procrels: procrels.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
hashpair: hashpair.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

hashpair.o: hashpair.h
merge_mono.o: merge_mono.h
markowitz.o: markowitz.h
dclist.o: dclist.h
swar.o: swar.h
report.o: report.h
sparse_mat.o: sparse_mat.h
mst.o: mst.h

freerel: freerel.o ../sieve/fb.o hashpair.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
duplicates: duplicates.o hashpair.o files.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp
purge: purge.o hashpair.o files.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
merge: merge.o markowitz.o swar.o mst.o report.o dclist.o merge_mono.o sparse_mat.o sparse.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lm
replay: replay.o sparse.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
checkdep: checkdep.o
	$(CC) $(MY_LFLAGS) -o $@ $^
linalg: linalg.o gauss.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
balance: balance.o readbuffer.o rowset_heap.o
	$(CXX) $(MY_LXXFLAGS) -o $@ $^ -lutils -lgmp -lm

gauss.o: gauss.c
	$(CC) $(MY_CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $< -c -o $@
gauss.o: gauss.h

transpose: transpose.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^

characters: characters.o files.o gauss.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

allsqrt: allsqrt.o files.o hashpair.o gzip.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp

apply_perm: apply_perm.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

# Don't forget the dependency on libutils !
procrels hashpair freerel duplicates purge merge replay linalg balance: $(TOP)/utils/libutils.a

######################################################################
## MPI section
######################################################################
# MPI_DIR should have been defined in Makefile.local
MPI_LIB=$(MPI_DIR)/lib
MPI_INC=-I$(MPI_DIR)/include
MPI_CC=mpicc

merge_mpi.o: merge_mpi.c merge_mono.c merge.c
	$(MPI_CC) -DUSE_MPI $(MPI_INC) $(MY_CFLAGS) merge_mono.c -c -o merge_mono.o
	$(MPI_CC) -DUSE_MPI $(MPI_INC) $(MY_CFLAGS) merge.c -c -o merge.o
	$(MPI_CC) -DUSE_MPI $(MPI_INC) $(MY_CFLAGS) $< -c -o $@

# for xterm debug, thanks ET: gore, isn't it? But we need to export somehow
# the path for the libraries
MPI_CC2:=$(shell $(MPI_CC) -showme | sed -e 's/-I/-isystem /g')
MPI_CC2:=$(MPI_CC2) -Wl,-rpath,$(MPI_LIB)
mpimerge: merge.o merge_mpi.o merge_mono.o sparse.o gzip.o
	$(MPI_CC2) -DUSE_MPI $(MPI_INC) $(MY_LFLAGS) -L$(MPI_LIB) -o $@ $^ -lutils -lm

hismerge: hismerge.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

######################################################################

$(binaries): $(TOP)/utils/libutils.a

clean:
	-rm -f $(binaries) *.o *~ *.bak
