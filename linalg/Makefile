TOP:=..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

binaries:=allsqrt replay merge purge freerel linalg characters duplicates \
	  transpose checkdep balance apply_perm mkbitstrings
all:	$(binaries)

active-targets:	\
	freerel duplicates purge merge replay transpose \
        balance characters allsqrt apply_perm mkbitstrings

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
duplicates: duplicates.o hashpair.o files.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp
purge: purge.o hashpair.o files.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm
merge: merge.o markowitz.o swar.o mst.o report.o dclist.o merge_mono.o sparse_mat.o sparse.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lm
postprocess: postprocess.o markowitz.o mst.o report.o merge_mono.o sparse_mat.o sparse.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lm
replay: replay.o sparse.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
checkdep: checkdep.o
	$(CC) $(MY_LFLAGS) -o $@ $^
linalg: linalg.o gauss.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils
balance: balance.o readbuffer.o rowset_heap.o
	$(CXX) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp -lm

mkbitstrings: mkbitstrings.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

gauss.o: gauss.c
	$(CC) $(MY_CFLAGS) -DMULTI_ROW=3 -DNO_MAIN $< -c -o $@
gauss.o: gauss.h

transpose: transpose.o readmat.o
	$(CC) $(MY_LFLAGS) -o $@ $^

characters: characters.o files.o gauss.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils -lgmp

allsqrt: allsqrt.o files.o hashpair.o
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

# Work around bugs that clutter output. Avoid forced optimization level.
mpicclibs:=$(shell $(MPI_CC) -show | sed -e 's/^.*-L\(\S*\)\>.*$$/\1/')
MPI_CC:=$(shell $(MPI_CC) -show | sed -e 's/-I/-isystem /g' -e 's/-O[0-9]//g')
MPI_CC:=$(MPI_CC) -DUSE_MPI
MPI_CC+=$(MPI_INC)
MPI_CC2:=$(MPI_CC) -Wl,-rpath,$(mpicclibs)
MPI_CC2+=-L$(MPI_LIB)

merge_mono_usempi.o: merge_mono.c
	$(MPI_CC) $(MY_CFLAGS) $< -c -o $@

merge_usempi.o: merge.c
	$(MPI_CC) $(MY_CFLAGS) $< -c -o $@

merge_mpi.o: merge_mpi.c
	$(MPI_CC) $(MY_CFLAGS) $< -c -o $@

mpimerge: merge_usempi.o merge_mpi.o merge_mono_usempi.o sparse.o
	$(MPI_CC2) $(MY_LFLAGS) -o $@ $^ -lutils -lm

hismerge: hismerge.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lutils

######################################################################

$(binaries): $(TOP)/utils/libutils.a

clean:
	-rm -f $(binaries) *.o *~ *.bak
