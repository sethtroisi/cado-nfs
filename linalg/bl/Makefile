
# For a static build, say LDFLAGS=-static
CC:=mpicc
CXX:=mpic++

TOP:=../..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

ifneq ($(FAKE),)
CFLAGS+=-DFAKEMPI
endif

binaries:=GenerateDataFiles GetDataFiles MainLanczos TestMult
all: $(binaries)

#Test: Test.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#Test2: test2.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#Test4: test4.o 
#	$(CC) $(MY_LFLAGS) -o $@ $^
#LanczosTest: LanczosTest.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#Lanczos2Core: Lanczos2Core.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#Lanczos: Lanczos.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#LanczosMPI: LanczosMPI.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#LanczosMPIOpt: LanczosMPIOpt.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
#LanczosMPIOptComplete: LanczosMPIOptComplete.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
#	$(CC) $(MY_LFLAGS) -o $@ $^
GenerateDataFiles: GenerateDataFiles.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
	$(CC) $(MY_LFLAGS) -o $@ $^
GetDataFiles: GetDataFiles.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o
	$(CC) $(MY_LFLAGS) -o $@ $^
MainLanczos: MainLanczos.o ProcLanczos.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
	$(CC) $(MY_LFLAGS) -o $@ $^
TestMult: TestMult.o ProcLanczos.o basic.o echelon.o mat_ops.o random.o alloc.o ReadWrite.o timing.o
	$(CC) $(MY_LFLAGS) -o $@ $^

alloc.o: alloc.c alloc.h 
basic.o: basic.c alloc.h random.h mat_ops.h struct.h \
  echelon.h basic.h
echelon.o: echelon.c alloc.h mat_ops.h struct.h \
  echelon.h
GenerateDataFiles.o: GenerateDataFiles.c \
  basic.h random.h alloc.h mat_ops.h struct.h echelon.h ReadWrite.h
GetDataFiles.o: GetDataFiles.c basic.h random.h \
  alloc.h mat_ops.h struct.h echelon.h ReadWrite.h
MainLanczos.o: MainLanczos.c mpi_select.h \
  basic.h random.h alloc.h mat_ops.h \
  struct.h echelon.h ReadWrite.h timing.h ProcLanczos.h
mat_ops.o: mat_ops.c mpi_select.h \
   mat_ops.h struct.h alloc.h \
  ReadWrite.h timing.h
ProcLanczos.o: ProcLanczos.c mpi_select.h \
   basic.h random.h alloc.h mat_ops.h \
  struct.h echelon.h ReadWrite.h timing.h
random.o: random.c random.h alloc.h 
ReadWrite.o: ReadWrite.c mpi_select.h \
   struct.h ReadWrite.h mat_ops.h
struct.o: struct.c
TestMultBit.o: TestMultBit.c mpi_select.h fakempi.h \
  basic.h random.h alloc.h mat_ops.h struct.h echelon.h ReadWrite.h \
  timing.h ProcLanczos.h
TestMult.o: TestMult.c mpi_select.h \
   basic.h random.h alloc.h mat_ops.h \
  struct.h echelon.h ReadWrite.h timing.h ProcLanczos.h
timing.o: timing.c timing.h

clean:
	-rm -f $(binaries) *.o *~ *.bak
