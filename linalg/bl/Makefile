# Default flags ; can be overriden by methods specified in main.mkf,
# like:
#
# STATIC=yes --> adds -static
# GMP=<path> --> uses GMP as in path
# MPFR=<path> --> uses MPFR as in path

TOPPATH:=.

CC:=mpicc
CXX:=mpicc    #g++
ARQCV:=ar qcv
OFLAGS:=
PFLAGS:=
IFLAGS:=
LFLAGS:=

include $(TOPPATH)/util/main.mki

######################################################################


#TEST_SOURCES:=Test.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
#TEST_LIBS:=
#-include .Test.mkf

#TEST2_SOURCES:=test2.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
#TEST2_LIBS:=
#-include .test2.mkf

#TEST4_SOURCES:=test4.c 
#TEST4_LIBS:=
#-include .test4.mkf

#LANCZOSTEST_SOURCES:=LanczosTest.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
#LANCZOSTEST_LIBS:=
#-include .LanczosTest.mkf


#LANCZOS2CORE_SOURCES:=Lanczos2Core.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c
#LANCZOS2CORE_LIBS:=
#-include .Lanczos2Core.mkf


#LANCZOS_SOURCES:=Lanczos.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
#LANCZOS_LIBS:=
#-include .Lanczos.mkf


GENERATEDATAFILES_SOURCES:=GenerateDataFiles.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
GENERATEDATAFILES_LIBS:=
-include .GenerateDataFiles.mkf


GETDATAFILES_SOURCES:=GetDataFiles.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c
GETDATAFILES_LIBS:=
-include .GetDataFiles.mkf

#LANCZOSMPI_SOURCES:=LanczosMPI.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c
#LANCZOSMPI_LIBS:=
#-include .LanczosMPI.mkf


#LANCZOSMPIOPT_SOURCES:=LanczosMPIOpt.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c
#LANCZOSMPIOPT_LIBS:=
#-include .LanczosMPIOpt.mkf

#LANCZOSMPIOPTCOMPLETE_SOURCES:=LanczosMPIOptComplete.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c
#LANCZOSMPIOPTCOMPLETE_LIBS:=
#-include .LanczosMPIOptComplete.mkf

MAINLANCZOS_SOURCES:=MainLanczos.c ProcLanczos.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c struct.c
LANCZOSMPIOPTCOMPLETE_LIBS:=
-include .MainLanczos.mkf

TESTMULT_SOURCES:=TestMult.c ProcLanczos.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c struct.c
LANCZOSMPIOPTCOMPLETE_LIBS:=
-include .TestMult.mkf

TESTMULTBIT_SOURCES:=TestMultBit.c ProcLanczos.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c struct.c
LANCZOSMPIOPTCOMPLETE_LIBS:=
-include .TestMultBit.mkf


#TESTMULTBIT_SOURCES:=toto.c 
#LANCZOSMPIOPTCOMPLETE_LIBS:=
#-include .toto.mkf


#MAINLANCZOSTESTS_SOURCES:=MainLanczosTestes.c ProcLanczosTestes.c basic.c echelon.c mat_ops.c random.c alloc.c ReadWrite.c timing.c struct.c
#MAINLANCZOSTESTS_LIBS:=
#-include .MainLanczosTestes.mkf

######################################################################

ifneq ($(ALLDEPS),)
-include $(ALLDEPS)
endif

late_all: $(BIN_TARGETS) $(LIB_TARGETS)

dep: $(ALLDEPS)
