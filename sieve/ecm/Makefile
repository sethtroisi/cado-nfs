TOP:=../..
include $(TOP)/Makefile.common

MY_CFLAGS+=-I$(TOP)/utils -I$(TOP)
MY_LFLAGS+=-L$(TOP)/utils

all: libfacul.a
more: find_primes_hit addchain testbench
tests: testbench

facul_objs:=ecm_ul.o ecm_15ul.o pm1_ul.o pm1_15ul.o pp1_ul.o pp1_15ul.o \
	facul_doit_ul.o facul_doit_15ul.o stage2.o prac_bc.o \
	makeplan.o facul.o 

libfacul.a: $(facul_objs)
	ar r $@ $^

find_primes_hit: find_primes_hit.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp
addchain: addchain.o prac_bc.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
testbench: testbench.o libfacul.a
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
ecm_check: ecm_check.o libfacul.a
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
ecm_check.o: facul.h ecm.h $(TOP)/utils/modredc_ul_default.h

stage2.o: stage2.h ../basicnt.h $(TOP)/utils/utils.h
prac_bc.o: prac_bc.h
addchain.o: $(TOP)/utils/modredc_ul.h
testbench.o: facul.h pm1.h pp1.h ecm.h $(TOP)/utils/modredc_ul.h
find_primes_hit.o: $(TOP)/utils/modredc_ul.h  $(TOP)/utils/modredc_ul.c
facul.o: facul.h pm1.h pp1.h ecm.h $(TOP)/utils/modredc_ul.h $(TOP)/utils/modredc_15ul.h
pm1_ul.o pp1_ul.o ecm_ul.o facul_doit_ul.o : $(TOP)/utils/modredc_ul.h $(TOP)/utils/modredc_ul_default.h
pm1_15ul.o pp1_15ul.o ecm_15ul.o facul_doit_15ul.o : $(TOP)/utils/modredc_15ul.h $(TOP)/utils/modredc_15ul_default.h 

clean:
	-rm -f $(facul_objs) libfacul.a addchain.o testbench.o
	-rm -f find_primes_hit addchain testbench
