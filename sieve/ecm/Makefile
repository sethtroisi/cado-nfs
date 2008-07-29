TOP:=../..
include $(TOP)/Makefile.common

MY_CFLAGS+=-I$(TOP)/utils -I$(TOP)
MY_LFLAGS+=-L$(TOP)/utils

PRIMEGEN=~/lib/primegen.a

all: libfacul.a
more: find_primes_hit addchain testbench test_pm1s2
tests: testbench test_pm1s2

facul_objs:=ecm.o pm1.o pp1.o stage2.o prac_bc.o facul.o

libfacul.a: ecm.o pm1.o pp1.o stage2.o prac_bc.o facul.o
	ar r $@ ecm.o pm1.o pp1.o stage2.o prac_bc.o facul.o

find_primes_hit: find_primes_hit.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp
addchain: addchain.o prac_bc.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
testbench: testbench.o pm1.o pp1.o ecm.o stage2.o prac_bc.o $(PRIMEGEN)
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
test_pm1s2: test_pm1s2.o pm1.o pp1.o stage2.o prac_bc.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm

pm1.o: pm1.h pp1.h stage2.h $(TOP)/utils/modredc_ul.h
pp1.o: pp1.h $(TOP)/utils/modredc_ul.h
ecm.o: ecm.h $(TOP)/utils/modredc_ul.h 
stage2.o: stage2.h ../basicnt.h $(TOP)/utils/utils.h
prac_bc.o: prac_bc.h
facul.o: facul.h pm1.h pp1.h ecm.h
addchain.o: $(TOP)/utils/modredc_ul.h
testbench.o: $(TOP)/utils/modredc_ul.h
test_pm1s2.o: $(TOP)/utils/modredc_ul.h
find_primes_hit.o: $(TOP)/utils/modredc_ul.h  $(TOP)/utils/modredc_ul.c

clean:
	-rm -f $(TARGETS) *.o libfacul.a
	-rm -f find_primes_hit addchain testbench test_pm1s2
