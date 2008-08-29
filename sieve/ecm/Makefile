TOP:=../..
include $(TOP)/Makefile.common

MY_CFLAGS+=-I$(TOP)/utils -I$(TOP)
MY_LFLAGS+=-L$(TOP)/utils

PRIMEGEN=~/lib/primegen.a

all: libfacul.a
more: find_primes_hit addchain testbench test_pm1s2
tests: testbench test_pm1s2

facul_objs:=ecm.o ecm_15.o pm1.o pm1_15.o pp1.o pp1_15.o stage2.o prac_bc.o \
	makeplan.o facul.o facul_doit_ul.o facul_doit_15ul.o

libfacul.a: $(facul_objs)
	ar r $@ $^

find_primes_hit: find_primes_hit.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp
addchain: addchain.o prac_bc.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
testbench: testbench.o libfacul.a $(PRIMEGEN)
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lfacul -lgmp -lm
test_pm1s2: test_pm1s2.o pm1.o pp1.o stage2.o prac_bc.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm
pm1.o: pm1.c pm1.h $(TOP)/utils/modredc_ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDCUL $<
pm1_15.o: pm1.c pm1.h $(TOP)/utils/modredc_15ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDC15UL $<
pp1.o: pp1.c pp1.h $(TOP)/utils/modredc_ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDCUL $<
pp1_15.o: pp1.c pp1.h $(TOP)/utils/modredc_15ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDC15UL $<
ecm.o: ecm.c ecm.h $(TOP)/utils/modredc_ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDCUL $<
ecm_15.o: ecm.c ecm.h $(TOP)/utils/modredc_15ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDC15UL $<
facul_doit_ul.o: facul_doit.c facul_doit.h facul.h $(TOP)/utils/modredc_ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDCUL $<
facul_doit_15ul.o: facul_doit.c facul_doit.h facul.h $(TOP)/utils/modredc_ul.h $(TOP)/utils/modredc_15ul.h
	$(CC) $(MY_CFLAGS) -c -o $@ -DMODREDC15UL $<

stage2.o: stage2.h ../basicnt.h $(TOP)/utils/utils.h
prac_bc.o: prac_bc.h
addchain.o: $(TOP)/utils/modredc_ul.h
testbench.o: $(TOP)/utils/modredc_ul.h
test_pm1s2.o: $(TOP)/utils/modredc_ul.h
find_primes_hit.o: $(TOP)/utils/modredc_ul.h  $(TOP)/utils/modredc_ul.c
facul.o: facul.h pm1.h pp1.h ecm.h $(TOP)/utils/modredc_ul.h $(TOP)/utils/modredc_15ul.h

clean:
	-rm -f $(facul_objs) libfacul.a
	-rm -f find_primes_hit addchain testbench test_pm1s2
