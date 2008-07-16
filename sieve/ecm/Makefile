GMPBUILD=
#NDEBUG=-DNDEBUG
#EXPENSIVE=-DWANT_ASSERT_EXPENSIVE
CFLAGS=-g -O2 $(NDEBUG) $(EXPENSIVE) -Wall -W # -Wconversion -pedantic
INCLUDE=-I ../../utils/ -I ../../
CC=gcc
TARGETS=addchain ecm.o pm1.o pp1.o stage2.o prac_bc.o # find_primes_hit 
TARGETS_TEST=testbench test_pm1s2
PRIMEGEN=~/lib/primegen.a

all: $(TARGETS)
tests: $(TARGETS_TEST)

pm1.o: pm1.c pm1.h pp1.h stage2.h ../../utils/modredc_ul.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

pp1.o: pp1.c pp1.h ../../utils/modredc_ul.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

ecm.o: ecm.c ecm.h ../../utils/modredc_ul.h 
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

stage2.o: stage2.c stage2.h ../basicnt.h ../../utils/utils.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

prac_bc.o: prac_bc.c prac_bc.h
	 $(CC) $(CFLAGS) -c -o $@ $<

find_primes_hit: find_primes_hit.c ../../utils/modredc_ul.h  ../../utils/modredc_ul.c
	$(CC) $(CFLAGS) $(GMPBUILD) $(INCLUDE) $< -lgmp ../../utils/modredc_ul.c -lm -o $@

addchain: addchain.c prac_bc.o ../../utils/modredc_ul.h ../../utils/modredc_ul.c
	$(CC) $(CFLAGS) $(INCLUDE) $< prac_bc.o ../../utils/modredc_ul.c -lm -o $@

testbench: testbench.c pm1.o pp1.o ecm.o stage2.o prac_bc.o ../../utils/modredc_ul.h ../../utils/modredc_ul.c $(PRIMEGEN)
	$(CC) $(CFLAGS) $(INCLUDE) $< ../../utils/modredc_ul.c pm1.o pp1.o ecm.o stage2.o prac_bc.o ../../utils/getprime.o -o $@ -lgmp $(PRIMEGEN) -lm

test_pm1s2: test_pm1s2.c pm1.o pp1.o stage2.o prac_bc.o ../../utils/modredc_ul.h ../../utils/modredc_ul.c ../../utils/getprime.o
	$(CC) $(CFLAGS) $(INCLUDE) $< ../../utils/modredc_ul.c ../../utils/getprime.o pm1.o pp1.o stage2.o prac_bc.o -o $@ -lgmp
clean:
	rm -f $(TARGETS) 
