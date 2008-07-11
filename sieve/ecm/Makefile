GMPBUILD=
#NDEBUG=-DNDEBUG
CFLAGS=-g -O2 $(NDEBUG) -Wall -W # -Wconversion -pedantic
INCLUDE=-I ../../utils/ -I ../../
CC=gcc
TARGETS=addchain ecm.o pm1.o pp1.o stage2.o prac_bc.o # find_primes_hit 
TARGETS_TEST=testbench test_pm1s2
PRIMEGEN=~/lib/primegen.a

all: $(TARGETS)
tests: $(TARGETS_TEST)

chain_50.c: addchain
	./addchain -p 50 > $@
chain_50_100.c: addchain
	./addchain -p 100 50 > $@
chain_100_150.c: addchain
	./addchain -p 150 100 > $@
chain_150_200.c: addchain
	./addchain -p 200 150 > $@
chain_200_300.c: addchain
	./addchain -p 300 200 > $@
chain_300_400.c: addchain
	./addchain -p 400 300 > $@
chain_400_500.c: addchain
	./addchain -p 500 400 > $@

bc_50.h: addchain
	./addchain -pb 50 | head -n 1 > $@
bc_100.h: addchain
	./addchain -pb 100 | head -n 1 > $@
bc_150.h: addchain
	./addchain -pb 150 | head -n 1 > $@
bc_200.h: addchain
	./addchain -pb 200 | head -n 1 > $@
bc_300.h: addchain
	./addchain -pb 300 | head -n 1 > $@
bc_400.h: addchain
	./addchain -pb 400 | head -n 1 > $@
bc_500.h: addchain
	./addchain -pb 500 | head -n 1 > $@

pm1.o: pm1.c pm1.h stage2.h ../../utils/modredc_ul.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

pp1.o: pp1.c pp1.h chain_50.c chain_50_100.c chain_100_150.c chain_150_200.c \
        chain_200_300.c chain_300_400.c chain_400_500.c ../../utils/modredc_ul.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

ecm.o: ecm.c ecm.h ../../utils/modredc_ul.h  ../../utils/modredc_ul.c bc_50.h bc_100.h \
        bc_150.h bc_200.h bc_300.h bc_400.h bc_500.h
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

test_pm1s2: test_pm1s2.c pm1.o stage2.o prac_bc.o ../../utils/modredc_ul.h ../../utils/modredc_ul.c ../../utils/getprime.o
	$(CC) $(CFLAGS) $(INCLUDE) $< ../../utils/modredc_ul.c ../../utils/getprime.o pm1.o stage2.o prac_bc.o -o $@ -lgmp
clean:
	rm -f $(TARGETS) 
