GMPBUILD=
#NDEBUG=-DNDEBUG
CFLAGS=-g -O0 $(NDEBUG) -Wall -W -Wconversion -pedantic
INCLUDE=-I ../../utils/
CC=gcc
TARGETS=find_primes_hit addchain pp1.o

all: $(TARGETS)

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

pp1.o : pp1.c pp1.h chain_50.c chain_50_100.c chain_100_150.c chain_150_200.c \
        chain_200_300.c chain_300_400.c chain_400_500.c ../../utils/mod_ul.h
	$(CC) $(CFLAGS) $(INCLUDE) -c -o $@ $<

find_primes_hit: find_primes_hit.c ../../utils/mod_ul.h  ../../utils/mod_ul.c
	$(CC) $(CFLAGS) $(GMPBUILD) $(INCLUDE) $< -lgmp ../../utils/mod_ul.c -lm -o $@

addchain: addchain.c ../../utils/mod_ul.h ../../utils/mod_ul.c
	$(CC) $(CFLAGS) $(INCLUDE) $< ../../utils/mod_ul.c -lm -o $@

clean:
	rm $(TARGETS) 

