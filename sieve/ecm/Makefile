GMPBUILD=
CFLAGS=-g -O2 -DNDEBUG -Wall -W -Wmissing-prototypes -Wconversion -pedantic
CC=gcc
TARGETS=find_primes_hit addchain

all: $(TARGETS)

find_primes_hit: find_primes_hit.c ../../utils/mod_ul.h  ../../utils/mod_ul.c
	$(CC) $(CFLAGS) $(GMPBUILD) -I ../../utils/ $< -lgmp ../../utils/mod_ul.c -lm -o $@

addchain: addchain.c ../../utils/mod_ul.h ../../utils/mod_ul.c
	$(CC) $(CFLAGS) -I ../../utils/ $< ../../utils/mod_ul.c -lm -o $@

clean:
	rm $(TARGETS) 

