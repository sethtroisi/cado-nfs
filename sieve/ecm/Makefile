GMPBUILD=
CFLAGS=-g -O2 -DNDEBUG -Wall -W -Wmissing-prototypes -Wconversion -pedantic
CC=gcc

find_primes_hit: find_primes_hit.c
	$(CC) $(CFLAGS) $(GMPBUILD) $< -lgmp -lm -o $@

clean:
	rm find_primes_hit

