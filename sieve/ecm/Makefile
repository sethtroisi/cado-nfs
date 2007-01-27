GMPBUILD=
CFLAGS=-g -O2 -Wall -W -Wmissing-prototypes -pedantic
CC=gcc

find_primes_hit: find_primes_hit.c
	$(CC) $(CFLAGS) -I$(GMPBUILD) $< $(GMPBUILD)/.libs/libgmp.a -lm -o $@

clean:
	rm find_primes_hit

