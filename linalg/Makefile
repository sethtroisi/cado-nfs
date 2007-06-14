CFLAGS=-O2 -g -DREV=\"`svnversion`\"
LDFLAGS=$(GMP)/lib/libgmp.a

procrels: procrels.c
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

linalg: linalg.c gauss.o
	gcc $(CFLAGS) gauss.o -I.. $< $(LDFLAGS) -o $@

gauss.o: gauss.c
	gcc -c -g -DMULTI_ROW=1 -DNO_MAIN -I.. $< -o $@

ratsqrt: ratsqrt.c
	gcc $(CFLAGS) $< $(LDFLAGS) -o $@

characters: characters.c gauss.o
	gcc $(CFLAGS) -I.. $< gauss.o $(LDFLAGS) -o $@

all: procrels linalg ratsqrt characters

clean:
	\rm procrels gauss.o linalg ratsqrt characters
