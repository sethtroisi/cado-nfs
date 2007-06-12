CFLAGS=-O2 -g
LDFLAGS=/usr/lib64/libgmp.a

procrels: procrels.c
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

linalg: linalg.c gauss.o
	gcc $(CFLAGS) gauss.o -I.. $< $(LDFLAGS) -o $@

gauss.o: gauss.c
	gcc -c -O4 -DNDEBUG -funroll-loops -DMULTI_ROW=1 -DNO_MAIN $< -o $@

ratsqrt: ratsqrt.c
	gcc $(CFLAGS) $< $(LDFLAGS) -o $@

characters: characters.c
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

all: procrels linalg ratsqrt characters

clean:
	\rm procrels gauss.o linalg ratsqrt characters
