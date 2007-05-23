CFLAGS=-O2 -g
LDFLAGS=/usr/lib64/libgmp.a

procrels: procrels.c
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

clean:
	\rm procrels
