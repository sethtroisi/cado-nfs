GMP=

CFLAGS=-O2 -g -I$(GMP)/include -DREV=\"`svnversion`\"
LDFLAGS=$(GMP)/lib/libgmp.a -lm

polyselect: polyselect.c ../cado.h
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

clean:
	rm polyselect *~
