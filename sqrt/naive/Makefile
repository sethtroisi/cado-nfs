#libgmp=$(GMP)/lib/libgmp.a
libgmp=-lgmp

CFLAGS=-g
LDFLAGS=../../utils/libutils.a $(libgmp)

algsqrt: algsqrt.c poly.o
	gcc $(CFLAGS) -I../.. -o $@ $< poly.o $(LDFLAGS) 

poly.o: poly.c
	gcc $(CLFAGS) -c $< -o $@

clean:
	rm -f poly.o algsqrt

