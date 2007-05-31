GMP=

polyselect: polyselect.c ../cado.h
	gcc -I.. -O2 -g $< -o $@ $(GMP)/lib/libgmp.a -lm

clean:
	rm polyselect *~
