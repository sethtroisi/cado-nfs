GMP?=/usr/local/
GMPINC?=$(GMP)/include
GMPLIB?=$(GMP)/lib

CC:=gcc
CFLAGS:=-W -Wall -O4 -DNDEBUG -std=c99 -I$(GMPINC)

all: cantor bench

cantor128.o: mpfq_2_128.h.32 mpfq_2_128.h.64
mpfq_2_128.o: mpfq_2_128.c.32 mpfq_2_128.c.64
mpfq_2_128.o: mpfq_2_128.h.32 mpfq_2_128.h.64
cantor_main.o:
c.o:
	$(CC) $(CFLAGS) -m sse2 -o $@ -c $<
bench: bench.o cantor128.o
	$(CC) $(FLAGS) -o $@ $^ -L$(GMPLIB) -lgmp
cantor: cantor_main.o cantor128.o
	$(CC) $(FLAGS) -msse2 -o $@ $^ -L$(GMPLIB) -lgmp
test: cantor
	./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
clean: 
	-rm -f cantor bench *.o *~
