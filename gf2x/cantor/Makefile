FLAGS=-O4 -DNDEBUG 
#FLAGS=-g


SRC_CANTOR=cantor128.c mpfq_2_128.h.32 mpfq_2_128.h.64 mpfq_2_128.c.32 mpfq_2_128.c.64

cantor128.o: $(SRC_CANTOR)
	gcc $(FLAGS) -c -msse2 cantor128.c mpfq_2_128.c

bench: cantor128.o bench.c
	gcc $(FLAGS) -o bench cantor128.o bench.c -lgmp

bench_notrunc: cantor128.o bench.c
	gcc $(FLAGS) -DMULCANTOR=mulCantor128Plain -o bench_notrunc cantor128.o bench.c -lgmp

cantor: $(SRC_CANTOR)
	gcc $(FLAGS) -DWANT_MAIN -msse2 -o cantor cantor128.c mpfq_2_128.c -lgmp

clean: 
	rm -f cantor bench cantor128.o mpfq_2_128.o bench_notrunc

all: bench bench_notrunc cantor
