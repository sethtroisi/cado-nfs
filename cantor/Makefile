
CC:=gcc

all: cantor bench cantor_gm bench_gm cantor_gmt bench_gmt

ifdef   WITHIN_CADO
TOP:=..
include $(TOP)/Makefile.common
else
include ./Makefile.common.gf2x
endif

cantor128.o: mpfq_2_128.h.32 mpfq_2_128.h.64 cantor128.h
cantor128_gm.o: mpfq_2_128.h.32 mpfq_2_128.h.64 cantor128.h cantor128.c
	$(CC) $(MY_CFLAGS) -DCANTOR_GM -o $@ -c cantor128.c
cantor128_gmt.o: mpfq_2_128.h.32 mpfq_2_128.h.64 cantor128.h cantor128.c
	$(CC) $(MY_CFLAGS) -DCANTOR_GM -DCANTOR_GM_TRUNCATE -o $@ -c cantor128.c

mpfq_2_128.o: mpfq_2_128.c.32 mpfq_2_128.c.64
mpfq_2_128.o: mpfq_2_128.h.32 mpfq_2_128.h.64

bench: bench.o cantor128.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
cantor: cantor_main.o cantor128.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
bench_gm: bench.o cantor128_gm.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
cantor_gm: cantor_main.o cantor128_gm.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
bench_gmt: bench.o cantor128_gmt.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
cantor_gmt: cantor_main.o cantor128_gmt.o ; $(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

test: cantor
	./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
test_gm: cantor_gm
	CANTOR=./cantor_gm ./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
test_gmt: cantor_gmt
	CANTOR=./cantor_gmt ./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`

clean: 
	-rm -f cantor bench *.o *~
