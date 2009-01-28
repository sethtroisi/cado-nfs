
CC:=gcc

CFLAGS:=-O4 -DNDEBUG -msse2
ALWAYS_CFLAGS:=-W -Wall -std=c99

# Start with defaults that should never be overridden
# These flags have late expansion, on purpose.
MY_CFLAGS=$(ALWAYS_CFLAGS)
MY_LFLAGS=

# Add the possibly overridden variables.
MY_CFLAGS+=$(CFLAGS)
MY_LFLAGS+=$(LDFLAGS)

# Support a library (GMP)
# Within Makefile.local, or on the command line, either put
# GMP=/some/prefix/, and/or GMP_LIBDIR=/some/directory and
# GMP_INCDIR=/some/directory
# Note that the proper -l must still be specified per-target.
ifneq ($(GMP),)
GMP_LIBDIR ?= $(GMP)/lib
GMP_INCDIR ?= $(GMP)/include
endif
ifneq ($(GMP_LIBDIR),)
GMP_LFLAGS ?= $(addprefix -L,$(GMP_LIBDIR))
endif
ifneq ($(GMP_INCDIR),)
GMP_CFLAGS ?= $(addprefix -I,$(GMP_INCDIR))
GMP_CXXFLAGS ?= $(addprefix -I,$(GMP_INCDIR))
endif
MY_CFLAGS+=$(GMP_CFLAGS)
MY_LFLAGS+=$(GMP_LFLAGS)

all: cantor bench cantor_gm bench_gm

.c.o:
	$(CC) $(MY_CFLAGS) -o $@ -c $<

cantor128.o: mpfq_2_128.h.32 mpfq_2_128.h.64
cantor128_gm.o: mpfq_2_128.h.32 mpfq_2_128.h.64
mpfq_2_128.o: mpfq_2_128.c.32 mpfq_2_128.c.64
mpfq_2_128.o: mpfq_2_128.h.32 mpfq_2_128.h.64

bench: bench.o cantor128.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
cantor: cantor_main.o cantor128.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
bench_gm: bench.o cantor128_gm.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
cantor_gm: cantor_main.o cantor128_gm.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp
test: cantor
	./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`
test_gm: cantor_gm
	CANTOR=./cantor_gm ./check_cantor.sh `seq 60 70` `seq 120 140` `seq 250 270`

clean: 
	-rm -f cantor bench *.o *~
