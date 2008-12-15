
# Sorry, this makefile is kind of ugly, so pending a choice for an
# appropriate building tool, we'll simply avoid make -j for this one.
.NOTPARALLEL:

# Default flags
CFLAGS:=-O2 -DNDEBUG
CXXFLAGS:=$(CFLAGS)

ALWAYS_CFLAGS:=-std=c99 -g -W -Wall
ALWAYS_CXXFLAGS:=-g -W -Wall

# Start with defaults that should never be overridden
MY_CFLAGS=$(ALWAYS_CFLAGS)
MY_CXXFLAGS=$(ALWAYS_CXXFLAGS)
MY_LFLAGS=-L.

# Add the possibly overridden variables.
MY_CFLAGS+=$(CFLAGS)
MY_CXXFLAGS+=$(CXXFLAGS)
MY_LFLAGS+=$(LDFLAGS)
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
MY_CXXFLAGS+=$(GMP_CXXFLAGS)
MY_LFLAGS+=$(GMP_LFLAGS)

.c.o:
	$(CC) $(MY_CFLAGS) $< -c -o $@

ifeq ($(shell uname -m),i686)
w:=32
CFLAGS+=-msse2
CXXFLAGS+=-msse2
endif

ifeq ($(shell uname -m),x86_64)
w:=64
endif
ifndef w
$(error "Define w")
endif

.PHONY: all tune1

all: tune1 libgf2x.a tunetoom

LIBGF2X_FILES:=mul-toom.c mul-gf2x.c mulfft-bit.c

libgf2x.a: $(patsubst %.c,%.o,$(LIBGF2X_FILES))
	rm -f $@
	ar qcv $@ $^

gen_bb_mul_code:	gen_bb_mul_code.c
	$(CC) $(MY_CFLAGS) $< -o $@

mul-inlines-k%.c:	gen_bb_mul_code
	@./gen_bb_mul_code $(w) $* > $@

tuneup-k%: mul-inlines-k%.c mul-small.c
	@$(CC) $(MY_CFLAGS) \
		-DTUNE_MUL1 -DWORDSIZE=$(w)				\
		-DINLINES_FILE=\"mul-inlines-k$*.c\"			\
		-Wno-unused						\
		tuneup.c -o $@

tune1-k%: tuneup-k%
	@echo $* " : "`./tuneup-k$* -nokara 1`

tune1_files:= $(patsubst %,mul-inlines-k%.c,1 2 3 4 5 6 7 8)
tune1_progs:= $(addprefix tuneup-k,1 2 3 4 5 6 7 8)
tune1_actions:= $(addprefix tune1-k,1 2 3 4 5 6 7 8)

mul-inlines.c:
	-rm -f tune1.log
	$(MAKE) $(tune1_files) $(tune1_progs)
	$(MAKE) $(tune1_actions) | tee tune1.log
	./best_k.sh $(w) tune1.log

tune1:	mul-inlines.c

TUNETOOM_FILES:=$(LIBGF2X_FILES) tunetoom.c timing.c tuning-common.c
tunetoom: $(patsubst %.c,%.o,$(TUNETOOM_FILES))
	$(CC) $(MY_CFLAGS) -o $@ $^
UPDATE_THRESHOLDS_FILES:=update-thresholds.c replace.c
update-thresholds: $(patsubst %.c,%.o,$(UPDATE_THRESHOLDS_FILES))
	$(CC) $(MY_CFLAGS) -o $@ $^
TUNEFFT_FILES:=$(LIBGF2X_FILES) tunefft.c timing.c mulfft-bit.c tuning-common.c
tunefft: $(patsubst %.c,%.o,$(TUNEFFT_FILES))
	$(CC) $(MY_CFLAGS) -o $@ $^
BENCH_FILES:=bench.c
bench: $(patsubst %.c,%.o,$(BENCH_FILES))
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgf2x

factor.o: factor.cpp
	$(CXX) $(MY_CXXFLAGS) -o $@ -c $<

factor: factor.o
	$(CXX) $(MY_LFLAGS) -o $@ $< -lntl -lgmp -lgf2x

factor: libgf2x.a

.PHONY: do-tune-toom do-tune-fft

recover-tune: update-thresholds
	cat tunetoom.res tunefft.res | ./update-thresholds -o thresholds.h

do-tune-toom: tunetoom update-thresholds
	./tunetoom -s 1.05 2048 -o tunetoom.res
	./update-thresholds -o thresholds.h < tunetoom.res

do-tune-fft: tunefft update-thresholds
	./tunefft -s 1.02 8000000 -o tunefft.res
	./update-thresholds -o thresholds.h < tunefft.res


tune-all: update-thresholds
	$(MAKE) tune1
	$(MAKE) do-tune-toom
	$(MAKE) do-tune-fft

do-check%.in: check%.in factor
	tail -n +2 $< | \
	./factor `head -1 $<` > check$*.res
	@diff -b check$*.res check$*.exp	&& echo "Test $* passed"

check: $(addprefix do-,$(wildcard check*.in))

clean:
	-rm -f gen_bb_mul_code mul-inlines-k*.c tuneup-k* tune1-k*
	-rm -f mul-inlines.c *~ tune1.log
	-rm -f bench update-thresholds tunetoom tunefft libgf2x.a
	-rm -f check*.res factor
	-rm -f *.o

gf2x.h:	mul-inlines.c mul-small.c mul-basecase.c thresholds.h
mul-kara.o: gf2x.h thresholds.h
mul-tc3.o: gf2x.h thresholds.h
mul-tc3u.o: gf2x.h thresholds.h
mul-tc3w.o: gf2x.h thresholds.h
mul-tc4.o: gf2x.h thresholds.h
mul-toom.o: gf2x.h thresholds.h
mul-small.o: mul2t.c thresholds.h
timing.o: timing.h thresholds.h
toomspace.o: gf2x.h thresholds.h
tunetoom.o: gf2x.h replace.h thresholds.h
tunefft.o: gf2x.h replace.h thresholds.h
test.o: gf2x.h thresholds.h
replace.o: replace.h
factor.o: HalfGCD.c gf2x.h thresholds.h
mulfft-bit.o: gf2x.h thresholds.h
mul-gf2x.o: gf2x.h thresholds.h
