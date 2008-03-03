
# Sorry, this makefile is kind of ugly, so pending a choice for an
# appropriate building tool, we'll simply avoid make -j for this one.
.NOTPARALLEL:

CFLAGS:=-g -O2 -W -Wall -std=c99

ifeq ($(shell arch),i686)
w:=32
OFLAGS:=-msse2
endif

ifeq ($(shell arch),x86_64)
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
	gcc $(CFLAGS) $< -o $@

mul-inlines-k%.c:	gen_bb_mul_code
	@./gen_bb_mul_code $(w) $* > $@

tuneup-k%: mul-inlines-k%.c mul-small.c
	@gcc $(CFLAGS) \
		-DTUNE_MUL1 -DWORDSIZE=$(w)				\
		-DINLINES_FILE=\"mul-inlines-k$*.c\"			\
		-Wno-unused						\
		tuneup.c -o $@

tune1-k%: tuneup-k%
	@echo -n $* " : " ; ./tuneup-k$* -nokara 1

tune1_files:= $(patsubst %,mul-inlines-k%.c,$(shell seq 1 8))
tune1_progs:= $(addprefix tuneup-k,$(shell seq 1 8))
tune1_actions:= $(addprefix tune1-k,$(shell seq 1 8))

mul-inlines.c:
	-rm -f tune1.log
	$(MAKE) $(tune1_files) $(tune1_progs)
	$(MAKE) $(tune1_actions) | tee tune1.log
	./best_k.sh $(w) tune1.log

tune1:	mul-inlines.c

TUNETOOM_FILES:=$(LIBGF2X_FILES) tunetoom.c timing.c replace.c
tunetoom: $(patsubst %.c,%.o,$(TUNETOOM_FILES))
	gcc $(CFLAGS) -o $@ $^
TUNEFFT_FILES:=$(LIBGF2X_FILES) tunefft.c timing.c mulfft-bit.c replace.c
tunefft: $(patsubst %.c,%.o,$(TUNEFFT_FILES))
	gcc $(CFLAGS) -o $@ $^
TEST_FILES:=$(LIBGF2X_FILES) test.c timing.c
test: $(patsubst %.c,%.o,$(TEST_FILES))
	gcc $(CFLAGS) -o $@ $^

factor.o: factor.c
	g++ -g -o $@ -c $<

factor: factor.o libgf2x.a
	g++ -g factor.o -L. -lntl -lgmp -lgf2x -o $@

tune-all:
	$(MAKE) tune1
	$(MAKE) tunetoom && ./tunetoom 2048
	$(MAKE) tunefft && ./tunefft 100000

do-check%.in: check%.in factor
	tail -n +2 $< | \
	./factor `head -1 $<` > check$*.res
	@diff -b check$*.res check$*.exp	&& echo "Test $* passed"

check: $(addprefix do-,$(wildcard check*.in))

clean:
	-rm -f gen_bb_mul_code mul-inlines-k*.c tuneup-k* tune1-k*
	-rm -f mul-inlines.c *~ tune1.log
	-rm -f test tunetoom tunefft libgf2x.a
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
replace.o: replace.h thresholds.h
factor.o: HalfGCD.c gf2x.h thresholds.h
