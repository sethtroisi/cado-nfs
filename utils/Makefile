TOP:=..
include $(TOP)/Makefile.common

# All parts of libutils.
names:=	cado_poly relation timing plain_poly \
	getprime mod_ul modredc_ul modredc_15ul modul_poly gmp_aux \
	fpoly rootfinder params gcd_int64 gcd_uint64 discriminant	\
	random_generation mpz_array

src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

bins:=cache cachesize check_rels
all: libutils.a $(bins)

cache: cache.o timing.o
	$(CC) $(MY_LFLAGS) $^ -o $@

cachesize: cachesize.o
	$(CC) $(MY_LFLAGS) $^ -o $@

check_rels: check_rels.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -L. -lutils -lgmp

test-rootfinder: test-rootfinder.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -L. -lutils -lgmp

test-rootfinder: libutils.a
check_rels: libutils.a
mod_ul.o: mod_ul.h mod_ul_common.c ularith.h
modredc_ul.o: modredc_ul.h mod_ul_common.c ularith.h
modredc_15ul.o: modredc_15ul.h ularith.h

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	-rm -f *.bak *.o libutils.a $(bins)

test: test-rootfinder
	./test-rootfinder-script.sh < testbase-rootfinder.quick



# ALL source files must come here. Using $(wildcard *.c) is discouraged
# for obvious reasons.
sources:=       \
    cache.c cachesize.c cado_poly.c check_rels.c discriminant.c fpoly.c	\
    gcd_int64.c gcd_uint64.c getprime.c gmp_aux.c mod_ul.c mod_ul_common.c \
    modredc_15ul.c modredc_ul.c modul_poly.c params.c plain_poly.c \
    random_generation.c relation.c rootfinder.c test-rootfinder.c timing.c

dep: $(patsubst %,.deps/%.d,$(sources))
-include $(patsubst %,.deps/%.d,$(sources))


