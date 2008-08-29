TOP:=..
include $(TOP)/Makefile.common

# All parts of libutils.
names:=	cado_poly relation timing plain_poly \
	getprime mod_ul modredc_ul modredc_15ul modul_poly gmp_aux \
	fpoly rootfinder params

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

check_rels: libutils.a

mod_ul.o: mod_ul.h mod_ul_common.c

modredc_ul.o: modredc_ul.h mod_ul_common.c

modredc_15ul.o: modredc_15ul.h

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	-rm -f *.bak *.o libutils.a $(bins)
