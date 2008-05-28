TOP:=..
include $(TOP)/Makefile.common

# All parts of libutils.
names:=	polyfile relation timing long_poly \
	getprime mod_ul modul_poly gmp_aux \
	fpoly

src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a cache cachesize check_rels

cache: cache.o timing.o
	$(CC) $(MY_LFLAGS) $^ -o $@

cachesize: cachesize.o
	$(CC) $(MY_LFLAGS) $^ -o $@

check_rels: check_rels.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -L. -lutils -lgmp

check_rels: libutils.a

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	-rm -f *~ *.bak *.o libutils.a cache cachesize


