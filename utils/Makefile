include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPINCLUDE ?= $(GMP)/include/
OPTCFLAGS ?= -g -O2 -DNDEBUG -W -Wall

CFLAGS = $(OPTCFLAGS)

names:=polyfile relation timing long_poly getprime mod_ul modul_poly gmp_aux \
       fpoly
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a cache cachesize check_rels

%.o: %.c
	$(CC) $(CFLAGS) -c $< -I.. -I$(GMPINCLUDE)

cache: cache.c timing.o
	$(CC) $(CFLAGS) $< timing.c -I.. -o $@

cachesize: cachesize.c
	$(CC) $(CFLAGS) $< -o $@

check_rels: check_rels.o libutils.a
	$(CC) $(CFLAGS) $^ -o $@ -lgmp

libutils.a: $(obj)
	ar r libutils.a $(obj)


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

clean:
	-rm -f *~ *.bak $(obj) libutils.a cache cachesize


