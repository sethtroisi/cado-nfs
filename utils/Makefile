include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPINCLUDE ?= $(GMP)/include/
OPTCFLAGS ?= -g -O2 -W -Wall

CFLAGS = $(OPTCFLAGS)

names:=polyfile relation timing long_poly getprime mod_ul gmp_aux
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a cache

%.o: %.c
	gcc $(CFLAGS) -c $< -I.. -I$(GMPINCLUDE)

cache: cache.c timing.o
	gcc $(CFLAGS) $< timing.c -I.. -o $@

libutils.a: $(obj)
	ar r libutils.a $(obj)


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

clean:
	-rm -f *~ *.bak $(obj) libutils.a cache


