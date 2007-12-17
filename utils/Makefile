include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPINCLUDE ?= $(GMP)/include/
CFLAGS=-g -O2 -W -Wall

names:=polyfile relation timing long_poly getprime mod_ul
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a

%.o: %.c
	gcc $(CFLAGS) -c $< -I.. -I$(GMPINCLUDE)

libutils.a: $(obj)
	ar r libutils.a $(obj)

Makefile.local:
	echo "Makefile.local does not exist. It can be used to store paths \
	to libraries etc.  Making an empty file to avoid warnings in \
	future invocations of make."
	touch Makefile.local

clean:
	-rm -f *~ *.bak $(obj) libutils.a


