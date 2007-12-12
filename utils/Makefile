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

clean:
	-rm -f *~ *.bak $(obj) libutils.a


