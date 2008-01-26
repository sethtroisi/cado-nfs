include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/

OPTCFLAGS ?= -O2 -g -W -Wall

CFLAGS=$(OPTCFLAGS) -I $(GMPINCLUDE) -DREV=\"`svnversion .`\"

LDFLAGS=../../utils/libutils.a $(GMPLIB)

all: fmalgsqrt algsqrt

Makefile.local:
	@echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../../Makefile.local .

.c.o:
	gcc $(CFLAGS) -I../.. -o $@ -c $<

algsqrt: algsqrt.o poly.o
	gcc -o $@ $^ $(LDFLAGS) 

fmalgsqrt: fmalgsqrt.o poly.o
	gcc -o $@ $^ $(LDFLAGS)

ratsqrt: ratsqrt.o
	gcc -o $@ $^ $(LDFLAGS)

clean:
	-rm -f poly algsqrt ratsqrt fmalgsqrt *.o *~

