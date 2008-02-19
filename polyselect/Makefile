include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= -L $(GMP)/lib -lgmp
GMPINCLUDE ?= $(GMP)/include/

OPTCFLAGS ?= -O2 -g -W -Wall

CFLAGS=$(OPTCFLAGS) -std=c99 -I$(GMPINCLUDE) -DREV=\"`svnversion .`\"
LDFLAGS=../utils/libutils.a $(GMPLIB) -lm

polyselect: polyselect.c ../cado.h
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

clean:
	-rm -f polyselect *~ *.o
