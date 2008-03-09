include Makefile.local

# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= -L$(GMP)/lib -lgmp
GMPINCLUDE ?= $(GMP)/include/

OPTCFLAGS ?= -O2 -DNDEBUG -g -W -Wall

CFLAGS=$(OPTCFLAGS) -std=c99 -I$(GMPINCLUDE) -DREV=\"`svnversion .`\"
LDFLAGS=-L../utils $(GMPLIB) -lutils -lm

all: polyselect

.c.o:
	$(CC) $(CFLAGS) -I.. $< -c -o $@

polyselect.o: polyselect.c ../cado.h

polyselect: polyselect.o
	$(CC) $^ $(LDFLAGS) -o $@


Makefile.local:
	echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../Makefile.local .

clean:
	-rm -f polyselect *~ *.o
