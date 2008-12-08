include Makefile.local
# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/
CFLAGS ?= -g -O2 -DNDEBUG -W -Wall

targets:=convert_poly convert_rels contains_zerochar 

all: $(targets)

.SUFFIXES: .c

.c:
	gcc $(CFLAGS) $< -I.. -I$(GMPINCLUDE) $(GMPLIB) -o $@

# the dependence on ../linalg/gzip.o is a dirty hack, once we switch to cmake
# this can be done in a clean way, or better the gzip stuff should go in utils
convert_rels: convert_rels.c
	gcc $(CFLAGS) $< ../linalg/gzip.o -I.. -I$(GMPINCLUDE) $(GMPLIB) -o $@

clean:
	\rm -f $(targets) *~
