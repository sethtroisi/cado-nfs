# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPINCLUDE ?= $(GMP)/include/
CFLAGS=-g -O2 -W -Wall

targets:=convert_poly convert_rels

all: $(targets)

.SUFFIXES: .c

.c:
	gcc $(CFLAGS) $< -I.. -I$(GMPINCLUDE) $(GMP)/lib/libgmp.a -o $@

clean:
	\rm -f $(targets) *~
