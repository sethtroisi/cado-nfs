include Makefile.local
# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/
CFLAGS ?= -g -O2 -W -Wall

targets:=convert_poly convert_rels

all: $(targets)

.SUFFIXES: .c

.c:
	gcc $(CFLAGS) $< -I.. -I$(GMPINCLUDE) $(GMPLIB) -o $@

clean:
	\rm -f $(targets) *~
