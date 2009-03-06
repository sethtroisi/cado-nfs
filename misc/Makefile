TOP:=..
include $(TOP)/Makefile.common

include Makefile.local
# Set default filenames/paths if not already set.
GMP ?= /usr/local/
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include/
CFLAGS ?= -g -O2 -DNDEBUG -W -Wall -I$(TOP)/utils

targets:=convert_poly convert_rels contains_zerochar 

all: $(targets)

.SUFFIXES: .c

.c:
	gcc $(CFLAGS) $< -I.. -I$(GMPINCLUDE) $(GMPLIB) -o $@

convert_rels: convert_rels.c
	gcc $(CFLAGS) $< -I.. -I$(GMPINCLUDE) $(GMPLIB) -o $@ -L$(TOP)/utils -lutils

clean:
	\rm -f $(targets) *~
