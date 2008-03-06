# Usage: make <target>
#
# Paths and anything specific to a host go in Makefile.local. It can specify
# the full path/filename to the ECM libarary in ECMLIB, and full path to ECM
# header file in ECMINCLUDE. If not specified, defaults are derived by the
# variable assignments below.
# Same thing for GMPLIB and GMPINCLUDE for the GMP library.

include Makefile.local

# Set default filenames/paths for GMP if not already set.
GMP ?= /usr/local
GMPLIB ?= $(GMP)/lib/libgmp.a
GMPINCLUDE ?= $(GMP)/include

# Set USE_TIFA_FROM_SVN to one to use TIFA headers/library from local SVN
# repository. Of course TIFA has to be compiled.
USE_TIFA_FROM_SVN ?= 1

# Set default filenames/paths for TIFA if not already set.
ifeq ($(USE_TIFA_FROM_SVN), 1)
    TIFA = ../tifa
    TIFALIB ?= $(TIFA)/build/lib/libtifa.a
    TIFAINCFLAGS = -I$(TIFA) -I$(TIFA)/lib -I$(TIFA)/lib/algo/include \
                   -I$(TIFA)/lib/data/include -I$(TIFA)/lib/utils/include
else
    TIFA ?= /usr/local
    TIFALIB ?= $(TIFA)/lib/libtifa.a
    TIFAINCFLAGS ?= $(TIFA)/include/tifa
endif

# Set default filenames/paths for utils if not already set.
UTILSLIB ?= ../../utils/libutils.a
UTILSINCLUDE ?= ../../utils

OPTFLAGS ?= -std=c99 -O2 -DNDEBUG -g -fkeep-inline-functions -Wall -W
CFLAGS=$(OPTFLAGS) -DREV=\"`svnversion .`\"
LDFLAGS ?= -lm
INCFLAGS ?= -I. -I.. -I../.. -I$(GMPINCLUDE) $(TIFAINCFLAGS) -I$(UTILSINCLUDE)

all: checknorms checknorms.2lp

Makefile.local:
	@echo "Makefile.local does not exist. Making a link to parent dir's one."
	ln -s ../../Makefile.local .

checknorms: checknorms.c checknorms.h Makefile.local $(TIFALIB) $(UTILSLIB)
	gcc $(INCFLAGS) $(CFLAGS) -o $@ $< $(UTILSLIB) $(TIFALIB) $(GMPLIB) $(LDFLAGS)

checknorms.2lp: checknorms.2lp.c Makefile.local $(TIFALIB) $(UTILSLIB)
	gcc $(INCFLAGS) $(CFLAGS) -o $@ $< $(UTILSLIB) $(TIFALIB) $(GMPLIB) $(LDFLAGS)

config.h: config_h.in
	cp config_h.in config.h

clean:
	rm -f checknorms
	rm -f checknorms.2lp

