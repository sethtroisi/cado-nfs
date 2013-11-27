BIN = sieve makefb normstats checkrels testirred descent_init
LIB = fppol

fppol_OBJ = f2pol.o f3pol.o fppol.o fppol_test.o fppol_set.o  \
	    fppol_arith.o fppol_mul.o fppol_div.o fppol_mod.o \
	    fppol_gcd.o fppol_io.o fppol_facttools.o fppol_internal.o \
	    fppol_conv.o polyfactor.o

sieve_OBJ = $(fppol_OBJ) ffspol.o fqpol.o ijvec.o params.o fb.o qlat.o    \
	    buckets.o norm.o latsieve.o polyfactor.o smoothness.o sieve.o \
	    memusage.o timing.o

makefb_OBJ = $(fppol_OBJ) ffspol.o fqpol.o polyfactor.o smoothness.o \
	     memusage.o timing.o params.o

normstats_OBJ = $(fppol_OBJ) ffspol.o fqpol.o ijvec.o params.o \
	        norm.o normstats.o memusage.o timing.o 

checkrels_OBJ = $(fppol_OBJ) ffspol.o memusage.o timing.o params.o norm.o 

freerels_OBJ = $(fppol_OBJ) ffspol.o memusage.o timing.o fqpol.o \
	       polyfactor.o

testirred_OBJ = $(fppol_OBJ) ffspol.o memusage.o timing.o params.o fqpol.o \
		polyfactor.o

descent_init_OBJ =  $(fppol_OBJ) ffspol.o params.o polyfactor.o smoothness.o \
	    memusage.o timing.o

# The binaries goes in CADO-NFS default build dirirectory.
# Can be overwritten in Makefile.local
HOST=$(shell hostname)
BUILDDIR = ../build/$(HOST)/ffs

# Use gf2x from CADO-NFS. Can also be defined or overwritten in Makefile.local
# GF2XPREFIX = ../build/$(HOST)/gf2x
OPTFLAGS  = 
WARNFLAGS = -Wall -Wextra -pedantic

FFSREV = `./version.sh`

CPPFLAGS    = -DFFSREV=\"$(FFSREV)\" -DSCALE=0 -I../utils -I..
CFLAGS      = $(OPTFLAGS) $(WARNFLAGS) -std=c99   -g
TARGET_ARCH = -mtune=native -march=native
LDFLAGS     = 
LDLIBS      = -lm -lpthread

#begin for freerels
#BIN += freerels
#LDFLAGS  += -L$(BUILDDIR)/../utils #for freerels
#LDLIBS   += -lutils -lm -lgmp #for freerels
#end for freerels

# Override variables in here!
# See Makefile.local.example for a small documentation about compile
# options.
-include Makefile.local

# Use GF2X?
ifdef GF2XPREFIX
  CPPFLAGS += -I$(GF2XPREFIX)/include -I$(GF2XPREFIX) -I../gf2x -I. -DHAVE_GF2X
  LDFLAGS  += -L$(GF2XPREFIX)/lib -L$(GF2XPREFIX)/.libs
  LDLIBS   += -lgf2x
endif

SRC = $(wildcard *.c)
INC = $(wildcard *.h)
OBJ = $(SRC:.c=.o)

FIELDS = f2 f3
f2_DEF = USE_F2
f3_DEF = USE_F3

.PHONY: all $(FIELDS) clean dirclean

all: $(FIELDS)

define BIN_template
  $(1)_$(2)_OBJ = $$(filter $$($(1)_OBJ), \
                     $$(addprefix $$($(1)_BUILDDIR)/, $$($(2)_OBJ)))
  $(1)_$(2)_BIN = $$($(1)_BUILDDIR)/$(2)

  $(1): $$($(1)_$(2)_BIN)

  $$($(1)_$(2)_BIN): $$($(1)_$(2)_OBJ)

  clean::
	rm -f $$($(1)_$(2)_BIN)

endef

define LIB_template
  $(1)_$(2)_OBJ = $$(filter $$($(1)_OBJ), \
                     $$(addprefix $$($(1)_BUILDDIR)/, $$($(2)_OBJ)))
  $(1)_$(2)_LIB = $$($(1)_BUILDDIR)/lib$(2)-$(1).a

  $(1): $$($(1)_$(2)_LIB)

  $$($(1)_$(2)_LIB): $$($(1)_$(2)_OBJ)
	$$(AR) rcs $$@ $$^

  clean::
	rm -f $$($(1)_$(2)_LIB)

endef

define FIELD_template
  $(1)_BUILDDIR = $$(BUILDDIR)/$(1)
  $(1)_DEP      = $$($(1)_BUILDDIR)/.depend
  $(1)_OBJ      = $$(addprefix $$($(1)_BUILDDIR)/,$$(filter-out $$(foreach f,\
                     $$(filter-out $(1),$$(FIELDS)),$$(f)pol.o), $$(OBJ)))
  $(1)_CPPFLAGS = -D$$($(1)_DEF)

  $$(foreach bin,$$(BIN),$$(eval $$(call BIN_template,$(1),$$(bin))))

  $$(foreach lib,$$(LIB),$$(eval $$(call LIB_template,$(1),$$(lib))))

  $$($(1)_BUILDDIR):
	mkdir -p $$@

  $$($(1)_OBJ): $$($(1)_BUILDDIR)/%.o: Makefile | $$($(1)_BUILDDIR)
	$$(CC) $$(CFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -c -o $$@ $$*.c

  $$($(1)_DEP): $$(SRC) $$(INC) Makefile | $$($(1)_BUILDDIR)
	$$(CC) $$(CFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -MM $$(SRC) | sed -e "s+^.*\.o:+$$($(1)_BUILDDIR)/&+" >$$@

  clean::
	rm -f $$($(1)_DEP) $$($(1)_OBJ) # $$($(1)_BIN)
	[ -z "`find $$($(1)_BUILDDIR) -maxdepth 0 -empty`" ] \
	  || rmdir $$($(1)_BUILDDIR)

  dirclean::
	rm -rf $$($(1)_BUILDDIR)

  -include $$($(1)_DEP)

endef

$(foreach f,$(FIELDS),$(eval $(call FIELD_template,$(f))))

clean dirclean::
	[ -z "`find $(BUILDDIR) -maxdepth 0 -empty`" ] || rmdir $(BUILDDIR)
