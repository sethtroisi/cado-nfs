SRC = $(wildcard *.c)
INC = $(wildcard *.h)
OBJ = $(SRC:.c=.o)
BIN = sieve

BUILDDIR = ./build

OPTFLAGS  = 
WARNFLAGS = -Wall -Wextra -pedantic

FFSREV = `./version.sh`

CPPFLAGS    = -DFFSREV=\"$(FFSREV)\"
CFLAGS      = $(OPTFLAGS) $(WARNFLAGS) -std=c99   -g
CXXFLAGS    = $(OPTFLAGS) $(WARNFLAGS) -std=c++0x -g
TARGET_ARCH = -mtune=native -march=native 
LDFLAGS     = 
LDLIBS      = -lntl -lstdc++

# Override variables in here!
-include Makefile.local

# Use GF2X?
ifdef GF2XPREFIX
  CPPFLAGS    += -I$(GF2XPREFIX)/include -DHAVE_GF2X
  TARGET_ARCH += -mpclmul
  LDFLAGS     += -L$(GF2XPREFIX)/lib
  LDLIBS      += -lgf2x
endif

FIELDS = f2 f3
f2_DEF = USE_F2
f3_DEF = USE_F3

.PHONY: all $(FIELDS) clean dirclean

all: $(FIELDS)

define BUILD_template
  $(1)_BUILDDIR = $$(BUILDDIR)/$(1)
  $(1)_DEP      = $$($(1)_BUILDDIR)/.depend
  $(1)_OBJ      = $$(addprefix $$($(1)_BUILDDIR)/,$$(filter-out $$(foreach f,\
                     $$(filter-out $(1),$$(FIELDS)),$$(f)pol.o), $$(OBJ)))
  $(1)_BIN      = $$(addprefix $$($(1)_BUILDDIR)/,$$(BIN))
  $(1)_CPPFLAGS = -D$$($(1)_DEF)

  $(1): $$($(1)_BIN)

  $$($(1)_BUILDDIR):
	mkdir -p $$@

  $$($(1)_OBJ): $$($(1)_BUILDDIR)/%.o: Makefile | $$($(1)_BUILDDIR)
	$$(CC) $$(CFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -c -o $$@ $$*.c

  $$($(1)_BUILDDIR)/cofactor.o: Makefile | $$($(1)_BUILDDIR)
	$$(CXX) $$(CXXFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -c -o $$@ cofactor.cc

  $$($(1)_BIN): $$($(1)_OBJ) $$($(1)_BUILDDIR)/cofactor.o

  $$($(1)_DEP): $$(SRC) $$(INC) cofactor.cc cofactor.hh Makefile | $$($(1)_BUILDDIR)
	$$(CC) $$(CFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -MM $$(SRC) cofactor.cc | sed -e "s+^.*\.o:+$$($(1)_BUILDDIR)/&+" >$$@

  clean::
	rm -f $$($(1)_DEP) $$($(1)_OBJ) $$($(1)_BUILDDIR)/cofactor.o $$($(1)_BIN)
	[ -z "`find $$($(1)_BUILDDIR) -maxdepth 0 -empty`" ] \
	  || rmdir $$($(1)_BUILDDIR)

  dirclean::
	rm -rf $$($(1)_BUILDDIR)

  -include $$($(1)_DEP)

endef

$(foreach f,$(FIELDS),$(eval $(call BUILD_template,$(f))))

clean dirclean::
	[ -z "`find $(BUILDDIR) -maxdepth 0 -empty`" ] || rmdir $(BUILDDIR)
