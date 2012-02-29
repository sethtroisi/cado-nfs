SRC = $(wildcard *.c)
INC = $(wildcard *.h)
OBJ = $(SRC:.c=.o)
BIN = sieve

#OPTFLAGS = -O3
OPTFLAGS =
WARNFLAGS = -Wall -Wextra -pedantic

CPPFLAGS    =
CFLAGS      = $(OPTFLAGS) $(WARNFLAGS) -std=c99 -g # -DTRACE_POS=49755
CXXFLAGS    = $(OPTFLAGS) $(WARNFLAGS) -std=c++0x -g
TARGET_ARCH = -mtune=native -march=native 
LDFLAGS     = 
LDLIBS      = -lntl -lstdc++

FIELDS = f2 f3
f2_DEF = USE_F2
f3_DEF = USE_F3

.PHONY: all clean

define BUILD_template
  $(1)_OBJ      = $$(filter-out $$(foreach f,$$(FIELDS),$$(f)pol.$(1).o), \
                     $$(OBJ:.o=.$(1).o)) $(1)pol.$(1).o
  $(1)_BIN      = $$(BIN:=.$(1))
  $(1)_CPPFLAGS = -D$$($(1)_DEF)

  all: $$($(1)_BIN)

  $$(filter fppol%.$(1).o, $$($(1)_OBJ)): $(1)pol.h

  $$($(1)_OBJ): %.$(1).o: Makefile
	$$(CC) $$(CFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -c -o $$@ $$*.c

  cofactor.$(1).o: Makefile
	$$(CXX) $$(CXXFLAGS) $$(CPPFLAGS) $$($(1)_CPPFLAGS) $$(TARGET_ARCH) \
	  -c -o $$@ cofactor.cc

  $$($(1)_BIN): $$($(1)_OBJ) cofactor.$(1).o

  clean::
	rm -f $$($(1)_OBJ) cofactor.$(1).o $$($(1)_BIN)
endef

$(foreach f,$(FIELDS),$(eval $(call BUILD_template,$(f))))

clean::
	rm -f .depend

.depend: $(SRC) $(INC) cofactor.cc cofactor.hh Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) -MM $(SRC) cofactor.cc           \
	  | sed "s+\(.*\)\.o+$(foreach f,$(FIELDS),\1.$(f).o )+" >$@

-include .depend
