.PHONY: clean all

TOP:=.
include $(TOP)/Makefile.common

export CC
export CXX
export SCONS
export GMP_LIBDIR
export GMP_INCDIR

# Pass down these overrides to scons
ifneq ($(CC),)
SCONS_FLAGS+=CC=$(CC)
endif
ifneq ($(GMP_LIBDIR),)
SCONS_FLAGS+=GMP_LIBDIR=$(GMP_LIBDIR)
endif
ifneq ($(GMP_INCDIR),)
SCONS_FLAGS+=GMP_INCDIR=$(GMP_INCDIR)
endif


all:
	$(MAKE) -C utils
	$(MAKE) -C polyselect
	(cd postsieve/tifa; $(SCONS) $(SCONS_FLAGS))
	$(MAKE) -C sieve
	$(MAKE) -C postsieve/checknorms
	$(MAKE) -C linalg
	$(MAKE) -C sqrt/naive
	$(MAKE) -C gf2x
	$(MAKE) -C gf2x/cantor
	$(MAKE) -C linalg/bw
	$(MAKE) -C linalg/bl	FAKE=1

clean:
	$(MAKE) -C utils		clean
	$(MAKE) -C polyselect		clean
	(cd postsieve/tifa; $(SCONS) $(SCONS_FLAGS) -c)
	$(MAKE) -C sieve		clean
	$(MAKE) -C postsieve/checknorms	clean
	$(MAKE) -C linalg		clean
	$(MAKE) -C sqrt/naive		clean
	$(MAKE) -C gf2x			clean
	$(MAKE) -C gf2x/cantor		clean
	$(MAKE) -C linalg/bw		clean
	$(MAKE) -C linalg/bl		clean

