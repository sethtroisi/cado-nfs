.PHONY: clean all

TOP:=.

-include Makefile.local

# Currently, this stuff is equivalent to forcing ``make WITHIN_CADO=1''
# for all sub-makes. This allows the makefiles to react differently
# depending on whether they're within the cado source tree, or within
# standalone packages (relevant to gf2x).
WITHIN_CADO:=1
export WITHIN_CADO

all:
	$(MAKE) -C utils libutils.a check_rels
	$(MAKE) -C polyselect
	$(MAKE) -C sieve/ecm libfacul.a
	$(MAKE) -C sieve makefb las
	$(MAKE) -C linalg active-targets
	$(MAKE) all-gf2x
	$(MAKE) -C cantor
	$(MAKE) -C linalg/bwc
	$(MAKE) -C linalg/bwc/lingen
	$(MAKE) -C sqrt/naive algsqrt

clean:
	$(MAKE) -C utils		clean
	$(MAKE) -C polyselect		clean
	$(MAKE) -C sieve/ecm		clean
	$(MAKE) -C sieve                clean
	$(MAKE) -C linalg               clean
	$(MAKE) clean-gf2x
	$(MAKE) -C cantor		clean
	$(MAKE) -C linalg/bwc		clean
	$(MAKE) -C linalg/bwc/lingen	clean
	$(MAKE) -C sqrt/naive		clean

# Because gf2x is a separate autotools project, we need to call it in a
# special way. We disable shared libs because it's irrelevant here.

all-gf2x:
	if [ ! -f gf2x/Makefile ] ; then \
		cd gf2x ;	\
		./configure --disable-shared $(CONFIGURE_EXTRA) ;	\
	else : ; fi
	$(MAKE) -C gf2x

clean-gf2x:
	if [ ! -f gf2x/Makefile ] ; then : ; else $(MAKE) -C gf2x clean ; fi

