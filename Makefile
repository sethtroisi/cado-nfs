.PHONY: clean all

TOP:=.

# It's useful to pass these down to sub-makefiles if we do have the
# information. Some of the code base does not (yet) explicitly load
# Makefile.local becase it's got a life outside cado
export CC
export CXX
export GMP_LIBDIR
export GMP_INCDIR

all:
	$(MAKE) -C utils libutils.a check_rels
	$(MAKE) -C polyselect
	$(MAKE) -C sieve/ecm libfacul.a
	$(MAKE) -C sieve makefb las
	$(MAKE) -C linalg freerel duplicates purge merge replay transpose \
                   balance characters allsqrt apply_perm
	$(MAKE) -C gf2x
	$(MAKE) -C gf2x/cantor
	$(MAKE) -C linalg/bw
	$(MAKE) -C sqrt/naive algsqrt

clean:
	$(MAKE) -C utils		clean
	$(MAKE) -C polyselect		clean
	$(MAKE) -C sieve/ecm		clean
	$(MAKE) -C sieve                clean
	$(MAKE) -C linalg               clean
	$(MAKE) -C gf2x			clean
	$(MAKE) -C gf2x/cantor		clean
	$(MAKE) -C linalg/bw		clean
	$(MAKE) -C sqrt/naive		clean

