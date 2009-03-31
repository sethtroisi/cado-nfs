.PHONY: clean all

TOP:=.

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

