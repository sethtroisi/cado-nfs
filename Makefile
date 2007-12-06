
.PHONY: clean all

all:
	$(MAKE) -C utils
	$(MAKE) -C polyselect
	$(MAKE) -C sieve
	$(MAKE) -C linalg
	$(MAKE) -C sqrt/naive

clean:
	$(MAKE) -C utils	clean
	$(MAKE) -C polyselect	clean
	$(MAKE) -C sieve	clean
	$(MAKE) -C linalg	clean
	$(MAKE) -C sqrt/naive	clean
