
.PHONY: clean all

all:
	make -C utils
	make -C polyselect
	make -C sieve
	make -C linalg
	make -C sqrt/naive

clean:
	make -C utils	clean
	make -C polyselect	clean
	make -C sieve	clean
	make -C linalg	clean
	make -C sqrt/naive	clean
