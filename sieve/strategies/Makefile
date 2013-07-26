all: primedb buildchain

FLAGS=-Wall -std=c99 -g

primedb: primedb.c
	gcc $(FLAGS) -o primedb primedb.c -lgmp

buildchain: methods.c methods.h buildchain.c
	gcc $(FLAGS) -o buildchain methods.c buildchain.c -lm
