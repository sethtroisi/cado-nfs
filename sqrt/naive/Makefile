TOP:=../..
include $(TOP)/Makefile.common

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

all: algsqrt

algsqrt: algsqrt.o poly.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp -lm

fmalgsqrt: fmalgsqrt.o poly.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

ratsqrt: ratsqrt.o
	$(CC) $(MY_LFLAGS) -o $@ $^ -lgmp

# Add some dependencies
algsqrt fmalgsqrt ratsqrt: $(TOP)/utils/libutils.a

clean:
	-rm -f poly algsqrt ratsqrt fmalgsqrt *.o *~

