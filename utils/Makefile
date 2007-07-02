CFLAGS=-g -O2

names:=polyfile relation cputime long_poly getprime
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a

%.o: %.c
	gcc $(CFLAGS) -c $< -I..

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	rm -f $(obj) libutils.a

