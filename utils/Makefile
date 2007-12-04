CFLAGS=-g -O2

names:=polyfile relation timing long_poly getprime mod_ul
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a

%.o: %.c
	gcc $(CFLAGS) -c $< -I..

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	-rm -f *~ *.bak $(obj) libutils.a

