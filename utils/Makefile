CFLAGS=-g

names:=polyfile relation
src:=$(patsubst %,%.c,$(names))
obj:=$(patsubst %,%.o,$(names))

all: libutils.a

%.o: %.c
	gcc $(CFLAGS) -c $< -I..

libutils.a: $(obj)
	ar r libutils.a $(obj)

clean:
	rm -f $(obj) libutils.a

