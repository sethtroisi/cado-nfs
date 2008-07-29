TOP:=..
include $(TOP)/Makefile.common


# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

all: polyselect kleinjung

polyselect: polyselect.o aux.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm

kleinjung: kleinjung.o aux.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -lutils -lgmp -lm

# Some more dependencies
polyselect.o: polyselect.c aux.o $(TOP)/cado.h

kleinjung.o: kleinjung.c aux.o $(TOP)/cado.h

aux.o: aux.c aux.h $(TOP)/cado.h

clean:
	-rm -f polyselect kleinjung *~ *.o
