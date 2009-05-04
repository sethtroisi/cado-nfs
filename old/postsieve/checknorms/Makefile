TOP:=../..

include $(TOP)/Makefile.common
include $(TOP)/Makefile.tifa-support

# We require the cado utils library
MY_LFLAGS+=-L$(TOP)/utils
MY_CFLAGS+=-I$(TOP)/utils

all: checknorms checknorms.2lp

checknorms: checknorms.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -ltifa -lutils -lgmp -lm
checknorms.2lp: checknorms.2lp.o
	$(CC) $(MY_LFLAGS) $^ -o $@ -ltifa -lutils -lgmp -lm

config.h: config_h.in
	cp config_h.in config.h

clean:
	-rm -f checknorms checknorms.2lp *~ *.o
