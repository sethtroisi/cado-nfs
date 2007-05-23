CFLAGS=-O2 -g
LDFLAGS=/usr/lib64/libgmp.a

procrels: procrels.c
	gcc $(CFLAGS) -I.. $< $(LDFLAGS) -o $@

gauss: gauss.c
	gcc -O4 -DNDEBUG -funroll-loops -DMULTI_ROW=3 -DNB_TEST=10 $< $(LDFLAGS) -o $@


clean:
	\rm procrels gauss
