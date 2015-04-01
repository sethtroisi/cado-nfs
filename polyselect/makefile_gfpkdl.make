test : gfpkdl_test
	echo 'test gfpkdl polyselect for a few small p.'
	./gfpkdl_test -p 10000000019 -k 2 -label test01
	./gfpkdl_test -p 10000000019 -k 2 -label test02
	./gfpkdl_test -p 10000000019 -k 2 -label test03
	echo 'list of .poly files:'
	ls -1 *test*.poly


gfpkdl_test : gfpkdlpolyselect.o gfpkdlpolyselect_test.o
	cc -lgmp -std=c99 -o gfpkdl_test gfpkdlpolyselect_test.o gfpkdlpolyselect.o -I . -I ../ -I ../utils -I ../build/titans -I ../build/cormoran/ -I-

gfpkdlpolyselect_test.o : gfpkdlpolyselect_test.c gfpkdlpolyselect.h
	cc -lgmp -std=c99 -c gfpkdlpolyselect_test.c -I . -I ../ -I ../utils -I ../build/titans -I ../build/cormoran/ -I-

gfpkdlpolyselect.o : gfpkdlpolyselect.c gfpkdlpolyselect.h  table_t_Py_f_deg4_type0_h1_t-200--200.c
	cc -lgmp -std=c99 -c gfpkdlpolyselect.c -I . -I ../ -I ../utils -I ../build/titans/ -I ../build/cormoran/ -I-

