gfpkdl_test : gfpkdlpolyselect.o gfpkdlpolyselect_test.o
	cc -lgmp -std=c99 -o gfpkdl_test gfpkdlpolyselect_test.o gfpkdlpolyselect.o

gfpkdlpolyselect_test.o : gfpkdlpolyselect_test.c gfpkdlpolyselect.h table_t_Py_f_deg4_type0_h1_t-200--200.c
	cc -lgmp -std=c99 -c gfpkdlpolyselect_test.c  table_t_Py_f_deg4_type0_h1_t-200--200.c -I . -I ../ -I ../utils -I ../build/titans/ -I-

gfpkdlpolyselect.o : gfpkdlpolyselect.c gfpkdlpolyselect.h table_t_Py_f_deg4_type0_h1_t-200--200.c
	cc -lgmp -std=c99 -c gfpkdlpolyselect.c  table_t_Py_f_deg4_type0_h1_t-200--200.c -I . -I ../ -I ../utils -I ../build/titans/ -I-

