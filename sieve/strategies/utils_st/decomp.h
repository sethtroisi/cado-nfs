#ifndef DECOMP_H
#define DECOMP_H

#include <stdio.h>

typedef struct decomp{
  int* tab;
  int len; //length of the decomposition
  double nb_elem; //number of elements which satisfy this decomposition
}decomp_t;


decomp_t*
decomp_create (int len, int* tab, double nb_elem);

void
decomp_free (decomp_t* t);

double
decomp_get_nb_elem (decomp_t* t);

int*
decomp_get_decomp (decomp_t* t);

int
decomp_get_len (decomp_t* t);

void
decomp_set_decomp (decomp_t* t, int* tab, int len);

void
decomp_set_nb_elem (decomp_t* t, double nb_elem);

void
decomp_print_file (decomp_t* t, FILE* output_file);


void
decomp_print (decomp_t* t);





#endif /* DECOMP_H */
