#ifndef SIEVING_INTERVAL_H
#define SIEVING_INTERVAL_H

#include <stdio.h>
#include <stdint.h>

typedef struct
{
  unsigned int * h;
  unsigned int t;
}  s_sieving_inverval_t;

typedef s_sieving_inverval_t sieving_interval_t[1];
typedef s_sieving_inverval_t * sieving_interval_ptr;
typedef const s_sieving_inverval_t * sieving_interval_srcptr;

/*
  Initialise a sieving interval.

  H: the sieving interval.
  t: dimension of the lattice.
*/
void sieving_interval_init(sieving_interval_ptr H, unsigned int t);

/*
  Set the ith value of a sieving interval.

  H: the sieving interval.
  i: index of the value.
  value: the value of h[i].
*/
void sieving_interval_set_hi(sieving_interval_ptr H, unsigned int i, unsigned
                             int value);

/*
  Find the number of elements contained in the sieving region defined by H and
   set the result in nb.

  nb: number of element in the sieving region.
  H: the sieving interval.
*/
void sieving_interval_number_element(uint64_t * nb, sieving_interval_srcptr H);

/*
  Delete a sieving interval.

  H: the sieving interval.
*/
void sieving_interval_clear(sieving_interval_ptr H);

/*
  To print a sieving interval.

  H: the sieving interval.
*/
void sieving_interval_printf(sieving_interval_srcptr H);

/*
  To write a sieving interval in a file.

  file: the file.
  H: the sieving interval.
*/
void sieving_interval_fprintf(FILE * file, sieving_interval_srcptr H);

#endif /* SIEVING_INTERVAL_H */
